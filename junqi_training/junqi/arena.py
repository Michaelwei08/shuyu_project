"""Reproducible, parallel evaluation.

Design rules, all of which the previous training loop broke:

* **Identical openings.** The opening is a pure function of the match seed, so
  every agent judged on seed *k* plays the same position.
* **Colour-swapped pairs.** Each seed is played twice with the sides exchanged,
  which removes first-move and layout advantage from the comparison.
* **Common random numbers.** A player's RNG seed comes from the *side* it plays
  and the match seed, never from which weights it carries -- so two candidates
  meet byte-identical sampled hidden states.
* **A denser objective than win/loss.** Flag capture, headquarters pressure and
  flag defence all enter the score, so a 200-game screen carries far more signal
  per game than win rate alone.
"""

from __future__ import annotations

import math
import os
import random
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from statistics import fmean, stdev

from .bot import BotWeights, _distance, _piece_value
from .deployment import strategic_deployment
from .game import Game
from .opponents import AgentSpec, Pool
from .types import Owner, Piece, PieceKind, Position

MAX_MOVES = 300
HORIZON = 12.0

# Shaping stays small next to the result so that wins always dominate; it exists
# to cut the variance of a fixed-size sample, not to change what "better" means.
SHAPE_FLAG_CAPTURE = 0.10
SHAPE_HQ_PRESSURE = 0.08
SHAPE_HQ_DEFENSE = 0.08
SHAPE_MATERIAL = 0.04


def make_opening(
    seed: int, screen_caps: dict[Owner, int | None] | None = None
) -> dict[Position, Piece]:
    """Both armies deployed from the seed alone, so openings are comparable.

    Arithmetic only -- `hash()` of anything containing a str is salted per
    process, which would make openings differ between worker processes and
    quietly destroy the paired design.

    ``screen_caps`` optionally overrides how many mines each side may put on
    its flag headquarters' neighbours. It is the one deployment knob the paired
    harness can vary, and it stays a plain int per side so a `Job` remains
    picklable and every worker rebuilds the identical board.
    """
    rng = random.Random(seed * 2_654_435_761 + 12_345)
    caps = screen_caps or {}
    # One stream per side. A shared stream would make the number of draws the
    # bot's deployment consumes depend on its screen cap, which would silently
    # change the *human's* army too -- and then a capped-vs-sealed comparison
    # would be measuring a different opponent as well as a different screen.
    streams = {owner: random.Random(rng.randrange(2**32)) for owner in Owner}
    board: dict[Position, Piece] = {}
    for owner in (Owner.BOT, Owner.HUMAN):
        board.update(
            strategic_deployment(owner, streams[owner], screen_cap=caps.get(owner))
        )
    return board


@lru_cache(maxsize=32)
def _load_weights(path: str) -> BotWeights:
    return BotWeights.load(path)


@dataclass(frozen=True)
class Job:
    weights: BotWeights
    opponent: AgentSpec
    seed: int
    subject_side: int  # Owner value the candidate plays
    #: Mines the *subject* may seal its flag headquarters with. `None` uses the
    #: default policy. Only the subject's side varies, so a candidate and an
    #: incumbent still face byte-identical opposing armies on the same seed and
    #: the comparison stays paired.
    subject_screen_cap: int | None = None


@dataclass(frozen=True)
class MatchResult:
    opponent: str
    seed: int
    subject_side: int
    result: float  # 1 win / .5 draw / 0 loss, from the candidate's view
    won_by_flag: bool
    lost_by_flag: bool
    hq_pressure: float  # mean closeness of our raiders to their headquarters
    hq_defense: float  # mean distance of their raiders from ours
    material: float  # final material margin, normalised
    plies: int

    @property
    def score(self) -> float:
        return (
            self.result
            + SHAPE_FLAG_CAPTURE * (float(self.won_by_flag) - float(self.lost_by_flag))
            + SHAPE_HQ_PRESSURE * self.hq_pressure
            + SHAPE_HQ_DEFENSE * self.hq_defense
            + SHAPE_MATERIAL * self.material
        )


def _closest_raider(game: Game, side: Owner, targets: list[Position]) -> float:
    if not targets:
        return HORIZON
    from .board import HEADQUARTERS

    return min(
        (
            float(_distance(position, target))
            for position, piece in game.board.items()
            if piece.owner == side
            and piece.kind.movable
            and position not in HEADQUARTERS
            for target in targets
        ),
        default=HORIZON,
    )


def _material_margin(game: Game, side: Owner) -> float:
    total = sum(
        _piece_value(piece.kind) * (1 if piece.owner == side else -1)
        for piece in game.board.values()
    )
    return max(-1.0, min(1.0, total / 120.0))


def play_match(job: Job) -> MatchResult:
    subject_side = Owner(job.subject_side)
    game = Game(
        board=make_opening(job.seed, {subject_side: job.subject_screen_cap}),
        turn=Owner(job.seed % 2),
    )
    opponent_weights = (
        _load_weights(job.opponent.weights_path)
        if job.opponent.weights_path
        else job.weights
    )
    from .search_bot import SearchBot

    subject = SearchBot(
        job.weights,
        seed=job.seed * 2 + int(subject_side),
        samples=3,
        beam_width=8,
        reply_width=4,
    )
    challenger = job.opponent.build(
        opponent_weights, seed=job.seed * 2 + int(subject_side.other)
    )
    players = {subject_side: subject, subject_side.other: challenger}

    pressure: list[float] = []
    defense: list[float] = []
    while not game.over and game.move_count < MAX_MOVES:
        game.apply(players[game.turn].choose_move(game))
        pressure.append(
            1.0
            - _closest_raider(
                game, subject_side, game.flag_candidates(subject_side.other)
            )
            / HORIZON
        )
        defense.append(
            _closest_raider(
                game, subject_side.other, game.flag_candidates(subject_side)
            )
            / HORIZON
        )

    flag_taken = _flag_was_captured(game)
    if game.winner == subject_side:
        result = 1.0
    elif game.winner is None:
        result = 0.5
    else:
        result = 0.0
    return MatchResult(
        opponent=job.opponent.name,
        seed=job.seed,
        subject_side=int(subject_side),
        result=result,
        won_by_flag=flag_taken and game.winner == subject_side,
        lost_by_flag=flag_taken and game.winner == subject_side.other,
        hq_pressure=fmean(pressure) if pressure else 0.0,
        hq_defense=fmean(defense) if defense else 0.0,
        material=_material_margin(game, subject_side),
        plies=game.move_count,
    )


def _flag_was_captured(game: Game) -> bool:
    if not game.records:
        return False
    last = game.records[-1]
    return last.defender is not None and last.defender.kind == PieceKind.FLAG


def build_jobs(
    weights: BotWeights,
    pool: Pool,
    seeds: list[int],
    subject_screen_cap: int | None = None,
) -> list[Job]:
    """Every seed against every opponent, both colours -- the paired design."""
    return [
        Job(weights, spec, seed, side, subject_screen_cap)
        for seed in seeds
        for spec in pool.specs
        for side in (int(Owner.HUMAN), int(Owner.BOT))
    ]


def seeds_for(games: int, pool: Pool, offset: int = 0) -> list[int]:
    per_seed = max(1, len(pool) * 2)
    count = max(1, math.ceil(games / per_seed))
    return [offset + index for index in range(count)]


def default_workers() -> int:
    """Deliberately leaves most of the machine alone.

    This is a desktop, not a cluster node: saturating every core makes the
    whole system unusable for the hours an evaluation run takes. Raise it with
    `--workers` when nobody is using the machine.
    """
    return max(1, (os.cpu_count() or 2) // 3)


def _deprioritise() -> None:
    """Drop worker priority so the UI keeps getting scheduled."""
    try:
        if sys.platform == "win32":
            import ctypes

            below_normal = 0x00004000
            kernel32 = ctypes.windll.kernel32
            kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), below_normal)
        else:
            os.nice(10)
    except Exception:  # pragma: no cover - best effort, never fatal
        pass


MAX_FAILURE_RATE = 0.02


def _report_progress(futures: list, workers: int) -> None:
    """Print a heartbeat as matches land.

    A 4000-game comparison takes minutes and used to print nothing until it
    finished, which looks indistinguishable from a hang.
    """
    total = len(futures)
    if total < 200:
        return
    step = max(1, total // 20)
    state = {"done": 0, "started": time.perf_counter()}
    lock = threading.Lock()

    def tick(_future) -> None:
        with lock:
            state["done"] += 1
            done = state["done"]
        if done % step and done != total:
            return
        elapsed = time.perf_counter() - state["started"]
        rate = done / elapsed if elapsed else 0.0
        remaining = (total - done) / rate if rate else 0.0
        print(
            f"    {done}/{total} matches ({done / total:.0%}) "
            f"{rate:.1f}/s, ~{remaining:.0f}s left",
            flush=True,
        )

    for future in futures:
        future.add_done_callback(tick)


def run_jobs(
    jobs: list[Job], workers: int | None = None
) -> list[MatchResult | None]:
    """Play every job, tolerating isolated failures.

    A single raised game must not throw away hours of evaluation, so a failed
    job becomes `None` and its pair is dropped downstream. Anything above
    `MAX_FAILURE_RATE` is a bug rather than a hiccup and raises.
    """
    count = default_workers() if workers is None else workers
    results: list[MatchResult | None] = []
    failures: list[str] = []

    if count <= 1:
        for job in jobs:
            try:
                results.append(play_match(job))
            except Exception as error:  # noqa: BLE001 - recorded, then re-raised in bulk
                failures.append(f"{job.opponent.name}/seed{job.seed}: {error!r}")
                results.append(None)
    else:
        with ProcessPoolExecutor(
            max_workers=count, initializer=_deprioritise
        ) as executor:
            futures = [executor.submit(play_match, job) for job in jobs]
            _report_progress(futures, count)
            for job, future in zip(jobs, futures, strict=True):
                try:
                    results.append(future.result())
                except Exception as error:  # noqa: BLE001
                    failures.append(f"{job.opponent.name}/seed{job.seed}: {error!r}")
                    results.append(None)

    if failures and len(failures) > MAX_FAILURE_RATE * max(1, len(jobs)):
        raise RuntimeError(
            f"{len(failures)}/{len(jobs)} matches failed, e.g. {failures[:3]}"
        )
    if failures:
        print(f"  warning: {len(failures)}/{len(jobs)} matches failed", flush=True)
    return results


@dataclass
class PoolReport:
    results: list[MatchResult]

    @property
    def games(self) -> int:
        return len(self.results)

    @property
    def win_rate(self) -> float:
        return fmean(result.result for result in self.results)

    @property
    def score(self) -> float:
        return fmean(result.score for result in self.results)

    def by_opponent(self) -> dict[str, tuple[int, float, float]]:
        buckets: dict[str, list[MatchResult]] = {}
        for result in self.results:
            buckets.setdefault(result.opponent, []).append(result)
        return {
            name: (
                len(group),
                fmean(item.result for item in group),
                fmean(item.score for item in group),
            )
            for name, group in sorted(buckets.items())
        }

    def format(self) -> str:
        lines = [
            f"{'opponent':<18} {'games':>6} {'win%':>8} {'score':>8}",
            "-" * 44,
        ]
        for name, (count, win, score) in self.by_opponent().items():
            lines.append(f"{name:<18} {count:>6} {win:>7.1%} {score:>8.3f}")
        lines.append("-" * 44)
        lines.append(
            f"{'AGGREGATE':<18} {self.games:>6} {self.win_rate:>7.1%} "
            f"{self.score:>8.3f}"
        )
        return "\n".join(lines)


def evaluate(
    weights: BotWeights,
    pool: Pool,
    seeds: list[int],
    workers: int | None = None,
    subject_screen_cap: int | None = None,
) -> PoolReport:
    played = run_jobs(build_jobs(weights, pool, seeds, subject_screen_cap), workers)
    return PoolReport([result for result in played if result is not None])


@dataclass
class Comparison:
    candidate: PoolReport
    incumbent: PoolReport
    mean_difference: float
    standard_error: float
    p_value: float

    @property
    def significant(self) -> bool:
        return self.p_value < 0.05 and self.mean_difference > 0

    @property
    def detectable(self) -> float:
        """Smallest true improvement this sample size could have found.

        Worth printing: if it is larger than any plausible gain, a rejection
        says nothing about the candidate, only about the sample size.
        """
        return 1.645 * self.standard_error

    def format(self) -> str:
        return (
            f"candidate {self.candidate.score:.3f} (win {self.candidate.win_rate:.1%}) "
            f"vs incumbent {self.incumbent.score:.3f} "
            f"(win {self.incumbent.win_rate:.1%})\n"
            f"paired difference {self.mean_difference:+.4f} "
            f"+/- {self.standard_error:.4f} (SE), p = {self.p_value:.4f}, "
            f"n = {self.candidate.games} games each\n"
            f"minimum detectable improvement at this n: {self.detectable:+.4f}"
        )


def compare(
    candidate: BotWeights,
    incumbent: BotWeights,
    pool: Pool,
    seeds: list[int],
    workers: int | None = None,
    candidate_screen_cap: int | None = None,
    incumbent_screen_cap: int | None = None,
) -> Comparison:
    """Run both weight sets over an identical job list and pair the outcomes.

    The two ``screen_cap`` arguments are the one way a *non-weight* change can
    enter this comparison: they alter only the subject's own deployment, so the
    opposing army on a given seed is unchanged and the pairing survives.
    """
    candidate_jobs = build_jobs(candidate, pool, seeds, candidate_screen_cap)
    incumbent_jobs = build_jobs(incumbent, pool, seeds, incumbent_screen_cap)
    combined = run_jobs(candidate_jobs + incumbent_jobs, workers)
    split = len(candidate_jobs)
    # Drop a pair entirely if either half failed -- a half-pair would break the
    # pairing that makes this comparison low-variance in the first place.
    pairs = [
        (a, b)
        for a, b in zip(combined[:split], combined[split:], strict=True)
        if a is not None and b is not None
    ]
    if not pairs:
        raise RuntimeError("every paired match failed")
    left = [a for a, _ in pairs]
    right = [b for _, b in pairs]
    differences = [a.score - b.score for a, b in pairs]
    mean = fmean(differences)
    if len(differences) > 1:
        deviation = stdev(differences)
        error = deviation / math.sqrt(len(differences))
    else:
        error = 0.0
    if error == 0.0:
        p_value = 0.0 if mean > 0 else 1.0
    else:
        z = mean / error
        p_value = 1.0 - 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
    return Comparison(PoolReport(left), PoolReport(right), mean, error, p_value)


def archive(weights: BotWeights, directory: str | Path, label: str) -> Path:
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=True)
    destination = folder / f"{label}.json"
    weights.save(destination)
    return destination
