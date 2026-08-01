"""Evolutionary weight search judged against an opponent pool.

The previous loop accepted a mutation when it scored >= 55% over six games
against the immediately previous version. Six games has a standard error near
20 points, so that gate accepted noise roughly as often as it accepted
improvements, and "better than last time" is not the same as "better".

This loop instead:

* screens a mutation cheaply against the whole pool,
* re-runs survivors at the acceptance sample size,
* accepts only on a *paired* difference that clears a significance test,
* archives every accepted model so it can rejoin the pool as a historical
  opponent.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import fields, replace
from pathlib import Path

from .arena import archive, compare, evaluate, seeds_for
from .bot import BotWeights
from .opponents import discover_history, standard_pool

DEFAULT_MODEL = Path("models/bot_weights.json")
DEFAULT_ARCHIVE = Path("models/history")

# Sample sizes. Screening is deliberately cheap and noisy; nothing is accepted
# on a screening result alone.
SCREEN_GAMES = 160
ACCEPT_GAMES = 600
RELEASE_GAMES = 2_000

# Coefficients whose sign is known a priori. Letting the search rediscover them
# wastes generations, and it previously converged on a negative camp bonus --
# a bot that avoided the safe squares.
FLOORS = {
    "flag_capture": 50.0,
    "mobility": 0.0,
    "camp": 0.0,
    "hq_pressure": 0.0,
    "hq_strike": 0.0,
    "mine_risk": 0.0,
    "engineer_mine": 0.0,
    "unknown_risk": 0.0,
    "belief_battle": 0.0,
    "eval_material": 0.1,
    "eval_mobility": 0.0,
    "eval_hq_attack": 0.0,
    "eval_hq_attack_certain": 0.0,
    "eval_hq_defense": 0.0,
    "eval_hq_defense_certain": 0.0,
}
CEILINGS = {"engineer_waste": 0.0}
FROZEN = {"noise", "eval_terminal"}


def mutate(
    weights: BotWeights,
    rng: random.Random,
    scale: float,
    coordinates: int = 4,
) -> BotWeights:
    """Perturb a few coefficients, not all of them.

    Two earlier mistakes, both of which made every candidate worse:

    * a 35% step on all 22 weights at once is a huge jump in 22 dimensions
      away from an already-decent point, so essentially nothing survived;
    * the step floor of `max(0.35, abs(value))` meant a coefficient of 0.12 got
      a step wider than itself, and the `FLOORS` clamp then pinned about half
      of those draws to exactly zero -- a systematic bias toward 0 for every
      small weight rather than a search.
    """
    names = [
        descriptor.name
        for descriptor in fields(weights)
        if descriptor.name not in FROZEN
    ]
    values: dict[str, float] = {}
    for name in rng.sample(names, min(coordinates, len(names))):
        value = getattr(weights, name)
        step = rng.gauss(0.0, scale * max(0.02, abs(value)))
        candidate = value + step
        if name in FLOORS:
            candidate = max(FLOORS[name], candidate)
        if name in CEILINGS:
            candidate = min(CEILINGS[name], candidate)
        values[name] = candidate
    values["noise"] = max(0.01, weights.noise * 0.98)
    return replace(weights, **values)


def train(
    generations: int,
    seed: int,
    output: Path,
    workers: int | None,
    screen_games: int,
    accept_games: int,
    archive_dir: Path,
    start: BotWeights | None = None,
) -> BotWeights:
    rng = random.Random(seed)
    incumbent = start or BotWeights()
    # Opponents must be identical for candidate and incumbent, so pin the
    # weight-driven ones to a fixed anchor rather than letting them be built
    # from whichever model is being measured.
    anchor = output.parent / "defaults.json"
    if not anchor.exists():
        BotWeights().save(anchor)
    pool = standard_pool(
        history=discover_history(archive_dir), anchor=str(anchor)
    )
    scale = 0.25
    accepted = 0

    for generation in range(1, generations + 1):
        candidate = mutate(incumbent, rng, scale)
        screen_seeds = seeds_for(screen_games, pool, offset=seed + generation * 1_000)
        screen = compare(candidate, incumbent, pool, screen_seeds, workers)
        print(
            f"gen {generation:02d}/{generations} screen: "
            f"{screen.mean_difference:+.4f} (p={screen.p_value:.3f}, "
            f"n={screen.candidate.games})",
            flush=True,
        )
        # Let anything not clearly worse through to the real test; screening is
        # a cost filter, not a decision.
        if screen.mean_difference <= -0.01:
            scale = max(0.04, scale * 0.9)
            continue

        accept_seeds = seeds_for(accept_games, pool, offset=seed + generation * 7_919)
        verdict = compare(candidate, incumbent, pool, accept_seeds, workers)
        print("  " + verdict.format().replace("\n", "\n  "), flush=True)
        if verdict.significant:
            incumbent = candidate
            accepted += 1
            label = f"gen{generation:03d}"
            archive(incumbent, archive_dir, label)
            incumbent.save(output)
            pool = standard_pool(
                history=discover_history(archive_dir), anchor=str(anchor)
            )
            print(f"  ACCEPTED -> {output} (archived as {label})", flush=True)
            # 1/5 success rule: widen after a hit, shrink while missing.
            scale = min(0.6, scale * 1.3)
        else:
            print("  rejected (not significant)", flush=True)
            scale = max(0.04, scale * 0.9)

    incumbent.save(output)
    print(f"\n{accepted}/{generations} accepted. Model saved: {output.resolve()}")
    return incumbent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="通过对手池自对弈训练军棋 Bot 权重"
    )
    parser.add_argument("--generations", type=int, default=12)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help=(
            "并行进程数。默认只用约三分之一核心并降低进程优先级，"
            "以免长时间评测把整台机器占死；无人使用时可手动调高"
        ),
    )
    parser.add_argument(
        "--screen-games", type=int, default=SCREEN_GAMES, help="快速筛选局数"
    )
    parser.add_argument(
        "--accept-games", type=int, default=ACCEPT_GAMES, help="接受前的验证局数"
    )
    parser.add_argument("--resume", type=Path, default=None)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    if arguments.generations < 1:
        raise SystemExit("generations 至少为 1")
    start = (
        BotWeights.load(arguments.resume)
        if arguments.resume is not None and arguments.resume.exists()
        else None
    )
    train(
        generations=arguments.generations,
        seed=arguments.seed,
        output=arguments.output,
        workers=arguments.workers,
        screen_games=arguments.screen_games,
        accept_games=arguments.accept_games,
        archive_dir=arguments.archive,
        start=start,
    )


if __name__ == "__main__":
    main()
