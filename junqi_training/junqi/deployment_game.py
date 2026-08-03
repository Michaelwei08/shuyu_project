"""The opening as a normal-form game, and its equilibrium.

Deployment is chosen simultaneously, once, before either side has seen
anything -- so it is a one-shot matrix game, and its solution is a
*distribution* over openings rather than a single best one. Everything else in
this project tunes a policy against a fixed pool; this is the one part with an
actual game-theoretic answer, and it is small enough to solve exactly in the
standard library.

Why it is worth doing at all: on 2026-08-01 a single bit of the deployment
generator (`SCREEN_MINE_CAP`, 3 vs 2) measured **-0.0985 +/- 0.0217 over 806
paired games** -- about as large as total omniscience, and larger than every
weight-space result this project has ever produced. Deployment is also the
least-explored thing here: `strategic_deployment` was a point generator with
hardcoded constants.

Two separate questions, which this module keeps apart on purpose:

* **Which family is best against a fixed opponent?** A pool measurement, and
  the answer is a *pure* strategy -- against an opponent that never adapts, no
  mixture can beat the best pure response. `junqi.benchmark --deployment` is
  the tool for that one.
* **What should you play against an opponent who is also choosing?** That is
  this module. The answer is a mixture, and it is the only thing here that
  survives an opponent who watches you.

The matrix is measured on **win rate, not `MatchResult.score`**. Score's
shaping terms are written from the subject's seat and are not recorded for the
other one, so it is not zero-sum and `A[j][i] = 1 - A[i][j]` would be false.

Antisymmetry is also why only the upper triangle is played. Cell (i, j) is the
seed played twice, once with the subject on each side; cell (j, i) is *the same
two games* with the seats read the other way round. Playing both would be
duplicated work reported as independent evidence. The diagonal is measured
anyway, even though it is 0.5 by construction, because a diagonal that comes
back off 0.5 means the harness is broken -- the same role
`test_identical_weights_produce_a_zero_paired_difference` plays for `compare`.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from .arena import DEFAULT_SEARCH, Job, run_jobs
from .bot import BotWeights
from .deployment import FAMILIES
from .opponents import AgentSpec
from .types import Owner

#: The subject's own policy on the other side of the board. A true mirror --
#: same weights, same search budget as `DEFAULT_SEARCH` -- so that the only
#: thing differing between the two seats is the opening, which is the whole
#: point of the matrix.
def mirror_spec(model: str) -> AgentSpec:
    samples, beam_width, reply_width = DEFAULT_SEARCH
    return AgentSpec(
        "mirror",
        "search",
        weights_path=model,
        samples=samples,
        beam_width=beam_width,
        reply_width=reply_width,
    )


def cell_jobs(
    weights: BotWeights,
    spec: AgentSpec,
    row: str,
    column: str,
    seeds: list[int],
) -> list[Job]:
    """One matrix cell: every seed, both colours, row family vs column family."""
    return [
        Job(
            weights,
            spec,
            seed,
            side,
            subject_deployment=row,
            opponent_deployment=column,
        )
        for seed in seeds
        for side in (int(Owner.HUMAN), int(Owner.BOT))
    ]


def measure_matrix(
    families: list[str],
    weights: BotWeights,
    model: str,
    games_per_cell: int,
    workers: int | None,
    offset: int = 900_000,
) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], int]]:
    """Win rate of every row family against every column family.

    Every cell uses the *same* seed list, so the openings underlying two cells
    differ only in the family that generated them -- common random numbers, the
    same discipline `compare()` uses.
    """
    spec = mirror_spec(model)
    seeds = [offset + index for index in range(max(1, games_per_cell // 2))]
    # A diagonal cell is the *same game* twice -- same board, same turn, same
    # per-seat RNG seeds, only the seat the subject is called from differs -- so
    # it returns exactly 0.5 however many seeds it is given. Four is enough to
    # assert that, and the rest of the budget goes where there is signal.
    check_seeds = seeds[:2]
    pairs = [
        (row, column)
        for index, row in enumerate(families)
        for column in families[index:]
    ]
    jobs: list[Job] = []
    for row, column in pairs:
        jobs.extend(
            cell_jobs(
                weights, spec, row, column, check_seeds if row == column else seeds
            )
        )

    off_diagonal = sum(1 for row, column in pairs if row != column)
    print(
        f"{off_diagonal} cells x {2 * len(seeds)} games "
        f"+ {len(families)} mirror checks = {len(jobs)} games "
        f"({len(families)} families, upper triangle only)",
        flush=True,
    )
    started = time.perf_counter()
    played = run_jobs(jobs, workers)
    elapsed = time.perf_counter() - started
    print(f"{len(jobs)} games in {elapsed:.0f}s", flush=True)

    matrix: dict[tuple[str, str], float] = {}
    counts: dict[tuple[str, str], int] = {}
    cursor = 0
    for row, column in pairs:
        size = 2 * len(check_seeds if row == column else seeds)
        window = played[cursor : cursor + size]
        cursor += size
        results = [item.result for item in window if item is not None]
        if not results:
            raise RuntimeError(f"cell {row} vs {column} produced no games")
        rate = sum(results) / len(results)
        if row == column and abs(rate - 0.5) > 1e-9:
            raise RuntimeError(
                f"mirror check failed: {row} against itself scored {rate:.4f}, "
                "not 0.5 -- the two seats are no longer playing the same game"
            )
        matrix[(row, column)] = rate
        counts[(row, column)] = len(results)
        if row != column:
            # The same games, read from the other seat.
            matrix[(column, row)] = 1.0 - rate
            counts[(column, row)] = len(results)
    return matrix, counts


def regret_matching(
    families: list[str],
    matrix: dict[tuple[str, str], float],
    iterations: int = 200_000,
) -> list[float]:
    """Nash of a symmetric zero-sum matrix game, in the standard library.

    Hart & Mas-Colell regret matching. For a symmetric zero-sum game the
    time-average of self-play under regret matching converges to a Nash
    equilibrium, so one strategy vector answers for both seats.
    """
    size = len(families)
    payoff = [
        [matrix[(families[i], families[j])] - 0.5 for j in range(size)]
        for i in range(size)
    ]
    regret = [0.0] * size
    total = [0.0] * size
    for _ in range(iterations):
        positive = [value if value > 0 else 0.0 for value in regret]
        mass = sum(positive)
        strategy = (
            [value / mass for value in positive]
            if mass > 0
            else [1.0 / size] * size
        )
        for index, share in enumerate(strategy):
            total[index] += share
        utility = [
            sum(payoff[i][j] * strategy[j] for j in range(size)) for i in range(size)
        ]
        value = sum(strategy[i] * utility[i] for i in range(size))
        for index in range(size):
            regret[index] += utility[index] - value
    mass = sum(total)
    return [share / mass for share in total]


def exploitability(
    families: list[str], matrix: dict[tuple[str, str], float], strategy: list[float]
) -> tuple[float, str]:
    """How much a best responder gains against `strategy`, and what it plays.

    Zero for an exact equilibrium. This is the number that says whether mixing
    was worth anything: a pure strategy in a game with no pure Nash has a
    positive best-response gap, and the mixture's job is to close it.
    """
    size = len(families)
    gains = [
        sum(
            (matrix[(families[i], families[j])] - 0.5) * strategy[j]
            for j in range(size)
        )
        for i in range(size)
    ]
    best = max(range(size), key=lambda index: gains[index])
    return gains[best], families[best]


#: Matches sharing one seed share one opening, so they are not independent.
#: Measured at 1.26 over 806 paired games on 2026-08-01; the same correction
#: applies here, and it enters the bootstrap as a reduction in effective n.
DESIGN_EFFECT = 1.26


def bootstrap_equilibrium(
    families: list[str],
    matrix: dict[tuple[str, str], float],
    counts: dict[tuple[str, str], int],
    draws: int = 400,
    seed: int = 20260803,
    design_effect: float = DESIGN_EFFECT,
) -> tuple[list[float], list[float], list[float]]:
    """Mean and 5th/95th percentile share for each family, over cell noise.

    An equilibrium read off a matrix whose cells carry a +/-3 point standard
    error can be an artefact of that error. Resampling each cell and re-solving
    says how much of the mixture is real. This is the same question the 12%
    design-effect correction asks of a p-value, applied to a solve instead of a
    difference -- and it carries the same correction, because the games in a
    cell share openings rather than being 200 independent draws.
    """
    import random as _random

    rng = _random.Random(seed)
    samples: list[list[float]] = []
    for _ in range(draws):
        resampled: dict[tuple[str, str], float] = {}
        for index, row in enumerate(families):
            for column in families[index:]:
                if row == column:
                    resampled[(row, column)] = 0.5
                    continue
                games = max(1, round(counts[(row, column)] / design_effect))
                rate = matrix[(row, column)]
                wins = sum(1 for _ in range(games) if rng.random() < rate)
                value = wins / games
                resampled[(row, column)] = value
                resampled[(column, row)] = 1.0 - value
        samples.append(regret_matching(families, resampled, iterations=20_000))

    size = len(families)
    means = [sum(sample[i] for sample in samples) / draws for i in range(size)]
    low, high = [], []
    for i in range(size):
        column = sorted(sample[i] for sample in samples)
        low.append(column[int(0.05 * draws)])
        high.append(column[min(draws - 1, int(0.95 * draws))])
    return means, low, high


def format_matrix(
    families: list[str],
    matrix: dict[tuple[str, str], float],
    counts: dict[tuple[str, str], int],
) -> str:
    width = max(len(name) for name in families) + 1
    header = " " * width + "".join(f"{name:>12}" for name in families)
    lines = [header, "-" * len(header)]
    for row in families:
        cells = "".join(f"{matrix[(row, column)]:>11.1%} " for column in families)
        lines.append(f"{row:<{width}}{cells}")
    lines.append("")
    off_diagonal = [
        count for (row, column), count in counts.items() if row != column
    ]
    sample = min(off_diagonal) if off_diagonal else 0
    lines.append(
        f"games per off-diagonal cell: {sample} "
        f"(upper triangle measured, lower mirrored; diagonal is 0.5 by construction)"
    )
    lines.append(f"standard error of one cell: +/-{0.5 / max(1, sample) ** 0.5:.3f}")
    return "\n".join(lines)


def report(
    families: list[str],
    matrix: dict[tuple[str, str], float],
    counts: dict[tuple[str, str], int],
    design_effect: float = DESIGN_EFFECT,
) -> tuple[list[float], float]:
    """Print the matrix, the equilibrium, and how much of it survives noise."""
    print(format_matrix(families, matrix, counts))

    strategy = regret_matching(families, matrix)
    gap, responder = exploitability(families, matrix, strategy)
    print("\nequilibrium mixture:")
    for name, share in sorted(zip(families, strategy), key=lambda item: -item[1]):
        if share > 1e-4:
            print(f"  {name:<12} {share:>7.1%}")
    print(f"\nexploitability of the mixture: {gap:+.4f} (best response: {responder})")
    print("how much a best responder gains against each pure family:")
    for name in families:
        pure = [1.0 if item == name else 0.0 for item in families]
        pure_gap, pure_responder = exploitability(families, matrix, pure)
        print(f"  pure {name:<12} {pure_gap:+.4f} (best response: {pure_responder})")

    means, low, high = bootstrap_equilibrium(
        families, matrix, counts, design_effect=design_effect
    )
    print(
        f"\nbootstrap over cell noise (400 resamples, design effect "
        f"{design_effect:.2f}, 5th-95th percentile):"
    )
    for name, mean, lower, upper in sorted(
        zip(families, means, low, high), key=lambda item: -item[1]
    ):
        print(f"  {name:<12} {mean:>7.1%}   [{lower:>6.1%}, {upper:>6.1%}]")

    mix = ",".join(
        f"{name}={share:.4f}"
        for name, share in zip(families, strategy)
        if share > 1e-4
    )
    print(f"\nvalidate it against the pool with:\n  --deployment mix:{mix}")
    return strategy, gap


def main() -> None:
    parser = argparse.ArgumentParser(description="布阵阶段的正规形博弈与均衡")
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument(
        "--games-per-cell",
        type=int,
        default=120,
        help="每个矩阵格的对局数（两种先后手各半）",
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument(
        "--families",
        nargs="+",
        default=None,
        help=f"参与博弈的布阵族，默认全部：{sorted(FAMILIES)}",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("models/deployment_matrix.json"),
        help="矩阵与均衡的落盘位置",
    )
    parser.add_argument(
        "--analyse",
        type=Path,
        default=None,
        help="不重新对局，直接从已保存的矩阵重算均衡与自助区间",
    )
    parser.add_argument(
        "--design-effect",
        type=float,
        default=DESIGN_EFFECT,
        help="同一开局下多局的相关性修正，用于自助区间（默认 1.26，实测值）",
    )
    arguments = parser.parse_args()

    if arguments.analyse is not None:
        saved = json.loads(arguments.analyse.read_text(encoding="utf-8"))
        families = saved["families"]
        matrix = {
            (key.split("|")[0], key.split("|")[1]): value
            for key, value in saved["matrix"].items()
        }
        counts = {
            key: saved["games_per_cell"] for key in matrix if key[0] != key[1]
        }
        print(f"loaded {arguments.analyse} ({saved['games_per_cell']} games per cell)")
        print()
        report(families, matrix, counts, arguments.design_effect)
        return

    families = arguments.families or list(FAMILIES)
    weights = (
        BotWeights.load(arguments.model)
        if arguments.model.exists()
        else BotWeights()
    )
    print(f"families: {families}")
    matrix, counts = measure_matrix(
        families,
        weights,
        str(arguments.model),
        arguments.games_per_cell,
        arguments.workers,
    )
    print()
    strategy, gap = report(families, matrix, counts, arguments.design_effect)

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    off_diagonal = [count for (row, column), count in counts.items() if row != column]
    arguments.out.write_text(
        json.dumps(
            {
                "families": families,
                "games_per_cell": min(off_diagonal),
                "matrix": {
                    f"{row}|{column}": value for (row, column), value in matrix.items()
                },
                "equilibrium": dict(zip(families, strategy)),
                "exploitability": gap,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\nwritten to {arguments.out}")


if __name__ == "__main__":
    main()
