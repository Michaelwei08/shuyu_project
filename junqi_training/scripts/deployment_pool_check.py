"""Is a deployment family worth points against the pool?

`junqi.deployment_game` answers a different question: what to deploy against an
opponent who is *also* choosing, solved as a normal-form game over mirror
matches. This script answers the practical one -- does any of that beat the
shipped generator against the actual judging pool.

The two must stay separate. Against a pool that never adapts, no mixture can
beat the best pure family, so a mixture doing nothing here is expected and is
not evidence against the equilibrium. What this measures is whether the
*shipped point* is the wrong point.

Selection here is honest by construction: the equilibrium is chosen on
self-play and judged on the pool, so the pool is genuinely held out from the
choice. But if you pick a family by reading the aggregate column below, you
have fitted to the pool -- which is why the held-out four are reported
separately, exactly as `generalisation_check.py` does. Adopt on the aggregate
(D012); use the held-out column to check the aggregate is not an artefact.

    python scripts/deployment_pool_check.py --games 1600 --workers 46
    python scripts/deployment_pool_check.py --games 2400 --workers 46 \\
        --from-matrix models/deployment_matrix.json

Design notes, all of which this shares with `arena.compare`:

* Every candidate plays the *same* job list as the baseline -- same seeds, same
  opponents, same colours, same hidden-state sampling -- differing only in the
  subject's own deployment. The opposing army on a given seed is untouched,
  because `make_opening` draws each side from its own stream.
* The baseline is played **once** and reused for every candidate. Running
  `compare()` per family would replay it each time, which costs a game list per
  candidate and buys nothing.
* Standard errors are reported both naive and clustered by seed. The naive one
  is optimistic by about 12% here, because the 32 matches sharing a seed share
  one opening -- measured design effect 1.26 on 2026-08-01.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from statistics import fmean, stdev

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from junqi.arena import build_jobs, run_jobs, seeds_for  # noqa: E402
from junqi.bot import BotWeights  # noqa: E402
from junqi.deployment import FAMILIES  # noqa: E402
from junqi.opponents import standard_pool  # noqa: E402

BASELINE = "standard"
#: The same four `generalisation_check.py` holds out, so the two scripts'
#: numbers are comparable. They are the pool's strongest non-cheating members.
HELD_OUT = ("search-deep", "search-mid", "selective", "selective-strict")


def paired(differences: list[float]) -> tuple[float, float, float]:
    """Mean, standard error and one-sided p for a list of paired differences."""
    mean = fmean(differences)
    if len(differences) < 2:
        return mean, 0.0, 1.0
    error = stdev(differences) / math.sqrt(len(differences))
    if error == 0.0:
        return mean, 0.0, 0.0 if mean > 0 else 1.0
    z = mean / error
    return mean, error, 1.0 - 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def clustered_error(groups: dict[int, list[float]]) -> float:
    """Standard error treating each seed, not each match, as the unit.

    The 32 matches sharing an opening are not 32 independent observations.
    """
    means = [fmean(values) for values in groups.values() if values]
    if len(means) < 2:
        return 0.0
    return stdev(means) / math.sqrt(len(means))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--anchor", type=Path, default=Path("models/defaults.json"))
    parser.add_argument("--games", type=int, default=1600, help="每个候选的对局数")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--offset", type=int, default=770_000)
    parser.add_argument(
        "--families",
        nargs="*",
        default=None,
        help=f"要检验的布阵族，默认除基线外全部：{sorted(FAMILIES)}",
    )
    parser.add_argument(
        "--from-matrix",
        type=Path,
        default=None,
        help="把 deployment_game 解出的均衡混合也加入候选",
    )
    arguments = parser.parse_args()
    # No `--no-history`: `standard_pool` takes no history unless it is given
    # some, and the archived models predate the current coefficients, so they
    # would not be the opponents they were trained as.

    candidates = list(
        arguments.families
        if arguments.families is not None
        else [name for name in FAMILIES if name != BASELINE]
    )
    if arguments.from_matrix is not None:
        saved = json.loads(arguments.from_matrix.read_text(encoding="utf-8"))
        mix = ",".join(
            f"{name}={share:.4f}"
            for name, share in saved["equilibrium"].items()
            if share > 1e-4
        )
        if len(saved["equilibrium"]) and mix:
            candidates.append(f"mix:{mix}")

    weights = (
        BotWeights.load(arguments.model)
        if arguments.model.exists()
        else BotWeights()
    )
    pool = standard_pool(anchor=str(arguments.anchor) if arguments.anchor else None)
    seeds = seeds_for(arguments.games, pool, offset=arguments.offset)
    print(f"anchor: {arguments.anchor}")
    print(f"pool:   {len(pool)} opponents")
    print(f"seeds:  {len(seeds)} x {2 * len(pool)} = {len(seeds) * 2 * len(pool)} games "
          f"per arm")
    print(f"arms:   {BASELINE} (baseline) + {len(candidates)} candidates")
    for name in candidates:
        print(f"          {name}")

    arms = [BASELINE, *candidates]
    jobs = []
    for name in arms:
        jobs.extend(build_jobs(weights, pool, seeds, subject_deployment=name))
    per_arm = len(jobs) // len(arms)
    print(f"total:  {len(jobs)} games\n", flush=True)

    started = time.perf_counter()
    played = run_jobs(jobs, arguments.workers)
    print(f"\n{len(jobs)} games in {time.perf_counter() - started:.0f}s\n")

    windows = {
        name: played[index * per_arm : (index + 1) * per_arm]
        for index, name in enumerate(arms)
    }
    base = windows[BASELINE]

    header = (
        f"{'family':<28} {'win%':>7} {'vs standard':>13} {'SE':>8} "
        f"{'clustered':>10} {'p':>8} {'held-out':>10}"
    )
    print(header)
    print("-" * len(header))
    baseline_win = fmean(item.result for item in base if item is not None)
    print(f"{BASELINE + ' (baseline)':<28} {baseline_win:>6.1%} {'--':>13} "
          f"{'--':>8} {'--':>10} {'--':>8} {'--':>10}")

    rows = []
    for name in candidates:
        window = windows[name]
        differences: list[float] = []
        by_seed: dict[int, list[float]] = defaultdict(list)
        held: list[float] = []
        wins: list[float] = []
        for candidate, incumbent in zip(window, base, strict=True):
            if candidate is None or incumbent is None:
                continue  # a half pair would break the pairing
            delta = candidate.score - incumbent.score
            differences.append(delta)
            by_seed[candidate.seed].append(delta)
            wins.append(candidate.result)
            if candidate.opponent in HELD_OUT:
                held.append(delta)
        if not differences:
            print(f"{name:<28} every paired match failed")
            continue
        mean, error, p_value = paired(differences)
        cluster = clustered_error(by_seed)
        held_mean = fmean(held) if held else float("nan")
        rows.append((name, mean, error, cluster, p_value, held_mean, len(differences)))
        print(
            f"{name:<28} {fmean(wins):>6.1%} {mean:>+13.4f} {error:>8.4f} "
            f"{cluster:>10.4f} {p_value:>8.4f} {held_mean:>+10.4f}"
        )

    if not rows:
        raise SystemExit("no candidate produced a usable pair")
    worst_cluster = max(row[3] for row in rows)
    worst_naive = max(row[2] for row in rows)
    # One seed cannot produce a between-seed standard deviation, and a run that
    # small has no business reporting a detectable effect at all.
    basis, error = (
        ("clustered SE", worst_cluster)
        if worst_cluster > 0
        else ("naive SE -- too few seeds to cluster", worst_naive)
    )
    print(
        f"\nn = {rows[0][6]} paired games per candidate over {len(seeds)} seeds "
        f"({sum(1 for spec in pool.specs if spec.name in HELD_OUT)} of "
        f"{len(pool)} opponents held out). "
        f"Minimum detectable improvement: {1.645 * error:+.4f} ({basis})."
    )
    print(
        "A positive aggregate that the held-out column does not echo is the "
        "pool grading its own homework; a held-out effect that the aggregate "
        "dilutes away is real but too small to adopt (D012)."
    )


if __name__ == "__main__":
    main()
