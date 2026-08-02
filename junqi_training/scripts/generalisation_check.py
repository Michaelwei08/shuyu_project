"""Does a fitted model beat opponents it never trained against?

Mandatory before adopting anything fit to self-play data. `junqi.training`'s
evolutionary search moves four coefficients per generation behind a paired
p-value, which is low-capacity enough that fitting the pool and being judged by
the pool is survivable. A 96-parameter value function is not: on 2026-08-01 the
first fitted evaluation scored +0.1158 (p = 0.001) against the pool it trained
on and **-0.0013 (p = 0.52)** against four opponents held out of its training,
with everything else -- weights, seeds, opponents, sample size -- identical. The
whole effect was the benchmark grading its own homework.

Usage: retrain with those opponents excluded, then measure against exactly them.

    python -m junqi.value_training --games 500 --workers N \\
        --exclude search-deep search-mid selective selective-strict
    python scripts/generalisation_check.py --workers N

A change that survives this is worth a full acceptance run. One that does not is
worth nothing, however good its in-pool number looks.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from junqi.arena import compare, seeds_for  # noqa: E402
from junqi.bot import BotWeights  # noqa: E402
from junqi.opponents import Pool, standard_pool  # noqa: E402

DEFAULT_HELD_OUT = ("search-deep", "search-mid", "selective", "selective-strict")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--anchor", type=Path, default=Path("models/defaults.json"))
    parser.add_argument("--games", type=int, default=400)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--offset", type=int, default=555_000)
    parser.add_argument(
        "--held-out", nargs="*", default=list(DEFAULT_HELD_OUT),
        help="必须与 value_training 的 --exclude 完全一致",
    )
    parser.add_argument(
        "--scales", type=float, nargs="*", default=[15.0, 35.0, 70.0],
        help="要检验的 eval_value_scale 取值",
    )
    arguments = parser.parse_args()

    held_out = set(arguments.held_out)
    full = standard_pool(anchor=str(arguments.anchor))
    pool = Pool([spec for spec in full.specs if spec.name in held_out])
    missing = held_out - {spec.name for spec in pool.specs}
    if missing:
        raise SystemExit(f"not in the pool: {sorted(missing)}")

    print(f"held-out opponents: {[spec.name for spec in pool.specs]}")
    print("the value model must NOT have trained against these\n")
    shipped = BotWeights.load(arguments.model)
    # Derive the baseline by zeroing the scale rather than using the model as
    # it sits on disk. Once the value term was adopted into bot_weights.json,
    # `shipped` already carried scale 15, so comparing against it pitted the
    # model against itself and printed "+0.0000 +/- 0.0000, p = 1.0000" -- which
    # reads exactly like a clean negative result and is not one.
    baseline = replace(shipped, eval_value_scale=0.0)
    seeds = seeds_for(arguments.games, pool, offset=arguments.offset)

    for scale in arguments.scales:
        if scale == 0.0:
            print(f"  eval_value_scale {scale:>6.1f}: skipped -- identical to "
                  f"the baseline, nothing to compare")
            continue
        verdict = compare(
            replace(shipped, eval_value_scale=scale),
            baseline,
            pool,
            seeds,
            arguments.workers,
        )
        if verdict.standard_error == 0.0:
            raise SystemExit(
                "paired difference has zero variance: the two sides are the "
                "same player. Check that --model and the scale actually differ."
            )
        flag = "PASS" if verdict.significant else "no effect"
        print(
            f"  eval_value_scale {scale:>6.1f}: {verdict.mean_difference:+.4f} "
            f"+/- {verdict.standard_error:.4f}, p = {verdict.p_value:.4f}, "
            f"n = {verdict.candidate.games}   [{flag}]"
        )


if __name__ == "__main__":
    main()
