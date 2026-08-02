"""How much does the learned evaluation improve with more self-play data?

    python scripts/value_scaling.py --sizes 1600 3200 6400 --workers 46

This exists because data, not capacity, turned out to be the binding constraint.
A 500-game fit scored **-0.0013 (p = 0.52)** against opponents held out of its
training; the same 96 parameters at 2000 games scored **+0.0631 (p = 0.019)**.
Meanwhile a capacity ladder on the same pipeline showed 30 -> 96 parameters
buying nothing on held-out AUC (0.776 -> 0.771), with train AUC *below* test AUC
at every rung -- underfitting, not overfitting. So the question worth compute is
where the data curve flattens, not how many parameters to add.

For each size this collects, fits, writes `models/value.json`, and plays the
held-out opponents. Sizes are run largest-first so a run killed part-way still
leaves the best model on disk.

The held-out opponents must be the same ones excluded from training, or the
whole point is lost -- both come from `--held-out` here, so they cannot drift
apart the way they could when the two commands were run separately.
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from junqi.arena import compare, seeds_for  # noqa: E402
from junqi.bot import BotWeights  # noqa: E402
from junqi.opponents import Pool, standard_pool  # noqa: E402
from junqi.value import auc, fit, reset_cache  # noqa: E402
from junqi.value_training import collect  # noqa: E402

DEFAULT_HELD_OUT = ("search-deep", "search-mid", "selective", "selective-strict")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[1600, 3200, 6400])
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--anchor", type=Path, default=Path("models/defaults.json"))
    parser.add_argument("--output", type=Path, default=Path("models/value.json"))
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument(
        "--stride",
        type=int,
        default=4,
        help="每隔几手取一个样本；6400 局全取会吃掉几个 GB 内存",
    )
    parser.add_argument(
        "--check-games",
        type=int,
        default=2400,
        help=(
            "每个规模用多少局做留出对手检验。这个默认值是算出来的："
            "2000 局训练出的模型在 n=1600 时实测 +0.0210 +/- 0.0153，"
            "要在 p<0.05 下检出 +0.021 需要 SE < 0.0128，即约 2400 局。"
            "用 400 局（SE 约 0.031）只会得到一串没有意义的 no effect"
        ),
    )
    parser.add_argument("--scale", type=float, default=15.0,
                        help="只检验一个 eval_value_scale。扫多个会抬高假阳性率")
    parser.add_argument("--held-out", nargs="*", default=list(DEFAULT_HELD_OUT))
    arguments = parser.parse_args()

    held_out = set(arguments.held_out)
    full = standard_pool(anchor=str(arguments.anchor))
    check_pool = Pool([spec for spec in full.specs if spec.name in held_out])
    missing = held_out - {spec.name for spec in check_pool.specs}
    if missing:
        raise SystemExit(f"not in the pool: {sorted(missing)}")

    shipped = BotWeights.load(arguments.model)
    candidate = replace(shipped, eval_value_scale=arguments.scale)
    check_seeds = seeds_for(arguments.check_games, check_pool, offset=555_000)

    print(f"held out of training AND used for the check: {sorted(held_out)}")
    print(f"eval_value_scale = {arguments.scale}, "
          f"{arguments.check_games} check games per size\n")

    results: list[tuple[int, int, float, float, float, float]] = []
    for size in sorted(arguments.sizes, reverse=True):
        started = time.perf_counter()
        print(f"=== {size} games ===", flush=True)
        games = [
            rows
            for rows in collect(
                size,
                arguments.model,
                arguments.anchor,
                arguments.workers,
                frozenset(held_out),
                arguments.stride,
            )
            if rows
        ]
        rows = [row for game in games for row in game]
        print(f"  {len(games)} decisive games, {len(rows)} positions", flush=True)

        order = list(range(len(games)))
        random.Random(11).shuffle(order)
        cut = int(len(order) * 0.7)
        train = [row for index in order[:cut] for row in games[index]]
        test = [row for index in order[cut:] for row in games[index]]
        model = fit([(v, y) for v, y, _, _ in train], epochs=arguments.epochs)

        learned = auc(
            [model.predict_features(v) for v, _, _, _ in test],
            [y for _, y, _, _ in test],
        )
        baseline = auc(
            [s for _, _, _, s in test], [y for _, y, _, _ in test]
        )
        model.save(arguments.output)
        reset_cache()  # this process must re-read the model it just wrote
        print(f"  held-out AUC: shipped {baseline:.3f} -> learned {learned:.3f}",
              flush=True)

        verdict = compare(
            candidate, shipped, check_pool, check_seeds, arguments.workers
        )
        print(
            f"  vs held-out opponents: {verdict.mean_difference:+.4f} "
            f"+/- {verdict.standard_error:.4f}, p = {verdict.p_value:.4f}, "
            f"n = {verdict.candidate.games}"
            f"   [{'PASS' if verdict.significant else 'no effect'}]"
        )
        print(f"  {time.perf_counter() - started:.0f}s\n", flush=True)
        results.append(
            (
                size,
                len(rows),
                baseline,
                learned,
                verdict.mean_difference,
                verdict.standard_error,
            )
        )

    print(f"{'games':>7}{'positions':>11}{'shipped':>10}{'learned':>10}"
          f"{'paired diff':>14}{'SE':>9}")
    print("-" * 61)
    for size, positions, base, learn, diff, error in sorted(results):
        print(f"{size:>7}{positions:>11}{base:>10.3f}{learn:>10.3f}"
              f"{diff:>+14.4f}{error:>9.4f}")
    print(
        "\nA curve still climbing means more games is the cheapest lever left.\n"
        "A flat one means this feature set is done, and the next question is\n"
        "whether the opponent pool -- not the model -- is what limits it."
    )


if __name__ == "__main__":
    main()
