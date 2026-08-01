"""Run a model against the opponent pool at a chosen sample size.

    python -m junqi.benchmark --games 200                 # dev screening
    python -m junqi.benchmark --games 600                 # acceptance
    python -m junqi.benchmark --games 2000 --seeds 3      # release gate

`--seeds N` repeats the whole pool over N disjoint seed blocks, so a release
number is not one lucky block of openings.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from .arena import PoolReport, compare, evaluate, seeds_for
from .bot import BotWeights
from .opponents import discover_history, standard_pool


DEFAULT_ANCHOR = Path("models/defaults.json")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="对手池基准评测")
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--games", type=int, default=200)
    parser.add_argument("--seeds", type=int, default=1, help="独立种子块数量")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help=(
            "并行进程数。默认只用约三分之一核心并降低进程优先级，"
            "以免长时间评测把整台机器占死；无人使用时可手动调高"
        ),
    )
    parser.add_argument("--archive", type=Path, default=Path("models/history"))
    parser.add_argument(
        "--no-history",
        action="store_true",
        help="不要把历史模型加入对手池",
    )
    parser.add_argument(
        "--anchor",
        type=Path,
        default=None,
        help=(
            "把对手池中所有依赖权重的对手固定到这个模型文件。"
            f"默认 {DEFAULT_ANCHOR}（存在时）。不固定的话，对手会跟着被测模型一起变，"
            "对手池就不再是一把尺子"
        ),
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=None,
        help=(
            "与另一个模型做配对比较（相同开局、相同暗子采样）。"
            "训练期间对手池会随接受的模型一起变强，因此代际分数看不出净进步；"
            "要判断净进步，请对固定对手池用 --baseline 加 --no-history"
        ),
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    weights = (
        BotWeights.load(arguments.model)
        if arguments.model.exists()
        else BotWeights()
    )
    history = [] if arguments.no_history else discover_history(arguments.archive)
    anchor = arguments.anchor
    if anchor is None and DEFAULT_ANCHOR.exists():
        anchor = DEFAULT_ANCHOR
    pool = standard_pool(history=history, anchor=str(anchor) if anchor else None)
    print(f"anchor: {anchor if anchor else 'NONE (opponents track the model)'}")
    print(f"pool: {len(pool)} opponents -> {[s.name for s in pool.specs]}")

    per_block = max(1, arguments.games // max(1, arguments.seeds))
    blocks = [
        seeds_for(per_block, pool, offset=100_000 * (block + 1))
        for block in range(arguments.seeds)
    ]
    started = time.perf_counter()

    if arguments.baseline is not None:
        baseline = BotWeights.load(arguments.baseline)
        every_seed = [seed for block in blocks for seed in block]
        verdict = compare(weights, baseline, pool, every_seed, arguments.workers)
        print(f"\n=== {arguments.model.name} vs {arguments.baseline.name} ===")
        print("candidate:")
        print(verdict.candidate.format())
        print("\nbaseline:")
        print(verdict.incumbent.format())
        print()
        print(verdict.format())
        elapsed = time.perf_counter() - started
        played = verdict.candidate.games * 2
        print(f"\n{played} games in {elapsed:.0f}s")
        return

    all_results = []
    for block, seeds in enumerate(blocks):
        report = evaluate(weights, pool, seeds, arguments.workers)
        all_results.extend(report.results)
        print(
            f"\n--- seed block {block + 1}/{arguments.seeds} "
            f"({report.games} games) ---"
        )
        print(report.format())
    elapsed = time.perf_counter() - started

    if arguments.seeds > 1:
        combined = PoolReport(all_results)
        print(f"\n=== combined ({combined.games} games) ===")
        print(combined.format())
    print(
        f"\n{len(all_results)} games in {elapsed:.0f}s "
        f"({elapsed / max(1, len(all_results)):.2f}s/game wall)"
    )


if __name__ == "__main__":
    main()
