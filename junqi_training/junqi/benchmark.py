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

from .arena import DEFAULT_SEARCH, PoolReport, compare, evaluate, seeds_for
from .bot import BotWeights
from .opponents import discover_exploiters, discover_history, standard_pool


DEFAULT_ANCHOR = Path("models/defaults.json")
DEFAULT_MODEL = Path("models/bot_weights.json")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="对手池基准评测")
    # A sentinel rather than the path itself, so an explicitly passed --model
    # that does not exist can be told apart from the default one being absent on
    # a fresh checkout. See the guard in `main`.
    parser.add_argument("--model", type=Path, default=None)
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
    parser.add_argument(
        "--screen-cap",
        type=int,
        default=None,
        help=(
            "被测方军旗大本营相邻格最多埋几颗雷（默认 3，即全封）。"
            "传 2 可复现 2026-08-01 那次被实测否决的松雷布局"
        ),
    )
    parser.add_argument(
        "--baseline-screen-cap",
        type=int,
        default=None,
        help=(
            "对照模型的同一参数。只改被测方自己的布阵，"
            "对手军队不变，所以配对比较依然成立"
        ),
    )
    parser.add_argument(
        "--deployment",
        type=str,
        default=None,
        help=(
            "被测方使用的布阵族（见 junqi.deployment.FAMILIES）。"
            "只改被测方自己的阵型，对手军队不变，配对依然成立。"
            "--screen-cap 是它的旧接口，等价于一个只改雷数的临时族"
        ),
    )
    parser.add_argument(
        "--baseline-deployment",
        type=str,
        default=None,
        help="对照模型的布阵族",
    )
    parser.add_argument(
        "--search",
        type=int,
        nargs=3,
        metavar=("SAMPLES", "BEAM", "REPLIES"),
        default=None,
        help=(
            f"被测方的搜索预算，默认 {' '.join(map(str, DEFAULT_SEARCH))}。"
            "所有权重都是在这个较浅的设置下调出来的，而网页“深思”档实际是 "
            "28 18 5 外加一层 Python 侧没有的延伸搜索。"
            "用 --search 28 18 5 可检验某个结论在实际部署深度下是否还成立"
            "（同时作用于候选与对照，保持配对）"
        ),
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    # A model path the user typed and that does not exist must be fatal. It used
    # to fall back to `BotWeights()`, which differs from the shipped model in 14
    # of 36 fields -- so a mistyped or unsynced `--model models/ab/foo.json`
    # quietly measured *default weights* against the baseline and printed a
    # perfectly ordinary-looking result. `models/` never travels with
    # `sync_remote.py push`, so an ablation that has not been regenerated on the
    # box is exactly how this happens.
    if arguments.model is not None and not arguments.model.exists():
        raise SystemExit(
            f"--model {arguments.model} does not exist.\n"
            "If this is an ablation, `models/` is not synced by "
            "scripts/sync_remote.py -- run `python scripts/make_ablations.py` "
            "first. Refusing to fall back to default weights, which would "
            "measure something you did not ask for."
        )
    model_path = arguments.model or DEFAULT_MODEL
    weights = (
        BotWeights.load(model_path) if model_path.exists() else BotWeights()
    )
    history = [] if arguments.no_history else discover_history(arguments.archive)
    anchor = arguments.anchor
    if anchor is None and DEFAULT_ANCHOR.exists():
        anchor = DEFAULT_ANCHOR
    exploiters = discover_exploiters(Path("models/exploiters"))
    pool = standard_pool(
        history=history,
        anchor=str(anchor) if anchor else None,
        exploiters=exploiters,
    )
    print(f"anchor: {anchor if anchor else 'NONE (opponents track the model)'}")
    print(f"pool: {len(pool)} opponents -> {[s.name for s in pool.specs]}")

    # A named family wins over the older bare screen cap; passing neither keeps
    # the shipped generator.
    deployment = arguments.deployment or arguments.screen_cap
    baseline_deployment = arguments.baseline_deployment or arguments.baseline_screen_cap
    if deployment is not None or baseline_deployment is not None:
        print(
            f"deployment: subject={deployment or 'standard'} "
            f"baseline={baseline_deployment or 'standard'}"
        )

    search = tuple(arguments.search) if arguments.search else None
    if search:
        print(f"search budget: samples={search[0]} beam={search[1]} replies={search[2]}"
              f"  (default {DEFAULT_SEARCH})")

    per_block = max(1, arguments.games // max(1, arguments.seeds))
    blocks = [
        seeds_for(per_block, pool, offset=100_000 * (block + 1))
        for block in range(arguments.seeds)
    ]
    started = time.perf_counter()

    if arguments.baseline is not None:
        baseline = BotWeights.load(arguments.baseline)
        every_seed = [seed for block in blocks for seed in block]
        verdict = compare(
            weights,
            baseline,
            pool,
            every_seed,
            arguments.workers,
            candidate_deployment=deployment,
            incumbent_deployment=baseline_deployment,
            search=search,
        )
        print(f"\n=== {model_path.name} vs {arguments.baseline.name} ===")
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
        report = evaluate(
            weights, pool, seeds, arguments.workers, deployment, search
        )
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
