"""Play an LLM against the search ladder and report its search-equivalent.

    python -m llmarena.run_matches --model claude-opus-5 --scaffold legal --seeds 3

Prints the plan and exits unless ``--confirm`` is passed, because this is the
script that spends money.

The headline number is **search-equivalent**: the rollout budget at which the
subject's win rate crosses 50%. "GPT-X with legal-move filtering plays like 40
rollouts" is interpretable and comparable across models in a way an Elo figure
pinned to an arbitrary anchor is not. Elo can still be fitted over the same
games afterwards; it is the display layer, not the measurement.

Two guarantees carried over from `junqi.arena`, and one that is lost:

* **Kept -- colour-swapped pairs.** Every seed is played twice with the sides
  exchanged, so first-move and layout advantage cancel.
* **Kept -- one opening per seed.** `make_opening` is arithmetic on the seed, so
  a given seed is the same position for every rung.
* **Lost -- common random numbers.** `arena.compare` is exact because two
  candidates meet byte-identical states; `temperature` no longer exists on this
  model family, so an LLM cannot be pinned to one sample. A warm prompt cache
  restores replay-level reproducibility but not first-run determinism. Standard
  errors here are therefore larger than the ones in `CLAUDE.md`, and matches
  sharing a seed are correlated -- that round measured a design effect of 1.26,
  so multiply any SE by ~1.12 before believing a marginal p-value.

Threads, not the arena's process pool: these matches are network-bound, and
`run_jobs` aborts a whole run above a 2% failure rate, which one rate-limit
storm would trip.
"""

from __future__ import annotations

import argparse
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import fmean

from junqi.arena import Job, play_match
from junqi.bot import BotWeights
from junqi.opponents import AgentSpec
from junqi.types import Owner

from .anthropic_completer import USAGE
from .cost import MODELS, estimate, estimate_tokens
from .view import RULES

#: `(name, samples, beam_width, reply_width)`. Strength is summarised as the
#: product -- the number of rollouts a move decision costs -- which is a
#: definition, not a measured equivalence, but it is monotone and spans
#: two orders of magnitude across the ladder.
LADDER: tuple[tuple[str, int, int, int], ...] = (
    ("search-shallow", 1, 4, 1),
    ("search-mid", 3, 8, 3),
    ("search-deep", 6, 10, 4),
)

MEAN_PLIES = 85


def rollouts(samples: int, beam: int, reply: int) -> int:
    return samples * beam * reply


def search_equivalent(points: list[tuple[int, float]]) -> str:
    """Interpolate the rollout budget where the win rate crosses 50%."""
    ordered = sorted(points)
    if all(rate >= 0.5 for _, rate in ordered):
        return f"> {ordered[-1][0]} rollouts (beat every rung)"
    if all(rate < 0.5 for _, rate in ordered):
        return f"< {ordered[0][0]} rollouts (lost to every rung)"
    for (low_x, low_y), (high_x, high_y) in zip(ordered, ordered[1:]):
        if low_y >= 0.5 > high_y:
            span = low_y - high_y
            share = (low_y - 0.5) / span if span else 0.0
            log_x = math.log10(low_x) + share * (
                math.log10(high_x) - math.log10(low_x)
            )
            return f"~{10 ** log_x:.0f} rollouts"
    return "non-monotone across the ladder -- report the per-rung table instead"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="claude-opus-5")
    parser.add_argument("--scaffold", default="legal", choices=("raw", "legal", "derived"))
    parser.add_argument("--effort", default="low")
    parser.add_argument("--thinking", default="adaptive", choices=("adaptive", "disabled"))
    parser.add_argument(
        "--seeds", type=int, default=3, help="seeds per rung; each is played twice"
    )
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument(
        "--rungs",
        nargs="*",
        default=[name for name, *_ in LADDER],
        help="subset of the ladder, or 'random' for the floor check",
    )
    parser.add_argument("--anchor", type=Path, default=Path("models/defaults.json"))
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    parser.add_argument("--confirm", action="store_true", help="actually spend money")
    args = parser.parse_args()

    subject = AgentSpec(
        f"llm-{args.model}",
        "external",
        builder="llmarena.agent:build_agent",
        options=(
            ("completer", "anthropic"),
            ("model", args.model),
            ("scaffold", args.scaffold),
            ("effort", args.effort),
            ("thinking", args.thinking),
            ("cache", str(args.cache)),
        ),
    )

    anchor = str(args.anchor) if args.anchor.exists() else None
    rungs: list[tuple[str, AgentSpec, int]] = []
    for name in args.rungs:
        if name == "random":
            rungs.append(("random", AgentSpec("random", "random"), 1))
            continue
        match = next((row for row in LADDER if row[0] == name), None)
        if match is None:
            raise SystemExit(f"unknown rung {name!r}; choose from {[r[0] for r in LADDER]} or 'random'")
        label, samples, beam, reply = match
        rungs.append(
            (
                label,
                AgentSpec(
                    label,
                    "search",
                    weights_path=anchor,
                    samples=samples,
                    beam_width=beam,
                    reply_width=reply,
                ),
                rollouts(samples, beam, reply),
            )
        )

    seeds = list(range(args.seed_offset, args.seed_offset + args.seeds))
    weights = BotWeights.load(anchor) if anchor else BotWeights()
    jobs = [
        (label, budget, Job(weights, spec, seed, side, subject_spec=subject))
        for label, spec, budget in rungs
        for seed in seeds
        for side in (int(Owner.HUMAN), int(Owner.BOT))
    ]

    calls = round(len(jobs) * MEAN_PLIES * 0.5)
    print(f"plan: {len(jobs)} games ({len(rungs)} rungs x {args.seeds} seeds x 2 colours)")
    print(f"      ~{calls:,} model calls at ~{MEAN_PLIES} plies/game")
    if args.model in MODELS:
        guess = estimate(
            args.model, calls, 1500, 600, cacheable_tokens=estimate_tokens(RULES)
        )
        print(f"      rough cost ~${guess.dollars:.2f} (character-based estimate)")
    print(f"      1 game in 12 hits the 300-ply cap -- budget headroom\n")
    if not args.confirm:
        print("dry run. re-run with --confirm to play.")
        return

    results: dict[str, list] = {label: [] for label, *_ in rungs}
    budgets = {label: budget for label, _, budget in rungs}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(play_match, job): label for label, _, job in jobs}
        done = 0
        for future in as_completed(futures):
            label = futures[future]
            done += 1
            try:
                results[label].append(future.result())
            except Exception as error:  # noqa: BLE001 - one bad game is not fatal
                print(f"  game failed ({label}): {error!r}", flush=True)
            if done % max(1, len(jobs) // 10) == 0:
                print(f"  {done}/{len(jobs)} games", flush=True)

    print()
    print(f"{'rung':<18} {'rollouts':>9} {'games':>6} {'LLM win%':>9} {'score':>8}")
    print("-" * 55)
    points: list[tuple[int, float]] = []
    for label, *_ in rungs:
        group = results[label]
        if not group:
            continue
        # MatchResult is from the *subject's* view, and the subject is the LLM.
        win = fmean(r.result for r in group)
        print(
            f"{label:<18} {budgets[label]:>9} {len(group):>6} "
            f"{win:>8.1%} {fmean(r.score for r in group):>8.3f}"
        )
        if label != "random":
            points.append((budgets[label], win))

    if points:
        print("-" * 55)
        print(f"search-equivalent: {search_equivalent(points)}")
    print(f"\nusage: {USAGE.format()}")


if __name__ == "__main__":
    main()
