"""Turn probe results into a finding.

    python -m llmarena.analyse --results data/probe-results-cli-full.jsonl

A headline accuracy is not a result -- ``belief 75%`` at n=12 has a 95%
interval running from 43% to 94%, which is consistent with almost anything.
This prints Wilson intervals so the width is visible, and then splits each
probe kind along the axis that would explain the errors:

* ``flag_candidates`` by how many headquarters are still standing. The
  one-candidate case is the decided one -- the deduction that actually hands
  the bot a target -- and it is much rarer than the two-candidate case, which
  is why :func:`llmarena.probes.balanced` stratifies on it.
* ``belief`` by the size of the deduced rank set. A singleton means the battle
  log pins the piece exactly; a set of eight means it barely constrains it. If
  accuracy falls as the set narrows, the model is failing at the *chain*, not
  at reading one battle.
* ``legal_moves`` by how many destinations the piece has. Small sets are road
  steps; large ones are rail slides, where the rules actually bite.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean

from .probes import PROBE_KINDS, read_jsonl


def wilson(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson score interval -- honest at small n, unlike normal approx."""
    if total == 0:
        return (0.0, 0.0)
    rate = successes / total
    denominator = 1 + z * z / total
    centre = (rate + z * z / (2 * total)) / denominator
    spread = (
        z
        * math.sqrt(rate * (1 - rate) / total + z * z / (4 * total * total))
        / denominator
    )
    return (max(0.0, centre - spread), min(1.0, centre + spread))


def _line(label: str, group: list[dict], width: int = 22) -> str:
    exact = sum(1 for r in group if r["exact"])
    low, high = wilson(exact, len(group))
    return (
        f"  {label:<{width}} {len(group):>5} "
        f"{exact / len(group):>7.1%}  [{low:.0%}, {high:.0%}]"
        f"  jaccard {fmean(r['jaccard'] for r in group):.3f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", type=Path, default=Path("data/probe-results-cli-full.jsonl")
    )
    parser.add_argument("--probes", type=Path, default=Path("data/probes.jsonl"))
    parser.add_argument("--show-failures", type=int, default=6)
    args = parser.parse_args()

    with args.results.open(encoding="utf-8") as stream:
        results = [json.loads(line) for line in stream if line.strip()]
    results = [r for r in results if "error" not in r]
    meta = {probe.probe_id: probe for probe in read_jsonl(args.probes)}

    print(f"{len(results)} graded probes\n")
    print(f"  {'group':<22} {'n':>5} {'exact':>7}  {'95% CI':<12} {'jaccard':>7}")
    print("  " + "-" * 62)
    for kind in PROBE_KINDS:
        group = [r for r in results if r["kind"] == kind]
        if group:
            print(_line(kind, group))
    print("  " + "-" * 62)
    print(_line("ALL", results))

    # --- flag_candidates: is the decided case harder than the ambiguous one?
    flags = [r for r in results if r["kind"] == "flag_candidates"]
    if flags:
        print("\nflag_candidates by headquarters still standing:")
        buckets: dict[int, list[dict]] = defaultdict(list)
        for result in flags:
            probe = meta.get(result["probe_id"])
            buckets[probe.meta.get("count", 0) if probe else 0].append(result)
        for count in sorted(buckets):
            label = "1 left (decided)" if count == 1 else f"{count} left"
            print(_line(label, buckets[count]))

    # --- belief: does accuracy track how tight the deduction is?
    beliefs = [r for r in results if r["kind"] == "belief"]
    if beliefs:
        print("\nbelief by size of the deduced rank set (1 = pinned exactly):")
        buckets = defaultdict(list)
        for result in beliefs:
            probe = meta.get(result["probe_id"])
            size = probe.meta.get("size", len(result["label"])) if probe else 0
            buckets[min(size, 6)].append(result)
        for size in sorted(buckets):
            label = f"{size} ranks" if size < 6 else "6+ ranks"
            print(_line(label, buckets[size]))

        scored = [r for r in beliefs if "correct" in r]
        if scored:
            correct = sum(1 for r in scored if r["correct"])
            low, high = wilson(correct, len(scored))
            tighter = sum(
                1 for r in scored if r["correct"] and r.get("tighter_by", 0) > 0
            )
            unsound = sum(1 for r in scored if not r["sound"])
            wide = sum(1 for r in scored if r["over_wide"])
            print(
                "\nbelief scored against the TRUE rank, not the engine's set"
                "\n(the engine does per-square battle constraints only -- no rank"
                "\n counting, no bomb exhaustion -- so its set is a loose upper"
                "\n bound and a better answer can be a strict subset):"
            )
            print(
                f"  correct                          {correct:>4}/{len(scored)}"
                f"  {correct / len(scored):.1%}  [{low:.0%}, {high:.0%}]"
            )
            print(f"  ...of which BEAT the engine      {tighter:>4}"
                  f"   (eliminated ranks the label still allowed)")
            print(f"  excluded the true rank           {unsound:>4}   <- real error")
            print(f"  kept a proven-impossible rank    {wide:>4}   <- real error")
        else:
            print(
                "\nbelief: no true_kind in the probe metadata -- regenerate the "
                "suite to score against ground truth instead of the engine's "
                "upper bound (see llmarena.probes.score)"
            )

    # --- tightest: how far below the engine's set, without losing the truth
    tight = [r for r in results if r["kind"] == "tightest" and "correct" in r]
    if tight:
        correct = [r for r in tight if r["correct"]]
        gains = [r["tighter_by"] for r in correct]
        print("\ntightest (ask for the smallest provable set):")
        low, high = wilson(len(correct), len(tight))
        print(
            f"  correct                          {len(correct):>4}/{len(tight)}"
            f"  {len(correct) / len(tight):.1%}  [{low:.0%}, {high:.0%}]"
        )
        print(f"  excluded the true rank           "
              f"{sum(1 for r in tight if not r['sound']):>4}   <- real error")
        print(f"  kept a proven-impossible rank    "
              f"{sum(1 for r in tight if r['over_wide']):>4}   <- real error")
        if gains:
            beat = sum(1 for g in gains if g > 0)
            print(
                f"  ranks eliminated beyond the engine: mean {fmean(gains):.2f}, "
                f"max {max(gains)}, better on {beat}/{len(correct)}"
            )
            print(
                "  (the engine's own score here is 0 by construction -- this is"
                "\n   the headroom `eliminate_dead_ranks` is failing to capture)"
            )

    # --- legal_moves: road steps versus rail slides
    legals = [r for r in results if r["kind"] == "legal_moves"]
    if legals:
        print("\nlegal_moves by number of destinations:")
        buckets = defaultdict(list)
        for result in legals:
            size = len(result["label"])
            buckets["5-7 (short)" if size <= 7 else "8+ (rail slide)"].append(result)
        for label in sorted(buckets):
            print(_line(label, buckets[label]))

    worst = sorted(results, key=lambda r: r["jaccard"])[: args.show_failures]
    if worst and worst[0]["jaccard"] < 1.0:
        print(f"\nworst {args.show_failures} answers:")
        for result in worst:
            if result["jaccard"] >= 1.0:
                continue
            print(f"  {result['probe_id']}  jaccard {result['jaccard']:.2f}")
            print(f"    label:     {' '.join(result['label'])}")
            print(f"    predicted: {' '.join(result['predicted']) or '(nothing)'}")


if __name__ == "__main__":
    main()
