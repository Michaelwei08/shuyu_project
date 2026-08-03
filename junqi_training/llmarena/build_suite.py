"""Generate the diagnostic battery.

    python -m llmarena.build_suite --games 400 --per-kind 150

Costs no API calls. Over-generates positions and keeps only the ones whose
label discriminates -- see :mod:`llmarena.probes` for why uniform sampling does
not work here.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from junqi.bot import BotWeights

from .probes import balanced, generate, write_jsonl


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--games", type=int, default=400)
    parser.add_argument("--per-kind", type=int, default=150)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--sample-every", type=int, default=6)
    parser.add_argument("--max-plies", type=int, default=120)
    parser.add_argument("--weights", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--out", type=Path, default=Path("data/probes.jsonl"))
    args = parser.parse_args()

    weights = (
        BotWeights.load(str(args.weights)) if args.weights.exists() else BotWeights()
    )
    seeds = range(args.seed_offset, args.seed_offset + args.games)
    harvested = list(
        generate(
            list(seeds),
            weights,
            max_plies=args.max_plies,
            sample_every=args.sample_every,
        )
    )
    available = Counter(probe.kind for probe in harvested)
    chosen = balanced(harvested, args.per_kind)
    written = Counter(probe.kind for probe in chosen)

    path = write_jsonl(chosen, args.out)
    print(f"harvested {len(harvested)} probes from {args.games} games")
    for kind in sorted(available):
        short = " (short)" if written[kind] < args.per_kind else ""
        print(f"  {kind:<16} available {available[kind]:>6}  kept {written[kind]:>5}{short}")
    print(f"wrote {len(chosen)} probes to {path}")


if __name__ == "__main__":
    main()
