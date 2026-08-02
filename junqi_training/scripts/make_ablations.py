"""Write the baseline models the current changes are measured against.

`compare()` varies *weights*, not code, so anything wanting a paired p-value has
to be reachable by turning a coefficient on or off. These files are the "off"
ends of that switch:

    models/ab/pre-2026-08-01.json   the old blind-attack pricing

Regenerate after retraining, so the ablations track the shipped model:

    python scripts/make_ablations.py

## What the 2026-08-01 round actually measured (806 paired games each)

    blind pricing       +0.0218 +/- 0.0226, p = 0.17   inconclusive, kept
    defender supply     +0.0008 +/- 0.0097, p = 0.47   dead, shipped at 0
    screen cap 2 vs 3   -0.0985 +/- 0.0217             reverted to 3

Only the first is still a live question, and 806 games cannot resolve it: the
minimum detectable improvement at that sample size is +0.037, larger than the
effect. It needs `--games 2400 --seeds 3`.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from junqi.bot import BotWeights  # noqa: E402

#: The value `unknown_risk` carried in the model that played the ten replayed
#: games (fingerprint 9fc6aa1f, engine tag "2026-07-31 move-distance").
LEGACY_UNKNOWN_RISK = 0.0907965881


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--out", type=Path, default=Path("models/ab"))
    arguments = parser.parse_args()

    shipped = BotWeights.load(arguments.model)
    variants = {
        # Every other 2026-08-01 coefficient now ships at its old value, so this
        # single file is the whole "before" side of the remaining comparison.
        "pre-2026-08-01": replace(
            shipped, blind_battle=0.0, unknown_risk=LEGACY_UNKNOWN_RISK
        ),
    }
    # Separating the ordering fix from the rescaling. At 9.0 the branch spans
    # 14.6 points across the rank order where the old term spanned 0.73, and the
    # 2418-game result splits cleanly by opponent class: +0.5..+4.5 against all
    # eight greedy opponents, -1.6..-5.4 against all three search opponents.
    # That is the signature of a louder heuristic drowning out the rollout, not
    # of better odds. These keep the ordering and walk the volume back down.
    for scale in (2.0, 4.0, 6.0):
        variants[f"blind-{scale:g}"] = replace(shipped, blind_battle=scale)
    # The learned leaf evaluation is now ADOPTED at scale 15, so the shipped
    # model is the "on" end and this is the "off" end -- the file to re-baseline
    # against, and the one to reach for if a later change makes the fitted
    # evaluation look suspect.
    variants["value-off"] = replace(shipped, eval_value_scale=0.0)
    # Scales either side of the adopted one, for a later power-matched sweep.
    # Only 15 has been measured; 35 and 70 passed once at n=400, which today's
    # record says means nothing.
    for scale in (35.0, 70.0):
        variants[f"value-{scale:g}"] = replace(shipped, eval_value_scale=scale)
    arguments.out.mkdir(parents=True, exist_ok=True)
    for name, weights in variants.items():
        destination = arguments.out / f"{name}.json"
        weights.save(destination)
        print(f"wrote {destination}")

    stale = {f"{name}.json" for name in variants}
    for leftover in sorted(arguments.out.glob("*.json")):
        if leftover.name not in stale:
            leftover.unlink()
            print(f"removed stale ablation {leftover}")


if __name__ == "__main__":
    main()
