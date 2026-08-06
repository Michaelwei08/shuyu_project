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
    # The learned leaf evaluation ships at 0: real, but too small to clear the
    # aggregate bar (+0.0093 +/- 0.0078 at n=4836, p = 0.11). The shipped model
    # is therefore the "off" end and this is the "on" end -- flip to it the
    # moment the browser budget stops binding or the bar changes.
    variants["value-on"] = replace(shipped, eval_value_scale=15.0)
    # A piece that wins a headquarters probe can never move again, and nothing
    # priced that: `hq_strike` is 14.0 and rank-independent, so the branch is
    # indifferent between probing with an engineer and probing with the
    # commander -- while the rollout prefers the commander, because it survives.
    # Measured over 320 pool games: COMMANDER was entombed 41 times to
    # ENGINEER's 1, mean frozen value 8.8, 70.9 plies each.
    #
    # Three scales because the penalty has to be read against `hq_strike`. At
    # 1.0 a commander probe costs 11 of that 14 and an engineer probe 3, so the
    # ordering flips without forbidding the probe outright; 0.5 barely tilts it,
    # 2.0 makes any expensive probe worse than not probing.
    for scale in (0.5, 1.0, 2.0):
        variants[f"entomb-{scale:g}"] = replace(shipped, eval_hq_entomb=scale)
    # `reply_insight` has shipped at 0 since it was written and has never been
    # measured. Unlike every other candidate here it is a *loud* change:
    # divergence_check puts it at 85-94% of games with an SD of 0.55, against
    # entombment's 11% and 0.15. So it is genuinely tried, and a 2400-game run
    # can resolve it.
    #
    # 0 assumes the replier is blind to our ranks, and that assumption is
    # measurably wrong: tracking `OpponentKnowledge` from the opponent's seat
    # over 128 pool games, the mean belief about a surviving attacker of ours is
    # ~4.4 of ~9 movable ranks and 13.1% are pinned to a single rank. That puts
    # the empirically implied insight near 0.13-0.3, so 0.25 is the principled
    # value and the other two bracket it.
    for scale in (0.25, 0.5, 1.0):
        variants[f"insight-{scale:g}"] = replace(shipped, reply_insight=scale)
    # 2026-08-04, all three from a human's reading of a lost game. Sealing the
    # flag doors (a deployment change) was worth +0.0384; these are the
    # evaluation-side versions of the same intuitions.
    #
    # 1. Guarding the flag should outweigh a capture. `eval_hq_guard` already
    #    exists at 5.5 per guard, against `capture` at 2.8 plus battle terms --
    #    which is why a major general walked off a flag door to take a 3-point
    #    engineer. No new coefficient needed, just volume.
    #    **Adopted at 11.0 on 2026-08-04**, doubled from 5.5, on two independent
    #    2400-game runs: +0.0087 +/- 0.0064 and +0.0099 +/- 0.0080 (clustered),
    #    point estimates 0.0012 apart, combining to +0.0092 +/- 0.0050,
    #    one-sided p = 0.033. Neither run cleared 0.05 alone (both p ~ 0.08); the
    #    decision rule was fixed before the second one ran.
    #
    #    Named by absolute value now that the base has moved. `guard-off` is 5.5,
    #    the pre-adoption value and the switched-off end; 22.0 asks whether more
    #    is better still, and is genuinely untested -- the local probe that made
    #    it look bad reported the mean column, which is not trustworthy.
    variants["guard-off"] = replace(shipped, eval_hq_guard=5.5)
    variants["guard-22"] = replace(shipped, eval_hq_guard=22.0)
    # Belief from *behaviour* rather than only from battles. `OpponentKnowledge`
    # learned nothing from a quiet move; a railway corner identifies the piece as
    # an engineer with certainty, measured at 1.67 per game in 89% of games and
    # 267/267 correct. It matters more since the flag doors became fully mined:
    # an engineer is the only rank that clears a mine and survives.
    #
    # This is the first candidate aimed at the axis where the headroom actually
    # is -- perfect information is worth +18 points and none of it has been taken.
    variants["engineer-deduction"] = replace(shipped, use_engineer_deduction=1.0)
    # 2. Being one square from their flag should pay more than being six. The
    #    defensive side already has that lump (`eval_hq_breach` 26.0); the
    #    offensive side had nothing, so closing 2->1 paid what 6->5 pays.
    #    Divergence, n=256: 8.0 fires on 55.1% of games for +0.0385 +/- 0.0205,
    #    20.0 on 64.8% for +0.0652 +/- 0.0226. **Rising with scale**, which is
    #    the opposite of every null this project has recorded -- entombment and
    #    guard both got quieter when turned up. 40.0 is here to find where it
    #    turns over, because a term that only ever improves has not been bounded.
    #    **Adopted at 20.0 on 2026-08-04**, and it is the first replicated
    #    positive here: +0.0178 +/- 0.0083 (clustered 0.0092, p = 0.016) and, on
    #    openings it was not measured on, +0.0168 +/- 0.0088 (clustered 0.0093,
    #    p = 0.029). Point estimates 0.001 apart. `storm-off` is now the
    #    switched-off end; 40.0 remains to ask whether 20 is the right value,
    #    since the scale was picked off a local estimate that did not reproduce.
    variants["storm-off"] = replace(shipped, eval_hq_storm=0.0)
    for scale in (8.0, 40.0):
        variants[f"storm-{scale:g}"] = replace(shipped, eval_hq_storm=scale)
    # 3. A camp cannot be captured into, so a piece in one is untouchable -- and
    #    depth should matter. The flat `camp` pays 1.1 for every camp equally.
    for scale in (0.5, 1.5):
        variants[f"campdepth-{scale:g}"] = replace(shipped, camp_depth=scale)
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
