"""Re-emit `models/defaults.json` as a complete, explicit weight vector.

The anchor pins every weight-driven opponent in the judging pool, so it is
supposed to be a fixed yardstick. It was not one.

`BotWeights.load` keeps only the keys it recognises and lets the dataclass
supply the rest -- deliberately, so an archived model written before a
coefficient existed still loads. Applied to the anchor that behaviour is a
silent trap: every coefficient added to `BotWeights` since the file was written
resolves to *today's* default, so the yardstick moves on every commit that adds
a field, with the file on disk untouched.

Measured on 2026-08-01: the anchor was missing seven keys and was running
`unknown_risk = 0.12` **and** `blind_battle = 9.0` at the same time -- both
blind-attack penalties at once, a configuration never shipped and never
intended (the `unknown_risk` docstring says the old behaviour is `unknown_risk
= 0.0908` AND `blind_battle = 0`). Ten of thirteen pool opponents were built
from it, so the pool systematically under-attacked unknown squares, which is
precisely the axis a human beat this bot on: 76% of self-initiated battles won
against the bot's 46%.

The fix is to write every field explicitly. The only *behavioural* change is
`blind_battle 9.0 -> 0.0`, which removes the double penalty and restores the
single coherent pricing the anchor had before that field existed. Everything
else is pinned at what the anchor was already doing, so the yardstick moves as
little as possible while becoming a yardstick again.

    python scripts/rebuild_anchor.py            # show the diff
    python scripts/rebuild_anchor.py --write

`test_the_anchor_pins_every_coefficient` fails if a new field is ever added
without rerunning this.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, fields, replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from junqi.bot import BotWeights  # noqa: E402

DEFAULT_ANCHOR = Path("models/defaults.json")

#: What each field should be when the anchor never carried it. These reproduce
#: the anchor's behaviour from *before* the field existed, not today's default.
BACKFILL = {
    # The branch did not exist, so no prior-odds term was applied. Leaving this
    # at the dataclass default is what put both blind penalties on at once.
    "blind_battle": 0.0,
    "eval_hq_supply": 0.0,
    "reply_insight": 0.0,
    "eval_value_scale": 0.0,
    # Python scored `base + rollout` before this became a coefficient.
    "search_base_weight": 1.0,
}


def rebuild(path: Path) -> tuple[BotWeights, dict[str, tuple[float, float]]]:
    """Return the repaired anchor and the fields whose behaviour changes."""
    raw = json.loads(path.read_text(encoding="utf-8"))
    current = BotWeights.load(path)
    overrides = {
        name: value for name, value in BACKFILL.items() if name not in raw
    }
    repaired = replace(current, **overrides)

    before, after = asdict(current), asdict(repaired)
    changed = {
        name: (before[name], after[name])
        for name in before
        if before[name] != after[name]
    }
    return repaired, changed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor", type=Path, default=DEFAULT_ANCHOR)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()

    raw = json.loads(arguments.anchor.read_text(encoding="utf-8"))
    absent = [
        descriptor.name
        for descriptor in fields(BotWeights)
        if descriptor.name not in raw
    ]
    repaired, changed = rebuild(arguments.anchor)

    print(f"{arguments.anchor}: {len(raw)} keys on disk, "
          f"{len(fields(BotWeights))} fields in BotWeights")
    if absent:
        print("\nabsent from the file, so silently taking today's default:")
        for name in absent:
            print(f"  {name:<22} -> {getattr(repaired, name)}")
    if changed:
        print("\nBEHAVIOUR CHANGES once pinned explicitly:")
        for name, (before, after) in changed.items():
            print(f"  {name:<22} {before} -> {after}")
    else:
        print("\nno behaviour change; the anchor was already complete")

    if not arguments.write:
        print("\n(dry run -- pass --write to apply)")
        return
    repaired.save(arguments.anchor)
    print(f"\nwrote {arguments.anchor} with all "
          f"{len(fields(BotWeights))} fields explicit")
    print("Numbers measured against the old anchor are not comparable: "
          "re-baseline before quoting them.")


if __name__ == "__main__":
    main()
