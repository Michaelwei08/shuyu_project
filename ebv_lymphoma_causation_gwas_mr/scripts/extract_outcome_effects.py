from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_by_variant(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["variant_id"]: row for row in csv.DictReader(handle)}


def extract(instruments_csv: Path, outcome_csv: Path, output_csv: Path, trait: str) -> None:
    instruments = read_by_variant(instruments_csv)
    outcomes = read_by_variant(outcome_csv)
    fields = ["variant_id", "trait", "effect_allele", "other_allele", "beta_outcome", "se_outcome", "p_outcome", "status"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for variant_id, instrument in sorted(instruments.items()):
            outcome = outcomes.get(variant_id)
            if outcome is None:
                writer.writerow(
                    {
                        "variant_id": variant_id,
                        "trait": trait,
                        "effect_allele": instrument.get("effect_allele", ""),
                        "other_allele": instrument.get("other_allele", ""),
                        "beta_outcome": "",
                        "se_outcome": "",
                        "p_outcome": "",
                        "status": "missing_from_outcome_sumstats",
                    }
                )
            else:
                writer.writerow(
                    {
                        "variant_id": variant_id,
                        "trait": trait,
                        "effect_allele": outcome["effect_allele"],
                        "other_allele": outcome["other_allele"],
                        "beta_outcome": outcome["beta_outcome"],
                        "se_outcome": outcome["se_outcome"],
                        "p_outcome": outcome.get("p_outcome", ""),
                        "status": "found",
                    }
                )


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract outcome effect sizes for selected MR instruments.")
    parser.add_argument("instruments_csv", type=Path)
    parser.add_argument("outcome_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--trait", required=True)
    args = parser.parse_args()
    extract(args.instruments_csv, args.outcome_csv, args.output_csv, args.trait)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
