from __future__ import annotations

import argparse
import csv
from pathlib import Path


PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def read_by_variant(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return {row["variant_id"]: row for row in reader}


def harmonize(exposure_csv: Path, outcome_csv: Path, output_csv: Path) -> None:
    exposure = read_by_variant(exposure_csv)
    outcome = read_by_variant(outcome_csv)
    fields = [
        "variant_id",
        "effect_allele",
        "other_allele",
        "beta_exposure",
        "se_exposure",
        "beta_outcome",
        "se_outcome",
        "outcome_trait",
        "harmonization_status",
        "harmonization_note",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for variant_id, exp in sorted(exposure.items()):
            out = outcome.get(variant_id)
            if out is None:
                writer.writerow(
                    {
                        "variant_id": variant_id,
                        "effect_allele": exp.get("effect_allele", ""),
                        "other_allele": exp.get("other_allele", ""),
                        "beta_exposure": exp.get("beta_exposure", ""),
                        "se_exposure": exp.get("se_exposure", ""),
                        "beta_outcome": "",
                        "se_outcome": "",
                        "outcome_trait": "",
                        "harmonization_status": "drop",
                        "harmonization_note": "variant absent from outcome summary statistics",
                    }
                )
                continue
            ea_exp = exp["effect_allele"].upper()
            oa_exp = exp["other_allele"].upper()
            ea_out = out["effect_allele"].upper()
            oa_out = out["other_allele"].upper()
            status = "keep"
            note = "alleles aligned"
            beta_outcome = float(out["beta_outcome"])
            if (ea_exp, oa_exp) in PALINDROMIC:
                status = "review"
                note = "palindromic allele pair; verify allele frequencies before final MR"
            if ea_exp == oa_out and oa_exp == ea_out:
                beta_outcome *= -1
                note = "outcome effect flipped to exposure effect allele"
            elif not (ea_exp == ea_out and oa_exp == oa_out):
                status = "drop"
                note = "allele mismatch"
            writer.writerow(
                {
                    "variant_id": variant_id,
                    "effect_allele": ea_exp,
                    "other_allele": oa_exp,
                    "beta_exposure": exp["beta_exposure"],
                    "se_exposure": exp["se_exposure"],
                    "beta_outcome": f"{beta_outcome:.10g}" if status != "drop" else "",
                    "se_outcome": out.get("se_outcome", "") if status != "drop" else "",
                    "outcome_trait": out.get("outcome_trait", ""),
                    "harmonization_status": status,
                    "harmonization_note": note,
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Harmonize exposure and outcome summary statistics for MR.")
    parser.add_argument("exposure_instruments_csv", type=Path)
    parser.add_argument("outcome_sumstats_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    harmonize(args.exposure_instruments_csv, args.outcome_sumstats_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
