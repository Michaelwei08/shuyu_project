from __future__ import annotations

import argparse
import csv
from pathlib import Path


SUBTYPES = {
    "V1_pan_viral": "v1_pan_viral",
    "V3_high_ebv_low_anello": "v3_high_ebv_low_anello",
}


def number(value: str) -> float:
    return float(value)


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1", "y"}


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def summarize(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    fields = [
        "scope",
        "subtype",
        "n_samples",
        "mean_ebv_load",
        "mean_anellovirus_load",
        "adverse_outcome_fraction",
        "immunosuppression_fraction",
        "interpretation_hint",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        scopes = ["all_diseases", "HL", "PTLD", "DLBCL"]
        for scope in scopes:
            scoped_rows = [row for row in rows if row["disease_label"] != "control"] if scope == "all_diseases" else [row for row in rows if row["disease_label"] == scope]
            for subtype, column in SUBTYPES.items():
                selected = [row for row in scoped_rows if truthy(row.get(column, ""))]
                if not selected:
                    writer.writerow(
                        {
                            "scope": scope,
                            "subtype": subtype,
                            "n_samples": 0,
                            "mean_ebv_load": "",
                            "mean_anellovirus_load": "",
                            "adverse_outcome_fraction": "",
                            "immunosuppression_fraction": "",
                            "interpretation_hint": "no subtype-positive samples in example input",
                        }
                    )
                    continue
                ebv = [number(row["ebv_load"]) for row in selected]
                anello = [number(row["anellovirus_load"]) for row in selected]
                adverse = [1.0 if row.get("clinical_outcome") == "adverse" else 0.0 for row in selected]
                immune = [0.0 if row.get("immunosuppression_status") in {"none", ""} else 1.0 for row in selected]
                if subtype.startswith("V3"):
                    hint = "candidate disease-associated EBV pattern if replicated in real data"
                else:
                    hint = "candidate immune-state marker if immunosuppression fraction remains high"
                writer.writerow(
                    {
                        "scope": scope,
                        "subtype": subtype,
                        "n_samples": len(selected),
                        "mean_ebv_load": f"{mean(ebv):.6f}",
                        "mean_anellovirus_load": f"{mean(anello):.6f}",
                        "adverse_outcome_fraction": f"{mean(adverse):.6f}",
                        "immunosuppression_fraction": f"{mean(immune):.6f}",
                        "interpretation_hint": hint,
                    }
                )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize V1/V3 subtype interpretation features by disease.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    summarize(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
