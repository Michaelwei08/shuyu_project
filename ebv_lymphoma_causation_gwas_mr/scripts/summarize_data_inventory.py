from __future__ import annotations

import argparse
import csv
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1", "y"}


def status(row: dict[str, str]) -> str:
    has_manifest = truthy(row.get("has_sample_manifest", ""))
    has_viral = truthy(row.get("has_viral_detection", ""))
    has_genotype = truthy(row.get("has_genotype", "")) or truthy(row.get("has_imputed_genotype", ""))
    has_pcs = truthy(row.get("has_ancestry_pcs", ""))
    has_depth = truthy(row.get("has_depth", ""))
    has_outcome = truthy(row.get("has_outcome", ""))
    if has_manifest and has_viral and has_genotype and has_pcs and has_depth:
        return "gwas_ready"
    if has_manifest and has_viral and has_depth and has_outcome:
        return "association_only"
    if has_manifest and has_viral:
        return "phenotype_only"
    return "not_feasible_until_data_arrives"


def summarize(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=[*fields, "computed_status", "completion_note"])
        writer.writeheader()
        for row in reader:
            computed = status(row)
            output = dict(row)
            output["computed_status"] = computed
            output["completion_note"] = "inventory evaluated; real-data tasks remain blocked where data are absent"
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize P3 data inventory readiness.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    summarize(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
