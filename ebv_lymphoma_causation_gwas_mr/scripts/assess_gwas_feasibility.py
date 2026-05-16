from __future__ import annotations

import argparse
import csv
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1", "y"}


def readiness(row: dict[str, str]) -> tuple[str, str]:
    has_viral = truthy(row.get("has_viral_detection", ""))
    has_genotype = truthy(row.get("has_genotype", ""))
    has_imputed = truthy(row.get("has_imputed_genotype", ""))
    has_pcs = truthy(row.get("has_ancestry_pcs", ""))
    outcome = truthy(row.get("outcome_available", ""))
    has_depth = row.get("sequencing_depth", "").strip() not in {"", "NA", "unknown"}
    has_genotype_source = has_genotype or has_imputed
    if has_viral and has_genotype_source and has_pcs and has_depth:
        return "gwas_ready", "viral phenotype, genotype/imputed genotype, ancestry PCs, and depth are present"
    if has_viral and outcome and has_depth:
        return "association_only", "viral phenotype, outcome, and depth are present; genotype inputs are incomplete"
    if has_viral:
        return "phenotype_only", "viral phenotype is present; association covariates or outcomes are incomplete"
    return "blocked", "viral phenotype is missing"


def assess(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "computed_analysis_status", "computed_readiness_reason"])
        writer.writeheader()
        for row in reader:
            output = dict(row)
            computed_status, reason = readiness(row)
            output["computed_analysis_status"] = computed_status
            output["computed_readiness_reason"] = reason
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Assess GWAS and association feasibility from a data inventory table.")
    parser.add_argument("inventory_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    assess(args.inventory_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
