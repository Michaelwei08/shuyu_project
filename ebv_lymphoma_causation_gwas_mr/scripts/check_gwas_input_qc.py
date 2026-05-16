from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_PHENOTYPE_COLUMNS = {
    "sample_id",
    "ebv_load",
    "anellovirus_load",
    "pan_viral_burden",
    "phenotype_version",
    "depth_pass",
    "qc_status",
}

REQUIRED_COVARIATE_COLUMNS = {
    "sample_id",
    "age",
    "sex",
    "ancestry_pc1",
    "ancestry_pc2",
    "ancestry_pc3",
    "cohort",
    "batch",
    "sample_type",
    "sequencing_depth",
}

REQUIRED_GENOTYPE_COLUMNS = {
    "dataset",
    "genotype_path",
    "format",
    "genome_build",
    "has_sample_qc",
    "has_variant_qc",
    "has_ancestry_pcs",
}


def is_missing(value: str) -> bool:
    return value.strip() in {"", "NA", "N/A", "unknown", "Unknown"}


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1", "y"}


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def validate_columns(observed: list[str], required: set[str], label: str) -> None:
    missing = sorted(required.difference(observed))
    if missing:
        joined = ", ".join(missing)
        raise SystemExit(f"{label} is missing required columns: {joined}")


def genotype_manifest_status(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "missing_manifest_rows"
    if any(is_missing(row.get("genotype_path", "")) for row in rows):
        return "missing_genotype_path"
    if not all(truthy(row.get("has_sample_qc", "")) for row in rows):
        return "needs_sample_qc"
    if not all(truthy(row.get("has_variant_qc", "")) for row in rows):
        return "needs_variant_qc"
    if not all(truthy(row.get("has_ancestry_pcs", "")) for row in rows):
        return "needs_ancestry_pcs"
    return "genotype_manifest_ready"


def sample_qc_status(
    in_phenotype: bool,
    in_covariates: bool,
    phenotype_row: dict[str, str] | None,
    covariate_row: dict[str, str] | None,
    genotype_status: str,
) -> tuple[str, str]:
    reasons: list[str] = []
    if not in_phenotype:
        reasons.append("missing_phenotype")
    if not in_covariates:
        reasons.append("missing_covariates")
    if phenotype_row and not truthy(phenotype_row.get("depth_pass", "")):
        reasons.append("depth_fail_or_unset")
    if covariate_row:
        required_covariates = ["age", "sex", "ancestry_pc1", "ancestry_pc2", "ancestry_pc3", "sequencing_depth"]
        missing_covariates = [field for field in required_covariates if is_missing(covariate_row.get(field, ""))]
        if missing_covariates:
            reasons.append("missing_covariates:" + "|".join(missing_covariates))
    if genotype_status != "genotype_manifest_ready":
        reasons.append(genotype_status)
    if reasons:
        return "fail", ";".join(reasons)
    return "pass", "ready_for_gwas_input_formatting"


def check_inputs(phenotype_csv: Path, covariate_csv: Path, genotype_manifest_csv: Path, output_csv: Path) -> None:
    phenotype_fields, phenotype_rows = read_table(phenotype_csv)
    covariate_fields, covariate_rows = read_table(covariate_csv)
    genotype_fields, genotype_rows = read_table(genotype_manifest_csv)
    validate_columns(phenotype_fields, REQUIRED_PHENOTYPE_COLUMNS, "Phenotype table")
    validate_columns(covariate_fields, REQUIRED_COVARIATE_COLUMNS, "Covariate table")
    validate_columns(genotype_fields, REQUIRED_GENOTYPE_COLUMNS, "Genotype manifest")

    phenotype_by_sample = {row["sample_id"]: row for row in phenotype_rows}
    covariate_by_sample = {row["sample_id"]: row for row in covariate_rows}
    genotype_status = genotype_manifest_status(genotype_rows)
    sample_ids = sorted(set(phenotype_by_sample).union(covariate_by_sample))

    output_fields = [
        "sample_id",
        "in_phenotype",
        "in_covariates",
        "genotype_manifest_status",
        "sample_qc_status",
        "sample_qc_reason",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        for sample_id in sample_ids:
            phenotype_row = phenotype_by_sample.get(sample_id)
            covariate_row = covariate_by_sample.get(sample_id)
            status, reason = sample_qc_status(
                phenotype_row is not None,
                covariate_row is not None,
                phenotype_row,
                covariate_row,
                genotype_status,
            )
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "in_phenotype": str(phenotype_row is not None).lower(),
                    "in_covariates": str(covariate_row is not None).lower(),
                    "genotype_manifest_status": genotype_status,
                    "sample_qc_status": status,
                    "sample_qc_reason": reason,
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Check GWAS phenotype, covariate, and genotype-manifest readiness.")
    parser.add_argument("phenotype_csv", type=Path)
    parser.add_argument("covariate_csv", type=Path)
    parser.add_argument("genotype_manifest_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    check_inputs(args.phenotype_csv, args.covariate_csv, args.genotype_manifest_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
