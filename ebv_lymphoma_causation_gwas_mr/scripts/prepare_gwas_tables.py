from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_by_sample(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError(f"{path} must contain sample_id")
        return {row["sample_id"]: row for row in reader}


def write_tables(phenotype_csv: Path, covariate_csv: Path, phenotype_name: str, output_prefix: Path) -> None:
    phenotypes = load_by_sample(phenotype_csv)
    covariates = load_by_sample(covariate_csv)
    shared = sorted(set(phenotypes) & set(covariates))
    phenotype_path = output_prefix.with_suffix(".phenotype.tsv")
    covariate_path = output_prefix.with_suffix(".covariates.tsv")

    with phenotype_path.open("w", newline="", encoding="utf-8") as pheno_out:
        writer = csv.DictWriter(pheno_out, fieldnames=["FID", "IID", phenotype_name], delimiter="\t")
        writer.writeheader()
        for sample_id in shared:
            value = phenotypes[sample_id].get(phenotype_name, "")
            if value in {"", "NA", "filtered"}:
                continue
            writer.writerow({"FID": sample_id, "IID": sample_id, phenotype_name: value})

    covariate_fields = ["age", "sex", "ancestry_pc1", "ancestry_pc2", "ancestry_pc3", "ancestry_pc4", "ancestry_pc5", "cohort", "batch", "sample_type", "sequencing_depth"]
    with covariate_path.open("w", newline="", encoding="utf-8") as covar_out:
        writer = csv.DictWriter(covar_out, fieldnames=["FID", "IID", *covariate_fields], delimiter="\t")
        writer.writeheader()
        for sample_id in shared:
            row = covariates[sample_id]
            writer.writerow({"FID": sample_id, "IID": sample_id, **{field: row.get(field, "") for field in covariate_fields}})


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare GWAS phenotype and covariate tables from project CSV templates.")
    parser.add_argument("phenotype_csv", type=Path)
    parser.add_argument("covariate_csv", type=Path)
    parser.add_argument("phenotype_name")
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()
    write_tables(args.phenotype_csv, args.covariate_csv, args.phenotype_name, args.output_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
