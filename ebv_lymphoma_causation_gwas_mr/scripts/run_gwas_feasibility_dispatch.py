from __future__ import annotations

import argparse
import csv
from pathlib import Path


PHENOTYPES = ["ebv_load", "anellovirus_load", "v1_pan_viral", "v3_high_ebv_low_anello"]


def read_qc_status(path: Path) -> tuple[bool, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return False, "no QC rows"
    failures = [row for row in rows if row.get("sample_qc_status") != "pass"]
    if failures:
        return False, f"{len(failures)} sample rows fail GWAS input QC"
    return True, "all sample rows pass GWAS input QC"


def dispatch(qc_csv: Path, output_csv: Path) -> None:
    ready, reason = read_qc_status(qc_csv)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["phenotype", "gwas_status", "command_or_reason"])
        writer.writeheader()
        for phenotype in PHENOTYPES:
            if ready:
                command = f"plink2 --pheno {phenotype}.phenotype.tsv --glm hide-covar --out {phenotype}_gwas"
                writer.writerow({"phenotype": phenotype, "gwas_status": "ready_to_run", "command_or_reason": command})
            else:
                writer.writerow({"phenotype": phenotype, "gwas_status": "not_feasible_until_data_arrives", "command_or_reason": reason})


def main() -> int:
    parser = argparse.ArgumentParser(description="Dispatch GWAS runs only when input QC passes.")
    parser.add_argument("qc_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    dispatch(args.qc_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
