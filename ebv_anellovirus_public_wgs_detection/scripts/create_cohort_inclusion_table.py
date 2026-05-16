from __future__ import annotations

import argparse
import csv
from pathlib import Path


def status(row: dict[str, str]) -> str:
    if row.get("immediately_runnable", "").lower() == "yes":
        return "ready_now"
    if row.get("access_status", "").lower() == "controlled":
        return "needs_access"
    if row.get("accession_or_id", "") == "TBD":
        return "needs_accession_confirmation"
    return "needs_review"


def create(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        columns = ["disease", "cohort_name", "accession_or_id", "sample_count", "sequencing_type", "access_status", "priority", "inclusion_status", "notes"]
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        for row in reader:
            writer.writerow(
                {
                    "disease": row.get("disease", ""),
                    "cohort_name": row.get("cohort_name", ""),
                    "accession_or_id": row.get("accession_or_id", ""),
                    "sample_count": row.get("sample_count", ""),
                    "sequencing_type": row.get("sequencing_type", ""),
                    "access_status": row.get("access_status", ""),
                    "priority": row.get("priority", ""),
                    "inclusion_status": status(row),
                    "notes": row.get("notes", ""),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Create public WGS cohort inclusion table.")
    parser.add_argument("inventory_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    create(args.inventory_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
