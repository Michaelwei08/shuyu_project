from __future__ import annotations

import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "disease",
    "cohort_name",
    "accession_or_id",
    "data_source",
    "source_url",
    "sequencing_type",
    "file_availability",
    "access_status",
    "priority",
}


def validate(path: Path) -> list[str]:
    errors: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            return [f"Missing required columns: {', '.join(sorted(missing))}"]
        for line_number, row in enumerate(reader, start=2):
            if row.get("disease", "").strip() not in {"HL", "DLBCL", "PTLD", "other"}:
                errors.append(f"Line {line_number}: unsupported disease label")
            if not row.get("cohort_name", "").strip():
                errors.append(f"Line {line_number}: cohort_name is required")
            if row.get("priority", "").strip() not in {"high", "medium", "low", "exclude"}:
                errors.append(f"Line {line_number}: priority must be high, medium, low, or exclude")
    return errors


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python validate_cohort_inventory.py <cohort_inventory.csv>")
        return 2
    errors = validate(Path(sys.argv[1]))
    if errors:
        print("Cohort inventory validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print("Cohort inventory validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
