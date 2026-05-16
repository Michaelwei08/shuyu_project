from __future__ import annotations

import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_id",
    "ebv_load",
    "anellovirus_load",
    "other_herpesvirus_load",
    "pan_viral_burden",
    "v1_pan_viral",
    "v3_high_ebv_low_anello",
    "phenotype_version",
}


def is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def validate(path: Path) -> list[str]:
    errors: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            return [f"Missing required columns: {', '.join(sorted(missing))}"]
        for line_number, row in enumerate(reader, start=2):
            if not row.get("sample_id", "").strip():
                errors.append(f"Line {line_number}: sample_id is required")
            for column in ("ebv_load", "anellovirus_load", "other_herpesvirus_load", "pan_viral_burden"):
                value = row.get(column, "").strip()
                if value and not is_number(value):
                    errors.append(f"Line {line_number}: {column} must be numeric")
            for column in ("v1_pan_viral", "v3_high_ebv_low_anello"):
                value = row.get(column, "").strip().lower()
                if value not in {"true", "false", "unknown", ""}:
                    errors.append(f"Line {line_number}: {column} must be true, false, unknown, or blank")
    return errors


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python validate_phenotype_table.py <phenotype.csv>")
        return 2
    errors = validate(Path(sys.argv[1]))
    if errors:
        print("Phenotype validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print("Phenotype validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
