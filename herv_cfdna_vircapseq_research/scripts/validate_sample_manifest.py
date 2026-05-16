from __future__ import annotations

import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_id",
    "cohort",
    "assay_type",
    "expected_group",
    "hiv_status",
    "htlv_status",
    "file_format",
}

VALID_ASSAYS = {"targeted", "WGS", "cfWGS", "unknown"}
VALID_STATUS = {"positive", "negative", "unknown", "not_applicable"}


def validate_manifest(path: Path) -> list[str]:
    errors: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            return [f"Missing required columns: {', '.join(sorted(missing))}"]

        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            sample_id = row.get("sample_id", "").strip()
            if not sample_id:
                errors.append(f"Line {line_number}: sample_id is required")
            elif sample_id in seen:
                errors.append(f"Line {line_number}: duplicate sample_id {sample_id}")
            seen.add(sample_id)

            assay_type = row.get("assay_type", "").strip()
            if assay_type not in VALID_ASSAYS:
                errors.append(f"Line {line_number}: assay_type must be one of {sorted(VALID_ASSAYS)}")

            for column in ("hiv_status", "htlv_status"):
                value = row.get(column, "").strip()
                if value not in VALID_STATUS:
                    errors.append(f"Line {line_number}: {column} must be one of {sorted(VALID_STATUS)}")

            file_format = row.get("file_format", "").strip()
            if file_format not in {"FASTQ", "BAM", "CRAM", "unknown"}:
                errors.append(f"Line {line_number}: file_format must be FASTQ, BAM, CRAM, or unknown")

    return errors


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python validate_sample_manifest.py <manifest.csv>")
        return 2
    errors = validate_manifest(Path(sys.argv[1]))
    if errors:
        print("Manifest validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print("Manifest validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
