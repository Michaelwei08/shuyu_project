from __future__ import annotations

import argparse
import csv
from pathlib import Path


MISSING_VALUES = {"", "NA", "N/A", "na", "n/a", "null", "None", "none"}


def check_missingness(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        counts = {field: 0 for field in fields}
        total = 0
        for row in reader:
            total += 1
            for field in fields:
                if row.get(field, "") in MISSING_VALUES:
                    counts[field] += 1

    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["column", "rows", "missing", "missing_fraction"])
        writer.writeheader()
        for field in fields:
            missing = counts[field]
            writer.writerow(
                {
                    "column": field,
                    "rows": total,
                    "missing": missing,
                    "missing_fraction": f"{missing / total:.6f}" if total else "",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize missingness for a CSV table.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    check_missingness(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
