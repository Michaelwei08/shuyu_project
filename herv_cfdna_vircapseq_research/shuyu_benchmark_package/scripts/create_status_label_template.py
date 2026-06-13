from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a fillable known-status template from a sample manifest.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "sample_id" not in rows[0]:
        raise ValueError("Manifest is empty or lacks sample_id")
    sample_ids = [row["sample_id"].strip() for row in rows]
    if any(not sample for sample in sample_ids):
        raise ValueError("Manifest contains an empty sample_id")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Manifest contains duplicate sample_id values")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "subject_id", "status", "notes"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sample_id": row["sample_id"],
                    "subject_id": row.get("subject_id", ""),
                    "status": "unknown",
                    "notes": "replace status with positive, negative, or keep unknown",
                }
            )
    print(f"Wrote {len(rows)} label rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
