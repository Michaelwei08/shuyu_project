from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_calls(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def compare(calls_csv: Path, output_csv: Path) -> None:
    rows = load_calls(calls_csv)
    by_subject: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_subject[row.get("subject_id", "")].append(row)

    fields = ["comparison_type", "subject_or_group", "sample_ids", "hiv_calls", "htlv_calls", "note"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for subject_id, subject_rows in sorted(by_subject.items()):
            assays = {row.get("assay_type") for row in subject_rows}
            if {"targeted", "WGS"}.issubset(assays):
                writer.writerow(
                    {
                        "comparison_type": "targeted_vs_wgs_overlap",
                        "subject_or_group": subject_id,
                        "sample_ids": "|".join(row["sample_id"] for row in subject_rows),
                        "hiv_calls": "|".join(f"{row['assay_type']}:{row['hiv_call']}" for row in subject_rows),
                        "htlv_calls": "|".join(f"{row['assay_type']}:{row['htlv_call']}" for row in subject_rows),
                        "note": "subject has both assay types in dry-run data",
                    }
                )
        for group in ["hiv_cohort", "local_60_wgs"]:
            group_rows = [row for row in rows if row.get("cohort") == group and row.get("assay_type") == "WGS"]
            if not group_rows:
                continue
            writer.writerow(
                {
                    "comparison_type": "wgs_group_summary",
                    "subject_or_group": group,
                    "sample_ids": "|".join(row["sample_id"] for row in group_rows),
                    "hiv_calls": "|".join(row["hiv_call"] for row in group_rows),
                    "htlv_calls": "|".join(row["htlv_call"] for row in group_rows),
                    "note": "WGS dry-run comparison across HIV cohort and HL/local WGS samples",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare targeted/WGS overlap and WGS groups from frozen calls.")
    parser.add_argument("calls_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    compare(args.calls_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
