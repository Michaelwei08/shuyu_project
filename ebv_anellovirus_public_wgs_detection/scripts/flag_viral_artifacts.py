from __future__ import annotations

import argparse
import csv
from pathlib import Path


def number(row: dict[str, str], column: str) -> float:
    value = row.get(column, "").strip()
    return float(value) if value else 0.0


def flag(row: dict[str, str], min_breadth: float, min_reads: float) -> str:
    total = number(row, "total_viral_reads")
    ebv = number(row, "ebv_reads")
    anello = number(row, "anellovirus_reads")
    breadth = number(row, "max_viral_breadth")
    if total == 0:
        return "no_viral_signal"
    if total < min_reads:
        return "low_read_count_review"
    if breadth and breadth < min_breadth:
        return "low_breadth_review"
    if ebv == 0 and anello == 0 and total > 0:
        return "non_target_viral_review"
    return "pass_first_review"


def run(input_csv: Path, output_csv: Path, min_breadth: float, min_reads: float) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "artifact_review_flag"])
        writer.writeheader()
        for row in reader:
            output = dict(row)
            output["artifact_review_flag"] = flag(row, min_breadth=min_breadth, min_reads=min_reads)
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Flag viral calls that need artifact review.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--min-breadth", type=float, default=0.01)
    parser.add_argument("--min-reads", type=float, default=2)
    args = parser.parse_args()
    run(args.input_csv, args.output_csv, args.min_breadth, args.min_reads)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
