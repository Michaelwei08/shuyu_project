from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def summarize(input_csv: Path, output_csv: Path, min_depth: int) -> None:
    covered = defaultdict(int)
    total = defaultdict(int)
    depth_sum = defaultdict(float)
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "virus", "position", "depth"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            key = (row["sample_id"], row["virus"])
            depth = float(row["depth"])
            total[key] += 1
            depth_sum[key] += depth
            if depth >= min_depth:
                covered[key] += 1

    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["sample_id", "virus", "positions_observed", "positions_covered", "breadth", "mean_depth"])
        writer.writeheader()
        for sample_id, virus in sorted(total):
            observed = total[(sample_id, virus)]
            covered_count = covered[(sample_id, virus)]
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "virus": virus,
                    "positions_observed": observed,
                    "positions_covered": covered_count,
                    "breadth": f"{covered_count / observed:.6f}" if observed else "",
                    "mean_depth": f"{depth_sum[(sample_id, virus)] / observed:.6f}" if observed else "",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize viral breadth of coverage from per-position depth.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--min-depth", type=int, default=1)
    args = parser.parse_args()
    summarize(args.input_csv, args.output_csv, args.min_depth)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
