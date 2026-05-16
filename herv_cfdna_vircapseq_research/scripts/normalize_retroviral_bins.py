from __future__ import annotations

import argparse
import csv
from pathlib import Path


COUNT_COLUMNS = [
    "unique_hiv",
    "unique_htlv",
    "unique_herv_line1",
    "ambiguous_retroviral",
    "other_viral",
    "human",
    "unassigned",
]


def load_depths(path: Path) -> dict[str, float]:
    depths: dict[str, float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "sample_id" not in (reader.fieldnames or []) or "usable_read_pairs" not in (reader.fieldnames or []):
            raise ValueError("Depth CSV must contain sample_id and usable_read_pairs columns")
        for row in reader:
            depths[row["sample_id"]] = float(row["usable_read_pairs"])
    return depths


def normalize(input_csv: Path, output_csv: Path, depth_csv: Path | None) -> None:
    depths = load_depths(depth_csv) if depth_csv else {}
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        normalized_columns = [f"{column}_per_million" for column in COUNT_COLUMNS]
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "normalization_denominator", *normalized_columns])
        writer.writeheader()
        for row in reader:
            sample_id = row["sample_id"]
            denominator = depths.get(sample_id)
            if denominator is None:
                denominator = float(row.get("total_binned_reads", 0) or 0)
            output = dict(row)
            output["normalization_denominator"] = f"{denominator:.0f}"
            for column in COUNT_COLUMNS:
                value = float(row.get(column, 0) or 0)
                output[f"{column}_per_million"] = "" if denominator <= 0 else f"{value / denominator * 1_000_000:.6f}"
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize retroviral bin counts per million usable read pairs.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--depth-csv", type=Path)
    args = parser.parse_args()
    normalize(args.input_csv, args.output_csv, args.depth_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
