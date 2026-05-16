from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def summarize(input_csv: Path, output_csv: Path) -> None:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "read_id", "virus", "host_chrom", "host_position", "junction_type"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        seen: set[tuple[str, str]] = set()
        for row in reader:
            sample_id = row["sample_id"]
            read_id = row["read_id"]
            if (sample_id, read_id) in seen:
                continue
            seen.add((sample_id, read_id))
            virus = row["virus"].lower()
            counts[sample_id]["total_junction_reads"] += 1
            if "hiv" in virus:
                counts[sample_id]["hiv_junction_reads"] += 1
            elif "htlv" in virus:
                counts[sample_id]["htlv_junction_reads"] += 1
            else:
                counts[sample_id]["other_viral_junction_reads"] += 1

    columns = ["sample_id", "total_junction_reads", "hiv_junction_reads", "htlv_junction_reads", "other_viral_junction_reads"]
    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        for sample_id in sorted(counts):
            row = {"sample_id": sample_id}
            row.update({column: counts[sample_id].get(column, 0) for column in columns if column != "sample_id"})
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize host-virus junction reads by sample.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    summarize(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
