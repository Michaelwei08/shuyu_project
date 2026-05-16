from __future__ import annotations

import argparse
import csv
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1", "y"}


def summarize(input_csv: Path, output_csv: Path, p_threshold: float) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"variant_id", "chrom", "position", "p"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        rows = [row for row in reader if float(row["p"]) <= p_threshold]

    fields = ["variant_id", "chrom", "position", "p", "in_mhc_region", "hit_class", "note"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in sorted(rows, key=lambda item: float(item["p"])):
            in_mhc = truthy(row.get("in_mhc_region", "false"))
            writer.writerow(
                {
                    "variant_id": row["variant_id"],
                    "chrom": row["chrom"],
                    "position": row["position"],
                    "p": row["p"],
                    "in_mhc_region": str(in_mhc).lower(),
                    "hit_class": "mhc" if in_mhc else "non_mhc",
                    "note": f"p<={p_threshold}",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize GWAS hits passing a P-value threshold.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    args = parser.parse_args()
    summarize(args.input_csv, args.output_csv, args.p_threshold)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
