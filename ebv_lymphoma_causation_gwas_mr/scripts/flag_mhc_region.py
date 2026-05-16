from __future__ import annotations

import argparse
import csv
from pathlib import Path


MHC_CHROM = {"6", "chr6"}
MHC_START = 25_000_000
MHC_END = 34_000_000


def flag(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        required = {"chrom", "position"}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "in_mhc_region"])
        writer.writeheader()
        for row in reader:
            chrom = row["chrom"]
            pos = int(float(row["position"]))
            output = dict(row)
            output["in_mhc_region"] = str(chrom in MHC_CHROM and MHC_START <= pos <= MHC_END).lower()
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Flag variants in the extended MHC region.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    flag(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
