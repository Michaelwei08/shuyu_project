from __future__ import annotations

import argparse
import csv
from pathlib import Path


def call_status(row: dict[str, str], threshold: float) -> str:
    value = row.get("ebv_reads_per_million", "")
    if not value:
        return "unknown"
    return "positive" if float(value) >= threshold else "negative"


def compare(input_csv: Path, output_csv: Path, threshold: float) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "ebv_call_status", "ebv_annotation_concordance"])
        writer.writeheader()
        for row in reader:
            output = dict(row)
            call = call_status(row, threshold)
            known = row.get("ebv_status_known", "unknown").strip().lower()
            if known not in {"positive", "negative"} or call == "unknown":
                concordance = "not_assessable"
            elif known == call:
                concordance = "concordant"
            else:
                concordance = "discordant"
            output["ebv_call_status"] = call
            output["ebv_annotation_concordance"] = concordance
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare WGS EBV calls against known EBV annotations.")
    parser.add_argument("joined_calls_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--threshold", type=float, default=1.0)
    args = parser.parse_args()
    compare(args.joined_calls_csv, args.output_csv, args.threshold)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
