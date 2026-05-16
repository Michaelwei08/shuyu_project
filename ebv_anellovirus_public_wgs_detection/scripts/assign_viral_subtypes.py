from __future__ import annotations

import argparse
import csv
from pathlib import Path


def number(row: dict[str, str], column: str) -> float:
    value = row.get(column, "").strip()
    return float(value) if value else 0.0


def subtype(row: dict[str, str], ebv_threshold: float, anello_threshold: float, pan_threshold: float) -> str:
    ebv = number(row, "ebv_reads_per_million")
    anello = number(row, "anellovirus_reads_per_million")
    herpes = number(row, "other_herpesvirus_reads_per_million")
    other = number(row, "other_viral_reads_per_million")
    viral_families_positive = sum(
        [
            ebv >= ebv_threshold,
            anello >= anello_threshold,
            herpes > 0,
            other > 0,
        ]
    )
    pan = ebv + anello + herpes + other
    if ebv >= ebv_threshold and anello < anello_threshold:
        return "V3_high_ebv_low_anello"
    if viral_families_positive >= 2 and pan >= pan_threshold:
        return "V1_pan_viral"
    if pan == 0:
        return "no_confident_viral_signal"
    return "unclassified_viral_signal"


def assign(input_csv: Path, output_csv: Path, ebv_threshold: float, anello_threshold: float, pan_threshold: float) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle, output_csv.open("w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        writer = csv.DictWriter(out, fieldnames=[*fieldnames, "draft_viral_subtype"])
        writer.writeheader()
        for row in reader:
            output = dict(row)
            output["draft_viral_subtype"] = subtype(row, ebv_threshold, anello_threshold, pan_threshold)
            writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Assign draft viral subtype labels from normalized viral features.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--ebv-threshold", type=float, default=1.0)
    parser.add_argument("--anello-threshold", type=float, default=1.0)
    parser.add_argument("--pan-threshold", type=float, default=2.0)
    args = parser.parse_args()
    assign(args.input_csv, args.output_csv, args.ebv_threshold, args.anello_threshold, args.pan_threshold)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
