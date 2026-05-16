from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_id",
    "read_id",
    "virus",
    "family",
    "mapq",
    "alignment_length",
    "is_duplicate",
}


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def summarize(path: Path, min_mapq: float, min_length: int) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    seen_reads: set[tuple[str, str, str]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            if truthy(row.get("is_duplicate", "")):
                continue
            if float(row["mapq"]) < min_mapq:
                continue
            if int(float(row["alignment_length"])) < min_length:
                continue
            sample_id = row["sample_id"].strip()
            read_id = row["read_id"].strip()
            virus = row["virus"].strip()
            family = row["family"].strip()
            key = (sample_id, read_id, virus)
            if key in seen_reads:
                continue
            seen_reads.add(key)
            counts[sample_id]["total_viral_reads"] += 1
            if virus.lower() in {"ebv", "epstein-barr virus", "human herpesvirus 4", "hhv-4"}:
                counts[sample_id]["ebv_reads"] += 1
            elif family == "Herpesviridae":
                counts[sample_id]["other_herpesvirus_reads"] += 1
            elif family == "Anelloviridae":
                counts[sample_id]["anellovirus_reads"] += 1
            else:
                counts[sample_id]["other_viral_reads"] += 1
    return counts


def write_summary(counts: dict[str, Counter[str]], output_path: Path) -> None:
    columns = [
        "sample_id",
        "total_viral_reads",
        "ebv_reads",
        "other_herpesvirus_reads",
        "anellovirus_reads",
        "other_viral_reads",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for sample_id in sorted(counts):
            row = {"sample_id": sample_id}
            row.update({column: counts[sample_id].get(column, 0) for column in columns if column != "sample_id"})
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize viral hits into EBV, other herpesvirus, and anellovirus bins.")
    parser.add_argument("hits_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--min-mapq", type=float, default=20)
    parser.add_argument("--min-length", type=int, default=40)
    args = parser.parse_args()

    counts = summarize(args.hits_csv, min_mapq=args.min_mapq, min_length=args.min_length)
    write_summary(counts, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
