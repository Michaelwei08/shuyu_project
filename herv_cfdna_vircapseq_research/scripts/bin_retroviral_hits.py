from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_id",
    "read_id",
    "reference_category",
    "mapq",
    "alignment_score",
    "alignment_length",
}

RETROVIRAL = {"hiv", "htlv", "herv", "line1"}


@dataclass(frozen=True)
class Hit:
    sample_id: str
    read_id: str
    category: str
    mapq: float
    score: float
    length: int


def parse_float(value: str) -> float:
    return float(value.strip()) if value.strip() else 0.0


def parse_int(value: str) -> int:
    return int(float(value.strip())) if value.strip() else 0


def load_hits(path: Path, min_mapq: float, min_length: int) -> dict[tuple[str, str], list[Hit]]:
    grouped: dict[tuple[str, str], list[Hit]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            hit = Hit(
                sample_id=row["sample_id"].strip(),
                read_id=row["read_id"].strip(),
                category=row["reference_category"].strip().lower(),
                mapq=parse_float(row["mapq"]),
                score=parse_float(row["alignment_score"]),
                length=parse_int(row["alignment_length"]),
            )
            if not hit.sample_id or not hit.read_id:
                continue
            if hit.mapq < min_mapq or hit.length < min_length:
                continue
            grouped[(hit.sample_id, hit.read_id)].append(hit)
    return grouped


def assign_bin(hits: list[Hit], score_delta: float) -> str:
    if not hits:
        return "unassigned"
    ranked = sorted(hits, key=lambda hit: (hit.score, hit.mapq, hit.length), reverse=True)
    best = ranked[0]
    tied = [hit for hit in ranked if best.score - hit.score <= score_delta]
    tied_categories = {hit.category for hit in tied}

    if len(tied_categories & RETROVIRAL) > 1:
        return "ambiguous_retroviral"
    if best.category == "hiv":
        return "unique_hiv"
    if best.category == "htlv":
        return "unique_htlv"
    if best.category in {"herv", "line1"}:
        return "unique_herv_line1"
    if best.category == "other_viral":
        return "other_viral"
    if best.category == "human":
        return "human"
    return "unassigned"


def summarize(grouped_hits: dict[tuple[str, str], list[Hit]], score_delta: float) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for (sample_id, _read_id), hits in grouped_hits.items():
        counts[sample_id][assign_bin(hits, score_delta)] += 1
    return counts


def write_summary(counts: dict[str, Counter[str]], output_path: Path) -> None:
    bins = [
        "unique_hiv",
        "unique_htlv",
        "unique_herv_line1",
        "ambiguous_retroviral",
        "other_viral",
        "human",
        "unassigned",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "total_binned_reads", *bins])
        writer.writeheader()
        for sample_id in sorted(counts):
            total = sum(counts[sample_id].values())
            row = {"sample_id": sample_id, "total_binned_reads": total}
            row.update({bin_name: counts[sample_id].get(bin_name, 0) for bin_name in bins})
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description="Bin competitive retroviral alignment hits into benchmark categories.")
    parser.add_argument("hits_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--min-mapq", type=float, default=20)
    parser.add_argument("--min-length", type=int, default=40)
    parser.add_argument("--score-delta", type=float, default=5)
    args = parser.parse_args()

    grouped = load_hits(args.hits_csv, min_mapq=args.min_mapq, min_length=args.min_length)
    counts = summarize(grouped, score_delta=args.score_delta)
    write_summary(counts, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
