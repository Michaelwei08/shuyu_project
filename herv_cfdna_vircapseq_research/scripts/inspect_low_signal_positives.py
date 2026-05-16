from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def top_categories(hits_csv: Path) -> dict[str, str]:
    categories: dict[str, set[str]] = defaultdict(set)
    with hits_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            categories[row["sample_id"]].add(row["reference_category"].lower())
    return {sample_id: "|".join(sorted(values)) for sample_id, values in categories.items()}


def inspect(manifest_csv: Path, calls_csv: Path, hits_csv: Path, output_csv: Path) -> None:
    calls = {row["sample_id"]: row for row in read_rows(calls_csv)}
    hit_categories = top_categories(hits_csv)
    fields = ["sample_id", "expected_group", "failed_expected_virus", "observed_call", "observed_hit_categories", "review_note"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in read_rows(manifest_csv):
            sample_id = row["sample_id"]
            call = calls.get(sample_id)
            if call is None:
                continue
            expected_group = row["expected_group"]
            if expected_group == "hiv_positive" and call["hiv_call"] != "positive":
                writer.writerow(
                    {
                        "sample_id": sample_id,
                        "expected_group": expected_group,
                        "failed_expected_virus": "hiv",
                        "observed_call": call["hiv_call"],
                        "observed_hit_categories": hit_categories.get(sample_id, ""),
                        "review_note": "low HIV signal after frozen filtering; check assay sensitivity, subtype/reference mismatch, and WGS depth",
                    }
                )
            if expected_group == "htlv_positive" and call["htlv_call"] != "positive":
                writer.writerow(
                    {
                        "sample_id": sample_id,
                        "expected_group": expected_group,
                        "failed_expected_virus": "htlv",
                        "observed_call": call["htlv_call"],
                        "observed_hit_categories": hit_categories.get(sample_id, ""),
                        "review_note": "low HTLV signal after frozen filtering; check HTLV type, reference mismatch, and target compatibility",
                    }
                )


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect known-positive samples missed by frozen HIV/HTLV calls.")
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("calls_csv", type=Path)
    parser.add_argument("hits_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    inspect(args.manifest_csv, args.calls_csv, args.hits_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
