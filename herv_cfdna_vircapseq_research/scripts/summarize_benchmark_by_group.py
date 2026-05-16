from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METRICS = [
    "unique_hiv_per_million",
    "unique_htlv_per_million",
    "unique_herv_line1_per_million",
    "ambiguous_retroviral_per_million",
]


def load_groups(manifest_csv: Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    with manifest_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "sample_id" not in (reader.fieldnames or []) or "expected_group" not in (reader.fieldnames or []):
            raise ValueError("Manifest must contain sample_id and expected_group")
        for row in reader:
            groups[row["sample_id"]] = row["expected_group"]
    return groups


def summarize(manifest_csv: Path, bins_csv: Path, output_csv: Path) -> None:
    groups = load_groups(manifest_csv)
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    with bins_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = row["sample_id"]
            group = groups.get(sample_id, "unknown")
            for metric in METRICS:
                value = row.get(metric, "")
                if value:
                    values[(group, metric)].append(float(value))

    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["expected_group", "metric", "n", "mean", "min", "max"])
        writer.writeheader()
        for group, metric in sorted(values):
            metric_values = values[(group, metric)]
            writer.writerow(
                {
                    "expected_group": group,
                    "metric": metric,
                    "n": len(metric_values),
                    "mean": f"{sum(metric_values) / len(metric_values):.6f}",
                    "min": f"{min(metric_values):.6f}",
                    "max": f"{max(metric_values):.6f}",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize normalized retroviral benchmark metrics by expected group.")
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("bins_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    summarize(args.manifest_csv, args.bins_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
