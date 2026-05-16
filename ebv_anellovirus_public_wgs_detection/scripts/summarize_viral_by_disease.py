from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METRICS = ["ebv_reads_per_million", "anellovirus_reads_per_million", "other_herpesvirus_reads_per_million"]


def summarize(input_csv: Path, output_csv: Path) -> None:
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            disease = row.get("disease", "unknown") or "unknown"
            for metric in METRICS:
                value = row.get(metric, "")
                if value:
                    values[(disease, metric)].append(float(value))
    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["disease", "metric", "n", "mean", "min", "max"])
        writer.writeheader()
        for disease, metric in sorted(values):
            metric_values = values[(disease, metric)]
            writer.writerow(
                {
                    "disease": disease,
                    "metric": metric,
                    "n": len(metric_values),
                    "mean": f"{sum(metric_values) / len(metric_values):.6f}",
                    "min": f"{min(metric_values):.6f}",
                    "max": f"{max(metric_values):.6f}",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize viral loads by disease.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    summarize(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
