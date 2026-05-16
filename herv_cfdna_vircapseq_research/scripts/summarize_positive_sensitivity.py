from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def summarize(calls_csv: Path, output_csv: Path) -> None:
    counts: dict[str, Counter[str]] = {"hiv_positive": Counter(), "htlv_positive": Counter(), "healthy_negative": Counter()}
    with calls_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            group = row.get("expected_group", "")
            if group not in counts:
                continue
            counts[group]["n"] += 1
            if row["hiv_call"] == "positive":
                counts[group]["hiv_positive_calls"] += 1
            if row["htlv_call"] == "positive":
                counts[group]["htlv_positive_calls"] += 1
            if row["hiv_call"] == "review" or row["htlv_call"] == "review":
                counts[group]["review_calls"] += 1

    fields = ["expected_group", "n", "hiv_positive_calls", "htlv_positive_calls", "review_calls", "hiv_sensitivity_or_fp_rate", "htlv_sensitivity_or_fp_rate"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for group, counter in counts.items():
            n = counter["n"]
            writer.writerow(
                {
                    "expected_group": group,
                    "n": n,
                    "hiv_positive_calls": counter["hiv_positive_calls"],
                    "htlv_positive_calls": counter["htlv_positive_calls"],
                    "review_calls": counter["review_calls"],
                    "hiv_sensitivity_or_fp_rate": f"{counter['hiv_positive_calls'] / n:.6f}" if n else "",
                    "htlv_sensitivity_or_fp_rate": f"{counter['htlv_positive_calls'] / n:.6f}" if n else "",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize frozen-call sensitivity and healthy-control false positives.")
    parser.add_argument("calls_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    summarize(args.calls_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
