from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def number(value: str) -> float | None:
    if value in {"", "NA", "N/A", "na", "null"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def variance(values: list[float]) -> float:
    if len(values) < 2:
        return float("nan")
    center = mean(values)
    return sum((value - center) ** 2 for value in values) / (len(values) - 1)


def welch_t(group_a: list[float], group_b: list[float]) -> tuple[float, float]:
    if len(group_a) < 2 or len(group_b) < 2:
        return float("nan"), float("nan")
    mean_a = mean(group_a)
    mean_b = mean(group_b)
    var_a = variance(group_a)
    var_b = variance(group_b)
    se = math.sqrt(var_a / len(group_a) + var_b / len(group_b))
    if se == 0:
        return mean_a - mean_b, float("nan")
    return mean_a - mean_b, (mean_a - mean_b) / se


def run(input_csv: Path, outcome_column: str, value_column: str, positive_label: str, output_csv: Path) -> None:
    positive: list[float] = []
    negative: list[float] = []
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            value = number(row.get(value_column, ""))
            if value is None:
                continue
            if row.get(outcome_column, "") == positive_label:
                positive.append(value)
            else:
                negative.append(value)
    effect, t_stat = welch_t(positive, negative)
    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "outcome_column",
                "positive_label",
                "value_column",
                "n_positive",
                "n_reference",
                "mean_positive",
                "mean_reference",
                "mean_difference",
                "welch_t_statistic",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "outcome_column": outcome_column,
                "positive_label": positive_label,
                "value_column": value_column,
                "n_positive": len(positive),
                "n_reference": len(negative),
                "mean_positive": f"{mean(positive):.6f}" if positive else "",
                "mean_reference": f"{mean(negative):.6f}" if negative else "",
                "mean_difference": f"{effect:.6f}" if not math.isnan(effect) else "",
                "welch_t_statistic": f"{t_stat:.6f}" if not math.isnan(t_stat) else "",
            }
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a simple two-group association summary for viral phenotypes.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--outcome-column", required=True)
    parser.add_argument("--positive-label", required=True)
    parser.add_argument("--value-column", required=True)
    args = parser.parse_args()
    run(args.input_csv, args.outcome_column, args.value_column, args.positive_label, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
