from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


CATEGORIES = ["HIV1", "HIV2", "HTLV1", "HTLV2", "HERV", "LINE1"]


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def count(row: dict[str, str], category: str) -> int:
    return int(float(row.get(category, "") or 0))


def index_counts(path: Path) -> dict[str, dict[str, str]]:
    rows = read_table(path)
    return {row["sample"]: row for row in rows}


def compare_pair(old_path: Path, new_path: Path, cohort: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    old = index_counts(old_path)
    new = index_counts(new_path)
    if set(old) != set(new):
        missing_new = sorted(set(old) - set(new))
        missing_old = sorted(set(new) - set(old))
        raise ValueError(f"Sample mismatch for {cohort}: missing_new={missing_new[:5]}, missing_old={missing_old[:5]}")

    samples: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    for category in CATEGORIES:
        old_values = [count(old[sample], category) for sample in sorted(old)]
        new_values = [count(new[sample], category) for sample in sorted(new)]
        for sample, old_value, new_value in zip(sorted(old), old_values, new_values):
            samples.append(
                {
                    "cohort": cohort,
                    "sample": sample,
                    "category": category,
                    "old_count": old_value,
                    "new_count": new_value,
                    "new_to_old_ratio": f"{new_value / old_value:.6f}" if old_value else "NA",
                }
            )
        summary.append(
            {
                "cohort": cohort,
                "category": category,
                "samples": len(old_values),
                "old_nonzero_samples": sum(value > 0 for value in old_values),
                "new_nonzero_samples": sum(value > 0 for value in new_values),
                "old_total": sum(old_values),
                "new_total": sum(new_values),
                "old_median": statistics.median(old_values),
                "new_median": statistics.median(new_values),
                "new_to_old_total_ratio": f"{sum(new_values) / sum(old_values):.6f}" if sum(old_values) else "NA",
            }
        )
    return samples, summary


def write_report(path: Path, rows: list[dict[str, object]]) -> None:
    lines = [
        "# Viral-Only vs hg38-Inclusive Method Comparison",
        "",
        "The new method uses hg38 competition, RefSeq retroviral references, category-specific MAPQ thresholds, unique-best filtering, and coordinate deduplication.",
        "",
        "| cohort | category | old nonzero | new nonzero | old total | new total | new/old total |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['cohort']} | {row['category']} | {row['old_nonzero_samples']} | "
            f"{row['new_nonzero_samples']} | {row['old_total']} | {row['new_total']} | "
            f"{row['new_to_old_total_ratio']} |"
        )
    lines.extend(
        [
            "",
            "Interpret exogenous-virus retention together with negative-control specificity. Loss of HERV/LINE1 after adding hg38 indicates reclassification of reads as human, not biological absence of those elements.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare viral-only and hg38-inclusive retrovirus count tables.")
    parser.add_argument("--old-targeted", type=Path, required=True)
    parser.add_argument("--new-targeted", type=Path, required=True)
    parser.add_argument("--old-wgs", type=Path, required=True)
    parser.add_argument("--new-wgs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    targeted_samples, targeted_summary = compare_pair(args.old_targeted, args.new_targeted, "targeted_HTLV")
    wgs_samples, wgs_summary = compare_pair(args.old_wgs, args.new_wgs, "HIV_HL_WGS")
    samples = targeted_samples + wgs_samples
    summary = targeted_summary + wgs_summary
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "method_comparison_by_sample.tsv", samples, list(samples[0].keys()))
    write_tsv(args.output_dir / "method_comparison_summary.tsv", summary, list(summary[0].keys()))
    write_report(args.output_dir / "method_comparison_report.md", summary)
    print(f"Wrote old-versus-new method comparison to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
