from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


SPECIES = ["HTLV1", "HTLV2", "HIV1", "HIV2", "HERV"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def count_value(row: dict[str, str]) -> int:
    return int(float(row["deduplicated_reads"] or 0))


def quantile(values: list[int], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    index = (len(ordered) - 1) * q
    lo = int(index)
    hi = min(lo + 1, len(ordered) - 1)
    frac = index - lo
    return ordered[lo] * (1 - frac) + ordered[hi] * frac


def summarize_figure2(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_sample: dict[tuple[str, str], dict[str, int]] = {}
    labels: dict[str, str] = {}
    for row in rows:
        cohort = row["cohort"]
        sample = row["sample"]
        labels[cohort] = row["cohort_label"]
        by_sample.setdefault((cohort, sample), {})[row["species"]] = count_value(row)

    sample_rows: list[dict[str, object]] = []
    for (cohort, sample), counts in sorted(by_sample.items()):
        herv = counts.get("HERV", 0)
        displayed_total = sum(counts.get(species, 0) for species in SPECIES)
        exogenous_total = sum(counts.get(species, 0) for species in ["HTLV1", "HTLV2", "HIV1", "HIV2"])
        sample_rows.append(
            {
                "cohort": cohort,
                "cohort_label": labels.get(cohort, cohort),
                "sample": sample,
                "HERV": herv,
                "exogenous_total": exogenous_total,
                "displayed_total": displayed_total,
                "herv_fraction_of_displayed": round(herv / displayed_total, 6) if displayed_total else 0,
                "herv_gt_exogenous": int(herv > exogenous_total),
                "herv_positive": int(herv > 0),
            }
        )

    cohort_rows: list[dict[str, object]] = []
    cohorts = sorted({row["cohort"] for row in rows})
    for cohort in cohorts:
        subset = [row for row in sample_rows if row["cohort"] == cohort]
        herv_values = [int(row["HERV"]) for row in subset]
        fractions = [float(row["herv_fraction_of_displayed"]) for row in subset]
        cohort_rows.append(
            {
                "cohort": cohort,
                "cohort_label": labels.get(cohort, cohort),
                "samples": len(subset),
                "herv_positive_samples": sum(1 for value in herv_values if value > 0),
                "herv_positive_fraction": round(sum(1 for value in herv_values if value > 0) / len(subset), 6)
                if subset
                else 0,
                "herv_min": min(herv_values) if herv_values else 0,
                "herv_q25": round(quantile(herv_values, 0.25), 3),
                "herv_median": round(statistics.median(herv_values), 3) if herv_values else 0,
                "herv_mean": round(statistics.mean(herv_values), 3) if herv_values else 0,
                "herv_q75": round(quantile(herv_values, 0.75), 3),
                "herv_max": max(herv_values) if herv_values else 0,
                "median_herv_fraction_of_displayed": round(statistics.median(fractions), 6) if fractions else 0,
                "herv_gt_exogenous_samples": sum(int(row["herv_gt_exogenous"]) for row in subset),
            }
        )
    return sample_rows, cohort_rows


def get_detection(rows: list[dict[str, str]], cohort: str, species: str) -> dict[str, str] | None:
    for row in rows:
        if row["cohort"] == cohort and row["species"] == species:
            return row
    return None


def write_report(path: Path, cohort_rows: list[dict[str, object]], figure3_rows: list[dict[str, str]]) -> None:
    lines = [
        "# HERV Background QC From Local Figure Tables",
        "",
        "## Scope",
        "",
        "- Input is `figure2.tsv` and `figure3.tsv`, so this checks the plotted result tables only.",
        "- LINE1 is not present in these downloaded figure tables; HERV-vs-LINE1 specificity must be tested from full `filtered_category_counts.tsv` on the remote run.",
        "",
        "## Findings",
        "",
    ]
    for row in cohort_rows:
        lines.append(
            "- "
            f"{row['cohort_label']}: HERV-positive {row['herv_positive_samples']}/{row['samples']}; "
            f"median={row['herv_median']}; mean={row['herv_mean']}; "
            f"range={row['herv_min']}-{row['herv_max']}; "
            f"median fraction of displayed retroviral signal={row['median_herv_fraction_of_displayed']}"
        )

    lines.extend(["", "## Detection Table Cross-Check", ""])
    for cohort_row in cohort_rows:
        cohort = str(cohort_row["cohort"])
        det = get_detection(figure3_rows, cohort, "HERV")
        if det:
            lines.append(
                f"- {cohort_row['cohort_label']}: figure3 HERV detection "
                f"{det['detected_samples']}/{det['total_samples']} at threshold {det['threshold']}."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- HERV in nearly every WGS sample is plausible because endogenous retroviral sequence is part of the human genome.",
            "- It is not a final HERV-specific biological call because the current plotted result used a compact HERV-K113 decoy rather than hg38-wide repeat-aware alignment.",
            "- The next specificity test should rerun with hg38 or a richer HERV/RepeatMasker/Dfam decoy set, stricter MAPQ, and coordinate/UMI-aware deduplication.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="QC ubiquitous HERV signal from local Figure 2/3 TSVs.")
    parser.add_argument("--figure2-data", type=Path, required=True)
    parser.add_argument("--figure3-data", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    figure2_rows = read_tsv(args.figure2_data)
    figure3_rows = read_tsv(args.figure3_data)
    sample_rows, cohort_rows = summarize_figure2(figure2_rows)

    write_tsv(
        args.output_dir / "herv_background_sample_qc.tsv",
        sample_rows,
        [
            "cohort",
            "cohort_label",
            "sample",
            "HERV",
            "exogenous_total",
            "displayed_total",
            "herv_fraction_of_displayed",
            "herv_gt_exogenous",
            "herv_positive",
        ],
    )
    write_tsv(
        args.output_dir / "herv_background_cohort_qc.tsv",
        cohort_rows,
        [
            "cohort",
            "cohort_label",
            "samples",
            "herv_positive_samples",
            "herv_positive_fraction",
            "herv_min",
            "herv_q25",
            "herv_median",
            "herv_mean",
            "herv_q75",
            "herv_max",
            "median_herv_fraction_of_displayed",
            "herv_gt_exogenous_samples",
        ],
    )
    write_report(args.output_dir / "herv_background_qc_report.md", cohort_rows, figure3_rows)
    print(f"Wrote HERV background QC outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
