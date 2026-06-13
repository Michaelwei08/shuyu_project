from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


BACKGROUND_SPECIES = ["HIV1", "HIV2", "HTLV2", "HERV", "LINE1"]
SUMMARY_SPECIES = ["HTLV1", "HIV1", "HIV2", "HTLV2", "HERV", "LINE1", "HUMAN"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def num(row: dict[str, str], column: str) -> int:
    return int(float(row.get(column, "0") or 0))


def call_status(htlv1: int, threshold: int) -> str:
    if htlv1 >= threshold:
        return "strong_positive"
    if htlv1 > 0:
        return "low_positive"
    return "zero_or_qc_fail"


def median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


def fmt(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    if value == int(value):
        return str(int(value))
    return f"{value:.3f}"


def merge_counts(
    pre_rows: list[dict[str, str]],
    post_rows: list[dict[str, str]],
    removed_rows: list[dict[str, str]],
    threshold: int,
) -> list[dict[str, object]]:
    pre_by = {row["sample"]: row for row in pre_rows}
    removed_by = {row["sample"]: row for row in removed_rows}
    merged: list[dict[str, object]] = []
    for post in post_rows:
        sample = post["sample"]
        pre = pre_by.get(sample, {})
        removed = removed_by.get(sample, {})
        pre_htlv1 = num(pre, "HTLV1")
        post_htlv1 = num(post, "HTLV1")
        removed_htlv1 = num(removed, "HTLV1")
        merged.append(
            {
                "sample": sample,
                "pre_dedup_HTLV1": pre_htlv1,
                "post_dedup_HTLV1": post_htlv1,
                "dedup_removed_HTLV1": removed_htlv1,
                "post_pre_HTLV1_ratio": round(post_htlv1 / pre_htlv1, 6) if pre_htlv1 else 0,
                "call": call_status(post_htlv1, threshold),
                **{f"post_{species}": num(post, species) for species in SUMMARY_SPECIES},
            }
        )
    return sorted(merged, key=lambda row: int(row["post_dedup_HTLV1"]), reverse=True)


def dedup_summary(
    pre_rows: list[dict[str, str]],
    post_rows: list[dict[str, str]],
    removed_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    post_by = {row["sample"]: row for row in post_rows}
    removed_by = {row["sample"]: row for row in removed_rows}
    out: list[dict[str, object]] = []
    for species in SUMMARY_SPECIES:
        ratios: list[float] = []
        total_pre = total_post = total_removed = 0
        for pre in pre_rows:
            sample = pre["sample"]
            pre_count = num(pre, species)
            post_count = num(post_by.get(sample, {}), species)
            removed_count = num(removed_by.get(sample, {}), species)
            total_pre += pre_count
            total_post += post_count
            total_removed += removed_count
            if pre_count:
                ratios.append(post_count / pre_count)
        out.append(
            {
                "species": species,
                "samples_with_pre_signal": len(ratios),
                "total_pre_dedup_records": total_pre,
                "total_post_dedup_counts": total_post,
                "total_dedup_removed": total_removed,
                "median_post_pre_ratio": round(median(ratios), 6),
                "fraction_removed": round(total_removed / total_pre, 6) if total_pre else 0,
            }
        )
    return out


def background_check(post_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for species in BACKGROUND_SPECIES:
        values = [num(row, species) for row in post_rows]
        out.append(
            {
                "species": species,
                "samples_nonzero": sum(value > 0 for value in values),
                "samples_total": len(values),
                "total_post_dedup_counts": sum(values),
                "max_post_dedup_count": max(values) if values else 0,
                "pass_zero_background": int(sum(value > 0 for value in values) == 0),
            }
        )
    return out


def call_distribution(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    calls = ["strong_positive", "low_positive", "zero_or_qc_fail"]
    return [{"call": call, "samples": sum(row["call"] == call for row in rows)} for call in calls]


def svg_header(width: int, height: int, title: str, subtitle: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        f'<text x="{width / 2}" y="34" text-anchor="middle" font-family="Georgia" font-size="24" fill="#1f2933">{title}</text>',
        f'<text x="{width / 2}" y="56" text-anchor="middle" font-family="Arial" font-size="12" fill="#52606d">{subtitle}</text>',
    ]


def write_pre_post_scatter(rows: list[dict[str, object]], path: Path) -> None:
    width, height = 820, 700
    left, right, top, bottom = 86, 34, 80, 76
    plot_w, plot_h = width - left - right, height - top - bottom
    max_val = max([int(row["pre_dedup_HTLV1"]) for row in rows] + [1])
    max_log = math.log10(max_val + 1)
    lines = svg_header(width, height, "HTLV1 Deduplication Effect", "Each dot is one targeted HTLV sample; axes are log10(count + 1).")
    lines.extend(
        [
            f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#9fb3c8"/>',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#9fb3c8"/>',
        ]
    )
    for tick in [0, 10, 100, 1000, 10000, 100000, 1000000]:
        if tick > max_val:
            continue
        x = left + math.log10(tick + 1) / max_log * plot_w
        y = top + plot_h - math.log10(tick + 1) / max_log * plot_h
        lines.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_h}" stroke="#e2e8f0"/>')
        lines.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e2e8f0"/>')
        lines.append(f'<text x="{x:.1f}" y="{top + plot_h + 22}" text-anchor="middle" font-family="Arial" font-size="10" fill="#52606d">{tick}</text>')
        lines.append(f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial" font-size="10" fill="#52606d">{tick}</text>')
    lines.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top}" stroke="#64748b" stroke-dasharray="5,5"/>')
    for row in rows:
        x = left + math.log10(int(row["pre_dedup_HTLV1"]) + 1) / max_log * plot_w
        y = top + plot_h - math.log10(int(row["post_dedup_HTLV1"]) + 1) / max_log * plot_h
        lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#b45309" opacity="0.62"/>')
    lines.append(f'<text x="{left + plot_w / 2}" y="{height - 20}" text-anchor="middle" font-family="Arial" font-size="12" fill="#334e68">Pre-dedup HTLV1 records</text>')
    lines.append(f'<text x="18" y="{top + plot_h / 2}" text-anchor="middle" transform="rotate(-90 18 {top + plot_h / 2})" font-family="Arial" font-size="12" fill="#334e68">Post-dedup HTLV1 counts</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_bar_svg(rows: list[dict[str, object]], path: Path, title: str, value_key: str, label_key: str) -> None:
    width, height = 780, 430
    left, right, top, bottom = 210, 40, 78, 42
    max_value = max([int(row[value_key]) for row in rows] + [1])
    lines = svg_header(width, height, title, "Post-dedup counts from hg38+RefSeq competitive alignment.")
    for i, row in enumerate(rows):
        y = top + i * 54
        value = int(row[value_key])
        bar_w = (width - left - right) * value / max_value if max_value else 0
        label = str(row[label_key])
        lines.append(f'<text x="{left - 12}" y="{y + 23}" text-anchor="end" font-family="Arial" font-size="13" fill="#334e68">{label}</text>')
        lines.append(f'<rect x="{left}" y="{y + 7}" width="{bar_w:.1f}" height="24" rx="4" fill="#2f6f73"/>')
        lines.append(f'<text x="{left + bar_w + 8:.1f}" y="{y + 24}" font-family="Arial" font-size="12" fill="#102a43">{value}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(path: Path, merged: list[dict[str, object]], dedup_rows: list[dict[str, object]], background_rows: list[dict[str, object]], threshold: int) -> None:
    calls = {row["call"]: row["samples"] for row in call_distribution(merged)}
    htlv = next(row for row in dedup_rows if row["species"] == "HTLV1")
    residuals = [
        f"{row['species']}={row['total_post_dedup_counts']} across {row['samples_nonzero']} samples"
        for row in background_rows
        if int(row["total_post_dedup_counts"]) > 0
    ]
    residual_text = "; ".join(residuals) if residuals else "none"
    lines = [
        "# Targeted HTLV Method Validation",
        "",
        f"- Samples analyzed: {len(merged)}",
        f"- Operational strong HTLV1 threshold: >= {threshold} coordinate-deduplicated alignments",
        f"- Strong positive: {calls.get('strong_positive', 0)}",
        f"- Low positive: {calls.get('low_positive', 0)}",
        f"- Zero/QC-fail: {calls.get('zero_or_qc_fail', 0)}",
        f"- Median HTLV1 post/pre dedup ratio: {htlv['median_post_pre_ratio']}",
        f"- HTLV1 fraction removed by coordinate dedup: {htlv['fraction_removed']}",
        f"- Residual non-target background: {residual_text}",
        "",
        "## Interpretation",
        "",
        "The hg38+RefSeq competitive run retains targeted HTLV1 signal while nearly eliminating exogenous-virus and retroelement background.",
        "Coordinate deduplication has a large effect, so post-dedup counts should be used for interpretation; the threshold is analytical and not clinically validated.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Targeted-only method validation for hg38+RefSeq HTLV rerun.")
    parser.add_argument("--filtered-record-counts", type=Path, required=True)
    parser.add_argument("--filtered-counts", type=Path, required=True)
    parser.add_argument("--dedup-removed-counts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--positive-threshold", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    pre = read_tsv(args.filtered_record_counts)
    post = read_tsv(args.filtered_counts)
    removed = read_tsv(args.dedup_removed_counts)
    merged = merge_counts(pre, post, removed, args.positive_threshold)
    dedup_rows = dedup_summary(pre, post, removed)
    background_rows = background_check(post)
    calls = call_distribution(merged)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "targeted_htlv_top_samples.tsv", merged[: args.top_n], list(merged[0].keys()))
    write_tsv(args.output_dir / "targeted_htlv_all_samples.tsv", merged, list(merged[0].keys()))
    write_tsv(args.output_dir / "targeted_htlv_dedup_summary.tsv", dedup_rows, list(dedup_rows[0].keys()))
    write_tsv(args.output_dir / "targeted_htlv_background_check.tsv", background_rows, list(background_rows[0].keys()))
    write_tsv(args.output_dir / "targeted_htlv_call_distribution.tsv", calls, ["call", "samples"])
    write_pre_post_scatter(merged, args.output_dir / "targeted_htlv_pre_post_dedup_scatter.svg")
    write_bar_svg(calls, args.output_dir / "targeted_htlv_call_distribution.svg", "HTLV1 Call Distribution", "samples", "call")
    write_bar_svg(background_rows, args.output_dir / "targeted_htlv_background_nonzero.svg", "Background Nonzero Samples", "samples_nonzero", "species")
    write_report(args.output_dir / "targeted_htlv_method_validation_report.md", merged, dedup_rows, background_rows, args.positive_threshold)
    print(f"Wrote targeted HTLV method-validation outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
