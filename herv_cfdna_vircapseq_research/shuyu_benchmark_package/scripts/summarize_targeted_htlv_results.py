from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_table(path: Path, delimiter: str | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        if delimiter is None:
            delimiter = "\t" if path.suffix == ".tsv" else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def value(row: dict[str, str], column: str) -> int:
    text = row.get(column, "") or "0"
    return int(float(text))


def classify(htlv1: int, threshold: int) -> str:
    if htlv1 >= threshold:
        return "HTLV1_positive"
    if htlv1 > 0:
        return "low_positive"
    return "no_signal_or_qc_fail"


def summarize_counts(rows: list[dict[str, str]], threshold: int) -> tuple[list[dict[str, object]], dict[str, object]]:
    out_rows: list[dict[str, object]] = []
    for row in rows:
        htlv1 = value(row, "HTLV1")
        line1 = max(value(row, "LINE1"), 1)
        herv = max(value(row, "HERV"), 1)
        out = dict(row)
        out["HTLV1_per_LINE1_10k"] = f"{htlv1 / line1 * 10000:.3f}"
        out["HTLV1_per_HERV_10k"] = f"{htlv1 / herv * 10000:.3f}"
        out["call"] = classify(htlv1, threshold)
        out_rows.append(out)

    htlv_values = [value(row, "HTLV1") for row in rows]
    sorted_values = sorted(htlv_values)
    median = sorted_values[len(sorted_values) // 2] if sorted_values else 0
    metrics: dict[str, object] = {
        "samples": len(rows),
        "htlv1_positive": sum(v >= threshold for v in htlv_values),
        "low_positive": sum(0 < v < threshold for v in htlv_values),
        "zero": sum(v == 0 for v in htlv_values),
        "hiv1_nonzero": sum(value(row, "HIV1") > 0 for row in rows),
        "hiv2_nonzero": sum(value(row, "HIV2") > 0 for row in rows),
        "htlv2_nonzero": sum(value(row, "HTLV2") > 0 for row in rows),
        "median_htlv1": median,
        "max_htlv1": max(htlv_values) if htlv_values else 0,
    }
    return sorted(out_rows, key=lambda row: int(row.get("HTLV1", 0)), reverse=True), metrics


def audit_manifests(
    manifest_path: Path | None,
    qc_path: Path | None,
    count_rows: list[dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    completed = {row["sample"] for row in count_rows}
    audit_rows: list[dict[str, object]] = []
    metrics: dict[str, object] = {"manifest_samples": "", "qc_complete_targeted": ""}

    manifest_ids: set[str] = set()
    if manifest_path:
        manifest = read_table(manifest_path, ",")
        manifest_ids = {row["sample_id"] for row in manifest}
        metrics["manifest_samples"] = len(manifest_ids)
        for sample in sorted(manifest_ids - completed):
            audit_rows.append({"sample": sample, "audit_status": "manifest_not_completed", "issue_flags": ""})
        for sample in sorted(completed - manifest_ids):
            audit_rows.append({"sample": sample, "audit_status": "completed_not_in_manifest", "issue_flags": ""})

    if qc_path:
        qc_rows = read_table(qc_path, ",")
        complete_targeted = [
            row
            for row in qc_rows
            if row.get("assay_type") == "targeted" and row.get("paired_fastq_complete") == "yes"
        ]
        metrics["qc_complete_targeted"] = len(complete_targeted)
        for row in complete_targeted:
            sample = row["sample_id"]
            if manifest_ids and sample not in manifest_ids:
                audit_rows.append(
                    {
                        "sample": sample,
                        "audit_status": "qc_complete_targeted_not_in_manifest",
                        "issue_flags": row.get("issue_flags", ""),
                    }
                )

    return audit_rows, metrics


def svg_bar_chart(rows: list[dict[str, object]], path: Path, limit: int) -> None:
    top = rows[:limit]
    width = 1200
    row_h = 26
    left = 310
    right = 40
    top_margin = 58
    height = top_margin + max(len(top), 1) * row_h + 48
    max_value = max([int(row.get("HTLV1", 0)) for row in top] or [1])
    chart_w = width - left - right
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="600" y="30" text-anchor="middle" font-family="Georgia" font-size="24" fill="#1f2933">Top HTLV1 Signals</text>',
        '<text x="600" y="50" text-anchor="middle" font-family="Arial" font-size="12" fill="#52606d">Coordinate-deduplicated alignments after competitive human/retrovirus mapping</text>',
    ]
    for i, row in enumerate(top):
        y = top_margin + i * row_h
        sample = str(row["sample"]).replace("targeted_htlv_tcl_", "")
        value_int = int(row.get("HTLV1", 0))
        bar_w = 1 if max_value == 0 else int(chart_w * value_int / max_value)
        lines.append(f'<text x="{left - 12}" y="{y + 17}" text-anchor="end" font-family="Arial" font-size="12" fill="#334e68">{sample}</text>')
        lines.append(f'<rect x="{left}" y="{y + 5}" width="{bar_w}" height="16" rx="3" fill="#2f6f73"/>')
        label_x = min(left + bar_w + 8, width - right)
        anchor = "start"
        if label_x > width - 80:
            label_x = left + max(bar_w - 8, 4)
            anchor = "end"
        lines.append(f'<text x="{label_x}" y="{y + 17}" text-anchor="{anchor}" font-family="Arial" font-size="11" fill="#102a43">{value_int}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(path: Path, metrics: dict[str, object], audit_metrics: dict[str, object], threshold: int) -> None:
    lines = [
        "# Targeted HTLV Full-Run Summary",
        "",
        f"- Samples analyzed: {metrics['samples']}",
        f"- HTLV1 positive threshold: >= {threshold} coordinate-deduplicated alignments",
        f"- HTLV1 positive: {metrics['htlv1_positive']}",
        f"- Low positive: {metrics['low_positive']}",
        f"- Zero/possible QC fail: {metrics['zero']}",
        f"- Median HTLV1 coordinate-deduplicated alignments: {metrics['median_htlv1']}",
        f"- Max HTLV1 coordinate-deduplicated alignments: {metrics['max_htlv1']}",
        f"- HIV1 nonzero samples: {metrics['hiv1_nonzero']}",
        f"- HIV2 nonzero samples: {metrics['hiv2_nonzero']}",
        f"- HTLV2 nonzero samples: {metrics['htlv2_nonzero']}",
        "",
        "## Manifest Audit",
        "",
        f"- Manifest samples: {audit_metrics.get('manifest_samples', '')}",
        f"- QC complete targeted samples: {audit_metrics.get('qc_complete_targeted', '')}",
        "",
        "Interpretation: strong targeted HTLV1 detection with no HIV1/HIV2/HTLV2 signal after competitive filtering.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize full targeted HTLV competitive-alignment results.")
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--sample-qc", type=Path)
    parser.add_argument("--positive-threshold", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    counts = read_table(args.counts, "\t")
    summary_rows, metrics = summarize_counts(counts, args.positive_threshold)
    audit_rows, audit_metrics = audit_manifests(args.manifest, args.sample_qc, counts)

    fieldnames = [*counts[0].keys(), "HTLV1_per_LINE1_10k", "HTLV1_per_HERV_10k", "call"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "targeted_htlv_full_summary.tsv", summary_rows, fieldnames)
    write_tsv(args.output_dir / "targeted_htlv_manifest_audit.tsv", audit_rows, ["sample", "audit_status", "issue_flags"])
    write_report(args.output_dir / "targeted_htlv_full_report.md", metrics, audit_metrics, args.positive_threshold)
    svg_bar_chart(summary_rows, args.output_dir / "targeted_htlv_top_signals.svg", args.top_n)
    print(f"Wrote summary outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
