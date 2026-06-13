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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def value(row: dict[str, str], column: str) -> int:
    return int(float(row.get(column, "") or 0))


def rate_per_10k(numerator: int, denominator: int, decimals: int = 6) -> str:
    return f"{numerator / denominator * 10000:.{decimals}f}" if denominator > 0 else "NA"


def group_from_manifest(row: dict[str, str]) -> str:
    if row.get("hiv_status") == "positive":
        return "HIV_WGS"
    if row.get("expected_group") == "hl_wgs_hiv_htlv_negative_control":
        return "HL_control"
    return row.get("expected_group") or "unknown"


def median(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def annotate_rows(counts: list[dict[str, str]], manifest: list[dict[str, str]]) -> list[dict[str, object]]:
    manifest_by_sample = {row["sample_id"]: row for row in manifest}
    rows: list[dict[str, object]] = []
    for row in counts:
        sample = row["sample"]
        meta = manifest_by_sample.get(sample, {})
        line1 = value(row, "LINE1")
        herv = value(row, "HERV")
        hiv_total = value(row, "HIV1") + value(row, "HIV2")
        htlv_total = value(row, "HTLV1") + value(row, "HTLV2")
        annotated: dict[str, object] = {
            "sample": sample,
            "subject_id": meta.get("subject_id", ""),
            "group": group_from_manifest(meta),
            "hiv_status": meta.get("hiv_status", ""),
            "expected_group": meta.get("expected_group", ""),
            "HERV": value(row, "HERV"),
            "HIV1": value(row, "HIV1"),
            "HIV2": value(row, "HIV2"),
            "HTLV1": value(row, "HTLV1"),
            "HTLV2": value(row, "HTLV2"),
            "LINE1": value(row, "LINE1"),
            "HIV_total": hiv_total,
            "HTLV_total": htlv_total,
            "HIV_per_LINE1_10k": rate_per_10k(hiv_total, line1),
            "HTLV_per_LINE1_10k": rate_per_10k(htlv_total, line1),
            "HERV_per_LINE1_10k": rate_per_10k(value(row, "HERV"), line1, decimals=3),
            "HIV_per_HERV_10k": rate_per_10k(hiv_total, herv),
        }
        rows.append(annotated)
    return rows


def summarize_groups(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups = sorted({str(row["group"]) for row in rows})
    summary: list[dict[str, object]] = []
    for group in groups:
        subset = [row for row in rows if row["group"] == group]
        summary.append(
            {
                "group": group,
                "samples": len(subset),
                "HIV_total_nonzero": sum(int(row["HIV_total"]) > 0 for row in subset),
                "HTLV_total_nonzero": sum(int(row["HTLV_total"]) > 0 for row in subset),
                "median_HIV_total": f"{median([int(row['HIV_total']) for row in subset]):.1f}",
                "median_HTLV_total": f"{median([int(row['HTLV_total']) for row in subset]):.1f}",
                "median_HERV": f"{median([int(row['HERV']) for row in subset]):.1f}",
                "median_LINE1": f"{median([int(row['LINE1']) for row in subset]):.1f}",
                "max_HIV_total": max([int(row["HIV_total"]) for row in subset] or [0]),
                "max_HTLV_total": max([int(row["HTLV_total"]) for row in subset] or [0]),
            }
        )
    return summary


def write_report(path: Path, rows: list[dict[str, object]], group_rows: list[dict[str, object]]) -> None:
    hiv_group = [row for row in rows if row["group"] == "HIV_WGS"]
    hl_group = [row for row in rows if row["group"] == "HL_control"]
    lines = [
        "# WGS HIV/HL Retrovirus Screen Summary",
        "",
        f"- Samples analyzed: {len(rows)}",
        f"- HIV WGS samples: {len(hiv_group)}",
        f"- HL negative-control samples: {len(hl_group)}",
        f"- HIV signal in HIV WGS samples: {sum(int(row['HIV_total']) > 0 for row in hiv_group)}",
        f"- HIV signal in HL controls: {sum(int(row['HIV_total']) > 0 for row in hl_group)}",
        f"- HTLV signal in HL controls: {sum(int(row['HTLV_total']) > 0 for row in hl_group)}",
        "",
        "## Group Summary",
        "",
        "| group | samples | HIV_total_nonzero | HTLV_total_nonzero | median_HERV | median_LINE1 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in group_rows:
        lines.append(
            f"| {row['group']} | {row['samples']} | {row['HIV_total_nonzero']} | "
            f"{row['HTLV_total_nonzero']} | {row['median_HERV']} | {row['median_LINE1']} |"
        )
    lines.extend(
        [
            "",
            "Interpretation guide: HL controls should remain HIV/HTLV zero. HIV WGS may remain zero at low subsampling depth; that supports specificity but does not prove HIV sensitivity.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def svg_scatter(rows: list[dict[str, object]], path: Path) -> None:
    width, height = 1000, 560
    left, right, top, bottom = 90, 40, 70, 80
    plot_w, plot_h = width - left - right, height - top - bottom
    max_line1 = max([int(row["LINE1"]) for row in rows] or [1])
    max_herv = max([int(row["HERV"]) for row in rows] or [1])
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="500" y="34" text-anchor="middle" font-family="Georgia" font-size="24" fill="#1f2933">WGS HERV/LINE1 Background</text>',
        '<text x="500" y="56" text-anchor="middle" font-family="Arial" font-size="12" fill="#52606d">Competitive reference screen; colored by expected WGS group</text>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#9fb3c8"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#9fb3c8"/>',
        f'<text x="{left + plot_w/2:.0f}" y="{height-28}" text-anchor="middle" font-family="Arial" font-size="13" fill="#334e68">LINE1 filtered read IDs</text>',
        f'<text x="22" y="{top + plot_h/2:.0f}" text-anchor="middle" transform="rotate(-90 22 {top + plot_h/2:.0f})" font-family="Arial" font-size="13" fill="#334e68">HERV filtered read IDs</text>',
        '<circle cx="785" cy="32" r="6" fill="#c2410c"/><text x="800" y="37" font-family="Arial" font-size="12">HIV WGS</text>',
        '<circle cx="875" cy="32" r="6" fill="#2563eb"/><text x="890" y="37" font-family="Arial" font-size="12">HL control</text>',
    ]
    for row in rows:
        x = left + (int(row["LINE1"]) / max_line1) * plot_w if max_line1 else left
        y = height - bottom - (int(row["HERV"]) / max_herv) * plot_h if max_herv else height - bottom
        color = "#c2410c" if row["group"] == "HIV_WGS" else "#2563eb"
        stroke = "#111827" if int(row["HIV_total"]) or int(row["HTLV_total"]) else "none"
        lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="{color}" stroke="{stroke}" stroke-width="1.5"/>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize WGS HIV/HL competitive retrovirus screen.")
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    counts = read_table(args.counts, "\t")
    manifest = read_table(args.manifest, ",")
    rows = annotate_rows(counts, manifest)
    group_rows = summarize_groups(rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.output_dir / "wgs_retro_sample_summary.tsv",
        rows,
        [
            "sample",
            "subject_id",
            "group",
            "hiv_status",
            "expected_group",
            "HERV",
            "HIV1",
            "HIV2",
            "HTLV1",
            "HTLV2",
            "LINE1",
            "HIV_total",
            "HTLV_total",
            "HIV_per_LINE1_10k",
            "HTLV_per_LINE1_10k",
            "HERV_per_LINE1_10k",
            "HIV_per_HERV_10k",
        ],
    )
    write_tsv(args.output_dir / "wgs_retro_group_summary.tsv", group_rows, list(group_rows[0].keys()) if group_rows else ["group"])
    write_report(args.output_dir / "wgs_retro_report.md", rows, group_rows)
    svg_scatter(rows, args.output_dir / "wgs_herv_line1_background.svg")
    print(f"Wrote WGS summary outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
