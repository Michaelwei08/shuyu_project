from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


SPECIES = ["HTLV1", "HTLV2", "HIV1", "HIV2", "HERV"]
COHORT_LABELS = {
    "targeted_htlv": "Targeted HTLV TCL",
    "HIV_WGS": "WGS HIV+ DLBCL/BL",
    "HL_control": "WGS HL controls",
}
COLORS = {
    "targeted_htlv": "#b45309",
    "HIV_WGS": "#c2410c",
    "HL_control": "#2563eb",
    "HTLV1": "#b45309",
    "HTLV2": "#92400e",
    "HIV1": "#c2410c",
    "HIV2": "#dc2626",
    "HERV": "#475569",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def num(row: dict[str, str], key: str) -> int:
    return int(float(row.get(key, "") or 0))


def load_plot_rows(targeted_summary: Path, wgs_summary: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(targeted_summary):
        for species in SPECIES:
            rows.append(
                {
                    "cohort": "targeted_htlv",
                    "cohort_label": COHORT_LABELS["targeted_htlv"],
                    "sample": row["sample"],
                    "species": species,
                    "deduplicated_reads": num(row, species),
                }
            )
    for row in read_tsv(wgs_summary):
        cohort = row["group"]
        for species in SPECIES:
            rows.append(
                {
                    "cohort": cohort,
                    "cohort_label": COHORT_LABELS.get(cohort, cohort),
                    "sample": row["sample"],
                    "species": species,
                    "deduplicated_reads": num(row, species),
                }
            )
    return rows


def median(values: list[int]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2


def stable_jitter(text: str, width: float) -> float:
    value = int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16) / 0xFFFFFFFF
    return (value - 0.5) * width


def esc(text: object) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def save_svg(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def figure1(path: Path) -> None:
    width, height = 1500, 760
    boxes = [
        ("Data inputs", "Targeted HTLV TCL FASTQ\nWGS HIV+/HL FASTQ", 60, 110),
        ("Manifest + QC", "prepare_shuyu_benchmark.py\nqc_shuyu_manifest.py", 300, 110),
        ("Competitive reference", "HIV1/HIV2/HTLV1/HTLV2\nHERV/LINE1 decoys\nbwa index", 540, 110),
        ("Alignment", "run_retro_pilot_alignment.py\nbwa mem | samtools sort\n--full-input --jobs", 780, 110),
        ("Filtering + counting", "MAPQ >= 20\naligned length >= 60 bp\ndeduplicate read_id + category", 1020, 110),
        ("Results + figures", "summarize_*.py\naudit_filtered_alignments.py\nplot_p1_figures.py", 1260, 110),
    ]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="750" y="48" text-anchor="middle" font-family="Georgia" font-size="30" fill="#1f2933">Figure 1. P1 Computational Workflow</text>',
        '<text x="750" y="78" text-anchor="middle" font-family="Arial" font-size="15" fill="#52606d">Competitive alignment separates exogenous retroviruses from HERV/LINE1 background</text>',
    ]
    for i, (title, body, x, y) in enumerate(boxes):
        lines.append(f'<rect x="{x}" y="{y}" width="190" height="170" rx="18" fill="#ffffff" stroke="#bcccdc" stroke-width="2"/>')
        lines.append(f'<text x="{x + 95}" y="{y + 32}" text-anchor="middle" font-family="Arial" font-size="17" font-weight="700" fill="#102a43">{esc(title)}</text>')
        for j, part in enumerate(body.split("\n")):
            lines.append(f'<text x="{x + 95}" y="{y + 68 + j * 24}" text-anchor="middle" font-family="Arial" font-size="13" fill="#334e68">{esc(part)}</text>')
        if i < len(boxes) - 1:
            x1 = x + 195
            x2 = boxes[i + 1][2] - 8
            lines.append(f'<line x1="{x1}" y1="{y + 85}" x2="{x2}" y2="{y + 85}" stroke="#486581" stroke-width="2"/>')
            lines.append(f'<polygon points="{x2},{y + 85} {x2 - 12},{y + 78} {x2 - 12},{y + 92}" fill="#486581"/>')
    metrics = [
        ("Targeted HTLV", "71 QC-pass samples; 67 strong HTLV1 positive; 0 HIV1/HIV2/HTLV2 off-target calls"),
        ("Full WGS HIV/HL", "60 samples; 2/37 HIV-labeled WGS with single-fragment HIV1; 0/23 HL controls with HIV/HTLV"),
        ("Audit rule", "No read sequence/quality in audit tables; report only compact alignment metadata"),
    ]
    y0 = 360
    for i, (title, body) in enumerate(metrics):
        y = y0 + i * 95
        lines.append(f'<rect x="120" y="{y}" width="1260" height="70" rx="14" fill="#f0f4f8" stroke="#d9e2ec"/>')
        lines.append(f'<text x="150" y="{y + 28}" font-family="Arial" font-size="17" font-weight="700" fill="#102a43">{esc(title)}</text>')
        lines.append(f'<text x="150" y="{y + 52}" font-family="Arial" font-size="14" fill="#334e68">{esc(body)}</text>')
    lines.append("</svg>")
    save_svg(path, lines)


def y_log(value: int, max_value: int, top: float, height: float) -> float:
    return top + height - (math.log10(value + 1) / math.log10(max_value + 1)) * height


def figure2(rows: list[dict[str, object]], path: Path) -> None:
    cohorts = ["targeted_htlv", "HIV_WGS", "HL_control"]
    width, panel_h = 1200, 355
    left, right, top_pad, bottom_pad = 95, 35, 70, 60
    panel_w = width - left - right
    height = 70 + len(cohorts) * panel_h
    max_value = max(int(row["deduplicated_reads"]) for row in rows)
    ticks = [0, 1, 10, 100, 1000, 10000, 100000, 1000000]
    ticks = [t for t in ticks if t <= max_value] + ([max_value] if max_value not in ticks else [])
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="600" y="36" text-anchor="middle" font-family="Georgia" font-size="28" fill="#1f2933">Figure 2. Deduplicated Viral/Retroelement Reads by Cohort</text>',
        '<text x="600" y="60" text-anchor="middle" font-family="Arial" font-size="13" fill="#52606d">Points are samples; translucent bars show cohort medians. Y-axis uses log10(count + 1) with count labels.</text>',
    ]
    for p, cohort in enumerate(cohorts):
        y0 = 75 + p * panel_h
        plot_top = y0 + top_pad
        plot_h = panel_h - top_pad - bottom_pad
        lines.append(f'<text x="{left}" y="{y0 + 28}" font-family="Arial" font-size="19" font-weight="700" fill="#102a43">{COHORT_LABELS[cohort]}</text>')
        lines.append(f'<line x1="{left}" y1="{plot_top + plot_h}" x2="{left + panel_w}" y2="{plot_top + plot_h}" stroke="#9fb3c8"/>')
        lines.append(f'<line x1="{left}" y1="{plot_top}" x2="{left}" y2="{plot_top + plot_h}" stroke="#9fb3c8"/>')
        for tick in ticks:
            y = y_log(tick, max_value, plot_top, plot_h)
            lines.append(f'<line x1="{left - 5}" y1="{y:.1f}" x2="{left + panel_w}" y2="{y:.1f}" stroke="#e2e8f0"/>')
            lines.append(f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial" font-size="10" fill="#52606d">{tick}</text>')
        for i, species in enumerate(SPECIES):
            x = left + (i + 0.5) * panel_w / len(SPECIES)
            subset = [r for r in rows if r["cohort"] == cohort and r["species"] == species]
            vals = [int(r["deduplicated_reads"]) for r in subset]
            med = median(vals)
            y_med = y_log(int(med), max_value, plot_top, plot_h)
            lines.append(f'<rect x="{x - 32}" y="{y_med:.1f}" width="64" height="{plot_top + plot_h - y_med:.1f}" fill="{COLORS[species]}" opacity="0.18"/>')
            lines.append(f'<line x1="{x - 35}" y1="{y_med:.1f}" x2="{x + 35}" y2="{y_med:.1f}" stroke="{COLORS[species]}" stroke-width="3"/>')
            for r in subset:
                val = int(r["deduplicated_reads"])
                px = x + stable_jitter(str(r["sample"]) + species, 42)
                py = y_log(val, max_value, plot_top, plot_h)
                lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3.2" fill="{COLORS[cohort]}" opacity="0.52"/>')
            lines.append(f'<text x="{x}" y="{plot_top + plot_h + 22}" text-anchor="middle" font-family="Arial" font-size="13" fill="#334e68">{species}</text>')
    lines.append(f'<text x="22" y="{height / 2:.1f}" text-anchor="middle" transform="rotate(-90 22 {height / 2:.1f})" font-family="Arial" font-size="14" fill="#334e68">Deduplicated read count</text>')
    lines.append("</svg>")
    save_svg(path, lines)


def detection_rows(rows: list[dict[str, object]], threshold: int) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for cohort in ["targeted_htlv", "HIV_WGS", "HL_control"]:
        samples = sorted({str(r["sample"]) for r in rows if r["cohort"] == cohort})
        for species in SPECIES:
            vals = [int(r["deduplicated_reads"]) for r in rows if r["cohort"] == cohort and r["species"] == species]
            detected = sum(v >= threshold for v in vals)
            out.append(
                {
                    "cohort": cohort,
                    "cohort_label": COHORT_LABELS[cohort],
                    "species": species,
                    "detected_samples": detected,
                    "total_samples": len(samples),
                    "proportion": detected / len(samples) if samples else 0,
                    "threshold": threshold,
                }
            )
    return out


def figure3(rows: list[dict[str, object]], path: Path) -> None:
    width, height = 1100, 620
    left, right, top, bottom = 90, 40, 80, 95
    plot_w, plot_h = width - left - right, height - top - bottom
    cohorts = ["targeted_htlv", "HIV_WGS", "HL_control"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="550" y="36" text-anchor="middle" font-family="Georgia" font-size="28" fill="#1f2933">Figure 3. Detection Proportions by Cohort</text>',
        '<text x="550" y="60" text-anchor="middle" font-family="Arial" font-size="13" fill="#52606d">Detection defined as deduplicated read count &gt; 0 unless a higher threshold is requested.</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#9fb3c8"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#9fb3c8"/>',
    ]
    for pct in [0, 25, 50, 75, 100]:
        y = top + plot_h - pct / 100 * plot_h
        lines.append(f'<line x1="{left - 5}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e2e8f0"/>')
        lines.append(f'<text x="{left - 12}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial" font-size="11" fill="#52606d">{pct}%</text>')
    group_w = plot_w / len(SPECIES)
    bar_w = group_w / 5
    for i, species in enumerate(SPECIES):
        cx = left + i * group_w + group_w / 2
        lines.append(f'<text x="{cx}" y="{top + plot_h + 28}" text-anchor="middle" font-family="Arial" font-size="13" fill="#334e68">{species}</text>')
        for j, cohort in enumerate(cohorts):
            rec = next(r for r in rows if r["cohort"] == cohort and r["species"] == species)
            prop = float(rec["proportion"])
            h = prop * plot_h
            x = cx - 1.5 * bar_w + j * bar_w
            y = top + plot_h - h
            lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w - 3:.1f}" height="{h:.1f}" fill="{COLORS[cohort]}" rx="3"/>')
            if prop > 0:
                lines.append(f'<text x="{x + (bar_w - 3)/2:.1f}" y="{y - 5:.1f}" text-anchor="middle" font-family="Arial" font-size="10" fill="#102a43">{prop * 100:.0f}%</text>')
    for j, cohort in enumerate(cohorts):
        x = 720 + j * 125
        lines.append(f'<rect x="{x}" y="84" width="12" height="12" fill="{COLORS[cohort]}"/>')
        lines.append(f'<text x="{x + 18}" y="95" font-family="Arial" font-size="12" fill="#334e68">{COHORT_LABELS[cohort]}</text>')
    lines.append(f'<text x="24" y="{top + plot_h / 2:.1f}" text-anchor="middle" transform="rotate(-90 24 {top + plot_h / 2:.1f})" font-family="Arial" font-size="14" fill="#334e68">Proportion detected</text>')
    lines.append("</svg>")
    save_svg(path, lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate P1 workflow, read-count, and detection-proportion figures.")
    parser.add_argument("--targeted-summary", type=Path, required=True)
    parser.add_argument("--wgs-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--detection-threshold", type=int, default=1)
    args = parser.parse_args()

    rows = load_plot_rows(args.targeted_summary, args.wgs_summary)
    det_rows = detection_rows(rows, args.detection_threshold)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "figure2_plot_data.tsv", rows, ["cohort", "cohort_label", "sample", "species", "deduplicated_reads"])
    write_tsv(args.output_dir / "figure3_detection_proportions.tsv", det_rows, ["cohort", "cohort_label", "species", "detected_samples", "total_samples", "proportion", "threshold"])
    figure1(args.output_dir / "figure1_workflow.svg")
    figure2(rows, args.output_dir / "figure2_deduplicated_reads_by_cohort.svg")
    figure3(det_rows, args.output_dir / "figure3_detection_proportions.svg")
    print(f"Wrote figures and plot data to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
