from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


SPECIES = ["HTLV1", "HTLV2", "HIV1", "HIV2", "HERV"]
COHORTS = ["targeted_htlv", "HIV_WGS", "HL_control"]
COHORT_COLORS = {
    "targeted_htlv": "#b45309",
    "HIV_WGS": "#c2410c",
    "HL_control": "#2563eb",
}
SPECIES_COLORS = {
    "HTLV1": "#b45309",
    "HTLV2": "#92400e",
    "HIV1": "#c2410c",
    "HIV2": "#dc2626",
    "HERV": "#475569",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_svg(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def text(value: object) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def count_value(row: dict[str, str]) -> int:
    return int(float(row["deduplicated_reads"] or 0))


def pct(value: float) -> str:
    if value == 0:
        return "0%"
    if value == 1:
        return "100%"
    return f"{value * 100:.1f}%"


def median(values: list[int]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2


def mean(values: list[int]) -> float:
    return sum(values) / len(values) if values else 0.0


def format_count(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:.0f}"


def stable_jitter(seed: str, width: float) -> float:
    raw = hashlib.sha256(seed.encode("utf-8")).hexdigest()[:8]
    return (int(raw, 16) / 0xFFFFFFFF - 0.5) * width


def log_y(value: int, max_value: int, top: float, height: float) -> float:
    return top + height - (math.log10(value + 1) / math.log10(max_value + 1)) * height


def figure2(rows: list[dict[str, str]], output: Path) -> None:
    cohort_labels = {row["cohort"]: row["cohort_label"] for row in rows}
    max_value = max(count_value(row) for row in rows)
    width = 1320
    panel_h = 392
    left = 104
    right = 42
    header = 86
    top_pad = 78
    bottom_pad = 70
    plot_w = width - left - right
    height = header + len(COHORTS) * panel_h + 28
    ticks = [0, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]
    ticks = [tick for tick in ticks if tick <= max_value]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="660" y="36" text-anchor="middle" font-family="Georgia" font-size="28" fill="#1f2933">Figure 2. Deduplicated Reads by Viral Species and Cohort</text>',
        '<text x="660" y="62" text-anchor="middle" font-family="Arial" font-size="13" fill="#52606d">Dots are samples; horizontal bars show medians; teal diamonds and labels show averages. Y axis is log10(count + 1).</text>',
        f'<text x="{left}" y="{header - 8}" font-family="Arial" font-size="12" font-weight="700" fill="#334e68">Deduplicated read count</text>',
    ]
    for panel, cohort in enumerate(COHORTS):
        y0 = header + panel * panel_h
        plot_top = y0 + top_pad
        plot_h = panel_h - top_pad - bottom_pad
        plot_bottom = plot_top + plot_h
        label = cohort_labels.get(cohort, cohort)
        cohort_rows = [row for row in rows if row["cohort"] == cohort]
        sample_count = len({row["sample"] for row in cohort_rows})
        lines.append(f'<text x="{left}" y="{y0 + 31}" font-family="Arial" font-size="19" font-weight="700" fill="#102a43">{text(label)}</text>')
        lines.append(f'<text x="{left}" y="{y0 + 54}" font-family="Arial" font-size="12" fill="#627d98">n = {sample_count} samples</text>')
        lines.append(f'<line x1="{left}" y1="{plot_bottom}" x2="{left + plot_w}" y2="{plot_bottom}" stroke="#9fb3c8"/>')
        lines.append(f'<line x1="{left}" y1="{plot_top}" x2="{left}" y2="{plot_bottom}" stroke="#9fb3c8"/>')
        for tick in ticks:
            y = log_y(tick, max_value, plot_top, plot_h)
            lines.append(f'<line x1="{left - 5}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e2e8f0"/>')
            lines.append(f'<text x="{left - 12}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial" font-size="10" fill="#52606d">{tick}</text>')
        for idx, species in enumerate(SPECIES):
            x = left + (idx + 0.5) * plot_w / len(SPECIES)
            subset = [row for row in cohort_rows if row["species"] == species]
            values = [count_value(row) for row in subset]
            med = median(values)
            avg = mean(values)
            y_med = log_y(int(med), max_value, plot_top, plot_h)
            lines.append(f'<rect x="{x - 40}" y="{y_med:.1f}" width="80" height="{plot_bottom - y_med:.1f}" rx="4" fill="{SPECIES_COLORS[species]}" opacity="0.16"/>')
            lines.append(f'<line x1="{x - 43}" y1="{y_med:.1f}" x2="{x + 43}" y2="{y_med:.1f}" stroke="{SPECIES_COLORS[species]}" stroke-width="3"/>')
            if med > 0:
                label_y = max(y_med - 8, plot_top + 12)
                lines.append(f'<text x="{x}" y="{label_y:.1f}" text-anchor="middle" font-family="Arial" font-size="10" fill="#102a43">median {med:.0f}</text>')
            for row in subset:
                val = count_value(row)
                px = x + stable_jitter(row["sample"] + species, 50)
                py = log_y(val, max_value, plot_top, plot_h)
                radius = 3.2 if val else 2.2
                opacity = 0.58 if val else 0.28
                lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{radius}" fill="{COHORT_COLORS[cohort]}" opacity="{opacity}"/>')
            y_avg = log_y(int(avg), max_value, plot_top, plot_h)
            diamond = f"{x:.1f},{y_avg - 6:.1f} {x + 6:.1f},{y_avg:.1f} {x:.1f},{y_avg + 6:.1f} {x - 6:.1f},{y_avg:.1f}"
            lines.append(f'<polygon points="{diamond}" fill="#0f766e" stroke="#064e3b" stroke-width="1.2"/>')
            lines.append(f'<text x="{x}" y="{plot_bottom + 28}" text-anchor="middle" font-family="Arial" font-size="13" fill="#334e68">{species}</text>')
            lines.append(f'<text x="{x}" y="{plot_bottom + 47}" text-anchor="middle" font-family="Arial" font-size="10" fill="#0f766e">avg {format_count(avg)}</text>')
    lines.append("</svg>")
    write_svg(output, lines)


def figure3(rows: list[dict[str, str]], output: Path) -> None:
    width = 1220
    height = 680
    left = 92
    right = 42
    top = 112
    bottom = 112
    plot_w = width - left - right
    plot_h = height - top - bottom
    cohort_labels = {row["cohort"]: row["cohort_label"] for row in rows}
    lookup = {(row["cohort"], row["species"]): row for row in rows}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf6"/>',
        '<text x="610" y="38" text-anchor="middle" font-family="Georgia" font-size="28" fill="#1f2933">Figure 3. Detection Proportions by Viral Species and Cohort</text>',
        '<text x="610" y="64" text-anchor="middle" font-family="Arial" font-size="13" fill="#52606d">Detection is deduplicated read count &gt;= 1. Labels show detected / total samples.</text>',
        f'<text x="{left}" y="{top - 24}" font-family="Arial" font-size="12" font-weight="700" fill="#334e68">Proportion detected</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#9fb3c8"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#9fb3c8"/>',
    ]
    for tick in [0, 25, 50, 75, 100]:
        y = top + plot_h - tick / 100 * plot_h
        lines.append(f'<line x1="{left - 6}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e2e8f0"/>')
        lines.append(f'<text x="{left - 14}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial" font-size="11" fill="#52606d">{tick}%</text>')
    group_w = plot_w / len(SPECIES)
    bar_w = group_w / 4.8
    for idx, species in enumerate(SPECIES):
        cx = left + idx * group_w + group_w / 2
        lines.append(f'<text x="{cx}" y="{top + plot_h + 32}" text-anchor="middle" font-family="Arial" font-size="13" fill="#334e68">{species}</text>')
        for j, cohort in enumerate(COHORTS):
            row = lookup[(cohort, species)]
            prop = float(row["proportion"])
            detected = int(row["detected_samples"])
            total = int(row["total_samples"])
            h = prop * plot_h
            x = cx - 1.5 * bar_w + j * bar_w
            y = top + plot_h - h
            lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w - 5:.1f}" height="{h:.1f}" rx="4" fill="{COHORT_COLORS[cohort]}"/>')
            label = f"{detected}/{total}"
            label_y = y - 7 if h > 22 else y - 9
            if h == 0:
                label_y = top + plot_h - 6
            lines.append(f'<text x="{x + (bar_w - 5) / 2:.1f}" y="{label_y:.1f}" text-anchor="middle" font-family="Arial" font-size="10" fill="#102a43">{label}</text>')
    legend_y = height - 52
    legend_x = 250
    for j, cohort in enumerate(COHORTS):
        x = legend_x + j * 250
        lines.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" rx="2" fill="{COHORT_COLORS[cohort]}"/>')
        lines.append(f'<text x="{x + 22}" y="{legend_y + 12}" font-family="Arial" font-size="12" fill="#334e68">{text(cohort_labels.get(cohort, cohort))}</text>')
    lines.append("</svg>")
    write_svg(output, lines)


def validate_svg(path: Path) -> None:
    data = path.read_text(encoding="utf-8")
    if data.count("<svg") != 1 or not data.rstrip().endswith("</svg>"):
        raise SystemExit(f"Invalid SVG structure: {path}")
    if "NaN" in data or "nan" in data:
        raise SystemExit(f"Invalid numeric value in SVG: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Render local Figure 2 and Figure 3 from downloaded plot-data TSVs.")
    parser.add_argument("--figure2-data", type=Path, required=True)
    parser.add_argument("--figure3-data", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    fig2_rows = read_tsv(args.figure2_data)
    fig3_rows = read_tsv(args.figure3_data)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig2 = args.output_dir / "figure2_deduplicated_reads_local.svg"
    fig3 = args.output_dir / "figure3_detection_proportions_local.svg"
    figure2(fig2_rows, fig2)
    figure3(fig3_rows, fig3)
    validate_svg(fig2)
    validate_svg(fig3)
    print(f"Wrote {fig2}")
    print(f"Wrote {fig3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
