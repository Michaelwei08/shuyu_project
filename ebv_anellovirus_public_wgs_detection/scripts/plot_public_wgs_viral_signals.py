from __future__ import annotations

import argparse
import csv
from pathlib import Path


PLOTS = {
    "ebv": "ebv_reads_per_million",
    "anellovirus": "anellovirus_reads_per_million",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def plot_by_disease(rows: list[dict[str, str]], metric: str, title: str, output: Path) -> None:
    grouped: dict[str, list[float]] = {}
    for row in rows:
        disease = row.get("disease", "unknown") or "unknown"
        value = row.get(metric, "")
        if value:
            grouped.setdefault(disease, []).append(float(value))
    width, height = 900, 420
    plot_w = 620
    x0 = 180
    y0 = 120
    gap = 60
    max_value = max([max(values) for values in grouped.values()] or [1])
    max_value = max(max_value, 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="900" height="420" fill="#f8fafc"/>',
        f'<text x="40" y="48" font-family="Segoe UI, Arial, sans-serif" font-size="26" font-weight="700" fill="#0f172a">{title}</text>',
    ]
    for idx, disease in enumerate(sorted(grouped)):
        values = grouped[disease]
        mean_value = sum(values) / len(values)
        y = y0 + idx * gap
        bar_w = plot_w * mean_value / max_value
        parts.append(f'<text x="40" y="{y + 20}" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="600" fill="#334155">{disease}</text>')
        parts.append(f'<rect x="{x0}" y="{y}" width="{plot_w}" height="28" rx="6" fill="#e2e8f0"/>')
        parts.append(f'<rect x="{x0}" y="{y}" width="{bar_w}" height="28" rx="6" fill="#0f766e"/>')
        parts.append(f'<text x="{x0 + plot_w + 18}" y="{y + 20}" font-family="Segoe UI, Arial, sans-serif" font-size="13" font-weight="700" fill="#0f172a">{mean_value:.3f}</text>')
    parts.append("</svg>")
    output.write_text("".join(parts), encoding="utf-8")


def plot_ebv_vs_anello(rows: list[dict[str, str]], output: Path) -> None:
    width, height = 620, 520
    plot_x, plot_y, plot_w, plot_h = 80, 90, 440, 320
    points = []
    for row in rows:
        ebv = row.get("ebv_reads_per_million", "")
        anello = row.get("anellovirus_reads_per_million", "")
        if ebv and anello:
            points.append((float(ebv), float(anello), row.get("disease", "unknown")))
    max_x = max([point[0] for point in points] or [1])
    max_y = max([point[1] for point in points] or [1])
    max_x = max(max_x, 1)
    max_y = max(max_y, 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="620" height="520" fill="#f8fafc"/>',
        '<text x="40" y="48" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">EBV vs anellovirus signal</text>',
        f'<line x1="{plot_x}" y1="{plot_y + plot_h}" x2="{plot_x + plot_w}" y2="{plot_y + plot_h}" stroke="#94a3b8" stroke-width="2"/>',
        f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y + plot_h}" stroke="#94a3b8" stroke-width="2"/>',
    ]
    for ebv, anello, disease in points:
        x = plot_x + ebv / max_x * plot_w
        y = plot_y + plot_h - anello / max_y * plot_h
        parts.append(f'<circle cx="{x}" cy="{y}" r="7" fill="#2563eb" opacity="0.85"><title>{disease}</title></circle>')
    parts.extend(
        [
            f'<text x="{plot_x + plot_w / 2}" y="{plot_y + plot_h + 48}" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155" text-anchor="middle">EBV reads per million</text>',
            f'<text x="24" y="{plot_y + plot_h / 2}" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">Anello</text>',
            "</svg>",
        ]
    )
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create simple SVG plots for public WGS viral signals.")
    parser.add_argument("joined_or_subtype_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(args.joined_or_subtype_csv)
    for name, metric in PLOTS.items():
        plot_by_disease(rows, metric, f"{name.upper()} signal by disease", args.output_dir / f"{name}_signal_by_disease.svg")
    plot_ebv_vs_anello(rows, args.output_dir / "ebv_vs_anellovirus.svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
