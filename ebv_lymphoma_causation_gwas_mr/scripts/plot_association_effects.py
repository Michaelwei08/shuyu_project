from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def plot(rows: list[dict[str, str]], output: Path) -> None:
    selected = [
        row
        for row in rows
        if row["test_name"] == "disease_vs_control"
        and row["model"] == "unadjusted"
        and row["phenotype"] in {"ebv_load", "anellovirus_load"}
        and row["effect_estimate"]
    ]
    width, height = 980, 520
    left, top, plot_w, row_h = 360, 90, 460, 34
    values = [float(row["effect_estimate"]) for row in selected]
    max_abs = max([abs(value) for value in values] or [1.0])
    max_abs = max(max_abs, 1.0)
    zero_x = left + plot_w / 2
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="980" height="520" fill="#f8fafc"/>',
        '<text x="490" y="42" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">Example viral association effect sizes</text>',
        '<text x="490" y="68" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Unadjusted disease-vs-control scaffold effects; replace example data before interpretation.</text>',
        f'<line x1="{zero_x}" y1="{top - 10}" x2="{zero_x}" y2="{top + row_h * max(len(selected), 1) + 10}" stroke="#94a3b8" stroke-width="2"/>',
    ]
    for idx, row in enumerate(selected):
        y = top + idx * row_h
        value = float(row["effect_estimate"])
        ci_low = float(row["ci95_low"]) if row["ci95_low"] else value
        ci_high = float(row["ci95_high"]) if row["ci95_high"] else value
        x = zero_x + (value / max_abs) * (plot_w / 2)
        x_low = zero_x + (ci_low / max_abs) * (plot_w / 2)
        x_high = zero_x + (ci_high / max_abs) * (plot_w / 2)
        label = f'{row["group"]} {row["phenotype"]}'
        color = "#2563eb" if row["phenotype"] == "ebv_load" else "#0f766e"
        parts.append(f'<text x="40" y="{y + 7}" font-family="Segoe UI, Arial, sans-serif" font-size="13" font-weight="600" fill="#334155">{label}</text>')
        parts.append(f'<line x1="{x_low:.2f}" y1="{y}" x2="{x_high:.2f}" y2="{y}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<circle cx="{x:.2f}" cy="{y}" r="6" fill="{color}"/>')
        parts.append(f'<text x="840" y="{y + 5}" font-family="Segoe UI, Arial, sans-serif" font-size="12" fill="#0f172a">{value:.2f}</text>')
    parts.append("</svg>")
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot scaffold association effect sizes from association-grid output.")
    parser.add_argument("association_csv", type=Path)
    parser.add_argument("output_svg", type=Path)
    args = parser.parse_args()
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    plot(load_rows(args.association_csv), args.output_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
