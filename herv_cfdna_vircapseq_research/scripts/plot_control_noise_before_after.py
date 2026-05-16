from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


ORDERED_METRICS = [
    "unique_hiv_per_million",
    "unique_htlv_per_million",
    "ambiguous_retroviral_per_million",
]

STAGE_COLORS = {
    "before_filtering": "#b91c1c",
    "after_filtering": "#0f766e",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"stage", "metric", "reads_per_million"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def label(metric: str) -> str:
    return metric.replace("_per_million", "").replace("_", " ").upper()


def render(rows: list[dict[str, str]], output: Path) -> None:
    values: dict[tuple[str, str], float] = {}
    for row in rows:
        values[(row["metric"], row["stage"])] = float(row["reads_per_million"])

    width, height = 980, 520
    left, top = 270, 120
    plot_w, bar_h, gap = 560, 24, 58
    max_value = max(values.values() or [1.0])
    max_value = max(max_value, 1.0)
    title = "Healthy-control retroviral noise before and after filtering"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="980" height="520" fill="#f8fafc"/>',
        f'<text x="490" y="46" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="25" font-weight="700" fill="#0f172a">{html.escape(title)}</text>',
        '<text x="490" y="76" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Example control-calibration view; replace with real healthy-control benchmark output before interpretation.</text>',
    ]
    for idx, metric in enumerate(ORDERED_METRICS):
        y = top + idx * (2 * bar_h + gap)
        parts.append(f'<text x="38" y="{y + 34}" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="700" fill="#334155">{html.escape(label(metric))}</text>')
        for stage_idx, stage in enumerate(["before_filtering", "after_filtering"]):
            value = values.get((metric, stage), 0.0)
            bar_w = plot_w * value / max_value
            bar_y = y + stage_idx * (bar_h + 8)
            text_y = bar_y + 17
            stage_text = "before" if stage == "before_filtering" else "after"
            parts.append(f'<text x="{left - 18}" y="{text_y}" text-anchor="end" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">{stage_text}</text>')
            parts.append(f'<rect x="{left}" y="{bar_y}" width="{plot_w}" height="{bar_h}" rx="6" fill="#e2e8f0"/>')
            parts.append(f'<rect x="{left}" y="{bar_y}" width="{bar_w:.2f}" height="{bar_h}" rx="6" fill="{STAGE_COLORS[stage]}"/>')
            parts.append(f'<text x="{left + plot_w + 18}" y="{text_y}" font-family="Segoe UI, Arial, sans-serif" font-size="13" font-weight="700" fill="#0f172a">{value:.2f}</text>')
    parts.extend(
        [
            '<text x="490" y="474" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Primary goal: drive exogenous HIV/HTLV-like signal toward zero in healthy controls while retaining positive-cohort sensitivity.</text>',
            "</svg>",
        ]
    )
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot healthy-control retroviral noise before and after filtering.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_svg", type=Path)
    args = parser.parse_args()
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    render(load_rows(args.input_csv), args.output_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
