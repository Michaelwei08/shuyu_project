from __future__ import annotations

import argparse
import csv
from pathlib import Path


METRICS = {
    "hiv": "unique_hiv_per_million",
    "htlv": "unique_htlv_per_million",
    "ambiguous": "ambiguous_retroviral_per_million",
}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def svg_bar_plot(rows: list[dict[str, str]], metric: str, title: str, output: Path) -> None:
    groups = []
    values = []
    for row in rows:
        if row["metric"] != metric:
            continue
        groups.append(row["expected_group"])
        values.append(float(row["mean"]))

    width = 860
    height = 420
    margin_left = 150
    plot_w = 620
    plot_h = 240
    max_value = max(values) if values else 1.0
    max_value = max(max_value, 1.0)
    row_gap = 58
    start_y = 120

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="860" height="420" fill="#f8fafc"/>',
        f'<text x="40" y="48" font-family="Segoe UI, Arial, sans-serif" font-size="26" font-weight="700" fill="#0f172a">{title}</text>',
        f'<text x="40" y="78" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Mean normalized signal by benchmark group. Replace example rows with real benchmark output before interpretation.</text>',
    ]
    for idx, (group, value) in enumerate(zip(groups, values)):
        y = start_y + idx * row_gap
        bar_w = plot_w * value / max_value if max_value else 0
        parts.append(f'<text x="40" y="{y + 20}" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="600" fill="#334155">{group}</text>')
        parts.append(f'<rect x="{margin_left}" y="{y}" width="{plot_w}" height="28" rx="6" fill="#e2e8f0"/>')
        parts.append(f'<rect x="{margin_left}" y="{y}" width="{bar_w}" height="28" rx="6" fill="#2563eb"/>')
        parts.append(f'<text x="{margin_left + plot_w + 18}" y="{y + 20}" font-family="Segoe UI, Arial, sans-serif" font-size="13" font-weight="700" fill="#0f172a">{value:.3f}</text>')
    parts.append("</svg>")
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot normalized retroviral benchmark metrics by group.")
    parser.add_argument("group_summary_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(args.group_summary_csv)
    for name, metric in METRICS.items():
        svg_bar_plot(rows, metric, f"{name.upper()} signal by benchmark group", args.output_dir / f"{name}_signal_by_group.svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
