from __future__ import annotations

import argparse
import csv
import html
import re
from collections import defaultdict
from pathlib import Path


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "unknown"


def load_coverage(path: Path) -> dict[tuple[str, str], list[tuple[int, float]]]:
    grouped: dict[tuple[str, str], list[tuple[int, float]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "virus", "position", "depth"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            grouped[(row["sample_id"], row["virus"])].append((int(row["position"]), float(row["depth"])))
    for values in grouped.values():
        values.sort()
    return grouped


def render_plot(sample_id: str, virus: str, values: list[tuple[int, float]], output: Path) -> None:
    width, height = 900, 420
    left, top, plot_w, plot_h = 90, 90, 720, 240
    max_position = max((position for position, _ in values), default=1)
    min_position = min((position for position, _ in values), default=1)
    max_depth = max((depth for _, depth in values), default=1.0)
    max_depth = max(max_depth, 1.0)
    span = max(max_position - min_position, 1)
    bar_w = max(2, min(18, plot_w / max(len(values), 1) * 0.65))

    title = f"{virus} coverage distribution: {sample_id}"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="900" height="420" fill="#f8fafc"/>',
        f'<text x="450" y="42" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">{html.escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#ffffff" stroke="#cbd5e1"/>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#64748b" stroke-width="2"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#64748b" stroke-width="2"/>',
    ]

    for position, depth in values:
        x = left + ((position - min_position) / span) * plot_w
        bar_h = (depth / max_depth) * plot_h
        y = top + plot_h - bar_h
        parts.append(
            f'<rect x="{x - bar_w / 2:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" fill="#0f766e" opacity="0.82">'
            f'<title>position {position}: depth {depth:g}</title></rect>'
        )

    tick_y = top + plot_h + 24
    parts.extend(
        [
            f'<text x="{left}" y="{tick_y}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">{min_position}</text>',
            f'<text x="{left + plot_w}" y="{tick_y}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">{max_position}</text>',
            f'<text x="{left + plot_w / 2}" y="{top + plot_h + 58}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="600" fill="#334155">Viral reference position</text>',
            f'<text x="{left - 18}" y="{top + 6}" text-anchor="end" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">{max_depth:g}</text>',
            f'<text x="{left - 18}" y="{top + plot_h + 4}" text-anchor="end" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">0</text>',
            f'<text x="28" y="{top + plot_h / 2}" transform="rotate(-90 28 {top + plot_h / 2})" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="14" font-weight="600" fill="#334155">Depth</text>',
            '<text x="450" y="388" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Use this plot to spot narrow pileups, broad viral coverage, or reference-specific artifacts.</text>',
            "</svg>",
        ]
    )
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot per-position viral coverage distributions as SVG files.")
    parser.add_argument("coverage_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    grouped = load_coverage(args.coverage_csv)
    for (sample_id, virus), values in sorted(grouped.items()):
        output = args.output_dir / f"coverage_{safe_name(virus)}_{safe_name(sample_id)}.svg"
        render_plot(sample_id, virus, values, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
