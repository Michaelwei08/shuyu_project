from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def plot(rows: list[dict[str, str]], output: Path) -> None:
    selected = [row for row in rows if row["scope"] == "all_diseases" and row["n_samples"] != "0"]
    width, height = 880, 460
    left, top = 120, 110
    panel_w, panel_h = 280, 240
    max_ebv = max([float(row["mean_ebv_load"]) for row in selected] or [1.0])
    max_anello = max([float(row["mean_anellovirus_load"]) for row in selected] or [1.0])
    max_ebv = max(max_ebv, 1.0)
    max_anello = max(max_anello, 1.0)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="880" height="460" fill="#f8fafc"/>',
        '<text x="440" y="42" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">Example V1/V3 interpretation features</text>',
        '<text x="440" y="70" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Toy data only; use real viral calls and clinical metadata before interpretation.</text>',
        f'<rect x="{left}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#cbd5e1"/>',
        f'<rect x="{left + 360}" y="{top}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#cbd5e1"/>',
        f'<text x="{left + panel_w / 2}" y="{top - 18}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="15" font-weight="700" fill="#334155">Mean EBV load</text>',
        f'<text x="{left + 360 + panel_w / 2}" y="{top - 18}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="15" font-weight="700" fill="#334155">Mean anellovirus load</text>',
    ]
    for idx, row in enumerate(selected):
        y = top + 55 + idx * 90
        ebv_w = panel_w * float(row["mean_ebv_load"]) / max_ebv
        anello_w = panel_w * float(row["mean_anellovirus_load"]) / max_anello
        label = row["subtype"].replace("_", " ")
        parts.append(f'<text x="24" y="{y + 18}" font-family="Segoe UI, Arial, sans-serif" font-size="13" font-weight="700" fill="#334155">{label}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{ebv_w:.2f}" height="24" rx="6" fill="#2563eb"/>')
        parts.append(f'<rect x="{left + 360}" y="{y}" width="{anello_w:.2f}" height="24" rx="6" fill="#0f766e"/>')
        parts.append(f'<text x="{left + ebv_w + 8:.2f}" y="{y + 17}" font-family="Segoe UI, Arial, sans-serif" font-size="12" fill="#0f172a">{float(row["mean_ebv_load"]):.2f}</text>')
        parts.append(f'<text x="{left + 360 + anello_w + 8:.2f}" y="{y + 17}" font-family="Segoe UI, Arial, sans-serif" font-size="12" fill="#0f172a">{float(row["mean_anellovirus_load"]):.2f}</text>')
    parts.append("</svg>")
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot subtype interpretation features from summary table.")
    parser.add_argument("subtype_summary_csv", type=Path)
    parser.add_argument("output_svg", type=Path)
    args = parser.parse_args()
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    plot(load_rows(args.subtype_summary_csv), args.output_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
