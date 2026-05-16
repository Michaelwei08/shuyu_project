from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_hits(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"variant_id", "chrom", "position", "p"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def chrom_sort(chrom: str) -> tuple[int, str]:
    clean = chrom.lower().replace("chr", "")
    if clean.isdigit():
        return int(clean), clean
    return 99, clean


def manhattan(rows: list[dict[str, str]], output: Path) -> None:
    width, height = 920, 460
    left, top, plot_w, plot_h = 80, 90, 760, 270
    ordered = sorted(rows, key=lambda row: (chrom_sort(row["chrom"]), float(row["position"])))
    max_logp = max([-math.log10(float(row["p"])) for row in ordered] or [1.0])
    max_logp = max(max_logp, 1.0)
    denom = max(len(ordered) - 1, 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="920" height="460" fill="#f8fafc"/>',
        '<text x="460" y="44" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">Example GWAS Manhattan plot</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#ffffff" stroke="#cbd5e1"/>',
    ]
    for idx, row in enumerate(ordered):
        x = left + idx / denom * plot_w
        y = top + plot_h - (-math.log10(float(row["p"])) / max_logp * plot_h)
        color = "#ef4444" if row.get("in_mhc_region", "false").lower() == "true" else "#2563eb"
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="6" fill="{color}" opacity="0.88"><title>{row["variant_id"]}</title></circle>')
    parts.extend(
        [
            f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#64748b" stroke-width="2"/>',
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#64748b" stroke-width="2"/>',
            f'<text x="{left + plot_w / 2}" y="{top + plot_h + 45}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">Variants ordered by chromosome and position</text>',
            f'<text x="28" y="{top + plot_h / 2}" transform="rotate(-90 28 {top + plot_h / 2})" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">-log10(P)</text>',
            '<text x="460" y="424" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Red points are MHC-region variants; exclude them from primary MR unless explicitly justified.</text>',
            "</svg>",
        ]
    )
    output.write_text("".join(parts), encoding="utf-8")


def qq(rows: list[dict[str, str]], output: Path) -> None:
    width, height = 520, 520
    left, top, plot_w, plot_h = 80, 80, 360, 320
    observed = sorted([-math.log10(float(row["p"])) for row in rows])
    n = len(observed)
    expected = [-math.log10((idx + 0.5) / max(n, 1)) for idx in range(n)]
    max_value = max(observed + expected + [1.0])
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="520" height="520" fill="#f8fafc"/>',
        '<text x="260" y="42" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="23" font-weight="700" fill="#0f172a">Example GWAS QQ plot</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#ffffff" stroke="#cbd5e1"/>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top}" stroke="#94a3b8" stroke-width="2" stroke-dasharray="6 6"/>',
    ]
    for exp, obs in zip(expected, observed):
        x = left + exp / max_value * plot_w
        y = top + plot_h - obs / max_value * plot_h
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5" fill="#0f766e" opacity="0.88"/>')
    parts.extend(
        [
            f'<text x="{left + plot_w / 2}" y="{top + plot_h + 42}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">Expected -log10(P)</text>',
            f'<text x="24" y="{top + plot_h / 2}" transform="rotate(-90 24 {top + plot_h / 2})" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#334155">Observed -log10(P)</text>',
            "</svg>",
        ]
    )
    output.write_text("".join(parts), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create example Manhattan and QQ SVG plots from GWAS hit tables.")
    parser.add_argument("gwas_hits_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_hits(args.gwas_hits_csv)
    manhattan(rows, args.output_dir / "example_gwas_manhattan.svg")
    qq(rows, args.output_dir / "example_gwas_qq.svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
