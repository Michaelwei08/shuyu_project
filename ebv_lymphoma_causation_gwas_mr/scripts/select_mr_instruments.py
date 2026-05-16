from __future__ import annotations

import argparse
import csv
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "yes", "1", "y"}


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def row_p(row: dict[str, str]) -> float:
    return float(row.get("p", row.get("p_exposure", "nan")))


def row_position(row: dict[str, str]) -> int | None:
    value = row.get("position", "")
    if not value:
        return None
    return int(float(value))


def clump_by_window(rows: list[dict[str, str]], window_bp: int) -> list[dict[str, str]]:
    kept: list[dict[str, str]] = []
    occupied: dict[str, list[int]] = {}
    for row in sorted(rows, key=row_p):
        chrom = row.get("chrom", "")
        pos = row_position(row)
        if pos is None:
            kept.append(row)
            continue
        if any(abs(pos - prior) <= window_bp for prior in occupied.get(chrom, [])):
            continue
        occupied.setdefault(chrom, []).append(pos)
        kept.append(row)
    return kept


def select(input_csv: Path, output_csv: Path, p_threshold: float, exclude_mhc: bool, window_bp: int) -> None:
    fields, rows = read_rows(input_csv)
    required = {"variant_id", "effect_allele", "other_allele", "beta", "se"}
    missing = required.difference(fields)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    selected: list[dict[str, str]] = []
    for row in rows:
        if row_p(row) > p_threshold:
            continue
        if exclude_mhc and truthy(row.get("in_mhc_region", "false")):
            continue
        if float(row["se"]) <= 0:
            continue
        selected.append(row)
    selected = clump_by_window(selected, window_bp)

    out_fields = [
        "variant_id",
        "chrom",
        "position",
        "effect_allele",
        "other_allele",
        "beta_exposure",
        "se_exposure",
        "p_exposure",
        "in_mhc_region",
        "instrument_source",
        "selection_note",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        for row in selected:
            writer.writerow(
                {
                    "variant_id": row["variant_id"],
                    "chrom": row.get("chrom", ""),
                    "position": row.get("position", ""),
                    "effect_allele": row["effect_allele"],
                    "other_allele": row["other_allele"],
                    "beta_exposure": row["beta"],
                    "se_exposure": row["se"],
                    "p_exposure": row.get("p", row.get("p_exposure", "")),
                    "in_mhc_region": row.get("in_mhc_region", "false"),
                    "instrument_source": input_csv.as_posix(),
                    "selection_note": f"p<={p_threshold}; exclude_mhc={exclude_mhc}; physical_window_bp={window_bp}",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Select candidate MR instruments from exposure GWAS hits.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--include-mhc", action="store_true")
    parser.add_argument("--physical-window-bp", type=int, default=1_000_000)
    args = parser.parse_args()
    select(args.input_csv, args.output_csv, args.p_threshold, not args.include_mhc, args.physical_window_bp)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
