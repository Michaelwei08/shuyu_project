from __future__ import annotations

import argparse
import csv
from pathlib import Path

from run_mr_summary import run


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_ivw_estimate(path: Path) -> str:
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["method"] == "IVW":
                return row["estimate"]
    return ""


def sensitivity(harmonized_csv: Path, output_csv: Path, temp_dir: Path) -> None:
    rows = read_rows(harmonized_csv)
    if len(rows) < 3:
        raise ValueError("Leave-one-out sensitivity needs at least three instruments for a useful scaffold run.")
    temp_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    out_fields = ["dropped_variant_id", "n_instruments", "ivw_estimate", "note"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        for idx, row in enumerate(rows):
            retained = [candidate for pos, candidate in enumerate(rows) if pos != idx]
            temp_input = temp_dir / f"loo_{idx}.csv"
            temp_output = temp_dir / f"loo_{idx}.mr.csv"
            write_rows(temp_input, retained, fieldnames)
            run(temp_input, temp_output)
            writer.writerow(
                {
                    "dropped_variant_id": row.get("variant_id", f"row_{idx}"),
                    "n_instruments": len(retained),
                    "ivw_estimate": load_ivw_estimate(temp_output),
                    "note": "leave-one-instrument-out IVW scaffold",
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run lightweight leave-one-out MR sensitivity summaries.")
    parser.add_argument("harmonized_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--temp-dir", type=Path, default=Path("results/mr_sensitivity_tmp"))
    args = parser.parse_args()
    sensitivity(args.harmonized_csv, args.output_csv, args.temp_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
