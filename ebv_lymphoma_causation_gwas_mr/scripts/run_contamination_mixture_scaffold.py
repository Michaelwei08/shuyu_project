from __future__ import annotations

import argparse
import csv
from pathlib import Path


def run(harmonized_csv: Path, output_csv: Path) -> None:
    ratios: list[float] = []
    with harmonized_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            exposure = float(row["beta_exposure"])
            outcome = float(row["beta_outcome"])
            if exposure != 0:
                ratios.append(outcome / exposure)
    ratios.sort()
    estimate = ratios[len(ratios) // 2] if ratios else float("nan")
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["method", "estimate", "n_instruments", "status", "note"])
        writer.writeheader()
        writer.writerow(
            {
                "method": "ConMix",
                "estimate": f"{estimate:.6f}" if ratios else "",
                "n_instruments": len(ratios),
                "status": "scaffold_only_not_validated",
                "note": "Median-ratio proxy only. Final ConMix must be rerun with a validated MR package before causal claims.",
            }
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Record a ConMix scaffold result and validation warning.")
    parser.add_argument("harmonized_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    run(args.harmonized_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
