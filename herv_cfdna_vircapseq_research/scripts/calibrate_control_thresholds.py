from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


VIRUSES = {
    "hiv": "unique_hiv_per_million",
    "htlv": "unique_htlv_per_million",
}


def load_groups(manifest_csv: Path) -> dict[str, str]:
    with manifest_csv.open(newline="", encoding="utf-8") as handle:
        return {row["sample_id"]: row["expected_group"] for row in csv.DictReader(handle)}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def calibrate(manifest_csv: Path, normalized_csv: Path, thresholds_json: Path, noise_csv: Path, margin: float) -> None:
    groups = load_groups(manifest_csv)
    rows = load_rows(normalized_csv)
    controls = [row for row in rows if groups.get(row["sample_id"]) == "healthy_negative"]
    thresholds: dict[str, object] = {
        "version": "example_v1",
        "training_group": "healthy_negative",
        "margin_per_million": margin,
        "minimum_breadth": {"hiv": 0.20, "htlv": 0.20},
        "minimum_junction_reads": {"hiv": 0, "htlv": 0},
        "notes": "Dry-run threshold set. Recalibrate using real healthy controls before final analysis.",
    }
    viral_thresholds: dict[str, float] = {}
    for virus, metric in VIRUSES.items():
        control_values = [float(row.get(metric, 0) or 0) for row in controls]
        viral_thresholds[virus] = (max(control_values) if control_values else 0.0) + margin
    thresholds["minimum_per_million"] = viral_thresholds
    thresholds_json.write_text(json.dumps(thresholds, indent=2), encoding="utf-8")

    with noise_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "metric", "reads_per_million"])
        writer.writeheader()
        for virus, metric in VIRUSES.items():
            before = max([float(row.get(metric, 0) or 0) for row in controls] or [0.0])
            after = 0.0 if before <= viral_thresholds[virus] else before
            writer.writerow({"stage": "before_filtering", "metric": metric, "reads_per_million": f"{before:.6f}"})
            writer.writerow({"stage": "after_filtering", "metric": metric, "reads_per_million": f"{after:.6f}"})
        ambiguous = max([float(row.get("ambiguous_retroviral_per_million", 0) or 0) for row in controls] or [0.0])
        writer.writerow({"stage": "before_filtering", "metric": "ambiguous_retroviral_per_million", "reads_per_million": f"{ambiguous:.6f}"})
        writer.writerow({"stage": "after_filtering", "metric": "ambiguous_retroviral_per_million", "reads_per_million": f"{ambiguous:.6f}"})


def main() -> int:
    parser = argparse.ArgumentParser(description="Calibrate dry-run HIV/HTLV thresholds from healthy controls.")
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("normalized_csv", type=Path)
    parser.add_argument("thresholds_json", type=Path)
    parser.add_argument("noise_csv", type=Path)
    parser.add_argument("--margin-per-million", type=float, default=0.1)
    args = parser.parse_args()
    calibrate(args.manifest_csv, args.normalized_csv, args.thresholds_json, args.noise_csv, args.margin_per_million)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
