from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def by_sample(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["sample_id"]: row for row in rows}


def coverage_by_sample(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in read_rows(path):
        result[(row["sample_id"], row["virus"].lower())] = row
    return result


def junction_by_sample(path: Path) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    for row in read_rows(path):
        result[row["sample_id"]] = {
            "hiv": int(row.get("hiv_junction_reads", 0) or 0),
            "htlv": int(row.get("htlv_junction_reads", 0) or 0),
        }
    return result


def call_virus(sample: dict[str, str], coverage: dict[str, str] | None, junctions: dict[str, int], thresholds: dict[str, object], virus: str) -> tuple[str, str]:
    metric = f"unique_{virus}_per_million"
    signal = float(sample.get(metric, 0) or 0)
    min_signal = float(thresholds["minimum_per_million"][virus])  # type: ignore[index]
    min_breadth = float(thresholds["minimum_breadth"][virus])  # type: ignore[index]
    breadth = float((coverage or {}).get("breadth", 0) or 0)
    junction_count = junctions.get(virus, 0)
    if signal <= min_signal:
        return "negative", f"{metric}={signal:.6f} <= threshold={min_signal:.6f}"
    if breadth < min_breadth and junction_count <= 0:
        return "review", f"signal passes threshold but breadth={breadth:.3f} and junction_reads={junction_count}"
    return "positive", f"signal={signal:.6f}; breadth={breadth:.3f}; junction_reads={junction_count}"


def call(manifest_csv: Path, normalized_csv: Path, coverage_csv: Path, junction_csv: Path, thresholds_json: Path, output_csv: Path) -> None:
    manifest = by_sample(read_rows(manifest_csv))
    normalized = read_rows(normalized_csv)
    coverage = coverage_by_sample(coverage_csv)
    junctions = junction_by_sample(junction_csv)
    thresholds = json.loads(thresholds_json.read_text(encoding="utf-8"))
    fields = [
        "sample_id",
        "subject_id",
        "cohort",
        "assay_type",
        "expected_group",
        "hiv_call",
        "hiv_reason",
        "htlv_call",
        "htlv_reason",
        "ambiguous_retroviral_per_million",
        "unique_herv_line1_per_million",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in normalized:
            sample_id = row["sample_id"]
            meta = manifest.get(sample_id, {})
            sample_junctions = junctions.get(sample_id, {"hiv": 0, "htlv": 0})
            hiv_call, hiv_reason = call_virus(row, coverage.get((sample_id, "hiv")), sample_junctions, thresholds, "hiv")
            htlv_call, htlv_reason = call_virus(row, coverage.get((sample_id, "htlv")), sample_junctions, thresholds, "htlv")
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "subject_id": meta.get("subject_id", ""),
                    "cohort": meta.get("cohort", ""),
                    "assay_type": meta.get("assay_type", ""),
                    "expected_group": meta.get("expected_group", "unknown"),
                    "hiv_call": hiv_call,
                    "hiv_reason": hiv_reason,
                    "htlv_call": htlv_call,
                    "htlv_reason": htlv_reason,
                    "ambiguous_retroviral_per_million": row.get("ambiguous_retroviral_per_million", ""),
                    "unique_herv_line1_per_million": row.get("unique_herv_line1_per_million", ""),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Apply frozen HIV/HTLV thresholds to normalized retroviral bins.")
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("normalized_csv", type=Path)
    parser.add_argument("coverage_csv", type=Path)
    parser.add_argument("junction_csv", type=Path)
    parser.add_argument("thresholds_json", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    call(args.manifest_csv, args.normalized_csv, args.coverage_csv, args.junction_csv, args.thresholds_json, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
