from __future__ import annotations

import argparse
import csv
from pathlib import Path

from prepare_shuyu_benchmark import MANIFEST_COLUMNS, clean_sample_id, discover_files, group_files


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_manifest(data_dir: Path, hiv_status: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if not data_dir.is_dir():
        raise FileNotFoundError(f"Targeted HIV data directory does not exist: {data_dir}")
    grouped = group_files(discover_files(data_dir))
    if not grouped:
        raise ValueError(f"No FASTQ files found under {data_dir}")

    manifest: list[dict[str, str]] = []
    issues: list[dict[str, str]] = []
    for sample_key, files in sorted(grouped.items()):
        r1 = [item for item in files if item.mate == "R1"]
        r2 = [item for item in files if item.mate == "R2"]
        if len(r1) != 1 or len(r2) != 1:
            issues.append(
                {
                    "sample_key": sample_key,
                    "r1_files": str(len(r1)),
                    "r2_files": str(len(r2)),
                    "issue": "expected_exactly_one_R1_and_one_R2",
                }
            )
            continue
        sample = f"targeted_hiv_{clean_sample_id(sample_key)}"
        row = {column: "" for column in MANIFEST_COLUMNS}
        row.update(
            {
                "sample_id": sample,
                "subject_id": clean_sample_id(sample_key),
                "cohort": "shuyu_targeted_hiv",
                "disease_group": "hiv_infection",
                "assay_type": "targeted",
                "expected_group": "hiv_targeted_expected_positive",
                "hiv_status": hiv_status,
                "htlv_status": "negative",
                "herv_analysis_scope": "family_plus_background",
                "file_format": "FASTQ",
                "read1_path": str(r1[0].path.resolve()),
                "read2_path": str(r2[0].path.resolve()),
                "reference_build": "hg38",
                "library_prep": "targeted_hiv_capture",
                "batch": "targeted_hiv",
                "notes": "Targeted HIV benchmark; confirm authoritative infection status separately.",
            }
        )
        manifest.append(row)
    return manifest, issues


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a strict paired-FASTQ manifest for targeted HIV data.")
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--issues-output", type=Path, required=True)
    parser.add_argument("--hiv-status", choices=["positive", "negative", "unknown"], default="unknown")
    args = parser.parse_args()

    manifest, issues = build_manifest(args.data_dir, args.hiv_status)
    if not manifest:
        raise ValueError("No complete R1/R2 pairs were found; inspect the issues table")
    write_csv(args.output, manifest, MANIFEST_COLUMNS)
    write_csv(args.issues_output, issues, ["sample_key", "r1_files", "r2_files", "issue"])
    print(f"Wrote {len(manifest)} complete targeted HIV samples to {args.output}")
    print(f"Wrote {len(issues)} pairing issues to {args.issues_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
