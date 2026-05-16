from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
ALIGNMENT_EXTENSIONS = (".bam", ".cram")
KNOWN_EXTENSIONS = FASTQ_EXTENSIONS + ALIGNMENT_EXTENSIONS

MANIFEST_COLUMNS = [
    "sample_id",
    "subject_id",
    "cohort",
    "disease_group",
    "assay_type",
    "expected_group",
    "hiv_status",
    "htlv_status",
    "herv_analysis_scope",
    "file_format",
    "file_path",
    "read1_path",
    "read2_path",
    "bam_or_cram_path",
    "reference_build",
    "estimated_depth",
    "library_prep",
    "batch",
    "notes",
]


@dataclass(frozen=True)
class SequencingFile:
    path: Path
    sample_key: str
    mate: str
    file_format: str


def strip_known_extension(name: str) -> str:
    lowered = name.lower()
    for extension in KNOWN_EXTENSIONS:
        if lowered.endswith(extension):
            return name[: -len(extension)]
    return Path(name).stem


def sample_key_from_name(name: str) -> tuple[str, str]:
    stem = strip_known_extension(name)
    mate = ""
    patterns = [
        (r"(.+?)(?:[_\.-])R?1(?:[_\.-].*)?$", "R1"),
        (r"(.+?)(?:[_\.-])R?2(?:[_\.-].*)?$", "R2"),
    ]
    for pattern, value in patterns:
        match = re.match(pattern, stem, flags=re.IGNORECASE)
        if match:
            return clean_sample_id(match.group(1)), value
    return clean_sample_id(stem), mate


def clean_sample_id(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "unknown_sample"


def infer_file_format(path: Path) -> str:
    lowered = path.name.lower()
    if lowered.endswith(FASTQ_EXTENSIONS):
        return "FASTQ"
    if lowered.endswith(".bam"):
        return "BAM"
    if lowered.endswith(".cram"):
        return "CRAM"
    return "unknown"


def discover_files(root: Path) -> list[SequencingFile]:
    if not root.exists():
        return []
    files: list[SequencingFile] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if not path.name.lower().endswith(KNOWN_EXTENSIONS):
            continue
        sample_key, mate = sample_key_from_name(path.name)
        files.append(SequencingFile(path=path, sample_key=sample_key, mate=mate, file_format=infer_file_format(path)))
    return files


def group_files(files: list[SequencingFile]) -> dict[str, list[SequencingFile]]:
    grouped: dict[str, list[SequencingFile]] = {}
    for item in files:
        grouped.setdefault(item.sample_key, []).append(item)
    return grouped


def infer_wgs_metadata(sample_id: str) -> dict[str, str]:
    upper = sample_id.upper()
    if upper.startswith("HIV"):
        return {
            "cohort": "shuyu_wgs_hiv_dlbcl_bl",
            "disease_group": "dlbcl_bl_hiv_infection",
            "expected_group": "hiv_wgs_expected_low_positive",
            "hiv_status": "positive",
            "htlv_status": "negative",
            "notes": "Filename starts with HIV; Shuyu says these are HIV-infected patients and may contain a few HIV reads.",
        }
    if upper.startswith("HL"):
        return {
            "cohort": "shuyu_wgs_hl_no_hiv",
            "disease_group": "hl_without_hiv_infection",
            "expected_group": "hl_wgs_hiv_htlv_negative_control",
            "hiv_status": "negative",
            "htlv_status": "negative",
            "notes": "Filename starts with HL; Shuyu says these are without HIV infection. HIV/HTLV signal should be treated as likely false positive or HERV/retroelement-related until proven otherwise.",
        }
    return {
        "cohort": "shuyu_wgs_unknown_prefix",
        "disease_group": "unknown",
        "expected_group": "unknown",
        "hiv_status": "unknown",
        "htlv_status": "unknown",
        "notes": "WGS file did not start with HIV or HL; confirm label with Shuyu.",
    }


def infer_htlv_targeted_metadata(sample_id: str) -> dict[str, str]:
    return {
        "cohort": "shuyu_targeted_htlv_tcl",
        "disease_group": "tcl_htlv_enriched",
        "expected_group": "htlv_targeted_enriched_expected_many_positive",
        "hiv_status": "negative",
        "htlv_status": "unknown",
        "notes": "Targeted HTLV TCL data. Shuyu says many samples have HTLV infection and targeted sequencing should enrich HTLV reads, but not necessarily every sample is positive.",
    }


def choose_paths(items: list[SequencingFile]) -> dict[str, str]:
    read1 = next((str(item.path) for item in items if item.mate == "R1"), "")
    read2 = next((str(item.path) for item in items if item.mate == "R2"), "")
    alignment = next((str(item.path) for item in items if item.file_format in {"BAM", "CRAM"}), "")
    single = ""
    if not read1 and not alignment and items:
        single = str(items[0].path)
    file_format = "unknown"
    if alignment:
        file_format = "CRAM" if alignment.lower().endswith(".cram") else "BAM"
    elif read1 or read2 or single:
        file_format = "FASTQ"
    return {
        "file_format": file_format,
        "file_path": single,
        "read1_path": read1,
        "read2_path": read2,
        "bam_or_cram_path": alignment,
    }


def build_rows(
    wgs_dir: Path,
    htlv_targeted_dir: Path,
    reference_build: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    access_rows: list[dict[str, str]] = []
    sources = [
        ("wgs_60samples_hiv_hl", wgs_dir, "WGS"),
        ("targeted_htlv_tcl", htlv_targeted_dir, "targeted"),
    ]
    for source_name, root, assay_type in sources:
        files = discover_files(root)
        access_rows.append(
            {
                "source_name": source_name,
                "path": str(root),
                "path_exists": "yes" if root.exists() else "no",
                "discovered_sequence_files": str(len(files)),
            }
        )
        for sample_key, items in sorted(group_files(files).items()):
            metadata = infer_wgs_metadata(sample_key) if assay_type == "WGS" else infer_htlv_targeted_metadata(sample_key)
            paths = choose_paths(items)
            sample_id = f"{source_name}_{sample_key}"
            rows.append(
                {
                    "sample_id": sample_id,
                    "subject_id": sample_key,
                    "cohort": metadata["cohort"],
                    "disease_group": metadata["disease_group"],
                    "assay_type": assay_type,
                    "expected_group": metadata["expected_group"],
                    "hiv_status": metadata["hiv_status"],
                    "htlv_status": metadata["htlv_status"],
                    "herv_analysis_scope": "family_plus_background",
                    "file_format": paths["file_format"],
                    "file_path": paths["file_path"],
                    "read1_path": paths["read1_path"],
                    "read2_path": paths["read2_path"],
                    "bam_or_cram_path": paths["bam_or_cram_path"],
                    "reference_build": reference_build,
                    "estimated_depth": "",
                    "library_prep": "WGS" if assay_type == "WGS" else "targeted_htlv_capture",
                    "batch": source_name,
                    "notes": metadata["notes"],
                }
            )
    return rows, access_rows


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def select_pilot(rows: list[dict[str, str]], hiv_wgs: int, hl_wgs: int, htlv_targeted: int) -> list[dict[str, str]]:
    buckets = [
        ("hiv_wgs", hiv_wgs, lambda row: row["assay_type"] == "WGS" and row["hiv_status"] == "positive"),
        ("hl_wgs_control", hl_wgs, lambda row: row["assay_type"] == "WGS" and row["expected_group"].startswith("hl_wgs")),
        ("htlv_targeted", htlv_targeted, lambda row: row["assay_type"] == "targeted"),
    ]
    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    for _name, limit, predicate in buckets:
        for row in [candidate for candidate in rows if predicate(candidate)][:limit]:
            if row["sample_id"] not in seen:
                selected.append(row)
                seen.add(row["sample_id"])
    return selected


def write_plan(output_dir: Path, rows: list[dict[str, str]], pilot: list[dict[str, str]], access_rows: list[dict[str, str]]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        key = f"{row['assay_type']}::{row['expected_group']}"
        counts[key] = counts.get(key, 0) + 1
    lines = [
        "# Shuyu Benchmark Plan",
        "",
        "## Data Access Check",
        "",
    ]
    for row in access_rows:
        lines.append(f"- `{row['source_name']}`: path `{row['path']}`, exists={row['path_exists']}, files={row['discovered_sequence_files']}.")
    lines.extend(["", "## Manifest Summary", ""])
    for key, value in sorted(counts.items()):
        lines.append(f"- `{key}`: {value} samples")
    lines.extend(
        [
            "",
            "## Pilot Samples",
            "",
            "Run the pilot first before scaling to all samples.",
        ]
    )
    for row in pilot:
        lines.append(f"- `{row['sample_id']}`: {row['expected_group']}, format={row['file_format']}")
    lines.extend(
        [
            "",
            "## Benchmark Logic",
            "",
            "- HIV-prefixed WGS samples: expected HIV-infected group; a few HIV reads may be detectable.",
            "- HL-prefixed WGS samples: HIV/HTLV-negative controls; HIV or HTLV calls are likely false positives, cross-mapping, or HERV/retroelement signal.",
            "- Targeted HTLV TCL samples: expected HTLV-enriched group; many but not all samples should be HTLV-positive.",
            "",
            "## Next Required Output From Alignment",
            "",
            "The existing P1 scripts expect a competitive hit table with columns:",
            "",
            "`sample_id,read_id,reference_category,mapq,alignment_score,alignment_length`",
            "",
            "Allowed `reference_category` values should include `human`, `hiv`, `htlv`, `herv`, `line1`, and `other_viral`.",
        ]
    )
    (output_dir / "benchmark_plan.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_command_templates(output_dir: Path) -> None:
    commands = """#!/usr/bin/env bash
set -euo pipefail

ROOT="herv_cfdna_vircapseq_research"
OUT="$ROOT/shuyu_benchmark_package/output"

python "$ROOT/scripts/validate_sample_manifest.py" "$OUT/sample_manifest.csv"

# After competitive alignment is complete, write:
#   $OUT/competitive_hits.csv
#   $OUT/depths.csv
#   $OUT/viral_coverage.csv
#   $OUT/junction_reads.csv

python "$ROOT/scripts/bin_retroviral_hits.py" "$OUT/competitive_hits.csv" "$OUT/retroviral_bins.csv"
python "$ROOT/scripts/normalize_retroviral_bins.py" "$OUT/retroviral_bins.csv" "$OUT/retroviral_bins_normalized.csv" --depth-csv "$OUT/depths.csv"
python "$ROOT/scripts/summarize_viral_coverage.py" "$OUT/viral_coverage.csv" "$OUT/viral_coverage_summary.csv"
python "$ROOT/scripts/summarize_junction_reads.py" "$OUT/junction_reads.csv" "$OUT/junction_summary.csv"
python "$ROOT/scripts/calibrate_control_thresholds.py" "$OUT/sample_manifest.csv" "$OUT/retroviral_bins_normalized.csv" "$OUT/frozen_thresholds.json" "$OUT/control_noise_before_after.csv"
python "$ROOT/scripts/call_retroviral_status.py" "$OUT/sample_manifest.csv" "$OUT/retroviral_bins_normalized.csv" "$OUT/viral_coverage_summary.csv" "$OUT/junction_summary.csv" "$OUT/frozen_thresholds.json" "$OUT/frozen_calls.csv"
python "$ROOT/scripts/inspect_control_false_signal.py" "$OUT/sample_manifest.csv" "$OUT/competitive_hits.csv" "$OUT/control_false_signal.csv"
python "$ROOT/scripts/summarize_positive_sensitivity.py" "$OUT/frozen_calls.csv" "$OUT/sensitivity_summary.csv"
python "$ROOT/scripts/inspect_low_signal_positives.py" "$OUT/sample_manifest.csv" "$OUT/frozen_calls.csv" "$OUT/competitive_hits.csv" "$OUT/low_signal_positive_review.csv"
python "$ROOT/scripts/compare_assays_and_wgs.py" "$OUT/frozen_calls.csv" "$OUT/assay_comparison.csv"
"""
    (output_dir / "run_after_alignment.sh").write_text(commands, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Shuyu HIV/HL WGS and HTLV targeted benchmark manifests.")
    parser.add_argument("--wgs-dir", type=Path, default=Path("/drive3/Shuyu_for_Michael/WGS_60samples_HIV_HL"))
    parser.add_argument("--htlv-targeted-dir", type=Path, default=Path("/drive3/Shuyu_for_Michael/Targeted_HTLV_TCLsamples"))
    parser.add_argument("--output-dir", type=Path, default=Path("herv_cfdna_vircapseq_research/shuyu_benchmark_package/output"))
    parser.add_argument("--reference-build", default="unknown")
    parser.add_argument("--pilot-hiv-wgs", type=int, default=2)
    parser.add_argument("--pilot-hl-wgs", type=int, default=2)
    parser.add_argument("--pilot-htlv-targeted", type=int, default=3)
    args = parser.parse_args()

    rows, access_rows = build_rows(args.wgs_dir, args.htlv_targeted_dir, args.reference_build)
    pilot = select_pilot(rows, args.pilot_hiv_wgs, args.pilot_hl_wgs, args.pilot_htlv_targeted)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "sample_manifest.csv", rows, MANIFEST_COLUMNS)
    write_csv(args.output_dir / "pilot_manifest.csv", pilot, MANIFEST_COLUMNS)
    write_csv(args.output_dir / "data_access_check.csv", access_rows, ["source_name", "path", "path_exists", "discovered_sequence_files"])
    write_plan(args.output_dir, rows, pilot, access_rows)
    write_command_templates(args.output_dir)

    print(f"Wrote {len(rows)} manifest rows to {args.output_dir / 'sample_manifest.csv'}")
    print(f"Wrote {len(pilot)} pilot rows to {args.output_dir / 'pilot_manifest.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
