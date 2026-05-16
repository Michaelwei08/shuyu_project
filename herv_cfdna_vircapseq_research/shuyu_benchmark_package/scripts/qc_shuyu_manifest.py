from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


MANIFEST_REQUIRED = {
    "sample_id",
    "subject_id",
    "cohort",
    "assay_type",
    "expected_group",
    "file_format",
    "file_path",
    "read1_path",
    "read2_path",
    "bam_or_cram_path",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = MANIFEST_REQUIRED.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest missing required columns: {', '.join(sorted(missing))}")
        return list(reader)


def file_size_gb(path_text: str) -> tuple[bool, float]:
    if not path_text:
        return False, 0.0
    path = Path(path_text)
    if not path.exists():
        return False, 0.0
    return True, path.stat().st_size / 1_000_000_000


def bool_text(value: bool) -> str:
    return "yes" if value else "no"


def pct(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator * 100:.1f}"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sample_qc(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    qc_rows: list[dict[str, object]] = []
    for row in rows:
        r1_exists, r1_gb = file_size_gb(row.get("read1_path", ""))
        r2_exists, r2_gb = file_size_gb(row.get("read2_path", ""))
        single_exists, single_gb = file_size_gb(row.get("file_path", ""))
        alignment_exists, alignment_gb = file_size_gb(row.get("bam_or_cram_path", ""))
        has_any_existing_file = any([r1_exists, r2_exists, single_exists, alignment_exists])
        paired_fastq_complete = row.get("file_format") == "FASTQ" and r1_exists and r2_exists
        size_ratio = ""
        if r1_gb > 0 and r2_gb > 0:
            size_ratio = f"{max(r1_gb, r2_gb) / min(r1_gb, r2_gb):.3f}"
        issue_flags: list[str] = []
        if not has_any_existing_file:
            issue_flags.append("no_existing_sequence_file")
        if row.get("file_format") == "FASTQ" and not paired_fastq_complete:
            issue_flags.append("incomplete_fastq_pair")
        if size_ratio and float(size_ratio) > 1.25:
            issue_flags.append("r1_r2_size_imbalance")
        qc_rows.append(
            {
                "sample_id": row["sample_id"],
                "subject_id": row["subject_id"],
                "cohort": row["cohort"],
                "assay_type": row["assay_type"],
                "expected_group": row["expected_group"],
                "file_format": row["file_format"],
                "read1_exists": bool_text(r1_exists),
                "read2_exists": bool_text(r2_exists),
                "single_file_exists": bool_text(single_exists),
                "alignment_exists": bool_text(alignment_exists),
                "paired_fastq_complete": bool_text(paired_fastq_complete),
                "read1_gb": f"{r1_gb:.3f}",
                "read2_gb": f"{r2_gb:.3f}",
                "single_file_gb": f"{single_gb:.3f}",
                "alignment_gb": f"{alignment_gb:.3f}",
                "total_sequence_gb": f"{r1_gb + r2_gb + single_gb + alignment_gb:.3f}",
                "r1_r2_size_ratio": size_ratio,
                "issue_flags": ";".join(issue_flags) if issue_flags else "ok",
            }
        )
    return qc_rows


def summary_rows(qc_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in qc_rows:
        groups[(str(row["cohort"]), str(row["assay_type"]), str(row["expected_group"]))].append(row)
    output: list[dict[str, object]] = []
    for (cohort, assay_type, expected_group), rows in sorted(groups.items()):
        total = len(rows)
        paired = sum(row["paired_fastq_complete"] == "yes" for row in rows)
        missing = sum("no_existing_sequence_file" in str(row["issue_flags"]) for row in rows)
        incomplete = sum("incomplete_fastq_pair" in str(row["issue_flags"]) for row in rows)
        imbalanced = sum("r1_r2_size_imbalance" in str(row["issue_flags"]) for row in rows)
        total_gb = sum(float(row["total_sequence_gb"]) for row in rows)
        output.append(
            {
                "cohort": cohort,
                "assay_type": assay_type,
                "expected_group": expected_group,
                "samples": total,
                "paired_fastq_complete": paired,
                "paired_fastq_complete_pct": pct(paired, total),
                "missing_all_sequence_files": missing,
                "incomplete_fastq_pairs": incomplete,
                "r1_r2_size_imbalanced": imbalanced,
                "total_sequence_gb": f"{total_gb:.3f}",
            }
        )
    return output


def duplicate_subject_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    by_subject: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        by_subject[row["subject_id"]].append(row["sample_id"])
    duplicates = []
    for subject_id, sample_ids in sorted(by_subject.items()):
        if len(sample_ids) > 1:
            duplicates.append(
                {
                    "subject_id": subject_id,
                    "sample_count": len(sample_ids),
                    "sample_ids": ";".join(sample_ids),
                }
            )
    return duplicates


def write_report(
    report_path: Path,
    manifest_path: Path,
    qc_rows: list[dict[str, object]],
    summaries: list[dict[str, object]],
    duplicates: list[dict[str, object]],
) -> None:
    issue_rows = [row for row in qc_rows if row["issue_flags"] != "ok"]
    lines = [
        "# Shuyu Manifest QC",
        "",
        f"Manifest: `{manifest_path}`",
        "",
        "## Bottom Line",
        "",
        f"- Total samples: {len(qc_rows)}",
        f"- Samples with issues: {len(issue_rows)}",
        f"- Duplicate subject IDs: {len(duplicates)}",
        "",
        "## Group Summary",
        "",
        "| Cohort | Assay | Expected group | Samples | Complete FASTQ pairs | Incomplete FASTQ pairs | Total GB |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            f"| {row['cohort']} | {row['assay_type']} | {row['expected_group']} | {row['samples']} | "
            f"{row['paired_fastq_complete']} ({row['paired_fastq_complete_pct']}%) | "
            f"{row['incomplete_fastq_pairs']} | {row['total_sequence_gb']} |"
        )
    lines.extend(["", "## Samples Needing Review", ""])
    if not issue_rows:
        lines.append("- None.")
    else:
        for row in issue_rows[:50]:
            lines.append(
                f"- `{row['sample_id']}`: {row['issue_flags']}; "
                f"R1={row['read1_exists']} ({row['read1_gb']} GB), "
                f"R2={row['read2_exists']} ({row['read2_gb']} GB), "
                f"single={row['single_file_exists']} ({row['single_file_gb']} GB)."
            )
        if len(issue_rows) > 50:
            lines.append(f"- ... {len(issue_rows) - 50} additional issue rows omitted from report; see CSV.")
    lines.extend(["", "## Recommendation", ""])
    if issue_rows:
        lines.append("Resolve incomplete FASTQ pairs or confirm they are intentionally single-end before alignment.")
    else:
        lines.append("Manifest is ready for pilot alignment.")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="QC Shuyu sample manifest without reading FASTQ contents.")
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    rows = read_rows(args.manifest_csv)
    qc_rows = sample_qc(rows)
    summaries = summary_rows(qc_rows)
    duplicates = duplicate_subject_rows(rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "sample_file_qc.csv", qc_rows, list(qc_rows[0].keys()) if qc_rows else ["sample_id"])
    write_csv(
        args.output_dir / "manifest_qc_summary.csv",
        summaries,
        [
            "cohort",
            "assay_type",
            "expected_group",
            "samples",
            "paired_fastq_complete",
            "paired_fastq_complete_pct",
            "missing_all_sequence_files",
            "incomplete_fastq_pairs",
            "r1_r2_size_imbalanced",
            "total_sequence_gb",
        ],
    )
    write_csv(
        args.output_dir / "duplicate_subject_ids.csv",
        duplicates,
        ["subject_id", "sample_count", "sample_ids"],
    )
    write_report(args.output_dir / "manifest_qc_report.md", args.manifest_csv, qc_rows, summaries, duplicates)
    print(f"Wrote QC outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
