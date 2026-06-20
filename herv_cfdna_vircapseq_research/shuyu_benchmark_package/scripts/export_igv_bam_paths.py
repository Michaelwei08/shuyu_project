from __future__ import annotations

import argparse
import csv
from pathlib import Path


COUNT_COLUMNS = ["HIV1", "HIV2", "HTLV1", "HTLV2", "HERV", "LINE1", "HUMAN"]


def read_counts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def value(row: dict[str, str], key: str) -> int:
    try:
        return int(row.get(key, "0") or 0)
    except ValueError:
        return 0


def add_rows(
    out_rows: list[dict[str, object]],
    cohort: str,
    counts_path: Path,
    work_dir: Path,
    nonzero_exogenous_only: bool,
) -> None:
    for row in read_counts(counts_path):
        sample = row["sample"]
        if nonzero_exogenous_only and not any(value(row, key) > 0 for key in ["HIV1", "HIV2", "HTLV1", "HTLV2", "HERV"]):
            continue
        bam = work_dir / "bam" / f"{sample}.retrovirus.bam"
        bai = bam.with_name(f"{bam.name}.bai")
        out_rows.append(
            {
                "cohort": cohort,
                "sample": sample,
                "bam": str(bam),
                "bai": str(bai),
                "bam_exists": "yes" if bam.exists() else "unknown_or_not_visible_here",
                "bai_exists": "yes" if bai.exists() else "unknown_or_not_visible_here",
                **{key: value(row, key) for key in COUNT_COLUMNS},
            }
        )


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["cohort", "sample", "bam", "bai", "bam_exists", "bai_exists", *COUNT_COLUMNS]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Export BAM/BAI paths with viral counts for IGV review.")
    parser.add_argument("--targeted-counts", type=Path, required=True)
    parser.add_argument("--targeted-work-dir", type=Path, required=True)
    parser.add_argument("--wgs-counts", type=Path, required=True)
    parser.add_argument("--wgs-work-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nonzero-exogenous-only", action="store_true")
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    add_rows(rows, "targeted_htlv", args.targeted_counts, args.targeted_work_dir, args.nonzero_exogenous_only)
    add_rows(rows, "wgs_hiv_hl", args.wgs_counts, args.wgs_work_dir, args.nonzero_exogenous_only)
    write_tsv(args.output, rows)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
