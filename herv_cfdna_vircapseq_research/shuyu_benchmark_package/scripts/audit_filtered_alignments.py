from __future__ import annotations

import argparse
import csv
import hashlib
import re
import subprocess
from pathlib import Path

from retro_alignment_filters import dedup_key, sam_tag


def read_table(path: Path, delimiter: str | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        if delimiter is None:
            delimiter = "\t" if path.suffix == ".tsv" else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def aligned_len(cigar: str) -> int:
    total = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in "M=X":
            total += int(n)
    return total


def tag_value(fields: list[str], tag: str) -> str:
    prefix = f"{tag}:"
    for field in fields[11:]:
        if field.startswith(prefix):
            return field.split(":", 2)[-1]
    return "NA"


def read_hash(read_id: str) -> str:
    return hashlib.sha256(read_id.encode("utf-8")).hexdigest()[:16]


def load_reference_map(path: Path) -> dict[str, str]:
    rows = read_table(path, ",")
    return {row["reference_id"]: row["category"] for row in rows}


def selected_samples(rows: list[dict[str, str]], categories: set[str], sample_ids: list[str]) -> list[str]:
    if sample_ids:
        return sample_ids
    selected: list[str] = []
    for row in rows:
        if any(int(float(row.get(category, "") or 0)) > 0 for category in categories):
            selected.append(row["sample"])
    return selected


def sam_records(
    samtools_exe: str,
    bam_path: Path,
    reference_ids: list[str],
) -> list[list[str]]:
    cmd = [samtools_exe, "view", "-F", "4", str(bam_path), *reference_ids]
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, check=True)
    return [line.split("\t") for line in result.stdout.splitlines() if line]


def audit_sample(
    sample: str,
    bam_dir: Path,
    reference_map: dict[str, str],
    categories: set[str],
    min_mapq: int,
    min_aligned_length: int,
    samtools_exe: str,
    include_read_id: bool,
    dedup_mode: str,
    require_unique_best: bool,
) -> list[dict[str, object]]:
    bam_path = bam_dir / f"{sample}.retrovirus.bam"
    if not bam_path.exists():
        raise FileNotFoundError(f"Missing BAM for {sample}: {bam_path}")
    refs = [ref for ref, category in reference_map.items() if category in categories]
    grouped: dict[tuple[str, str], list[list[str]]] = {}
    seen: set[tuple[str, ...]] = set()
    for fields in sam_records(samtools_exe, bam_path, refs):
        if len(fields) < 11:
            continue
        read_id, ref, mapq_text, cigar = fields[0], fields[2], fields[4], fields[5]
        category = reference_map.get(ref)
        if category not in categories:
            continue
        alen = aligned_len(cigar)
        if int(mapq_text) < min_mapq or alen < min_aligned_length:
            continue
        if require_unique_best:
            alignment_score = sam_tag(fields, "AS")
            suboptimal_score = sam_tag(fields, "XS")
            if alignment_score is not None and suboptimal_score is not None and alignment_score <= suboptimal_score:
                continue
        key = dedup_key(dedup_mode, fields, category, alen, None)
        if key is not None:
            key_text = tuple(str(item) for item in key)
            if key_text in seen:
                continue
            seen.add(key_text)
        grouped.setdefault((read_id, category), []).append(fields)

    rows: list[dict[str, object]] = []
    for (read_id, category), records in sorted(grouped.items(), key=lambda item: (item[0][1], item[0][0])):
        refs_seen = [record[2] for record in records]
        rows.append(
            {
                "sample": sample,
                "category": category,
                "read_id_hash": read_hash(read_id),
                "read_id": read_id if include_read_id else "",
                "record_count": len(records),
                "references": ",".join(refs_seen),
                "positions": ",".join(record[3] for record in records),
                "mapqs": ",".join(record[4] for record in records),
                "cigars": ",".join(record[5] for record in records),
                "aligned_lengths": ",".join(str(aligned_len(record[5])) for record in records),
                "NM": ",".join(tag_value(record, "NM") for record in records),
                "AS": ",".join(tag_value(record, "AS") for record in records),
                "XS": ",".join(tag_value(record, "XS") for record in records),
                "SA_present": any(tag_value(record, "SA") != "NA" for record in records),
                "XA_present": any(tag_value(record, "XA") != "NA" for record in records),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit filtered viral alignments without outputting sequences or qualities.")
    parser.add_argument("--sample-summary", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--reference-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--category", action="append")
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-aligned-length", type=int, default=60)
    parser.add_argument("--samtools-exe", default="samtools")
    parser.add_argument("--include-read-id", action="store_true")
    parser.add_argument("--dedup-mode", choices=["none", "read_id", "coordinate", "fragment"], default="none")
    parser.add_argument("--require-unique-best", action="store_true")
    args = parser.parse_args()

    categories = set(args.category or ["HIV1", "HIV2"])
    summary_rows = read_table(args.sample_summary, "\t")
    reference_map = load_reference_map(args.reference_map)
    samples = selected_samples(summary_rows, categories, args.sample_id)

    rows: list[dict[str, object]] = []
    for sample in samples:
        rows.extend(
            audit_sample(
                sample,
                args.bam_dir,
                reference_map,
                categories,
                args.min_mapq,
                args.min_aligned_length,
                args.samtools_exe,
                args.include_read_id,
                args.dedup_mode,
                args.require_unique_best,
            )
        )
    write_tsv(
        args.output,
        rows,
        [
            "sample",
            "category",
            "read_id_hash",
            "read_id",
            "record_count",
            "references",
            "positions",
            "mapqs",
            "cigars",
            "aligned_lengths",
            "NM",
            "AS",
            "XS",
            "SA_present",
            "XA_present",
        ],
    )
    print(f"Wrote {len(rows)} audited fragments to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
