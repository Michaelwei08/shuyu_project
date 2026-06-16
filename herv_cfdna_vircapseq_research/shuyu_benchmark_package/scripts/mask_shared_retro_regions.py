from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


DNA = set("ACGT")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_header_id(header: str) -> str:
    return header[1:].strip().split()[0]


def read_reference_map(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "reference_id" not in rows[0] or "category" not in rows[0]:
        raise ValueError("reference map must contain reference_id and category columns")
    return {row["reference_id"]: row["category"] for row in rows}


def iter_fasta(path: Path):
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks).upper()


def write_record(handle, header: str, sequence: str, width: int = 80) -> None:
    handle.write(header + "\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def revcomp(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def valid_kmer(kmer: str) -> bool:
    return set(kmer) <= DNA


def category_ref_ids(reference_map: dict[str, str], categories: set[str]) -> set[str]:
    refs = {ref for ref, category in reference_map.items() if category in categories}
    missing = sorted(categories - set(reference_map.values()))
    if missing:
        raise ValueError(f"categories absent from reference map: {', '.join(missing)}")
    if not refs:
        raise ValueError(f"no references found for categories: {', '.join(sorted(categories))}")
    return refs


def load_sequences(path: Path, ref_ids: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for header, sequence in iter_fasta(path):
        ref_id = parse_header_id(header)
        if ref_id in ref_ids:
            sequences[ref_id] = sequence
    missing = sorted(ref_ids - set(sequences))
    if missing:
        raise ValueError(f"reference IDs absent from FASTA: {', '.join(missing[:10])}")
    return sequences


def build_kmer_index(sequence: str, k: int) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for start in range(0, len(sequence) - k + 1):
        kmer = sequence[start : start + k]
        if not valid_kmer(kmer):
            continue
        index.setdefault(kmer, []).append(start)
        rc = revcomp(kmer)
        if rc != kmer:
            index.setdefault(rc, []).append(start)
    return index


def mark(mask: bytearray, start: int, end: int) -> None:
    for idx in range(max(start, 0), min(end, len(mask))):
        mask[idx] = 1


def mask_intervals(mask: bytearray) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(mask):
        if value and start is None:
            start = idx
        elif not value and start is not None:
            intervals.append((start, idx))
            start = None
    if start is not None:
        intervals.append((start, len(mask)))
    return intervals


def pct(count: int, total: int) -> str:
    return f"{count / total * 100:.6f}" if total else "0.000000"


def analyze_pair(query_seq: str, against_seq: str, k: int) -> tuple[bytearray, bytearray, int]:
    query_mask = bytearray(len(query_seq))
    against_mask = bytearray(len(against_seq))
    index = build_kmer_index(against_seq, k)
    shared = 0
    for start in range(0, len(query_seq) - k + 1):
        kmer = query_seq[start : start + k]
        if not valid_kmer(kmer):
            continue
        hits = index.get(kmer)
        if not hits:
            continue
        shared += 1
        mark(query_mask, start, start + k)
        for hit in hits:
            mark(against_mask, hit, hit + k)
    return query_mask, against_mask, shared


def apply_masks(input_fasta: Path, output_fasta: Path, masks: dict[str, bytearray]) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as out:
        for header, sequence in iter_fasta(input_fasta):
            ref_id = parse_header_id(header)
            mask = masks.get(ref_id)
            if mask:
                sequence = masked_sequence(sequence, mask)
            write_record(out, header, sequence)


def masked_sequence(sequence: str, mask: bytearray) -> str:
    chars = list(sequence)
    for start, end in mask_intervals(mask):
        chars[start:end] = "N" * (end - start)
    return "".join(chars)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_bed(path: Path, masks: dict[str, bytearray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for ref_id in sorted(masks):
            for start, end in mask_intervals(masks[ref_id]):
                handle.write(f"{ref_id}\t{start}\t{end}\tshared_kmer_mask\n")


def write_report(
    path: Path,
    k: int,
    mask_categories: list[str],
    against_categories: list[str],
    summary_rows: list[dict[str, object]],
) -> None:
    lines = [
        "# Shared Retroelement Mask Report",
        "",
        f"- K-mer size: {k}",
        f"- Masked categories: {', '.join(mask_categories)}",
        f"- Compared against categories: {', '.join(against_categories)}",
        "",
        "| reference | category | length | masked bp | masked % | retained % |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['reference_id']} | {row['category']} | {row['length']} | "
            f"{row['masked_bases']} | {row['masked_pct']} | {row['retained_pct']} |"
        )
    lines.extend(
        [
            "",
            "Pairwise before/after exact-kmer similarity metrics are written to the pair-summary TSV.",
            "",
            "This is an exact shared-kmer mask. If VirCAPP-seq uses a different production masking rule, use that rule for the final benchmark.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Mask HIV/HTLV regions sharing exact k-mers with HERV references.")
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--reference-map", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--mask-bed", type=Path, required=True)
    parser.add_argument("--summary-tsv", type=Path, required=True)
    parser.add_argument("--pair-summary-tsv", type=Path, required=True)
    parser.add_argument("--report-md", type=Path, required=True)
    parser.add_argument("--output-map", type=Path, help="Optional copy of the reference map next to the masked FASTA.")
    parser.add_argument("--mask-category", action="append", required=True)
    parser.add_argument("--against-category", action="append", required=True)
    parser.add_argument("--kmer-size", type=int, default=40)
    args = parser.parse_args()
    if args.kmer_size < 15:
        parser.error("--kmer-size should be >= 15 for a specific retroviral mask")

    reference_map = read_reference_map(args.reference_map)
    mask_categories = set(args.mask_category)
    against_categories = set(args.against_category)
    mask_refs = category_ref_ids(reference_map, mask_categories)
    against_refs = category_ref_ids(reference_map, against_categories)
    mask_sequences = load_sequences(args.input_fasta, mask_refs)
    against_sequences = load_sequences(args.input_fasta, against_refs)

    aggregate_masks = {ref_id: bytearray(len(seq)) for ref_id, seq in mask_sequences.items()}
    pair_rows: list[dict[str, object]] = []
    for query_ref, query_seq in mask_sequences.items():
        for against_ref, against_seq in against_sequences.items():
            query_mask, against_mask, shared = analyze_pair(query_seq, against_seq, args.kmer_size)
            for start, end in mask_intervals(query_mask):
                mark(aggregate_masks[query_ref], start, end)
            pair_rows.append(
                {
                    "query_reference_id": query_ref,
                    "query_category": reference_map[query_ref],
                    "against_reference_id": against_ref,
                    "against_category": reference_map[against_ref],
                    "kmer_size": args.kmer_size,
                    "shared_query_kmers": shared,
                    "query_length": len(query_seq),
                    "query_similar_bases": sum(query_mask),
                    "query_similar_pct": pct(sum(query_mask), len(query_seq)),
                    "against_length": len(against_seq),
                    "against_similar_bases": sum(against_mask),
                    "against_similar_pct": pct(sum(against_mask), len(against_seq)),
                    "post_mask_shared_query_kmers": 0,
                    "post_mask_query_similar_bases": 0,
                    "post_mask_query_similar_pct": "0.000000",
                }
            )

    for row in pair_rows:
        query_ref = str(row["query_reference_id"])
        against_ref = str(row["against_reference_id"])
        masked_query = masked_sequence(mask_sequences[query_ref], aggregate_masks[query_ref])
        post_query_mask, _, post_shared = analyze_pair(masked_query, against_sequences[against_ref], args.kmer_size)
        post_similar = sum(post_query_mask)
        row["post_mask_shared_query_kmers"] = post_shared
        row["post_mask_query_similar_bases"] = post_similar
        row["post_mask_query_similar_pct"] = pct(post_similar, len(masked_query))

    summary_rows: list[dict[str, object]] = []
    for ref_id, sequence in sorted(mask_sequences.items()):
        masked = sum(aggregate_masks[ref_id])
        summary_rows.append(
            {
                "reference_id": ref_id,
                "category": reference_map[ref_id],
                "length": len(sequence),
                "masked_bases": masked,
                "masked_pct": pct(masked, len(sequence)),
                "retained_bases": len(sequence) - masked,
                "retained_pct": pct(len(sequence) - masked, len(sequence)),
            }
        )

    apply_masks(args.input_fasta, args.output_fasta, aggregate_masks)
    write_bed(args.mask_bed, aggregate_masks)
    write_tsv(args.summary_tsv, summary_rows, list(summary_rows[0].keys()) if summary_rows else [])
    write_tsv(args.pair_summary_tsv, pair_rows, list(pair_rows[0].keys()) if pair_rows else [])
    write_report(args.report_md, args.kmer_size, args.mask_category, args.against_category, summary_rows)
    if args.output_map:
        args.output_map.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.reference_map, args.output_map)
    print(f"Wrote masked FASTA to {args.output_fasta}")
    print(f"Wrote mask BED to {args.mask_bed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
