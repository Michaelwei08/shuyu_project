from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


RETRO_CATEGORIES = ["HIV1", "HIV2", "HTLV1", "HTLV2", "HERV"]


def read_table(path: Path, delimiter: str | None = None) -> list[dict[str, str]]:
    if delimiter is None:
        delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def as_int(row: dict[str, str], key: str) -> int:
    return int(float(row.get(key, "") or 0))


def pct(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator * 100:.2f}" if denominator else "NA"


def is_hiv_pretreatment(row: dict[str, str]) -> bool:
    sample = row.get("sample_id", "")
    subject = row.get("subject_id", "")
    return row.get("hiv_status") == "positive" and ("-P1" in sample or "-P1" in subject)


def detection_rows(
    counts: list[dict[str, str]],
    manifest: list[dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    counts_by_sample = {row["sample"]: row for row in counts}
    rows: list[dict[str, object]] = []
    for meta in manifest:
        if not is_hiv_pretreatment(meta):
            continue
        count = counts_by_sample.get(meta["sample_id"], {})
        hiv1 = as_int(count, "HIV1")
        hiv2 = as_int(count, "HIV2")
        htlv1 = as_int(count, "HTLV1")
        htlv2 = as_int(count, "HTLV2")
        rows.append(
            {
                "sample": meta["sample_id"],
                "subject_id": meta.get("subject_id", ""),
                "HIV1": hiv1,
                "HIV2": hiv2,
                "HIV_total": hiv1 + hiv2,
                "HTLV1": htlv1,
                "HTLV2": htlv2,
                "detected_hiv": int(hiv1 + hiv2 > 0),
            }
        )
    detected = sum(int(row["detected_hiv"]) for row in rows)
    summary = {
        "group": "HIV pretreatment WGS (-P1)",
        "samples": len(rows),
        "hiv_detected": detected,
        "hiv_detection_pct": pct(detected, len(rows)),
        "hiv1_detected": sum(int(row["HIV1"]) > 0 for row in rows),
        "hiv2_detected": sum(int(row["HIV2"]) > 0 for row in rows),
        "htlv_detected": sum(int(row["HTLV1"]) + int(row["HTLV2"]) > 0 for row in rows),
    }
    return rows, summary


def reference_rows(reference_map: list[dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in reference_map:
        category = row.get("category", "")
        if category not in RETRO_CATEGORIES:
            continue
        rows.append(
            {
                "category": category,
                "reference_id": row.get("reference_id", ""),
                "source": row.get("source", ""),
                "description": row.get("description", ""),
            }
        )
    return sorted(rows, key=lambda row: (str(row["category"]), str(row["reference_id"])))


def mask_rows(
    mask_summary: list[dict[str, str]],
    pair_summary: list[dict[str, str]],
) -> list[dict[str, object]]:
    post_pct_by_ref: dict[str, float] = defaultdict(float)
    kmer_by_ref: dict[str, str] = {}
    against_by_ref: dict[str, set[str]] = defaultdict(set)
    for row in pair_summary:
        ref_id = row.get("query_reference_id", "")
        if row.get("against_category") != "HERV":
            continue
        post_pct_by_ref[ref_id] = max(post_pct_by_ref[ref_id], float(row.get("post_mask_query_similar_pct") or 0))
        kmer_by_ref[ref_id] = row.get("kmer_size", "")
        against_by_ref[ref_id].add(row.get("against_reference_id", ""))

    rows: list[dict[str, object]] = []
    for row in mask_summary:
        category = row.get("category", "")
        if category not in {"HIV1", "HTLV1"}:
            continue
        ref_id = row["reference_id"]
        before_pct = row.get("masked_pct", "0")
        rows.append(
            {
                "query_category": category,
                "query_reference_id": ref_id,
                "against_category": "HERV",
                "against_reference_ids": ",".join(sorted(against_by_ref.get(ref_id, set()))),
                "kmer_size": kmer_by_ref.get(ref_id, ""),
                "query_length_bp": row.get("length", ""),
                "before_mask_similar_bp": row.get("masked_bases", ""),
                "before_mask_similar_pct": before_pct,
                "after_mask_similar_pct": f"{post_pct_by_ref.get(ref_id, 0.0):.6f}",
                "retained_after_mask_pct": row.get("retained_pct", ""),
            }
        )
    return sorted(rows, key=lambda row: str(row["query_category"]))


def optional_masked_comparison(
    unmasked: list[dict[str, str]],
    masked: list[dict[str, str]] | None,
) -> list[dict[str, object]]:
    if masked is None:
        return []
    masked_by_sample = {row["sample"]: row for row in masked}
    rows: list[dict[str, object]] = []
    for row in unmasked:
        sample = row["sample"]
        masked_row = masked_by_sample.get(sample, {})
        rows.append(
            {
                "sample": sample,
                "unmasked_HIV1": as_int(row, "HIV1"),
                "masked_HIV1": as_int(masked_row, "HIV1"),
                "unmasked_HTLV1": as_int(row, "HTLV1"),
                "masked_HTLV1": as_int(masked_row, "HTLV1"),
                "unmasked_HERV": as_int(row, "HERV"),
                "masked_HERV": as_int(masked_row, "HERV"),
            }
        )
    return rows


def write_report(
    path: Path,
    detection_summary: dict[str, object],
    detected_rows: list[dict[str, object]],
    refs: list[dict[str, object]],
    masks: list[dict[str, object]],
    has_masked_counts: bool,
) -> None:
    detected_samples = [str(row["sample"]) for row in detected_rows if int(row["detected_hiv"])]
    lines = [
        "# Shuyu Reply Metrics",
        "",
        "## HIV pretreatment detection",
        "",
        f"- Denominator: {detection_summary['samples']} HIV-positive WGS pretreatment samples labeled `-P1`.",
        f"- HIV detected: {detection_summary['hiv_detected']}/{detection_summary['samples']} ({detection_summary['hiv_detection_pct']}%).",
        f"- HIV1 detected: {detection_summary['hiv1_detected']}; HIV2 detected: {detection_summary['hiv2_detected']}.",
        f"- HTLV detected in these HIV pretreatment samples: {detection_summary['htlv_detected']}.",
        f"- Detected sample IDs: {', '.join(detected_samples) if detected_samples else 'none'}.",
        "",
        "## Viral reference IDs",
        "",
        "| category | reference_id | source |",
        "|---|---|---|",
    ]
    for row in refs:
        lines.append(f"| {row['category']} | {row['reference_id']} | {row['source']} |")

    lines.extend(
        [
            "",
            "## HERV/HIV/HTLV masking status",
            "",
            "- Current primary-only counts are from competitive hg38 + retroviral alignment with secondary/supplementary reads excluded.",
        ]
    )
    if has_masked_counts:
        lines.append("- A masked-reference count table was also provided and compared in `masked_vs_unmasked_counts.tsv`.")
    else:
        lines.append("- No masked-reference count table was provided to this summary; masking below is a reference-space mask calculation, not a masked-result rerun.")

    lines.extend(
        [
            "- Similarity below is exact shared-kmer similarity against HERV used by `mask_shared_retro_regions.py`; it is not the VirCAPP production mask unless Shuyu's provided masked panel is used directly.",
            "",
            "| query | reference | before mask similar % | after mask similar % | retained % |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in masks:
        lines.append(
            f"| {row['query_category']} | {row['query_reference_id']} | "
            f"{row['before_mask_similar_pct']} | {row['after_mask_similar_pct']} | "
            f"{row['retained_after_mask_pct']} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create Shuyu reply metrics from primary-only P1 outputs.")
    parser.add_argument("--wgs-counts", type=Path, required=True)
    parser.add_argument("--wgs-manifest", type=Path, required=True)
    parser.add_argument("--reference-map", type=Path, required=True)
    parser.add_argument("--mask-summary", type=Path, required=True)
    parser.add_argument("--pair-summary", type=Path, required=True)
    parser.add_argument("--masked-wgs-counts", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    wgs_counts = read_table(args.wgs_counts, "\t")
    wgs_manifest = read_table(args.wgs_manifest, ",")
    refs = reference_rows(read_table(args.reference_map, ","))
    masks = mask_rows(read_table(args.mask_summary, "\t"), read_table(args.pair_summary, "\t"))
    masked_counts = read_table(args.masked_wgs_counts, "\t") if args.masked_wgs_counts else None

    p1_rows, p1_summary = detection_rows(wgs_counts, wgs_manifest)
    comparison = optional_masked_comparison(wgs_counts, masked_counts)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "hiv_p1_detection.tsv", p1_rows, list(p1_rows[0].keys()) if p1_rows else ["sample"])
    write_tsv(args.output_dir / "hiv_p1_detection_summary.tsv", [p1_summary], list(p1_summary.keys()))
    write_tsv(args.output_dir / "viral_reference_ids.tsv", refs, ["category", "reference_id", "source", "description"])
    write_tsv(args.output_dir / "mask_similarity_summary.tsv", masks, list(masks[0].keys()) if masks else ["query_category"])
    if comparison:
        write_tsv(args.output_dir / "masked_vs_unmasked_counts.tsv", comparison, list(comparison[0].keys()))
    write_report(args.output_dir / "shuyu_reply_metrics.md", p1_summary, p1_rows, refs, masks, bool(comparison))
    print(f"Wrote Shuyu reply metrics to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
