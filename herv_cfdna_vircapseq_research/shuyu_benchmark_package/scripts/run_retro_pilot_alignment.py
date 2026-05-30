from __future__ import annotations

import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from retro_alignment_filters import filtered_counts, parse_category_min_mapq

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

def load_reference_map(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    mapping: dict[str, str] = {}
    for row in rows:
        mapping[row["reference_id"]] = row["category"]
    return mapping

def run(cmd: list[str], log_path: Path | None = None, stdout_path: Path | None = None) -> None:
    log_handle = log_path.open("w", encoding="utf-8") if log_path else None
    stdout_handle = stdout_path.open("w", encoding="utf-8") if stdout_path else None
    try:
        subprocess.run(cmd, stderr=log_handle, stdout=stdout_handle, check=True)
    finally:
        if log_handle:
            log_handle.close()
        if stdout_handle:
            stdout_handle.close()

def subsample_pair(row: dict[str, str], outdir: Path, pairs: int, seed: int, force: bool) -> tuple[Path, Path]:
    sample = row["sample_id"]
    out1 = outdir / f"{sample}.R1.fastq"
    out2 = outdir / f"{sample}.R2.fastq"
    if not force and out1.exists() and out2.exists() and out1.stat().st_size > 0 and out2.stat().st_size > 0:
        return out1, out2
    run(["seqtk", "sample", f"-s{seed}", row["read1_path"], str(pairs)], stdout_path=out1)
    run(["seqtk", "sample", f"-s{seed}", row["read2_path"], str(pairs)], stdout_path=out2)
    return out1, out2

def original_pair(row: dict[str, str]) -> tuple[Path, Path]:
    read1 = Path(row["read1_path"])
    read2 = Path(row["read2_path"])
    if not read1.exists() or read1.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty R1 for {row['sample_id']}: {read1}")
    if not read2.exists() or read2.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty R2 for {row['sample_id']}: {read2}")
    return read1, read2

def align_pair(
    sample: str,
    ref_fasta: Path,
    read1: Path,
    read2: Path,
    bam_path: Path,
    log_path: Path,
    samtools_exe: str,
    threads: int,
    sort_threads: int,
    sort_tmp_dir: Path | None,
    keep_unmapped: bool,
) -> None:
    bwa_cmd = ["bwa", "mem", "-t", str(threads), str(ref_fasta), str(read1), str(read2)]
    view_cmd = [samtools_exe, "view", "-bS", "-"]
    if not keep_unmapped:
        view_cmd[2:2] = ["-F", "4"]
    sort_cmd = [samtools_exe, "sort", "-@", str(sort_threads)]
    if sort_tmp_dir is not None:
        sort_tmp_dir.mkdir(parents=True, exist_ok=True)
        sort_cmd.extend(["-T", str(sort_tmp_dir / sample)])
    sort_cmd.extend(["-o", str(bam_path)])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=log_handle)
        view = subprocess.Popen(view_cmd, stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=log_handle)
        if bwa.stdout:
            bwa.stdout.close()
        sort = subprocess.Popen(sort_cmd, stdin=view.stdout, stderr=log_handle)
        if view.stdout:
            view.stdout.close()
        sort_rc = sort.wait()
        view_rc = view.wait()
        bwa_rc = bwa.wait()
    if bwa_rc != 0 or view_rc != 0 or sort_rc != 0:
        raise subprocess.CalledProcessError(sort_rc or view_rc or bwa_rc, f"alignment pipeline for {sample}")
    run([samtools_exe, "index", str(bam_path)])

def bam_is_usable(bam_path: Path, samtools_exe: str) -> bool:
    bai_path = bam_path.with_name(f"{bam_path.name}.bai")
    if not bam_path.exists() or bam_path.stat().st_size == 0 or not bai_path.exists():
        return False
    result = subprocess.run([samtools_exe, "quickcheck", str(bam_path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return result.returncode == 0

def write_idxstats(bam_path: Path, idxstats_path: Path, samtools_exe: str) -> None:
    run([samtools_exe, "idxstats", str(bam_path)], stdout_path=idxstats_path)

def raw_idxstats_counts(idxstats_path: Path, reference_map: dict[str, str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    with idxstats_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            ref, mapped = fields[0], int(fields[2])
            category = reference_map.get(ref, ref)
            counts[category] = counts.get(category, 0) + mapped
    return counts

def select_rows(rows: list[dict[str, str]], sample_ids: list[str]) -> list[dict[str, str]]:
    if not sample_ids:
        return rows
    requested = set(sample_ids)
    selected = [row for row in rows if row["sample_id"] in requested]
    missing = requested - {row["sample_id"] for row in selected}
    if missing:
        raise SystemExit(f"Requested sample_id not found in manifest: {', '.join(sorted(missing))}")
    return selected

def process_sample(
    index: int,
    row: dict[str, str],
    args: argparse.Namespace,
    reference_map: dict[str, str],
    categories: list[str],
    fastq_dir: Path,
    bam_dir: Path,
    result_dir: Path,
    log_dir: Path,
) -> tuple[int, dict[str, object], dict[str, object], dict[str, object], dict[str, object], dict[str, object]]:
    sample = row["sample_id"]
    if args.full_input:
        read1, read2 = original_pair(row)
    else:
        read1, read2 = subsample_pair(row, fastq_dir, args.pairs, args.seed, args.force_subsample)
    bam = bam_dir / f"{sample}.retrovirus.bam"
    idxstats = result_dir / f"{sample}.idxstats.tsv"
    status = "reused_bam" if args.resume and bam_is_usable(bam, args.samtools_exe) else "aligned"
    if status == "aligned":
        align_pair(
            sample,
            args.reference_fasta,
            read1,
            read2,
            bam,
            log_dir / f"{sample}.bwa.log",
            args.samtools_exe,
            args.threads,
            args.sort_threads,
            args.sort_tmp_dir,
            args.keep_unmapped,
        )
    if not args.resume or not idxstats.exists() or idxstats.stat().st_size == 0:
        write_idxstats(bam, idxstats, args.samtools_exe)

    raw = raw_idxstats_counts(idxstats, reference_map)
    filt, filtered_records, dedup_removed = filtered_counts(
        bam,
        reference_map,
        args.min_mapq,
        args.category_min_mapq_map,
        args.min_aligned_length,
        args.samtools_exe,
        args.dedup_mode,
        args.umi_regex,
        args.require_unique_best,
    )
    raw_row = {"sample": sample, **{category: raw.get(category, 0) for category in categories}}
    filtered_row = {"sample": sample, **{category: filt.get(category, 0) for category in categories}}
    filtered_record_row = {"sample": sample, **{category: filtered_records.get(category, 0) for category in categories}}
    dedup_removed_row = {"sample": sample, **{category: dedup_removed.get(category, 0) for category in categories}}
    run_row = {
        "sample": sample,
        "status": status,
        "read1": str(read1),
        "read2": str(read2),
        "bam": str(bam),
        "idxstats": str(idxstats),
        "pairs": "full" if args.full_input else args.pairs,
        "jobs": args.jobs,
        "threads": args.threads,
        "sort_threads": args.sort_threads,
        "min_mapq": args.min_mapq,
        "category_min_mapq": args.category_min_mapq_text,
        "min_aligned_length": args.min_aligned_length,
        "dedup_mode": args.dedup_mode,
        "umi_regex": args.umi_regex or "",
        "require_unique_best": args.require_unique_best,
        "sort_tmp_dir": str(args.sort_tmp_dir or ""),
        "keep_unmapped": args.keep_unmapped,
    }
    return index, raw_row, filtered_row, filtered_record_row, dedup_removed_row, run_row

def main() -> int:
    parser = argparse.ArgumentParser(description="Run subsampled BWA retrovirus pilot and filtered category counts.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-map", type=Path, required=True)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--pairs", type=int, default=1_000_000)
    parser.add_argument("--full-input", action="store_true", help="Use original FASTQs directly instead of subsampling.")
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--sort-threads", type=int, default=4)
    parser.add_argument("--samtools-exe", default="samtools")
    parser.add_argument("--jobs", type=int, default=1, help="Samples to process concurrently.")
    parser.add_argument("--sort-tmp-dir", type=Path, help="Directory for samtools sort temporary files.")
    parser.add_argument("--keep-unmapped", action="store_true", help="Keep unmapped reads in output BAMs.")
    parser.add_argument("--resume", action="store_true", help="Reuse existing valid BAM/index outputs.")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--category-min-mapq", action="append", default=[], help="Override MAPQ by category, e.g. HUMAN:60.")
    parser.add_argument("--min-aligned-length", type=int, default=60)
    parser.add_argument("--dedup-mode", choices=["none", "read_id", "coordinate", "fragment", "umi"], default="read_id")
    parser.add_argument("--umi-regex", help="Regex used with --dedup-mode umi.")
    parser.add_argument("--require-unique-best", action="store_true", help="Drop records with AS <= XS when both tags exist.")
    parser.add_argument("--force-subsample", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if args.dedup_mode == "umi" and not args.umi_regex:
        raise SystemExit("--dedup-mode umi requires --umi-regex")

    rows = select_rows(read_csv(args.manifest), args.sample_id)
    reference_map = load_reference_map(args.reference_map)
    categories = sorted(set(reference_map.values()))
    args.category_min_mapq_map = parse_category_min_mapq(args.category_min_mapq)
    unknown_categories = sorted(set(args.category_min_mapq_map) - set(categories))
    if unknown_categories:
        raise SystemExit(f"--category-min-mapq uses categories not found in reference map: {', '.join(unknown_categories)}")
    args.category_min_mapq_text = ",".join(f"{category}:{args.category_min_mapq_map[category]}" for category in sorted(args.category_min_mapq_map))

    fastq_dir = args.work_dir / ("fastq_full_inputs" if args.full_input else f"fastq_{args.pairs}")
    bam_dir = args.work_dir / "bam"
    result_dir = args.work_dir / "results"
    log_dir = args.work_dir / "logs"
    for path in (fastq_dir, bam_dir, result_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)

    results: list[tuple[int, dict[str, object], dict[str, object], dict[str, object], dict[str, object], dict[str, object]]] = []
    if args.jobs == 1:
        for index, row in enumerate(rows):
            print(f"Running {row['sample_id']}")
            result = process_sample(index, row, args, reference_map, categories, fastq_dir, bam_dir, result_dir, log_dir)
            results.append(result)
            print(f"Completed {row['sample_id']} ({result[5]['status']})")
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(process_sample, index, row, args, reference_map, categories, fastq_dir, bam_dir, result_dir, log_dir): row["sample_id"]
                for index, row in enumerate(rows)
            }
            for completed, future in enumerate(as_completed(futures), start=1):
                sample = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    raise SystemExit(f"Sample failed: {sample}: {exc}") from exc
                results.append(result)
                print(f"Completed {completed}/{len(futures)} {sample} ({result[5]['status']})")

    ordered = sorted(results, key=lambda item: item[0])
    raw_rows = [item[1] for item in ordered]
    filtered_rows = [item[2] for item in ordered]
    filtered_record_rows = [item[3] for item in ordered]
    dedup_removed_rows = [item[4] for item in ordered]
    run_rows = [item[5] for item in ordered]

    write_csv(result_dir / "raw_idxstats_category_counts.tsv", raw_rows, ["sample", *categories])
    write_csv(result_dir / "filtered_record_category_counts.tsv", filtered_record_rows, ["sample", *categories])
    write_csv(result_dir / "filtered_category_counts.tsv", filtered_rows, ["sample", *categories])
    write_csv(result_dir / "dedup_removed_category_counts.tsv", dedup_removed_rows, ["sample", *categories])
    write_csv(result_dir / "run_manifest.tsv", run_rows, list(run_rows[0].keys()) if run_rows else ["sample"])
    print(f"Wrote {result_dir / 'filtered_category_counts.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
