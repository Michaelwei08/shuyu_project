from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


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


def aligned_len(cigar: str) -> int:
    total = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in "M=X":
            total += int(n)
    return total


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
    threads: int,
    sort_tmp_dir: Path | None,
) -> None:
    bwa_cmd = ["bwa", "mem", "-t", str(threads), str(ref_fasta), str(read1), str(read2)]
    view_cmd = ["samtools", "view", "-bS", "-"]
    sort_cmd = ["samtools", "sort", "-@", "4"]
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
    run(["samtools", "index", str(bam_path)])


def write_idxstats(bam_path: Path, idxstats_path: Path) -> None:
    run(["samtools", "idxstats", str(bam_path)], stdout_path=idxstats_path)


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


def filtered_counts(
    bam_path: Path,
    reference_map: dict[str, str],
    min_mapq: int,
    min_aligned_length: int,
) -> dict[str, int]:
    counts: dict[str, int] = {}
    seen: set[tuple[str, str]] = set()
    result = subprocess.run(["samtools", "view", "-F", "4", str(bam_path)], text=True, stdout=subprocess.PIPE, check=True)
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        read_id, ref, mapq_text, cigar = fields[0], fields[2], fields[4], fields[5]
        if ref not in reference_map:
            continue
        if int(mapq_text) < min_mapq:
            continue
        if aligned_len(cigar) < min_aligned_length:
            continue
        category = reference_map[ref]
        key = (read_id, category)
        if key in seen:
            continue
        seen.add(key)
        counts[category] = counts.get(category, 0) + 1
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


def main() -> int:
    parser = argparse.ArgumentParser(description="Run subsampled BWA retrovirus pilot and filtered category counts.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-map", type=Path, required=True)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--pairs", type=int, default=1_000_000)
    parser.add_argument(
        "--full-input",
        action="store_true",
        help="Use original FASTQs directly instead of subsampling with seqtk. Intended for full targeted runs.",
    )
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--sort-tmp-dir",
        type=Path,
        help="Directory for samtools sort temporary files. Use a local filesystem such as /tmp/cpwei_sort if network storage gives Illegal seek.",
    )
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-aligned-length", type=int, default=60)
    parser.add_argument("--force-subsample", action="store_true")
    args = parser.parse_args()

    rows = select_rows(read_csv(args.manifest), args.sample_id)
    reference_map = load_reference_map(args.reference_map)
    categories = sorted(set(reference_map.values()))

    fastq_dir = args.work_dir / ("fastq_full_inputs" if args.full_input else f"fastq_{args.pairs}")
    bam_dir = args.work_dir / "bam"
    result_dir = args.work_dir / "results"
    log_dir = args.work_dir / "logs"
    for path in (fastq_dir, bam_dir, result_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)

    raw_rows: list[dict[str, object]] = []
    filtered_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []

    for row in rows:
        sample = row["sample_id"]
        print(f"Running {sample}")
        if args.full_input:
            read1, read2 = original_pair(row)
        else:
            read1, read2 = subsample_pair(row, fastq_dir, args.pairs, args.seed, args.force_subsample)
        bam = bam_dir / f"{sample}.retrovirus.bam"
        idxstats = result_dir / f"{sample}.idxstats.tsv"
        align_pair(
            sample,
            args.reference_fasta,
            read1,
            read2,
            bam,
            log_dir / f"{sample}.bwa.log",
            args.threads,
            args.sort_tmp_dir,
        )
        write_idxstats(bam, idxstats)

        raw = raw_idxstats_counts(idxstats, reference_map)
        filt = filtered_counts(bam, reference_map, args.min_mapq, args.min_aligned_length)
        raw_rows.append({"sample": sample, **{category: raw.get(category, 0) for category in categories}})
        filtered_rows.append({"sample": sample, **{category: filt.get(category, 0) for category in categories}})
        run_rows.append(
            {
                "sample": sample,
                "read1": str(read1),
                "read2": str(read2),
                "bam": str(bam),
                "idxstats": str(idxstats),
                "pairs": "full" if args.full_input else args.pairs,
                "min_mapq": args.min_mapq,
                "min_aligned_length": args.min_aligned_length,
                "sort_tmp_dir": str(args.sort_tmp_dir or ""),
            }
        )

    write_csv(result_dir / "raw_idxstats_category_counts.tsv", raw_rows, ["sample", *categories])
    write_csv(result_dir / "filtered_category_counts.tsv", filtered_rows, ["sample", *categories])
    write_csv(result_dir / "run_manifest.tsv", run_rows, list(run_rows[0].keys()) if run_rows else ["sample"])
    print(f"Wrote {result_dir / 'filtered_category_counts.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
