# Remote Final Two Pre-Full-Run Steps

These commands assume the `retro_qc` conda environment is active and that `bwa`, `samtools`, `seqtk`, and `efetch` are available.

## 1. Build Competitive HIV/HTLV/HERV/LINE1 Reference

```bash
cd /home/alizadehlab/cpwei/shuyu_project
WORK=/home/alizadehlab/cpwei/shuyu_project/local_work/retro_competitive
mkdir -p "$WORK/ref"

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/make_retro_competitive_reference.py \
  --from-efetch-defaults \
  --output-fasta "$WORK/ref/retro_competitive.fa" \
  --output-map "$WORK/ref/retro_competitive_reference_map.csv" \
  --index

cat "$WORK/ref/retro_competitive_reference_map.csv"
```

Default records:

- HIV-1 HXB2: `K03455`
- HIV-2 ROD: `M15390`
- HTLV-1 ATK: `J02029`
- HTLV-2 Mo: `M10060`
- HERV-K113: `AY037928`
- LINE1 L1.3: `L19088`

If better local HERV/LINE1 consensus FASTAs exist, append them with repeated `--extra-fasta CATEGORY:/path/to/file.fa`.

## 2. Re-run 7-Sample Pilot With Competitive Reference

```bash
cd /home/alizadehlab/cpwei/shuyu_project
REFWORK=/home/alizadehlab/cpwei/shuyu_project/local_work/retro_competitive
PILOTWORK=/home/alizadehlab/cpwei/shuyu_project/local_work/shuyu_pilot7_competitive

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_retro_pilot_alignment.py \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/pilot_manifest.csv \
  --work-dir "$PILOTWORK" \
  --reference-fasta "$REFWORK/ref/retro_competitive.fa" \
  --reference-map "$REFWORK/ref/retro_competitive_reference_map.csv" \
  --pairs 1000000 \
  --threads 8 \
  --min-mapq 20 \
  --min-aligned-length 60

cat "$PILOTWORK/results/filtered_category_counts.tsv"
```

Expected gate:

- HL WGS controls stay zero or near-zero for HIV and HTLV.
- `TCL201`, `TCL202`, and `TCL203` retain HTLV1 signal after HERV/LINE1 decoys are present.
- Spurious short HIV signal stays removed.

## 3. Larger WGS HIV-vs-HL Check

Run 5 million pairs for one HIV WGS and one HL WGS control:

```bash
cd /home/alizadehlab/cpwei/shuyu_project
REFWORK=/home/alizadehlab/cpwei/shuyu_project/local_work/retro_competitive
WGSWORK=/home/alizadehlab/cpwei/shuyu_project/local_work/shuyu_wgs_5m_competitive

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_retro_pilot_alignment.py \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/sample_manifest.csv \
  --work-dir "$WGSWORK" \
  --reference-fasta "$REFWORK/ref/retro_competitive.fa" \
  --reference-map "$REFWORK/ref/retro_competitive_reference_map.csv" \
  --sample-id wgs_60samples_hiv_hl_HIV001-P1_S01 \
  --sample-id wgs_60samples_hiv_hl_HL350-P1_S02 \
  --pairs 5000000 \
  --threads 8 \
  --min-mapq 20 \
  --min-aligned-length 60

cat "$WGSWORK/results/filtered_category_counts.tsv"
```

Expected gate:

- HL control should remain zero for HIV/HTLV.
- HIV WGS may remain zero at 5M pairs; that does not fail the method because WGS HIV signal may be extremely sparse.
- If HIV appears in HIV WGS but not HL control, proceed to larger HIV WGS sampling or full WGS scan.

## Full-Round Decision

Proceed to full targeted HTLV first if:

- 7-sample competitive pilot keeps HTLV targeted samples positive.
- HL controls remain negative.
- HERV/LINE1 decoys do not absorb all HTLV signal.

## 4. Full Targeted HTLV Run

Only run complete targeted FASTQ pairs. This uses original FASTQs directly and does not copy/subsample reads:

```bash
cd /home/alizadehlab/cpwei/shuyu_project
REFWORK=/home/alizadehlab/cpwei/shuyu_project/local_work/retro_competitive
FULLHTLV=/drive3/cpwei/shuyu_runs/targeted_htlv_full_competitive
SORTTMP=/drive3/cpwei/tmp/samtools_sort
mkdir -p "$FULLHTLV" "$SORTTMP"

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_retro_pilot_alignment.py \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/targeted_htlv_complete_manifest.csv \
  --work-dir "$FULLHTLV" \
  --reference-fasta "$REFWORK/ref/retro_competitive.fa" \
  --reference-map "$REFWORK/ref/retro_competitive_reference_map.csv" \
  --full-input \
  --jobs 4 \
  --threads 16 \
  --sort-threads 2 \
  --sort-tmp-dir "$SORTTMP" \
  --resume \
  --min-mapq 20 \
  --min-aligned-length 60

cat "$FULLHTLV/results/filtered_category_counts.tsv"
```

By default, this stores only mapped reads in the BAMs. Use `--keep-unmapped` only for debugging, not for full FASTQ runs, because all-read BAMs can exceed home quotas quickly.

`--jobs` controls sample-level parallelism. `--threads` controls BWA threads per active sample. Total CPU pressure is approximately `jobs * (threads + sort_threads)`, plus some single-threaded counting overhead. On a large shared server, start with `--jobs 4 --threads 16 --sort-threads 2`; increase only if CPU is high and `/drive3` I/O is not the bottleneck.

Monitor progress from another terminal:

```bash
FULLHTLV=/drive3/cpwei/shuyu_runs/targeted_htlv_full_competitive
ls "$FULLHTLV/results"/*.idxstats.tsv 2>/dev/null | wc -l
du -sh "$FULLHTLV"
ps -u cpwei -o pid,etime,pcpu,pmem,cmd | grep -E 'bwa|samtools|run_retro' | grep -v grep
```

Proceed to full WGS only after choosing whether WGS should use:

- larger subsampling,
- all-read viral alignment,
- or human-first unmapped/poorly mapped read extraction.

## 5. WGS HIV-vs-HL 5M-Pair Screen

Run all 60 WGS samples at 5 million pairs per sample before considering full WGS. This verifies control specificity without committing to the full WGS compute/storage cost.

```bash
cd /home/alizadehlab/cpwei/shuyu_project

python - <<'PY'
import csv
inp = "herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/sample_manifest.csv"
out = "herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/wgs_complete_manifest.csv"
rows = [
    row for row in csv.DictReader(open(inp))
    if row["assay_type"] == "WGS" and row["read1_path"] and row["read2_path"]
]
with open(out, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"wrote {len(rows)} WGS rows to {out}")
PY

REFWORK=/drive3/cpwei/shuyu_runs/archive_before_home_cleanup_2026-05-17/retro_competitive
WGSWORK=/drive3/cpwei/shuyu_runs/wgs_hiv_hl_5m_competitive
SORTTMP=/drive3/cpwei/tmp/samtools_sort
mkdir -p "$WGSWORK" "$SORTTMP"

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_retro_pilot_alignment.py \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/wgs_complete_manifest.csv \
  --work-dir "$WGSWORK" \
  --reference-fasta "$REFWORK/ref/retro_competitive.fa" \
  --reference-map "$REFWORK/ref/retro_competitive_reference_map.csv" \
  --pairs 5000000 \
  --jobs 4 \
  --threads 16 \
  --sort-threads 2 \
  --sort-tmp-dir "$SORTTMP" \
  --resume \
  --min-mapq 20 \
  --min-aligned-length 60

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/summarize_wgs_retro_results.py \
  --counts "$WGSWORK/results/filtered_category_counts.tsv" \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/wgs_complete_manifest.csv \
  --output-dir "$WGSWORK/results/final_summary"

cat "$WGSWORK/results/final_summary/wgs_retro_report.md"
```

Monitor progress:

```bash
WGSWORK=/drive3/cpwei/shuyu_runs/wgs_hiv_hl_5m_competitive
find "$WGSWORK/results" -name "*.idxstats.tsv" | wc -l
ps -u cpwei -o pid,etime,pcpu,pmem,cmd | grep -E 'bwa|samtools|run_retro' | grep -v grep
du -sh "$WGSWORK"
```
