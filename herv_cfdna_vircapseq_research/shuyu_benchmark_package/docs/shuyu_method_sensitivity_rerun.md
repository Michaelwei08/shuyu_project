# Shuyu Method Sensitivity Rerun

This rerun is designed to test whether HIV/HTLV calls are stable and whether ubiquitous HERV counts are biological background or a reference/mapping artifact.

## 1. Build Alternative Reference Panels

Current isolate panel:

```bash
REF_CURRENT=/drive3/cpwei/shuyu_runs/retro_reference_current
mkdir -p "$REF_CURRENT/ref"

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/make_retro_competitive_reference.py \
  --from-efetch-defaults \
  --panel current \
  --output-fasta "$REF_CURRENT/ref/retro_competitive.current.fa" \
  --output-map "$REF_CURRENT/ref/retro_competitive.current.reference_map.csv" \
  --index
```

RefSeq panel using Shuyu's requested `NC_*` accessions:

```bash
REF_REFSEQ=/drive3/cpwei/shuyu_runs/retro_reference_refseq
mkdir -p "$REF_REFSEQ/ref"

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/make_retro_competitive_reference.py \
  --from-efetch-defaults \
  --panel refseq \
  --output-fasta "$REF_REFSEQ/ref/retro_competitive.refseq.fa" \
  --output-map "$REF_REFSEQ/ref/retro_competitive.refseq.reference_map.csv" \
  --index
```

Full human-plus-virus competitive reference, if an hg38 FASTA is available:

```bash
HG38=/path/to/hg38.fa
REF_HG38=/drive3/cpwei/shuyu_runs/retro_reference_hg38_refseq
mkdir -p "$REF_HG38/ref"

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/make_retro_competitive_reference.py \
  --from-efetch-defaults \
  --panel refseq \
  --human-fasta "$HG38" \
  --output-fasta "$REF_HG38/ref/hg38_plus_retro.refseq.fa" \
  --output-map "$REF_HG38/ref/hg38_plus_retro.refseq.reference_map.csv" \
  --index
```

## 2. Run Stricter Alignment/Counting Conditions

Example targeted HTLV rerun with stricter viral MAPQ and coordinate deduplication:

```bash
RUN=/drive3/cpwei/shuyu_runs/targeted_htlv_refseq_mapq40_coord_dedup
SORTTMP=/drive3/cpwei/tmp/samtools_sort

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_retro_pilot_alignment.py \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/targeted_htlv_complete_manifest.csv \
  --work-dir "$RUN" \
  --reference-fasta "$REF_REFSEQ/ref/retro_competitive.refseq.fa" \
  --reference-map "$REF_REFSEQ/ref/retro_competitive.refseq.reference_map.csv" \
  --full-input \
  --jobs 10 \
  --threads 16 \
  --sort-threads 2 \
  --sort-tmp-dir "$SORTTMP" \
  --resume \
  --min-mapq 40 \
  --min-aligned-length 60 \
  --dedup-mode coordinate \
  --require-unique-best
```

Example hg38-plus-virus run with category-specific MAPQ thresholds:

```bash
RUN=/drive3/cpwei/shuyu_runs/wgs_hg38_refseq_mapq_human60_viral40
SORTTMP=/drive3/cpwei/tmp/samtools_sort

python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_retro_pilot_alignment.py \
  --manifest herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/wgs_complete_manifest.csv \
  --work-dir "$RUN" \
  --reference-fasta "$REF_HG38/ref/hg38_plus_retro.refseq.fa" \
  --reference-map "$REF_HG38/ref/hg38_plus_retro.refseq.reference_map.csv" \
  --full-input \
  --jobs 4 \
  --threads 16 \
  --sort-threads 2 \
  --sort-tmp-dir "$SORTTMP" \
  --resume \
  --min-mapq 40 \
  --category-min-mapq HUMAN:60 \
  --min-aligned-length 60 \
  --dedup-mode coordinate \
  --require-unique-best
```

## 3. Outputs To Compare

Each rerun writes:

```text
results/raw_idxstats_category_counts.tsv
results/filtered_record_category_counts.tsv
results/filtered_category_counts.tsv
results/dedup_removed_category_counts.tsv
results/run_manifest.tsv
```

Use `filtered_record_category_counts.tsv` versus `filtered_category_counts.tsv` to measure whether deduplication changed the counts. Use `dedup_removed_category_counts.tsv` to quantify what was removed.

The key HERV specificity question is whether HERV remains ubiquitous after:

- RefSeq viral accessions.
- hg38-plus-retroelement competitive reference.
- MAPQ 40 viral filtering and HUMAN MAPQ 60 filtering.
- coordinate or UMI-aware deduplication.
- `--require-unique-best` filtering.
