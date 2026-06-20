# Remaining P1 Validation Steps

Run from the remote repository after activating `retro_qc`.

## 1. Targeted HIV Data

When Shuyu provides the paired FASTQs:

```bash
export TARGETED_HIV_DIR=/drive3/Shuyu_for_Michael/Targeted_HIV_samples
export JOBS=8 THREADS=8 SORT_THREADS=2
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh targeted-hiv
```

The script fails if FASTQ mates are missing or duplicated. Pairing issues are written to
`shuyu_benchmark_package/output/targeted_hiv_pairing_issues.csv`.

## 2. Final Tables and Figures

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh finalize
```

Figures are written to `shuyu_benchmark_package/output/final_rerun_figures/`.

## 2b. Primary-Only Sensitivity Check

This creates the second version Shuyu requested. It reuses the existing BAMs and writes
`primary_only_*` count tables after excluding secondary (`0x100`) and supplementary (`0x800`)
SAM alignments.

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh primary-only
```

Key outputs:

- `$FULLHTLV/results/primary_only_filtered_category_counts.tsv`
- `$FULLHTLV/results/final_summary_primary_only/`
- `$WGSFULL/results/primary_only_filtered_category_counts.tsv`
- `$WGSFULL/results/final_summary_primary_only/`

## 3. HERV-Shared Region Calculation and Masking

This calculates exact shared-kmer similarity between HIV1/HERV and HTLV1/HERV, masks the
shared HIV1/HTLV1 positions, and indexes a new masked reference. The default k-mer size is 40.

```bash
export MASK_KMER=40
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh mask-reference
```

Key outputs:

- `$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.pair_similarity.tsv`
- `$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.mask_summary.tsv`
- `$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.mask.bed`
- `$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.mask_report.md`
- `$MASKED_REFERENCE_FASTA`

The pair-summary TSV contains both pre-mask and post-mask exact-kmer similarity:

- `query_similar_pct`: pre-mask percent of HIV1/HTLV1 bases covered by HERV-shared k-mers.
- `post_mask_query_similar_pct`: residual percent after applying the mask.
- `retained_pct` in the mask-summary TSV: genome fraction retained after masking.

`NC_007605.1_masked.fa` is an EBV masked reference, so it is not sufficient for this
HIV/HTLV/HERV question unless Shuyu provides the matching VirCAPP masking method or masked
HIV/HTLV/HERV reference files.

## 4. Shuyu Masked Viral Panel Follow-Up

Shuyu provided the panel path:

```text
/drive3/shuyu/references/HIV1_masked/viral_sel_v1_MASKED_HIV1masked.fa
```

First export the BAM paths from the current hg38 rerun so Shuyu can inspect them in IGV:

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_shuyu_masked_panel_validation.sh export-current-igv
```

Key outputs:

- `shuyu_benchmark_package/output/igv_bam_paths/current_primary_only_all_bams.tsv`
- `shuyu_benchmark_package/output/igv_bam_paths/current_primary_only_nonzero_exogenous_bams.tsv`

To test Shuyu's 180-virus masked panel directly, build a new reference containing hg38/HERV/LINE1
from the current base reference plus Shuyu's masked viral panel:

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_shuyu_masked_panel_validation.sh build-reference
```

Then rerun full targeted HTLV and WGS with secondary/supplementary alignments excluded:

```bash
export JOBS=8 THREADS=8 SORT_THREADS=2
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_shuyu_masked_panel_validation.sh rerun
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_shuyu_masked_panel_validation.sh summarize
```

The Shuyu-panel rerun outputs are written to:

- `/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_shuyu_masked_panel_primary_only`
- `/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_shuyu_masked_panel_primary_only`
- `shuyu_benchmark_package/output/igv_bam_paths/shuyu_panel_all_bams.tsv`

## 5. Exact-Kmer Masked Full Rerun

After `mask-reference` finishes, rerun targeted HTLV and HIV/HL WGS against the masked reference:

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh masked-rerun
```

Key outputs:

- `$FULLHTLV_MASKED/results/filtered_category_counts.tsv`
- `$FULLHTLV_MASKED/results/final_summary/`
- `$WGSFULL_MASKED/results/filtered_category_counts.tsv`
- `$WGSFULL_MASKED/results/final_summary/`

This rerun also excludes secondary and supplementary alignments, applies human MAPQ >=60,
viral MAPQ >=40, aligned length >=60 bp, unique-best filtering, and coordinate deduplication.

## 6. Old-versus-New Comparison

The defaults use the existing viral-only and hg38-inclusive run directories. Override
`OLDHTLV`, `OLDWGS`, `FULLHTLV`, or `WGSFULL` if the paths differ.

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh compare
cat herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/method_comparison/method_comparison_report.md
```
