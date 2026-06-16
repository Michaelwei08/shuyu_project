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

## 4. Masked Full Rerun

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

## 5. Old-versus-New Comparison

The defaults use the existing viral-only and hg38-inclusive run directories. Override
`OLDHTLV`, `OLDWGS`, `FULLHTLV`, or `WGSFULL` if the paths differ.

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh compare
cat herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/method_comparison/method_comparison_report.md
```
