# Final Four P1 Validation Steps

Run from the remote repository after activating `retro_qc`.

## 1. Known HTLV Status Labels

The evaluator requires `sample_id` and `status`; additional metadata columns are allowed:

```text
sample_id<TAB>status
targeted_htlv_tcl_TCL201-C1_S07<TAB>positive
```

Allowed statuses are `positive`, `negative`, and `unknown`. Then run:

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh label-template
```

This creates `output/targeted_htlv_status_to_fill.tsv` with all exact sample IDs. Have
Shuyu fill the `status` column, move the completed file to protected metadata storage,
and then run:

```bash
export HTLV_LABELS=/drive3/cpwei/private_metadata/targeted_htlv_status.tsv
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh labels
```

## 2. Targeted HIV Data

When Shuyu provides the paired FASTQs:

```bash
export TARGETED_HIV_DIR=/drive3/Shuyu_for_Michael/Targeted_HIV_samples
export JOBS=8 THREADS=8 SORT_THREADS=2
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh targeted-hiv
```

The script fails if FASTQ mates are missing or duplicated. Pairing issues are written to
`shuyu_benchmark_package/output/targeted_hiv_pairing_issues.csv`. It also creates
`targeted_hiv_status_to_fill.tsv`. After authoritative labels are filled, rerun with:

```bash
export HIV_LABELS=/drive3/cpwei/private_metadata/targeted_hiv_status.tsv
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh targeted-hiv
```

`--resume` reuses the completed BAMs, so the second run performs post-processing rather than repeating alignment.

## 3. Final Tables and Figures

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh finalize
```

Figures are written to `shuyu_benchmark_package/output/final_rerun_figures/`.

## 4. Old-versus-New Comparison

The defaults use the existing viral-only and hg38-inclusive run directories. Override
`OLDHTLV`, `OLDWGS`, `FULLHTLV`, or `WGSFULL` if the paths differ.

```bash
bash herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/run_p1_final_validation.sh compare
cat herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/method_comparison/method_comparison_report.md
```
