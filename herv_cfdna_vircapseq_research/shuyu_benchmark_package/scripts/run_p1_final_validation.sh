#!/usr/bin/env bash
set -euo pipefail

STEP="${1:-help}"
PROJECT="${PROJECT:-/home/alizadehlab/cpwei/shuyu_project}"
PACKAGE="$PROJECT/herv_cfdna_vircapseq_research/shuyu_benchmark_package"
SCRIPTS="$PACKAGE/scripts"
OUT="${OUT:-$PACKAGE/output}"

FULLHTLV="${FULLHTLV:-/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord}"
WGSFULL="${WGSFULL:-/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord}"
OLDHTLV="${OLDHTLV:-/drive3/cpwei/shuyu_runs/targeted_htlv_full_competitive}"
OLDWGS="${OLDWGS:-/drive3/cpwei/shuyu_runs/wgs_hiv_hl_full_competitive}"
REFDIR="${REFDIR:-/drive3/cpwei/shuyu_runs/retro_reference_hg38_refseq/ref}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REFDIR/hg38_plus_retro.refseq.fa}"
REFERENCE_MAP="${REFERENCE_MAP:-$REFDIR/hg38_plus_retro.refseq.reference_map.csv}"
SORTTMP="${SORTTMP:-/drive3/cpwei/tmp/samtools_sort}"
TARGETED_HIV_WORK="${TARGETED_HIV_WORK:-/drive3/cpwei/shuyu_runs/targeted_hiv_hg38_refseq_mapq_human60_viral40_coord}"
MASK_KMER="${MASK_KMER:-40}"
MASKDIR="${MASKDIR:-/drive3/cpwei/shuyu_runs/retro_reference_hg38_refseq_masked_hiv1_htlv1_vs_herv_k${MASK_KMER}/ref}"
MASKED_REFERENCE_FASTA="${MASKED_REFERENCE_FASTA:-$MASKDIR/hg38_plus_retro.refseq.masked_hiv1_htlv1_vs_herv.k${MASK_KMER}.fa}"
MASKED_REFERENCE_MAP="${MASKED_REFERENCE_MAP:-$MASKDIR/hg38_plus_retro.refseq.masked_hiv1_htlv1_vs_herv.k${MASK_KMER}.reference_map.csv}"
FULLHTLV_MASKED="${FULLHTLV_MASKED:-/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_refseq_masked_k${MASK_KMER}_mapq_human60_viral40_coord}"
WGSFULL_MASKED="${WGSFULL_MASKED:-/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_refseq_masked_k${MASK_KMER}_mapq_human60_viral40_coord}"

require_file() {
  [[ -f "$1" ]] || { echo "Required file not found: $1" >&2; exit 2; }
}

require_dir() {
  [[ -d "$1" ]] || { echo "Required directory not found: $1" >&2; exit 2; }
}

step_targeted_hiv() {
  [[ -n "${TARGETED_HIV_DIR:-}" ]] || { echo "Required environment variable is unset: TARGETED_HIV_DIR" >&2; exit 2; }
  require_dir "$TARGETED_HIV_DIR"
  require_file "$REFERENCE_FASTA"
  require_file "$REFERENCE_MAP"
  mkdir -p "$OUT" "$SORTTMP" "$TARGETED_HIV_WORK"

  python "$SCRIPTS/prepare_targeted_hiv_manifest.py" \
    --data-dir "$TARGETED_HIV_DIR" \
    --output "$OUT/targeted_hiv_complete_manifest.csv" \
    --issues-output "$OUT/targeted_hiv_pairing_issues.csv" \
    --hiv-status unknown

  python "$SCRIPTS/run_retro_pilot_alignment.py" \
    --manifest "$OUT/targeted_hiv_complete_manifest.csv" \
    --work-dir "$TARGETED_HIV_WORK" \
    --reference-fasta "$REFERENCE_FASTA" \
    --reference-map "$REFERENCE_MAP" \
    --full-input \
    --jobs "${JOBS:-8}" \
    --threads "${THREADS:-8}" \
    --sort-threads "${SORT_THREADS:-2}" \
    --sort-tmp-dir "$SORTTMP" \
    --resume \
    --min-mapq 40 \
    --category-min-mapq HUMAN:60 \
    --min-aligned-length 60 \
    --dedup-mode coordinate \
    --require-unique-best

}

step_finalize() {
  require_file "$FULLHTLV/results/filtered_category_counts.tsv"
  require_file "$WGSFULL/results/filtered_category_counts.tsv"

  python "$SCRIPTS/summarize_targeted_htlv_results.py" \
    --counts "$FULLHTLV/results/filtered_category_counts.tsv" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --sample-qc "$OUT/manifest_qc/sample_file_qc.csv" \
    --output-dir "$FULLHTLV/results/final_summary" \
    --positive-threshold "${HTLV_THRESHOLD:-100}"

  python "$SCRIPTS/analyze_targeted_htlv_method_validation.py" \
    --filtered-record-counts "$FULLHTLV/results/filtered_record_category_counts.tsv" \
    --filtered-counts "$FULLHTLV/results/filtered_category_counts.tsv" \
    --dedup-removed-counts "$FULLHTLV/results/dedup_removed_category_counts.tsv" \
    --output-dir "$FULLHTLV/results/targeted_method_validation" \
    --positive-threshold "${HTLV_THRESHOLD:-100}"

  python "$SCRIPTS/summarize_wgs_retro_results.py" \
    --counts "$WGSFULL/results/filtered_category_counts.tsv" \
    --manifest "$OUT/wgs_complete_manifest.csv" \
    --output-dir "$WGSFULL/results/final_summary"

  python "$SCRIPTS/plot_p1_figures.py" \
    --targeted-summary "$FULLHTLV/results/final_summary/targeted_htlv_full_summary.tsv" \
    --wgs-summary "$WGSFULL/results/final_summary/wgs_retro_sample_summary.tsv" \
    --output-dir "$PACKAGE/output/final_rerun_figures" \
    --detection-threshold 1
}

step_compare() {
  require_file "$OLDHTLV/results/filtered_category_counts.tsv"
  require_file "$OLDWGS/results/filtered_category_counts.tsv"
  require_file "$FULLHTLV/results/filtered_category_counts.tsv"
  require_file "$WGSFULL/results/filtered_category_counts.tsv"
  python "$SCRIPTS/compare_retro_method_runs.py" \
    --old-targeted "$OLDHTLV/results/filtered_category_counts.tsv" \
    --new-targeted "$FULLHTLV/results/filtered_category_counts.tsv" \
    --old-wgs "$OLDWGS/results/filtered_category_counts.tsv" \
    --new-wgs "$WGSFULL/results/filtered_category_counts.tsv" \
    --output-dir "$PACKAGE/output/method_comparison"
}

step_primary_only() {
  require_file "$REFERENCE_FASTA"
  require_file "$REFERENCE_MAP"
  require_file "$OUT/targeted_htlv_complete_manifest.csv"
  require_file "$OUT/wgs_complete_manifest.csv"

  python "$SCRIPTS/run_retro_pilot_alignment.py" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --work-dir "$FULLHTLV" \
    --reference-fasta "$REFERENCE_FASTA" \
    --reference-map "$REFERENCE_MAP" \
    --full-input \
    --jobs "${JOBS:-8}" \
    --threads "${THREADS:-8}" \
    --sort-threads "${SORT_THREADS:-2}" \
    --sort-tmp-dir "$SORTTMP" \
    --resume \
    --min-mapq 40 \
    --category-min-mapq HUMAN:60 \
    --min-aligned-length 60 \
    --dedup-mode coordinate \
    --require-unique-best \
    --exclude-secondary-supplementary \
    --result-prefix primary_only_

  python "$SCRIPTS/run_retro_pilot_alignment.py" \
    --manifest "$OUT/wgs_complete_manifest.csv" \
    --work-dir "$WGSFULL" \
    --reference-fasta "$REFERENCE_FASTA" \
    --reference-map "$REFERENCE_MAP" \
    --full-input \
    --jobs "${JOBS:-8}" \
    --threads "${THREADS:-8}" \
    --sort-threads "${SORT_THREADS:-2}" \
    --sort-tmp-dir "$SORTTMP" \
    --resume \
    --min-mapq 40 \
    --category-min-mapq HUMAN:60 \
    --min-aligned-length 60 \
    --dedup-mode coordinate \
    --require-unique-best \
    --exclude-secondary-supplementary \
    --result-prefix primary_only_

  python "$SCRIPTS/summarize_targeted_htlv_results.py" \
    --counts "$FULLHTLV/results/primary_only_filtered_category_counts.tsv" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --sample-qc "$OUT/manifest_qc/sample_file_qc.csv" \
    --output-dir "$FULLHTLV/results/final_summary_primary_only" \
    --positive-threshold "${HTLV_THRESHOLD:-100}"

  python "$SCRIPTS/analyze_targeted_htlv_method_validation.py" \
    --filtered-record-counts "$FULLHTLV/results/primary_only_filtered_record_category_counts.tsv" \
    --filtered-counts "$FULLHTLV/results/primary_only_filtered_category_counts.tsv" \
    --dedup-removed-counts "$FULLHTLV/results/primary_only_dedup_removed_category_counts.tsv" \
    --output-dir "$FULLHTLV/results/targeted_method_validation_primary_only" \
    --positive-threshold "${HTLV_THRESHOLD:-100}"

  python "$SCRIPTS/summarize_wgs_retro_results.py" \
    --counts "$WGSFULL/results/primary_only_filtered_category_counts.tsv" \
    --manifest "$OUT/wgs_complete_manifest.csv" \
    --output-dir "$WGSFULL/results/final_summary_primary_only"

  python "$SCRIPTS/audit_filtered_alignments.py" \
    --sample-summary "$WGSFULL/results/final_summary_primary_only/wgs_retro_sample_summary.tsv" \
    --bam-dir "$WGSFULL/bam" \
    --reference-map "$REFERENCE_MAP" \
    --output "$WGSFULL/results/final_summary_primary_only/wgs_exogenous_call_audit.tsv" \
    --category HIV1 --category HIV2 --category HTLV1 --category HTLV2 \
    --min-mapq 40 \
    --min-aligned-length 60 \
    --dedup-mode coordinate \
    --require-unique-best \
    --exclude-secondary-supplementary
}

step_mask_reference() {
  require_file "$REFERENCE_FASTA"
  require_file "$REFERENCE_MAP"
  mkdir -p "$MASKDIR"

  python "$SCRIPTS/mask_shared_retro_regions.py" \
    --input-fasta "$REFERENCE_FASTA" \
    --reference-map "$REFERENCE_MAP" \
    --output-fasta "$MASKED_REFERENCE_FASTA" \
    --output-map "$MASKED_REFERENCE_MAP" \
    --mask-bed "$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.mask.bed" \
    --summary-tsv "$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.mask_summary.tsv" \
    --pair-summary-tsv "$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.pair_similarity.tsv" \
    --report-md "$MASKDIR/hiv1_htlv1_vs_herv.k${MASK_KMER}.mask_report.md" \
    --mask-category HIV1 \
    --mask-category HTLV1 \
    --against-category HERV \
    --kmer-size "$MASK_KMER"

  bwa index "$MASKED_REFERENCE_FASTA"
}

step_masked_rerun() {
  require_file "$MASKED_REFERENCE_FASTA"
  require_file "$MASKED_REFERENCE_MAP"
  require_file "$OUT/targeted_htlv_complete_manifest.csv"
  require_file "$OUT/wgs_complete_manifest.csv"
  mkdir -p "$FULLHTLV_MASKED" "$WGSFULL_MASKED" "$SORTTMP"

  python "$SCRIPTS/run_retro_pilot_alignment.py" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --work-dir "$FULLHTLV_MASKED" \
    --reference-fasta "$MASKED_REFERENCE_FASTA" \
    --reference-map "$MASKED_REFERENCE_MAP" \
    --full-input \
    --jobs "${JOBS:-8}" \
    --threads "${THREADS:-8}" \
    --sort-threads "${SORT_THREADS:-2}" \
    --sort-tmp-dir "$SORTTMP" \
    --resume \
    --min-mapq 40 \
    --category-min-mapq HUMAN:60 \
    --min-aligned-length 60 \
    --dedup-mode coordinate \
    --require-unique-best \
    --exclude-secondary-supplementary

  python "$SCRIPTS/run_retro_pilot_alignment.py" \
    --manifest "$OUT/wgs_complete_manifest.csv" \
    --work-dir "$WGSFULL_MASKED" \
    --reference-fasta "$MASKED_REFERENCE_FASTA" \
    --reference-map "$MASKED_REFERENCE_MAP" \
    --full-input \
    --jobs "${JOBS:-8}" \
    --threads "${THREADS:-8}" \
    --sort-threads "${SORT_THREADS:-2}" \
    --sort-tmp-dir "$SORTTMP" \
    --resume \
    --min-mapq 40 \
    --category-min-mapq HUMAN:60 \
    --min-aligned-length 60 \
    --dedup-mode coordinate \
    --require-unique-best \
    --exclude-secondary-supplementary

  python "$SCRIPTS/summarize_targeted_htlv_results.py" \
    --counts "$FULLHTLV_MASKED/results/filtered_category_counts.tsv" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --sample-qc "$OUT/manifest_qc/sample_file_qc.csv" \
    --output-dir "$FULLHTLV_MASKED/results/final_summary" \
    --positive-threshold "${HTLV_THRESHOLD:-100}"

  python "$SCRIPTS/summarize_wgs_retro_results.py" \
    --counts "$WGSFULL_MASKED/results/filtered_category_counts.tsv" \
    --manifest "$OUT/wgs_complete_manifest.csv" \
    --output-dir "$WGSFULL_MASKED/results/final_summary"
}

usage() {
  cat <<'EOF'
Usage: bash run_p1_final_validation.sh STEP

Steps:
  targeted-hiv  Build manifest and align targeted HIV data.
                Requires TARGETED_HIV_DIR.
  finalize      Regenerate current targeted/WGS summaries and figures.
  primary-only  Recount existing BAMs after removing secondary/supplementary alignments.
  mask-reference Calculate HIV1/HTLV1-vs-HERV shared-kmer masks and index the masked FASTA.
  masked-rerun   Rerun targeted HTLV and WGS against the masked reference.
  compare       Compare old viral-only and new hg38-inclusive runs.
  all           Run targeted-hiv, finalize, primary-only, and compare.

Optional resource variables: JOBS, THREADS, SORT_THREADS, SORTTMP.
Operational threshold: HTLV_THRESHOLD (default 100).
EOF
}

cd "$PROJECT"
case "$STEP" in
  targeted-hiv) step_targeted_hiv ;;
  finalize) step_finalize ;;
  primary-only) step_primary_only ;;
  mask-reference) step_mask_reference ;;
  masked-rerun) step_masked_rerun ;;
  compare) step_compare ;;
  all) step_targeted_hiv; step_finalize; step_primary_only; step_compare ;;
  help|-h|--help) usage ;;
  *) echo "Unknown step: $STEP" >&2; usage >&2; exit 2 ;;
esac
