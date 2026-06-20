#!/usr/bin/env bash
set -euo pipefail

STEP="${1:-help}"
PROJECT="${PROJECT:-/home/alizadehlab/cpwei/shuyu_project}"
PACKAGE="$PROJECT/herv_cfdna_vircapseq_research/shuyu_benchmark_package"
SCRIPTS="$PACKAGE/scripts"
OUT="${OUT:-$PACKAGE/output}"

BASE_REFDIR="${BASE_REFDIR:-/drive3/cpwei/shuyu_runs/retro_reference_hg38_refseq/ref}"
BASE_REFERENCE_FASTA="${BASE_REFERENCE_FASTA:-$BASE_REFDIR/hg38_plus_retro.refseq.fa}"
BASE_REFERENCE_MAP="${BASE_REFERENCE_MAP:-$BASE_REFDIR/hg38_plus_retro.refseq.reference_map.csv}"
SHUYU_PANEL_FASTA="${SHUYU_PANEL_FASTA:-/drive3/shuyu/references/HIV1_masked/viral_sel_v1_MASKED_HIV1masked.fa}"
SHUYU_PANEL_REFDIR="${SHUYU_PANEL_REFDIR:-/drive3/cpwei/shuyu_runs/shuyu_masked_panel_hg38_herv_line1/ref}"
SHUYU_PANEL_REFERENCE_FASTA="${SHUYU_PANEL_REFERENCE_FASTA:-$SHUYU_PANEL_REFDIR/hg38_herv_line1_plus_shuyu_masked_panel.fa}"
SHUYU_PANEL_REFERENCE_MAP="${SHUYU_PANEL_REFERENCE_MAP:-$SHUYU_PANEL_REFDIR/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv}"
SHUYU_PANEL_INVENTORY="${SHUYU_PANEL_INVENTORY:-$SHUYU_PANEL_REFDIR/shuyu_masked_panel_inventory.csv}"
SORTTMP="${SORTTMP:-/drive3/cpwei/tmp/samtools_sort}"

FULLHTLV_CURRENT="${FULLHTLV_CURRENT:-/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord}"
WGSFULL_CURRENT="${WGSFULL_CURRENT:-/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord}"
FULLHTLV_PANEL="${FULLHTLV_PANEL:-/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_shuyu_masked_panel_primary_only}"
WGSFULL_PANEL="${WGSFULL_PANEL:-/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_shuyu_masked_panel_primary_only}"
PANEL_FILTER_ARGS=(--filter-category HERV --filter-category HIV1 --filter-category HIV2 --filter-category HTLV1 --filter-category HTLV2 --filter-category LINE1 --filter-category OTHER_VIRAL)

require_file() {
  [[ -f "$1" ]] || { echo "Required file not found: $1" >&2; exit 2; }
}

run_alignment() {
  local manifest="$1"
  local work_dir="$2"
  python "$SCRIPTS/run_retro_pilot_alignment.py" \
    --manifest "$manifest" \
    --work-dir "$work_dir" \
    --reference-fasta "$SHUYU_PANEL_REFERENCE_FASTA" \
    --reference-map "$SHUYU_PANEL_REFERENCE_MAP" \
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
    "${PANEL_FILTER_ARGS[@]}"
}

step_build_reference() {
  require_file "$BASE_REFERENCE_FASTA"
  require_file "$BASE_REFERENCE_MAP"
  require_file "$SHUYU_PANEL_FASTA"
  mkdir -p "$SHUYU_PANEL_REFDIR"
  python "$SCRIPTS/make_shuyu_masked_panel_reference.py" \
    --base-fasta "$BASE_REFERENCE_FASTA" \
    --base-map "$BASE_REFERENCE_MAP" \
    --panel-fasta "$SHUYU_PANEL_FASTA" \
    --output-fasta "$SHUYU_PANEL_REFERENCE_FASTA" \
    --output-map "$SHUYU_PANEL_REFERENCE_MAP" \
    --panel-inventory "$SHUYU_PANEL_INVENTORY" \
    --index
}

step_rerun() {
  require_file "$SHUYU_PANEL_REFERENCE_FASTA"
  require_file "$SHUYU_PANEL_REFERENCE_MAP"
  require_file "$OUT/targeted_htlv_complete_manifest.csv"
  require_file "$OUT/wgs_complete_manifest.csv"
  mkdir -p "$FULLHTLV_PANEL" "$WGSFULL_PANEL" "$SORTTMP"
  run_alignment "$OUT/targeted_htlv_complete_manifest.csv" "$FULLHTLV_PANEL"
  run_alignment "$OUT/wgs_complete_manifest.csv" "$WGSFULL_PANEL"
}

step_summarize() {
  require_file "$FULLHTLV_PANEL/results/filtered_category_counts.tsv"
  require_file "$WGSFULL_PANEL/results/filtered_category_counts.tsv"
  python "$SCRIPTS/summarize_targeted_htlv_results.py" \
    --counts "$FULLHTLV_PANEL/results/filtered_category_counts.tsv" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --sample-qc "$OUT/manifest_qc/sample_file_qc.csv" \
    --output-dir "$FULLHTLV_PANEL/results/final_summary" \
    --positive-threshold "${HTLV_THRESHOLD:-100}"
  python "$SCRIPTS/summarize_wgs_retro_results.py" \
    --counts "$WGSFULL_PANEL/results/filtered_category_counts.tsv" \
    --manifest "$OUT/wgs_complete_manifest.csv" \
    --output-dir "$WGSFULL_PANEL/results/final_summary"
  python "$SCRIPTS/export_igv_bam_paths.py" \
    --targeted-counts "$FULLHTLV_PANEL/results/filtered_category_counts.tsv" \
    --targeted-work-dir "$FULLHTLV_PANEL" \
    --wgs-counts "$WGSFULL_PANEL/results/filtered_category_counts.tsv" \
    --wgs-work-dir "$WGSFULL_PANEL" \
    --output "$PACKAGE/output/igv_bam_paths/shuyu_panel_all_bams.tsv"
}

step_export_current_igv() {
  require_file "$FULLHTLV_CURRENT/results/primary_only_filtered_category_counts.tsv"
  require_file "$WGSFULL_CURRENT/results/primary_only_filtered_category_counts.tsv"
  python "$SCRIPTS/export_igv_bam_paths.py" \
    --targeted-counts "$FULLHTLV_CURRENT/results/primary_only_filtered_category_counts.tsv" \
    --targeted-work-dir "$FULLHTLV_CURRENT" \
    --wgs-counts "$WGSFULL_CURRENT/results/primary_only_filtered_category_counts.tsv" \
    --wgs-work-dir "$WGSFULL_CURRENT" \
    --output "$PACKAGE/output/igv_bam_paths/current_primary_only_all_bams.tsv"
  python "$SCRIPTS/export_igv_bam_paths.py" \
    --targeted-counts "$FULLHTLV_CURRENT/results/primary_only_filtered_category_counts.tsv" \
    --targeted-work-dir "$FULLHTLV_CURRENT" \
    --wgs-counts "$WGSFULL_CURRENT/results/primary_only_filtered_category_counts.tsv" \
    --wgs-work-dir "$WGSFULL_CURRENT" \
    --output "$PACKAGE/output/igv_bam_paths/current_primary_only_nonzero_exogenous_bams.tsv" \
    --nonzero-exogenous-only
}

usage() {
  cat <<'EOF'
Usage: bash run_shuyu_masked_panel_validation.sh STEP

Steps:
  build-reference     Build hg38/HERV/LINE1 + Shuyu masked 180-virus panel reference.
  rerun               Rerun targeted HTLV and WGS with primary-only filtering.
  summarize           Summarize Shuyu-panel rerun and export IGV BAM paths.
  export-current-igv  Export current hg38 rerun BAM paths for Shuyu IGV review.
  all                 Run build-reference, rerun, summarize.

Key override variables:
  SHUYU_PANEL_FASTA, JOBS, THREADS, SORT_THREADS, SORTTMP
EOF
}

cd "$PROJECT"
case "$STEP" in
  build-reference) step_build_reference ;;
  rerun) step_rerun ;;
  summarize) step_summarize ;;
  export-current-igv) step_export_current_igv ;;
  all) step_build_reference; step_rerun; step_summarize ;;
  help|-h|--help) usage ;;
  *) echo "Unknown step: $STEP" >&2; usage >&2; exit 2 ;;
esac
