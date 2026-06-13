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

require_file() {
  [[ -f "$1" ]] || { echo "Required file not found: $1" >&2; exit 2; }
}

require_dir() {
  [[ -d "$1" ]] || { echo "Required directory not found: $1" >&2; exit 2; }
}

require_env() {
  [[ -n "${!1:-}" ]] || { echo "Required environment variable is unset: $1" >&2; exit 2; }
}

step_labels() {
  require_env HTLV_LABELS
  require_file "$HTLV_LABELS"
  require_file "$FULLHTLV/results/filtered_category_counts.tsv"
  python "$SCRIPTS/evaluate_known_status_performance.py" \
    --counts "$FULLHTLV/results/filtered_category_counts.tsv" \
    --labels "$HTLV_LABELS" \
    --category HTLV1 \
    --threshold "${HTLV_THRESHOLD:-100}" \
    --output-dir "$FULLHTLV/results/known_status_performance"
}

step_label_template() {
  require_file "$OUT/targeted_htlv_complete_manifest.csv"
  LABEL_TEMPLATE="${LABEL_TEMPLATE:-$OUT/targeted_htlv_status_to_fill.tsv}"
  python "$SCRIPTS/create_status_label_template.py" \
    --manifest "$OUT/targeted_htlv_complete_manifest.csv" \
    --output "$LABEL_TEMPLATE"
  echo "Fill the status column, store the completed file outside Git, then run the labels step."
}

step_targeted_hiv() {
  require_env TARGETED_HIV_DIR
  require_dir "$TARGETED_HIV_DIR"
  require_file "$REFERENCE_FASTA"
  require_file "$REFERENCE_MAP"
  mkdir -p "$OUT" "$SORTTMP" "$TARGETED_HIV_WORK"

  python "$SCRIPTS/prepare_targeted_hiv_manifest.py" \
    --data-dir "$TARGETED_HIV_DIR" \
    --output "$OUT/targeted_hiv_complete_manifest.csv" \
    --issues-output "$OUT/targeted_hiv_pairing_issues.csv" \
    --hiv-status unknown

  python "$SCRIPTS/create_status_label_template.py" \
    --manifest "$OUT/targeted_hiv_complete_manifest.csv" \
    --output "$OUT/targeted_hiv_status_to_fill.tsv"

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

  if [[ -n "${HIV_LABELS:-}" ]]; then
    require_file "$HIV_LABELS"
    python "$SCRIPTS/evaluate_known_status_performance.py" \
      --counts "$TARGETED_HIV_WORK/results/filtered_category_counts.tsv" \
      --labels "$HIV_LABELS" \
      --category HIV1 \
      --threshold "${HIV_THRESHOLD:-100}" \
      --output-dir "$TARGETED_HIV_WORK/results/known_status_performance"
  else
    echo "HIV_LABELS is unset; alignment completed but performance evaluation was skipped."
    echo "Fill $OUT/targeted_hiv_status_to_fill.tsv and rerun with HIV_LABELS set."
  fi
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

usage() {
  cat <<'EOF'
Usage: bash run_p1_final_validation.sh STEP

Steps:
  label-template Create a 71-sample HTLV status sheet for Shuyu to fill.
  labels        Evaluate HTLV1 against HTLV_LABELS (sample_id,status TSV/CSV).
  targeted-hiv  Build manifest and align targeted HIV data.
                Requires TARGETED_HIV_DIR; evaluates performance when HIV_LABELS is set.
  finalize      Regenerate current targeted/WGS summaries and figures.
  compare       Compare old viral-only and new hg38-inclusive runs.
  all           Run labels, targeted-hiv, finalize, and compare.

Optional resource variables: JOBS, THREADS, SORT_THREADS, SORTTMP.
Operational thresholds: HTLV_THRESHOLD (default 100), HIV_THRESHOLD (default 100).
EOF
}

cd "$PROJECT"
case "$STEP" in
  label-template) step_label_template ;;
  labels) step_labels ;;
  targeted-hiv) step_targeted_hiv ;;
  finalize) step_finalize ;;
  compare) step_compare ;;
  all) step_labels; step_targeted_hiv; step_finalize; step_compare ;;
  help|-h|--help) usage ;;
  *) echo "Unknown step: $STEP" >&2; usage >&2; exit 2 ;;
esac
