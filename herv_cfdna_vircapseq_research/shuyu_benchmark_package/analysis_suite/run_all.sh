#!/usr/bin/env bash
#
# run_all.sh -- run the whole analysis suite (a1..a7) into one --outdir on the
# cluster, then render the figures (a8) over those outputs.
#
# Every step gets its documented default input paths. Steps are independent:
# one failing step does not stop the others, and the failures are listed at the
# end. Exit status is 0 only when every step succeeded.
#
# a9 and a10 are OPTIONAL follow-ups and are not part of the default run. Each is
# gated on an input the default run does not have, and is reported as SKIPPED --
# never as a failure -- when that input is absent:
#   a9   needs a clinical CD4 table   -> set CD4_TABLE=<file>
#   a10  needs BAMs                   -> needs <A10_RUN>/bam/*.bam to exist
#
# USAGE
#   bash run_all.sh                       # documented defaults (a1..a8)
#   bash run_all.sh /path/to/outdir       # same, different output directory
#   OUTDIR=... RUNS_ROOT=... bash run_all.sh
#   CD4_TABLE=/path/to/cd4.tsv bash run_all.sh          # also runs a9
#
# ENVIRONMENT OVERRIDES (all optional)
#   OUTDIR      where every .tsv lands            (default below)
#   RUNS_ROOT   root holding the run directories  (/path/to/runs)
#   PYTHON      python interpreter                (python3)
#   SAMTOOLS    samtools executable               (samtools)
#   CD4_TABLE   clinical CD4 table for a9         (unset -> a9 is SKIPPED)
#   A10_RUN     run directory a10 audits          (default: the WGS panel run;
#                                                  no bam/*.bam -> a10 SKIPPED)
#   ONLY        space-separated step ids to run    e.g. ONLY="a1 a8"
#
# IDENTIFIERS
#   Only the *_sample_key.tsv files carry real sample names. Every other output,
#   including this script's per-step logs, uses the anonymous ids S01..Snn.
#   Do not commit or email the key files or the logs directory without checking.
#
# Date: 2026-07-26
set -euo pipefail

# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #
RUNS_ROOT="${RUNS_ROOT:-/path/to/runs}"
OUTDIR="${OUTDIR:-${1:-${RUNS_ROOT}/panel_report_20260725/suite_out}}"
SAMTOOLS="${SAMTOOLS:-samtools}"
ONLY="${ONLY:-}"

# --------------------------------------------------------------------------- #
# interpreter selection
#
# The modules need Python >= 3.7 ("from __future__ import annotations"). Several
# clusters still ship a 3.6 as plain "python3", which fails at import time in
# every module at once. If PYTHON is set we honour it but still check it; if not,
# we probe the usual candidates and take the first that is new enough.
# --------------------------------------------------------------------------- #
py_ok () {  # $1 = interpreter; true if it exists and is >= 3.7
    command -v "$1" >/dev/null 2>&1 || return 1
    "$1" -c 'import sys; sys.exit(0 if sys.version_info[:2] >= (3, 7) else 1)' >/dev/null 2>&1
}

if [ -n "${PYTHON:-}" ]; then
    if ! py_ok "$PYTHON"; then
        echo "ERROR: PYTHON=$PYTHON is missing or older than 3.7." >&2
        "$PYTHON" --version >&2 2>/dev/null || true
        echo "       The modules need >= 3.7 (from __future__ import annotations)." >&2
        exit 2
    fi
else
    PYTHON=""
    for _cand in python3.12 python3.11 python3.10 python3.9 python3.8 python3 python; do
        if py_ok "$_cand"; then PYTHON="$_cand"; break; fi
    done
    if [ -z "$PYTHON" ]; then
        echo "ERROR: no Python >= 3.7 found on PATH (tried python3.12 .. python)." >&2
        echo "       Set one explicitly, e.g. PYTHON=/usr/bin/python3.11 bash run_all.sh" >&2
        exit 2
    fi
fi
export PYTHON

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIGDIR="${OUTDIR}/figures"
LOGDIR="${OUTDIR}/logs"

# documented default inputs (see each module's --help)
TARGETED_CURRENT="${RUNS_ROOT}/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord"
WGS_CURRENT="${RUNS_ROOT}/wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord"
WGS_PANEL="${RUNS_ROOT}/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only"
TARGETED_PANEL="${RUNS_ROOT}/targeted_htlv_hg38_shuyu_masked_panel_refixed_primary_only"
PANEL_REFMAP="${RUNS_ROOT}/shuyu_masked_panel_hg38_herv_line1_refixed/ref/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv"
BASE_REFMAP="${RUNS_ROOT}/retro_reference_hg38_refseq/ref/hg38_plus_retro.refseq.reference_map.csv"
LADDER_SUMMARY="${RUNS_ROOT}/reply_to_shuyu_primary_only/kmer_ladder_summary.tsv"
LADDER_DIR="${RUNS_ROOT}/reply_to_shuyu_primary_only/kmer_ladder"
MASK_METRICS_DIR="${RUNS_ROOT}/retro_reference_hg38_refseq_mask_metrics_k40"
MASKED_BUILD_DIR="${RUNS_ROOT}/retro_reference_hg38_refseq_masked_hiv1_htlv1_vs_herv_k40"

DETECTION_THRESHOLD=100

# optional-step inputs (see the gates further down)
CD4_TABLE="${CD4_TABLE:-}"
A10_RUN="${A10_RUN:-$WGS_PANEL}"

if ! mkdir -p "$OUTDIR" "$FIGDIR" "$LOGDIR"; then
    echo "FATAL: cannot create the output directory ${OUTDIR}" >&2
    exit 2
fi

STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S')"
OK_STEPS=()
FAILED_STEPS=()
SKIPPED_STEPS=()

# --------------------------------------------------------------------------- #
# step runner
# --------------------------------------------------------------------------- #
selected() {
    [ -z "$ONLY" ] && return 0
    local want
    for want in $ONLY; do
        [ "$want" = "$1" ] && return 0
    done
    return 1
}

run_step() {
    local id="$1"; shift
    local question="$1"; shift

    if ! selected "$id"; then
        SKIPPED_STEPS+=("${id} (not in ONLY)")
        return 0
    fi

    local log="${LOGDIR}/${id}.log"
    echo ""
    echo "=============================================================================="
    echo "[${id}]  ${question}"
    echo "  module : $(basename "$1")"
    echo "  outdir : ${OUTDIR}"
    echo "  log    : ${log}"
    echo "  started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=============================================================================="

    # pipefail makes this see the module's status, not tee's
    if "$PYTHON" "$@" 2>&1 | tee "$log"; then
        echo "[${id}] OK"
        OK_STEPS+=("$id")
    else
        local rc=$?
        echo "[${id}] FAILED (exit ${rc}) -- continuing; see ${log}"
        FAILED_STEPS+=("${id} (exit ${rc})")
    fi
    return 0
}

# Gated step. $2 is either the literal "ok" or the reason the step cannot run;
# a step that cannot run is SKIPPED, so the default run's exit status does not
# change just because an optional input is missing.
maybe_step() {
    local id="$1"; shift
    local gate="$1"; shift

    if ! selected "$id"; then
        SKIPPED_STEPS+=("${id} (not in ONLY)")
        return 0
    fi
    if [ "$gate" != "ok" ]; then
        SKIPPED_STEPS+=("${id} (${gate})")
        echo ""
        echo "=============================================================================="
        echo "[${id}]  SKIPPED -- ${gate}"
        echo "=============================================================================="
        return 0
    fi
    run_step "$id" "$@"
    return 0
}

echo "=============================================================================="
echo "VIRAL SEQUENCING ANALYSIS SUITE -- run_all.sh"
echo "=============================================================================="
echo "  started    : ${STARTED_AT}"
echo "  suite dir  : ${HERE}"
echo "  runs root  : ${RUNS_ROOT}"
echo "  outdir     : ${OUTDIR}"
echo "  figures    : ${FIGDIR}"
echo "  logs       : ${LOGDIR}"
echo "  python     : ${PYTHON}  ($("${PYTHON}" --version 2>&1))"
echo "  samtools   : ${SAMTOOLS}"
if [ -n "$ONLY" ]; then
    echo "  ONLY       : ${ONLY}"
fi
if [ ! -d "$RUNS_ROOT" ]; then
    echo ""
    echo "  NOTE: runs root ${RUNS_ROOT} does not exist here. Every module reports"
    echo "        its missing inputs with a WARN and exits 0, so this run will"
    echo "        produce empty or header-only tables rather than failing."
fi

# --------------------------------------------------------------------------- #
# a1 .. a7 -- all into one --outdir
# --------------------------------------------------------------------------- #
run_step a1 "How well does the panel detect HTLV-1 / HIV-1?" \
    "${HERE}/a1_detection_performance.py" \
    --runs-root "$RUNS_ROOT" \
    --run "$TARGETED_CURRENT" \
    --run "$WGS_CURRENT" \
    --targets HTLV1,HIV1 \
    --threshold "$DETECTION_THRESHOLD" \
    --outdir "$OUTDIR" \
    --prefix detection \
    --samtools "$SAMTOOLS"

run_step a2 "What does adding the human genome as a competitor change?" \
    "${HERE}/a2_reference_comparison.py" \
    --runs-root "$RUNS_ROOT" \
    --base-refmap "$BASE_REFMAP" \
    --panel-refmap "$PANEL_REFMAP" \
    --outdir "$OUTDIR" \
    --prefix a2 \
    --samtools "$SAMTOOLS"

run_step a3 "What does k cost and buy in the HERV masking ladder?" \
    "${HERE}/a3_kmer_ladder.py" \
    --ladder-summary "$LADDER_SUMMARY" \
    --ladder-dir "$LADDER_DIR" \
    --mask-metrics-dir "$MASK_METRICS_DIR" \
    --masked-build-dir "$MASKED_BUILD_DIR" \
    --outdir "$OUTDIR" \
    --prefix kmer_ladder

run_step a4 "Which detections survive subsampling to 5M reads?" \
    "${HERE}/a4_depth_sensitivity.py" \
    --runs-root "$RUNS_ROOT" \
    --sub-run wgs_hiv_hl_5m_competitive \
    --full-run wgs_hiv_hl_full_competitive \
    --outdir "$OUTDIR" \
    --prefix depth_sensitivity \
    --samtools "$SAMTOOLS"

run_step a5 "Is a per-reference signal a real genome or a pile-up?" \
    "${HERE}/a5_reference_depth_profiles.py" \
    --run "$WGS_PANEL" \
    --run "$TARGETED_PANEL" \
    --ref ebv1 --ref ebv2 --ref hhv6b \
    --refmap "$PANEL_REFMAP" \
    --outdir "$OUTDIR" \
    --prefix refprofile \
    --samtools "$SAMTOOLS"

run_step a6 "Is there read-level evidence of HTLV-1 integration?" \
    "${HERE}/a6_htlv_junctions.py" \
    --run-dir "$TARGETED_CURRENT" \
    --refmap "$BASE_REFMAP" \
    --outdir "$OUTDIR" \
    --prefix htlv_junction \
    --samtools "$SAMTOOLS"

run_step a7 "What is the anellovirus burden and coinfection structure?" \
    "${HERE}/a7_virome_structure.py" \
    --runs "$WGS_PANEL" "$TARGETED_PANEL" \
    --refmap "$PANEL_REFMAP" \
    --base-refmap "$BASE_REFMAP" \
    --outdir "$OUTDIR" \
    --prefix a7_virome \
    --samtools "$SAMTOOLS"

# --------------------------------------------------------------------------- #
# a9, a10 -- OPTIONAL follow-ups, gated on inputs the default run does not have.
# Both read a7's outputs out of "$OUTDIR", so they must come after a7.
# --------------------------------------------------------------------------- #
A9_GATE="ok"
if [ -z "$CD4_TABLE" ]; then
    A9_GATE="no CD4 table; set CD4_TABLE=<file> to enable"
elif [ ! -f "$CD4_TABLE" ]; then
    A9_GATE="CD4 table not found at ${CD4_TABLE}"
fi

maybe_step a9 "$A9_GATE" \
    "Does the anellovirus burden rise as CD4 falls?" \
    "${HERE}/a9_cd4_correlation.py" \
    --indir "$OUTDIR" \
    --outdir "$OUTDIR" \
    --cd4 "$CD4_TABLE" \
    --a7-prefix a7_virome

A10_GATE="ok"
if [ ! -d "$A10_RUN" ]; then
    A10_GATE="run directory ${A10_RUN} is not present"
else
    # an unmatched glob stays literal, so [0] always exists and -e settles it
    A10_BAMS=( "${A10_RUN}"/bam/*.bam )
    if [ ! -e "${A10_BAMS[0]}" ]; then
        A10_GATE="no BAM at ${A10_RUN}/bam/*.bam"
    fi
fi

maybe_step a10 "$A10_GATE" \
    "Is the low-count anellovirus signal real virus or cross-mapping?" \
    "${HERE}/a10_anello_read_audit.py" \
    --run "$A10_RUN" \
    --refmap "$PANEL_REFMAP" \
    --indir "$OUTDIR" \
    --outdir "$OUTDIR" \
    --prefix anello_read_audit \
    --samtools "$SAMTOOLS"

# --------------------------------------------------------------------------- #
# a8 -- figures over whatever a1..a7 actually wrote
# --------------------------------------------------------------------------- #
run_step a8 "Render the figures from the tables a1..a7 wrote" \
    "${HERE}/a8_figures.py" \
    --indir "$OUTDIR" \
    --outdir "$FIGDIR" \
    --prefix a8 \
    --default-threshold "$DETECTION_THRESHOLD"

# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
echo ""
echo "=============================================================================="
echo "SUMMARY"
echo "=============================================================================="
echo "  started  : ${STARTED_AT}"
echo "  finished : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  succeeded: ${#OK_STEPS[@]}  ${OK_STEPS[*]:-(none)}"
echo "  failed   : ${#FAILED_STEPS[@]}"
if [ "${#FAILED_STEPS[@]}" -gt 0 ]; then
    for step in "${FAILED_STEPS[@]}"; do
        echo "             ${step}"
    done
fi
if [ "${#SKIPPED_STEPS[@]}" -gt 0 ]; then
    echo "  skipped  : ${#SKIPPED_STEPS[@]}"
    for step in "${SKIPPED_STEPS[@]}"; do
        echo "             ${step}"
    done
fi

echo ""
echo "  tables   : ${OUTDIR}"
echo "  figures  : ${FIGDIR}"
echo "  logs     : ${LOGDIR}"

echo ""
echo "  IDENTIFIERS -- these files carry real sample names, nothing else does:"
found_key=0
for key in "${OUTDIR}"/*_sample_key.tsv "${FIGDIR}"/*_sample_key.tsv; do
    if [ -f "$key" ]; then
        echo "    ${key}"
        found_key=1
    fi
done
if [ "$found_key" -eq 0 ]; then
    echo "    (none written)"
fi
echo "  Do NOT commit or email the *_sample_key.tsv files."
echo "  Per-step logs are anonymised, but skim them before sharing."

if [ "${#FAILED_STEPS[@]}" -gt 0 ]; then
    echo ""
    echo "Finished with ${#FAILED_STEPS[@]} failed step(s)."
    exit 1
fi
echo ""
echo "All selected steps succeeded."
exit 0
