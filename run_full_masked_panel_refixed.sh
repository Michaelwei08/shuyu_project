#!/usr/bin/env bash
set -eo pipefail
umask 077

PROJECT=/home/alizadehlab/cpwei/shuyu_project
PACKAGE=$PROJECT/herv_cfdna_vircapseq_research/shuyu_benchmark_package
SCRIPTS=$PACKAGE/scripts
OUT=$PACKAGE/output

export PROJECT PACKAGE SCRIPTS OUT

export SHUYU_PANEL_FASTA=/drive3/shuyu/references/HIV1_masked/viral_sel_v1_MASKED_HIV1masked.fa
export SHUYU_PANEL_REFDIR=/drive3/cpwei/shuyu_runs/shuyu_masked_panel_hg38_herv_line1_refixed/ref
export SHUYU_PANEL_REFERENCE_FASTA=$SHUYU_PANEL_REFDIR/hg38_herv_line1_plus_shuyu_masked_panel.fa
export SHUYU_PANEL_REFERENCE_MAP=$SHUYU_PANEL_REFDIR/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv
export SHUYU_PANEL_INVENTORY=$SHUYU_PANEL_REFDIR/shuyu_masked_panel_inventory.csv

export FULLHTLV_CURRENT=/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord
export WGSFULL_CURRENT=/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord

export FULLHTLV_PANEL=/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_shuyu_masked_panel_refixed_primary_only
export WGSFULL_PANEL=/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only

export SORTTMP=/drive3/cpwei/tmp/samtools_sort
export JOBS=8
export THREADS=12
export SORT_THREADS=2
export HTLV_THRESHOLD=100

LOGDIR=/drive3/cpwei/shuyu_runs/logs
mkdir -p "$LOGDIR" "$SORTTMP"
LOG=$LOGDIR/shuyu_masked_panel_refixed_full_$(date +%Y%m%d_%H%M%S).log
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date -Is)"
echo "Log: $LOG"

test -f "$SHUYU_PANEL_REFERENCE_FASTA"
test -f "$SHUYU_PANEL_REFERENCE_MAP"
test -f "$SHUYU_PANEL_REFERENCE_FASTA.bwt"
test -f "$OUT/targeted_htlv_complete_manifest.csv"
test -f "$OUT/wgs_complete_manifest.csv"

python - <<'PY'
import csv, collections, os, sys

p = os.environ["SHUYU_PANEL_REFERENCE_MAP"]
c = collections.Counter(r["category"] for r in csv.DictReader(open(p)))
print("Reference category counts:", dict(c))

required = {"HUMAN", "HERV", "LINE1", "HIV1", "HIV2", "HTLV1", "HTLV2"}
missing = sorted(required - set(c))
if missing:
    sys.exit("Missing required categories in reference map: " + ",".join(missing))
PY

python -B -m py_compile \
  "$SCRIPTS/run_retro_pilot_alignment.py" \
  "$SCRIPTS/summarize_targeted_htlv_results.py" \
  "$SCRIPTS/summarize_wgs_retro_results.py" \
  "$SCRIPTS/export_igv_bam_paths.py"

bash "$SCRIPTS/run_shuyu_masked_panel_validation.sh" rerun
bash "$SCRIPTS/run_shuyu_masked_panel_validation.sh" summarize

mkdir -p "$PACKAGE/output/igv_bam_paths"

python "$SCRIPTS/export_igv_bam_paths.py" \
  --targeted-counts "$FULLHTLV_PANEL/results/filtered_category_counts.tsv" \
  --targeted-work-dir "$FULLHTLV_PANEL" \
  --wgs-counts "$WGSFULL_PANEL/results/filtered_category_counts.tsv" \
  --wgs-work-dir "$WGSFULL_PANEL" \
  --output "$PACKAGE/output/igv_bam_paths/shuyu_panel_refixed_all_bams.tsv"

python "$SCRIPTS/export_igv_bam_paths.py" \
  --targeted-counts "$FULLHTLV_PANEL/results/filtered_category_counts.tsv" \
  --targeted-work-dir "$FULLHTLV_PANEL" \
  --wgs-counts "$WGSFULL_PANEL/results/filtered_category_counts.tsv" \
  --wgs-work-dir "$WGSFULL_PANEL" \
  --output "$PACKAGE/output/igv_bam_paths/shuyu_panel_refixed_nonzero_exogenous_bams.tsv" \
  --nonzero-exogenous-only

python - <<'PY'
from pathlib import Path
import csv, os

package = Path(os.environ["PACKAGE"])
pairs = [
    (
        "targeted_htlv",
        Path(os.environ["FULLHTLV_CURRENT"]) / "results" / "primary_only_filtered_category_counts.tsv",
        Path(os.environ["FULLHTLV_PANEL"]) / "results" / "filtered_category_counts.tsv",
    ),
    (
        "wgs_hiv_hl",
        Path(os.environ["WGSFULL_CURRENT"]) / "results" / "primary_only_filtered_category_counts.tsv",
        Path(os.environ["WGSFULL_PANEL"]) / "results" / "filtered_category_counts.tsv",
    ),
]

species = ["HERV", "HIV1", "HIV2", "HTLV1", "HTLV2", "LINE1", "OTHER_VIRAL", "HUMAN"]
rows = []

for cohort, current_path, panel_path in pairs:
    if not current_path.exists():
        raise SystemExit(f"Missing current counts: {current_path}")
    if not panel_path.exists():
        raise SystemExit(f"Missing panel counts: {panel_path}")

    current = {r["sample"]: r for r in csv.DictReader(open(current_path), delimiter="\t")}
    panel = {r["sample"]: r for r in csv.DictReader(open(panel_path), delimiter="\t")}

    for sample in sorted(set(current) | set(panel)):
        row = {"cohort": cohort, "sample": sample}
        for sp in species:
            row[f"current_{sp}"] = int(float(current.get(sample, {}).get(sp, 0) or 0))
            row[f"panel_{sp}"] = int(float(panel.get(sample, {}).get(sp, 0) or 0))
        rows.append(row)

out = package / "output" / "shuyu_panel_refixed_vs_current_counts.tsv"
cols = ["cohort", "sample"] + [f"{prefix}_{sp}" for sp in species for prefix in ("current", "panel")]
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

print("Wrote comparison:", out)
PY

python - <<'PY'
import csv, os
from pathlib import Path

checks = [
    ("targeted", Path(os.environ["OUT"]) / "targeted_htlv_complete_manifest.csv", Path(os.environ["FULLHTLV_PANEL"]) / "results" / "filtered_category_counts.tsv"),
    ("WGS", Path(os.environ["OUT"]) / "wgs_complete_manifest.csv", Path(os.environ["WGSFULL_PANEL"]) / "results" / "filtered_category_counts.tsv"),
]

for name, manifest, counts in checks:
    expected = {r["sample_id"] for r in csv.DictReader(open(manifest))}
    observed = {r["sample"] for r in csv.DictReader(open(counts), delimiter="\t")}
    print(f"{name}: expected={len(expected)} completed={len(observed)} missing={len(expected - observed)} unexpected={len(observed - expected)}")
PY

echo "Targeted report:"
cat "$FULLHTLV_PANEL/results/final_summary/targeted_htlv_full_report.md"

echo "WGS report:"
cat "$WGSFULL_PANEL/results/final_summary/wgs_retro_report.md"

echo "Key outputs:"
ls -lh "$FULLHTLV_PANEL/results/final_summary"
ls -lh "$WGSFULL_PANEL/results/final_summary"
ls -lh "$PACKAGE/output/igv_bam_paths/shuyu_panel_refixed_"*.tsv
ls -lh "$PACKAGE/output/shuyu_panel_refixed_vs_current_counts.tsv"

echo "Finished: $(date -Is)"
