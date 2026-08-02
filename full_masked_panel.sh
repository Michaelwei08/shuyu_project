#!/usr/bin/env bash
set -euo pipefail

cd /home/alizadehlab/cpwei/shuyu_project
# conda activate retro_qc

PROJECT=$PWD
PACKAGE=$PROJECT/herv_cfdna_vircapseq_research/shuyu_benchmark_package
SCRIPTS=$PACKAGE/scripts
OUT=$PACKAGE/output

export SHUYU_PANEL_FASTA=/drive3/shuyu/references/HIV1_masked/viral_sel_v1_MASKED_HIV1masked.fa
export BASE_REFDIR=/drive3/cpwei/shuyu_runs/retro_reference_hg38_refseq/ref
export BASE_REFERENCE_FASTA=$BASE_REFDIR/hg38_plus_retro.refseq.fa
export BASE_REFERENCE_MAP=$BASE_REFDIR/hg38_plus_retro.refseq.reference_map.csv

export SHUYU_PANEL_REFDIR=/drive3/cpwei/shuyu_runs/shuyu_masked_panel_hg38_herv_line1_refixed/ref
export SHUYU_PANEL_REFERENCE_FASTA=$SHUYU_PANEL_REFDIR/hg38_herv_line1_plus_shuyu_masked_panel.fa
export SHUYU_PANEL_REFERENCE_MAP=$SHUYU_PANEL_REFDIR/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv
export SHUYU_PANEL_INVENTORY=$SHUYU_PANEL_REFDIR/shuyu_masked_panel_inventory.csv

export FULLHTLV_CURRENT=/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord
export WGSFULL_CURRENT=/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord

export FULLHTLV_PANEL=/drive3/cpwei/shuyu_runs/targeted_htlv_hg38_shuyu_masked_panel_primary_only
export WGSFULL_PANEL=/drive3/cpwei/shuyu_runs/wgs_hiv_hl_hg38_shuyu_masked_panel_primary_only

export SORTTMP=/drive3/cpwei/tmp/samtools_sort
export JOBS=4
export THREADS=16
export SORT_THREADS=2
export HTLV_THRESHOLD=100

LOGDIR=/drive3/cpwei/shuyu_runs/logs
mkdir -p "$LOGDIR" "$SORTTMP"
LOG=$LOGDIR/shuyu_masked_panel_full_$(date +%Y%m%d_%H%M%S).log
exec > >(tee -a "$LOG") 2>&1

echo "Started: $(date -Is)"
echo "Log: $LOG"

test -f "$SHUYU_PANEL_FASTA"
test -f "$BASE_REFERENCE_FASTA"
test -f "$BASE_REFERENCE_MAP"
test -f "$OUT/targeted_htlv_complete_manifest.csv"
test -f "$OUT/wgs_complete_manifest.csv"

python -B -m py_compile \
  "$SCRIPTS/make_shuyu_masked_panel_reference.py" \
  "$SCRIPTS/run_retro_pilot_alignment.py" \
  "$SCRIPTS/summarize_targeted_htlv_results.py" \
  "$SCRIPTS/summarize_wgs_retro_results.py" \
  "$SCRIPTS/export_igv_bam_paths.py"

if [[ ! -s "$SHUYU_PANEL_REFERENCE_FASTA.bwt" ]]; then
  bash "$SCRIPTS/run_shuyu_masked_panel_validation.sh" build-reference
else
  echo "Reusing existing BWA index: $SHUYU_PANEL_REFERENCE_FASTA"
fi

bash "$SCRIPTS/run_shuyu_masked_panel_validation.sh" rerun
bash "$SCRIPTS/run_shuyu_masked_panel_validation.sh" summarize

python "$SCRIPTS/export_igv_bam_paths.py" \
  --targeted-counts "$FULLHTLV_PANEL/results/filtered_category_counts.tsv" \
  --targeted-work-dir "$FULLHTLV_PANEL" \
  --wgs-counts "$WGSFULL_PANEL/results/filtered_category_counts.tsv" \
  --wgs-work-dir "$WGSFULL_PANEL" \
  --output "$PACKAGE/output/igv_bam_paths/shuyu_panel_nonzero_exogenous_bams.tsv" \
  --nonzero-exogenous-only

python - <<'PY'
from pathlib import Path
import csv, os

pairs = [
    ("targeted_htlv", os.environ["FULLHTLV_CURRENT"] + "/results/primary_only_filtered_category_counts.tsv",
     os.environ["FULLHTLV_PANEL"] + "/results/filtered_category_counts.tsv"),
    ("wgs_hiv_hl", os.environ["WGSFULL_CURRENT"] + "/results/primary_only_filtered_category_counts.tsv",
     os.environ["WGSFULL_PANEL"] + "/results/filtered_category_counts.tsv"),
]
out = Path(os.environ["PACKAGE"]) / "output" / "shuyu_panel_vs_current_counts.tsv"
species = ["HERV","HIV1","HIV2","HTLV1","HTLV2","LINE1","OTHER_VIRAL"]

rows = []
for cohort, current_path, panel_path in pairs:
    current = {r["sample"]: r for r in csv.DictReader(open(current_path), delimiter="\t")}
    panel = {r["sample"]: r for r in csv.DictReader(open(panel_path), delimiter="\t")}
    for sample in sorted(set(current) | set(panel)):
        row = {"cohort": cohort, "sample": sample}
        for sp in species:
            row[f"current_{sp}"] = int(float(current.get(sample, {}).get(sp, 0) or 0))
            row[f"panel_{sp}"] = int(float(panel.get(sample, {}).get(sp, 0) or 0))
        rows.append(row)

cols = ["cohort","sample"] + [f"{prefix}_{sp}" for sp in species for prefix in ("current","panel")]
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(out)
PY

echo "Targeted panel report:"
cat "$FULLHTLV_PANEL/results/final_summary/targeted_htlv_full_report.md"

echo "WGS panel report:"
cat "$WGSFULL_PANEL/results/final_summary/wgs_retro_report.md"

echo "Completion check:"
echo -n "Targeted idxstats: "
find "$FULLHTLV_PANEL/results" -name "*.idxstats.tsv" | wc -l
echo -n "WGS idxstats: "
find "$WGSFULL_PANEL/results" -name "*.idxstats.tsv" | wc -l

echo "Key outputs:"
ls -lh "$FULLHTLV_PANEL/results/final_summary"
ls -lh "$WGSFULL_PANEL/results/final_summary"
ls -lh "$PACKAGE/output/igv_bam_paths/"*shuyu_panel*".tsv"
ls -lh "$PACKAGE/output/shuyu_panel_vs_current_counts.tsv"

echo "Finished: $(date -Is)"
