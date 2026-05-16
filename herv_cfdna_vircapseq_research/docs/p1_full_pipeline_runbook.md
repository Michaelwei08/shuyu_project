# P1 Full Pipeline Runbook

## Goal

Separate confident exogenous HIV/HTLV signal from endogenous HERV/LINE1 and ambiguous retroviral background in targeted sequencing and WGS.

## Inputs

- Sample manifest with `sample_id`, `subject_id`, `cohort`, `assay_type`, `expected_group`, infection status, and file paths.
- Competitive alignment hit table with `sample_id`, `read_id`, `reference_category`, `mapq`, `alignment_score`, and `alignment_length`.
- Usable read-pair depth table.
- Optional viral coverage table.
- Optional host-virus junction table.

## Dry-Run Command Sequence

```powershell
python scripts/bin_retroviral_hits.py data/manifests/example_full_competitive_hits.csv results/example_full_retroviral_bins.csv
python scripts/normalize_retroviral_bins.py results/example_full_retroviral_bins.csv results/example_full_retroviral_bins_normalized.csv --depth-csv data/manifests/example_full_depths.csv
python scripts/summarize_viral_coverage.py data/manifests/example_full_viral_coverage.csv results/example_full_viral_coverage_summary.csv
python scripts/summarize_junction_reads.py data/manifests/example_full_junction_reads.csv results/example_full_junction_summary.csv
python scripts/calibrate_control_thresholds.py data/manifests/example_full_sample_manifest.csv results/example_full_retroviral_bins_normalized.csv results/example_full_frozen_thresholds.json results/example_full_control_noise_before_after.csv
python scripts/call_retroviral_status.py data/manifests/example_full_sample_manifest.csv results/example_full_retroviral_bins_normalized.csv results/example_full_viral_coverage_summary.csv results/example_full_junction_summary.csv results/example_full_frozen_thresholds.json results/example_full_frozen_calls.csv
python scripts/inspect_control_false_signal.py data/manifests/example_full_sample_manifest.csv data/manifests/example_full_competitive_hits.csv results/example_full_control_false_signal.csv
python scripts/summarize_positive_sensitivity.py results/example_full_frozen_calls.csv results/example_full_sensitivity_summary.csv
python scripts/inspect_low_signal_positives.py data/manifests/example_full_sample_manifest.csv results/example_full_frozen_calls.csv data/manifests/example_full_competitive_hits.csv results/example_full_low_signal_positive_review.csv
python scripts/compare_assays_and_wgs.py results/example_full_frozen_calls.csv results/example_full_assay_comparison.csv
```

## Interpretation Rules

- Healthy controls define the first specificity threshold.
- HIV and HTLV positives are evaluated only after thresholds are frozen.
- Reads close to HERV/LINE1 are not promoted to HIV/HTLV calls.
- WGS and targeted sequencing are compared separately because their sensitivity and background differ.
- Dry-run example outputs are not biological evidence; replace all `example_full_*` files with Shuyu/Takeshi/local WGS data before interpretation.

## Current Dry-Run Behavior

- Healthy control has zero frozen HIV/HTLV positive calls.
- Targeted HIV and targeted HTLV examples remain positive after filtering.
- HIV WGS overlap is below the dry-run threshold and is flagged for review.
- HL WGS control remains negative for HIV/HTLV and retains endogenous HERV/LINE1 background.
