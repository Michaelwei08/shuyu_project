# HIV / HTLV / HERV Benchmark Memo

## Objective

Benchmark retrovirus-aware detection on targeted sequencing and WGS using known HIV-positive cases, Takeshi's HTLV cohort, healthy controls, and the 60 local WGS samples.

## Current Inputs

- Sample manifest: pending Shuyu.
- Sequencing files: pending Shuyu.
- Combined reference plan: complete.
- First-pass binning scripts: complete.
- Normalization and coverage scripts: complete.

## Planned Read Bins

- `unique_hiv`
- `unique_htlv`
- `unique_herv_line1`
- `ambiguous_retroviral`
- `other_viral`
- `human`

## Primary Benchmark Questions

- Do healthy controls have near-zero high-confidence HIV and HTLV signal after filtering?
- Do known HIV-positive samples retain HIV signal?
- Do known or suspected HTLV-positive samples show HTLV signal?
- How much signal remains ambiguous after HERV and LINE1 decoy competition?

## Failure Modes To Track

- Conserved retroviral sequence misassigned to HIV or HTLV.
- HERV or LINE1 background inflated by short or low-complexity reads.
- Targeted-capture off-target bias.
- Reference mismatch for HIV subtype or HTLV subtype.
- Duplicate-driven false signal.
- Low WGS depth preventing locus-level HERV interpretation.

## Next Improvements

- Populate real sample manifest.
- Run healthy controls first to calibrate specificity.
- Freeze thresholds before final positive-cohort sensitivity analysis.
- Add subtype-specific HIV/HTLV references if sensitivity is poor.
