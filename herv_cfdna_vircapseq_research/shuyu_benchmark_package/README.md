# Shuyu HIV/HTLV/HERV Benchmark Package

## Purpose

Prepare the benchmark inputs for Shuyu's shared data without copying protected sequencing data into this workspace.

Expected remote paths:

- `/drive3/Shuyu_for_Michael/WGS_60samples_HIV_HL`
- `/drive3/Shuyu_for_Michael/Targeted_HTLV_TCLsamples`

The package creates:

- a full sample manifest inferred from filenames,
- a small pilot manifest,
- a data-access check table,
- a benchmark plan,
- command templates for the existing P1 retrovirus/HERV scripts.

## Quick Start On The Machine With `/drive3`

Run from the repo root:

```bash
python herv_cfdna_vircapseq_research/shuyu_benchmark_package/scripts/prepare_shuyu_benchmark.py \
  --wgs-dir /drive3/Shuyu_for_Michael/WGS_60samples_HIV_HL \
  --htlv-targeted-dir /drive3/Shuyu_for_Michael/Targeted_HTLV_TCLsamples \
  --output-dir herv_cfdna_vircapseq_research/shuyu_benchmark_package/output
```

Then validate the generated manifest:

```bash
python herv_cfdna_vircapseq_research/scripts/validate_sample_manifest.py \
  herv_cfdna_vircapseq_research/shuyu_benchmark_package/output/sample_manifest.csv
```

## Interpretation Of Inferred Labels

- WGS files starting with `HIV`: expected HIV-infected patient group. HIV reads may be present at low level.
- WGS files starting with `HL`: expected HIV/HTLV-negative control group. HIV or HTLV calls should be treated as false positive, cross-mapping, or HERV/retroelement signal unless independently validated.
- Targeted TCL HTLV files: HTLV-enriched group. HTLV should be detectable in many but not necessarily all samples.

## What This Package Does Not Do

- It does not download, copy, or modify raw sequencing data.
- It does not run alignment because reference paths and compute environment are not yet known.
- It does not claim clinical association evidence.

## Next Analytical Step

After the manifest is generated, run a small pilot first:

- 2 HIV-prefixed WGS samples.
- 2 HL-prefixed WGS controls.
- 3 targeted HTLV TCL samples.

The initial success criterion is specificity-first:

- HL WGS controls should have near-zero HIV and HTLV signal after filtering.
- HTLV targeted samples should show strong HTLV signal in many samples.
- HIV WGS samples may have low but plausible HIV signal.
