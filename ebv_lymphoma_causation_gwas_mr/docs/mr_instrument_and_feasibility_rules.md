# MR Instrument And Feasibility Rules

## Primary Instrument Selection

- Select candidate EBV-load instruments from exposure GWAS results or validated published exposure summary statistics.
- Exclude extended MHC-region variants from the primary MR instrument set because of pleiotropy and linkage-disequilibrium risk.
- Use genome-wide significance for final instruments where possible.
- If genome-wide significant non-MHC instruments are absent, record the relaxed-threshold exploratory instrument set separately and do not present it as primary causal evidence.

## Clumping

The current scaffold supports physical-window pruning for toy data. Final analysis should use LD clumping against an ancestry-matched reference panel before any causal interpretation.

## Harmonization

- Match exposure and outcome variants by `variant_id`.
- Align outcome effects to the exposure effect allele.
- Flip outcome beta when alleles are reversed.
- Flag palindromic SNPs for allele-frequency review.
- Drop variants with unresolved allele mismatches.

## Method Feasibility

- IVW, weighted median, and MR-Egger are available as lightweight scaffold summaries in `scripts/run_mr_summary.py`.
- MR-Egger intercept is emitted, but final pleiotropy testing should include a standard error and P value from a validated MR implementation.
- Leave-one-instrument-out sensitivity is scaffolded in `scripts/run_mr_sensitivity.py`.
- ConMix is not fully implemented in this workspace. The current output is a placeholder proxy and must not be used for final claims.

## Current Example Outcome

In the toy GWAS example, excluding the MHC leaves one relaxed-threshold non-MHC instrument at `p <= 0.01` and zero non-MHC instruments at `p <= 5e-8`. This demonstrates the pipeline behavior but does not support any biological inference.
