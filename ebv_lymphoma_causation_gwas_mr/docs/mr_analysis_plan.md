# Mendelian Randomization Analysis Plan

## Nature 2026 Methods To Reproduce

The Nature 2026 virome paper used four Mendelian randomization methods:

- Weighted median
- Inverse-variance weighted (IVW)
- MR-Egger
- Contamination mixture (ConMix)

The primary MR instrument set should exclude the MHC region to reduce pleiotropy risk from linkage disequilibrium.

## Primary Exposure

- EBV DNA load or a validated EBV-load proxy from WGS.

## Candidate Outcomes

- Hodgkin lymphoma.
- Diffuse large B-cell lymphoma.
- PTLD, if usable outcome GWAS summary statistics exist.

## Sensitivity Checks

- MR-Egger intercept for directional pleiotropy.
- Leave-one-instrument-out analysis if enough instruments exist.
- Alternative instrument significance thresholds if the primary set is small.
- Separate analysis excluding known immune-response loci if pleiotropy is suspected.

## Interpretation Categories

- Supports EBV causation.
- Supports immune-state marker.
- Inconclusive.
- Not feasible with available instruments or outcome statistics.
