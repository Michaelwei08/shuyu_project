# Phenotype Freeze Record

## Freeze Status

Phenotype definitions are frozen as `v0_scaffold` for dry-run and script-validation work only.

Final biological analyses must create a new freeze record after real viral detection outputs are reviewed. The current freeze does not authorize final association, GWAS, or MR claims from example data.

## Frozen Scaffold Definitions

- `ebv_load`: normalized EBV read-pair count after mapping, duplicate, breadth, and depth filters.
- `anellovirus_load`: normalized aggregate anellovirus read-pair count.
- `other_herpesvirus_load`: normalized non-EBV herpesvirus signal.
- `pan_viral_burden`: aggregate normalized viral signal across included viral families.
- `v1_pan_viral`: broad viral signal across multiple viral families.
- `v3_high_ebv_low_anello`: high EBV with low anellovirus signal.

## Required Real-Data Refreeze

Before final analysis, update this file with:

- Input viral-call table path and checksum or timestamp.
- Final thresholds for EBV, anellovirus, and pan-viral burden.
- Depth and batch adjustment rules.
- Inclusion/exclusion criteria.
- Samples included by cohort and disease group.
