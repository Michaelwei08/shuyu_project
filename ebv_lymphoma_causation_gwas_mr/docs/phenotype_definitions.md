# Phenotype Definitions

## Viral Load Phenotypes

- `ebv_load`: normalized EBV read-pair count, with breadth and mapping-quality filters applied.
- `anellovirus_load`: normalized aggregate anellovirus read-pair count across selected representative references.
- `other_herpesvirus_load`: normalized aggregate read-pair count for non-EBV herpesviruses.
- `pan_viral_burden`: normalized aggregate signal across EBV, non-EBV herpesviruses, anelloviruses, and any other viruses included in the final panel.

## Subtype Phenotypes

- `v1_pan_viral`: subtype with broad viral signal across multiple viral families. Working interpretation: may reflect immune suppression or pan-viral reactivation.
- `v3_high_ebv_low_anello`: subtype with high EBV and low anellovirus signal. Working interpretation: stronger candidate for EBV disease contribution.

## Missingness Rules

- Missing sequencing data is `missing`, not zero.
- Viral reads below quality thresholds are `filtered`, not zero.
- A sample with adequate depth and no confident viral reads can be coded as zero for that viral phenotype.

## Freeze Rule

Phenotype definitions must be frozen before clinical association testing. Any post-hoc threshold changes should create a new phenotype version.

Current status: definitions exist as version `v0`, but they are not frozen for final analysis until the first real viral-detection output table is reviewed.
