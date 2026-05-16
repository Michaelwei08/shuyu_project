# GWAS And MR Methods Note

## Viral Phenotypes

Viral-load phenotypes will be taken from the public/local WGS viral-detection workflow after read-level filtering, duplicate handling, breadth checks, and depth normalization. Primary phenotypes are EBV load, anellovirus load, pan-viral burden, V1 pan-viral subtype, and V3 high-EBV low-anellovirus subtype.

## Association Models

Initial clinical association tests should compare viral phenotypes across HL, PTLD, DLBCL, and controls. Primary adjusted models should include sequencing depth, batch, sample type, cohort, and disease-specific covariates when available. PTLD models should explicitly include transplant and immunosuppression variables when metadata permit.

## GWAS

GWAS should be run only for datasets with viral phenotypes, genotype or imputed genotype data, ancestry principal components, sequencing depth, and core covariates. The local 60 WGS samples are expected to support phenotype validation and association checks, not discovery GWAS, unless combined with larger genotype-ready cohorts.

## MR

The MR workflow follows the Nature 2026 virome analysis pattern: weighted median, IVW, MR-Egger, and ConMix. The primary instrument set should exclude the extended MHC region. This workspace currently contains scaffold implementations for weighted median, IVW, MR-Egger, allele harmonization, leave-one-out sensitivity, and relaxed-threshold instrument exploration. Final ConMix should be run in a validated MR package.

## Reporting Constraint

Do not claim EBV causation from association alone. Label evidence as causal, immune-state-marker, or inconclusive only after subtype-specific association, instrument quality, MR sensitivity, and metadata completeness are reviewed.
