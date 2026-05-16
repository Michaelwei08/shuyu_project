# Missing Data Rules

## Viral Phenotypes

- Missing sequencing data is encoded as missing.
- Failed viral-detection QC is encoded as `filtered`.
- Adequate sequencing with no confident viral reads can be encoded as zero.

## Covariates

- Age: keep missing as `NA`; do not impute for primary analysis unless required by model software.
- Sex: keep missing as `NA`; exclude from sex-adjusted models if unavailable.
- Ancestry PCs: required for GWAS; samples without PCs are not GWAS-ready.
- Batch: required when multiple sequencing batches are combined.
- Sample type: required when tumor, blood, saliva, and normal samples are mixed.
- Sequencing depth: required for viral-load association models.
- Tumor purity: include when available for tumor WGS; do not require for blood WGS.
- Transplant and immunosuppression variables: required for PTLD interpretation when available.

## Analysis Status Labels

- `gwas_ready`: genotype data, ancestry PCs, phenotype, and required covariates are present.
- `association_only`: viral phenotype and clinical metadata are present, but GWAS inputs are incomplete.
- `phenotype_only`: viral phenotype can be summarized, but association covariates are incomplete.
- `blocked`: required data are unavailable or access is not approved.
