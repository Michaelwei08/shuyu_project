# GWAS Execution Plan

## First Phenotype

Use `ebv_load` as the first GWAS phenotype only if sample size, genotype data, ancestry PCs, and core covariates are available.

## Input Files

- Viral phenotype table.
- Covariate table.
- Genotype manifest.
- Sample inclusion list.

## Basic QC

- Confirm sample IDs match across phenotype, covariate, and genotype inputs.
- Confirm ancestry PCs are available.
- Confirm phenotype missingness is acceptable.
- Confirm batch and sequencing depth are available for adjustment.
- Flag MHC-region associations separately.

Run the scaffolded input check before any GWAS:

```powershell
python scripts/check_gwas_input_qc.py data/templates/viral_phenotype_template.csv data/templates/covariate_template.csv data/templates/genotype_manifest_template.csv results/example_gwas_input_qc.csv
```

The example template is expected to fail because real viral phenotypes, covariates, and genotype QC have not been loaded.

## Command Template

The exact GWAS command depends on genotype format and available software. For BGEN/PLINK-style workflows, prepare:

```powershell
plink2 --bgen input.bgen ref-first --sample input.sample --pheno gwas_phenotype.tsv --covar gwas_covariates.tsv --glm hide-covar --out ebv_load_gwas
```

Do not run final GWAS until phenotype version and sample inclusion criteria are frozen.
