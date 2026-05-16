# GWAS Feasibility Notes

## Local 60 WGS Samples

The local 60-sample WGS set is likely useful for viral phenotype validation and clinical association sanity checks. It is unlikely to be adequately powered for discovery GWAS by itself unless it is combined with a larger cohort or used only for targeted genotype checks.

## Public WGS Samples

Public WGS datasets can support viral phenotype detection if raw reads or BAM/CRAM files are accessible. They may not automatically support GWAS unless genotype calls, imputation, ancestry PCs, and consent-compatible phenotype metadata are available.

## External Outcome GWAS Candidates

- Hodgkin lymphoma: published GWAS meta-analysis source identified.
- DLBCL: published European and East Asian GWAS sources identified.
- PTLD: no clearly MR-ready GWAS summary statistics identified yet.
- Immune-mediated comparator traits: candidate OpenGWAS endpoints identified for multiple sclerosis, systemic lupus erythematosus, rheumatoid arthritis, and inflammatory bowel disease. These are comparator outcomes only; use them after allele harmonization and ancestry checks.

## Analysis Status Logic

- `gwas_ready`: genotype data, ancestry PCs, phenotype, and required covariates are present.
- `association_only`: viral phenotype and clinical metadata are present, but genotype inputs are incomplete.
- `phenotype_only`: viral phenotype can be summarized, but association covariates are incomplete.
- `blocked`: required data are unavailable or access is not approved.
