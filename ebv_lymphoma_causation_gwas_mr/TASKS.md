# TASKS.md

## Problem 3: EBV Causation in HL, PTLD, and DLBCL

## Small Task List

### A. Phenotype Definition

- [x] P3.01 Create `docs/` for phenotype definitions.
- [x] P3.02 Define the EBV viral-load phenotype.
- [x] P3.03 Define the anellovirus viral-load phenotype.
- [x] P3.04 Define the other-herpesvirus viral-load phenotype.
- [x] P3.05 Define the pan-viral burden phenotype.
- [x] P3.06 Define the V1 pan-viral subtype.
- [x] P3.07 Define the V3 high-EBV low-anellovirus subtype.
- [x] P3.08 Define how low-depth samples will be handled.
- [x] P3.09 Define how missing viral calls will be handled.
- [x] P3.10 Freeze phenotype definitions before association testing.

### B. Data Inventory

- [x] P3.11 List the 60 local WGS samples.
- [x] P3.12 Record which local samples are HIV cohort samples.
- [x] P3.13 Record which local samples are HL cohort samples.
- [x] P3.14 Link local WGS samples to viral detection outputs.
- [x] P3.15 Link public WGS samples to viral detection outputs.
- [x] P3.16 Record disease labels: HL, PTLD, DLBCL, control, or other.
- [x] P3.17 Record sample source: tumor, blood, saliva, or normal.
- [x] P3.18 Record sequencing depth for each sample.
- [x] P3.19 Record batch or sequencing center if available.
- [x] P3.20 Record tumor purity if available.
- [x] P3.21 Record transplant status if available.
- [x] P3.22 Record immunosuppression status if available.
- [x] P3.23 Record treatment variables if available.
- [x] P3.24 Record outcome variables if available.

### C. GWAS Feasibility

- [x] P3.25 Determine whether each dataset has genotype data.
- [x] P3.26 Determine whether each dataset has imputed genotypes.
- [x] P3.27 Determine whether each dataset has ancestry principal components.
- [x] P3.28 Determine whether sample size is adequate for GWAS.
- [x] P3.29 Mark datasets as GWAS-ready, association-only, or phenotype-only.
- [x] P3.30 Document why the 60 local WGS samples are or are not GWAS-ready.
- [x] P3.31 Identify any external GWAS summary statistics for HL.
- [x] P3.32 Identify any external GWAS summary statistics for DLBCL.
- [x] P3.33 Identify any external GWAS summary statistics for PTLD if available.
- [x] P3.34 Identify immune-mediated comparator traits if needed.

### D. Covariate Construction

- [x] P3.35 Create a covariate table template.
- [x] P3.36 Add age to the covariate table.
- [x] P3.37 Add sex to the covariate table.
- [x] P3.38 Add ancestry PCs to the covariate table.
- [x] P3.39 Add cohort to the covariate table.
- [x] P3.40 Add sequencing batch to the covariate table.
- [x] P3.41 Add sample type to the covariate table.
- [x] P3.42 Add sequencing depth to the covariate table.
- [x] P3.43 Add tumor purity if available.
- [x] P3.44 Add transplant status if available.
- [x] P3.45 Add immunosuppression variables if available.
- [x] P3.46 Check covariate missingness.
- [x] P3.47 Decide missing-data handling for each covariate.

### E. Clinical Association Analyses

- [x] P3.48 Test EBV load association with HL status.
- [x] P3.49 Test EBV load association with PTLD status.
- [x] P3.50 Test EBV load association with DLBCL status.
- [x] P3.51 Test anellovirus load association with HL status.
- [x] P3.52 Test anellovirus load association with PTLD status.
- [x] P3.53 Test anellovirus load association with DLBCL status.
- [x] P3.54 Test V1 subtype association with clinical outcomes.
- [x] P3.55 Test V3 subtype association with clinical outcomes.
- [x] P3.56 Repeat association tests after depth adjustment.
- [x] P3.57 Repeat association tests after batch adjustment.
- [x] P3.58 Repeat association tests within each disease group.
- [x] P3.59 Summarize association effect sizes and confidence intervals.

### F. GWAS Execution

- [x] P3.60 Select phenotype for first GWAS.
- [x] P3.61 Prepare GWAS phenotype file.
- [x] P3.62 Prepare GWAS covariate file.
- [x] P3.63 Prepare genotype input paths.
- [x] P3.64 Run basic genotype/sample QC.
- [x] P3.65 Run GWAS for EBV load if feasible.
- [x] P3.66 Run GWAS for anellovirus load if feasible.
- [x] P3.67 Run GWAS for V1 subtype if feasible.
- [x] P3.68 Run GWAS for V3 subtype if feasible.
- [x] P3.69 Generate Manhattan and QQ plots for each GWAS.
- [x] P3.70 Summarize genome-wide significant loci.
- [x] P3.71 Flag MHC-region associations separately.

### G. Mendelian Randomization

- [x] P3.72 Select candidate EBV-load instruments from GWAS or published summary statistics.
- [x] P3.73 Remove MHC-region instruments for the primary MR analysis.
- [x] P3.74 Clump instruments for linkage disequilibrium.
- [x] P3.75 Record instrument effect sizes on EBV load.
- [x] P3.76 Obtain outcome effect sizes for HL.
- [x] P3.77 Obtain outcome effect sizes for DLBCL.
- [x] P3.78 Obtain outcome effect sizes for PTLD if available.
- [x] P3.79 Harmonize exposure and outcome alleles.
- [x] P3.80 Run weighted median MR.
- [x] P3.81 Run inverse-variance weighted MR.
- [x] P3.82 Run MR-Egger.
- [x] P3.83 Run contamination mixture MR.
- [x] P3.84 Record MR effect sizes, standard errors, and P values.
- [x] P3.85 Test MR-Egger intercept for directional pleiotropy.
- [x] P3.86 Run leave-one-instrument-out sensitivity analysis if enough instruments exist.
- [x] P3.87 Repeat MR with alternative instrument thresholds if feasible.
- [x] P3.88 Document any MR method that is infeasible and why.

### H. Subtype Interpretation

- [x] P3.89 Compare EBV effects in V1 versus V3.
- [x] P3.90 Compare anellovirus effects in V1 versus V3.
- [x] P3.91 Test whether V3 is enriched for disease-associated EBV signal.
- [x] P3.92 Test whether V1 is enriched for immunosuppression markers.
- [x] P3.93 Analyze HL separately.
- [x] P3.94 Analyze PTLD separately.
- [x] P3.95 Analyze DLBCL separately.
- [x] P3.96 Write a causal-interpretation table: supports causation, supports immune-state marker, inconclusive.
- [x] P3.97 Write a short methods note for GWAS and MR.
- [x] P3.98 Write a short results memo for Ash and Shuyu.

## Success Criteria

- Viral subtype definitions are frozen before clinical association testing.
- GWAS and MR inputs are reproducible and provenance-tracked.
- Four MR methods from the Nature 2026 paper are run or explicitly marked infeasible with reasons.
- HL, PTLD, and DLBCL are analyzed separately before any pooled interpretation.
- V1 and V3 are evaluated as potentially different biological states.

## Risks

- The 60 local WGS samples are likely underpowered for GWAS and may be better suited for phenotype validation.
- PTLD causation may be confounded by transplant-related immunosuppression.
- EBV tumor presence, EBV blood DNA load, and EBV immune response are related but not interchangeable.
- MHC instruments can create pleiotropy risk in MR.
- Public WGS metadata may be insufficient for transplant status, immune suppression, or treatment adjustment.
