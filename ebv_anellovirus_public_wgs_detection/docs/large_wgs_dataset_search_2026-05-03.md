# Large WGS Dataset Search: Second Pass

## What Changed

This pass prioritized large, true WGS-scale datasets. I used the practical heuristic that normal human WGS BAM/FASTQ/CRAM data are often tens of GB per sample and treated files greater than or equal to 40 GB as strong evidence of true WGS-scale data rather than WES, panel, viral-only, or low-depth metagenomics.

## Strongest NIH/GDC Leads

1. `CGCI-HTMCP-DLBCL` is now the top NIH/GDC candidate. It is HIV-associated DLBCL, directly relevant to both the HIV and lymphoma-virus questions. The local GDC API manifest shows 116 WGS BAM files greater than or equal to 40 GB across 66 cases, with 38.7 TB of WGS-related files.

2. `TCGA-DLBC` is a strong DLBCL comparator. The local GDC API manifest shows 98 WGS BAM files greater than or equal to 40 GB across 42 cases, with 26.7 TB of WGS-related files.

3. `CGCI-BLGSP` is not DLBCL/HL/PTLD, but it is EBV-relevant lymphoma. The local GDC API manifest shows 463 WGS BAM files greater than or equal to 40 GB across 252 cases, with 154.4 TB of WGS-related files.

The downloaded GDC metadata summaries are in:

- `data/external_manifests/gdc_lymphoma_wgs_files_expanded_2026-05-03.csv`
- `data/external_manifests/gdc_lymphoma_wgs_project_summary_2026-05-03.csv`

## Strongest Disease-Specific Controlled Leads

- HL: EGA `EGAD00001009818` / `EGAS00001006884`, 53 BAM files and 9.1 TB.
- EBV-positive DLBCL: EGA `EGAD00001006858` / `EGAS00001004941`, 8 matched tumor-normal WGS patients plus targeted resequencing.
- DLBCL/FL: ICGC MMML-Seq `EGAS00001002199`, a large controlled mature B-cell lymphoma WGS resource.

## Recent Nature/Nature Genetics Datasets

These are not lymphoma WGS cohorts, but they are directly relevant as method and causal-inference anchors:

- Kamitaki et al., Nature 2026: UK Biobank, All of Us, and SPARK high-coverage WGS for EBV, HHV-6B, HHV-7, and anellovirus detection.
- Nyeo et al., Nature 2026: UK Biobank and All of Us EBV DNAemia GWAS with public summary statistics.
- Schmidt et al., Nature 2026: UK Biobank and All of Us EBVread-positive GWAS; GWAS Catalog IDs include `GCST90809298-GCST90809306`.
- Sasa et al., Nature Genetics 2025: 6,321 Japanese WGS samples for eHHV-6 and anellovirus associations, with JGA/HGVD access paths.

## Practical Recommendation

Submit controlled-access requests first for `CGCI-HTMCP-DLBCL`, `TCGA-DLBC`, `CGCI-BLGSP`, and EGA `EGAD00001009818`. `CGCI-HTMCP-DLBCL` should be prioritized because it is the cleanest intersection of HIV, DLBCL, WGS, and NIH/GDC access.

For immediate open-data work, use `PRJNA772081` only as a smaller smoke test. Its WGS runs are open and real, but the compressed run sizes are mostly below the 40 GB heuristic and the cohort is not a clean HL/DLBCL/PTLD public-WGS validation set.

PTLD remains a public-data gap. I found PTLD EBV genome sequencing and PTLD metagenomic publications, but no large public human WGS cohort suitable for EBV/anellovirus detection.

## Files Added

- `data/cohorts/large_wgs_dataset_candidates_2026-05-03.csv`
- `data/external_manifests/gdc_lymphoma_wgs_files_expanded_2026-05-03.csv`
- `data/external_manifests/gdc_lymphoma_wgs_project_summary_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJNA772081_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJNA707099_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJNA551423_2026-05-03.csv`
- `data/external_manifests/ncbi_runinfo_PRJDB18873_2026-05-03.csv`
