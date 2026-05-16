# TASKS.md

## Problem 2: EBV, Other Herpesviruses, and Anelloviruses from Public WGS

## Small Task List

### A. Public Cohort Discovery

- [x] P2.01 Create `data/cohorts/`.
- [x] P2.02 Create a public cohort inventory table.
- [x] P2.03 Search for public HL WGS cohorts.
- [x] P2.04 Search for public DLBCL WGS cohorts.
- [x] P2.05 Search for public PTLD WGS cohorts.
- [x] P2.06 Record accession IDs for each candidate cohort.
- [x] P2.07 Record publication or data-source links for each candidate cohort.
- [x] P2.08 Record whether raw reads are available.
- [x] P2.09 Record whether BAM or CRAM files are available.
- [x] P2.10 Record whether access is open or controlled.
- [x] P2.11 Record expected sample count for each cohort.
- [x] P2.12 Record sample type for each cohort: tumor, normal, blood, saliva, or cell line.
- [x] P2.13 Record whether clinical metadata are available.
- [x] P2.14 Record whether EBV status is already annotated.
- [x] P2.15 Mark cohorts that are immediately runnable.

### B. Prioritization

- [x] P2.16 Rank cohorts by data accessibility.
- [x] P2.17 Rank cohorts by disease relevance.
- [x] P2.18 Rank cohorts by sample count.
- [x] P2.19 Rank cohorts by metadata completeness.
- [x] P2.20 Select one pilot HL cohort if available.
- [x] P2.21 Select one pilot DLBCL cohort if available.
- [x] P2.22 Select one pilot PTLD cohort if available.
- [x] P2.23 Document why any disease group lacks an immediate pilot cohort.

### C. Viral Reference Panel

- [x] P2.24 Create `data/references/`.
- [x] P2.25 List EBV reference genome candidates.
- [x] P2.26 List HHV-6A reference genome candidates.
- [x] P2.27 List HHV-6B reference genome candidates.
- [x] P2.28 List HHV-7 reference genome candidates.
- [x] P2.29 List CMV reference genome candidates.
- [x] P2.30 List KSHV reference genome candidates.
- [x] P2.31 List HSV reference genome candidates.
- [x] P2.32 List representative anellovirus references.
- [x] P2.33 Add reference source and version metadata.
- [x] P2.34 Build a viral reference manifest.
- [x] P2.35 Decide whether to include the 31-virus panel from the Nature 2026 paper.
- [x] P2.36 Decide whether to include a broader RefSeq panel as a sensitivity run.

### D. Detection Workflow Design

- [x] P2.37 Define input mode for FASTQ data.
- [x] P2.38 Define input mode for BAM or CRAM data.
- [x] P2.39 Define how to extract unmapped or poorly mapped reads.
- [x] P2.40 Define whether all reads should also be screened.
- [x] P2.41 Define the viral alignment command template.
- [x] P2.42 Define minimum mapping quality for viral reads.
- [x] P2.43 Define minimum alignment length for viral reads.
- [x] P2.44 Define duplicate-handling rules.
- [x] P2.45 Define sample-level normalization denominator.
- [x] P2.46 Define viral breadth-of-coverage metric.
- [x] P2.47 Define EBV-positive call threshold for the pilot.
- [x] P2.48 Define herpesvirus-positive call threshold for the pilot.
- [x] P2.49 Define anellovirus-positive call threshold for the pilot.

### E. Pilot Execution

- [ ] P2.50 Download or stage a small pilot subset.
- [ ] P2.51 Run the workflow on one known or suspected EBV-positive sample if available.
- [ ] P2.52 Run the workflow on one expected negative or low-virus sample if available.
- [ ] P2.53 Run the workflow on one anellovirus-suspected sample if available.
- [x] P2.54 Generate per-sample viral read-count tables.
- [x] P2.55 Generate per-sample viral breadth tables.
- [x] P2.56 Generate per-sample normalized viral-load tables.
- [x] P2.57 Inspect EBV read distribution across the EBV genome.
- [x] P2.58 Inspect anellovirus read distribution across reference genomes.
- [x] P2.59 Flag samples with likely mapping artifacts.
- [x] P2.60 Document pilot runtime and storage requirements.

### F. Subtype And Association Preparation

- [x] P2.61 Define EBV feature columns.
- [x] P2.62 Define other-herpesvirus feature columns.
- [x] P2.63 Define anellovirus feature columns.
- [x] P2.64 Define pan-viral burden features.
- [x] P2.65 Map WGS features to existing viral subtype labels.
- [x] P2.66 Create a draft subtype assignment table.
- [x] P2.67 Identify clinical metadata needed for association testing.
- [x] P2.68 Join viral calls to available clinical metadata.
- [x] P2.69 Check missingness in clinical metadata.
- [x] P2.70 Decide which associations are feasible in each cohort.

### G. Reporting

- [x] P2.71 Create a cohort inclusion table.
- [x] P2.72 Create a viral detection summary table by disease.
- [x] P2.73 Plot EBV signal by disease and cohort.
- [x] P2.74 Plot anellovirus signal by disease and cohort.
- [x] P2.75 Plot EBV versus anellovirus signal.
- [x] P2.76 Compare WGS viral calls against known EBV annotations if available.
- [x] P2.77 Write a short public-WGS feasibility memo.
- [x] P2.78 List cohorts ready for inclusion in the viral subtype paper.

## Success Criteria

- EBV/herpesvirus and anellovirus detection is reproducible from public WGS.
- Viral calls include enough quality metrics to distinguish true signal from mapping artifacts.
- At least one public HL, DLBCL, or PTLD WGS cohort can be added to the viral subtype analysis.
- Clinical association results can be repeated or clearly marked as underpowered or metadata-limited.

## Risks

- Some public WGS cohorts may be controlled access or lack raw reads.
- Tumor-only WGS may complicate interpretation of systemic anellovirus signal.
- EBV reads may be sparse in blood-derived WGS and stronger in tumor WGS, depending on sample source.
- Batch effects and library-prep differences may dominate low-count viral signals.
