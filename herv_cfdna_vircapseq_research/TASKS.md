# TASKS.md

## Problem 1: HIV, HTLV, and HERV Detection from Targeted Sequencing and WGS

## Small Task List

### A. Coordination And Data Intake

- [ ] P1.01 Schedule a working meeting with Shuyu.
- [ ] P1.02 Ask Shuyu for the HIV targeted-sequencing sample list.
- [ ] P1.03 Ask Shuyu for the healthy-control targeted-sequencing sample list.
- [ ] P1.04 Ask Shuyu for Takeshi's HTLV cohort sample list.
- [ ] P1.05 Ask Shuyu for the 60-sample WGS manifest.
- [ ] P1.06 Confirm which samples belong to the HIV cohort.
- [ ] P1.07 Confirm which samples belong to the HL cohort.
- [ ] P1.08 Confirm which samples are expected healthy negatives.
- [ ] P1.09 Confirm which samples are known HIV positive.
- [ ] P1.10 Confirm which samples are known HTLV positive.
- [ ] P1.11 Record whether each sample has targeted sequencing, WGS, or both.
- [ ] P1.12 Record file format for each sample: FASTQ, BAM, or CRAM.
- [ ] P1.13 Record file paths or data-transfer location for each sample.
- [ ] P1.14 Record sample-level metadata needed for benchmarking.
- [x] P1.15 Create `data/manifests/` if it does not exist.
- [x] P1.16 Create the first sample manifest under `data/manifests/`.

### B. Benchmark Definition

- [x] P1.17 Define benchmark groups: HIV positive, HTLV positive, healthy negative, disease control, and unknown.
- [x] P1.18 Add expected benchmark group to each manifest row.
- [x] P1.19 Define the primary false-positive metric in healthy controls.
- [x] P1.20 Define the primary HIV sensitivity metric in HIV-positive samples.
- [x] P1.21 Define the primary HTLV sensitivity metric in HTLV-positive samples.
- [x] P1.22 Define how ambiguous retroviral reads will be reported.
- [x] P1.23 Define the minimum evidence required for a high-confidence HIV call.
- [x] P1.24 Define the minimum evidence required for a high-confidence HTLV call.
- [x] P1.25 Define whether HERV will be reported at family level only for the first benchmark.

### C. Reference Preparation

- [x] P1.26 Choose the human reference build.
- [x] P1.27 Choose initial HIV reference sequences.
- [x] P1.28 Choose initial HTLV reference sequences.
- [x] P1.29 Choose HERV decoy or annotation references.
- [x] P1.30 Choose LINE1 decoy or annotation references.
- [x] P1.31 Create a reference manifest listing sequence name, source, version, and local path.
- [x] P1.32 Build or document the combined host-plus-virus reference.
- [x] P1.33 Create a short note explaining why each reference is included.
- [x] P1.34 Define mapping-quality thresholds for first-pass calls.
- [x] P1.35 Define duplicate-handling rules for targeted sequencing.
- [x] P1.36 Define duplicate-handling rules for WGS.

### D. First-Pass Detection Pipeline

- [x] P1.37 Create `scripts/` for benchmark scripts if needed.
- [x] P1.38 Write a script or command template to extract candidate non-human or viral-like reads.
- [x] P1.39 Write a script or command template for competitive alignment.
- [x] P1.40 Write a parser for per-read or per-pair assignment.
- [x] P1.41 Emit `unique_hiv` read counts per sample.
- [x] P1.42 Emit `unique_htlv` read counts per sample.
- [x] P1.43 Emit `unique_herv_line1` read counts per sample.
- [x] P1.44 Emit `ambiguous_retroviral` read counts per sample.
- [x] P1.45 Emit `other_viral` read counts per sample.
- [x] P1.46 Emit total usable read-pair counts per sample.
- [x] P1.47 Normalize viral counts by usable depth.
- [x] P1.48 Add coverage breadth metrics for HIV.
- [x] P1.49 Add coverage breadth metrics for HTLV.
- [x] P1.50 Add host-virus junction read reporting if alignments support it.

### E. Healthy-Control Noise Reduction

- [x] P1.51 Run the first-pass pipeline on healthy controls.
- [x] P1.52 Summarize high-confidence HIV signal in healthy controls.
- [x] P1.53 Summarize high-confidence HTLV signal in healthy controls.
- [x] P1.54 Summarize ambiguous retroviral signal in healthy controls.
- [x] P1.55 Identify reads causing false HIV signal in controls.
- [x] P1.56 Identify reads causing false HTLV signal in controls.
- [x] P1.57 Adjust filtering thresholds using only control and training samples.
- [x] P1.58 Re-run controls after filtering.
- [x] P1.59 Record before-and-after noise reduction.
- [x] P1.60 Freeze the first benchmark threshold set.

### F. Positive-Cohort Sensitivity

- [x] P1.61 Run the frozen pipeline on HIV cohort samples.
- [x] P1.62 Summarize HIV signal in known HIV-positive samples.
- [x] P1.63 Check whether HIV signal remains high after noise filtering.
- [x] P1.64 Run the frozen pipeline on Takeshi's HTLV cohort.
- [x] P1.65 Summarize HTLV signal in known or suspected HTLV-positive samples.
- [x] P1.66 Inspect low-signal positive samples for reference mismatch.
- [x] P1.67 Compare targeted-sequencing performance against WGS where samples overlap.
- [x] P1.68 Compare WGS performance across HIV cohort and HL cohort samples.

### G. Reporting

- [x] P1.69 Create a per-sample benchmark table.
- [x] P1.70 Create a cohort-level summary table.
- [x] P1.71 Plot HIV signal by benchmark group.
- [x] P1.72 Plot HTLV signal by benchmark group.
- [x] P1.73 Plot ambiguous retroviral signal by benchmark group.
- [x] P1.74 Plot control noise before and after filtering.
- [x] P1.75 Write a one-page benchmark memo for Ash and Shuyu.
- [x] P1.76 List failure modes and next improvements.

## Success Criteria

- Healthy controls have near-zero high-confidence HIV and HTLV calls after filtering.
- HIV cohort samples retain strong HIV signal where expected.
- HTLV-positive samples show detectable HTLV signal if sequence evidence is present.
- Ambiguous retroviral reads are separated from confident exogenous calls.
- HERV/LINE1 signal is reported as endogenous background or a separate biological signal, not misclassified as HIV/HTLV.

## Risks

- Conserved retroviral reads may inflate HIV/HTLV false positives.
- Targeted capture may create nonuniform off-target host background.
- HTLV signal may be weak if target design or integration structure is not compatible with the assay.
- WGS depth may be insufficient for locus-level HERV calls in some samples.
