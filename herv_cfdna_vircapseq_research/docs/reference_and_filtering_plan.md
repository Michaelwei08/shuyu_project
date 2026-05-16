# Reference And Filtering Plan

## Initial Reference Choices

- Human reference: `GRCh38/hg38`, matching incoming aligned BAM/CRAM files whenever possible.
- HIV reference: HIV-1 HXB2 RefSeq `NC_001802.1` plus subtype-aware HIV references if Shuyu's cohort subtype information is available.
- HTLV reference: HTLV-1 RefSeq `NC_001436.1` and HTLV-2 RefSeq `NC_001488.1`.
- HERV decoys: HERV consensus or curated HERV sequences from Dfam/RepeatMasker resources.
- LINE1 decoys: LINE1/L1HS consensus or Dfam/RepeatMasker LINE1 resources.

## Combined Reference Policy

The first benchmark should align or classify against a combined host-plus-virus reference:

- human genome
- HIV references
- HTLV references
- HERV decoys
- LINE1 decoys
- optional other viral references if the assay uses a broader viral panel

The purpose is to prevent conserved retroviral reads from being forced into HIV or HTLV calls.

## First-Pass Thresholds

- Minimum mapping quality: `20`
- Minimum alignment length: `40 bp`
- Ambiguous assignment rule: if top retroviral hits are within `5 alignment-score units`, report `ambiguous_retroviral`.
- High-confidence HIV/HTLV calls require unique viral assignment after HERV and LINE1 decoy competition.
- Highest-confidence HIV/HTLV calls require either coverage over discriminative viral regions or host-virus junction support.

## Duplicate Handling

- Targeted sequencing: mark or collapse duplicate read pairs before final sample-level counting. Keep pre-deduplication counts for QC.
- WGS: mark duplicates and report both raw and duplicate-filtered viral counts when possible.

## HERV Reporting Scope

The first benchmark should report HERV/LINE1 at family level. Locus-level HERV calls should be postponed until the WGS depth and read structure are reviewed.
