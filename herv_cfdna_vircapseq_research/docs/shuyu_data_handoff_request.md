# Shuyu Data Handoff Request

## Purpose

To run the P1 HIV/HTLV/HERV benchmark on real data, we need the sample manifest and file locations for targeted sequencing and WGS.

## Requested Tables

Please provide one manifest table with these fields:

- `sample_id`
- `subject_id`
- `cohort`: `hiv_cohort`, `healthy_control`, `takeshi_htlv`, `local_60_wgs`, or other
- `disease_group`: HIV, HL, healthy, PTLD, DLBCL, or other
- `assay_type`: `targeted` or `WGS`
- `expected_group`: `hiv_positive`, `htlv_positive`, `healthy_negative`, `disease_control`, or `unknown`
- `hiv_status`: positive, negative, or unknown
- `htlv_status`: positive, negative, or unknown
- `file_format`: FASTQ, BAM, or CRAM
- `read1_path`, `read2_path`, or `bam_or_cram_path`
- `reference_build`
- `estimated_depth`
- `library_prep`
- `batch`
- `notes`

## Cohorts Needed

- HIV targeted-sequencing cohort.
- Healthy-control targeted-sequencing cohort.
- Takeshi's HTLV cohort.
- 60 local WGS samples from HIV and HL cohorts.

## Benchmark Questions

- Which healthy controls should be used for threshold calibration?
- Which HIV samples are known positives and should be held out for sensitivity testing?
- Which HTLV samples are known positives, and are they HTLV-1, HTLV-2, or unknown type?
- Which samples have both targeted sequencing and WGS for overlap comparison?
- Are file paths local, shared-drive, cloud bucket, or controlled-access handoff paths?

## Ready-To-Run Scaffold

Once the manifest and competitive alignment outputs are available, run the command sequence in `docs/p1_full_pipeline_runbook.md`.
