# Public WGS Viral Detection Workflow

## Input Modes

- FASTQ: align reads to human and viral references or run a direct host-plus-virus competitive workflow.
- BAM/CRAM: extract unmapped, poorly mapped, and optionally all reads for viral screening.
- Tumor-only WGS: report viral calls but interpret anellovirus signal cautiously.
- Tumor-normal WGS: compare tumor and normal viral signal where both are available.

## Initial Extraction Policy

- For BAM/CRAM, extract unmapped reads, supplementary alignments, soft-clipped reads if practical, and low mapping-quality reads for first-pass viral screening.
- For FASTQ, run direct viral screening and host-plus-virus competitive alignment when compute allows.
- For pilot work, screen all reads only on small subsets or after downsampling, then decide whether unmapped-only extraction loses sensitivity.

## Initial Thresholds

- Minimum viral mapping quality: `20`
- Minimum viral alignment length: `40 bp`
- Duplicate handling: remove duplicate read-pairs for final viral load, but retain duplicate-aware QC counts.
- Normalization denominator: usable nonduplicate read pairs, or total screened read pairs if usable read pairs are unavailable.
- Breadth metric: fraction of observed viral reference positions with depth at least `1`.
- EBV-positive pilot threshold: nonzero EBV reads plus breadth or multiple independent read-pairs; exact threshold must be calibrated per dataset.
- Herpesvirus-positive pilot threshold: nonzero family-specific reads plus mapping-quality review.
- Anellovirus-positive pilot threshold: nonzero anellovirus reads after duplicate removal, interpreted as a continuous burden unless strong coverage supports a binary call.

## First-Pass Outputs

- Total usable read pairs.
- Viral read pairs by virus.
- Viral read pairs by family.
- Breadth of coverage by virus.
- Mean depth over viral genome by virus.
- Mapping-quality distribution.
- Duplicate-aware normalized viral load.

## Initial Call Categories

- `ebv_detected`
- `other_herpesvirus_detected`
- `anellovirus_detected`
- `pan_viral_high`
- `no_confident_viral_signal`
- `artifact_review_required`

## Pilot Strategy

Start with the smallest accessible subset from each disease group. Run EBV and anellovirus detection separately, then add subtype assignment only after basic viral calls pass mapping-quality and coverage checks.
