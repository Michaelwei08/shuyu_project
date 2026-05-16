# Viral Feature And Subtype Definitions

## Feature Columns

- `ebv_reads`: duplicate-filtered EBV read-pair count.
- `ebv_per_million`: EBV reads per million usable read pairs.
- `other_herpesvirus_reads`: duplicate-filtered non-EBV herpesvirus read-pair count.
- `other_herpesvirus_per_million`: non-EBV herpesvirus reads per million usable read pairs.
- `anellovirus_reads`: duplicate-filtered anellovirus read-pair count.
- `anellovirus_per_million`: anellovirus reads per million usable read pairs.
- `pan_viral_reads`: EBV plus other herpesvirus plus anellovirus plus other viral reads.
- `pan_viral_per_million`: pan-viral reads per million usable read pairs.

## Draft Subtype Logic

- `V1_pan_viral`: broad viral signal across at least two viral families, especially when EBV and anellovirus are both elevated.
- `V3_high_ebv_low_anello`: EBV signal elevated while anellovirus signal remains low or absent.
- `no_confident_viral_signal`: no viral family passes the pilot threshold.
- `artifact_review_required`: viral signal dominated by low-complexity, duplicate, or low-breadth evidence.

## Notes

Subtype thresholds must be calibrated after the first real public WGS pilot. The definitions here specify columns and logic, not final biological cut points.
