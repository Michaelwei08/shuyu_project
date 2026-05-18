# P1 Interim Results: Targeted HTLV and WGS HIV/HL Retrovirus Screen

Date: 2026-05-18

## Scope

This note summarizes the current P1 benchmark for distinguishing exogenous retroviruses from HERV/LINE1 background using Shuyu's shared sequencing data.

Assay tracks completed so far:

- Targeted HTLV sequencing from TCL samples.
- WGS from HIV-labeled DLBCL/BL samples and HL HIV-negative controls.

The current implementation uses a competitive HIV/HTLV/HERV/LINE1 reference, BWA-MEM alignment, mapped-read BAM output, and filtered read counting. The current filter requires MAPQ >= 20 and aligned length >= 60 bp, with read IDs deduplicated within each reference category.

## Targeted HTLV Result

The targeted HTLV benchmark is a strong positive-control result.

Summary:

- Samples analyzed: 71.
- HTLV1 positive threshold: >= 100 filtered read IDs.
- HTLV1 positive: 67/71.
- Low positive: 3/71.
- Zero/possible QC failure: 1/71.
- Median HTLV1 filtered read IDs: 8,361.
- Maximum HTLV1 filtered read IDs: 218,993.
- HIV1 nonzero samples: 0/71.
- HIV2 nonzero samples: 0/71.
- HTLV2 nonzero samples: 0/71.

Manifest/QC audit:

- The targeted QC table contained 75 paired-complete targeted samples.
- The final complete manifest used for the full run contained 71 samples.
- The 4 omitted samples were excluded because of R1/R2 size imbalance.

Interpretation:

The competitive-filtered workflow robustly detects HTLV1 in targeted HTLV-enriched samples while suppressing spurious HIV1, HIV2, and HTLV2 signal. This supports the method as a working positive-control benchmark for exogenous retrovirus detection in targeted sequencing.

## WGS HIV/HL Result

The WGS benchmark is a strong specificity result and a weak-but-plausible sensitivity result.

Full WGS summary:

- Samples analyzed: 60.
- HIV-labeled WGS samples: 37.
- HL negative-control WGS samples: 23.
- HIV signal in HIV-labeled WGS: 2/37.
- HIV signal in HL controls: 0/23.
- HTLV signal in HL controls: 0/23.
- Median HERV filtered read IDs:
  - HIV-labeled WGS: 52,957.
  - HL controls: 52,447.
- Median LINE1 filtered read IDs:
  - HIV-labeled WGS: 4,212,518.
  - HL controls: 4,073,247.

The two WGS HIV-positive calls were each single-fragment HIV1 detections. Compact alignment inspection showed:

- HIV1 reference: K03455.1.
- MAPQ: 60.
- Aligned length: >120 bp.
- Secondary alignment score: XS:i:0.
- Count after deduplication: 1 HIV1 fragment per positive sample.

Interpretation:

The full WGS screen produced no HIV or HTLV signal in HL controls, supporting high specificity under the current competitive-filtered approach. Rare HIV1 signal was found only in HIV-labeled WGS samples, but each positive sample had only one deduplicated HIV1 fragment. This should be described as rare low-level HIV1 detectability, not as robust HIV quantification or viral-load measurement.

The similar HERV/LINE1 background between HIV-labeled WGS and HL controls argues against the HIV calls being explained by a broad retroelement-background shift.

## Current Conclusion

The P1 workflow has met the first practical benchmark:

- It detects abundant HTLV1 signal in targeted HTLV-positive/enriched data.
- It keeps healthy/negative-control-like WGS samples essentially free of HIV/HTLV false positives.
- It can detect rare HIV1-like WGS fragments in a small subset of HIV-labeled WGS samples, but WGS sensitivity is low.

The current evidence supports using this pipeline as a specificity-focused exogenous retrovirus screen, with targeted sequencing as the stronger assay for positive detection.

## Recommended Next Steps

1. Freeze the current counting rule as the first-pass benchmark:
   - Competitive HIV/HTLV/HERV/LINE1 reference.
   - MAPQ >= 20.
   - Aligned length >= 60 bp.
   - Deduplicate by read ID and reference category.

2. Add a compact fragment-level audit table for the two WGS HIV1-positive calls:
   - sample group, reference, position, MAPQ, CIGAR, aligned length, NM, AS, XS.
   - no read sequence or quality string.

3. Ask Shuyu for the targeted HIV cohort data when ready:
   - This is the correct sensitivity benchmark for HIV detection.
   - WGS appears too sparse for reliable HIV detection except rare single-fragment calls.

4. Run the same workflow on the HIV targeted cohort plus healthy controls:
   - Primary success criterion: HIV cohort retains high HIV signal.
   - Primary specificity criterion: healthy controls remain near zero.

5. Decide whether HERV enumeration is a separate second-pass analysis:
   - Current HERV counts are background/decoy counts, not locus-specific HERV enumeration.
   - Locus-level HERV insertion/activity analysis will require a different workflow.

6. Prepare a short update for Ash/Shuyu:
   - HTLV targeted benchmark succeeded.
   - WGS controls are clean.
   - Full WGS has only rare single-fragment HIV1 detections, so targeted HIV data is the necessary next sensitivity test.
