# Failure Modes And Next Improvements

## Failure Modes

- Conserved `gag`, `pol`, or `LTR` reads can inflate HIV or HTLV calls unless HERV and LINE1 decoys are included.
- Targeted sequencing can enrich off-target host or repetitive regions unevenly.
- Duplicate read pairs can create false viral load if not collapsed or marked.
- Short WGS fragments can support family-level signal but fail locus-level HERV interpretation.
- HTLV detection may fail if the cohort subtype is not represented in the initial reference panel.
- A low number of viral reads can be real, but without breadth or junction support it should remain low-confidence.

## Next Improvements

- Add cohort-specific HIV and HTLV references after Shuyu confirms subtype information.
- Add per-virus breadth and discriminative-region coverage thresholds.
- Add host-virus junction review for high-confidence proviral evidence.
- Add a manual review table for reads driving healthy-control false positives.
- Add final plots after real healthy-control and positive-cohort runs.
