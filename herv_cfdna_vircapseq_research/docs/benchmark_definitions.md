# Benchmark Definitions

## Read Bins

- `unique_hiv`: read pairs that align uniquely to HIV references and not to human, HTLV, HERV, or LINE1 decoys under the selected thresholds.
- `unique_htlv`: read pairs that align uniquely to HTLV references and not to human, HIV, HERV, or LINE1 decoys under the selected thresholds.
- `unique_herv_line1`: read pairs best explained by endogenous HERV or LINE1 references.
- `ambiguous_retroviral`: read pairs aligning to conserved retroviral sequence without enough evidence to distinguish HIV, HTLV, HERV, or LINE1.
- `other_viral`: read pairs assigned to non-HIV and non-HTLV viral references if a broader panel is used.
- `human`: read pairs best explained by the human reference.

## Primary Metrics

- Healthy-control false-positive rate: high-confidence HIV or HTLV calls per million usable read pairs in expected negative controls.
- HIV sensitivity proxy: retained high-confidence HIV signal in known HIV-positive samples after filtering.
- HTLV sensitivity proxy: retained high-confidence HTLV signal in known HTLV-positive samples after filtering.
- Ambiguity burden: ambiguous retroviral reads per million usable read pairs.

## Initial High-Confidence Call Requirements

- Unique viral assignment after competitive alignment.
- Mapping quality above the frozen benchmark threshold.
- Alignment length above the frozen benchmark threshold.
- No stronger alignment to HERV or LINE1 decoys.
- For the highest confidence tier, either coverage across discriminative viral regions or host-virus junction evidence.

## First Benchmark Policy

The first benchmark should optimize specificity using healthy controls, then evaluate retained sensitivity in known HIV-positive and HTLV-positive samples. Thresholds should be frozen before final positive-cohort evaluation.
