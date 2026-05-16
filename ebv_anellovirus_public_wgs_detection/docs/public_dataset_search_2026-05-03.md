# Public Dataset Search: EBV/Anellovirus Detection In HL, DLBCL, And PTLD

## Bottom Line

The best true human-WGS candidates are controlled-access, not immediately runnable:

- HL: EGA `EGAD00001009818` / `EGAS00001006884`.
- EBV-positive DLBCL: EGA `EGAD00001006858` / `EGAS00001004941`.
- Large DLBCL/FL WGS: EGA `EGAD00001003274` / `EGAS00001002199`.

The best immediately runnable raw FASTQ candidates are not ideal human WGS:

- Open EBV viral-genome WGS: `PRJDB18873` / `DRP013409`, including EBV samples from PTLD, DLBCL, and HL.
- Open DLBCL WES FASTQs on Dryad: `10.5061/dryad.612jm647m`.

## Practical Recommendation

1. Submit EGA access requests first for `EGAD00001009818` and `EGAD00001006858`.
2. Use `PRJDB18873` immediately to validate EBV-reference handling and EBV variant workflows, but do not use it for anellovirus or host-WGS claims.
3. Use the Dryad DLBCL WES FASTQs only as an open-code smoke test for viral-screening mechanics; exome capture biases make it unsuitable as primary WGS evidence.
4. Treat PTLD as the weakest public-WGS gap. The Dharnidharka PTLD metagenomic study is directly relevant for EBV/anellovirus, but I did not verify a public raw-data accession.

## Files Updated

- `data/cohorts/validated_public_dataset_candidates_2026-05-03.csv`
- `../ebv_lymphoma_causation_gwas_mr/data/summary_stats/validated_gwas_mr_candidate_datasets_2026-05-03.csv`
