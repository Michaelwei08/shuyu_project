# Competitive Alignment Input Schema

`scripts/bin_retroviral_hits.py` expects a CSV file with one row per read-pair/reference hit.

## Required Columns

- `sample_id`
- `read_id`
- `reference_category`
- `mapq`
- `alignment_score`
- `alignment_length`

## Allowed Reference Categories

- `human`
- `hiv`
- `htlv`
- `herv`
- `line1`
- `other_viral`

## Output Bins

- `unique_hiv`
- `unique_htlv`
- `unique_herv_line1`
- `ambiguous_retroviral`
- `other_viral`
- `human`
- `unassigned`

## Notes

The script is intentionally downstream of alignment. It does not choose the aligner. It standardizes how read-level competitive hits are converted into benchmarkable sample-level bins.
