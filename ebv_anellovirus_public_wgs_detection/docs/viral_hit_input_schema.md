# Viral Hit Input Schema

`scripts/summarize_viral_hits.py` expects one row per read-pair/virus hit.

## Required Columns

- `sample_id`
- `read_id`
- `virus`
- `family`
- `mapq`
- `alignment_length`
- `is_duplicate`

## Expected Families

- `Herpesviridae`
- `Anelloviridae`
- `Other`

## First-Pass Outputs

- Total nonduplicate viral read pairs.
- EBV read pairs.
- Other herpesvirus read pairs.
- Anellovirus read pairs.
- Per-family normalized counts if usable read-pair counts are supplied later.
