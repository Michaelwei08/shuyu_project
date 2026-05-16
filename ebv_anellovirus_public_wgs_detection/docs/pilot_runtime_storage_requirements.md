# Pilot Runtime And Storage Requirements

## Scope

This note covers a small public-WGS pilot only. It assumes BAM/CRAM inputs, extraction of unmapped or poorly mapped reads, viral alignment against the selected herpesvirus/anellovirus panel, and downstream summary table generation.

## Pilot Subset

- Start with 1 to 3 samples per disease group when access is available.
- Prefer one EBV-annotated positive sample, one expected low-virus sample, and one sample expected to carry anellovirus signal.
- Do not stage controlled-access raw reads inside this project folder. Keep only manifests, command logs, derived nonidentifying summaries, and plots here.

## Compute Estimate

- Unmapped-read extraction from one 30x WGS BAM/CRAM: usually hours, not minutes, depending on storage throughput.
- Viral alignment of extracted reads: usually shorter than human alignment because the reference is small.
- Optional all-read screening is the expensive path and should be downsampled first.
- Use 8 to 16 threads for extraction and viral alignment during the pilot.

## Storage Estimate

- Temporary FASTQ files can be large even for unmapped-only extraction; reserve tens of GB per sample for pilot work.
- Optional all-read FASTQ extraction can approach the original WGS file size and should not be done for many samples at once.
- Final derived outputs should be small: viral hit tables, coverage summaries, subtype assignments, artifact flags, and plots.

## Runtime Log Fields

Record these fields for each pilot run:

- `sample_id`
- `input_file`
- `input_format`
- `input_size_gb`
- `screening_mode`: `unmapped_only`, `low_mapq_softclip`, `all_reads_downsampled`, or `all_reads`
- `threads`
- `wall_clock_extract_minutes`
- `wall_clock_align_minutes`
- `temporary_storage_peak_gb`
- `final_output_size_mb`
- `notes`

## Decision Rule

After the first pilot subset, compare unmapped-only screening against any downsampled all-read screen. If unmapped-only extraction misses expected EBV or anellovirus signal, expand to low-MAPQ/soft-clipped reads before attempting full all-read screening at scale.
