# Command Templates

These templates define the expected workflow shape for public WGS. Replace paths and parameters before running on real data.

## Extract Unmapped Reads From BAM/CRAM

```powershell
samtools view -@ 8 -b -f 4 input.cram > sample.unmapped.bam
samtools fastq -@ 8 sample.unmapped.bam -1 sample.unmapped_R1.fastq.gz -2 sample.unmapped_R2.fastq.gz -0 sample.unmapped_singletons.fastq.gz -s sample.unmapped_singletons.fastq.gz
```

## Optional Pilot Screen On All Reads

```powershell
samtools fastq -@ 8 input.cram -1 sample.all_R1.fastq.gz -2 sample.all_R2.fastq.gz -0 sample.all_singletons.fastq.gz -s sample.all_singletons.fastq.gz
```

## Viral Alignment

```powershell
bwa mem -t 16 viral_reference_panel.fa sample.unmapped_R1.fastq.gz sample.unmapped_R2.fastq.gz > sample.viral.sam
samtools sort -@ 8 -o sample.viral.bam sample.viral.sam
samtools index sample.viral.bam
```

## Downstream Summary

Convert viral alignments into the schema documented in `docs/viral_hit_input_schema.md`, then run:

```powershell
python scripts/summarize_viral_hits.py viral_hits.csv viral_summary.csv
python scripts/normalize_viral_summary.py viral_summary.csv sample_depths.csv viral_summary_normalized.csv
python scripts/assign_viral_subtypes.py viral_summary_normalized.csv viral_subtypes.csv
```
