# Command Templates

These are templates. Replace paths and thread counts before running on real data.

## Extract Candidate Non-Human Reads From BAM/CRAM

```powershell
samtools view -@ 8 -b -f 4 input.cram > unmapped.bam
samtools view -@ 8 -b -q 0 input.cram > all_mapped_for_screening.bam
```

For a conservative pilot, screen both unmapped reads and all reads on a small subset, then decide whether unmapped-only screening is acceptable.

## FASTQ Competitive Alignment

```powershell
bwa mem -t 16 combined_host_virus_decoy.fa sample_R1.fastq.gz sample_R2.fastq.gz > sample.competitive.sam
samtools sort -@ 8 -o sample.competitive.bam sample.competitive.sam
samtools index sample.competitive.bam
```

## BAM/CRAM Competitive Realignment

```powershell
samtools fastq -@ 8 input_candidate_reads.bam -1 sample.candidate_R1.fastq.gz -2 sample.candidate_R2.fastq.gz -0 sample.candidate_singletons.fastq.gz -s sample.candidate_singletons.fastq.gz
bwa mem -t 16 combined_host_virus_decoy.fa sample.candidate_R1.fastq.gz sample.candidate_R2.fastq.gz > sample.candidate.competitive.sam
```

## Downstream Binning

Convert aligner output into the schema documented in `docs/competitive_alignment_input_schema.md`, then run:

```powershell
python scripts/bin_retroviral_hits.py competitive_hits.csv retroviral_bins.csv
python scripts/normalize_retroviral_bins.py retroviral_bins.csv retroviral_bins_normalized.csv --depth-csv sample_depths.csv
```
