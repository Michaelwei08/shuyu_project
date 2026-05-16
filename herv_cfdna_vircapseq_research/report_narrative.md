# Distinguishing HIV-1 and HTLV-1 from HERVs and LINE1 in cfDNA, WGS, and VirCapSeq Data

Date: 2026-04-17

### 1. Can exogenous retroviruses such as HIV-1 and HTLV-1 be distinguished from HERVs and LINE1?

Yes, but only if the analysis explicitly models ambiguity instead of forcing every retrovirus-like read into a viral or endogenous bin.

- align against a combined reference containing `hg38`, HIV-1, HTLV-1, and endogenous decoys for major `HERV` and `LINE1` families
- call `exogenous` only from unique viral regions or host-virus junction reads
- call `endogenous` only when reads support human-flanked HERV or `LINE1` loci
- keep an explicit `ambiguous retroviral` bin for reads that match conserved `LTR/gag/pol` sequence but are not uniquely assignable

### 2. Can existing HERV enumeration pipelines be used or adapted for cfDNA?

Yes, but only partially.

Existing HERV DNA callers are usable for `selective insertion calling` from sufficiently deep `cfWGS`, but they are not validated drop-in `cfDNA` pipelines. The safest transfer into plasma DNA is:

- family-level `HERV/LINE1` counting
- whitelist-restricted locus analysis
- selective insertion calling only in deeper `cfWGS` samples with strong breakpoint support

### 3. How should this be handled differently for WGS versus VirCapSeq off-target reads?

They should be treated as different substrates, not as interchangeable inputs.

- `cfWGS` should be the primary dataset for endogenous `HERV/LINE1` work
- `VirCapSeq` should be used primarily for exogenous HIV-1 and HTLV-1 detection
- `VirCapSeq` off-target host reads should be treated only as a secondary exploratory substrate for family-level summaries, not as the main discovery dataset for new endogenous insertions


## Main Findings

### 1. Distinguishing exogenous retroviruses from HERVs is feasible, but only with explicit ambiguity handling

![Discrimination logic](visuals/narrative_discrimination_logic.svg)

- Virus capture approaches such as `ViroCap` sharply enrich viral reads, but are designed around viral reference genomes and probe specificity. In the Genome Research paper, targeted capture increased viral read fraction and breadth substantially, and the authors explicitly removed probe sequences with high similarity to the human genome during panel design.
- A low-stringency retroviral target-enrichment study reported that only a very small fraction of probes showed appreciable similarity to `HERV-K113`, and reads that matched the exogenous porcine retrovirus references did not remap to annotated genomic HERV loci, arguing that well-designed competitive analyses can separate exogenous from endogenous signal.
- Dedicated proviral capture approaches for HIV-1 and HTLV-1 show that capture-based sequencing can recover integrated proviral sequence, structure, and integration-site information when appropriate target molecules are present.

Implications:

- Reads mapping only to conserved retroviral motifs such as `LTR`, `gag`, or `pol` should be treated as potentially ambiguous.
- Reads mapping uniquely to HIV-1 or HTLV-1 specific regions, or reads spanning host-virus junctions, are much stronger evidence for true exogenous virus.
- A combined reference should include `hg38`, HIV-1 references, HTLV-1 references, and optionally decoy consensuses for HERV and `L1HS` to prevent forced assignment of ambiguous reads.
- A conservative operational rule is to bin retroviral reads into three classes:
  - `unique exogenous`
  - `unique endogenous`
  - `ambiguous retroviral`
- This is not an off-the-shelf tool output; it is a recommended analysis design to reduce false positives from homology.

### 2. Existing HERV DNA callers can be adapted to cfDNA, but they were built for WGS-like data, not fragmented plasma DNA

![Downloaded HERV benchmark plot](visuals/herv_tool_benchmark.svg)

#### ERVcaller

- Purpose: non-reference TE and ERV insertion detection and genotyping from short-read DNA sequencing.
- Why it matters: it is still one of the core HERV-aware DNA callers and supports custom TE reference libraries.
- Strengths:
  - Designed for DNA data rather than RNA.
  - Handles TE insertion discovery and genotyping.
  - Suitable when the goal is non-reference or polymorphic HERV discovery.
- Limitations for this project:
  - Like other short-read insertion callers, it relies on split or clipped and discordant evidence that is degraded by short cfDNA fragments and low effective molecule counts.
  - False positives remain a concern in repetitive contexts, especially without orthogonal confirmation.
- Recommendation: use as a `sensitive first-pass caller` on sufficiently deep `cfWGS`, not as a sole truth set.

#### RetroSnake / RetroSeq+

- Purpose: end-to-end `HERV-K` insertion pipeline on short-read genome data.
- Evidence: the iScience paper describes RetroSnake as an end-to-end, modular pipeline, computationally efficient at about four hours per genome, and built on a wet-lab validated protocol.
- Strengths:
  - Strong specificity-oriented design for `HERV-K`.
  - Easier operationalization than assembling a custom RetroSeq workflow from scratch.
- Limitations:
  - More HERV-K-focused than pan-HERV or pan-retroelement.
  - Still subject to the same cfDNA breakpoint-support problem as other short-read callers.
- Recommendation: use as a `specificity companion` to ERVcaller or xTea, especially if `HERV-K/HML-2` is the primary endogenous target family.

#### xTea

- Purpose: broad TE insertion detection from multiple sequencing technologies, including HERV and LINE1.
- Evidence: the Nature Communications paper states that xTea detects a wide range of retroelement insertions, including `HERV` insertions, `LINE1`, processed pseudogenes, and somatic or germline events.
- Strengths:
  - Best fit if the project needs `HERV + LINE1` in one framework.
  - Uses TE-type-specific filters and supports both germline and somatic settings.
  - Explicitly recognizes the difficulty of insertions landing inside repetitive regions and same-family repeats.
- Limitations:
  - The same paper also emphasizes that short reads miss insertions in repetitive and low-mappability regions, with long reads performing substantially better.
  - cfDNA fragmentation is likely to exacerbate exactly that failure mode.
- Recommendation: make xTea the `primary all-in-one TE caller` if LINE1 is as important as HERV.

#### STEAK and MELT

- `STEAK` is useful historically and can recover reference HERV contexts, but it is not the strongest first choice for de novo cfDNA insertion discovery.
- `MELT` is mature for population-scale TE work and is relevant for `LINE1`, `Alu`, and `SVA`, but it is less central than xTea for a combined HERV plus LINE1 question.
- Recommendation: secondary comparators only.

### 3. Activity-oriented extensions are secondary to the main enumeration question

#### Telescope

- Telescope is a locus-specific RNA-seq quantifier that reassigns ambiguous fragments with a Bayesian model and reports expression at specific TE insertions rather than only at subfamily level.
- This is directly relevant if the project later adds `cfRNA`, PBMC RNA-seq, or tissue RNA-seq.
- It is not a DNA caller and should not be presented as a cfDNA HERV enumeration tool.

#### ERVmap

- ERVmap focuses on locus-specific proviral ERV transcription analysis from RNA-seq using a curated proviral annotation set.
- It is useful as a stringent locus whitelist or annotation resource, not as a direct cfDNA solution.

Inference:

- `ERVmap` and `Telescope` annotations can be repurposed as candidate locus coordinate sets for cfDNA methylation or fragmentomic analysis at HERV promoter or `LTR` regions.
- That is an adaptation idea, not something already validated as a standard cfDNA pipeline.

### 4. WGS is preferable to VirCapSeq off-target reads for endogenous retroelement work

Evidence-backed points:

- Virus-capture sequencing was developed to enrich viral sequence, not to provide unbiased host-genome representation.
- In `ViroCap`, probes with high similarity to human sequence were removed during design, and capture strongly increased viral read fraction. This is desirable for virome analysis but means the host background is shaped by capture design and not equivalent to `WGS`.
- Off-target reads from hybrid capture panels can support broad genome-wide analyses such as copy-number reconstruction. Tools and studies in that space show the concept is valid for coarse signals, not that off-target reads are equivalent to unbiased WGS for repetitive-locus insertion discovery.

Implications:

- `VirCapSeq on-target reads`:
  - Good for HIV-1 or HTLV-1 detection if the panel captures them.
  - Potentially useful for proviral characterization or integration evidence if host-virus chimeric reads are recovered.
- `VirCapSeq off-target host reads`:
  - Reasonable only for exploratory family-level repeat burden, broad fragmentomics, or very coarse locus-of-interest coverage.
  - Poor substrate for de novo HERV insertion calling because coverage is sparse, nonuniform, and not designed around repetitive human loci.
- If the goal is endogenous retroelement discovery or locus enumeration, start from `cfWGS`.
- If the goal is exogenous HIV-1 or HTLV-1 detection, `VirCapSeq` may still be very useful, but it answers a different part of the problem.

### 5. The most realistic cfDNA adaptation is not full locus discovery everywhere, but tiered analysis

Recommended tiering:

1. `Exogenous virus detection`
   - Competitive alignment against `hg38 + HIV-1 + HTLV-1`.
   - Require unique viral evidence or host-virus junction reads for high-confidence calls.
   - Keep an explicit ambiguous-retroviral bin.

2. `Family-level endogenous retroelement burden`
   - Use `RepeatMasker` or `Dfam` annotations to summarize reads overlapping `LTR/ERV`, `HERV-K`, and `LINE1` categories after strict preprocessing.
   - This is the most plausible entry point for both `cfWGS` and exploratory off-target host reads.

3. `Selective insertion calling`
   - Run `xTea` and either `ERVcaller` or `RetroSnake` on deeper `cfWGS` samples only.
   - Restrict interpretation to candidates with strong breakpoint support, known polymorphic loci, or orthogonal confirmation.

4. `Proxy activity analysis`
   - If biological interest is activation rather than just abundance, add:
     - methylation analysis at `LTR` or `LINE1` loci
     - cfDNA fragmentomics at HERV `LTR` promoter regions using region-of-interest frameworks such as `LIQUORICE` or `Griffin`
     - or direct `cfRNA` or RNA-seq with `Telescope` or `ERVmap`

## Proposed Analysis Blueprint

### Control design using the known infected samples

The known `HIV-1+` and `HTLV-1+` samples should be used as more than just illustrative examples. They should anchor the pilot design.

![Control calibration design](visuals/narrative_control_calibration.svg)

1. Use them to verify that the combined reference correctly assigns reads into `unique exogenous`, `unique endogenous`, and `ambiguous retroviral` bins.
2. Use them to identify which HIV-1 and HTLV-1 regions remain most discriminative after short-fragment sequencing and alignment filtering.
3. Use them to measure false assignment into `HERV` or `LINE1` decoys under the exact preprocessing and alignment settings planned for cfDNA.
4. Use them as positive controls when deciding whether any apparent endogenous signal in `VirCapSeq` off-target reads is real biology or capture-bias background.

![Recommended workflow](visuals/analysis_workflow.svg)

### A. Primary path from cfWGS

1. Preprocess
   - adapter trimming
   - duplicate marking or collapse
   - strict contamination control
   - alignment to a combined host plus virus reference

2. Exogenous retrovirus track
   - call HIV-1 and HTLV-1 from unique viral alignments
   - inspect coverage breadth across discriminative viral regions
   - recover host-virus chimeric reads where present

3. Endogenous retroelement track
   - summarize `HERV`, `LTR`, and `LINE1` family burden with `RepeatMasker` or `Dfam`
   - run `xTea` for broad TE discovery
   - run `ERVcaller` or `RetroSnake` as a HERV-focused companion

4. Candidate prioritization
   - prefer loci with unique flanking support
   - cross-check against known HERV-K catalogs and dbRIP-style references where available
   - deprioritize calls supported only by conserved retroviral sequence without flank support

5. Validation
  - targeted PCR, ddPCR, or locus-specific capture
  - matched positive controls from known HIV-1 or HTLV-1 infected subjects
  - negative plasma controls and no-template controls

### B. Secondary path from VirCapSeq data

1. Use on-target reads for exogenous virus detection and characterization.
2. Quantify the deduplicated off-target host yield before committing to endogenous analysis.
3. If off-target depth is modest, restrict endogenous analysis to:
   - family-level repeat burden
   - coarse region-of-interest coverage
   - exploratory fragmentomic summaries
4. Do not use off-target VirCapSeq as the primary discovery dataset for new HERV insertions.

## Recommended Tool Roles

| Tool or resource | Best use here | Fit for cfDNA | Fit for VirCapSeq off-target | Recommendation |
| --- | --- | --- | --- | --- |
| `xTea` | Combined `HERV + LINE1` insertion calling from WGS | Moderate if cfWGS is deep enough | Low | Primary TE caller |
| `ERVcaller` | Sensitive HERV-aware DNA insertion calling | Moderate with caution | Low | Sensitive companion caller |
| `RetroSnake` | Specific `HERV-K` insertion workflow | Moderate with caution | Low | Specificity companion |
| `MELT` | Secondary `LINE1/Alu/SVA` comparison | Moderate | Low | Optional comparator |
| `STEAK` | Reference-context HERV support | Low to moderate | Low | Not first-line |
| `RepeatMasker` / `Dfam` | Family and locus annotation framework | High | Moderate for exploratory summaries | Core annotation resource |
| `Telescope` | RNA-based locus-specific HERV expression | Not applicable to DNA | Not applicable | Reserve for cfRNA or RNA-seq |
| `ERVmap` | RNA-based stringent proviral locus annotation | Not applicable to DNA directly | Not applicable | Use as locus whitelist only |
| `LIQUORICE` / `Griffin` | Region-of-interest cfDNA fragmentomics | Moderate | Moderate if off-target depth allows | Exploratory activity proxy |
| `VirusPredictor` | Sequence-level discrimination of infectious virus vs human ERV | Exploratory | Exploratory | Optional triage, not primary caller |

## Decision Summary

- Best assay for endogenous retroelement discovery: `cfWGS`
- Best assay for exogenous virus detection: `VirCapSeq` or related capture, if HIV-1 and HTLV-1 are in panel scope
- Best combined DNA caller when LINE1 matters: `xTea`
- Best HERV-focused companion: `ERVcaller` or `RetroSnake`
- Best way to study activity: add `cfRNA`, methylation, or fragmentomics rather than relying on cfDNA sequence alone

## Immediate Practical Next Steps

1. Assemble a pilot cohort with known `HIV-1+`, `HTLV-1+`, and negative plasma samples, and use it as a calibration set rather than only as a validation set.
2. Build a combined reference containing `hg38`, HIV-1 references, HTLV-1 references, and decoy consensuses for major HERV and `LINE1` families.
3. Define a three-bin classification rule up front: `unique exogenous`, `unique endogenous`, and `ambiguous retroviral`.
4. Run the exogenous discrimination track first on both `cfWGS` and any available `VirCapSeq` data to quantify how much viral-versus-endogenous misassignment remains after competitive alignment.
5. Run `xTea` plus `ERVcaller` or `RetroSnake` on the deeper `cfWGS` samples only, and treat those calls as the primary endogenous enumeration output.
6. Use `VirCapSeq` off-target host reads only if deduplicated host yield is sufficient for coarse endogenous summaries; do not force de novo HERV insertion discovery from that assay.

## Key Sources

- ERVcaller paper: https://academic.oup.com/bioinformatics/article/35/20/3913/5416145
- ERVcaller repository: https://github.com/xunchen85/ERVcaller
- RetroSnake paper: https://www.sciencedirect.com/science/article/pii/S2589004222015619
- RetroSnake repository: https://github.com/KHP-Informatics/RetroSnake
- xTea paper: https://www.nature.com/articles/s41467-021-24041-8
- xTea repository: https://github.com/parklab/xTea
- Telescope paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453
- Telescope repository: https://github.com/mlbendall/telescope
- ERVmap project: https://mtokuyama.github.io/ERVmap/
- Dfam: https://dfam.org/home
- RepeatMasker: https://www.repeatmasker.org/
- Enhanced virome sequencing through solution-based capture (ViroCap): https://genome.cshlp.org/content/early/2015/09/22/gr.191049.115.full.pdf
- Investigation of Human Cancers for Retrovirus by Low-Stringency Target Enrichment and High-Throughput Sequencing: https://pmc.ncbi.nlm.nih.gov/articles/PMC4541070/
- HIV-1 DNA-capture-seq paper: https://www.nature.com/articles/s41598-019-48681-5
- HTLV-1 viral DNA-capture-seq paper: https://doi.org/10.1016/j.celrep.2019.09.016
- SavvyCNV off-target copy-number analysis: https://pmc.ncbi.nlm.nih.gov/articles/PMC8959187/
- LIQUORICE note: framework mentioned in the literature review, but I did not validate a stable primary paper link in this pass
- Griffin paper: https://www.nature.com/articles/s41467-022-35076-w
- VirusPredictor paper: https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btae192/7643508
