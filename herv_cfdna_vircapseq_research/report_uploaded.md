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
