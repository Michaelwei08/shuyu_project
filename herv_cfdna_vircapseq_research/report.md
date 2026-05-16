# Distinguishing exogenous retroviruses from HERVs and LINE1 in cfDNA, WGS, and VirCapSeq data

Generated from the deep-research JSON outputs using the `research-report` workflow and revised to match the email's decision-focused framing.

## Direct Answer to the Email

This report is organized around the three concrete questions in the email, with the item-by-item deep-research summaries moved below as supporting material.

### 1. Can HIV-1 and HTLV-1 be distinguished from HERVs and LINE1?

Yes, but only with an explicit ambiguity-aware discrimination framework. The recommended rule is to align against a combined reference containing `hg38`, HIV-1, HTLV-1, and endogenous `HERV/LINE1` decoys, call `exogenous` only from unique viral regions or host-virus junctions, call `endogenous` only from human-flanked loci, and keep a third `ambiguous retroviral` bin for conserved retroviral sequence.

### 2. Can existing HERV enumeration pipelines be adapted to cfDNA?

Yes, but only partially. Existing HERV DNA callers can support selective insertion calling from sufficiently deep `cfWGS`, but the safest cfDNA transfer is family-level counting plus whitelist-restricted locus analysis. They should not be described as validated drop-in cfDNA pipelines.

### 3. How should this be handled for WGS versus VirCapSeq off-target reads?

Treat them as different substrates. `cfWGS` should be the primary dataset for endogenous `HERV/LINE1` work. `VirCapSeq` should be used primarily for exogenous HIV-1 and HTLV-1 detection, while off-target host reads should be limited to secondary exploratory summaries rather than used as the main discovery substrate for new endogenous insertions.

## Recommended Pilot Plan

1. Use the known `HIV-1+`, `HTLV-1+`, and negative samples as a calibration cohort, not just as retrospective validation material.
2. Build a combined host-plus-virus reference and predefine the three-bin classification rule: `unique exogenous`, `unique endogenous`, and `ambiguous retroviral`.
3. Run the exogenous discrimination track first on all available `cfWGS` and `VirCapSeq` data to quantify viral-versus-endogenous misassignment under the real assay conditions.
4. Run endogenous enumeration primarily on `cfWGS`, using family-level repeat counting plus `xTea` and a HERV-focused companion caller such as `ERVcaller` or `RetroSnake` on the deeper samples.
5. Use `VirCapSeq` off-target host reads only if deduplicated host yield is adequate for coarse endogenous summaries; otherwise stop at the exogenous-retrovirus result.

The downloaded benchmark supports the tool-specific emphasis in this plan: `ERVcaller` has the highest sensitivity in the extracted table at `0.77`, which is why the report favors `ERVcaller`-style sensitivity plus companion-caller confirmation over a single rigid caller choice.

`cfRNA`, methylation, and fragmentomics remain reasonable follow-on extensions if the question later shifts from enumeration to activity, but they are secondary to the present email's main request.

## Validation Snapshot

- Sources audited: 53
- HTTP OK: 49
- Sources with warnings: 16
- Narrative report: [report_narrative.md](/D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/report_narrative.md)
- Claim validation: [claim_validation.md](/D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/validation/claim_validation.md)
- Source audit: [source_validation.md](/D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/validation/source_validation.md)

## Visual Overview

### Assay Fit by Question

![Assay fit heatmap](visuals/assay_fit_heatmap.svg)

### Recommended Workflow

![Analysis workflow](visuals/analysis_workflow.svg)

### Validation Dashboard

![Validation dashboard](visuals/validation_dashboard.svg)

## Downloaded Benchmark Support

I downloaded the benchmark repository linked to the HERV tool-assessment paper and extracted Supplementary Table 3 into [herv_tool_benchmark_table_3.csv](/D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/data/herv_tool_benchmark_table_3.csv).

This benchmark supports the report's tool-specific ranking, especially the decision to emphasize `ERVcaller` and to treat `STEAK` as a secondary comparator rather than a sole discovery tool.

![Downloaded HERV benchmark plot](visuals/herv_tool_benchmark.svg)

- Highest sensitivity in the downloaded table: `ERVcaller` at `0.77`
- Highest precision in the downloaded table: `Steak` at `0.90`
- Highest long-read confirmation among single-value rows: `Retroseq+` at `78%`
- Class averages are mixed: HERV-specific mean sensitivity `0.48` vs generalist `0.54`, largely because `STEAK` is precision-heavy but sensitivity-poor
- HERV-specific mean precision remains slightly higher: `0.81` vs generalist `0.79`

The strongest support here is tool-specific rather than class-average: `ERVcaller` has the best balance in this benchmark, `STEAK` fits the report's characterization as a secondary high-precision comparator, and `Retroseq+` contributes confirmation-style value rather than leading on sensitivity or precision.

This downloaded benchmark does not by itself establish cfDNA suitability, which remains an inferential extension from WGS-focused benchmarking plus cfDNA assay constraints.

## Table of Contents

1. [cfDNA adaptation paths](#cfdna-adaptation-paths) - Category: Method question | Recommended Role: Tiered secondary analysis path: start with family-level counting, then move to whi...
2. [Exogenous retrovirus discrimination strategy](#exogenous-retrovirus-discrimination-strategy) - Category: Analysis question | Recommended Role: Primary decision framework for classifying retrovirus-like reads in legacy cfDNA o...
3. [HERV DNA insertion callers](#herv-dna-insertion-callers) - Category: Tool class | Recommended Role: Use xTea as the primary combined HERV plus LINE1 DNA caller when both classes matt...
4. [HERV locus and annotation resources](#herv-locus-and-annotation-resources) - Category: Reference resource | Recommended Role: Use RepeatMasker or Dfam as the core annotation layer, add curated HERV-K catalogs...
5. [LINE1 and broader TE callers](#line1-and-broader-te-callers) - Category: Tool class | Recommended Role: xTea: primary TE caller for the combined HERV plus LINE1 analysis. MELT: secondary...
6. [Validation strategy](#validation-strategy) - Category: Experimental design | Recommended Role: Primary post-caller validation framework for all candidate retrovirus-like, HERV,...
7. [VirCapSeq off-target reuse](#vircapseq-off-target-reuse) - Category: Method question | Recommended Role: Exploratory secondary analysis only: use for family-level repeat burden or backgro...

## cfDNA adaptation paths

### Framing

- **Name:** cfDNA adaptation paths
- **Category:** Method question
- **Purpose:** > Assess which parts of existing HERV and TE pipelines can be repurposed for cfDNA analysis, with emphasis on family-level counting, locus whitelists, methylation, and fragmentomics rather than de novo insertion discovery.
### Relevant Assays

- cfWGS
- low-pass cfWGS
- VirCapSeq off-target host reads
- cfDNA methylation sequencing
- cfDNA fragmentomics
- RNA-seq or cfRNA for whitelist derivation and orthogonal validation

### Technical Evaluation

### Input Requirements

- Aligned cfDNA reads with duplicate handling and contamination control
- RepeatMasker or Dfam family annotations, or equivalent repeat-family coordinate tables
- Curated HERV, LTR, and LINE1 locus whitelist sets when moving beyond family-level summaries
- Enough depth for family-level aggregation or ROI-level fragmentomic analysis
- Bisulfite, enzymatic conversion, or equivalent methylation-readout data if methylation is part of the question
- Fragment-size, end-motif, and genomic-position summaries if fragmentomics is part of the question

- **Output Resolution:** > Best at family-level abundance, locus-whitelist-restricted signal, methylation state, and fragmentomic enrichment; not reliable for genome-wide de novo insertion discovery in fragmented plasma DNA.
- **Strengths:** > Family-level counting transfers cleanly because multi-mapping reads can be aggregated instead of forced into a single locus, and repeat-landscape work in cfDNA shows that repeat families can be measured in plasma. RepeatMasker and Dfam already provide repeat-family coordinate systems, while Telescope and ERVmap annotations can be reused as curated locus whitelists. cfDNA methylation and fragmentomics add orthogonal proxy signals when breakpoint evidence is weak, which is a better fit for plasma DNA than full insertion calling.
- **Limitations:** > De novo insertion discovery is still weak in cfDNA because split, clipped, and discordant evidence is sparse. Locus-level claims need a whitelist and unique flanks, otherwise ambiguous retroviral mapping dominates. Methylation and fragmentomics are indirect activity proxies rather than direct insertion or transcription evidence. VirCapSeq off-target host reads are usually too sparse and capture-biased for broad discovery.

### Recommendation

- **Recommended Role:** > Tiered secondary analysis path: start with family-level counting, then move to whitelist-restricted loci, then methylation or fragmentomics, while reserving de novo insertion calling for deep cfWGS only.
- **Validation Needs:** > Use targeted PCR or ddPCR for any claimed new insertion, confirm methylation with bisulfite, enzymatic conversion, or cfTAPS-style chemistry on matched samples, benchmark whitelist loci against known HERV and LINE1 catalogs, and verify fragmentomic signals in independent cfWGS or matched orthogonal datasets.
### Key Sources

- https://academic.oup.com/bioinformatics/article/31/22/3593/240793
- https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6294949/
- https://www.repeatmasker.org/dev/
- https://www.repeatmasker.org/libraries/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4702899/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11323656/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6117241/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11684646/
- https://www.nature.com/articles/s41467-021-23445-w
- https://www.nature.com/articles/s41467-022-35076-w
- https://academic.oup.com/bioinformatics/article/35/20/3913/5416145
- https://www.nature.com/articles/s41467-021-24041-8

- **Notes:** > The safest cfDNA transfer is tiered. Family-level burden is the most defensible adaptation, locus whitelists are best treated as curated candidate sets rather than direct cfDNA outputs, and methylation or fragmentomics should be treated as proxy activity layers. RNA-centric TE tools such as Telescope and ERVmap are useful mainly as annotation sources and whitelist generators, not as direct cfDNA callers.

## Exogenous retrovirus discrimination strategy

### Framing

- **Name:** Exogenous retrovirus discrimination strategy
- **Category:** Analysis question
- **Purpose:** > Build a comparative bioinformatics workflow that separates reads from known exogenous human retroviruses from human endogenous retroviruses and LINE1-derived reads in existing cfDNA or capture-based sequencing datasets.
- **Relevant Assays:** cfWGS, VirCapSeq, targeted viral capture sequencing, shotgun WGS

### Technical Evaluation

### Input Requirements

- Short-read FASTQ, BAM, or CRAM, ideally with paired-end read information preserved
- A curated exogenous retrovirus reference panel that includes representative HIV, HTLV, and other known human-infective retrovirus genomes
- A human reference genome plus repeat-aware background references for HERV and LINE1 loci or consensus sequences
- Optional local assemblies or realigned candidate contigs for ambiguous retrovirus-like reads

- **Output Resolution:** > Read-level taxonomic assignment, contig-level confirmation, and optional locus-level integration assessment when human-viral junctions are present.
### Strengths

- Fits legacy cfDNA and capture datasets because it relies on comparative alignment and classification rather than de novo discovery.
- Uses explicit HERV and LINE1 decoy references so exogenous calls are made against the real endogenous background instead of by viral similarity alone.
- Can exploit diagnostic viral regions such as env, LTR structure, and other clade-specific sequence outside conserved reverse-transcriptase motifs.
- Supports both free viral fragments and proviral junction evidence if those reads are present in the dataset.
- Can be implemented with multiple independent classifiers and a consensus rule set to reduce single-tool bias.

### Limitations

- Short cfDNA fragments and capture bias often leave only conserved retroviral segments, which can remain ambiguous between exogenous retroviruses and HERVs.
- LINE1-derived reads can look retrovirus-like in reverse-transcriptase-rich regions unless non-LTR architecture and flanking context are checked.
- Highly repetitive HERV loci and multi-mapping reads limit strain-level claims from short-read data alone.
- Capture-based datasets may enrich target viral reads while leaving too little host background for confident endogenous-versus-exogenous locus interpretation.

### Recommendation

- **Recommended Role:** > Primary decision framework for classifying retrovirus-like reads in legacy cfDNA or capture data; use as a conservative filter that reports only exogenous calls with multi-evidence support.
### Validation Needs

- Re-align candidates against a combined database containing human, HERV or LTR, LINE1, and exogenous retrovirus references.
- Inspect mapping quality, multi-mapping behavior, and unique k-mer support for each candidate read or contig.
- Require confirmation from virus-specific regions or phylogenetic placement, not only reverse-transcriptase or protease fragments.
- Downgrade any call that is equally compatible with a known HERV locus or a LINE1-derived sequence.

### Key Sources

- https://pmc.ncbi.nlm.nih.gov/articles/PMC8669383/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC7093845/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9605557/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9650475/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4665012/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4541070/
- https://pubmed.ncbi.nlm.nih.gov/17636050/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC5131823/

- **Notes:** > Retroviruses and HERVs share gag-pol-env and LTR architecture, so the main discriminants are comparative: exact placement against curated human repeat panels, diagnostic sequence outside conserved pol, and locus context. LINE1 is non-LTR and is better excluded by its 5' UTR-ORF1-ORF2-3' UTR-poly(A) architecture, plus target-site duplication and flanking evidence. For conservative reporting, label unresolved short-pol fragments as retrovirus-like rather than exogenous.

## HERV DNA insertion callers

### Framing

- **Name:** HERV DNA insertion callers
- **Category:** Tool class
- **Purpose:** > Compare the main short-read DNA callers relevant to non-reference or polymorphic HERV detection and assign practical roles for cfWGS and related DNA sequencing studies.
### Relevant Assays

- cfWGS
- tissue WGS
- high-depth short-read DNA sequencing
- targeted DNA sequencing with recoverable host flanks

### Technical Evaluation

### Input Requirements

- Short-read FASTQ, BAM, or CRAM with split, soft-clipped, or discordant evidence preserved
- Human reference genome aligned to the same build used by downstream annotations
- TE or HERV reference library, often including HERV-K or broader repeat consensuses
- Enough depth to recover breakpoint-supporting evidence, especially for cfDNA

- **Output Resolution:** > Insertion-level or locus-level calls for non-reference or polymorphic HERV events, often with breakpoint refinement, genotype estimation, and supporting-read evidence.
### Strengths

- This caller class is designed for DNA rather than RNA and therefore maps cleanly onto the insertion-enumeration question.
- ERVcaller remains a strong HERV-aware sensitive caller with custom-library support and genotype output.
- RetroSnake adds an end-to-end, HERV-K-focused, specificity-oriented workflow with annotation and reporting around RetroSeq-style discovery.
- xTea is the broadest option when HERV and LINE1 need to be studied together in one pipeline.
- STEAK remains useful as a historical comparator and for some reference-context HERV analyses.

### Limitations

- All tools in this class depend on breakpoint-supporting split, clipped, or discordant reads, which are weakened by cfDNA fragmentation and low molecule counts.
- Sparse off-target host reads from VirCapSeq are generally a poor substrate for de novo HERV insertion discovery.
- Short-read mapping remains fragile in highly repetitive or same-family insertion contexts.
- RetroSnake is more HERV-K-focused than pan-HERV, and STEAK is generally less sensitive than the strongest current alternatives.

- **Cfdna Fit:** > Moderate for deep cfWGS with adequate breakpoint evidence; low for sparse plasma datasets or low-depth sequencing.
- **Distinction Value:** > Moderate, because these tools define endogenous retroviral background and candidate HERV loci but do not directly classify exogenous retrovirus reads.

### Recommendation

- **Recommended Role:** > Use xTea as the primary combined HERV plus LINE1 DNA caller when both classes matter, pair it with ERVcaller or RetroSnake for HERV-focused confirmation, and treat STEAK as a secondary comparator rather than a first-line discovery tool.
### Validation Needs

- Orthogonal PCR or ddPCR for candidate loci
- Manual review of breakpoint-supporting reads and flanking uniqueness
- Matched positive and negative controls
- Escalation to targeted resequencing or long-read confirmation for ambiguous repetitive loci

### Key Sources

- https://academic.oup.com/bioinformatics/article/35/20/3913/5416145
- https://github.com/xunchen85/ERVcaller
- https://www.sciencedirect.com/science/article/pii/S2589004222015619
- https://github.com/KHP-Informatics/RetroSnake
- https://www.nature.com/articles/s41467-021-24041-8
- https://github.com/parklab/xTea
- https://pmc.ncbi.nlm.nih.gov/articles/PMC5597868/
- https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2022.1062328/full

- **Notes:** > Practical ranking for this project: xTea is best when LINE1 is co-primary with HERV; ERVcaller is a sensitive HERV-aware companion; RetroSnake is useful when HERV-K specificity and pipeline completeness matter; STEAK is better treated as a benchmark comparator than a main caller. None of these should be used as a direct exogenous-virus classifier.

## HERV locus and annotation resources

### Framing

- **Name:** HERV locus and annotation resources
- **Category:** Reference resource
- **Purpose:** > Compare the main coordinate and annotation systems that can anchor HERV-related cfDNA analyses, including family-level repeat tracks, curated HERV-K catalogs, and RNA-derived locus whitelists.
- **Relevant Assays:** cfWGS, bulk WGS, tissue WGS, VirCapSeq off-target host reads, RNA-seq, cfRNA

### Technical Evaluation

### Input Requirements

- A consistent human reference genome build such as hg38 or T2T-CHM13
- RepeatMasker or Dfam annotations for family and locus overlap analysis
- Curated HERV-K or HML-2 locus catalogs when specific loci are in scope
- Optional ERVmap or Telescope annotation files when building conservative locus whitelists

- **Output Resolution:** > Ranges from family-level repeat annotation to curated locus-level whitelist coordinates, depending on the resource used; these resources annotate existing loci rather than discover new insertions by themselves.
### Strengths

- RepeatMasker and Dfam provide the most mature base coordinate layer for family-level HERV, LTR, and LINE1 overlap analysis.
- Curated HERV-K catalogs are the best option when the question narrows to known polymorphic or recently integrated HML-2 loci.
- ERVmap offers a conservative proviral whitelist that can be repurposed for high-specificity locus restriction.
- Telescope annotations are valuable if the project later adds cfRNA or tissue RNA-seq and wants locus-specific expression follow-up.

### Limitations

- None of these resources is a direct cfDNA discovery caller.
- RNA-derived resources such as ERVmap and Telescope are being repurposed when used on cfDNA, which is an inference rather than a validated standard.
- Catalog completeness varies across genome builds, populations, and HERV definitions.
- Sparse off-target host reads usually limit these resources to coarse overlap summaries rather than strong locus claims.

- **Cfdna Fit:** > High for RepeatMasker or Dfam family-level annotation, moderate for curated HERV-K locus whitelists, and low to moderate for repurposed RNA-derived whitelists such as ERVmap or Telescope.
- **Distinction Value:** > Moderate to high, because these resources define endogenous retroelement background and help prevent ambiguous retroviral reads from being overcalled as exogenous signal.

### Recommendation

- **Recommended Role:** > Use RepeatMasker or Dfam as the core annotation layer, add curated HERV-K catalogs for known polymorphic-locus tracking, and use ERVmap or Telescope only as conservative whitelist resources rather than primary cfDNA coordinate systems.
### Validation Needs

- Use the same genome build across alignment, annotation, and locus interpretation
- Cross-check catalog coordinates against the specific source version
- Require unique flanking evidence or orthogonal assays for any locus-level claim
- Treat RNA-derived whitelists as exploratory on cfDNA unless independently confirmed

### Key Sources

- https://www.repeatmasker.org/RepeatMasker/
- https://www.dfam.org/home
- https://pubmed.ncbi.nlm.nih.gov/23203985/
- https://pubmed.ncbi.nlm.nih.gov/26612867/
- https://pubmed.ncbi.nlm.nih.gov/32375827/
- https://mtokuyama.github.io/ERVmap/
- https://www.ervmap.com/
- https://pubmed.ncbi.nlm.nih.gov/30455304/
- https://github.com/mlbendall/telescope
- https://github.com/mlbendall/telescope_annotation_db

- **Notes:** > Resource ranking for this project: RepeatMasker or Dfam should anchor all family-level summaries; curated HERV-K catalogs are the best locus-specific whitelist for endogenous HML-2 questions; ERVmap is a conservative proviral whitelist; Telescope annotations are most defensible for RNA-based follow-up and only secondarily as exploratory cfDNA whitelists.

## LINE1 and broader TE callers

### Framing

- **Name:** LINE1 and broader TE callers
- **Category:** Tool class
- **Purpose:** > Evaluate whether xTea and MELT are suitable for combined HERV plus LINE1 insertion discovery, and assign each tool a role in a host-genome transposable element analysis workflow.
- **Relevant Assays:** > cfWGS, bulk WGS, and optionally high-depth WES for known MEI comparison; not a primary fit for VirCapSeq off-target reads.

### Technical Evaluation

- **Input Requirements:** > xTea needs aligned short-read or long-read sequencing, with WGS, WES, or hybrid data and repeat annotations; MELT primarily needs Illumina paired-end WGS BAM or CRAM with discordant-read and split-read evidence, and is most mature for population-scale cohorts.
- **Output Resolution:** > Locus-level insertion calls with breakpoint evidence and genotypes; xTea can annotate HERV, LINE1, Alu, SVA, processed pseudogenes, and other TE insertions, while MELT focuses on Alu, LINE1, and SVA MEIs.
- **Strengths:** > xTea is the broader caller: it supports multiple sequencing technologies, can detect HERV and LINE1 in one framework, uses TE-type-specific filters, and can genotype calls. MELT is mature, scalable, and well established for population-scale LINE1, Alu, and SVA calling, with strong use of discordant-read and split-read evidence.
- **Limitations:** > xTea still depends on breakpoint-supporting reads and therefore loses sensitivity in low-input, fragmented, or highly repetitive regions; cfDNA is a hard use case. MELT is narrower, with no core HERV calling support, and is optimized for Illumina paired-end WGS rather than fragmented cfDNA or sparse capture off-target data.

### Recommendation

- **Recommended Role:** > xTea: primary TE caller for the combined HERV plus LINE1 analysis. MELT: secondary comparator and known LINE1, Alu, and SVA genotyper, not the primary HERV caller.
- **Validation Needs:** > Confirm candidates with targeted PCR or ddPCR, inspect breakpoint-supporting reads and flanking uniqueness, require orthogonal support for repetitive loci, and compare calls against known positive controls and matched negatives.
### Key Sources

- https://github.com/parklab/xTea
- https://www.nature.com/articles/s41467-021-24041-8
- https://melt.igs.umaryland.edu/
- https://genome.cshlp.org/content/early/2017/08/30/gr.218032.116.full.pdf
- https://www.nature.com/articles/s41431-023-01478-7

- **Notes:** > If the project later restricts to LINE1-only analysis, MELT becomes more useful as a comparator; if HERV remains in scope, xTea should stay the main caller.

## Validation strategy

### Framing

- **Name:** Validation strategy
- **Category:** Experimental design
- **Purpose:** > Define conservative analytical controls and orthogonal confirmation standards for candidate exogenous retrovirus, HERV, and LINE1 findings in sequencing-based studies so that homology, capture bias, library artifacts, and multi-mapping are not mistaken for biological signal.
- **Relevant Assays:** cfWGS, VirCapSeq, shotgun WGS, targeted viral capture sequencing, long-read WGS, RNA-seq

### Technical Evaluation

### Input Requirements

- Candidate calls with supporting read-level evidence, locus coordinates, and mapping-quality context
- A combined human, exogenous retrovirus, HERV or LTR, and LINE1 reference or decoy panel
- Matched negative controls, blank libraries, and batch-matched non-candidate samples
- Optional orthogonal data such as targeted PCR, ddPCR, Sanger sequencing, or long-read confirmation
- Sample and library metadata sufficient to track batch effects, index hopping, and contamination risk

- **Output Resolution:** > Candidate-level accept or reject triage with locus-level confidence tiers and a documented orthogonal confirmation status.
### Strengths

- Separates true biological signal from capture, mapping, and library-preparation artifacts before interpretation.
- Uses evidence classes tailored to the target, including unique viral sequence, host-virus junctions, insertion breakpoints, and flanking host context.
- Applies across cfWGS and capture-based datasets because the logic is assay-agnostic.
- Creates an explicit hierarchy of negative controls, technical replicates, and decoy-aware re-alignment checks.
- Reduces overcalling by forcing different confirmation standards for exogenous retroviruses, HERVs, and LINE1.

### Limitations

- Cannot recover biological signal that the assay never captured, so poor substrate remains poor.
- Strict criteria can increase false negatives in low-input cfDNA or sparse off-target capture data.
- Orthogonal confirmation may be infeasible for rare mosaic events or highly repetitive loci.
- Family-level standards are not sufficient for locus-specific claims, and locus-specific standards do not scale to genome-wide screening.
- Performance depends on the completeness of the reference and decoy set used for adjudication.

### Recommendation

- **Recommended Role:** > Primary post-caller validation framework for all candidate retrovirus-like, HERV, and LINE1 calls before any biological interpretation or downstream reporting.
### Validation Needs

- Require at least one orthogonal confirmation channel for high-impact calls, such as locus-specific PCR, ddPCR, Sanger sequencing, or long-read sequencing.
- Retain matched negative controls, blank libraries, and batch-matched non-candidate samples to estimate assay-specific false-positive rates.
- Demand concordance across technical replicates or independent library preparations for borderline calls.
- For exogenous retroviruses, require unique virus-specific sequence or host-virus junction evidence rather than conserved pol or LTR fragments alone.
- For HERV and LINE1, require insertion hallmarks such as breakpoint support, flanking host sequence, target-site duplication, or polyA signal when applicable.
- Use manual review of read placement, multi-mapping, soft clips, and decoy competition before accepting a candidate as confirmed.

### Key Sources

- https://www.nature.com/articles/s41467-021-24041-8
- https://academic.oup.com/bioinformatics/article/35/20/3913/5416145
- https://www.nature.com/articles/s41598-019-48681-5.pdf
- https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02307-0
- https://doi.org/10.1016/j.celrep.2019.09.016
- https://link.springer.com/article/10.1186/s12977-020-00519-z
- https://www.sciencedirect.com/science/article/pii/S2589004222015619
- https://www.nature.com/articles/srep33598
- https://pubmed.ncbi.nlm.nih.gov/28655292/

- **Notes:** > The validation logic should be evidence-stratified rather than tool-specific. Unique viral sequence and host-virus junctions support exogenous retrovirus calls; HERV and LINE1 calls need locus-specific flanks and breakpoint hallmarks; shared retroviral motifs alone should be treated as ambiguous. This item is a validation layer, not a discovery caller.

## VirCapSeq off-target reuse

### Framing

- **Name:** VirCapSeq off-target reuse
- **Category:** Method question
- **Purpose:** > Assess whether off-target host reads from virus-capture sequencing can support HERV and LINE1 analysis, and define the realistic resolution ceiling given capture bias and sparse host recovery.
- **Relevant Assays:** VirCapSeq, targeted viral capture sequencing, cfDNA capture sequencing, hybrid capture metagenomics

### Technical Evaluation

### Input Requirements

- Aligned reads or BAM/CRAM with on-target and off-target reads retained
- Panel metadata showing bait design and any human-similarity filtering
- A human reference genome plus RepeatMasker or Dfam repeat annotations
- Curated HERV and LINE1 locus or family reference sets if locus-restricted analysis is attempted
- Enough off-target depth and library complexity to support family-level counting or region-of-interest summaries

- **Output Resolution:** > Best at family-level HERV and LINE1 abundance, coarse region-of-interest coverage, and at most curated whitelist-restricted locus summaries in unusually host-rich datasets; not reliable for genome-wide de novo insertion discovery or single-locus breakpoint resolution.
### Strengths

- Off-target reads can sometimes be recycled for repeat-family quantification, so the data are not automatically wasted.
- Family-level aggregation tolerates multi-mapping better than forcing ambiguous reads into a single HERV or LINE1 locus.
- Capture bycatch can still support broad background estimates, which is useful for filtering retroviral-like reads against endogenous repeat context.
- If host yield is unexpectedly high, the data may support exploratory region-of-interest or whitelist-restricted analyses.

### Limitations

- VirCapSeq panels are designed to enrich viral nucleic acid, and many probe sets intentionally avoid human-similar sequence, so host bycatch is sparse and panel-dependent.
- Off-target host coverage is nonuniform and usually too shallow for unbiased HERV or LINE1 insertion calling.
- Repetitive loci are exactly where short-read off-target data are most ambiguous, so locus-level claims are fragile without unique flanks.
- Resolution varies strongly with capture chemistry, sample input, and deduplication, so results are not directly comparable across panels.

### Recommendation

- **Recommended Role:** > Exploratory secondary analysis only: use for family-level repeat burden or background-aware filtering, not as a primary HERV or LINE1 discovery dataset.
### Validation Needs

- Measure off-target host depth and breadth first and stop short of locus-level claims if coverage is sparse.
- Benchmark any family-level signal against matched shotgun WGS or cfWGS from the same samples if available.
- Report multi-mapping reads only at family or class level unless unique flanking support exists.
- Confirm any locus-level candidate with orthogonal targeted sequencing or PCR/ddPCR.

### Key Sources

- https://genome.cshlp.org/content/early/2015/11/06/gr.191049.115.full.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4541070/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8959187/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4396974/
- https://www.nature.com/articles/s41467-021-24041-8
- https://academic.oup.com/bioinformatics/article/35/20/3913/5416145

- **Notes:** > Off-target host reads are a secondary byproduct of a virus-enrichment design, not a substitute for shotgun WGS. The realistic ceiling is coarse repeat-family burden, perhaps a small curated locus whitelist in unusually host-rich libraries, and almost never unbiased de novo insertion discovery. If the question is truly locus-level HERV or LINE1 resolution, cfWGS is the defensible primary assay.
