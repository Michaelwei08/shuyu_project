# Claim Validation

Date: 2026-04-15

This file validates the research conclusions in [report.md](</D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/report.md>) against the audited source set in [source_validation.md](</D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/validation/source_validation.md>).

## Method

- I ran the automated validator in [validate_sources.py](</D:/Stanford/research/Ash/herv_cfdna_vircapseq_research/scripts/validate_sources.py>) across all cited URLs in the project, excluding the validator's own output directory.
- The validator fetched live metadata where possible and supplemented blocked PMC pages with PubMed metadata and GitHub repositories with GitHub API metadata.
- I then reviewed the main project claims and classified each as:
  - `Direct`: explicitly supported by a primary source.
  - `Supported + inference`: supported by primary sources, but the project recommendation extends beyond what any one source directly claims.
  - `Search-based absence claim`: not contradicted in the validation pass, but inherently weaker because it depends on what was not found.

## Source Corrections Made

- Replaced the incorrect `PMC7801709` HIV citation with the correct Scientific Reports article page:
  - https://www.nature.com/articles/s41598-019-48681-5
- Replaced the incorrect `PMC9191799` HTLV citation with the correct DOI:
  - https://doi.org/10.1016/j.celrep.2019.09.016
- Replaced the incorrect Griffin link with:
  - https://www.nature.com/articles/s41467-022-35076-w
- Removed the broken LIQUORICE paper link from the source list because I did not validate a stable primary paper URL in this pass.

## Claim Matrix

| ID | Claim | Verdict | Basis | Primary sources |
| --- | --- | --- | --- | --- |
| C1 | Existing HERV tools split into DNA insertion callers versus RNA expression quantifiers. | Direct | ERVcaller, RetroSnake, xTea, and STEAK are DNA insertion or TE callers; Telescope and ERVmap are RNA-focused locus-expression resources. | https://academic.oup.com/bioinformatics/article/35/20/3913/5416145 ; https://www.sciencedirect.com/science/article/pii/S2589004222015619 ; https://www.nature.com/articles/s41467-021-24041-8 ; https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453 ; https://pubmed.ncbi.nlm.nih.gov/30455304/ |
| C2 | I did not identify a validated cfDNA-specific HERV enumeration pipeline. | Search-based absence claim | The validated source set contains HERV DNA callers, RNA quantifiers, repeat annotations, and cfDNA fragmentomics papers, but no source directly describes a mature cfDNA-native HERV insertion pipeline. This remains a literature-survey conclusion rather than a direct proof. | Entire audited source set; especially HERV caller papers, cfDNA reviews, and repeat-landscape papers |
| C3 | ViroCap or VirCapSeq enriches viral sequence and deliberately minimizes human-like probe content. | Direct | The ViroCap paper directly supports both strong viral enrichment and human-similarity filtering during probe design. | https://genome.cshlp.org/content/early/2015/11/06/gr.191049.115.full.pdf |
| C4 | Distinguishing exogenous retroviruses from HERV or LINE1 requires explicit ambiguity handling and competitive alignment. | Supported + inference | Primary sources support the underlying homology problem and the existence of capture-based proviral characterization, but the three-bin decision rule and combined decoy-reference strategy are project-specific analysis design choices. | https://pmc.ncbi.nlm.nih.gov/articles/PMC4541070/ ; https://www.nature.com/articles/s41598-019-48681-5 ; https://doi.org/10.1016/j.celrep.2019.09.016 |
| C5 | xTea is the strongest combined HERV plus LINE1 DNA caller for this project. | Supported + inference | xTea is directly validated as a broad TE caller that includes HERV and LINE1. The stronger project claim is the ranking decision that it is the best single framework when both endogenous classes matter. | https://www.nature.com/articles/s41467-021-24041-8 ; https://github.com/parklab/xTea |
| C6 | ERVcaller and RetroSnake are appropriate HERV-focused companion tools. | Supported + inference | Both tools are directly validated as HERV-aware DNA pipelines. Their recommended role as companion callers for a cfWGS pilot is a project decision, not a direct claim from the papers. | https://academic.oup.com/bioinformatics/article/35/20/3913/5416145 ; https://github.com/xunchen85/ERVcaller ; https://www.sciencedirect.com/science/article/pii/S2589004222015619 ; https://github.com/KHP-Informatics/RetroSnake |
| C7 | Telescope and ERVmap should not be presented as cfDNA DNA callers. | Direct | Telescope is explicitly an RNA-seq tool for locus-specific TE expression. ERVmap is explicitly a proviral ERV transcription resource. Neither source describes DNA insertion calling from cfDNA. | https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453 ; https://pubmed.ncbi.nlm.nih.gov/30455304/ |
| C8 | ERVmap and Telescope annotations can be repurposed as cfDNA locus whitelists. | Supported + inference | The annotation assets are real and useful, but their use as cfDNA locus whitelists is an adaptation proposed here rather than a validated standard workflow. | https://pubmed.ncbi.nlm.nih.gov/30455304/ ; https://github.com/mlbendall/telescope_annotation_db |
| C9 | VirCapSeq off-target host reads are suitable only for coarse endogenous summaries, not unbiased HERV insertion discovery. | Supported + inference | ViroCap directly supports that the assay is virus-enriching and human-depleted by design. Off-target CNV papers support the idea that bycatch can carry coarse host-genome signal. The specific conclusion about HERV insertion-discovery unsuitability is a conservative project inference. | https://genome.cshlp.org/content/early/2015/11/06/gr.191049.115.full.pdf ; https://pmc.ncbi.nlm.nih.gov/articles/PMC8959187/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4396974/ |
| C10 | Exogenous-virus calls should require unique viral sequence or host-virus junctions; endogenous HERV or LINE1 calls should require flanking or breakpoint hallmarks. | Supported + inference | The required evidence classes are consistent with capture-sequencing and TE-caller literature, but the exact validation rubric is a best-practice synthesis rather than a single-source prescription. | https://www.nature.com/articles/s41598-019-48681-5 ; https://doi.org/10.1016/j.celrep.2019.09.016 ; https://academic.oup.com/bioinformatics/article/35/20/3913/5416145 ; https://www.nature.com/articles/s41467-021-24041-8 |
| C11 | Griffin is a valid cfDNA fragmentomics framework for region-of-interest style analysis; LIQUORICE was not fully validated in this pass. | Direct for Griffin; needs caution for LIQUORICE | The Griffin citation is now corrected and validated. The project mention of LIQUORICE remains plausible, but I did not validate a stable primary source URL for it here, so it should be treated as lower-confidence until a primary citation is added. | https://www.nature.com/articles/s41467-022-35076-w |

## Bottom Line

- The main project architecture survives validation.
- The strongest direct support is for the existence and scope of the individual tools and resources.
- The weakest parts are not factual errors about tools; they are project-level inferences about how to combine those tools for cfDNA and VirCapSeq.
- The biggest concrete issues I found were citation hygiene problems, not major conceptual failures:
  - two incorrect capture-sequencing source links
  - one incorrect Griffin link
  - one unvalidated LIQUORICE source link

## Recommended Interpretation Standard

- Treat tool existence, scope, and assay fit statements as `validated`.
- Treat workflow ranking, assay preference, and combined-reference strategy statements as `validated but inferential`.
- Treat the absence of a dedicated cfDNA HERV pipeline as `current best search result`, not a proof of nonexistence.
