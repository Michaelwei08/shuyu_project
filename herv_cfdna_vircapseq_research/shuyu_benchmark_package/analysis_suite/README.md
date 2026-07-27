# Viral sequencing analysis suite (a1 .. a12)

Twelve standalone modules that turn the existing per-sample count tables, idxstats
and BAMs of the shuyu runs into the tables and figures behind the panel report.
Nothing is re-aligned here: a1 .. a4 and a7 read text tables, a5 .. a7 and a10 .. a12
may stream BAMs through `samtools view`, a8 only reshapes the `.tsv` files the others
wrote, and a9 reshapes a7's burden table against a clinical CD4 column.

**a1 .. a8 are the default run. a9 .. a12 are optional follow-ups** and are
skipped by `run_all.sh` unless their extra input exists (a9 needs a CD4 table;
a10, a11 and a12 each need BAMs).

Python standard library only, except a8, which imports matplotlib (Agg backend).
No network access. Every module writes tab-separated, pure-ASCII output into
`--outdir`, and a missing input is reported as
`WARN: <what> missing at <path>, skipping` with exit code 0 rather than a crash.

## Quickstart

All input paths in this suite are **placeholders** (`/path/to/runs`, ...): the real
controlled-data locations are deliberately not committed. Point the suite at the
real run root with `RUNS_ROOT=`, or pass `--runs-root` / `--run` explicitly per
module.

```sh
cd <this directory>
RUNS_ROOT=<real run root> OUTDIR=<scratch>/suite_out bash run_all.sh
```

`run_all.sh` runs a1 .. a7 with their documented default input paths into that one
`--outdir`, then runs a8 over those outputs into `<outdir>/figures`. Each step is
banner-logged and tee'd to `<outdir>/logs/<step>.log`; one failing step does not
stop the others, and the failures are listed in the closing summary (the script
exits 1 if any step failed). Override with `RUNS_ROOT=`, `PYTHON=`, `SAMTOOLS=`,
or run a subset with `ONLY="a1 a8"`.

a9 .. a12 are gated, not run by default:

```sh
# a9 runs only when a CD4 table is pointed at
CD4_TABLE=/path/to/clinical/cd4.tsv bash run_all.sh
# a10, a11 and a12 run only when <run>/bam/*.bam exists. A10_RUN defaults to the
# WGS panel run; A11_RUN and A12_RUN default to A10_RUN, and all three use the
# same gate.
A10_RUN=/path/to/runs/<a run with bam/> bash run_all.sh
```

A gated step that cannot run is reported as `SKIPPED`, never as a failure, so the
script's exit status is unchanged when the extra inputs are absent.

Order matters for the gated steps: a11 and a12 join a10's tables out of the same
`--outdir` (a11 takes a10's per-pair verdict, a12 takes a10's pooled hotspots) and
a12 also reads a7's burden table, so `run_all.sh` runs them after a7 and a10.
Each of those joins is optional - a missing upstream table is a `WARN`, and the
module still writes its own tables and exits 0.

## Modules

| Module | What it answers | Outputs (into `--outdir`) |
|---|---|---|
| `a1_detection_performance.py` | How well do target reads separate expected-positive from expected-negative samples? Confusion matrix at a read threshold (default >= 100), rank ROC AUC, and a threshold sweep. Without `--labels` it uses a **cohort-level proxy**, not per-sample truth, and says so in every row. | `detection_confusion.tsv`, `detection_threshold_sweep.tsv`, `detection_label_template.tsv`, `detection_sample_key.tsv` |
| `a2_reference_comparison.py` | What does adding hg38 as a competitor remove? Viral-only vs hg38-inclusive run per cohort, plus all-alignments vs primary-only inflation: reads, samples called, and pos -> neg flips per category. | `a2_reference_comparison_by_category.tsv`, `a2_reference_comparison_by_sample.tsv`, `a2_reference_comparison_flips.tsv`, `a2_sample_key.tsv` |
| `a3_kmer_ladder.py` | What does the k-mer size cost and buy in the HERV mask? Parses the already-computed ladder tables (schema is sniffed, not assumed) into long form and pivots mask-cost against residual HERV cross-mapping per k. | `kmer_ladder_long.tsv`, `kmer_ladder_summary_report.txt`, `kmer_ladder_sample_key.tsv` (only if a sample-level field was found) |
| `a4_depth_sensitivity.py` | Which detections survive subsampling to 5M reads? Per category: read ratio, detection retention vs full depth, per-sample yield vs the expected depth fraction, and a SATURATED / PARTIAL / DEPTH_LIMITED verdict. | `depth_sensitivity_by_category.tsv`, `depth_sensitivity_by_sample.tsv`, `depth_sensitivity_sample_key.tsv` |
| `a5_reference_depth_profiles.py` | Is a per-reference signal a real genome or one pile-up? Unique-best (AS > XS) coverage on named panel references: breadth, mean/median depth, binned profile, depth CV, max-bin fraction, plus a ciHHV-6 / pileup_like / active_like call. | `refprofile_summary.tsv`, `refprofile_bins.tsv`, `refprofile_sample_key.tsv` |
| `a6_htlv_junctions.py` | Is there read-level evidence that an HTLV-1 call is a real infection? Candidate host-virus junctions from discordant human-mate pairs and soft clips, clustered per sample, with a per-sample clonality proxy. **Candidates only - IGV review required.** | `htlv_junction_candidates.tsv`, `htlv_junction_per_sample.tsv`, `htlv_junction_sample_key.tsv` |
| `a7_virome_structure.py` | What is the anellovirus burden (immunocompetence proxy) and the coinfection structure? Per-sample TTV/TTMV richness, Shannon and RPM with a standard-library Mann-Whitney HIV vs HL test, plus a virus-group presence matrix and pairwise Jaccard. | `a7_virome_anellovirus_burden.tsv`, `a7_virome_anellovirus_group_test.tsv`, `a7_virome_coinfection_matrix.tsv`, `a7_virome_coinfection_pairs.tsv`, `a7_virome_virus_group_refs.tsv`, `a7_virome_sample_key.tsv` |
| `a8_figures.py` | Renders one figure per upstream table, with the finding as the title. Reads only `.tsv`, never a BAM. A missing input is a WARN and a `skipped_missing_input` row in the index. Bar charts are never drawn on a log axis - log/symlog panels use dots and lollipops, and every panel carries direct value labels. | `figures/fig_detection_threshold_sweep`, `fig_reference_comparison_by_category`, `fig_depth_sensitivity_by_category`, `fig_refprofile_coverage_tracks`, `fig_anellovirus_burden`, `fig_coinfection_pairs` (each `.png` + `.svg`), `figures/a8_figure_index.tsv`, `figures/a8_sample_key.tsv` (only if an input carried a real name) |
| `a9_cd4_correlation.py` | **Optional (needs `--cd4`).** Does the anellovirus burden rise as CD4 falls? Re-expresses a7's burden table against a user-supplied CD4 count: tie-corrected Spearman rho and Kendall tau-b per cohort x metric, clinical CD4 bands (<200 / 200-499 / >=500) with a tie-corrected Kruskal-Wallis, and a Mann-Whitney on CD4 by a7's richness>=3 cut. Standard library only, no figures. Without `--cd4` it writes only the input template and exits 0. **The pooled HIV+HL row largely restates a7's group difference - the HIV-only row is the informative one.** | `cd4_anello_correlation.tsv`, `cd4_anello_strata.tsv`, `cd4_anello_richness_contrast.tsv`, `cd4_anello_joined.tsv`, `cd4_input_template.tsv` (only when no usable CD4 table), `a9_cd4_sample_key.tsv` (only when the CD4 table used a real name) |
| `a10_anello_read_audit.py` | **Optional (needs BAMs).** Is the low-count anellovirus signal real virus or cross-mapping? Streams each BAM over its anellovirus references and scores per (sample, reference): pile-up (`max_window_fraction` over a sliding `--window`), duplicate POS+CIGAR fraction, breadth, MAPQ, soft clip, low complexity, then a `too_few_reads` / `pileup_like` / `duplicate_like` / `real_like` / `indeterminate` verdict. Adds a pooled relative-position distribution per group, a shared-hotspot table per reference, and per-sample Mann-Whitney + Fisher HIV vs HL. Chimpanzee-isolate references are a built-in negative control. **Below `--min-reads` nothing is called either way.** | `anello_read_audit_by_pair.tsv`, `anello_read_audit_by_group.tsv`, `anello_read_audit_pooled_positions.tsv`, `anello_read_audit_sample_key.tsv` (only when a BAM was found) |
| `a11_anello_read_forensics.py` | **Optional (needs BAMs).** a10 showed *that* the anellovirus reads behave like artefact; a11 asks *what they are*, from the BAM alone. Streams every read on an anellovirus reference and tests five mechanisms: soft clips matched against an embedded Illumina adapter/primer table (TruSeq R1/R2 and their stem, Nextera/Tn5 mosaic end, P5/P7, reverse complements, polyG/A/C/T), 3-mer entropy + homopolymer + DUST complexity of the aligned segment vs the clipped segment separately, mate really in the human genome, `AS == XS` equally-good alternatives (and `XA` hits), and alignment inside the `--utr-lo`/`--utr-hi` window. Collapses clips to a boundary-anchored key and reports the commonest, with how many distinct samples *and* references each recurs in. Chimpanzee-isolate references are the paired negative control (Wilcoxon), with Mann-Whitney/Fisher HIV vs HL beside it. **A recurrent clip key is a library/probe signature, not patient virus - BLAST it before naming it.** | `a11_forensics_by_pair.tsv`, `a11_forensics_by_group.tsv`, `a11_clip_sequences.tsv`, `a11_forensics_by_read.tsv` (only with `--emit-reads`), `a11_forensics_sample_key.tsv` (only when a BAM was found) |
| `a12_anello_utr_exclusion.py` | **Optional (needs BAMs).** Does a7's HIV vs HL anellovirus difference survive when every read that could come from the conserved terminal UTR is removed? Re-counts the same BAMs under a **cumulative** six-rung filter - rung 0 all reads, 1 drop anything overlapping the relative `--utr-lo`..`--utr-hi` window, 2 one read per distinct POS+CIGAR, 3 k-mer Shannon entropy >= `--min-entropy`, 4 `AS > XS` (a tie is dropped), 5 require >= `--min-distinct-positions` non-overlapping segments per (sample, reference) - where **each rung filters only what the previous rung left, so the reported counts are monotonically non-increasing**. Re-runs a7's headline Fisher (richness >= `--richness-cut`) and Mann-Whitney (burden) at every rung, sweeps the UTR window (`--sweep-lo`, plus disabled), and carries the chimpanzee arm as the residual false-positive rate. **Surviving the filter is necessary, not sufficient: nothing is realigned, and rung 1 also removes a true infection represented only through its UTR.** | `a12_utr_exclusion_ladder.tsv`, `a12_utr_exclusion_by_sample.tsv`, `a12_utr_exclusion_group_test.tsv`, `a12_utr_window_sweep.tsv`, `a12_utr_exclusion_sample_key.tsv` (only when a BAM was found) |

## Anonymisation rule

**Every sample is written as `S01 .. Snn`. The only file that may contain a real
sample identifier is `*_sample_key.tsv`.**

- Ids are assigned by **sorted real sample name** across everything one invocation
  processes, so the mapping is stable within a run. Width grows past 99 samples
  (`S001 ..`). a8 passes through labels that are already `S<digits>` and numbers
  any new real names after them. Two exceptions: **a9 does not mint ids at all** -
  it reuses a7's, joining through `a7_virome_sample_key.tsv` - and **a10, a11 and
  a12 mint their own three-digit `S001 .. Snnn`**, which therefore do NOT equal
  a7's two-digit ids; a10 and a12 carry an `a7_sample_anon` column and a11 carries
  an `a10_sample_anon` column so the tables can still be joined. a10, a11 and a12
  use the same sorted-name rule, so their ids agree whenever they process the same
  BAM set - which is what `run_all.sh` gives them.
- Each key file starts with the header comment
  `# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL`
  followed by the generating module and date. **Do not commit or email these files.**
  A key file is written only when a real name was actually handled, so an empty
  one never appears in the identifier list at the end of a run.
- Cohort group labels are derived from the real sample name and then the name is
  discarded. The rule is identical in all ten modules. The **primary** test is
  case-sensitive, `(?:^|_)(HIV|HL)[0-9]` (so `HIV<ID>` / `_HL<ID>` match but the
  lower-case cohort prefix `wgs_60samples_hiv_hl_` does not, which is what keeps
  the HL group from being swallowed by HIV). Only if that misses does the
  case-insensitive fallback apply: `_hiv` -> `HIV`, `_hl` -> `HL`, `tcl` or
  `targeted_htlv` -> `TCL`, otherwise `NA`. a2, a5, a6, a7, a10, a11 and a12 fall
  back to the run directory name for `TCL` when the sample name is uninformative.
- Warnings printed to stdout (and therefore to `<outdir>/logs/`) are anonymised
  too: BAM paths, idxstats paths and samtools error text are masked to the anon id
  or to a literal `<sample>` placeholder before printing, so a run log can be
  shared. In a11 and a12, samtools' own stderr additionally has every
  path-shaped token replaced with a literal `<path>` before it is printed, because
  samtools echoes the BAM path and masking the sample name alone would still leak
  the controlled-data mount. Skim the log anyway before sending it outside the
  group.
- No output is derived from PHI, and no controlled-access data belongs in this
  directory. See the workspace `CLAUDE.md`.

## Notes and caveats

- **a1 without `--labels` is not ground truth.** The cohort-as-label fallback
  measures agreement with a cohort assignment, not with a clinical assay. Fill in
  `detection_label_template.tsv` and rerun with `--labels` for real metrics.
- **a2 compares two different experiments**, not two measurements of one sample.
  Read the delta as the human-competitor effect, not as a change in viral load.
  Categories present in only one arm's reference are flagged `structural=yes`; a8
  plots them but never quotes them as the headline finding.
- **a3 reports the ladder, it does not choose k.** Exact shared-kmer masking is not
  the VirCAPP production rule, so the numbers are not vendor-comparable.
- **a6 emits candidates.** Targeted capture makes chimeric junction artefacts
  genuinely likely; confirm every cluster in IGV before calling integration.
- **a9 is a re-expression of a7, not new evidence.** The outcome is mostly zero,
  so every rho sits on a large tied block; CD4 and HIV status are collinear here,
  so the pooled HIV+HL rows largely restate a7's group difference; and the
  richness>=3 contrast uses the same anellovirus numbers as the correlation. Read
  the HIV-only rows, with `n_zero_metric` beside them. Cross-sectional throughout:
  no direction of causation.
- **a10 refuses to call the low-count pairs.** Below `--min-reads` (default 5) no
  method separates real virus from cross-mapping, so those pairs stay
  `too_few_reads` and are unresolved, *not* negative. The pooled 20-bin position
  distribution is descriptive only - reads within one sample are not independent -
  so only the per-sample Mann-Whitney and Fisher rows are tests. A shared hotspot
  argues artefact only when the agreeing pairs are also concentrated; read
  `median_max_window_fraction_at_hotspot` with it.
- **a11 names nothing and realigns nothing.** A clip is called an adapter only by
  its match to the embedded table; BLAST it before calling it a capture probe. A
  read that fails all five tests is *not* thereby real virus - it only means none
  of the five tested mechanisms fired. `--utr-lo`/`--utr-hi` is a coordinate
  window reported as such, not an annotated UTR. Pooled per-group fractions are
  weighted by whichever sample carries the most reads and are descriptive only;
  no p value is corrected for multiple testing.
- **a12 tests one mechanism, and its ladder is one-directional.** Surviving the
  filter is necessary but not sufficient to call real virus - a survivor can
  still be cross-mapped from an anellovirus the panel does not carry, another
  small circular DNA virus, or an unmasked human repeat. Rung 1 is conservative
  *against* detection by construction: a true infection represented in the panel
  only through its conserved UTR is removed too, which is what the window sweep
  and the chimpanzee arm are there to quantify. The rungs are cumulative, so read
  `reads_dropped_vs_previous_rung` as the marginal cost of that rung given the
  ones above it, not as an independent count.
- **The small-n tests in a9 .. a12 are normal approximations.** The
  standard-library Mann-Whitney is tie-corrected with a continuity correction,
  not exact; at n = 3 vs 3 with complete separation it reports p ~ 0.047 where
  the exact two-sided permutation p is 0.10. Fisher is exact. Read the
  Mann-Whitney p as approximate whenever a group has fewer than ~8 samples.
- Runs may legitimately produce header-only tables (no samtools, no BAM index,
  missing run directory). That is a WARN plus exit 0, by design, so the
  orchestrator can keep going.
