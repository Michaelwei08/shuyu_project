# Viral sequencing analysis suite (a1 .. a8)

Eight standalone modules that turn the existing per-sample count tables, idxstats
and BAMs of the shuyu runs into the tables and figures behind the panel report.
Nothing is re-aligned here: a1 .. a4 and a7 read text tables, a5 .. a7 may stream
BAMs through `samtools view`, and a8 only reshapes the `.tsv` files the others wrote.

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

## Anonymisation rule

**Every sample is written as `S01 .. Snn`. The only file that may contain a real
sample identifier is `*_sample_key.tsv`.**

- Ids are assigned by **sorted real sample name** across everything one invocation
  processes, so the mapping is stable within a run. Width grows past 99 samples
  (`S001 ..`). a8 passes through labels that are already `S<digits>` and numbers
  any new real names after them.
- Each key file starts with the header comment
  `# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL`
  followed by the generating module and date. **Do not commit or email these files.**
- Cohort group labels are derived from the real sample name and then the name is
  discarded. The rule is identical in all eight modules and is matched
  case-insensitively: `_HIV` -> `HIV`, `_HL` -> `HL`, `TCL` or `targeted_htlv` ->
  `TCL`, otherwise `NA`. a2, a5, a6 and a7 fall back to the run directory name for
  `TCL` when the sample name is uninformative.
- Warnings printed to stdout (and therefore to `<outdir>/logs/`) are anonymised
  too: BAM paths, idxstats paths and samtools error text are masked to the anon id
  or to a literal `<sample>` placeholder before printing, so a run log can be
  shared. Skim it anyway before sending it outside the group.
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
- Runs may legitimately produce header-only tables (no samtools, no BAM index,
  missing run directory). That is a WARN plus exit 0, by design, so the
  orchestrator can keep going.
