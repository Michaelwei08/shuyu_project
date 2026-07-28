#!/usr/bin/env python3
"""a11_anello_read_forensics.py -- what ARE the suspect anellovirus reads?

THE QUESTION
  a7 reported an anellovirus burden higher in HIV+ than in HL in the WGS cohort
  (37 vs 23 samples; richness >= 3 in 18/37 HIV+ and 0/23 HL, Fisher
  p = 2.3e-05). a10 then audited those reads and found the signal dominated by
  artefact: median max_window_fraction 0.729, median
  duplicate_position_fraction 0.552, median fraction soft-clipped >= 10 bp
  0.960, median breadth 0.068; the CHIMPANZEE anellovirus references - a
  negative control, since no human carries chimpanzee TTV - showed the same
  profile on 1952 reads in the HIV group; a hotspot at relative position
  0.81-0.93 recurred across unrelated references and independent samples
  (NC_014075.1 17/19 pairs at rel 0.813, NC_002076.2 17/24 at 0.846,
  NC_014078.1 13/20 at 0.904); and requiring a read-validated pair dropped the
  group difference to 8/37 vs 1/23, p = 0.134.

  a10 established THAT those reads behave like artefact. It did not establish
  WHAT they are. This module asks that from the BAM alone. In a ~3.7 kb
  anellovirus genome the 0.75-1.00 region is the conserved UTR, i.e. the
  expected cross-mapping magnet, so five mechanisms are tested directly:
  library adapter / primer sequence hanging off the alignment, low-complexity
  sequence, a mate that really sits in the human genome, an alignment that is
  equally good on several anellovirus references, and an alignment that sits in
  the conserved UTR.

WHAT IT COMPUTES
  Every read aligned to an anellovirus reference is streamed with
  samtools view -F <--exclude-flags> (default 0x904: unmapped, secondary,
  supplementary) and -q <--min-mapq>. MAPQ is NOT filtered by default and must
  not be: multi-mapping reads ARE the artefact under investigation, and a MAPQ
  filter would delete the evidence before it is measured. Duplicates are kept
  for the same reason. Per read:

  (1) SOFT CLIP ORIGIN. The CIGAR splits the read into an aligned segment, a
      left clip and a right clip (leading / trailing S, H skipped; SEQ excludes
      hard-clipped bases). Every clip >= --min-clip (default 10 bp) is tested
      against an embedded table of standard Illumina adapter / primer
      sequences: TruSeq universal (R1) and indexed (R2) adapters and their
      shared 13 bp stem, the Nextera / Tn5 mosaic end and its read-through
      complement, the P5 and P7 flowcell primers, the reverse complement of
      each of those, and polyG / polyA (plus their complements polyC / polyT)
      runs. A right clip is compared 5'->3' against the adapter prefix, a left
      clip is compared against the adapter suffix, over the overlapping
      length, allowing --adapter-mismatch mismatches (default 1); the longest
      overlap with the fewest mismatches wins and its name is reported.
      Clip sequences are then collapsed - anchored at the alignment boundary
      and truncated to --clip-key-len bases (default 20, 0 = whole clip) so
      that clips of different lengths from one adapter or probe collapse to one
      key - and the --top-clips (default 30) commonest keys are written out
      with their counts, the number of DISTINCT SAMPLES and DISTINCT REFERENCES
      they occur in, their own complexity, and any adapter match. ONE CLIP
      SEQUENCE RECURRING ACROSS MANY UNRELATED SAMPLES AND REFERENCES IS A
      PROBE / ADAPTER / CONTAMINANT SIGNATURE, and it is the key output of this
      module. A key seen in >= --recurrent-min-samples samples AND
      >= --recurrent-min-refs references is counted as recurrent.

  (2) COMPLEXITY, computed separately for the ALIGNED segment and for the
      CLIPPED segments of the same read, because the discriminating question is
      whether only the clipped part is junk (adapter read-through on an
      otherwise real alignment) or the whole read is (a low-complexity fragment
      that lands anywhere). Per segment: Shannon entropy over overlapping
      3-mers in bits, the longest homopolymer run as a fraction of segment
      length, GC fraction over ACGT, and a DUST-like triplet score
      sum(c*(c-1)/2) / (t-1) over t = L-2 triplets, normalised by its
      homopolymer maximum t/2 so that 0 = complex and 1.0 = every triplet
      identical regardless of length. A segment is flagged low complexity when
      entropy3 < --min-entropy (default 1.2 bits) or the homopolymer fraction
      is > --max-homopolymer-frac (default 0.5). Note that entropy over 3-mers
      is bounded by log2(L-2), so a 12 bp clip cannot exceed 3.3 bits; the
      1.2-bit flag is therefore a homopolymer / dinucleotide-run detector on
      short clips, not a general complexity test.

  (3) MATE ORIGIN, from FLAG, RNEXT and PNEXT against the reference map:
      same_anello_ref, other_anello_ref, chimp_anello_ref, HUMAN,
      other_category, unmapped, no_mate. A read whose mate sits reliably in the
      human genome is host-derived or a chimera, not virus. References with no
      reference-map entry fall into other_category and are also counted in
      n_mate_unknown_ref, so that bucket is never read as evidence.

  (4) ALIGNMENT AMBIGUITY. AS:i and XS:i are parsed and AS == XS is flagged
      (an equally good alignment exists elsewhere); when an XA:Z tag is present
      its alternate hits are counted and split into anellovirus and chimpanzee
      anellovirus hits. This is the cheap stand-in for a unique-k-mer test: a
      read with several equal-scoring anellovirus hits cannot support a
      species-level call, whatever its MAPQ.

  (5) UTR OVERLAP. The midpoint of the aligned segment is expressed as a
      relative position along its own reference and flagged when it falls in
      --utr-lo .. --utr-hi (default 0.75 .. 1.00).

  Everything is aggregated per (sample, reference) and per group (HIV / HL /
  TCL / ALL), for two reference sets: human anellovirus references and the
  chimpanzee-isolate references (--chimp-accessions), which are the built-in
  NEGATIVE CONTROL. Every metric is then contrasted human vs chimp with a
  paired Wilcoxon signed-rank test over the samples that carry both, and the
  two cohorts are compared HIV vs HL with Mann-Whitney U and Fisher exact.
  Every p value is printed with the n of each side beside it. IF THE HUMAN
  REFERENCES MATCH THE CHIMPANZEE CONTROL ON THESE METRICS, THEY ARE BEHAVING
  LIKE THE CONTROL, i.e. like cross-mapping.

WHAT THIS CAN AND CANNOT SUPPORT
  CAN: identify clip sequence that is library adapter, primer or polyG, and
  identify a clip sequence recurring across unrelated samples and references -
  which no per-sample biological model explains.
  CAN: show that the aligned part of the read is itself low complexity, that
  the mates sit in the human genome, or that the alignments are ambiguous
  between anellovirus references. Any one of those is sufficient to refuse a
  species-level anellovirus call.
  CAN: show that the human references and the chimpanzee negative control are
  indistinguishable, which is the strongest statement available here.
  CANNOT: prove the true origin of a read. Nothing is re-aligned and no
  external database is consulted (no network), so "not adapter, not low
  complexity, mate not human" is NOT evidence of real virus - it only means
  none of the five tested mechanisms fired.
  CANNOT: a single short adapter hit is weak. A 10 bp overlap with one
  mismatch matches a random clip about 3e-5 of the time per adapter entry, so
  roughly 1 in 1000 reads will match something by chance; recurrence across
  samples, not a single hit, is the evidence.
  CANNOT: name the UTR. --utr-lo/--utr-hi is a coordinate window in the
  reference's own orientation, reported and never assumed; check it against the
  reference annotation before calling it the UTR.
  CANNOT: treat reads as independent. Pooled per-group fractions are
  descriptive, weighted by whichever sample carries the most reads. Only the
  per-sample Wilcoxon, Mann-Whitney and Fisher rows are tests, and their n is
  small.
  CANNOT: distinguish PCR duplicates from independent evidence; that is a10's
  duplicate_position_fraction, joined here as context rather than recomputed.

WHAT IT WRITES (tab separated, pure ASCII, into --outdir)
  a11_forensics_by_pair.tsv    one row per (sample, anellovirus reference)
  a11_forensics_by_group.tsv   per (group, reference set) summaries, the
                               human-vs-chimp contrast on every metric, and the
                               HIV vs HL tests
  a11_clip_sequences.tsv       the top --top-clips collapsed clip sequences
                               with sample / reference recurrence
  a11_forensics_by_read.tsv    per-read detail; written ONLY with --emit-reads
                               because it is large. The read name is never
                               written - QNAME can embed a sample or run
                               identifier - only a truncated SHA1 of it.
  a11_forensics_sample_key.tsv real -> anonymous mapping. THE ONLY file that
                               contains real sample identifiers, written only
                               when a real name was actually read. Do not
                               commit or email it.
  Samples are anonymised to S001..Snnn by sorted real name, the same three-digit
  scheme a10 uses; when a10's own key is found in --indir its ids and verdicts
  are joined in as context (a10_sample_anon, a10_verdict) so the two tables can
  be lined up without either of them carrying a real name.

  Standard library only. No figures, so matplotlib is not imported. No network.
  Any missing input prints "WARN: <what> missing at <path>, skipping" and the
  module still writes headed tables and exits 0.

EXAMPLE
  python3 a11_anello_read_forensics.py \
      --run /path/to/runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only \
      --refmap /path/to/runs/shuyu_masked_panel_hg38_herv_line1_refixed/ref/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv \
      --indir /path/to/runs/panel_report_20260725/suite_out \
      --outdir /path/to/runs/panel_report_20260725/suite_out \
      --min-clip 10 --top-clips 30

  # RUNS_ROOT is honoured for the default input paths, as in run_all.sh:
  RUNS_ROOT=/real/run/root python3 a11_anello_read_forensics.py --outdir <scratch>

Written 2026-07-27.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import hashlib
import math
import os
import re
import subprocess
import sys
import tempfile

# Placeholder root: the real controlled-data location is deliberately not
# committed. run_all.sh exports RUNS_ROOT, so honour it here too.
RUNS_ROOT = (os.environ.get("SHUYU_RUNS_ROOT")
             or os.environ.get("RUNS_ROOT")
             or "/path/to/runs")

DEFAULT_RUNS = [RUNS_ROOT + "/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only"]
DEFAULT_REFMAP = (RUNS_ROOT + "/shuyu_masked_panel_hg38_herv_line1_refixed/ref/"
                  "hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv")
DEFAULT_INDIR = RUNS_ROOT + "/panel_report_20260725/suite_out"
DEFAULT_OUTDIR = RUNS_ROOT + "/panel_report_20260725/suite_out"
DEFAULT_PREFIX = "a11_forensics"
DEFAULT_ANELLO_ACC_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "anello_accessions.txt")
DEFAULT_CHIMP_ACC = "NC_014069.1,NC_014077.1,NC_014480.2"

# a10 tables looked for in --indir. The suffixes deliberately carry a10's own
# prefix: a bare "_by_pair.tsv" would also match THIS module's output sitting in
# the same directory, and a bare "_sample_key.tsv" would match a7's key.
A10_PAIR_NAME = "anello_read_audit_by_pair.tsv"
A10_PAIR_SUFFIX = "read_audit_by_pair.tsv"
A10_KEY_NAME = "anello_read_audit_sample_key.tsv"
A10_KEY_SUFFIX = "read_audit_sample_key.tsv"

ANELLO_KEYWORDS = [
    "anello", "torque teno", "torque-teno", "torquetenovirus",
    "transfusion transmitted virus", "transfusion-transmitted virus",
    "small anellovirus", "tt virus", "ttv", "ttmv", "ttmdv", "sen virus",
]
# Short keywords that must not match inside a longer word.
SHORT_KEYWORDS = set(["ttv", "ttmv", "ttmdv", "tt virus"])

ACC_RE = re.compile(
    r"(?<![A-Za-z0-9])([A-Z]{2}_\d{5,8}|[A-Z]{1,2}\d{5,6})(?:\.\d+)?(?![A-Za-z0-9])")

# Fallback human-reference recognition, used ONLY when the reference map has no
# entry for a mate's reference. chrEBV is in the hg38 analysis set and is viral,
# so it is excluded explicitly.
HUMAN_NAME_RE = re.compile(r"^(chr(\d{1,2}|X|Y|M)(_[A-Za-z0-9]+)?|"
                           r"(GL|KI|KZ|ML|JH|KN|KQ|KV)\d{4,7}(\.\d+)?)$")

GROUP_ORDER = ["HIV", "HL", "TCL", "NA"]

MATE_CLASSES = ["same_anello_ref", "other_anello_ref", "chimp_anello_ref",
                "HUMAN", "other_category", "unmapped", "no_mate"]

# Published Illumina adapter / primer stems. Reverse complements are generated
# at run time, so read orientation does not matter.
BASE_ADAPTERS = [
    ("TruSeq_R1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),
    ("TruSeq_R2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
    ("TruSeq_common_stem", "AGATCGGAAGAGC"),
    ("Nextera_Tn5_mosaic_end", "AGATGTGTATAAGAGACAG"),
    ("Illumina_P5_primer", "AATGATACGGCGACCACCGAGATCTACAC"),
    ("Illumina_P7_primer", "CAAGCAGAAGACGGCATACGAGAT"),
]
# Homopolymers are their own entries (polyC / polyT are the reverse complements
# of the polyG / polyA artefacts) and are not reverse-complemented again.
HOMOPOLYMER_ADAPTERS = [
    ("polyG", "G" * 32),
    ("polyA", "A" * 32),
    ("polyC", "C" * 32),
    ("polyT", "T" * 32),
]

TODAY = datetime.date.today().isoformat()
SCRIPT = os.path.basename(__file__)

# Accumulator layout. Everything is a plain counter or a list, so one (sample,
# reference) pair, one sample and one whole group are summarised by exactly the
# same merge + derive pair of functions.
COUNT_KEYS = [
    "n_reads", "n_seq", "n_mapq0",
    "n_reads_clipped", "n_clips", "n_clips_adapter", "n_reads_adapter",
    "aln_n", "aln_lowcomp", "clip_n", "clip_lowcomp",
    "n_with_as", "n_with_xs", "n_as_eq_xs", "n_xs_gt_as",
    "n_with_xa", "n_xa_hits", "n_xa_anello", "n_xa_chimp",
    "n_utr", "n_mate_unknown", "n_detail", "n_detail_skipped",
] + ["mate_" + c for c in MATE_CLASSES]

LIST_KEYS = [
    "mapqs", "clip_lens",
    "aln_len", "aln_ent", "aln_hp", "aln_gc", "aln_dust",
    "clip_len_detail", "clip_ent", "clip_hp", "clip_gc", "clip_dust",
    "xa_hits", "xa_anello_hits",
]

PAIR_COLUMNS = [
    "run", "sample_anon", "group", "reference_id", "ref_label", "ref_len",
    "chimp_flagged", "idxstats_mapped", "n_reads", "n_reads_with_seq",
    "median_mapq", "frac_mapq0",
    "n_reads_clipped", "frac_reads_clipped", "n_clips", "median_clip_len",
    "max_clip_len", "n_clips_adapter", "frac_clips_adapter",
    "n_reads_adapter", "frac_reads_adapter", "top_adapter",
    "top_adapter_clips", "n_distinct_clip_keys", "top_clip_key",
    "top_clip_key_clips", "n_recurrent_clip_keys", "n_clips_recurrent",
    "aln_median_len", "aln_median_entropy3", "aln_median_homopolymer_frac",
    "aln_median_gc", "aln_median_dust", "aln_frac_low_complexity",
    "clip_median_len", "clip_median_entropy3", "clip_median_homopolymer_frac",
    "clip_median_gc", "clip_median_dust", "clip_frac_low_complexity",
    "n_mate_same_anello_ref", "n_mate_other_anello_ref",
    "n_mate_chimp_anello_ref", "n_mate_HUMAN", "n_mate_other_category",
    "n_mate_unmapped", "n_mate_no_mate", "n_mate_unknown_ref",
    "frac_mate_same_anello_ref", "frac_mate_other_anello_ref",
    "frac_mate_HUMAN",
    "n_with_as", "n_with_xs", "n_as_eq_xs", "frac_as_eq_xs",
    "n_xs_gt_as", "frac_xs_gt_as",
    "n_with_xa", "frac_with_xa", "n_xa_hits", "n_xa_hits_anello",
    "frac_xa_hits_anello", "median_xa_hits", "median_xa_anello_hits",
    "n_utr_reads", "frac_utr_reads", "utr_lo", "utr_hi",
    "detail_capped", "a10_sample_anon", "a10_verdict",
    "a10_max_window_fraction", "flags",
]

GROUP_COLUMNS = [
    "row_type", "label", "ref_set", "n_samples", "n_pairs", "n_reads",
    "n_reads_with_seq", "frac_reads_clipped", "n_clips",
    "frac_clips_adapter", "frac_reads_adapter", "top_adapter",
    "aln_median_entropy3", "aln_frac_low_complexity",
    "clip_median_entropy3", "clip_frac_low_complexity",
    "frac_mate_same_anello_ref", "frac_mate_other_anello_ref",
    "frac_mate_chimp_anello_ref", "frac_mate_HUMAN",
    "frac_mate_other_category", "frac_mate_unmapped", "frac_mate_no_mate",
    "frac_as_eq_xs", "frac_with_xa", "frac_xa_hits_anello", "frac_utr_reads",
    "n_clip_keys", "n_recurrent_clip_keys", "top_clip_key",
    "top_clip_key_clips", "top_clip_key_samples", "top_clip_key_refs",
    "group1", "n1", "value_group1", "group2", "n2", "value_group2",
    "statistic", "statistic_value", "p_value", "effect", "note",
]

CLIP_COLUMNS = [
    "rank", "clip_key", "key_len", "n_clips", "n_reads", "n_samples",
    "n_references", "n_groups", "groups", "n_human_references",
    "n_chimp_references", "sides", "recurrent", "adapter", "adapter_overlap",
    "adapter_mismatches", "median_clip_len", "entropy3", "homopolymer_frac",
    "gc_frac", "dust_score", "low_complexity", "example_reference", "note",
]

READ_COLUMNS = [
    "run", "sample_anon", "group", "reference_id", "chimp_flagged",
    "read_key", "flag", "mapq", "pos", "aln_ref_span", "rel_midpoint",
    "in_utr", "cigar", "left_clip_len", "right_clip_len",
    "clip_adapter", "clip_adapter_side", "clip_adapter_overlap",
    "clip_adapter_mismatches", "left_clip_key", "right_clip_key",
    "aln_len", "aln_entropy3", "aln_homopolymer_frac", "aln_gc", "aln_dust",
    "aln_low_complexity", "clip_len", "clip_entropy3",
    "clip_homopolymer_frac", "clip_gc", "clip_dust", "clip_low_complexity",
    "mate_class", "mate_ref", "mate_pos", "AS", "XS", "as_eq_xs",
    "n_xa_hits", "n_xa_anello", "n_xa_chimp",
]

# (column label, key in the derived-metric dict). Used for BOTH the HIV vs HL
# test rows and the human-vs-chimpanzee contrast, so "every metric" means the
# same list in both places.
CONTRAST_METRICS = [
    ("frac_reads_clipped", "frac_reads_clipped"),
    ("frac_clips_adapter", "frac_clips_adapter"),
    ("frac_reads_adapter", "frac_reads_adapter"),
    ("aln_median_entropy3", "aln_median_entropy3"),
    ("aln_frac_low_complexity", "aln_frac_low_complexity"),
    ("clip_median_entropy3", "clip_median_entropy3"),
    ("clip_frac_low_complexity", "clip_frac_low_complexity"),
    ("clip_median_dust", "clip_median_dust"),
    ("frac_mate_same_anello_ref", "frac_mate_same_anello_ref"),
    ("frac_mate_other_anello_ref", "frac_mate_other_anello_ref"),
    ("frac_mate_HUMAN", "frac_mate_HUMAN"),
    ("frac_mate_unmapped", "frac_mate_unmapped"),
    ("frac_as_eq_xs", "frac_as_eq_xs"),
    ("frac_with_xa", "frac_with_xa"),
    ("frac_xa_hits_anello", "frac_xa_hits_anello"),
    ("frac_utr_reads", "frac_utr_reads"),
    ("frac_mapq0", "frac_mapq0"),
    ("median_mapq", "median_mapq"),
]


# --------------------------------------------------------------------------- #
# small utilities
# --------------------------------------------------------------------------- #
def warn_missing(what, path):
    print("WARN: %s missing at %s, skipping" % (what, path))


def to_ascii(text):
    if text is None:
        return ""
    out = []
    for ch in str(text):
        out.append(ch if 32 <= ord(ch) < 127 else "?")
    return "".join(out).replace("\t", " ").strip()


def fnum(value, digits=4):
    if value is None:
        return "NA"
    try:
        if isinstance(value, int):
            return str(value)
        if value != value or value in (float("inf"), float("-inf")):
            return "NA"
        return ("%." + str(digits) + "f") % value
    except (TypeError, ValueError):
        return "NA"


def fp(value):
    if value is None:
        return "NA"
    return "%.6g" % value


def median(values):
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return float(s[mid])
    return (float(s[mid - 1]) + float(s[mid])) / 2.0


def mean(values):
    if not values:
        return None
    return float(sum(values)) / float(len(values))


def safe_frac(num, den):
    if not den:
        return None
    return float(num) / float(den)


def redact(text, sample, anon):
    """Replace a real sample name inside a string with its anonymous ID."""
    if not sample:
        return text
    return str(text).replace(sample, anon)


# Any whitespace-delimited token carrying a path separator. samtools echoes the
# BAM path in most of its error messages, and that path names the
# controlled-data mount, so it must never reach stdout even after the sample
# name itself has been masked out.
PATH_TOKEN_RE = re.compile(r"[^\s'\"]*[/\\][^\s'\"]*")


def strip_paths(text):
    """Replace every filesystem-path-looking token with a literal <path>."""
    return PATH_TOKEN_RE.sub("<path>", str(text or ""))


def out_name(prefix, base):
    return (prefix + "_" + base) if prefix else base


def clip_out_name(prefix):
    """Name of the clip table.

    The suite writes <prefix>_<base>, but this table is the module's headline
    output and is specified as a11_clip_sequences.tsv, so a trailing
    "_forensics" is dropped from the prefix. Any other --prefix is used as is.
    """
    base = prefix or ""
    if base.endswith("_forensics"):
        base = base[:-len("_forensics")]
    return out_name(base, "clip_sequences.tsv")


def write_tsv(path, comments, header, rows):
    with open(path, "w", encoding="ascii", errors="replace", newline="") as fh:
        for line in comments:
            fh.write("# " + to_ascii(line) + "\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(to_ascii(c) for c in row) + "\n")


# --------------------------------------------------------------------------- #
# sample naming (identical rule to a5 / a7 / a10)
# --------------------------------------------------------------------------- #
def group_of(sample, run_name, use_run_name=True):
    """HIV / HL / TCL / NA from the real sample name (spec-defined rules)."""
    # The sample-specific label is HIV/HL immediately followed by a digit
    # (HIV<ID>, HL<ID>). Match that first and case-sensitively: the WGS cohort
    # prefix "wgs_60samples_hiv_hl_" is lowercase, and an upper()/lower()
    # test for "_HIV" would otherwise label every WGS sample HIV and leave
    # the HL group empty.
    _m = re.search(r"(?:^|_)(HIV|HL)[0-9]", sample or "")
    if _m:
        return _m.group(1)
    low = (sample or "").lower()
    if "_hiv" in low:
        return "HIV"
    if "_hl" in low:
        return "HL"
    if "targeted_htlv" in low or "tcl" in low:
        return "TCL"
    if use_run_name:
        rlow = (run_name or "").lower()
        if "targeted_htlv" in rlow or "tcl" in rlow:
            return "TCL"
    return "NA"


def anonymise(real_names):
    """real sample name -> S001..Snnn, ordered by sorted real name."""
    uniq = sorted(set(real_names))
    width = max(3, len(str(len(uniq))))
    return dict((name, "S" + str(i + 1).zfill(width))
                for i, name in enumerate(uniq))


# --------------------------------------------------------------------------- #
# statistics (standard library only)
# --------------------------------------------------------------------------- #
def mann_whitney_u(x, y):
    """Two-sided Mann-Whitney U with tie correction and normal approximation.

    Returns a dict with U for each sample, z, p, and the rank-biserial effect
    size (positive means group x tends to be larger). None if either group is
    empty.
    """
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return None
    marked = [(float(v), 0) for v in x] + [(float(v), 1) for v in y]
    marked.sort(key=lambda t: t[0])
    n = n1 + n2
    ranks = [0.0] * n
    tie_sizes = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and marked[j + 1][0] == marked[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = avg_rank
        tie_sizes.append(j - i + 1)
        i = j + 1
    r1 = 0.0
    for rank, item in zip(ranks, marked):
        if item[1] == 0:
            r1 += rank
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = float(n1) * float(n2) - u1
    mu = float(n1) * float(n2) / 2.0
    tie_term = sum(t ** 3 - t for t in tie_sizes)
    if n > 1:
        var = (float(n1) * float(n2) / 12.0) * (
            (n + 1) - tie_term / float(n * (n - 1)))
    else:
        var = 0.0
    if var <= 0:
        z, p = 0.0, 1.0
    else:
        sd = math.sqrt(var)
        if u1 > mu:
            z = (u1 - mu - 0.5) / sd
        elif u1 < mu:
            z = (u1 - mu + 0.5) / sd
        else:
            z = 0.0
        p = math.erfc(abs(z) / math.sqrt(2.0))
        p = min(1.0, max(0.0, p))
    return {
        "n1": n1, "n2": n2, "U1": u1, "U2": u2, "z": z, "p": p,
        "median1": median(x), "median2": median(y),
        "mean1": mean(x), "mean2": mean(y),
        "effect_r": 2.0 * u1 / (float(n1) * float(n2)) - 1.0,
    }


def wilcoxon_signed_rank(pairs):
    """Two-sided paired Wilcoxon signed-rank, normal approximation.

    pairs is a list of (x, y) from the SAME sample - here the value on the
    human anellovirus references and on the chimpanzee control references. A
    paired test is used because the two sides are not independent groups: one
    sample contributes to both. Zero differences are dropped (Wilcoxon's
    original treatment), ties take average ranks with the standard tie
    correction, and a 0.5 continuity correction is applied. Returns None when
    no non-zero difference is left.
    """
    usable = [(float(x), float(y)) for x, y in pairs
              if x is not None and y is not None]
    diffs = [x - y for x, y in usable]
    nonzero = [d for d in diffs if d != 0.0]
    n = len(nonzero)
    if n == 0:
        return None
    order = sorted(range(n), key=lambda i: abs(nonzero[i]))
    ranks = [0.0] * n
    tie_sizes = []
    i = 0
    while i < n:
        j = i
        while (j + 1 < n
               and abs(nonzero[order[j + 1]]) == abs(nonzero[order[i]])):
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg_rank
        tie_sizes.append(j - i + 1)
        i = j + 1
    w_plus = sum(r for r, d in zip(ranks, nonzero) if d > 0)
    w_minus = sum(r for r, d in zip(ranks, nonzero) if d < 0)
    mu = n * (n + 1) / 4.0
    var = n * (n + 1) * (2 * n + 1) / 24.0
    var -= sum(t ** 3 - t for t in tie_sizes) / 48.0
    if var <= 0:
        z, p = 0.0, 1.0
    else:
        sd = math.sqrt(var)
        if w_plus > mu:
            z = (w_plus - mu - 0.5) / sd
        elif w_plus < mu:
            z = (w_plus - mu + 0.5) / sd
        else:
            z = 0.0
        p = min(1.0, max(0.0, math.erfc(abs(z) / math.sqrt(2.0))))
    return {
        "n_pairs": len(usable), "n_nonzero": n,
        "n_zero": len(usable) - n,
        "W_plus": w_plus, "W_minus": w_minus, "z": z, "p": p,
        "median1": median([x for x, _y in usable]),
        "median2": median([y for _x, y in usable]),
        "median_diff": median(diffs),
    }


def _log_choose(n, k):
    if k < 0 or k > n:
        return None
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def fisher_exact_2x2(a, b, c, d):
    """Two-sided Fisher exact test on [[a, b], [c, d]], standard library only.

    a, b = successes / failures in group 1; c, d = the same in group 2. The
    two-sided p sums the hypergeometric probability of every table at least as
    unlikely as the observed one. Returns None when the table is degenerate.
    """
    for v in (a, b, c, d):
        if v is None or v < 0:
            return None
    r1, r2 = a + b, c + d
    c1, c2 = a + c, b + d
    n = r1 + r2
    if n <= 0 or r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
        return None
    log_den = _log_choose(n, c1)
    probs = {}
    lo = max(0, c1 - r2)
    hi = min(r1, c1)
    for k in range(lo, hi + 1):
        lp = _log_choose(r1, k) + _log_choose(r2, c1 - k) - log_den
        probs[k] = math.exp(lp)
    total = sum(probs.values())
    if total <= 0:
        return None
    p_obs = probs.get(a, 0.0)
    tol = 1.0 + 1e-9
    p_two = sum(v for v in probs.values() if v <= p_obs * tol) / total
    p_two = min(1.0, max(0.0, p_two))
    odds = None
    if b > 0 and c > 0:
        odds = (float(a) * float(d)) / (float(b) * float(c))
    elif a > 0 and d > 0:
        odds = float("inf")
    return {
        "a": a, "b": b, "c": c, "d": d, "n1": r1, "n2": r2,
        "p": p_two, "odds_ratio": odds,
        "prop1": (float(a) / float(r1)) if r1 else None,
        "prop2": (float(c) / float(r2)) if r2 else None,
    }


# --------------------------------------------------------------------------- #
# reference identification (same rules as a7 / a10)
# --------------------------------------------------------------------------- #
def accessions_in(text):
    return set(m.group(1).upper() for m in ACC_RE.finditer(text or ""))


def norm_acc(token):
    token = (token or "").strip().upper()
    if not token:
        return ""
    return token.split(".")[0]


def keyword_hit(text_low, keyword):
    if keyword in SHORT_KEYWORDS:
        return re.search(r"(?<![a-z0-9])" + re.escape(keyword) + r"(?![a-z0-9])",
                         text_low) is not None
    return keyword in text_low


def read_accession_list(path):
    """One accession or reference_id per line; '#' starts a comment."""
    items = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    items.append(line)
    except OSError:
        return None
    return items


def load_refmap(path):
    """reference_id -> dict(category, description, text, accs). {} if unusable."""
    refs = {}
    if not path:
        return refs
    if not os.path.exists(path):
        warn_missing("reference map", path)
        return refs
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh):
                rid = (row.get("reference_id") or "").strip()
                if not rid:
                    continue
                parts = [rid]
                for key in ("category", "source", "description",
                            "original_reference_id", "original_header"):
                    parts.append((row.get(key) or "").strip())
                text = to_ascii(" ".join(p for p in parts if p))
                refs[rid] = {
                    "category": (row.get("category") or "").strip().upper(),
                    "description": to_ascii(row.get("description") or ""),
                    "text": text,
                    "accs": set(norm_acc(a) for a in accessions_in(text)),
                }
    except (OSError, csv.Error) as exc:
        print("WARN: could not parse reference map at %s (%s), continuing "
              "without it" % (path, exc.__class__.__name__))
        return {}
    return refs


def ref_accs(reference_id, refmap):
    accs = set(norm_acc(a) for a in accessions_in(reference_id))
    info = refmap.get(reference_id)
    if info:
        accs |= info["accs"]
    return accs


def is_anello_ref(reference_id, refmap, anello_norm, anello_ids):
    if reference_id in anello_ids:
        return True
    if ref_accs(reference_id, refmap) & anello_norm:
        return True
    info = refmap.get(reference_id)
    if not info:
        return False
    if info["category"] in ("HUMAN", "HERV", "LINE1"):
        return False
    low = info["text"].lower()
    for kw in ANELLO_KEYWORDS:
        if keyword_hit(low, kw):
            return True
    return False


def is_chimp_ref(reference_id, refmap, chimp_norm, chimp_ids):
    if reference_id in chimp_ids:
        return True
    if ref_accs(reference_id, refmap) & chimp_norm:
        return True
    info = refmap.get(reference_id)
    if info:
        low = info["text"].lower()
        if "chimpanzee" in low or "pan troglodytes" in low:
            return True
    return False


def sanitize_label(text):
    out = re.sub(r"[^0-9A-Za-z]+", "_", str(text)).strip("_")
    return out[:48] if out else "NA"


def ref_label(reference_id, refmap):
    info = refmap.get(reference_id)
    desc = info["description"] if info else ""
    return sanitize_label(desc) if desc else sanitize_label(reference_id)


class RefClassifier(object):
    """reference name -> ANELLO / ANELLO_CHIMP / HUMAN / OTHER / UNKNOWN.

    Cached, because it is asked once per read for the mate's reference.
    UNKNOWN means the reference map had no entry and the name did not look like
    a human contig; it is reported separately so it is never read as evidence.
    """

    def __init__(self, refmap, anello_norm, anello_ids, chimp_norm, chimp_ids):
        self.refmap = refmap
        self.anello_norm = anello_norm
        self.anello_ids = anello_ids
        self.chimp_norm = chimp_norm
        self.chimp_ids = chimp_ids
        self.cache = {}

    def classify(self, name):
        hit = self.cache.get(name)
        if hit is not None:
            return hit
        out = self._classify(name)
        self.cache[name] = out
        return out

    def _classify(self, name):
        if is_anello_ref(name, self.refmap, self.anello_norm, self.anello_ids):
            if is_chimp_ref(name, self.refmap, self.chimp_norm, self.chimp_ids):
                return "ANELLO_CHIMP"
            return "ANELLO"
        info = self.refmap.get(name)
        if info:
            if info["category"] == "HUMAN":
                return "HUMAN"
            return "OTHER"
        if name.lower().startswith("chrebv"):
            return "OTHER"
        if HUMAN_NAME_RE.match(name or ""):
            return "HUMAN"
        return "UNKNOWN"


# --------------------------------------------------------------------------- #
# run inputs
# --------------------------------------------------------------------------- #
def list_bams(run_dir, bam_glob):
    return sorted(p for p in glob.glob(os.path.join(run_dir, bam_glob))
                  if os.path.isfile(p))


# Suffixes some pipelines add to the BAM basename but not to the sibling
# idxstats / count-table sample name. Stripped when deriving a sample id.
BAM_NAME_SUFFIXES = (".retrovirus", ".retro", ".markdup", ".dedup",
                     ".sorted", ".filtered")


def sample_of_bam(bam_path):
    base = os.path.basename(bam_path)
    if base.endswith(".bam"):
        base = base[:-4]
    # Pipeline BAMs are named <sample>.retrovirus.bam while the idxstats
    # beside them are <sample>.idxstats.tsv. Leaving the suffix on makes
    # every BAM-to-idxstats join fail silently (empty RPM denominators).
    for _suf in BAM_NAME_SUFFIXES:
        if base.endswith(_suf):
            base = base[: -len(_suf)]
            break
    return base


def bam_is_indexed(bam):
    return (os.path.exists(bam + ".bai")
            or os.path.exists(os.path.splitext(bam)[0] + ".bai")
            or os.path.exists(bam + ".csi"))


def idxstats_lengths(run_dir, sample):
    """{refname: (seqlen, mapped)} from results/*<sample>.idxstats.tsv, or {}."""
    res = os.path.join(run_dir, "results")
    tail = sample + ".idxstats.tsv"
    cands = [os.path.join(res, tail)]
    # tolerate a filename prefix such as "primary_only_" but never match a
    # different sample whose name merely ends with this one
    for hit in sorted(glob.glob(os.path.join(res, "*" + tail))):
        if os.path.basename(hit).endswith("_" + tail):
            cands.append(hit)
    path = next((c for c in cands if os.path.exists(c)), None)
    if path is None:
        return {}
    out = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4 or p[0] == "*":
                    continue
                try:
                    out[p[0]] = (int(p[1]), int(p[2]))
                except ValueError:
                    continue
    except OSError:
        return {}
    return out


def header_lengths(samtools, bam_path):
    """{refname: (seqlen, None)} from the BAM header (fallback for idxstats)."""
    try:
        proc = subprocess.run([samtools, "view", "-H", bam_path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
    except OSError:
        return {}
    if proc.returncode != 0:
        return {}
    out = {}
    for line in proc.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        name = length = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                try:
                    length = int(field[3:])
                except ValueError:
                    length = None
        if name and length:
            out[name] = (length, None)
    return out


def check_samtools(samtools):
    try:
        proc = subprocess.run([samtools, "--version"], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        warn_missing("samtools executable", samtools)
        return None
    first = (proc.stdout or proc.stderr or "").splitlines()
    return first[0].strip() if first else "samtools (version unknown)"


# --------------------------------------------------------------------------- #
# a10 context (optional)
# --------------------------------------------------------------------------- #
def resolve_in_indir(indir, exact_name, suffix):
    """Exact filename in indir first, then any file ending with suffix."""
    if not indir or not os.path.isdir(indir):
        return None
    cand = os.path.join(indir, exact_name)
    if os.path.exists(cand):
        return cand
    try:
        names = sorted(os.listdir(indir))
    except OSError:
        return None
    for name in names:
        if name.endswith(suffix) and os.path.isfile(os.path.join(indir, name)):
            return os.path.join(indir, name)
    return None


def read_commented_tsv(path):
    """Rows of a '#'-commented TSV as dicts. None when unreadable."""
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            lines = [ln for ln in fh if not ln.startswith("#")]
    except OSError:
        return None
    if not lines:
        return []
    try:
        return list(csv.DictReader(lines, delimiter="\t"))
    except csv.Error:
        return None


def load_a10_context(args):
    """(real sample -> a10 anon, (a10 anon, reference) -> a10 pair context).

    a10 mints its own S001..Snnn from whatever IT processed, so its ids are only
    safe to join through its key file. Without that key the a10_* columns stay
    NA rather than being guessed from an id that may not line up. Never fatal.
    """
    key_map = {}
    pair_ctx = {}
    key_path = args.a10_key or resolve_in_indir(args.indir, A10_KEY_NAME,
                                                A10_KEY_SUFFIX)
    if not key_path or not os.path.exists(key_path):
        warn_missing("a10 sample key", key_path or os.path.join(
            args.indir or "<indir>", A10_KEY_NAME))
    else:
        rows = read_commented_tsv(key_path)
        if rows is None:
            warn_missing("readable a10 sample key", key_path)
        else:
            for row in rows:
                real = (row.get("sample_real") or "").strip()
                anon = (row.get("sample_anon") or "").strip()
                if real and anon:
                    key_map[real] = anon

    pair_path = args.a10_pairs or resolve_in_indir(args.indir, A10_PAIR_NAME,
                                                   A10_PAIR_SUFFIX)
    if not pair_path or not os.path.exists(pair_path):
        warn_missing("a10 by-pair audit table", pair_path or os.path.join(
            args.indir or "<indir>", A10_PAIR_NAME))
    elif not key_map:
        print("NOTE: a10's by-pair table was found but its sample key was not, "
              "so a10 verdicts cannot be joined; a10_* columns stay NA")
    else:
        rows = read_commented_tsv(pair_path)
        if rows is None:
            warn_missing("readable a10 by-pair audit table", pair_path)
        else:
            for row in rows:
                anon = (row.get("sample_anon") or "").strip()
                rid = (row.get("reference_id") or "").strip()
                if not anon or not rid:
                    continue
                pair_ctx[(anon, rid)] = {
                    "verdict": (row.get("verdict") or "NA").strip(),
                    "mwf": (row.get("max_window_fraction") or "NA").strip(),
                }
    return key_map, pair_ctx


# --------------------------------------------------------------------------- #
# sequence forensics
# --------------------------------------------------------------------------- #
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(seq):
    return "".join(COMPLEMENT.get(ch, "N") for ch in reversed(seq.upper()))


def build_adapters():
    """[(name, sequence)] with reverse complements, deduplicated by sequence."""
    out = []
    seen = set()
    for name, seq in BASE_ADAPTERS:
        for nm, sq in ((name, seq), (name + "_rc", revcomp(seq))):
            if sq in seen:
                continue
            seen.add(sq)
            out.append((nm, sq))
    for name, seq in HOMOPOLYMER_ADAPTERS:
        if seq in seen:
            continue
        seen.add(seq)
        out.append((name, seq))
    return out


def adapter_hit(clip, side, adapters, max_mm, min_overlap):
    """Best adapter / primer match for one soft clip.

    A right clip is read-through starting at the alignment boundary, so it is
    compared 5'->3' against the adapter prefix; a left clip abuts the boundary
    from the other side, so its last bases are compared against the adapter
    suffix. Comparison runs over the overlapping length only.
    Returns (name, overlap, mismatches); (None, 0, None) when nothing matched.
    """
    best_name, best_ov, best_mm = None, 0, None
    if not clip:
        return best_name, best_ov, best_mm
    up = clip.upper()
    for name, seq in adapters:
        overlap = min(len(up), len(seq))
        if overlap < min_overlap:
            continue
        if side == "R":
            a, b = up[:overlap], seq[:overlap]
        else:
            a, b = up[len(up) - overlap:], seq[len(seq) - overlap:]
        mm = 0
        ok = True
        for i in range(overlap):
            if a[i] != b[i]:
                mm += 1
                if mm > max_mm:
                    ok = False
                    break
        if not ok:
            continue
        if overlap > best_ov or (overlap == best_ov and best_mm is not None
                                 and mm < best_mm):
            best_name, best_ov, best_mm = name, overlap, mm
    return best_name, best_ov, best_mm


def entropy3(seq):
    """Shannon entropy in bits over overlapping 3-mers. None below 3 bases.

    The ceiling is log2(len - 2), so short clips cannot score high; compare
    within a length, not across lengths.
    """
    n = len(seq)
    if n < 3:
        return None
    counts = {}
    for i in range(n - 2):
        k = seq[i:i + 3]
        counts[k] = counts.get(k, 0) + 1
    total = float(n - 2)
    h = 0.0
    for c in counts.values():
        p = c / total
        h -= p * math.log(p, 2)
    return h


def longest_homopolymer(seq):
    best = 0
    run = 0
    prev = ""
    for ch in seq:
        if ch == prev:
            run += 1
        else:
            run = 1
            prev = ch
        if run > best:
            best = run
    return best


def gc_fraction(seq):
    gc = 0
    acgt = 0
    for ch in seq:
        if ch in "GCgc":
            gc += 1
            acgt += 1
        elif ch in "ATat":
            acgt += 1
    return safe_frac(gc, acgt)


def dust_score(seq):
    """SDUST-style triplet score normalised to 0-1 (1.0 = one repeated triplet).

    raw = sum(c * (c - 1) / 2) / (t - 1) over t = len - 2 triplets, which grows
    with length; dividing by its homopolymer maximum t/2 makes clips of
    different lengths comparable.
    """
    n = len(seq)
    t = n - 2
    if t < 2:
        return None
    counts = {}
    for i in range(t):
        k = seq[i:i + 3]
        counts[k] = counts.get(k, 0) + 1
    s = sum(c * (c - 1) / 2.0 for c in counts.values())
    raw = s / float(t - 1)
    top = t / 2.0
    if top <= 0:
        return None
    return min(1.0, raw / top)


def segment_stats(seq, min_entropy, max_hp_frac):
    """Complexity of one segment: length, entropy3, homopolymer, GC, DUST."""
    n = len(seq)
    if n <= 0:
        return None
    up = seq.upper()
    ent = entropy3(up)
    hp = float(longest_homopolymer(up)) / float(n)
    low = (ent is not None and ent < min_entropy) or (hp > max_hp_frac)
    return {
        "len": n,
        "entropy3": ent,
        "homopolymer_frac": hp,
        "gc": gc_fraction(up),
        "dust": dust_score(up),
        "low_complexity": bool(low),
    }


def clip_key_of(clip, side, key_len):
    """Clip key, anchored at the alignment boundary.

    Right clips are truncated to their first key_len bases and left clips to
    their last key_len bases, so the same adapter read-through collapses to one
    key whatever the clip length. key_len 0 keeps the whole clip.
    """
    up = clip.upper()
    if key_len <= 0 or len(up) <= key_len:
        return up
    return up[:key_len] if side == "R" else up[len(up) - key_len:]


# --------------------------------------------------------------------------- #
# CIGAR
# --------------------------------------------------------------------------- #
def parse_cigar(cigar):
    """[(length, op), ...]; [] for '*' or a malformed string."""
    if not cigar or cigar == "*":
        return []
    ops = []
    num = 0
    seen = False
    for ch in cigar:
        if "0" <= ch <= "9":
            num = num * 10 + (ord(ch) - 48)
            seen = True
            continue
        if not seen:
            return []
        ops.append((num, ch))
        num = 0
        seen = False
    if seen:
        return []
    return ops


def cigar_parts(ops):
    """(reference bases consumed, leading soft clip, trailing soft clip)."""
    span = 0
    for num, ch in ops:
        if ch in "M=XDN":
            span += num
    lead = 0
    for num, ch in ops:
        if ch == "H":
            continue
        if ch == "S":
            lead = num
        break
    trail = 0
    for num, ch in reversed(ops):
        if ch == "H":
            continue
        if ch == "S":
            trail = num
        break
    return span, lead, trail


def parse_tags(fields):
    """(AS, XS, XA list of alternate reference names). Missing tags -> None/[]."""
    a_score = x_score = None
    xa = []
    for tag in fields[11:]:
        if tag.startswith("AS:i:"):
            try:
                a_score = int(tag[5:].strip())
            except ValueError:
                a_score = None
        elif tag.startswith("XS:i:"):
            try:
                x_score = int(tag[5:].strip())
            except ValueError:
                x_score = None
        elif tag.startswith("XA:Z:"):
            for hit in tag[5:].split(";"):
                hit = hit.strip()
                if not hit:
                    continue
                xa.append(hit.split(",")[0])
    return a_score, x_score, xa


# --------------------------------------------------------------------------- #
# accumulators
# --------------------------------------------------------------------------- #
def new_acc(ref_len=None):
    acc = {"ref_len": ref_len, "adapters": {}, "clip_keys": {}}
    for key in COUNT_KEYS:
        acc[key] = 0
    for key in LIST_KEYS:
        acc[key] = []
    return acc


def merge_acc(dst, src):
    for key in COUNT_KEYS:
        dst[key] += src[key]
    for key in LIST_KEYS:
        dst[key].extend(src[key])
    for key, val in src["adapters"].items():
        dst["adapters"][key] = dst["adapters"].get(key, 0) + val
    for key, val in src["clip_keys"].items():
        dst["clip_keys"][key] = dst["clip_keys"].get(key, 0) + val
    return dst


def merged(accs):
    out = new_acc(None)
    for acc in accs:
        merge_acc(out, acc)
    return out


def top_of(counter):
    """(key, count) of the commonest entry, ties broken by key. (None, 0) if empty."""
    if not counter:
        return None, 0
    key = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))[0]
    return key[0], key[1]


def derive(acc):
    """Every reported metric from one accumulator (pair, sample or group)."""
    n = acc["n_reads"]
    d = {}
    d["ref_len"] = acc.get("ref_len")
    d["n_reads"] = n
    d["n_reads_with_seq"] = acc["n_seq"]
    d["median_mapq"] = median(acc["mapqs"])
    d["frac_mapq0"] = safe_frac(acc["n_mapq0"], n)

    d["n_reads_clipped"] = acc["n_reads_clipped"]
    d["frac_reads_clipped"] = safe_frac(acc["n_reads_clipped"], n)
    d["n_clips"] = acc["n_clips"]
    d["median_clip_len"] = median(acc["clip_lens"])
    d["max_clip_len"] = max(acc["clip_lens"]) if acc["clip_lens"] else 0
    d["n_clips_adapter"] = acc["n_clips_adapter"]
    d["frac_clips_adapter"] = safe_frac(acc["n_clips_adapter"], acc["n_clips"])
    d["n_reads_adapter"] = acc["n_reads_adapter"]
    d["frac_reads_adapter"] = safe_frac(acc["n_reads_adapter"], n)
    name, count = top_of(acc["adapters"])
    d["top_adapter"] = name or "none"
    d["top_adapter_clips"] = count
    key, key_count = top_of(acc["clip_keys"])
    d["n_distinct_clip_keys"] = len(acc["clip_keys"])
    d["top_clip_key"] = key or "NA"
    d["top_clip_key_clips"] = key_count

    d["aln_median_len"] = median(acc["aln_len"])
    d["aln_median_entropy3"] = median(acc["aln_ent"])
    d["aln_median_homopolymer_frac"] = median(acc["aln_hp"])
    d["aln_median_gc"] = median(acc["aln_gc"])
    d["aln_median_dust"] = median(acc["aln_dust"])
    d["aln_frac_low_complexity"] = safe_frac(acc["aln_lowcomp"], acc["aln_n"])
    d["clip_median_len"] = median(acc["clip_len_detail"])
    d["clip_median_entropy3"] = median(acc["clip_ent"])
    d["clip_median_homopolymer_frac"] = median(acc["clip_hp"])
    d["clip_median_gc"] = median(acc["clip_gc"])
    d["clip_median_dust"] = median(acc["clip_dust"])
    d["clip_frac_low_complexity"] = safe_frac(acc["clip_lowcomp"], acc["clip_n"])

    for cls in MATE_CLASSES:
        d["n_mate_" + cls] = acc["mate_" + cls]
        d["frac_mate_" + cls] = safe_frac(acc["mate_" + cls], n)
    d["n_mate_unknown_ref"] = acc["n_mate_unknown"]

    d["n_with_as"] = acc["n_with_as"]
    d["n_with_xs"] = acc["n_with_xs"]
    d["n_as_eq_xs"] = acc["n_as_eq_xs"]
    d["frac_as_eq_xs"] = safe_frac(acc["n_as_eq_xs"], n)
    d["n_xs_gt_as"] = acc["n_xs_gt_as"]
    d["frac_xs_gt_as"] = safe_frac(acc["n_xs_gt_as"], n)
    d["n_with_xa"] = acc["n_with_xa"]
    d["frac_with_xa"] = safe_frac(acc["n_with_xa"], n)
    d["n_xa_hits"] = acc["n_xa_hits"]
    d["n_xa_hits_anello"] = acc["n_xa_anello"]
    d["frac_xa_hits_anello"] = safe_frac(acc["n_xa_anello"], acc["n_xa_hits"])
    d["median_xa_hits"] = median(acc["xa_hits"])
    d["median_xa_anello_hits"] = median(acc["xa_anello_hits"])

    d["n_utr_reads"] = acc["n_utr"]
    d["frac_utr_reads"] = safe_frac(acc["n_utr"], n)
    d["detail_capped"] = "1" if acc["n_detail_skipped"] else "0"
    return d


def pair_flags(d, args, n_recurrent_keys):
    """Human-readable flags; every criterion that fired, not just the top one."""
    flags = []

    def fired(value):
        return value is not None and value >= args.flag_frac

    if fired(d["frac_clips_adapter"]) and d["n_clips"] > 0:
        flags.append("adapter_clips")
    if fired(d["clip_frac_low_complexity"]):
        flags.append("clip_low_complexity")
    if fired(d["aln_frac_low_complexity"]):
        flags.append("aligned_low_complexity")
    if fired(d["frac_mate_HUMAN"]):
        flags.append("mate_human")
    if fired(d["frac_mate_other_anello_ref"]):
        flags.append("mate_other_anello_ref")
    if fired(d["frac_mate_unmapped"]):
        flags.append("mate_unmapped")
    if fired(d["frac_as_eq_xs"]):
        flags.append("equal_scoring_alternative")
    if fired(d["frac_with_xa"]):
        flags.append("multi_hit_xa")
    if fired(d["frac_utr_reads"]):
        flags.append("utr_window")
    if n_recurrent_keys > 0:
        flags.append("recurrent_clip_shared")
    return ";".join(flags) if flags else "none"


# --------------------------------------------------------------------------- #
# per-read streaming
# --------------------------------------------------------------------------- #
def read_key_of(qname):
    """Truncated SHA1 of QNAME.

    The read name itself is never written: it can embed a sample, run or flow
    cell identifier. The hash is only there so two rows can be recognised as
    the same read pair.
    """
    if not qname or qname == "*":
        return "NA"
    return hashlib.sha1(qname.encode("utf-8", "replace")).hexdigest()[:12]


def stream_reference_set(args, bam_path, ref_lens, chimp_of, classifier,
                         adapters, clip_index, read_fh, ctx):
    """Stream the human and the chimpanzee references of one BAM.

    They are streamed as two passes only so the chimp flag can be attached to
    every record without a per-read lookup; the accumulators come back in one
    dict keyed by reference, exactly as if it had been one call. Returns
    (accumulators, error_string), accumulators None when samtools failed.
    """
    out = {}
    for chimp_flag in (False, True):
        subset = dict((r, length) for r, length in ref_lens.items()
                      if bool(chimp_of.get(r, False)) == chimp_flag)
        if not subset:
            continue
        accs, err = stream_sample(args.samtools, bam_path, subset, chimp_flag,
                                  classifier, adapters, clip_index, read_fh,
                                  ctx, args)
        if accs is None:
            return None, err
        out.update(accs)
    return out, ""


def stream_sample(samtools, bam, ref_lens, chimp_flag, classifier, adapters,
                  clip_index, read_fh, ctx, args, chunk=200):
    """Stream one BAM over every wanted anellovirus reference.

    Returns (accumulators keyed by reference_id, error_string); the
    accumulators are None when samtools could not be run. ctx carries the run
    label, the anonymous id, the group and the real sample name; the real name
    is used ONLY to mask samtools' error text before it is printed.
    """
    sample = ctx["sample"]
    anon = ctx["sample_anon"]
    accs = {}
    ref_list = sorted(ref_lens)
    for start in range(0, len(ref_list), chunk):
        block = ref_list[start:start + chunk]
        cmd = [samtools, "view", "-F", str(args.exclude_flags),
               "-q", str(args.min_mapq)]
        if args.samtools_threads > 1:
            cmd += ["-@", str(args.samtools_threads)]
        cmd += [bam] + block
        # stderr goes to a temp file, not a pipe: it is only drained after the
        # stdout loop finishes, and a pipe that filled up first would deadlock.
        # Binary mode on purpose: TemporaryFile only grew an "errors" keyword in
        # Python 3.8, and this suite supports 3.7, so decode by hand below.
        err_fh = tempfile.TemporaryFile(mode="w+b")
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err_fh,
                                    universal_newlines=True, errors="replace")
        except OSError as exc:
            err_fh.close()
            return None, "samtools not runnable (%s)" % strip_paths(
                redact(exc, sample, anon))
        for line in proc.stdout:
            if not line or line[0] == "@":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            rname = f[2]
            ref_len = ref_lens.get(rname)
            if not ref_len or ref_len < 1:
                continue
            try:
                flag = int(f[1])
                pos0 = int(f[3]) - 1
                mapq = int(f[4])
            except ValueError:
                continue
            if pos0 < 0:
                pos0 = 0
            acc = accs.get(rname)
            if acc is None:
                acc = new_acc(ref_len)
                accs[rname] = acc
            handle_read(f, flag, pos0, mapq, rname, ref_len, acc, chimp_flag,
                        classifier, adapters, clip_index, read_fh, ctx, args)
        proc.stdout.close()
        rc = proc.wait()
        try:
            err_fh.seek(0)
            err = err_fh.read().decode("utf-8", "replace")
        except (IOError, OSError, ValueError, UnicodeError):
            err = ""
        err_fh.close()
        if rc != 0:
            # Redact FIRST, on the raw text: samtools echoes the BAM path, and
            # either truncating it or running it through to_ascii() (which maps
            # anything outside printable ASCII to "?") could break the real name
            # into a form a later replace() would no longer match. Then drop the
            # path itself, which names the controlled-data mount and is not the
            # sample name, so redact() alone would let it through.
            err = strip_paths(redact(err, sample, anon))
            err = " ".join(to_ascii(err.replace("\r", " ").replace("\n", " ")).split())
            return None, "samtools view rc=%d %s" % (rc, err[:160] or "no stderr")
    return accs, ""


def handle_read(f, flag, pos0, mapq, rname, ref_len, acc, chimp_flag,
                classifier, adapters, clip_index, read_fh, ctx, args):
    """All five forensic measurements for one SAM record."""
    acc["n_reads"] += 1
    acc["mapqs"].append(mapq)
    if mapq == 0:
        acc["n_mapq0"] += 1

    cigar = f[5]
    ops = parse_cigar(cigar)
    span, lead, trail = cigar_parts(ops)
    seq = f[9]
    has_seq = bool(seq) and seq != "*"
    if has_seq:
        acc["n_seq"] += 1

    # ---- (5) UTR overlap --------------------------------------------------- #
    rel = None
    if span > 0:
        rel = (float(pos0) + float(span) / 2.0) / float(ref_len)
    else:
        rel = float(pos0) / float(ref_len)
    rel = min(1.0, max(0.0, rel))
    in_utr = args.utr_lo <= rel <= args.utr_hi
    if in_utr:
        acc["n_utr"] += 1

    # ---- (3) mate origin --------------------------------------------------- #
    rnext = f[6]
    try:
        pnext = int(f[7])
    except (ValueError, IndexError):
        pnext = 0
    if not (flag & 0x1):
        mate_class = "no_mate"
        mate_ref = "NA"
    elif (flag & 0x8) or rnext in ("*", ""):
        mate_class = "unmapped"
        mate_ref = "NA"
    elif rnext == "=" or rnext == rname:
        mate_class = "same_anello_ref"
        mate_ref = rname
    else:
        mate_ref = rnext
        cls = classifier.classify(rnext)
        if cls == "ANELLO_CHIMP":
            mate_class = "chimp_anello_ref"
        elif cls == "ANELLO":
            mate_class = "other_anello_ref"
        elif cls == "HUMAN":
            mate_class = "HUMAN"
        else:
            mate_class = "other_category"
            if cls == "UNKNOWN":
                acc["n_mate_unknown"] += 1
    acc["mate_" + mate_class] += 1

    # ---- (4) alignment ambiguity ------------------------------------------- #
    a_score, x_score, xa = parse_tags(f)
    if a_score is not None:
        acc["n_with_as"] += 1
    if x_score is not None:
        acc["n_with_xs"] += 1
    as_eq_xs = (a_score is not None and x_score is not None
                and a_score == x_score)
    if as_eq_xs:
        acc["n_as_eq_xs"] += 1
    # strictly better elsewhere: this read is on this reference only because
    # the aligner had to choose one, which is worse than a tie
    if a_score is not None and x_score is not None and x_score > a_score:
        acc["n_xs_gt_as"] += 1
    n_xa_anello = 0
    n_xa_chimp = 0
    if xa:
        acc["n_with_xa"] += 1
        acc["n_xa_hits"] += len(xa)
        for name in xa:
            cls = classifier.classify(name)
            if cls in ("ANELLO", "ANELLO_CHIMP"):
                n_xa_anello += 1
                if cls == "ANELLO_CHIMP":
                    n_xa_chimp += 1
        acc["n_xa_anello"] += n_xa_anello
        acc["n_xa_chimp"] += n_xa_chimp

    # ---- (1) soft clips and (2) complexity --------------------------------- #
    detail_room = acc["n_detail"] < args.max_detail_reads
    if detail_room:
        acc["n_detail"] += 1
        if xa:
            # medians over the reads that actually carry an XA tag; the reads
            # without one are in frac_with_xa, not silently counted as zero
            acc["xa_hits"].append(len(xa))
            acc["xa_anello_hits"].append(n_xa_anello)
    else:
        acc["n_detail_skipped"] += 1

    aln_stats = None
    clip_stats = None
    left_key = right_key = "NA"
    hit_name = hit_side = "NA"
    hit_ov = hit_mm = None
    clip_total = 0

    if has_seq and lead + trail <= len(seq):
        aln_seq = seq[lead:len(seq) - trail] if (lead or trail) else seq
        if aln_seq:
            aln_stats = segment_stats(aln_seq, args.min_entropy,
                                      args.max_homopolymer_frac)
        if aln_stats:
            acc["aln_n"] += 1
            if aln_stats["low_complexity"]:
                acc["aln_lowcomp"] += 1
            if detail_room:
                acc["aln_len"].append(aln_stats["len"])
                if aln_stats["entropy3"] is not None:
                    acc["aln_ent"].append(aln_stats["entropy3"])
                acc["aln_hp"].append(aln_stats["homopolymer_frac"])
                if aln_stats["gc"] is not None:
                    acc["aln_gc"].append(aln_stats["gc"])
                if aln_stats["dust"] is not None:
                    acc["aln_dust"].append(aln_stats["dust"])

        clips = []
        if lead >= args.min_clip:
            clips.append(("L", seq[:lead]))
        if trail >= args.min_clip:
            clips.append(("R", seq[len(seq) - trail:]))
        if clips:
            acc["n_reads_clipped"] += 1
        keys_this_read = set()
        read_had_adapter = False
        for side, clip in clips:
            acc["n_clips"] += 1
            acc["clip_lens"].append(len(clip))
            clip_total += len(clip)
            name, ov, mm = adapter_hit(clip, side, adapters,
                                       args.adapter_mismatch,
                                       args.min_adapter_overlap)
            if name:
                acc["n_clips_adapter"] += 1
                acc["adapters"][name] = acc["adapters"].get(name, 0) + 1
                read_had_adapter = True
                if hit_name == "NA" or (hit_ov is not None and ov > hit_ov):
                    hit_name, hit_side, hit_ov, hit_mm = name, side, ov, mm
            key = clip_key_of(clip, side, args.clip_key_len)
            if side == "L":
                left_key = key
            else:
                right_key = key
            acc["clip_keys"][key] = acc["clip_keys"].get(key, 0) + 1
            index_clip(clip_index, key, side, len(clip), name, ov, mm, ctx,
                       rname, chimp_flag, key in keys_this_read, args)
            keys_this_read.add(key)
        if read_had_adapter:
            acc["n_reads_adapter"] += 1
        if clips:
            joined = "".join(c for _s, c in clips)
            clip_stats = segment_stats(joined, args.min_entropy,
                                       args.max_homopolymer_frac)
            if clip_stats:
                acc["clip_n"] += 1
                if clip_stats["low_complexity"]:
                    acc["clip_lowcomp"] += 1
                if detail_room:
                    acc["clip_len_detail"].append(clip_stats["len"])
                    if clip_stats["entropy3"] is not None:
                        acc["clip_ent"].append(clip_stats["entropy3"])
                    acc["clip_hp"].append(clip_stats["homopolymer_frac"])
                    if clip_stats["gc"] is not None:
                        acc["clip_gc"].append(clip_stats["gc"])
                    if clip_stats["dust"] is not None:
                        acc["clip_dust"].append(clip_stats["dust"])

    if read_fh is not None:
        write_read_row(read_fh, f, flag, pos0, mapq, rname, span, rel, in_utr,
                       cigar, lead, trail, hit_name, hit_side, hit_ov, hit_mm,
                       left_key, right_key, aln_stats, clip_stats, clip_total,
                       mate_class, mate_ref, pnext, a_score, x_score, as_eq_xs,
                       xa, n_xa_anello, n_xa_chimp, chimp_flag, ctx)


def index_clip(clip_index, key, side, clip_len, adapter, overlap, mismatches,
               ctx, rname, chimp_flag, seen_in_read, args):
    """Add one clip to the cohort-wide clip index."""
    entry = clip_index["keys"].get(key)
    if entry is None:
        if len(clip_index["keys"]) >= args.max_clip_keys:
            clip_index["dropped"] += 1
            return
        entry = {
            "n_clips": 0, "n_reads": 0, "samples": set(), "refs": set(),
            "groups": set(), "chimp_refs": set(), "sides": set(),
            "lens": [], "adapter": (adapter, overlap, mismatches),
            "example_ref": rname,
        }
        clip_index["keys"][key] = entry
    entry["n_clips"] += 1
    if not seen_in_read:
        entry["n_reads"] += 1
    entry["samples"].add(ctx["sample_anon"])
    entry["refs"].add(rname)
    entry["groups"].add(ctx["group"])
    if chimp_flag:
        entry["chimp_refs"].add(rname)
    entry["sides"].add(side)
    if len(entry["lens"]) < 10000:
        entry["lens"].append(clip_len)
    if adapter and entry["adapter"][0] is None:
        entry["adapter"] = (adapter, overlap, mismatches)


def write_read_row(fh, f, flag, pos0, mapq, rname, span, rel, in_utr, cigar,
                   lead, trail, hit_name, hit_side, hit_ov, hit_mm, left_key,
                   right_key, aln_stats, clip_stats, clip_total, mate_class,
                   mate_ref, pnext, a_score, x_score, as_eq_xs, xa,
                   n_xa_anello, n_xa_chimp, chimp_flag, ctx):
    def seg(stats, field, digits=4):
        if not stats or stats.get(field) is None:
            return "NA"
        value = stats[field]
        if isinstance(value, bool):
            return "1" if value else "0"
        if isinstance(value, int):
            return str(value)
        return fnum(value, digits)

    row = [
        ctx["run"], ctx["sample_anon"], ctx["group"], rname,
        "1" if chimp_flag else "0", read_key_of(f[0]), str(flag), str(mapq),
        str(pos0 + 1), str(span), fnum(rel, 4), "1" if in_utr else "0", cigar,
        str(lead), str(trail), hit_name, hit_side,
        "NA" if hit_ov is None else str(hit_ov),
        "NA" if hit_mm is None else str(hit_mm), left_key, right_key,
        seg(aln_stats, "len"), seg(aln_stats, "entropy3", 3),
        seg(aln_stats, "homopolymer_frac", 3), seg(aln_stats, "gc", 3),
        seg(aln_stats, "dust", 3), seg(aln_stats, "low_complexity"),
        str(clip_total), seg(clip_stats, "entropy3", 3),
        seg(clip_stats, "homopolymer_frac", 3), seg(clip_stats, "gc", 3),
        seg(clip_stats, "dust", 3), seg(clip_stats, "low_complexity"),
        mate_class, mate_ref, str(pnext),
        "NA" if a_score is None else str(a_score),
        "NA" if x_score is None else str(x_score),
        "1" if as_eq_xs else "0", str(len(xa)), str(n_xa_anello),
        str(n_xa_chimp),
    ]
    fh.write("\t".join(to_ascii(c) for c in row) + "\n")


# --------------------------------------------------------------------------- #
# comments shared by every output
# --------------------------------------------------------------------------- #
def common_comments(args, extra=()):
    lines = [
        "%s  generated %s" % (SCRIPT, TODAY),
        "no real sample identifiers in this file; ids are anonymous S001..Snnn "
        "(mapping in %s)" % out_name(args.prefix, "sample_key.tsv"),
        "read filter: samtools view -F %s -q %d (MAPQ is deliberately NOT "
        "filtered by default: multi-mapping reads are the artefact under "
        "investigation)" % (args.exclude_flags, args.min_mapq),
        "clips >= %d bp; adapter match allows %d mismatch(es) over an overlap "
        "of >= %d bp; clip keys are anchored at the alignment boundary and cut "
        "to %s"
        % (args.min_clip, args.adapter_mismatch, args.min_adapter_overlap,
           ("%d bp" % args.clip_key_len) if args.clip_key_len > 0
           else "their full length"),
        "low complexity = 3-mer entropy < %.2f bits or homopolymer fraction > "
        "%.2f; entropy over 3-mers is capped at log2(len-2), so short clips "
        "cannot score high - compare within a length"
        % (args.min_entropy, args.max_homopolymer_frac),
        "UTR window = relative %.2f-%.2f of each reference in its own "
        "orientation; reported, never assumed - check the annotation before "
        "calling it the UTR" % (args.utr_lo, args.utr_hi),
        "nothing is re-aligned and no database is consulted: a read that fails "
        "every test here is NOT thereby real virus, it only means none of the "
        "five tested mechanisms fired",
    ]
    lines.extend(extra)
    return lines


def write_sample_key(path, samples):
    """samples: list of (anon, real, group, runs, a10_anon)."""
    lines = [
        "CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
        "generated %s by %s" % (TODAY, SCRIPT),
    ]
    header = ["sample_anon", "sample_real", "group", "runs", "a10_sample_anon"]
    write_tsv(path, lines, header, samples)


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def write_by_pair(path, pairs, args):
    rows = []
    for p in sorted(pairs, key=lambda r: (r["group"], r["sample_anon"],
                                          r["reference_id"])):
        d = p["m"]
        rows.append([
            p["run"], p["sample_anon"], p["group"], p["reference_id"],
            p["ref_label"], str(d["ref_len"]), "1" if p["chimp"] else "0",
            "NA" if p["idxstats_mapped"] is None else str(p["idxstats_mapped"]),
            str(d["n_reads"]), str(d["n_reads_with_seq"]),
            fnum(d["median_mapq"], 1), fnum(d["frac_mapq0"], 4),
            str(d["n_reads_clipped"]), fnum(d["frac_reads_clipped"], 4),
            str(d["n_clips"]), fnum(d["median_clip_len"], 1),
            str(d["max_clip_len"]), str(d["n_clips_adapter"]),
            fnum(d["frac_clips_adapter"], 4), str(d["n_reads_adapter"]),
            fnum(d["frac_reads_adapter"], 4), d["top_adapter"],
            str(d["top_adapter_clips"]), str(d["n_distinct_clip_keys"]),
            d["top_clip_key"], str(d["top_clip_key_clips"]),
            str(p["n_recurrent_keys"]), str(p["n_clips_recurrent"]),
            fnum(d["aln_median_len"], 1), fnum(d["aln_median_entropy3"], 3),
            fnum(d["aln_median_homopolymer_frac"], 3),
            fnum(d["aln_median_gc"], 3), fnum(d["aln_median_dust"], 3),
            fnum(d["aln_frac_low_complexity"], 4),
            fnum(d["clip_median_len"], 1), fnum(d["clip_median_entropy3"], 3),
            fnum(d["clip_median_homopolymer_frac"], 3),
            fnum(d["clip_median_gc"], 3), fnum(d["clip_median_dust"], 3),
            fnum(d["clip_frac_low_complexity"], 4),
            str(d["n_mate_same_anello_ref"]), str(d["n_mate_other_anello_ref"]),
            str(d["n_mate_chimp_anello_ref"]), str(d["n_mate_HUMAN"]),
            str(d["n_mate_other_category"]), str(d["n_mate_unmapped"]),
            str(d["n_mate_no_mate"]), str(d["n_mate_unknown_ref"]),
            fnum(d["frac_mate_same_anello_ref"], 4),
            fnum(d["frac_mate_other_anello_ref"], 4),
            fnum(d["frac_mate_HUMAN"], 4),
            str(d["n_with_as"]), str(d["n_with_xs"]), str(d["n_as_eq_xs"]),
            fnum(d["frac_as_eq_xs"], 4), str(d["n_xs_gt_as"]),
            fnum(d["frac_xs_gt_as"], 4), str(d["n_with_xa"]),
            fnum(d["frac_with_xa"], 4), str(d["n_xa_hits"]),
            str(d["n_xa_hits_anello"]), fnum(d["frac_xa_hits_anello"], 4),
            fnum(d["median_xa_hits"], 2), fnum(d["median_xa_anello_hits"], 2),
            str(d["n_utr_reads"]), fnum(d["frac_utr_reads"], 4),
            fnum(args.utr_lo, 2), fnum(args.utr_hi, 2),
            d["detail_capped"], p["a10_anon"], p["a10_verdict"], p["a10_mwf"],
            p["flags"],
        ])
    write_tsv(path, common_comments(args, [
        "one row per (sample, anellovirus reference) carrying at least one read",
        "chimp_flagged=1 rows are the negative control: no human sample can "
        "carry chimpanzee TTV, so whatever they collect is cross-mapping and "
        "their metrics are the artefact baseline the human rows are read "
        "against",
        "frac_* columns are over that pair's reads; medians over the per-read "
        "values of that pair (the first %d reads in coordinate order if "
        "detail_capped=1)" % args.max_detail_reads,
        "a read with no SEQ (n_reads_with_seq below n_reads) cannot be clip- "
        "or complexity-analysed; it still counts in the mate, ambiguity, MAPQ "
        "and UTR columns",
        "as_eq_xs counts reads with an EQUALLY good alignment elsewhere; "
        "xs_gt_as counts reads with a BETTER one, i.e. reads sitting on this "
        "reference only because the aligner had to pick a primary",
        "a10_* columns are context copied from a10's by-pair table through "
        "a10's own sample key; they are not recomputed here",
        "this table cannot say where a read really came from - nothing is "
        "re-aligned - only which artefact mechanisms it is consistent with",
    ]), PAIR_COLUMNS, rows)


def blank_group_row():
    return dict((c, "NA") for c in GROUP_COLUMNS)


def fill_summary_row(row, d, label, ref_set, n_samples, n_pairs, clip_stats):
    row["row_type"] = "group_summary"
    row["label"] = label
    row["ref_set"] = ref_set
    row["n_samples"] = str(n_samples)
    row["n_pairs"] = str(n_pairs)
    row["n_reads"] = str(d["n_reads"])
    row["n_reads_with_seq"] = str(d["n_reads_with_seq"])
    row["frac_reads_clipped"] = fnum(d["frac_reads_clipped"], 4)
    row["n_clips"] = str(d["n_clips"])
    row["frac_clips_adapter"] = fnum(d["frac_clips_adapter"], 4)
    row["frac_reads_adapter"] = fnum(d["frac_reads_adapter"], 4)
    row["top_adapter"] = d["top_adapter"]
    row["aln_median_entropy3"] = fnum(d["aln_median_entropy3"], 3)
    row["aln_frac_low_complexity"] = fnum(d["aln_frac_low_complexity"], 4)
    row["clip_median_entropy3"] = fnum(d["clip_median_entropy3"], 3)
    row["clip_frac_low_complexity"] = fnum(d["clip_frac_low_complexity"], 4)
    for cls in MATE_CLASSES:
        row["frac_mate_" + cls] = fnum(d["frac_mate_" + cls], 4)
    row["frac_as_eq_xs"] = fnum(d["frac_as_eq_xs"], 4)
    row["frac_with_xa"] = fnum(d["frac_with_xa"], 4)
    row["frac_xa_hits_anello"] = fnum(d["frac_xa_hits_anello"], 4)
    row["frac_utr_reads"] = fnum(d["frac_utr_reads"], 4)
    row["n_clip_keys"] = str(clip_stats["n_keys"])
    row["n_recurrent_clip_keys"] = str(clip_stats["n_recurrent"])
    row["top_clip_key"] = clip_stats["top_key"]
    row["top_clip_key_clips"] = str(clip_stats["top_clips"])
    row["top_clip_key_samples"] = str(clip_stats["top_samples"])
    row["top_clip_key_refs"] = str(clip_stats["top_refs"])
    return row


def scope_clip_stats(acc, clip_index, args):
    """Clip-key summary for one scope, resolved against the cohort-wide index."""
    keys = acc["clip_keys"]
    n_recurrent = 0
    for key in keys:
        entry = clip_index["keys"].get(key)
        if entry and is_recurrent(entry, args):
            n_recurrent += 1
    top_key, top_clips = top_of(keys)
    entry = clip_index["keys"].get(top_key) if top_key else None
    return {
        "n_keys": len(keys),
        "n_recurrent": n_recurrent,
        "top_key": top_key or "NA",
        "top_clips": top_clips,
        "top_samples": len(entry["samples"]) if entry else 0,
        "top_refs": len(entry["refs"]) if entry else 0,
    }


def is_recurrent(entry, args):
    return (len(entry["samples"]) >= args.recurrent_min_samples
            and len(entry["refs"]) >= args.recurrent_min_refs)


def build_group_rows(pairs_by_set, clip_index, args, g1, g2):
    """(rows for by_group.tsv, per-scope summaries, tests for the stdout block)."""
    rows = []
    summaries = {}
    per_sample = {}          # ref_set -> sample_anon -> {"group", "d"}
    for ref_set in ("human_anello", "chimp_flagged"):
        pairs = pairs_by_set.get(ref_set, [])
        if not pairs:
            continue
        by_sample = {}
        for p in pairs:
            by_sample.setdefault(p["sample_anon"], []).append(p)
        per_sample[ref_set] = {}
        for anon, plist in by_sample.items():
            acc = merged([p["acc"] for p in plist])
            per_sample[ref_set][anon] = {"group": plist[0]["group"],
                                         "d": derive(acc)}
        by_group = {}
        for p in pairs:
            by_group.setdefault(p["group"], []).append(p)
        by_group["ALL"] = list(pairs)
        labels = [g for g in GROUP_ORDER if g in by_group]
        labels += [g for g in sorted(by_group) if g not in labels and g != "ALL"]
        labels.append("ALL")
        for label in labels:
            plist = by_group[label]
            acc = merged([p["acc"] for p in plist])
            d = derive(acc)
            n_samples = len(set(p["sample_anon"] for p in plist))
            row = fill_summary_row(blank_group_row(), d, label, ref_set,
                                   n_samples, len(plist),
                                   scope_clip_stats(acc, clip_index, args))
            row["note"] = ("fractions are pooled over reads in this scope, so "
                           "the sample with the most reads dominates; the "
                           "tests below use one value per sample")
            rows.append([row[c] for c in GROUP_COLUMNS])
            summaries[(ref_set, label)] = dict(row)

    tests = []

    # ---- human anellovirus references vs the chimpanzee control ------------ #
    human = per_sample.get("human_anello", {})
    chimp = per_sample.get("chimp_flagged", {})
    shared = sorted(set(human) & set(chimp))
    if not shared:
        # No control to contrast against: either --include-chimp folded the
        # chimpanzee references in, or no sample carried reads on both sets.
        # One explanatory row beats 18 rows of NA.
        row = blank_group_row()
        row["row_type"] = "refset_contrast"
        row["label"] = "not_run"
        row["ref_set"] = "human_anello_vs_chimp_flagged"
        row["n1"] = "0"
        row["n2"] = "0"
        row["note"] = ("no sample carried reads on both reference sets, so the "
                       "negative-control contrast was not run"
                       + (" (--include-chimp folded the chimpanzee references "
                          "into the human aggregates)" if args.include_chimp
                          else ""))
        rows.append([row[c] for c in GROUP_COLUMNS])
    for name, field in (CONTRAST_METRICS if shared else []):
        values = [(human[a]["d"].get(field), chimp[a]["d"].get(field))
                  for a in shared]
        res = wilcoxon_signed_rank(values)
        row = blank_group_row()
        row["row_type"] = "refset_contrast"
        row["label"] = name
        row["ref_set"] = "human_anello_vs_chimp_flagged"
        row["group1"] = "human_anello"
        row["group2"] = "chimp_flagged"
        row["statistic"] = "wilcoxon_signed_rank_W_plus"
        usable = [v for v in values if v[0] is not None and v[1] is not None]
        row["n1"] = str(len(usable))
        row["n2"] = str(len(usable))
        # The medians are filled in even when the test cannot run, because
        # "no non-zero difference" means the two reference sets were IDENTICAL
        # on this metric - the strongest form of the negative-control reading -
        # and an NA there would read as missing data instead.
        med1 = median([v[0] for v in usable])
        med2 = median([v[1] for v in usable])
        row["value_group1"] = fnum(med1, 4)
        row["value_group2"] = fnum(med2, 4)
        if res is None:
            row["p_value"] = "NA"
            row["effect"] = fnum(0.0, 4) if usable else "NA"
            row["note"] = ("not tested: %d sample(s) carried the metric on "
                           "both reference sets and %s; a test needs at least "
                           "one non-zero difference"
                           % (len(usable),
                              "every one of them gave the identical value on "
                              "the human and the chimpanzee references"
                              if usable else "none had it on both"))
        else:
            row["value_group1"] = fnum(res["median1"], 4)
            row["value_group2"] = fnum(res["median2"], 4)
            row["statistic_value"] = fnum(res["W_plus"], 1)
            row["p_value"] = fp(res["p"])
            row["effect"] = fnum(res["median_diff"], 4)
            note = ("paired over the %d sample(s) carrying reads on both "
                    "reference sets, %d zero difference(s) dropped; a LARGE p "
                    "means the human references are indistinguishable from the "
                    "chimpanzee negative control on this metric"
                    % (res["n_pairs"], res["n_zero"]))
            if res["n_nonzero"] < 6:
                note += ("; n=%d non-zero pairs, the normal approximation is "
                         "unreliable this small" % res["n_nonzero"])
            row["note"] = note
        rows.append([row[c] for c in GROUP_COLUMNS])
        tests.append(("contrast", name, res, len(usable), len(usable),
                      (med1, med2)))

    # ---- HIV vs HL on the human references --------------------------------- #
    s1 = [v for v in human.values() if v["group"] == g1]
    s2 = [v for v in human.values() if v["group"] == g2]
    for name, field in CONTRAST_METRICS:
        x = [v["d"][field] for v in s1 if v["d"].get(field) is not None]
        y = [v["d"][field] for v in s2 if v["d"].get(field) is not None]
        res = mann_whitney_u(x, y)
        row = blank_group_row()
        row["row_type"] = "group_test"
        row["label"] = name
        row["ref_set"] = "human_anello"
        row["group1"] = g1
        row["group2"] = g2
        row["n1"] = str(len(x))
        row["n2"] = str(len(y))
        row["statistic"] = "mann_whitney_U_group1"
        if res is None:
            row["p_value"] = "NA"
            row["note"] = ("not tested: %d %s sample(s) and %d %s sample(s) "
                           "had this metric" % (len(x), g1, len(y), g2))
        else:
            row["value_group1"] = fnum(res["median1"], 4)
            row["value_group2"] = fnum(res["median2"], 4)
            row["statistic_value"] = fnum(res["U1"], 1)
            row["p_value"] = fp(res["p"])
            row["effect"] = fnum(res["effect_r"], 4)
            note = ("two-sided normal approximation, tie- and "
                    "continuity-corrected, over one value per sample")
            if min(res["n1"], res["n2"]) < 5:
                note += "; n<5 in one group, p is approximate and underpowered"
            row["note"] = note
        rows.append([row[c] for c in GROUP_COLUMNS])
        tests.append(("group", name, res, len(x), len(y), None))

    # ---- HIV vs HL, samples carrying any of each artefact signature -------- #
    binary = [
        ("any_adapter_matched_clip", lambda d: (d["n_clips_adapter"] or 0) > 0),
        ("any_mate_in_human_genome", lambda d: (d["n_mate_HUMAN"] or 0) > 0),
        ("any_equal_scoring_alternative", lambda d: (d["n_as_eq_xs"] or 0) > 0),
        ("any_read_in_utr_window", lambda d: (d["n_utr_reads"] or 0) > 0),
    ]
    for name, test in binary:
        a = sum(1 for v in s1 if test(v["d"]))
        b = len(s1) - a
        c = sum(1 for v in s2 if test(v["d"]))
        e = len(s2) - c
        res = fisher_exact_2x2(a, b, c, e)
        row = blank_group_row()
        row["row_type"] = "group_test"
        row["label"] = name
        row["ref_set"] = "human_anello"
        row["group1"] = g1
        row["group2"] = g2
        row["n1"] = str(len(s1))
        row["n2"] = str(len(s2))
        row["value_group1"] = "%d/%d" % (a, len(s1))
        row["value_group2"] = "%d/%d" % (c, len(s2))
        row["statistic"] = "fisher_exact_2x2_odds_ratio"
        if res is None:
            row["p_value"] = "NA"
            row["note"] = ("not tested: degenerate 2x2 table [[%d,%d],[%d,%d]]"
                           % (a, b, c, e))
        else:
            row["statistic_value"] = ("inf" if res["odds_ratio"] == float("inf")
                                      else fnum(res["odds_ratio"], 4))
            row["p_value"] = fp(res["p"])
            row["note"] = ("two-sided Fisher exact on [[%d,%d],[%d,%d]]; "
                           "samples counted once each, not pairs"
                           % (a, b, c, e))
        rows.append([row[c] for c in GROUP_COLUMNS])
        tests.append(("fisher", name, res, len(s1), len(s2), None))

    return rows, summaries, tests, per_sample


def write_by_group(path, rows, args, g1, g2):
    write_tsv(path, common_comments(args, [
        "row_type=group_summary: one row per (group, reference set); "
        "ref_set=chimp_flagged is the chimpanzee negative control",
        "row_type=refset_contrast: human anellovirus references vs the "
        "chimpanzee control on EVERY metric, paired Wilcoxon signed-rank over "
        "the samples carrying both. A LARGE p is the finding here: it means "
        "the human references behave like the negative control",
        "row_type=group_test: %s vs %s, Mann-Whitney over per-sample values "
        "and Fisher exact over samples; n1/n2 sit beside every p value"
        % (g1, g2),
        "no p value here is corrected for multiple testing: %d metrics are "
        "reported twice over plus %d Fisher rows, so read them as a profile, "
        "not as %d independent tests" % (len(CONTRAST_METRICS), 4,
                                         2 * len(CONTRAST_METRICS) + 4),
        "pooled group fractions are descriptive: reads within one sample are "
        "not independent observations",
    ]), GROUP_COLUMNS, rows)


def build_clip_rows(clip_index, args):
    """(rows for the clip table, ordered entries, n recurrent keys overall)."""
    items = sorted(clip_index["keys"].items(),
                   key=lambda kv: (-kv[1]["n_clips"], -len(kv[1]["samples"]),
                                   kv[0]))
    n_recurrent = sum(1 for _k, e in items if is_recurrent(e, args))
    rows = []
    shown = items[:args.top_clips] if args.top_clips > 0 else items
    for rank, (key, e) in enumerate(shown, start=1):
        stats = segment_stats(key, args.min_entropy, args.max_homopolymer_frac)
        adapter, overlap, mismatches = e["adapter"]
        recurrent = is_recurrent(e, args)
        note = []
        if recurrent:
            note.append("recurrent across %d samples and %d references"
                        % (len(e["samples"]), len(e["refs"])))
        if adapter:
            note.append("matches %s" % adapter)
        if stats and stats["low_complexity"]:
            note.append("low complexity")
        if not note:
            note.append("no artefact signature fired for this key")
        rows.append([
            str(rank), key, str(len(key)), str(e["n_clips"]), str(e["n_reads"]),
            str(len(e["samples"])), str(len(e["refs"])), str(len(e["groups"])),
            ";".join(sorted(e["groups"])) or "NA",
            str(len(e["refs"]) - len(e["chimp_refs"])),
            str(len(e["chimp_refs"])), ";".join(sorted(e["sides"])),
            "1" if recurrent else "0", adapter or "none",
            "NA" if overlap in (None, 0) else str(overlap),
            "NA" if mismatches is None else str(mismatches),
            fnum(median(e["lens"]), 1),
            fnum(stats["entropy3"], 3) if stats else "NA",
            fnum(stats["homopolymer_frac"], 3) if stats else "NA",
            fnum(stats["gc"], 3) if stats else "NA",
            fnum(stats["dust"], 3) if stats else "NA",
            ("1" if stats["low_complexity"] else "0") if stats else "NA",
            e["example_ref"], "; ".join(note),
        ])
    return rows, items, n_recurrent


def write_clip_table(path, rows, args, n_keys, n_recurrent, dropped):
    extra = [
        "%s collapsed soft-clip sequences by clip count, out of %d distinct "
        "keys; %d key(s) are recurrent (>= %d samples AND >= %d references)"
        % (("the top %d" % args.top_clips) if args.top_clips > 0
           else "every one of the", n_keys, n_recurrent,
           args.recurrent_min_samples, args.recurrent_min_refs),
        "ONE CLIP SEQUENCE RECURRING ACROSS MANY UNRELATED SAMPLES AND "
        "REFERENCES IS A PROBE / ADAPTER / CONTAMINANT SIGNATURE. A biological "
        "anellovirus infection does not put the identical clipped sequence in "
        "unrelated patients on unrelated references",
        "clip_key is the clip truncated at the alignment boundary, so its "
        "sequence is short library sequence, not a genotype and not an "
        "identifier",
        "adapter columns: the best match over the overlapping length; a single "
        "short hit is weak (about 1 read in 1000 matches something by chance "
        "at %d bp with %d mismatch), recurrence is the evidence"
        % (args.min_adapter_overlap, args.adapter_mismatch),
        "this table cannot tell an adapter from a capture probe or from a "
        "conserved genomic motif: it reports the sequence and its spread, and "
        "the sequence must be BLASTed elsewhere before it is named",
    ]
    if dropped:
        extra.append("WARNING: %d clip(s) were not indexed because "
                     "--max-clip-keys %d was reached; counts below are a lower "
                     "bound" % (dropped, args.max_clip_keys))
    write_tsv(path, common_comments(args, extra), CLIP_COLUMNS, rows)


def open_read_table(path, args):
    """Open the per-read table and write its comment block and header."""
    fh = open(path, "w", encoding="ascii", errors="replace", newline="")
    for line in common_comments(args, [
            "one row per read aligned to an anellovirus reference (--emit-reads)",
            "read_key is a truncated SHA1 of QNAME: the read name itself is "
            "never written because it can embed a sample, run or flow-cell "
            "identifier",
            "left/right_clip_key are the boundary-anchored clip keys used in "
            "the clip table; aln_* describe the aligned segment and clip_* the "
            "concatenated clips >= %d bp of the SAME read" % args.min_clip,
            "this table is evidence for inspection, not a result: no row here "
            "is a conclusion about one read",
    ]):
        fh.write("# " + to_ascii(line) + "\n")
    fh.write("\t".join(READ_COLUMNS) + "\n")
    return fh


def write_headers_only(paths, args):
    note = ["no (sample, anellovirus reference) pair could be examined"]
    write_tsv(paths["pair"], common_comments(args, note), PAIR_COLUMNS, [])
    write_tsv(paths["group"], common_comments(args, note), GROUP_COLUMNS, [])
    write_tsv(paths["clips"], common_comments(args, note), CLIP_COLUMNS, [])


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Read-level forensics on the anellovirus signal: clip "
                    "origin (adapter / probe), complexity of the aligned vs "
                    "clipped part, mate origin, alignment ambiguity and UTR "
                    "overlap, against a chimpanzee-reference negative control. "
                    "Standard library only, no figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--run", action="append", default=None, metavar="DIR",
                   help="run directory with bam/ and results/ (repeatable). "
                        "A bare name is resolved under --runs-root. "
                        "Default: the WGS panel run.")
    p.add_argument("--runs-root", default=RUNS_ROOT,
                   help="root holding the run directories (env RUNS_ROOT / "
                        "SHUYU_RUNS_ROOT is honoured)")
    p.add_argument("--refmap", default=DEFAULT_REFMAP,
                   help="panel reference map CSV (reference labels, categories "
                        "and anellovirus keywords)")
    p.add_argument("--indir", default=DEFAULT_INDIR,
                   help="directory holding a10's outputs, used for context")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR,
                   help="output directory (created if absent)")
    p.add_argument("--prefix", default=DEFAULT_PREFIX,
                   help="filename prefix for the output files")
    p.add_argument("--a10-pairs", default=None,
                   help="explicit path to anello_read_audit_by_pair.tsv")
    p.add_argument("--a10-key", default=None,
                   help="explicit path to anello_read_audit_sample_key.tsv")
    p.add_argument("--anello-accessions", default=DEFAULT_ANELLO_ACC_FILE,
                   help="file of anellovirus accessions / reference ids")
    p.add_argument("--chimp-accessions", default=DEFAULT_CHIMP_ACC,
                   help="comma-separated chimpanzee-isolate accessions; these "
                        "are the negative control")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--samtools-threads", type=int, default=1,
                   help="samtools -@ value")
    p.add_argument("--bam-glob", default="bam/*.bam",
                   help="glob for BAMs relative to a run directory")
    p.add_argument("--exclude-flags", default="0x904",
                   help="samtools view -F value (unmapped, secondary, "
                        "supplementary); duplicates are kept on purpose")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="samtools -q value; 0 on purpose, multi-mapping reads "
                        "are the artefact this module characterises")
    p.add_argument("--max-ref-len", type=int, default=100000,
                   help="skip any 'anellovirus' reference longer than this; "
                        "anellovirus genomes are ~3.7 kb, so a longer hit is a "
                        "misclassified reference")
    p.add_argument("--min-clip", type=int, default=10,
                   help="soft-clip length in bp counted as a clip")
    p.add_argument("--adapter-mismatch", type=int, default=1,
                   help="mismatches allowed over the adapter overlap")
    p.add_argument("--min-adapter-overlap", type=int, default=10,
                   help="bases that must overlap before an adapter is called")
    p.add_argument("--clip-key-len", type=int, default=20,
                   help="bases kept from the alignment boundary when clip "
                        "sequences are collapsed (0 = the whole clip)")
    p.add_argument("--top-clips", type=int, default=30,
                   help="clip sequences written to the clip table (0 = all)")
    p.add_argument("--recurrent-min-samples", type=int, default=2,
                   help="distinct samples a clip key needs to count as "
                        "recurrent")
    p.add_argument("--recurrent-min-refs", type=int, default=2,
                   help="distinct references a clip key needs to count as "
                        "recurrent")
    p.add_argument("--max-clip-keys", type=int, default=500000,
                   help="cap on distinct clip keys held in memory")
    p.add_argument("--min-entropy", type=float, default=1.2,
                   help="3-mer entropy in bits below which a segment is called "
                        "low complexity")
    p.add_argument("--max-homopolymer-frac", type=float, default=0.5,
                   help="longest-homopolymer fraction above which a segment is "
                        "called low complexity")
    p.add_argument("--utr-lo", type=float, default=0.75,
                   help="start of the conserved-UTR window, relative position")
    p.add_argument("--utr-hi", type=float, default=1.00,
                   help="end of the conserved-UTR window, relative position")
    p.add_argument("--flag-frac", type=float, default=0.50,
                   help="fraction of a pair's reads at which a per-pair flag "
                        "fires")
    p.add_argument("--max-detail-reads", type=int, default=500000,
                   help="per (sample, reference) cap on reads whose per-read "
                        "metric values are held for the medians; counts are "
                        "never capped, and detail_capped marks a pair that hit "
                        "it")
    p.add_argument("--emit-reads", action="store_true",
                   help="also write the per-read table (large)")
    p.add_argument("--include-chimp", action="store_true",
                   help="fold the chimpanzee references into the human "
                        "aggregates instead of keeping them as the control")
    p.add_argument("--test-groups", default="HIV,HL",
                   help="the two group labels compared")
    p.add_argument("--limit", type=int, default=0,
                   help="examine only the first N BAMs per run (debugging)")
    p.add_argument("--no-run-name-group", action="store_true",
                   help="do not fall back to the run directory name for TCL")
    return p


def sanitise(args):
    if args.min_clip < 1:
        args.min_clip = 1
    if args.adapter_mismatch < 0:
        args.adapter_mismatch = 0
    if args.min_adapter_overlap < 4:
        args.min_adapter_overlap = 4
    if args.clip_key_len < 0:
        args.clip_key_len = 0
    if args.top_clips < 0:
        args.top_clips = 0
    if args.max_clip_keys < 1:
        args.max_clip_keys = 1
    if args.max_detail_reads < 1:
        args.max_detail_reads = 1
    if args.recurrent_min_samples < 1:
        args.recurrent_min_samples = 1
    if args.recurrent_min_refs < 1:
        args.recurrent_min_refs = 1
    if args.utr_lo < 0.0:
        args.utr_lo = 0.0
    if args.utr_hi > 1.0:
        args.utr_hi = 1.0
    if args.utr_hi < args.utr_lo:
        print("WARN: --utr-hi %.2f is below --utr-lo %.2f; using 0.75-1.00"
              % (args.utr_hi, args.utr_lo))
        args.utr_lo, args.utr_hi = 0.75, 1.00
    return args


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = sanitise(build_parser().parse_args(argv))

    test_groups = [g.strip() for g in args.test_groups.split(",") if g.strip()]
    if len(test_groups) != 2:
        print("WARN: --test-groups needs exactly two labels, got %r; using "
              "HIV,HL" % args.test_groups)
        test_groups = ["HIV", "HL"]
    g1, g2 = test_groups

    runs = args.run if args.run else list(DEFAULT_RUNS)
    runs = [r if os.path.isabs(r) else os.path.join(args.runs_root, r)
            for r in runs]

    print("")
    print("%s  %s" % (SCRIPT, TODAY))
    print("runs     %s" % ", ".join(os.path.basename(r.rstrip("/\\")) or r
                                    for r in runs))
    print("filter   samtools view -F %s -q %d (MAPQ not filtered by default)"
          % (args.exclude_flags, args.min_mapq))
    print("clips    >= %d bp | adapter mismatch %d over >= %d bp | clip key %s"
          % (args.min_clip, args.adapter_mismatch, args.min_adapter_overlap,
             ("%d bp" % args.clip_key_len) if args.clip_key_len else "full"))
    print("complex  entropy3 < %.2f bits or homopolymer > %.2f | UTR window "
          "%.2f-%.2f" % (args.min_entropy, args.max_homopolymer_frac,
                         args.utr_lo, args.utr_hi))
    print("")

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            warn_missing("writable output directory", args.outdir)
            return 0
    paths = {
        "pair": os.path.join(args.outdir, out_name(args.prefix, "by_pair.tsv")),
        "group": os.path.join(args.outdir, out_name(args.prefix,
                                                    "by_group.tsv")),
        "clips": os.path.join(args.outdir, clip_out_name(args.prefix)),
        "reads": os.path.join(args.outdir, out_name(args.prefix,
                                                    "by_read.tsv")),
        "key": os.path.join(args.outdir, out_name(args.prefix,
                                                  "sample_key.tsv")),
    }

    anello_raw = read_accession_list(args.anello_accessions)
    if anello_raw is None:
        warn_missing("anellovirus accession list", args.anello_accessions)
        anello_raw = []
    anello_norm = set(norm_acc(a) for a in anello_raw if norm_acc(a))
    anello_ids = set(a.strip() for a in anello_raw)
    chimp_raw = [a.strip() for a in args.chimp_accessions.split(",") if a.strip()]
    chimp_norm = set(norm_acc(a) for a in chimp_raw if norm_acc(a))
    chimp_ids = set(chimp_raw)

    refmap = load_refmap(args.refmap)
    classifier = RefClassifier(refmap, anello_norm, anello_ids, chimp_norm,
                               chimp_ids)
    adapters = build_adapters()

    # ---- pass 1: discover BAMs, anonymise ---------------------------------- #
    run_bams = []
    for run_dir in runs:
        run_dir = run_dir.rstrip("/\\")
        base = os.path.basename(run_dir) or run_dir
        if not os.path.isdir(run_dir):
            warn_missing("run directory", run_dir)
            continue
        if not os.path.isdir(os.path.join(run_dir, "results")):
            warn_missing("results directory for run " + base,
                         os.path.join(run_dir, "results"))
        bams = list_bams(run_dir, args.bam_glob)
        if not bams:
            warn_missing("BAM files for run " + base,
                         os.path.join(run_dir, args.bam_glob))
            continue
        if args.limit > 0:
            bams = bams[:args.limit]
        run_bams.append((run_dir, base, bams))

    real_names = set()
    for _run_dir, _base, bams in run_bams:
        for b in bams:
            real_names.add(sample_of_bam(b))
    anon_of = anonymise(sorted(real_names))

    a10_key, a10_pairs = load_a10_context(args)

    key_rows = []
    runs_of = {}
    group_of_sample = {}
    for _run_dir, base, bams in run_bams:
        for b in bams:
            s = sample_of_bam(b)
            runs_of.setdefault(s, [])
            if base not in runs_of[s]:
                runs_of[s].append(base)
            if group_of_sample.get(s) in (None, "NA"):
                group_of_sample[s] = group_of(s, base, not args.no_run_name_group)
    for name in sorted(real_names):
        key_rows.append([anon_of[name], name, group_of_sample.get(name, "NA"),
                         ";".join(runs_of.get(name, [])),
                         a10_key.get(name, "NA")])
    # Only write the key when there is an identifier to record: an empty file
    # headed "CONTAINS IDENTIFIERS" is misleading in a directory listing.
    if key_rows:
        write_sample_key(paths["key"], key_rows)
        ordered = sorted(real_names)
        print("samples  %d unique real names -> %s..%s"
              % (len(ordered), anon_of[ordered[0]], anon_of[ordered[-1]]))
        print("key      %s (identifiers; do not commit or email)" % paths["key"])
    else:
        print("key      not written: no sample name was read")

    if not run_bams:
        print("No usable run directory; writing headed, empty tables.")
        write_headers_only(paths, args)
        return 0

    ver = check_samtools(args.samtools)
    if ver is None:
        print("Nothing can be computed without samtools; writing headed, "
              "empty tables.")
        write_headers_only(paths, args)
        return 0
    print("samtools %s" % ver)
    print("adapters %d entries (%d stems plus reverse complements and "
          "homopolymers)" % (len(adapters), len(BASE_ADAPTERS)))
    print("")

    # ---- pass 2: read forensics -------------------------------------------- #
    clip_index = {"keys": {}, "dropped": 0}
    read_fh = open_read_table(paths["reads"], args) if args.emit_reads else None
    pairs = []
    oversize = set()
    n_no_index = 0
    n_failed = 0
    n_no_anello_ref = 0
    try:
        for run_dir, base, bams in run_bams:
            done = 0
            for bam_path in bams:
                sample = sample_of_bam(bam_path)
                anon = anon_of[sample]
                group = group_of_sample.get(sample, "NA")
                if not bam_is_indexed(bam_path):
                    n_no_index += 1
                    if n_no_index == 1:
                        warn_missing("BAM index for " + anon,
                                     os.path.join(run_dir, "bam",
                                                  "<sample>.bam.bai"))
                    continue
                lens = idxstats_lengths(run_dir, sample)
                from_idxstats = bool(lens)
                if not lens:
                    lens = header_lengths(args.samtools, bam_path)
                if not lens:
                    warn_missing("idxstats and BAM header for " + anon,
                                 os.path.join(run_dir, "results",
                                              "<sample>.idxstats.tsv"))
                    continue
                anello_refs = {}
                for rname, (seqlen, mapped) in lens.items():
                    if seqlen is None or seqlen < 1:
                        continue
                    if not is_anello_ref(rname, refmap, anello_norm, anello_ids):
                        continue
                    if seqlen > args.max_ref_len:
                        if rname not in oversize:
                            oversize.add(rname)
                            print("WARN: reference %s is %d bp, longer than "
                                  "--max-ref-len %d; it cannot be an "
                                  "anellovirus and is skipped"
                                  % (rname, seqlen, args.max_ref_len))
                        continue
                    if from_idxstats and not mapped:
                        continue
                    anello_refs[rname] = (seqlen, mapped)
                if not anello_refs:
                    n_no_anello_ref += 1
                    continue
                ref_lens = dict((r, v[0]) for r, v in anello_refs.items())
                chimp_of = dict(
                    (r, is_chimp_ref(r, refmap, chimp_norm, chimp_ids))
                    for r in anello_refs)
                ctx = {"run": base, "sample": sample, "sample_anon": anon,
                       "group": group}
                accs, err = stream_reference_set(
                    args, bam_path, ref_lens, chimp_of, classifier, adapters,
                    clip_index, read_fh, ctx)
                if accs is None:
                    n_failed += 1
                    if n_failed == 1:
                        # err was already masked inside stream_sample
                        print("WARN: samtools view failed for %s (%s); that "
                              "sample is skipped" % (anon, err))
                    continue
                done += 1
                a10_anon = a10_key.get(sample, "NA")
                for rname, (seqlen, mapped) in sorted(anello_refs.items()):
                    acc = accs.get(rname)
                    if acc is None or acc["n_reads"] == 0:
                        continue
                    ctx10 = a10_pairs.get((a10_anon, rname), {})
                    pairs.append({
                        "run": base,
                        "sample_anon": anon,
                        "group": group,
                        "reference_id": rname,
                        "ref_label": ref_label(rname, refmap),
                        "chimp": chimp_of.get(rname, False),
                        "idxstats_mapped": mapped,
                        "acc": acc,
                        "a10_anon": a10_anon,
                        "a10_verdict": ctx10.get("verdict", "NA"),
                        "a10_mwf": ctx10.get("mwf", "NA"),
                    })
            print("  %-46s %3d/%3d samples examined" % (base[:46], done,
                                                        len(bams)))
    finally:
        if read_fh is not None:
            read_fh.close()

    if n_no_index:
        print("WARN: %d BAM(s) without an index were skipped" % n_no_index)
    if n_no_anello_ref:
        print("NOTE: %d sample(s) had no anellovirus reference with a mapped "
              "read" % n_no_anello_ref)
    if clip_index["dropped"]:
        print("WARN: --max-clip-keys %d reached, %d clip(s) were not indexed; "
              "clip counts are a lower bound"
              % (args.max_clip_keys, clip_index["dropped"]))

    if not pairs:
        print("")
        print("No (sample, anellovirus reference) pair carried a read; writing "
              "headed, empty tables.")
        write_headers_only(paths, args)
        print("wrote: %s, %s, %s"
              % (paths["pair"], paths["group"], paths["clips"]))
        return 0

    # ---- derive, flag, write ------------------------------------------------ #
    for p in pairs:
        p["m"] = derive(p["acc"])
        n_rec = 0
        n_clips_rec = 0
        for key, count in p["acc"]["clip_keys"].items():
            entry = clip_index["keys"].get(key)
            if entry and is_recurrent(entry, args):
                n_rec += 1
                n_clips_rec += count
        p["n_recurrent_keys"] = n_rec
        p["n_clips_recurrent"] = n_clips_rec
        p["flags"] = pair_flags(p["m"], args, n_rec)

    human_pairs = [p for p in pairs if args.include_chimp or not p["chimp"]]
    chimp_pairs = [p for p in pairs if p["chimp"] and not args.include_chimp]
    pairs_by_set = {"human_anello": human_pairs}
    if chimp_pairs:
        pairs_by_set["chimp_flagged"] = chimp_pairs

    write_by_pair(paths["pair"], pairs, args)
    clip_rows, clip_items, n_recurrent = build_clip_rows(clip_index, args)
    write_clip_table(paths["clips"], clip_rows, args, len(clip_index["keys"]),
                     n_recurrent, clip_index["dropped"])
    group_rows, summaries, tests, per_sample = build_group_rows(
        pairs_by_set, clip_index, args, g1, g2)
    write_by_group(paths["group"], group_rows, args, g1, g2)

    # ---- stdout summary ----------------------------------------------------- #
    print_summary(pairs, human_pairs, chimp_pairs, clip_items,
                  n_recurrent, summaries, tests, args, g1, g2)

    print("")
    print("wrote:")
    for key in ("pair", "group", "clips"):
        print("  %s" % paths[key])
    if args.emit_reads:
        print("  %s" % paths["reads"])
    if key_rows:
        print("  %s" % paths["key"])
        print("REMINDER: %s contains real sample identifiers - do not commit "
              "or email it." % os.path.basename(paths["key"]))
    return 0


def print_summary(pairs, human_pairs, chimp_pairs, clip_items,
                  n_recurrent, summaries, tests, args, g1, g2):
    print("")
    print("-- pairs examined --")
    print("   %d (sample, reference) pair(s), %d human-anellovirus + %d "
          "chimpanzee-flagged control; %d read(s) in total"
          % (len(pairs), len(human_pairs), len(chimp_pairs),
             sum(p["m"]["n_reads"] for p in pairs)))
    flagged = {}
    for p in pairs:
        if p["flags"] == "none":
            continue
        for flag in p["flags"].split(";"):
            flagged[flag] = flagged.get(flag, 0) + 1
    if flagged:
        print("   pair flags: " + ", ".join(
            "%s=%d" % (k, flagged[k]) for k in sorted(flagged)))
    else:
        print("   pair flags: none fired at --flag-frac %.2f" % args.flag_frac)

    print("")
    print("-- top recurrent soft-clip sequences (the probe / adapter test) --")
    shown = 0
    for key, e in clip_items:
        if shown >= 8:
            break
        shown += 1
        print("   %-24s %5d clip(s) in %2d sample(s), %2d reference(s)  %s"
              % (key[:24], e["n_clips"], len(e["samples"]), len(e["refs"]),
                 (e["adapter"][0] or "no adapter match")))
    if not clip_items:
        print("   none: no clip of >= %d bp was seen" % args.min_clip)
    else:
        print("   %d of %d distinct clip key(s) recur in >= %d samples AND >= "
              "%d references." % (n_recurrent, len(clip_items),
                                  args.recurrent_min_samples,
                                  args.recurrent_min_refs))
        print("   A sequence in that list is library / probe / contaminant "
              "sequence, not patient virus.")

    print("")
    print("-- per group and reference set (pooled over reads, descriptive) --")
    print("   %-14s %-5s %7s %7s %9s %8s %8s %8s %8s %8s"
          % ("ref_set", "group", "samples", "pairs", "reads", "clip_ad",
             "clip_lc", "mate_hu", "AS=XS", "in_utr"))
    for ref_set in ("human_anello", "chimp_flagged"):
        for label in GROUP_ORDER + ["ALL"]:
            row = summaries.get((ref_set, label))
            if not row:
                continue
            print("   %-14s %-5s %7s %7s %9s %8s %8s %8s %8s %8s"
                  % (ref_set[:14], label, row["n_samples"], row["n_pairs"],
                     row["n_reads"], row["frac_clips_adapter"],
                     row["clip_frac_low_complexity"], row["frac_mate_HUMAN"],
                     row["frac_as_eq_xs"], row["frac_utr_reads"]))
    print("   clip_ad = clips matching an adapter, clip_lc = clipped segments "
          "of low complexity,")
    print("   mate_hu = reads whose mate is in the human genome, AS=XS = "
          "equally good alignment elsewhere.")

    print("")
    print("-- human references vs the CHIMPANZEE negative control --")
    print("   (paired Wilcoxon over samples carrying both; a LARGE p means the "
          "human references")
    print("    are behaving like the control, which is the artefact reading)")
    print("   %-28s %10s %10s %7s %10s" % ("metric", "human", "chimp",
                                           "n_pairs", "p"))
    any_contrast = False
    for kind, name, res, n1, _n2, extra in tests:
        if kind != "contrast":
            continue
        any_contrast = True
        if res is None:
            # no non-zero difference: the two reference sets were identical on
            # this metric, which is a result, not missing data
            med1, med2 = extra if extra else (None, None)
            print("   %-28s %10s %10s %7d %10s"
                  % (name[:28], fnum(med1, 3), fnum(med2, 3), n1,
                     "identical" if n1 else "no data"))
            continue
        print("   %-28s %10s %10s %7d %10s"
              % (name[:28], fnum(res["median1"], 3), fnum(res["median2"], 3),
                 res["n_pairs"], fp(res["p"])))
    if not any_contrast:
        if args.include_chimp:
            print("   not run: --include-chimp folded the chimpanzee "
                  "references into the human aggregates, so there is no "
                  "held-out negative control left")
        else:
            print("   not run: no sample carried reads on both the human and "
                  "the chimpanzee reference set")

    print("")
    print("-- %s vs %s on the human references (n beside every p) --"
          % (g1, g2))
    for kind, name, res, n1, n2, _extra in tests:
        if kind == "contrast":
            continue
        if res is None:
            print("   %-30s not tested (n %s=%d, %s=%d)"
                  % (name[:30], g1, n1, g2, n2))
            continue
        if kind == "group":
            print("   %-30s U=%s  n %s=%d %s=%d  p=%s  (medians %s vs %s)"
                  % (name[:30], fnum(res["U1"], 1), g1, n1, g2, n2,
                     fp(res["p"]), fnum(res["median1"], 3),
                     fnum(res["median2"], 3)))
        else:
            odds = ("inf" if res["odds_ratio"] == float("inf")
                    else fnum(res["odds_ratio"], 3))
            print("   %-30s OR=%s  n %s=%d %s=%d  p=%s  (%d/%d vs %d/%d)"
                  % (name[:30], odds, g1, n1, g2, n2, fp(res["p"]),
                     res["a"], res["n1"], res["c"], res["n2"]))

    print("")
    print("READ THIS BEFORE QUOTING ANY OF IT:")
    print("   Nothing here is re-aligned and no database is consulted. A read "
          "that fails every")
    print("   test above is NOT thereby real virus - it only means none of the "
          "five tested")
    print("   mechanisms fired. A clip sequence is named only by its match to "
          "the embedded")
    print("   adapter table; BLAST it before calling it a capture probe. And "
          "no p value here")
    print("   is corrected for multiple testing: read the profile, not any one "
          "row.")


if __name__ == "__main__":
    sys.exit(main())
