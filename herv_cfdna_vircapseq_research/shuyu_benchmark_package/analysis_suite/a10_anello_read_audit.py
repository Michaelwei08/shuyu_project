#!/usr/bin/env python3
"""a10_anello_read_audit.py -- is the low-count anellovirus signal real virus or
a mapping artefact?

THE QUESTION
  a7 reports an anellovirus (Torque teno virus) burden that is higher in HIV+
  than in HL within the WGS cohort. The HIV+ side is supported by samples with
  thousands of reads spread over many references; the HL side is supported by a
  handful of samples with 2-144 reads. A handful of reads is exactly what
  cross-mapping produces. Anellovirus genomes are only ~3.7 kb, their ORF1 is
  hypervariable but their UTR is highly conserved ACROSS anellovirus species, so
  the conserved UTR is the obvious cross-mapping magnet: reads from one
  anellovirus (or from an unrelated GC-rich fragment) land on every reference at
  the same conserved spot. This module tests that specific hypothesis at read
  level and refuses to call it either way when the read count is too small.

WHAT IT COMPUTES
  (1) PER (SAMPLE, REFERENCE), for every anellovirus reference with any mapped
      reads in the sample's idxstats. The BAM is streamed once per sample with
      samtools view -F <--exclude-flags> -q <--min-mapq> <bam> <ref> ... and the
      records are dispatched by RNAME, so one samtools call covers every
      anellovirus reference of that sample.
        n_reads, n_distinct_start_positions, n_distinct_pos_cigar
        duplicate_position_fraction  fraction of reads sharing an identical
                                     POS+CIGAR with at least one other read.
                                     Identical starts mean PCR duplicates or a
                                     single locus copied over, not a genome.
        breadth_1x                   fraction of the reference covered >= 1x,
                                     rebuilt from POS + CIGAR (M/=/X/D cover,
                                     N skips, I/S/H/P do not consume reference)
        position_span_frac           (max POS - min POS) / reference length
        max_window_fraction          THE pile-up detector: the fraction of all
                                     reads whose start falls inside the single
                                     most-occupied --window bp window. The
                                     window slides over the observed starts, it
                                     is not a fixed grid, so a pile-up straddling
                                     a bin edge is still caught.
        mean_mapq, median_mapq, frac_mapq0, frac_reads_softclip (a leading or
        trailing S of >= --min-clip bp), frac_low_complexity (SEQ whose two
        commonest bases cover >= --lowcomp-frac of the read: a crude flag, not a
        dust/entropy filter)
      MAPQ is NOT filtered by default (--min-mapq 0) on purpose: multi-mapping
      reads are the artefact, and filtering them out first would hide it.

  (2) THE DECISIVE CROSS-SAMPLE TEST. Read start positions are pooled per virus
      GROUP (HIV, HL) across every sample and every anellovirus reference, each
      start normalised to a relative position along its own reference (0-1), and
      binned into --pooled-bins (default 20) bins per group. Reported per group:
      the whole 20-bin distribution, the fraction of pooled reads in the single
      hottest bin, and how many distinct samples and references feed each bin. If
      HL reads collapse into one bin while HIV reads are spread across the
      genome, the HL signal is positional, i.e. artefactual.

  (3) SHARED HOTSPOTS PER REFERENCE. For each reference, the relative position of
      the hottest window is compared across samples; the module reports how many
      independent samples put their hotspot within --hotspot-tol (relative units)
      of the same place. Independent samples agreeing on one coordinate is the
      strongest artefact evidence available without re-alignment, and a hotspot
      sitting in the conserved terminal UTR is the expected signature. The
      coordinate is reported, never assumed - check it against the reference
      annotation before quoting it as "the UTR".

  (4) VERDICTS. Per (sample, reference), in this precedence:
        too_few_reads   n_reads < --min-reads (default 5)
        pileup_like     max_window_fraction >= --pileup-frac (default 0.5)
        duplicate_like  duplicate_position_fraction >= --dup-frac (default 0.5)
        real_like       breadth_1x >= --min-breadth AND max_window_fraction below
                        --pileup-frac: dispersed over the genome
        indeterminate   anything left over (e.g. dispersed but very low breadth).
                        It exists so nothing is forced into a verdict it does not
                        earn.
      A separate flags column lists every criterion that fired, so a verdict
      taken by precedence never hides the other evidence. Verdicts are then
      aggregated into one row per group so the HIV vs HL comparison is one line,
      with a standard-library Mann-Whitney U (per-sample metrics) and a
      standard-library Fisher exact 2x2 (samples with any real_like / any
      pileup_like). Every p value is printed with the n of each group beside it.

  Chimpanzee-isolate anellovirus references (--chimp-accessions, the same three
  a7 flags) are audited and written out with chimp_flagged=1 but are excluded
  from the pooled distribution, the group aggregates and the tests unless
  --include-chimp is given. They are a built-in negative control: a human sample
  has no chimpanzee TTV, so whatever those references collect is cross-mapping,
  and its positional signature is what a real signal must NOT look like.

WHAT THIS CAN AND CANNOT SUPPORT
  CAN: it can show that reads on a reference are positionally concentrated,
  duplicated, low-MAPQ, soft-clipped or low-complexity, and that independent
  samples share one hotspot. Those are sufficient grounds to call a specific
  (sample, reference) signal an artefact and to remove it from a burden metric.
  CAN: with enough reads, it can show the opposite - dispersed, broad, unique
  coverage that no cross-mapping model explains, i.e. real low-level infection.
  CANNOT: at 2-10 reads NO method distinguishes real virus from artefact. Five
  reads can be dispersed by chance and a real infection can be caught only at its
  conserved region. Those pairs stay too_few_reads and must be reported as
  unresolved, not silently pushed into "artefact" or "real". Since the HL group
  in this cohort carries 2-144 anellovirus reads in total, expect most HL pairs
  to land in too_few_reads: that is the honest answer, and it already means the
  HL side of the HIV vs HL contrast cannot be validated read by read.
  CANNOT: it does not re-align anything, so it cannot say where a cross-mapping
  read really came from, and it cannot rule out a true infection by a divergent
  anellovirus that the panel represents only through its conserved UTR.
  CANNOT: the pooled 20-bin distribution is a description, not a test - reads
  within one sample are not independent observations. Only the per-sample
  Mann-Whitney and the per-sample Fisher rows are tests, and their n is small.

WHAT IT WRITES (tab separated, pure ASCII, into --outdir, --prefix prepended)
  <prefix>_by_pair.tsv            one row per (sample, reference): every metric,
                                  the flags and the verdict
  <prefix>_by_group.tsv           one row per (group, reference set) with the
                                  verdict counts and medians, plus the
                                  Mann-Whitney and Fisher test rows
  <prefix>_pooled_positions.tsv   the pooled relative-position distribution per
                                  group (row_type=group_bin) and the per
                                  reference shared-hotspot summary
                                  (row_type=reference_hotspot)
  <prefix>_sample_key.tsv         real -> anonymous mapping. THE ONLY file that
                                  contains real sample identifiers; its first
                                  line says so. Do not commit or email it.
  With the default prefix: anello_read_audit_by_pair.tsv,
  anello_read_audit_by_group.tsv, anello_read_audit_pooled_positions.tsv,
  anello_read_audit_sample_key.tsv. Samples are anonymised to S001..Snnn by
  sorted real name. Those ids are three digits wide and therefore do NOT equal
  a7's two-digit ids; the a7_sample_anon column (filled from a7's key in
  --indir) is there so the two tables can still be joined.

  Standard library only. No figures, so matplotlib is not imported. No network.
  Any missing input prints "WARN: <what> missing at <path>, skipping" and the
  module still writes headed tables and exits 0.

EXAMPLE
  python3 a10_anello_read_audit.py \
      --run /path/to/runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only \
      --refmap /path/to/runs/shuyu_masked_panel_hg38_herv_line1_refixed/ref/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv \
      --indir /path/to/runs/panel_report_20260725/suite_out \
      --outdir /path/to/runs/panel_report_20260725/suite_out \
      --window 200 --min-reads 5 --pileup-frac 0.5

  # RUNS_ROOT is honoured for the default input paths, as in run_all.sh:
  RUNS_ROOT=/real/run/root python3 a10_anello_read_audit.py --outdir <scratch>

Written 2026-07-26.
"""
from __future__ import annotations

import argparse
import array
import csv
import datetime
import glob
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
DEFAULT_PREFIX = "anello_read_audit"
DEFAULT_ANELLO_ACC_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "anello_accessions.txt")
DEFAULT_CHIMP_ACC = "NC_014069.1,NC_014077.1,NC_014480.2"

# a7 tables looked for in --indir (suffix match, so a prefix is tolerated).
A7_KEY_SUFFIX = "_sample_key.tsv"
A7_KEY_NAME = "a7_virome_sample_key.tsv"
A7_BURDEN_SUFFIX = "anellovirus_burden.tsv"
A7_BURDEN_NAME = "a7_virome_anellovirus_burden.tsv"

ANELLO_KEYWORDS = [
    "anello", "torque teno", "torque-teno", "torquetenovirus",
    "transfusion transmitted virus", "transfusion-transmitted virus",
    "small anellovirus", "tt virus", "ttv", "ttmv", "ttmdv", "sen virus",
]
# Short keywords that must not match inside a longer word.
SHORT_KEYWORDS = set(["ttv", "ttmv", "ttmdv", "tt virus"])

ACC_RE = re.compile(
    r"(?<![A-Za-z0-9])([A-Z]{2}_\d{5,8}|[A-Z]{1,2}\d{5,6})(?:\.\d+)?(?![A-Za-z0-9])")

VERDICTS = ["too_few_reads", "pileup_like", "duplicate_like", "real_like",
            "indeterminate"]
GROUP_ORDER = ["HIV", "HL", "TCL", "NA"]

TODAY = datetime.date.today().isoformat()
SCRIPT = os.path.basename(__file__)

PAIR_COLUMNS = [
    "run", "sample_anon", "group", "reference_id", "ref_label", "ref_len",
    "chimp_flagged", "idxstats_mapped", "n_reads", "n_reads_with_seq",
    "n_flagged_duplicate", "n_distinct_start_positions", "n_distinct_pos_cigar",
    "max_identical_pos_cigar", "duplicate_position_fraction",
    "covered_bases_1x", "breadth_1x", "position_span_bp", "position_span_frac",
    "window_bp", "max_window_reads", "max_window_fraction",
    "max_window_start", "max_window_end", "max_window_rel_pos",
    "mean_mapq", "median_mapq", "frac_mapq0",
    "softclip_min_bp", "frac_reads_softclip", "frac_low_complexity",
    "a7_sample_anon", "a7_anello_reads_total", "a7_anello_richness",
    "flags", "verdict", "verdict_metrics",
]

GROUP_COLUMNS = [
    "row_type", "label", "ref_set",
    "n_samples", "n_pairs", "n_pairs_auditable", "n_reads",
    "n_too_few_reads", "n_pileup_like", "n_duplicate_like", "n_real_like",
    "n_indeterminate", "n_samples_any_real_like", "n_samples_any_pileup_like",
    "n_samples_all_too_few_reads",
    "median_reads_per_pair", "median_breadth", "median_max_window_fraction",
    "median_duplicate_position_fraction", "median_mapq", "median_frac_mapq0",
    "median_frac_softclip",
    "pooled_reads", "pooled_hottest_bin", "pooled_hottest_bin_rel_start",
    "pooled_hottest_bin_fraction",
    "group1", "n1", "value_group1", "group2", "n2", "value_group2",
    "statistic", "statistic_value", "p_value", "effect", "note",
]

POOLED_COLUMNS = [
    "row_type", "scope", "ref_set", "reference_id", "ref_label", "ref_len",
    "n_bins", "bin_index", "rel_start", "rel_end", "n_reads",
    "frac_of_scope_reads", "n_samples", "n_references",
    "hotspot_rel_pos", "n_pairs", "n_pairs_at_hotspot",
    "frac_pairs_at_hotspot", "median_max_window_fraction_at_hotspot",
    "hotspot_tolerance_rel", "note",
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


def redact(text, sample, anon):
    """Replace a real sample name inside a string with its anonymous ID."""
    if not sample:
        return text
    return str(text).replace(sample, anon)


# samtools echoes the full BAM path on failure, which names the controlled-data
# mount. Redacting the sample name alone is not enough; scrub path-like tokens
# too. Same helper as a11/a12.
PATH_TOKEN_RE = re.compile(r"[^\s'\"]*[/\\][^\s'\"]*")


def strip_paths(text):
    """Replace every filesystem-path-looking token with a literal <path>."""
    return PATH_TOKEN_RE.sub("<path>", str(text or ""))


def out_name(prefix, base):
    return (prefix + "_" + base) if prefix else base


def write_tsv(path, comments, header, rows):
    with open(path, "w", encoding="ascii", errors="replace", newline="") as fh:
        for line in comments:
            fh.write("# " + to_ascii(line) + "\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(to_ascii(c) for c in row) + "\n")


# --------------------------------------------------------------------------- #
# sample naming (identical rule to a5 / a7)
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


def _log_choose(n, k):
    if k < 0 or k > n:
        return None
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def fisher_exact_2x2(a, b, c, d):
    """Two-sided Fisher exact test on [[a, b], [c, d]], standard library only.

    a, b = successes / failures in group 1; c, d = the same in group 2. The
    two-sided p sums the hypergeometric probability of every table at least as
    unlikely as the observed one (the same convention scipy calls 'two-sided').
    Returns None when the table is degenerate (an empty row or column).
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
# reference identification
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
# a7 context (optional)
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


def load_a7_context(args):
    """(real_sample -> a7 anon, a7 anon -> {reads, richness}). Never fatal."""
    key_map = {}
    burden = {}
    key_path = args.a7_key or resolve_in_indir(args.indir, A7_KEY_NAME,
                                               A7_KEY_SUFFIX)
    if not key_path or not os.path.exists(key_path):
        warn_missing("a7 sample key", key_path or os.path.join(
            args.indir or "<indir>", A7_KEY_NAME))
    else:
        rows = read_commented_tsv(key_path)
        if rows is None:
            warn_missing("readable a7 sample key", key_path)
        else:
            for row in rows:
                real = (row.get("real_sample") or "").strip()
                anon = (row.get("anon_sample") or "").strip()
                if real and anon:
                    key_map[real] = anon

    burden_path = args.a7_burden or resolve_in_indir(args.indir, A7_BURDEN_NAME,
                                                     A7_BURDEN_SUFFIX)
    if not burden_path or not os.path.exists(burden_path):
        warn_missing("a7 anellovirus burden table", burden_path or os.path.join(
            args.indir or "<indir>", A7_BURDEN_NAME))
    else:
        rows = read_commented_tsv(burden_path)
        if rows is None:
            warn_missing("readable a7 anellovirus burden table", burden_path)
        else:
            for row in rows:
                anon = (row.get("sample") or "").strip()
                if not anon:
                    continue
                burden[anon] = {
                    "reads": (row.get("anello_reads_human_total") or "NA").strip(),
                    "richness": (row.get("anello_richness_human") or "NA").strip(),
                }
    return key_map, burden


# --------------------------------------------------------------------------- #
# read-level audit
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


def is_low_complexity(seq, frac):
    """True when the two commonest bases cover >= frac of the read."""
    n = len(seq)
    if n <= 0:
        return False
    counts = {}
    for ch in seq.upper():
        counts[ch] = counts.get(ch, 0) + 1
    top = sorted(counts.values(), reverse=True)[:2]
    return (float(sum(top)) / float(n)) >= frac


def new_acc(ref_len):
    return {
        "ref_len": ref_len,
        "n_reads": 0,
        "n_dup_flag": 0,
        "n_mapq0": 0,
        "n_softclip": 0,
        "n_seq": 0,
        "n_lowcomp": 0,
        "starts": [],
        "poscigar": {},
        "mapqs": [],
        "delta": array.array("i", [0]) * (ref_len + 2),
    }


def stream_sample(samtools, bam, regions, ref_lens, args, sample, anon,
                  chunk=200):
    """Stream one BAM over every wanted anellovirus reference.

    Returns (accumulators keyed by reference_id, error_string). The
    accumulators are None when samtools could not be run at all. sample / anon
    are used only to mask the real sample name out of samtools' error text
    before it is truncated, so a truncated path cannot leak half a name.
    """
    accs = {}
    ref_list = sorted(regions)
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
            return None, "samtools not runnable (%s)" % redact(exc, sample, anon)
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
            acc["n_reads"] += 1
            if flag & 0x400:
                acc["n_dup_flag"] += 1
            if mapq == 0:
                acc["n_mapq0"] += 1
            acc["mapqs"].append(mapq)
            acc["starts"].append(pos0)
            cigar = f[5]
            key = f[3] + ":" + cigar
            acc["poscigar"][key] = acc["poscigar"].get(key, 0) + 1
            seq = f[9]
            if seq and seq != "*":
                acc["n_seq"] += 1
                if is_low_complexity(seq, args.lowcomp_frac):
                    acc["n_lowcomp"] += 1
            ops = parse_cigar(cigar)
            if not ops:
                continue
            if ((ops[0][1] == "S" and ops[0][0] >= args.min_clip)
                    or (ops[-1][1] == "S" and ops[-1][0] >= args.min_clip)):
                acc["n_softclip"] += 1
            pos = pos0
            delta = acc["delta"]
            for num, ch in ops:
                if ch in "M=XD":
                    s = pos if pos > 0 else 0
                    e = pos + num
                    if e > ref_len:
                        e = ref_len
                    if e > s:
                        delta[s] += 1
                        delta[e] -= 1
                    pos += num
                elif ch == "N":
                    pos += num
                # I, S, H, P consume query or nothing on the reference
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
            # into a form a later replace() would no longer match.
            err = strip_paths(redact(err, sample, anon))
            err = " ".join(to_ascii(err.replace("\r", " ").replace("\n", " ")).split())
            return None, "samtools view rc=%d %s" % (rc, err[:160] or "no stderr")
    return accs, ""


def max_window(starts_sorted, window):
    """(reads in the fullest sliding window, window start, window end index).

    The window slides over the observed start positions, so a pile-up that
    straddles a fixed bin edge is still counted in one piece.
    """
    n = len(starts_sorted)
    if n == 0:
        return 0, None, None
    best, best_i, best_j = 0, 0, 0
    j = 0
    for i in range(n):
        if j < i:
            j = i
        while j + 1 < n and starts_sorted[j + 1] < starts_sorted[i] + window:
            j += 1
        count = j - i + 1
        if count > best:
            best, best_i, best_j = count, i, j
    return best, best_i, best_j


def pair_metrics(acc, args):
    """Every per-(sample, reference) metric from one accumulator."""
    n = acc["n_reads"]
    ref_len = acc["ref_len"]
    starts = sorted(acc["starts"])
    groups = acc["poscigar"]
    n_dup_reads = sum(c for c in groups.values() if c >= 2)
    max_ident = max(groups.values()) if groups else 0

    covered = 0
    depth = 0
    delta = acc["delta"]
    for i in range(ref_len):
        depth += delta[i]
        if depth:
            covered += 1

    span = (starts[-1] - starts[0]) if starts else 0
    best, bi, bj = max_window(starts, args.window)
    if best > 0:
        win_start = starts[bi] + 1
        win_end = min(starts[bi] + args.window, ref_len)
        rel = ((starts[bi] + starts[bj]) / 2.0) / float(ref_len)
        rel = min(1.0, max(0.0, rel))
    else:
        win_start = win_end = rel = None

    return {
        "ref_len": ref_len,
        "n_reads": n,
        "n_reads_with_seq": acc["n_seq"],
        "n_flagged_duplicate": acc["n_dup_flag"],
        "n_distinct_starts": len(set(starts)),
        "n_distinct_pos_cigar": len(groups),
        "max_identical_pos_cigar": max_ident,
        "dup_frac": (float(n_dup_reads) / float(n)) if n else None,
        "covered_bases": covered,
        "breadth": (float(covered) / float(ref_len)) if ref_len else None,
        "span_bp": span,
        "span_frac": (float(span) / float(ref_len)) if ref_len else None,
        "max_window_reads": best,
        "max_window_frac": (float(best) / float(n)) if n else None,
        "max_window_start": win_start,
        "max_window_end": win_end,
        "max_window_rel": rel,
        "mean_mapq": mean(acc["mapqs"]),
        "median_mapq": median(acc["mapqs"]),
        "frac_mapq0": (float(acc["n_mapq0"]) / float(n)) if n else None,
        "frac_softclip": (float(acc["n_softclip"]) / float(n)) if n else None,
        "frac_lowcomp": ((float(acc["n_lowcomp"]) / float(acc["n_seq"]))
                         if acc["n_seq"] else None),
        "starts": starts,
    }


def verdict_of(m, args):
    """(verdict, flags string, metrics string) for one (sample, reference)."""
    flags = []
    pileup = (m["max_window_frac"] is not None
              and m["max_window_frac"] >= args.pileup_frac)
    duplicated = (m["dup_frac"] is not None and m["dup_frac"] >= args.dup_frac)
    broad = (m["breadth"] is not None and m["breadth"] >= args.min_breadth)
    if pileup:
        flags.append("window_pileup")
    if duplicated:
        flags.append("duplicate_positions")
    if not broad:
        flags.append("low_breadth")
    if m["frac_mapq0"] is not None and m["frac_mapq0"] >= 0.5:
        flags.append("mapq0_majority")
    if m["frac_softclip"] is not None and m["frac_softclip"] >= 0.5:
        flags.append("softclip_majority")
    if m["frac_lowcomp"] is not None and m["frac_lowcomp"] >= 0.5:
        flags.append("low_complexity_majority")
    if (m["n_reads"] >= 2 and m["n_distinct_starts"] == 1):
        flags.append("single_start_position")

    if m["n_reads"] < args.min_reads:
        verdict = "too_few_reads"
    elif pileup:
        verdict = "pileup_like"
    elif duplicated:
        verdict = "duplicate_like"
    elif broad:
        verdict = "real_like"
    else:
        verdict = "indeterminate"

    metrics = ("reads=%d;distinct_starts=%d;dup_frac=%s;breadth=%s;"
               "max_window_frac=%s;span_frac=%s;median_mapq=%s"
               % (m["n_reads"], m["n_distinct_starts"], fnum(m["dup_frac"], 3),
                  fnum(m["breadth"], 4), fnum(m["max_window_frac"], 3),
                  fnum(m["span_frac"], 3), fnum(m["median_mapq"], 1)))
    return verdict, (";".join(flags) if flags else "none"), metrics


# --------------------------------------------------------------------------- #
# aggregation
# --------------------------------------------------------------------------- #
def weighted_mean(pairs, field):
    """Read-weighted mean of a per-pair field over auditable pairs."""
    num = 0.0
    den = 0.0
    for p in pairs:
        val = p["m"][field]
        if val is None:
            continue
        num += float(val) * float(p["m"]["n_reads"])
        den += float(p["m"]["n_reads"])
    return (num / den) if den > 0 else None


def per_sample_metrics(pairs, args):
    """sample_anon -> {group, weighted metrics, verdict flags} over one ref set."""
    by_sample = {}
    for p in pairs:
        by_sample.setdefault(p["sample_anon"], []).append(p)
    out = {}
    for anon, plist in by_sample.items():
        auditable = [p for p in plist if p["m"]["n_reads"] >= args.min_reads]
        out[anon] = {
            "group": plist[0]["group"],
            "n_pairs": len(plist),
            "n_auditable": len(auditable),
            "n_reads": sum(p["m"]["n_reads"] for p in plist),
            "max_window_frac": weighted_mean(auditable, "max_window_frac"),
            "breadth": weighted_mean(auditable, "breadth"),
            "dup_frac": weighted_mean(auditable, "dup_frac"),
            "frac_mapq0": weighted_mean(auditable, "frac_mapq0"),
            "any_real_like": any(p["verdict"] == "real_like" for p in plist),
            "any_pileup_like": any(p["verdict"] == "pileup_like" for p in plist),
            "all_too_few": all(p["verdict"] == "too_few_reads" for p in plist),
        }
    return out


MW_SPECS = [
    ("max_window_fraction_weighted", "max_window_frac"),
    ("breadth_1x_weighted", "breadth"),
    ("duplicate_position_fraction_weighted", "dup_frac"),
    ("frac_mapq0_weighted", "frac_mapq0"),
]


def pooled_bins(pairs, nbins):
    """scope -> {counts, samples per bin, refs per bin, total}."""
    scopes = {}
    for p in pairs:
        ref_len = p["m"]["ref_len"]
        if not ref_len:
            continue
        for scope in ("ALL", p["group"]):
            st = scopes.setdefault(scope, {
                "counts": [0] * nbins,
                "samples": [set() for _ in range(nbins)],
                "refs": [set() for _ in range(nbins)],
                "total": 0,
                "all_samples": set(),
                "all_refs": set(),
            })
            st["all_samples"].add(p["sample_anon"])
            st["all_refs"].add(p["reference_id"])
            for pos in p["m"]["starts"]:
                rel = float(pos) / float(ref_len)
                idx = int(rel * nbins)
                if idx >= nbins:
                    idx = nbins - 1
                if idx < 0:
                    idx = 0
                st["counts"][idx] += 1
                st["samples"][idx].add(p["sample_anon"])
                st["refs"][idx].add(p["reference_id"])
                st["total"] += 1
    return scopes


def hottest_bin(state):
    if not state or state["total"] <= 0:
        return None, None
    counts = state["counts"]
    best = max(range(len(counts)), key=lambda i: counts[i])
    return best, float(counts[best]) / float(state["total"])


def shared_hotspot(positions, tol):
    """(rel position of the densest cluster, n within tol) over relative posns."""
    if not positions:
        return None, 0
    best_pos, best_n = None, 0
    for anchor in sorted(positions):
        near = [q for q in positions if abs(q - anchor) <= tol]
        if len(near) > best_n:
            best_n = len(near)
            best_pos = median(near)
    return best_pos, best_n


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def common_comments(args, extra=()):
    lines = [
        "%s  generated %s" % (SCRIPT, TODAY),
        "no real sample identifiers in this file; ids are anonymous S001..Snnn "
        "(mapping in %s)" % out_name(args.prefix, "sample_key.tsv"),
        "read filter: samtools view -F %s -q %d (MAPQ is deliberately NOT "
        "filtered by default: multi-mapping reads are the artefact)"
        % (args.exclude_flags, args.min_mapq),
        "window %d bp, min-reads %d, pileup-frac %.2f, dup-frac %.2f, "
        "min-breadth %.2f, min-clip %d bp"
        % (args.window, args.min_reads, args.pileup_frac, args.dup_frac,
           args.min_breadth, args.min_clip),
    ]
    if args.min_reads > 1:
        lines.append(
            "at %s reads no method separates real virus from artefact; those "
            "pairs stay too_few_reads and are NOT evidence either way"
            % ("1" if args.min_reads == 2 else "1-%d" % (args.min_reads - 1)))
    else:
        lines.append(
            "--min-reads is 1, so nothing is held back as too_few_reads; at a "
            "few reads no method separates real virus from artefact, so read "
            "every verdict below %d reads as unresolved" % 5)
    lines.extend(extra)
    return lines


def write_sample_key(path, samples):
    """samples: list of (anon, real, group, runs, a7_anon)."""
    lines = [
        "CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
        "generated %s by %s" % (TODAY, SCRIPT),
    ]
    header = ["sample_anon", "sample_real", "group", "runs", "a7_sample_anon"]
    write_tsv(path, lines, header, samples)


def write_by_pair(path, pairs, args):
    rows = []
    for p in sorted(pairs, key=lambda r: (r["group"], r["sample_anon"],
                                          r["reference_id"])):
        m = p["m"]
        rows.append([
            p["run"], p["sample_anon"], p["group"], p["reference_id"],
            p["ref_label"], str(m["ref_len"]), "1" if p["chimp"] else "0",
            "NA" if p["idxstats_mapped"] is None else str(p["idxstats_mapped"]),
            str(m["n_reads"]), str(m["n_reads_with_seq"]),
            str(m["n_flagged_duplicate"]), str(m["n_distinct_starts"]),
            str(m["n_distinct_pos_cigar"]), str(m["max_identical_pos_cigar"]),
            fnum(m["dup_frac"], 4), str(m["covered_bases"]),
            fnum(m["breadth"], 6), str(m["span_bp"]), fnum(m["span_frac"], 4),
            str(args.window), str(m["max_window_reads"]),
            fnum(m["max_window_frac"], 4),
            "NA" if m["max_window_start"] is None else str(m["max_window_start"]),
            "NA" if m["max_window_end"] is None else str(m["max_window_end"]),
            fnum(m["max_window_rel"], 4),
            fnum(m["mean_mapq"], 2), fnum(m["median_mapq"], 1),
            fnum(m["frac_mapq0"], 4), str(args.min_clip),
            fnum(m["frac_softclip"], 4), fnum(m["frac_lowcomp"], 4),
            p["a7_anon"], p["a7_reads"], p["a7_richness"],
            p["flags"], p["verdict"], p["metrics"],
        ])
    write_tsv(path, common_comments(args, [
        "one row per (sample, anellovirus reference) with any mapped read",
        "max_window_fraction is the share of reads whose start falls in the "
        "fullest sliding %d bp window; >= %.2f is called pileup_like"
        % (args.window, args.pileup_frac),
        "chimp_flagged=1 rows are the negative control and are excluded from "
        "the group aggregates unless --include-chimp was given",
        "a7_* columns are context copied from the a7 burden table, joined "
        "through the real name; they are not recomputed here",
    ]), PAIR_COLUMNS, rows)


def blank_group_row():
    return dict((c, "NA") for c in GROUP_COLUMNS)


def build_group_rows(pairs_by_set, args, g1, g2):
    """(rows for by_group.tsv, list of (name, result) for the stdout summary)."""
    rows = []
    summaries = {}
    for ref_set in ("human_anello", "chimp_flagged"):
        pairs = pairs_by_set.get(ref_set, [])
        if not pairs:
            continue
        nbins = args.pooled_bins
        scopes = pooled_bins(pairs, nbins)
        by_group = {}
        for p in pairs:
            by_group.setdefault(p["group"], []).append(p)
        by_group["ALL"] = list(pairs)
        labels = [g for g in GROUP_ORDER if g in by_group]
        labels += [g for g in sorted(by_group) if g not in labels and g != "ALL"]
        labels.append("ALL")
        for label in labels:
            plist = by_group[label]
            auditable = [p for p in plist if p["m"]["n_reads"] >= args.min_reads]
            per_sample = per_sample_metrics(plist, args)
            state = scopes.get(label)
            hb, hbfrac = hottest_bin(state)
            row = blank_group_row()
            row["row_type"] = "group_summary"
            row["label"] = label
            row["ref_set"] = ref_set
            row["n_samples"] = str(len(per_sample))
            row["n_pairs"] = str(len(plist))
            row["n_pairs_auditable"] = str(len(auditable))
            row["n_reads"] = str(sum(p["m"]["n_reads"] for p in plist))
            for verdict in VERDICTS:
                row["n_" + verdict] = str(
                    sum(1 for p in plist if p["verdict"] == verdict))
            row["n_samples_any_real_like"] = str(
                sum(1 for v in per_sample.values() if v["any_real_like"]))
            row["n_samples_any_pileup_like"] = str(
                sum(1 for v in per_sample.values() if v["any_pileup_like"]))
            row["n_samples_all_too_few_reads"] = str(
                sum(1 for v in per_sample.values() if v["all_too_few"]))
            row["median_reads_per_pair"] = fnum(
                median([p["m"]["n_reads"] for p in plist]), 1)
            row["median_breadth"] = fnum(median(
                [p["m"]["breadth"] for p in auditable
                 if p["m"]["breadth"] is not None]), 4)
            row["median_max_window_fraction"] = fnum(median(
                [p["m"]["max_window_frac"] for p in auditable
                 if p["m"]["max_window_frac"] is not None]), 4)
            row["median_duplicate_position_fraction"] = fnum(median(
                [p["m"]["dup_frac"] for p in auditable
                 if p["m"]["dup_frac"] is not None]), 4)
            row["median_mapq"] = fnum(median(
                [p["m"]["median_mapq"] for p in auditable
                 if p["m"]["median_mapq"] is not None]), 1)
            row["median_frac_mapq0"] = fnum(median(
                [p["m"]["frac_mapq0"] for p in auditable
                 if p["m"]["frac_mapq0"] is not None]), 4)
            row["median_frac_softclip"] = fnum(median(
                [p["m"]["frac_softclip"] for p in auditable
                 if p["m"]["frac_softclip"] is not None]), 4)
            row["pooled_reads"] = str(state["total"]) if state else "0"
            row["pooled_hottest_bin"] = "NA" if hb is None else str(hb)
            row["pooled_hottest_bin_rel_start"] = (
                "NA" if hb is None else fnum(float(hb) / float(nbins), 4))
            row["pooled_hottest_bin_fraction"] = fnum(hbfrac, 4)
            row["note"] = ("medians are over the %d pair(s) with >= %d reads"
                           % (len(auditable), args.min_reads))
            rows.append([row[c] for c in GROUP_COLUMNS])
            if ref_set == "human_anello":
                summaries[label] = dict(row)
                summaries[label]["_per_sample"] = per_sample

    # ---- tests, human anellovirus references only ------------------------- #
    tests = []
    pairs = pairs_by_set.get("human_anello", [])
    per_sample = per_sample_metrics(pairs, args) if pairs else {}
    s1 = [v for v in per_sample.values() if v["group"] == g1]
    s2 = [v for v in per_sample.values() if v["group"] == g2]

    for name, field in MW_SPECS:
        x = [v[field] for v in s1 if v[field] is not None]
        y = [v[field] for v in s2 if v[field] is not None]
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
            row["note"] = ("not tested: %d %s sample(s) and %d %s sample(s) had "
                           "an auditable pair (>= %d reads)"
                           % (len(x), g1, len(y), g2, args.min_reads))
        else:
            row["value_group1"] = fnum(res["median1"], 4)
            row["value_group2"] = fnum(res["median2"], 4)
            row["statistic_value"] = fnum(res["U1"], 1)
            row["p_value"] = fp(res["p"])
            row["effect"] = fnum(res["effect_r"], 4)
            note = ("two-sided normal approximation, tie- and "
                    "continuity-corrected; values are group medians of the "
                    "read-weighted per-sample metric")
            if min(res["n1"], res["n2"]) < 5:
                note += "; n<5 in one group, p is approximate and underpowered"
            row["note"] = note
        rows.append([row[c] for c in GROUP_COLUMNS])
        tests.append((name, res, len(x), len(y)))

    for name, field in (("any_real_like_pair", "any_real_like"),
                        ("any_pileup_like_pair", "any_pileup_like")):
        a = sum(1 for v in s1 if v[field])
        b = len(s1) - a
        c = sum(1 for v in s2 if v[field])
        d = len(s2) - c
        res = fisher_exact_2x2(a, b, c, d)
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
            row["note"] = ("not tested: degenerate 2x2 table [[%d,%d],[%d,%d]] "
                           "(an empty row or column)" % (a, b, c, d))
        else:
            row["statistic_value"] = (
                "inf" if res["odds_ratio"] == float("inf")
                else fnum(res["odds_ratio"], 4))
            row["p_value"] = fp(res["p"])
            row["note"] = ("two-sided Fisher exact on [[%d,%d],[%d,%d]]; "
                           "samples counted once each, not pairs" % (a, b, c, d))
        rows.append([row[c] for c in GROUP_COLUMNS])
        tests.append((name, res, len(s1), len(s2)))

    return rows, summaries, tests


def write_by_group(path, rows, args, g1, g2):
    write_tsv(path, common_comments(args, [
        "row_type=group_summary: one row per (group, reference set), so the "
        "%s vs %s comparison is one line each" % (g1, g2),
        "row_type=group_test: the %s vs %s tests; n1/n2 are the group sizes "
        "that entered the test and sit beside every p value" % (g1, g2),
        "Mann-Whitney is over per-sample read-weighted metrics, Fisher is over "
        "samples; neither pools reads, because reads within a sample are not "
        "independent observations",
        "ref_set=chimp_flagged is the negative control (no human sample has "
        "chimpanzee TTV); it is summarised but never tested",
    ]), GROUP_COLUMNS, rows)


def build_pooled_rows(pairs, args):
    """(rows for pooled_positions.tsv, scopes state, hotspot records)."""
    nbins = args.pooled_bins
    scopes = pooled_bins(pairs, nbins)
    rows = []
    order = [s for s in ["ALL"] + GROUP_ORDER if s in scopes]
    order += [s for s in sorted(scopes) if s not in order]
    for scope in order:
        st = scopes[scope]
        for b in range(nbins):
            row = dict((c, "NA") for c in POOLED_COLUMNS)
            row["row_type"] = "group_bin"
            row["scope"] = scope
            row["ref_set"] = "human_anello"
            row["reference_id"] = "ALL_ANELLO_REFS"
            row["n_bins"] = str(nbins)
            row["bin_index"] = str(b)
            row["rel_start"] = fnum(float(b) / float(nbins), 4)
            row["rel_end"] = fnum(float(b + 1) / float(nbins), 4)
            row["n_reads"] = str(st["counts"][b])
            row["frac_of_scope_reads"] = fnum(
                (float(st["counts"][b]) / float(st["total"]))
                if st["total"] else None, 4)
            row["n_samples"] = str(len(st["samples"][b]))
            row["n_references"] = str(len(st["refs"][b]))
            row["note"] = ("pooled read starts, relative position along each "
                           "reference; %d read(s) from %d sample(s) in scope"
                           % (st["total"], len(st["all_samples"])))
            rows.append([row[c] for c in POOLED_COLUMNS])

    # ---- per-reference shared hotspots ------------------------------------ #
    hotspots = []
    by_ref = {}
    for p in pairs:
        if p["m"]["n_reads"] < args.min_reads:
            continue
        if p["m"]["max_window_rel"] is None:
            continue
        by_ref.setdefault(p["reference_id"], []).append(p)
    for rid in sorted(by_ref):
        plist = by_ref[rid]
        for scope in ["ALL"] + GROUP_ORDER:
            sub = (plist if scope == "ALL"
                   else [p for p in plist if p["group"] == scope])
            if len(sub) < 2:
                continue
            positions = [p["m"]["max_window_rel"] for p in sub]
            pos, n_at = shared_hotspot(positions, args.hotspot_tol)
            frac = float(n_at) / float(len(sub))
            # A shared coordinate only means something if the pairs sitting on
            # it are actually concentrated. Dispersed profiles still have a
            # fullest window somewhere, and with even coverage the tie is broken
            # arbitrarily, so carry the concentration of the agreeing pairs.
            at_hot = [p for p in sub
                      if pos is not None
                      and abs(p["m"]["max_window_rel"] - pos) <= args.hotspot_tol]
            med_conc = median([p["m"]["max_window_frac"] for p in at_hot
                               if p["m"]["max_window_frac"] is not None])
            row = dict((c, "NA") for c in POOLED_COLUMNS)
            row["row_type"] = "reference_hotspot"
            row["scope"] = scope
            row["ref_set"] = "human_anello"
            row["reference_id"] = rid
            row["ref_label"] = sub[0]["ref_label"]
            row["ref_len"] = str(sub[0]["m"]["ref_len"])
            row["n_reads"] = str(sum(p["m"]["n_reads"] for p in sub))
            row["n_samples"] = str(len(set(p["sample_anon"] for p in sub)))
            row["hotspot_rel_pos"] = fnum(pos, 4)
            row["n_pairs"] = str(len(sub))
            row["n_pairs_at_hotspot"] = str(n_at)
            row["frac_pairs_at_hotspot"] = fnum(frac, 4)
            row["median_max_window_fraction_at_hotspot"] = fnum(med_conc, 4)
            row["hotspot_tolerance_rel"] = fnum(args.hotspot_tol, 4)
            row["note"] = ("independent samples agreeing on one hotspot is "
                           "artefact evidence, but only when they are "
                           "concentrated: read "
                           "median_max_window_fraction_at_hotspot with it; "
                           "approximate bp = rel * ref_len")
            rows.append([row[c] for c in POOLED_COLUMNS])
            hotspots.append({
                "scope": scope, "reference_id": rid,
                "label": sub[0]["ref_label"],
                "ref_len": sub[0]["m"]["ref_len"], "pos": pos,
                "n_at": n_at, "n_pairs": len(sub), "frac": frac,
                "conc": med_conc,
            })
    return rows, scopes, hotspots


def write_pooled(path, rows, args):
    write_tsv(path, common_comments(args, [
        "row_type=group_bin: pooled read start positions per group, every "
        "start normalised to 0-1 along its own reference, %d bins"
        % args.pooled_bins,
        "row_type=reference_hotspot: per reference, how many independent "
        "(sample, reference) pairs put their fullest %d bp window within "
        "%.3f relative units of the same position" % (args.window,
                                                      args.hotspot_tol),
        "the pooled distribution is DESCRIPTIVE, not a test: reads within one "
        "sample are not independent. The tests are in %s"
        % out_name(args.prefix, "by_group.tsv"),
        "a shared hotspot in the conserved terminal UTR is the expected "
        "cross-mapping signature; the coordinate is reported, not assumed - "
        "check it against the reference annotation before naming it",
    ]), POOLED_COLUMNS, rows)


def write_headers_only(paths, args):
    write_tsv(paths["pair"], common_comments(args, ["no pair could be audited"]),
              PAIR_COLUMNS, [])
    write_tsv(paths["group"], common_comments(args, ["no pair could be audited"]),
              GROUP_COLUMNS, [])
    write_tsv(paths["pooled"], common_comments(args, ["no pair could be audited"]),
              POOLED_COLUMNS, [])


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Read-level audit of the anellovirus signal: pile-up, "
                    "duplicate, MAPQ and shared-hotspot evidence for real "
                    "virus vs cross-mapping. Standard library only, no figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--run", action="append", default=None, metavar="DIR",
                   help="run directory with bam/ and results/ (repeatable). "
                        "A bare name is resolved under --runs-root. "
                        "Default: the WGS panel run.")
    p.add_argument("--runs-root", default=RUNS_ROOT,
                   help="root holding the run directories (env RUNS_ROOT / "
                        "SHUYU_RUNS_ROOT is honoured)")
    p.add_argument("--refmap", default=DEFAULT_REFMAP,
                   help="panel reference map CSV (labels + anellovirus keywords)")
    p.add_argument("--indir", default=DEFAULT_INDIR,
                   help="directory holding the a7 outputs used for context")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR,
                   help="output directory (created if absent)")
    p.add_argument("--prefix", default=DEFAULT_PREFIX,
                   help="filename prefix for every output file")
    p.add_argument("--a7-burden", default=None,
                   help="explicit path to a7_virome_anellovirus_burden.tsv")
    p.add_argument("--a7-key", default=None,
                   help="explicit path to a7_virome_sample_key.tsv")
    p.add_argument("--anello-accessions", default=DEFAULT_ANELLO_ACC_FILE,
                   help="file of anellovirus accessions / reference ids")
    p.add_argument("--chimp-accessions", default=DEFAULT_CHIMP_ACC,
                   help="comma-separated chimpanzee-isolate accessions to flag")
    p.add_argument("--include-chimp", action="store_true",
                   help="fold the chimpanzee-isolate references into the group "
                        "aggregates instead of keeping them as a control")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--samtools-threads", type=int, default=1,
                   help="samtools -@ value")
    p.add_argument("--bam-glob", default="bam/*.bam",
                   help="glob for BAMs relative to a run directory")
    p.add_argument("--exclude-flags", default="0x904",
                   help="samtools view -F value (unmapped, secondary, "
                        "supplementary); duplicates are kept on purpose")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="samtools -q value; 0 on purpose, a MAPQ filter would "
                        "hide the cross-mapping this module looks for")
    p.add_argument("--window", type=int, default=200,
                   help="sliding window in bp for the pile-up detector")
    p.add_argument("--max-ref-len", type=int, default=100000,
                   help="skip any 'anellovirus' reference longer than this. "
                        "Anellovirus genomes are ~3.7 kb, so a longer hit is a "
                        "misclassified reference, and per-base coverage on a "
                        "chromosome-sized one would exhaust memory")
    p.add_argument("--pooled-bins", type=int, default=20,
                   help="bins for the pooled relative-position distribution")
    p.add_argument("--min-reads", type=int, default=5,
                   help="reads below which a pair stays too_few_reads")
    p.add_argument("--min-breadth", type=float, default=0.10,
                   help="breadth needed for a real_like verdict")
    p.add_argument("--pileup-frac", type=float, default=0.50,
                   help="max_window_fraction at or above which a pair is "
                        "pileup_like")
    p.add_argument("--dup-frac", type=float, default=0.50,
                   help="duplicate_position_fraction at or above which a pair "
                        "is duplicate_like")
    p.add_argument("--hotspot-tol", type=float, default=0.05,
                   help="relative distance within which two samples count as "
                        "sharing one hotspot")
    p.add_argument("--min-clip", type=int, default=10,
                   help="soft clip length in bp counted as clipped")
    p.add_argument("--lowcomp-frac", type=float, default=0.80,
                   help="two-base share at or above which a read is called "
                        "low complexity")
    p.add_argument("--test-groups", default="HIV,HL",
                   help="the two group labels compared")
    p.add_argument("--limit", type=int, default=0,
                   help="audit only the first N BAMs per run (debugging)")
    p.add_argument("--no-run-name-group", action="store_true",
                   help="do not fall back to the run directory name for TCL")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.window < 1:
        args.window = 200
    if args.pooled_bins < 1:
        args.pooled_bins = 20
    if args.min_reads < 1:
        args.min_reads = 1

    test_groups = [g.strip() for g in args.test_groups.split(",") if g.strip()]
    if len(test_groups) != 2:
        print("WARN: --test-groups needs exactly two labels, got %r; using HIV,HL"
              % args.test_groups)
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
    print("window   %d bp | min-reads %d | pileup-frac %.2f | dup-frac %.2f | "
          "min-breadth %.2f" % (args.window, args.min_reads, args.pileup_frac,
                                args.dup_frac, args.min_breadth))
    print("")

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            warn_missing("writable output directory", args.outdir)
            return 0
    paths = {
        "pair": os.path.join(args.outdir, out_name(args.prefix, "by_pair.tsv")),
        "group": os.path.join(args.outdir, out_name(args.prefix, "by_group.tsv")),
        "pooled": os.path.join(args.outdir,
                               out_name(args.prefix, "pooled_positions.tsv")),
        "key": os.path.join(args.outdir, out_name(args.prefix, "sample_key.tsv")),
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

    a7_key, a7_burden = load_a7_context(args)

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
                         a7_key.get(name, "NA")])
    # Only write the key when there is an identifier to record: an empty file
    # headed "CONTAINS IDENTIFIERS" is misleading in a directory listing and in
    # run_all.sh's closing identifier summary.
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
    print("")

    # ---- pass 2: audit ----------------------------------------------------- #
    pairs = []
    oversize = set()
    n_no_index = 0
    n_failed = 0
    n_no_anello_ref = 0
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
                                 os.path.join(run_dir, "bam", "<sample>.bam.bai"))
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
                              "--max-ref-len %d; it cannot be an anellovirus "
                              "and is skipped" % (rname, seqlen,
                                                  args.max_ref_len))
                    continue
                if from_idxstats and not mapped:
                    continue
                anello_refs[rname] = (seqlen, mapped)
            if not anello_refs:
                n_no_anello_ref += 1
                continue
            ref_lens = dict((r, v[0]) for r, v in anello_refs.items())
            accs, err = stream_sample(args.samtools, bam_path,
                                      list(anello_refs), ref_lens, args,
                                      sample, anon)
            if accs is None:
                n_failed += 1
                if n_failed == 1:
                    # err was already masked inside stream_sample
                    print("WARN: samtools view failed for %s (%s); that sample "
                          "is skipped" % (anon, err))
                continue
            done += 1
            for rname, (seqlen, mapped) in sorted(anello_refs.items()):
                acc = accs.get(rname)
                if acc is None:
                    acc = new_acc(seqlen)
                m = pair_metrics(acc, args)
                verdict, flags, metrics = verdict_of(m, args)
                a7_anon = a7_key.get(sample, "NA")
                ctx = a7_burden.get(a7_anon, {})
                pairs.append({
                    "run": base,
                    "sample_anon": anon,
                    "group": group,
                    "reference_id": rname,
                    "ref_label": ref_label(rname, refmap),
                    "chimp": is_chimp_ref(rname, refmap, chimp_norm, chimp_ids),
                    "idxstats_mapped": mapped,
                    "m": m,
                    "verdict": verdict,
                    "flags": flags,
                    "metrics": metrics,
                    "a7_anon": a7_anon,
                    "a7_reads": ctx.get("reads", "NA"),
                    "a7_richness": ctx.get("richness", "NA"),
                })
        print("  %-46s %3d/%3d samples audited" % (base[:46], done, len(bams)))
    if n_no_index:
        print("WARN: %d BAM(s) without an index were skipped" % n_no_index)
    if n_no_anello_ref:
        print("NOTE: %d sample(s) had no anellovirus reference with a mapped "
              "read" % n_no_anello_ref)

    if not pairs:
        print("")
        print("No (sample, anellovirus reference) pair carried a read; writing "
              "headed, empty tables.")
        write_headers_only(paths, args)
        print("wrote: %s, %s, %s"
              % (paths["pair"], paths["group"], paths["pooled"]))
        return 0

    human_pairs = [p for p in pairs if args.include_chimp or not p["chimp"]]
    chimp_pairs = [p for p in pairs if p["chimp"] and not args.include_chimp]
    pairs_by_set = {"human_anello": human_pairs}
    if chimp_pairs:
        pairs_by_set["chimp_flagged"] = chimp_pairs

    write_by_pair(paths["pair"], pairs, args)
    group_rows, summaries, tests = build_group_rows(pairs_by_set, args, g1, g2)
    write_by_group(paths["group"], group_rows, args, g1, g2)
    pooled_rows, scopes, hotspots = build_pooled_rows(human_pairs, args)
    write_pooled(paths["pooled"], pooled_rows, args)

    # ---- stdout summary ---------------------------------------------------- #
    print("")
    print("-- pairs audited --")
    n_chimp = sum(1 for p in pairs if p["chimp"])
    if args.include_chimp:
        print("   %d (sample, reference) pair(s); --include-chimp, so the %d "
              "chimpanzee-flagged pair(s) are folded in, not held out as a "
              "control" % (len(pairs), n_chimp))
    else:
        print("   %d (sample, reference) pair(s), %d human-anellovirus + %d "
              "chimpanzee-flagged control"
              % (len(pairs), len(human_pairs), len(chimp_pairs)))
    counts = dict((v, sum(1 for p in pairs if p["verdict"] == v))
                  for v in VERDICTS)
    print("   verdicts: " + ", ".join("%s=%d" % (v, counts[v]) for v in VERDICTS))

    print("")
    print("-- per group (human anellovirus references) --")
    print("   %-5s %7s %6s %8s %9s %8s %7s %7s %9s"
          % ("group", "samples", "pairs", "reads", "too_few", "pileup",
             "dup", "real", "hot_bin"))
    for label in [g for g in GROUP_ORDER + ["ALL"] if g in summaries]:
        row = summaries[label]
        print("   %-5s %7s %6s %8s %9s %8s %7s %7s %9s"
              % (label, row["n_samples"], row["n_pairs"], row["n_reads"],
                 row["n_too_few_reads"], row["n_pileup_like"],
                 row["n_duplicate_like"], row["n_real_like"],
                 row["pooled_hottest_bin_fraction"]))
    print("   hot_bin = share of that group's pooled read starts in the single "
          "hottest of %d relative-position bins" % args.pooled_bins)

    for label in (g1, g2):
        st = scopes.get(label)
        if not st or not st["total"]:
            print("   %s: no pooled anellovirus read start to distribute" % label)
            continue
        hb, frac = hottest_bin(st)
        nz = sum(1 for c in st["counts"] if c)
        print("   %s: %d pooled read start(s) from %d sample(s) over %d/%d "
              "bins; hottest bin %d (rel %.2f-%.2f) holds %.1f%%"
              % (label, st["total"], len(st["all_samples"]), nz,
                 args.pooled_bins, hb, float(hb) / args.pooled_bins,
                 float(hb + 1) / args.pooled_bins, 100.0 * frac))

    print("")
    print("-- shared hotspots across independent samples (tol %.3f rel) --"
          % args.hotspot_tol)
    wanted_scopes = ["ALL", g1, g2]
    shared = [h for h in hotspots if h["n_at"] >= 2 and h["scope"] in wanted_scopes]
    shared.sort(key=lambda h: (wanted_scopes.index(h["scope"]), -h["n_at"],
                               -h["frac"], h["reference_id"]))
    if not shared:
        print("   none: no reference had two samples agreeing on a hotspot")
    for h in shared[:12]:
        approx = (h["pos"] * h["ref_len"]) if h["pos"] is not None else None
        print("   %-5s %-30s %2d/%2d pair(s) at rel %s (~%s bp of %d), "
              "median max_window_frac there %s"
              % (h["scope"], h["reference_id"][:30], h["n_at"], h["n_pairs"],
                 fnum(h["pos"], 3),
                 "NA" if approx is None else str(int(approx)), h["ref_len"],
                 fnum(h["conc"], 2)))
    if shared:
        print("   a shared hotspot only argues artefact when the pairs on it "
              "are concentrated:")
        print("   a high n at a low median max_window_frac is just where "
              "dispersed coverage peaked.")

    print("")
    print("-- %s vs %s tests (standard library; n is beside every p) --" % (g1, g2))
    for name, res, n1, n2 in tests:
        if res is None:
            print("   %-38s not tested (n %s=%d, %s=%d)"
                  % (name, g1, n1, g2, n2))
            continue
        if "p" in res and "U1" in res:
            print("   %-38s U=%s  n %s=%d %s=%d  p=%s  (medians %s vs %s)"
                  % (name, fnum(res["U1"], 1), g1, n1, g2, n2, fp(res["p"]),
                     fnum(res["median1"], 3), fnum(res["median2"], 3)))
        else:
            odds = ("inf" if res["odds_ratio"] == float("inf")
                    else fnum(res["odds_ratio"], 3))
            print("   %-38s OR=%s  n %s=%d %s=%d  p=%s  (%d/%d vs %d/%d)"
                  % (name, odds, g1, n1, g2, n2, fp(res["p"]),
                     res["a"], res["n1"], res["c"], res["n2"]))

    n_too_few = counts["too_few_reads"]
    print("")
    print("READ THIS BEFORE QUOTING ANY OF IT:")
    print("   %d/%d pair(s) carry fewer than %d reads and stay too_few_reads. "
          "At that depth" % (n_too_few, len(pairs), args.min_reads))
    print("   no method separates real virus from cross-mapping, so those pairs "
          "are unresolved,")
    print("   not negative. A group whose anellovirus signal is mostly "
          "too_few_reads has not been")
    print("   validated read by read, in either direction.")

    print("")
    print("wrote:")
    for key in ("pair", "group", "pooled"):
        print("  %s" % paths[key])
    if key_rows:
        print("  %s" % paths["key"])
        print("REMINDER: %s contains real sample identifiers - do not commit or "
              "email it." % os.path.basename(paths["key"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
