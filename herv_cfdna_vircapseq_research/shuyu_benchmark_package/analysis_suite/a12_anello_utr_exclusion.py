#!/usr/bin/env python3
"""a12_anello_utr_exclusion.py -- does the HIV+ vs HL anellovirus difference
survive when every read that could have come from the conserved UTR is removed?

THE QUESTION
  a7 reported an anellovirus burden higher in HIV+ than in HL inside the WGS
  cohort: richness >= 3 in 18/37 HIV+ against 0/23 HL, Fisher p = 2.3e-05.
  a10 then audited those reads and found the signal is dominated by artefact -
  median max_window_fraction 0.729, median duplicate_position_fraction 0.552,
  96% of reads soft-clipped >= 10 bp, median breadth 0.068, the CHIMPANZEE
  references (which no human carries) showing the same profile with 1952 reads
  in the HIV group, and one shared hotspot at relative position 0.81-0.93
  recurring across unrelated references and independent samples (NC_014075.1
  17/19 pairs at rel 0.813, NC_002076.2 17/24 at 0.846, NC_014078.1 13/20 at
  0.904). In a ~3.7 kb anellovirus genome relative 0.75-1.0 is the conserved
  UTR: the expected cross-mapping magnet. a10 also showed that demanding a
  read-validated pair drops the contrast to 8/37 vs 1/23, p = 0.134.

  This module tests that mechanism directly and decisively. It does NOT
  re-align anything. It re-counts from the existing BAMs under a cumulative
  evidence filter, then re-runs a7's two headline tests on whatever survives,
  and prints a7's published numbers beside the recomputed ones so the change is
  explicit and quotable.

THE FILTER LADDER (cumulative, so the reader can see what each rung costs)
  0 all_reads              every read on an anellovirus reference that passes
                           samtools view -F <--exclude-flags> -q <--min-mapq>.
                           MAPQ is NOT filtered by default: multi-mapping is
                           the artefact and filtering it first would hide it.
  1 utr_excluded           drop any read whose ALIGNED SEGMENT (POS through the
                           reference bases consumed by M/=/X/D/N) overlaps the
                           UTR window --utr-lo .. --utr-hi in relative
                           coordinates along that reference (default 0.75-1.00).
                           Overlap, not containment: a read that merely clips
                           into the window is still a UTR-compatible read.
  2 deduplicated           one read per distinct POS+CIGAR. Identical starts are
                           PCR duplicates or one locus copied over, not a genome.
  3 entropy_filtered       drop reads whose --entropy-k mer Shannon entropy
                           (default 3-mer, log base 2) is below --min-entropy
                           (default 1.2 bits). A read with no SEQ cannot be
                           judged and is KEPT, counted separately.
  4 unambiguous            drop reads whose best and second-best alignment score
                           tie (AS == XS; the degenerate AS < XS is treated the
                           same). No AS tag means the aligner reported no score
                           and the read is KEPT, counted separately.
  5 multi_position         require >= --min-distinct-positions (default 2)
                           mutually NON-OVERLAPPING aligned segments per
                           (sample, reference). One pile-up position is not
                           evidence of a genome, however clean the reads are.
                           This is the only pair-level rung: it keeps or drops
                           the whole (sample, reference) pair.

WHAT IT THEN COMPUTES ON THE SURVIVORS
  Per sample and per rung: total surviving reads (burden), the same per million
  (--norm-source denominator, default the sample's own idxstats mapped total,
  identical to a7's default), and richness = the number of references carrying
  >= --min-reads surviving reads (default 5).
  Per rung it re-runs a7's two headline tests with the same standard-library
  implementations a7 and a10 use: Mann-Whitney U on burden (raw and per
  million) and Fisher exact 2x2 on richness >= --richness-cut (default 3).
  Every p is written and printed with the n of each group beside it.

  UTR-WINDOW SWEEP. The whole ladder is repeated with --utr-lo at each value of
  --sweep-lo (default 0.70, 0.75, 0.80, 0.85) and once with the UTR rung
  disabled entirely, and the same before/after statistics are reported per
  window, so the conclusion cannot be an artefact of where the boundary was
  drawn. When a10's pooled-position table is in --indir, each window also
  carries how many of a10's observed shared hotspots it actually covers.

  CHIMPANZEE CONTROL. The same pipeline runs on the chimpanzee-isolate
  references (--chimp-accessions; NC_014069.1, NC_014077.1, NC_014480.2) and is
  reported alongside every rung and every window. No human sample carries
  chimpanzee TTV, so whatever survives there is this filter's residual
  false-positive rate, measured rather than assumed.

WHAT THIS CAN AND CANNOT SUPPORT
  THIS TESTS ONE MECHANISM: cross-mapping onto the conserved terminal UTR.
  CAN: it can show that the group difference is carried by UTR-compatible,
  duplicated, low-complexity or ambiguously-placed reads, i.e. that it does not
  survive the removal of that one mechanism. That is a strong negative result.
  CAN: it can bound the residual false-positive rate with the chimpanzee arm.
  CANNOT: SURVIVING THE FILTER IS NECESSARY BUT NOT SUFFICIENT TO CALL REAL
  VIRUS. Reads outside the UTR can still be cross-mapped from an anellovirus
  the panel does not contain, from another small circular DNA virus, or from an
  unmasked human repeat. The confirmatory test is a COMPETITIVE REALIGNMENT of
  the surviving reads against a fuller reference set (all Anelloviridae plus a
  human and vector decoy), which this module deliberately does not attempt.
  CANNOT: it cannot rescue a true infection whose only represented sequence is
  the conserved UTR - such a read is removed by construction, so rung 1 is
  conservative against detection, and the sweep exists to show by how much.
  CANNOT: at a handful of reads no method separates real virus from artefact.
  A pair dropped at rung 5 is unresolved, not proven absent.

WHAT IT WRITES (tab separated, pure ASCII, into --outdir, --prefix prepended)
  <prefix>_utr_exclusion_ladder.tsv      sample, pair and read counts surviving
                                         at every rung, per group and per
                                         reference set (human + chimp control)
  <prefix>_utr_exclusion_by_sample.tsv   one row per sample: reads at each rung,
                                         richness and burden before and after,
                                         raw and per million, chimp residual
  <prefix>_utr_exclusion_group_test.tsv  the BEFORE / AFTER block: a7's
                                         published numbers, a7's own table if
                                         it is in --indir, and the recomputed
                                         Mann-Whitney and Fisher rows per rung
  <prefix>_utr_window_sweep.tsv          the same statistics for every UTR
                                         window and for the filter disabled
  <prefix>_utr_exclusion_sample_key.tsv  real -> anonymous mapping. THE ONLY
                                         file with real sample identifiers; its
                                         first line says so, and it is written
                                         only when a real name was actually
                                         read. Do not commit or email it.
  With the default prefix: a12_utr_exclusion_ladder.tsv,
  a12_utr_exclusion_by_sample.tsv, a12_utr_exclusion_group_test.tsv,
  a12_utr_window_sweep.tsv, a12_utr_exclusion_sample_key.tsv. Samples are
  anonymised to S001..Snnn by sorted real name, the same three-digit rule a10
  uses, so a10's ids match when the same BAM set is processed; the
  a7_sample_anon column (filled from a7's key in --indir) joins to a7's
  two-digit ids.

  Standard library only. No figures, so matplotlib is not imported. No network.
  Any missing input prints "WARN: <what> missing at <path>, skipping" and the
  module still writes headed tables and exits 0.

EXAMPLE
  python3 a12_anello_utr_exclusion.py \
      --run /path/to/runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only \
      --refmap /path/to/runs/shuyu_masked_panel_hg38_herv_line1_refixed/ref/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv \
      --indir /path/to/runs/panel_report_20260725/suite_out \
      --outdir /path/to/runs/panel_report_20260725/suite_out \
      --utr-lo 0.75 --utr-hi 1.00 --min-reads 5 --richness-cut 3

  # RUNS_ROOT is honoured for the default input paths, as in run_all.sh:
  RUNS_ROOT=/real/run/root python3 a12_anello_utr_exclusion.py --outdir <scratch>

Written 2026-07-27.
"""
from __future__ import annotations

import argparse
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
DEFAULT_PREFIX = "a12"
DEFAULT_ANELLO_ACC_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "anello_accessions.txt")
DEFAULT_CHIMP_ACC = "NC_014069.1,NC_014077.1,NC_014480.2"
DEFAULT_SWEEP_LO = "0.70,0.75,0.80,0.85"

# Upstream tables looked for in --indir (exact name first, then a suffix match).
A7_KEY_NAME = "a7_virome_sample_key.tsv"
A7_KEY_SUFFIX = "_sample_key.tsv"
A7_BURDEN_NAME = "a7_virome_anellovirus_burden.tsv"
A7_BURDEN_SUFFIX = "anellovirus_burden.tsv"
A10_POOLED_NAME = "anello_read_audit_pooled_positions.tsv"
A10_POOLED_SUFFIX = "pooled_positions.tsv"

# a7's published headline, quoted so the change is explicit. Overridable.
DEFAULT_A7_PUB_G1 = "18/37"
DEFAULT_A7_PUB_G2 = "0/23"
DEFAULT_A7_PUB_P = "2.3e-05"

ANELLO_KEYWORDS = [
    "anello", "torque teno", "torque-teno", "torquetenovirus",
    "transfusion transmitted virus", "transfusion-transmitted virus",
    "small anellovirus", "tt virus", "ttv", "ttmv", "ttmdv", "sen virus",
]
# Short keywords that must not match inside a longer word.
SHORT_KEYWORDS = set(["ttv", "ttmv", "ttmdv", "tt virus"])

ACC_RE = re.compile(
    r"(?<![A-Za-z0-9])([A-Z]{2}_\d{5,8}|[A-Z]{1,2}\d{5,6})(?:\.\d+)?(?![A-Za-z0-9])")

GROUP_ORDER = ["HIV", "HL", "TCL", "NA"]
REF_SETS = ["human_anello", "chimp_flagged"]

RUNG_NAMES = [
    "all_reads",
    "utr_excluded",
    "deduplicated",
    "entropy_filtered",
    "unambiguous",
    "multi_position",
]
N_RUNGS = len(RUNG_NAMES)
FINAL_RUNG = N_RUNGS - 1
# Derived once so the per-rung column names and the writer can never diverge.
RUNG_READ_COLUMNS = ["reads_rung%d_%s" % (i, name)
                     for i, name in enumerate(RUNG_NAMES)]

TODAY = datetime.date.today().isoformat()
SCRIPT = os.path.basename(__file__)

LADDER_COLUMNS = [
    "row_type", "ref_set", "scope", "window_label", "utr_lo", "utr_hi",
    "utr_filter_applied", "rung", "rung_name", "rung_filter",
    "n_samples_in_scope", "n_samples_with_surviving_read",
    "n_samples_with_detected_reference", "n_pairs_with_surviving_read",
    "n_pairs_at_min_reads", "n_reads", "reads_dropped_vs_previous_rung",
    "frac_reads_of_rung0", "frac_pairs_of_rung0",
    "median_reads_per_sample", "max_reads_per_sample",
    "min_reads_per_reference", "richness_cut", "n_samples_richness_ge_cut",
    "frac_samples_richness_ge_cut", "note",
]

BY_SAMPLE_COLUMNS = [
    "sample_anon", "group", "run", "window_label", "utr_lo", "utr_hi",
] + RUNG_READ_COLUMNS + [
    "n_refs_with_read_rung0", "n_refs_with_read_final",
    "min_reads_per_reference", "richness_rung0", "richness_final",
    "richness_cut", "richness_ge_cut_rung0", "richness_ge_cut_final",
    "norm_source", "norm_denominator", "burden_rpm_rung0", "burden_rpm_final",
    "top_reference_final", "top_reference_reads_final",
    "max_nonoverlapping_positions_final",
    "chimp_reads_rung0", "chimp_reads_final", "chimp_richness_final",
    "a7_sample_anon", "a7_anello_reads_total", "a7_anello_richness", "note",
]

TEST_COLUMNS = [
    "row_type", "source", "ref_set", "window_label", "utr_lo", "utr_hi",
    "rung", "rung_name", "metric", "test", "group1", "n1", "value_group1",
    "group2", "n2", "value_group2", "statistic", "statistic_value",
    "p_value", "effect", "note",
]

SWEEP_COLUMNS = [
    "row_type", "window_label", "utr_lo", "utr_hi", "utr_filter_applied",
    "is_primary_window", "ref_set", "rung", "rung_name",
    "n_samples_in_scope", "n_samples_with_surviving_read",
    "n_pairs_with_surviving_read", "n_pairs_at_min_reads", "n_reads",
    "frac_reads_of_rung0", "richness_cut", "min_reads_per_reference",
    "group1", "n1", "richness_ge_cut_group1", "group2", "n2",
    "richness_ge_cut_group2", "fisher_odds_ratio", "fisher_p",
    "burden_median_group1", "burden_median_group2", "mann_whitney_U_group1",
    "mann_whitney_p", "chimp_reads_surviving", "chimp_samples_with_read",
    "chimp_pairs_at_min_reads", "a10_shared_hotspots_total",
    "a10_shared_hotspots_inside_window", "note",
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


def write_tsv(path, comments, header, rows):
    with open(path, "w", encoding="ascii", errors="replace", newline="") as fh:
        for line in comments:
            fh.write("# " + to_ascii(line) + "\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(to_ascii(c) for c in row) + "\n")


def parse_fraction_pair(text, fallback):
    """'18/37' -> (18, 37). fallback on anything unparseable."""
    try:
        left, right = str(text).split("/", 1)
        a = int(left.strip())
        n = int(right.strip())
        if a < 0 or n < 0 or a > n:
            return fallback
        return a, n
    except (ValueError, AttributeError):
        return fallback


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
# statistics (standard library only; identical implementations to a7 / a10)
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


def odds_text(res):
    if res is None or res.get("odds_ratio") is None:
        return "NA"
    if res["odds_ratio"] == float("inf"):
        return "inf"
    return fnum(res["odds_ratio"], 4)


# --------------------------------------------------------------------------- #
# reference identification (identical rules to a10)
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


def sample_of_bam(bam_path):
    base = os.path.basename(bam_path)
    return base[:-4] if base.endswith(".bam") else base


def bam_is_indexed(bam):
    return (os.path.exists(bam + ".bai")
            or os.path.exists(os.path.splitext(bam)[0] + ".bai")
            or os.path.exists(bam + ".csi"))


def idxstats_table(run_dir, sample):
    """({refname: (seqlen, mapped)}, total_mapped, total_unmapped) or ({},0,0)."""
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
        return {}, 0, 0
    out = {}
    total_mapped = 0
    total_unmapped = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                try:
                    seqlen = int(p[1])
                    mapped = int(p[2])
                    unmapped = int(p[3])
                except ValueError:
                    continue
                if p[0] == "*":
                    total_unmapped += unmapped
                    continue
                total_mapped += mapped
                total_unmapped += unmapped
                out[p[0]] = (seqlen, mapped)
    except OSError:
        return {}, 0, 0
    return out, total_mapped, total_unmapped


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


def find_filtered_category_counts(results_dir):
    """Headline deduplicated category table, with or without a filename prefix."""
    hits = sorted(glob.glob(os.path.join(results_dir,
                                         "*filtered_category_counts.tsv")))
    hits = [h for h in hits if "record_category_counts" not in os.path.basename(h)]
    if hits:
        return hits[0]
    hits = sorted(glob.glob(os.path.join(results_dir, "*_category_counts.tsv")))
    hits = [h for h in hits if "record_category_counts" not in os.path.basename(h)]
    return hits[0] if hits else None


def parse_category_counts(path):
    """-> dict sample -> total of the category columns. None if unreadable."""
    table = {}
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = None
            for row in reader:
                if not row:
                    continue
                if header is None:
                    header = [c.strip() for c in row]
                    continue
                total = 0
                for idx in range(1, min(len(row), len(header))):
                    try:
                        total += int(float(row[idx]))
                    except (ValueError, TypeError):
                        continue
                table[row[0].strip()] = total
    except (OSError, csv.Error):
        return None
    return table


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
# upstream context (all optional, never fatal)
# --------------------------------------------------------------------------- #
def resolve_in_indir(indir, exact_name, suffix, must_contain=None,
                     skip_prefix=None):
    """Exact filename in indir first, then any file ending with suffix.

    must_contain / skip_prefix keep a generic suffix such as "_sample_key.tsv"
    from resolving to this module's own output when --indir == --outdir.
    """
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
        if not name.endswith(suffix):
            continue
        if must_contain and must_contain not in name:
            continue
        if skip_prefix and name.startswith(skip_prefix):
            continue
        full = os.path.join(indir, name)
        if os.path.isfile(full):
            return full
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


def load_a7_key(args):
    """real sample name -> a7 anonymous id. {} when the key is absent."""
    key_map = {}
    path = args.a7_key or resolve_in_indir(args.indir, A7_KEY_NAME,
                                           A7_KEY_SUFFIX, must_contain="a7",
                                           skip_prefix=args.prefix)
    if not path or not os.path.exists(path):
        warn_missing("a7 sample key", path or os.path.join(
            args.indir or "<indir>", A7_KEY_NAME))
        return key_map
    rows = read_commented_tsv(path)
    if rows is None:
        warn_missing("readable a7 sample key", path)
        return key_map
    for row in rows:
        real = (row.get("real_sample") or "").strip()
        anon = (row.get("anon_sample") or "").strip()
        if real and anon:
            key_map[real] = anon
    return key_map


def load_a7_burden(args):
    """a7 anon id -> {group, reads, richness}. {} when the table is absent.

    The table carries a7's own anonymous ids and its own group column, so a7's
    published contrast can be recomputed from it without any identifier join.
    """
    burden = {}
    path = args.a7_burden or resolve_in_indir(args.indir, A7_BURDEN_NAME,
                                              A7_BURDEN_SUFFIX,
                                              skip_prefix=args.prefix)
    if not path or not os.path.exists(path):
        warn_missing("a7 anellovirus burden table", path or os.path.join(
            args.indir or "<indir>", A7_BURDEN_NAME))
        return burden
    rows = read_commented_tsv(path)
    if rows is None:
        warn_missing("readable a7 anellovirus burden table", path)
        return burden
    for row in rows:
        anon = (row.get("sample") or "").strip()
        if not anon:
            continue
        reads = (row.get("anello_reads_human_total") or "NA").strip()
        rich = (row.get("anello_richness_human") or "NA").strip()
        try:
            rich_i = int(float(rich))
        except (ValueError, TypeError):
            rich_i = None
        try:
            reads_i = int(float(reads))
        except (ValueError, TypeError):
            reads_i = None
        burden[anon] = {
            "group": (row.get("group") or "NA").strip(),
            "reads": reads,
            "richness": rich,
            "reads_int": reads_i,
            "richness_int": rich_i,
        }
    return burden


def load_a10_hotspots(args):
    """[(reference_id, rel_pos, n_pairs_at_hotspot), ...] from a10, or None.

    Only the ALL-scope reference_hotspot rows with at least two agreeing pairs
    are kept: those are the shared hotspots a10 called artefact evidence. No
    sample-level join is involved, so no identifier is touched.
    """
    path = args.a10_pooled or resolve_in_indir(args.indir, A10_POOLED_NAME,
                                               A10_POOLED_SUFFIX,
                                               skip_prefix=args.prefix)
    if not path or not os.path.exists(path):
        warn_missing("a10 pooled-position table", path or os.path.join(
            args.indir or "<indir>", A10_POOLED_NAME))
        return None
    rows = read_commented_tsv(path)
    if rows is None:
        warn_missing("readable a10 pooled-position table", path)
        return None
    out = []
    for row in rows:
        if (row.get("row_type") or "").strip() != "reference_hotspot":
            continue
        if (row.get("scope") or "").strip() != "ALL":
            continue
        try:
            pos = float(row.get("hotspot_rel_pos"))
        except (TypeError, ValueError):
            continue
        try:
            n_at = int(float(row.get("n_pairs_at_hotspot") or 0))
        except (TypeError, ValueError):
            n_at = 0
        if n_at < 2:
            continue
        out.append(((row.get("reference_id") or "NA").strip(), pos, n_at))
    return out


# --------------------------------------------------------------------------- #
# read-level parsing
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


def ref_consumed(ops):
    """Reference bases consumed by a CIGAR (M/=/X/D/N)."""
    total = 0
    for num, ch in ops:
        if ch in "M=XDN":
            total += num
    return total


def kmer_entropy(seq, k):
    """Shannon entropy (log base 2) over the overlapping k-mers of a read.

    None when the read has no SEQ or is shorter than k: such a read cannot be
    judged and is kept by the entropy rung rather than silently dropped.
    """
    if not seq or seq == "*":
        return None
    s = seq.upper()
    n = len(s) - k + 1
    if n < 1:
        return None
    counts = {}
    for i in range(n):
        km = s[i:i + k]
        counts[km] = counts.get(km, 0) + 1
    total = float(n)
    h = 0.0
    for c in counts.values():
        p = float(c) / total
        h -= p * math.log(p, 2)
    return h


def max_nonoverlapping(reads):
    """Largest set of aligned segments that pairwise do not overlap.

    Classic interval-scheduling greedy (earliest end first), so the answer is
    the maximum, not an arbitrary chain. One pile-up position gives 1 however
    many reads sit on it; two clearly separated fragments give 2.
    """
    if not reads:
        return 0
    segs = sorted((r[1], r[0]) for r in reads)
    count = 0
    last_end = None
    for end, start in segs:
        if last_end is None or start >= last_end:
            count += 1
            last_end = end
    return count


def stream_sample(samtools, bam, refs, ref_lens, args, sample, anon, chunk=200):
    """Stream one BAM over its anellovirus references, keeping per-read facts.

    Returns (reference_id -> list of read tuples, error_string, counters). The
    dict is None when samtools could not be run at all. Each read tuple is
    (start0, end0, pos_cigar_key, kmer_entropy_or_None, ambiguous_bool).
    sample / anon are used only to mask the real sample name out of samtools'
    error text before it is truncated, so a truncated path cannot leak half a
    name.
    """
    reads = {}
    counters = {"no_seq": 0, "no_as_tag": 0, "no_cigar": 0}
    ref_list = sorted(refs)
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
                redact(exc, sample, anon)), counters
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
                pos0 = int(f[3]) - 1
            except ValueError:
                continue
            if pos0 < 0:
                pos0 = 0
            cigar = f[5]
            ops = parse_cigar(cigar)
            seq = f[9]
            if ops:
                span = ref_consumed(ops)
            else:
                counters["no_cigar"] += 1
                span = len(seq) if seq and seq != "*" else 1
            if span < 1:
                span = 1
            end0 = pos0 + span
            if end0 > ref_len:
                end0 = ref_len
            if end0 <= pos0:
                end0 = pos0 + 1
            if seq and seq != "*":
                ent = kmer_entropy(seq, args.entropy_k)
            else:
                counters["no_seq"] += 1
                ent = None
            a_score = None
            x_score = None
            for tag in f[11:]:
                if tag.startswith("AS:i:"):
                    try:
                        a_score = int(tag[5:])
                    except ValueError:
                        a_score = None
                elif tag.startswith("XS:i:"):
                    try:
                        x_score = int(tag[5:])
                    except ValueError:
                        x_score = None
            if a_score is None:
                counters["no_as_tag"] += 1
                ambiguous = False
            else:
                ambiguous = (x_score is not None and a_score <= x_score)
            reads.setdefault(rname, []).append(
                (pos0, end0, f[3] + ":" + cigar, ent, ambiguous))
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
            return None, ("samtools view rc=%d %s"
                          % (rc, err[:160] or "no stderr")), counters
    return reads, "", counters


# --------------------------------------------------------------------------- #
# the filter ladder
# --------------------------------------------------------------------------- #
def rung_filters(cfg, args):
    """Human-readable description of every rung under one window config."""
    if cfg["on"]:
        utr = ("aligned segment does not overlap relative %.2f-%.2f of the "
               "reference" % (cfg["lo"], cfg["hi"]))
    else:
        utr = "UTR exclusion DISABLED for this window (rung 1 == rung 0)"
    return [
        "every read on an anellovirus reference passing samtools view -F %s "
        "-q %d" % (args.exclude_flags, args.min_mapq),
        utr,
        "one read per distinct POS+CIGAR",
        "%d-mer Shannon entropy >= %.2f bits (a read with no SEQ is kept)"
        % (args.entropy_k, args.min_entropy),
        "AS > XS, or no XS tag; a tie AS == XS is dropped (no AS tag is kept)",
        ">= %d non-overlapping aligned segments per (sample, reference); the "
        "whole pair is dropped otherwise" % args.min_distinct_positions,
    ]


def apply_ladder(pair, cfg, args):
    """(reads surviving each rung, max non-overlapping positions, distinct starts).

    Cumulative: each rung filters what the previous one left. Rung 5 is the only
    pair-level rung, so it either keeps everything rung 4 left or nothing.
    """
    reads = pair["reads"]
    counts = [0] * N_RUNGS
    counts[0] = len(reads)

    if cfg["on"]:
        ref_len = float(pair["ref_len"])
        lo_bp = cfg["lo"] * ref_len
        hi_bp = cfg["hi"] * ref_len
        r1 = [r for r in reads if not (r[1] > lo_bp and r[0] < hi_bp)]
    else:
        r1 = list(reads)
    counts[1] = len(r1)

    seen = set()
    r2 = []
    for r in r1:
        if r[2] in seen:
            continue
        seen.add(r[2])
        r2.append(r)
    counts[2] = len(r2)

    r3 = [r for r in r2 if r[3] is None or r[3] >= args.min_entropy]
    counts[3] = len(r3)

    r4 = [r for r in r3 if not r[4]]
    counts[4] = len(r4)

    nonov = max_nonoverlapping(r4)
    ndistinct = len(set(r[0] for r in r4))
    counts[5] = len(r4) if nonov >= args.min_distinct_positions else 0
    return counts, nonov, ndistinct


def evaluate_config(pairs, cfg, universe, args):
    """Run the whole ladder under one UTR-window configuration.

    universe is every sample that was successfully streamed, including the ones
    with no anellovirus read at all: those samples are part of the cohort and
    must sit in the denominator of every proportion, or 18/37 could never be
    reproduced.
    """
    samples = {}
    for ref_set in REF_SETS:
        samples[ref_set] = {}
        for anon in sorted(universe):
            info = universe[anon]
            samples[ref_set][anon] = {
                "group": info["group"],
                "run": info["run"],
                "denom": info["denom"],
                "reads": [0] * N_RUNGS,
                "nrefs": [0] * N_RUNGS,
                "detected": [0] * N_RUNGS,
                "top_ref": "NA",
                "top_reads": 0,
                "nonov": 0,
            }
    per_pair = []
    for pair in pairs:
        counts, nonov, ndistinct = apply_ladder(pair, cfg, args)
        per_pair.append({"pair": pair, "counts": counts, "nonov": nonov,
                         "ndistinct": ndistinct})
        rec = samples.get(pair["ref_set"], {}).get(pair["sample_anon"])
        if rec is None:
            continue
        for i in range(N_RUNGS):
            rec["reads"][i] += counts[i]
            if counts[i] > 0:
                rec["nrefs"][i] += 1
            if counts[i] >= args.min_reads:
                rec["detected"][i] += 1
        if counts[FINAL_RUNG] > rec["top_reads"]:
            rec["top_reads"] = counts[FINAL_RUNG]
            rec["top_ref"] = pair["reference_id"]
        if nonov > rec["nonov"]:
            rec["nonov"] = nonov
    return {"cfg": cfg, "per_pair": per_pair, "samples": samples}


# --------------------------------------------------------------------------- #
# aggregation and testing
# --------------------------------------------------------------------------- #
def scope_pairs(ev, ref_set, scope):
    out = []
    for item in ev["per_pair"]:
        p = item["pair"]
        if p["ref_set"] != ref_set:
            continue
        if scope != "ALL" and p["group"] != scope:
            continue
        out.append(item)
    return out


def scope_samples(ev, ref_set, scope):
    recs = ev["samples"].get(ref_set, {})
    if scope == "ALL":
        return dict(recs)
    return dict((k, v) for k, v in recs.items() if v["group"] == scope)


def rpm_of(rec, rung):
    denom = rec["denom"]
    if not denom:
        return None
    return 1e6 * float(rec["reads"][rung]) / float(denom)


def run_tests(ev, ref_set, rung, args, g1, g2):
    """a7's headline tests, recomputed on what survives at one rung."""
    recs = ev["samples"].get(ref_set, {})
    s1 = [v for v in recs.values() if v["group"] == g1]
    s2 = [v for v in recs.values() if v["group"] == g2]

    raw1 = [float(v["reads"][rung]) for v in s1]
    raw2 = [float(v["reads"][rung]) for v in s2]
    rpm1 = [rpm_of(v, rung) for v in s1]
    rpm1 = [v for v in rpm1 if v is not None]
    rpm2 = [rpm_of(v, rung) for v in s2]
    rpm2 = [v for v in rpm2 if v is not None]
    rich1 = [float(v["detected"][rung]) for v in s1]
    rich2 = [float(v["detected"][rung]) for v in s2]

    a = sum(1 for v in s1 if v["detected"][rung] >= args.richness_cut)
    c = sum(1 for v in s2 if v["detected"][rung] >= args.richness_cut)
    any_a = sum(1 for v in s1 if v["detected"][rung] >= 1)
    any_c = sum(1 for v in s2 if v["detected"][rung] >= 1)

    return {
        "n1": len(s1), "n2": len(s2),
        "mw_raw": mann_whitney_u(raw1, raw2),
        "mw_rpm": mann_whitney_u(rpm1, rpm2),
        "mw_rich": mann_whitney_u(rich1, rich2),
        "n_rpm1": len(rpm1), "n_rpm2": len(rpm2),
        "fisher_rich": fisher_exact_2x2(a, len(s1) - a, c, len(s2) - c),
        "rich_a": a, "rich_c": c,
        "fisher_any": fisher_exact_2x2(any_a, len(s1) - any_a,
                                       any_c, len(s2) - any_c),
        "any_a": any_a, "any_c": any_c,
    }


def hotspots_in_window(hotspots, cfg):
    """(total shared hotspots a10 reported, how many the window covers)."""
    if hotspots is None:
        return None, None
    total = len(hotspots)
    if not cfg["on"]:
        return total, 0
    inside = sum(1 for (_rid, pos, _n) in hotspots
                 if cfg["lo"] <= pos <= cfg["hi"])
    return total, inside


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def common_comments(args, cfg, extra=()):
    if cfg["on"]:
        window = ("UTR window excluded: relative %.2f-%.2f of each reference"
                  % (cfg["lo"], cfg["hi"]))
    else:
        window = "UTR window exclusion DISABLED for this table"
    lines = [
        "%s  generated %s" % (SCRIPT, TODAY),
        "no real sample identifiers in this file; ids are anonymous S001..Snnn "
        "(mapping in %s)" % out_name(args.prefix, "utr_exclusion_sample_key.tsv"),
        "QUESTION: does the HIV vs HL anellovirus difference survive when every "
        "read that could come from the conserved UTR is removed?",
        "reads are RE-COUNTED from the existing BAMs; nothing is re-aligned",
        "read filter: samtools view -F %s -q %d (MAPQ is deliberately NOT "
        "filtered: multi-mapping reads are the artefact under test)"
        % (args.exclude_flags, args.min_mapq),
        window,
        "cumulative ladder: 0 all reads | 1 UTR-overlapping reads removed | "
        "2 POS+CIGAR deduplicated | 3 %d-mer entropy >= %.2f bits | "
        "4 AS != XS | 5 >= %d non-overlapping positions per (sample, reference)"
        % (args.entropy_k, args.min_entropy, args.min_distinct_positions),
        "a reference counts as detected at >= %d surviving reads; richness is "
        "the number of detected references; the headline cut is richness >= %d"
        % (args.min_reads, args.richness_cut),
        "normalisation: %s (per million = 1e6 * reads / denominator)"
        % args.norm_source,
        "THIS TESTS ONE MECHANISM: conserved-UTR cross-mapping. Surviving the "
        "filter is NECESSARY BUT NOT SUFFICIENT to call real virus - a "
        "competitive realignment against a fuller reference set is the "
        "confirmatory test and is not attempted here.",
        "ref_set=chimp_flagged is the negative control: no human sample carries "
        "chimpanzee TTV, so whatever survives there is the residual "
        "false-positive rate of this filter, not a detection",
    ]
    lines.extend(extra)
    return lines


def write_sample_key(path, rows):
    lines = [
        "CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
        "generated %s by %s" % (TODAY, SCRIPT),
    ]
    header = ["sample_anon", "sample_real", "group", "runs", "a7_sample_anon"]
    write_tsv(path, lines, header, rows)


def ladder_rows(ev, args, scopes):
    """Rows for the ladder table: rung x reference set x group scope."""
    cfg = ev["cfg"]
    filters = rung_filters(cfg, args)
    rows = []
    for ref_set in REF_SETS:
        if not any(p["pair"]["ref_set"] == ref_set for p in ev["per_pair"]):
            continue
        for scope in scopes:
            items = scope_pairs(ev, ref_set, scope)
            recs = scope_samples(ev, ref_set, scope)
            if not recs and not items:
                continue
            base_reads = sum(i["counts"][0] for i in items)
            base_pairs = sum(1 for i in items if i["counts"][0] > 0)
            prev_reads = None
            for rung in range(N_RUNGS):
                n_reads = sum(i["counts"][rung] for i in items)
                n_pairs = sum(1 for i in items if i["counts"][rung] > 0)
                n_pairs_min = sum(1 for i in items
                                  if i["counts"][rung] >= args.min_reads)
                with_read = sum(1 for v in recs.values() if v["reads"][rung] > 0)
                with_det = sum(1 for v in recs.values()
                               if v["detected"][rung] > 0)
                ge_cut = sum(1 for v in recs.values()
                             if v["detected"][rung] >= args.richness_cut)
                per_sample = [v["reads"][rung] for v in recs.values()]
                row = dict((c, "NA") for c in LADDER_COLUMNS)
                row["row_type"] = "ladder_rung"
                row["ref_set"] = ref_set
                row["scope"] = scope
                row["window_label"] = cfg["label"]
                row["utr_lo"] = fnum(cfg["lo"], 2) if cfg["on"] else "NA"
                row["utr_hi"] = fnum(cfg["hi"], 2) if cfg["on"] else "NA"
                row["utr_filter_applied"] = "1" if cfg["on"] else "0"
                row["rung"] = str(rung)
                row["rung_name"] = RUNG_NAMES[rung]
                row["rung_filter"] = filters[rung]
                row["n_samples_in_scope"] = str(len(recs))
                row["n_samples_with_surviving_read"] = str(with_read)
                row["n_samples_with_detected_reference"] = str(with_det)
                row["n_pairs_with_surviving_read"] = str(n_pairs)
                row["n_pairs_at_min_reads"] = str(n_pairs_min)
                row["n_reads"] = str(n_reads)
                row["reads_dropped_vs_previous_rung"] = (
                    "NA" if prev_reads is None else str(prev_reads - n_reads))
                row["frac_reads_of_rung0"] = fnum(
                    (float(n_reads) / float(base_reads)) if base_reads else None, 4)
                row["frac_pairs_of_rung0"] = fnum(
                    (float(n_pairs) / float(base_pairs)) if base_pairs else None, 4)
                row["median_reads_per_sample"] = fnum(median(per_sample), 2)
                row["max_reads_per_sample"] = str(max(per_sample) if per_sample
                                                  else 0)
                row["min_reads_per_reference"] = str(args.min_reads)
                row["richness_cut"] = str(args.richness_cut)
                row["n_samples_richness_ge_cut"] = str(ge_cut)
                row["frac_samples_richness_ge_cut"] = fnum(
                    (float(ge_cut) / float(len(recs))) if recs else None, 4)
                row["note"] = ("cumulative: rung %d is applied to what rung %d "
                               "left" % (rung, max(0, rung - 1)))
                rows.append([row[c] for c in LADDER_COLUMNS])
                prev_reads = n_reads
    return rows


def write_ladder(path, rows, args, cfg):
    write_tsv(path, common_comments(args, cfg, [
        "one row per (reference set, group scope, ladder rung): the sample, "
        "pair and read counts that survive that rung",
        "scope=ALL pools every group; the per-group rows show what each rung "
        "costs each arm, which is the point of the ladder",
        "n_samples_in_scope counts every successfully streamed sample, "
        "including samples with no anellovirus read at all - they belong in "
        "the denominator of every proportion",
        "CANNOT support: a rung that removes a signal has not proven the "
        "signal false, only that it is compatible with the mechanism that rung "
        "models; and a pair dropped at rung %d is unresolved, not negative"
        % FINAL_RUNG,
    ]), LADDER_COLUMNS, rows)


def by_sample_rows(ev, args, a7_key_of, a7_burden, universe):
    cfg = ev["cfg"]
    human = ev["samples"].get("human_anello", {})
    chimp = ev["samples"].get("chimp_flagged", {})
    rows = []
    for anon in sorted(human):
        rec = human[anon]
        crec = chimp.get(anon)
        a7_anon = a7_key_of.get(anon, "NA")
        ctx = a7_burden.get(a7_anon, {})
        row = dict((c, "NA") for c in BY_SAMPLE_COLUMNS)
        row["sample_anon"] = anon
        row["group"] = rec["group"]
        row["run"] = rec["run"]
        row["window_label"] = cfg["label"]
        row["utr_lo"] = fnum(cfg["lo"], 2) if cfg["on"] else "NA"
        row["utr_hi"] = fnum(cfg["hi"], 2) if cfg["on"] else "NA"
        for i, col in enumerate(RUNG_READ_COLUMNS):
            row[col] = str(rec["reads"][i])
        row["n_refs_with_read_rung0"] = str(rec["nrefs"][0])
        row["n_refs_with_read_final"] = str(rec["nrefs"][FINAL_RUNG])
        row["min_reads_per_reference"] = str(args.min_reads)
        row["richness_rung0"] = str(rec["detected"][0])
        row["richness_final"] = str(rec["detected"][FINAL_RUNG])
        row["richness_cut"] = str(args.richness_cut)
        row["richness_ge_cut_rung0"] = (
            "1" if rec["detected"][0] >= args.richness_cut else "0")
        row["richness_ge_cut_final"] = (
            "1" if rec["detected"][FINAL_RUNG] >= args.richness_cut else "0")
        row["norm_source"] = args.norm_source
        row["norm_denominator"] = ("NA" if not rec["denom"]
                                   else str(int(rec["denom"])))
        row["burden_rpm_rung0"] = fnum(rpm_of(rec, 0), 3)
        row["burden_rpm_final"] = fnum(rpm_of(rec, FINAL_RUNG), 3)
        row["top_reference_final"] = rec["top_ref"]
        row["top_reference_reads_final"] = str(rec["top_reads"])
        row["max_nonoverlapping_positions_final"] = str(rec["nonov"])
        row["chimp_reads_rung0"] = str(crec["reads"][0]) if crec else "0"
        row["chimp_reads_final"] = (str(crec["reads"][FINAL_RUNG]) if crec
                                    else "0")
        row["chimp_richness_final"] = (str(crec["detected"][FINAL_RUNG])
                                       if crec else "0")
        row["a7_sample_anon"] = a7_anon
        row["a7_anello_reads_total"] = ctx.get("reads", "NA")
        row["a7_anello_richness"] = ctx.get("richness", "NA")
        row["note"] = ("a7_* columns are context copied from a7's burden "
                       "table, not recomputed here")
        rows.append([row[c] for c in BY_SAMPLE_COLUMNS])
    return rows


def write_by_sample(path, rows, args, cfg):
    write_tsv(path, common_comments(args, cfg, [
        "one row per successfully streamed sample, at the primary UTR window "
        "only; the sweep table carries the other windows",
        "reads_rung0_all .. reads_rung5_multi_position are cumulative survivor "
        "counts summed over that sample's human anellovirus references",
        "max_nonoverlapping_positions_final is the largest number of mutually "
        "non-overlapping aligned segments on any single reference: 1 means one "
        "pile-up position and nothing else",
        "chimp_* columns are the negative control for the same sample",
        "CANNOT support: a non-zero final burden is not a virus call. It means "
        "only that those reads are not explained by the conserved-UTR "
        "mechanism, duplication, low complexity or an alignment tie",
    ]), BY_SAMPLE_COLUMNS, rows)


def blank_test_row():
    return dict((c, "NA") for c in TEST_COLUMNS)


def build_test_rows(ev, args, g1, g2, a7_burden, pub):
    """(rows, per-rung results) for the BEFORE / AFTER table."""
    cfg = ev["cfg"]
    rows = []

    # ---- BEFORE, as published by a7 (quoted, never recomputed) ------------- #
    (pa, pn1) = pub["g1"]
    (pc, pn2) = pub["g2"]
    row = blank_test_row()
    row["row_type"] = "before_published"
    row["source"] = "a7_as_published"
    row["ref_set"] = "human_anello"
    row["window_label"] = "not_applicable"
    row["rung"] = "NA"
    row["rung_name"] = "a7_published"
    row["metric"] = "richness_ge_%d" % args.richness_cut
    row["test"] = "fisher_exact_2x2"
    row["group1"] = g1
    row["n1"] = str(pn1)
    row["value_group1"] = "%d/%d" % (pa, pn1)
    row["group2"] = g2
    row["n2"] = str(pn2)
    row["value_group2"] = "%d/%d" % (pc, pn2)
    row["statistic"] = "p_as_published"
    row["p_value"] = pub["p"]
    row["note"] = ("quoted from a7's published result, NOT recomputed; a7 "
                   "counted a reference as present at its own --min-reads "
                   "(default 10) on idxstats mapped counts, so it is not "
                   "expected to equal this module's rung 0 exactly")
    rows.append([row[c] for c in TEST_COLUMNS])

    # ---- BEFORE, recomputed from a7's own burden table if it is present ---- #
    if a7_burden:
        b1 = [v for v in a7_burden.values() if v["group"] == g1]
        b2 = [v for v in a7_burden.values() if v["group"] == g2]
        a = sum(1 for v in b1 if v["richness_int"] is not None
                and v["richness_int"] >= args.richness_cut)
        c = sum(1 for v in b2 if v["richness_int"] is not None
                and v["richness_int"] >= args.richness_cut)
        res = fisher_exact_2x2(a, len(b1) - a, c, len(b2) - c)
        row = blank_test_row()
        row["row_type"] = "before_a7_table"
        row["source"] = "a7_burden_table_in_indir"
        row["ref_set"] = "human_anello"
        row["window_label"] = "not_applicable"
        row["rung"] = "NA"
        row["rung_name"] = "a7_table"
        row["metric"] = "richness_ge_%d" % args.richness_cut
        row["test"] = "fisher_exact_2x2"
        row["group1"] = g1
        row["n1"] = str(len(b1))
        row["value_group1"] = "%d/%d" % (a, len(b1))
        row["group2"] = g2
        row["n2"] = str(len(b2))
        row["value_group2"] = "%d/%d" % (c, len(b2))
        row["statistic"] = "fisher_exact_2x2_odds_ratio"
        row["statistic_value"] = odds_text(res)
        row["p_value"] = "NA" if res is None else fp(res["p"])
        row["note"] = ("recomputed here from a7's own anello_richness_human "
                       "column, so it uses a7's read threshold, not this "
                       "module's --min-reads")
        rows.append([row[c] for c in TEST_COLUMNS])

        x = [float(v["reads_int"]) for v in b1 if v["reads_int"] is not None]
        y = [float(v["reads_int"]) for v in b2 if v["reads_int"] is not None]
        mres = mann_whitney_u(x, y)
        row = blank_test_row()
        row["row_type"] = "before_a7_table"
        row["source"] = "a7_burden_table_in_indir"
        row["ref_set"] = "human_anello"
        row["window_label"] = "not_applicable"
        row["rung"] = "NA"
        row["rung_name"] = "a7_table"
        row["metric"] = "burden_reads_raw"
        row["test"] = "mann_whitney_u"
        row["group1"] = g1
        row["n1"] = str(len(x))
        row["group2"] = g2
        row["n2"] = str(len(y))
        if mres is not None:
            row["value_group1"] = fnum(mres["median1"], 4)
            row["value_group2"] = fnum(mres["median2"], 4)
            row["statistic"] = "mann_whitney_U_group1"
            row["statistic_value"] = fnum(mres["U1"], 1)
            row["p_value"] = fp(mres["p"])
            row["effect"] = fnum(mres["effect_r"], 4)
        row["note"] = ("a7's own anello_reads_human_total, medians reported; "
                       "this is the BEFORE burden test")
        rows.append([row[c] for c in TEST_COLUMNS])

    # ---- AFTER, recomputed per rung ---------------------------------------- #
    filters = rung_filters(cfg, args)
    per_rung = {}
    for rung in range(N_RUNGS):
        res = run_tests(ev, "human_anello", rung, args, g1, g2)
        per_rung[rung] = res
        specs = [
            ("burden_reads_raw", "mann_whitney_u", res["mw_raw"],
             res["n1"], res["n2"], None),
            ("burden_reads_per_million", "mann_whitney_u", res["mw_rpm"],
             res["n_rpm1"], res["n_rpm2"],
             "samples without a usable %s denominator are dropped from this "
             "row only" % args.norm_source),
            ("richness_references_detected", "mann_whitney_u", res["mw_rich"],
             res["n1"], res["n2"], None),
        ]
        for metric, test, mres, n1, n2, extra_note in specs:
            row = blank_test_row()
            row["row_type"] = "after_recomputed"
            row["source"] = "a12_rung%d" % rung
            row["ref_set"] = "human_anello"
            row["window_label"] = cfg["label"]
            row["utr_lo"] = fnum(cfg["lo"], 2) if cfg["on"] else "NA"
            row["utr_hi"] = fnum(cfg["hi"], 2) if cfg["on"] else "NA"
            row["rung"] = str(rung)
            row["rung_name"] = RUNG_NAMES[rung]
            row["metric"] = metric
            row["test"] = test
            row["group1"] = g1
            row["n1"] = str(n1)
            row["group2"] = g2
            row["n2"] = str(n2)
            note = "filter at this rung: %s" % filters[rung]
            if mres is None:
                row["p_value"] = "NA"
                note += ("; not tested: %d %s and %d %s sample(s)"
                         % (n1, g1, n2, g2))
            else:
                row["value_group1"] = fnum(mres["median1"], 4)
                row["value_group2"] = fnum(mres["median2"], 4)
                row["statistic"] = "mann_whitney_U_group1"
                row["statistic_value"] = fnum(mres["U1"], 1)
                row["p_value"] = fp(mres["p"])
                row["effect"] = fnum(mres["effect_r"], 4)
                note += ("; two-sided normal approximation, tie- and "
                         "continuity-corrected; values are group medians")
                if min(mres["n1"], mres["n2"]) < 5:
                    note += "; n<5 in one group, p is approximate"
            if extra_note:
                note += "; " + extra_note
            row["note"] = note
            rows.append([row[c] for c in TEST_COLUMNS])

        fspecs = [
            ("richness_ge_%d" % args.richness_cut, res["fisher_rich"],
             res["rich_a"], res["rich_c"],
             "a7's headline cut, recomputed on the survivors"),
            ("any_reference_detected", res["fisher_any"],
             res["any_a"], res["any_c"],
             "at least one reference with >= %d surviving reads"
             % args.min_reads),
        ]
        for metric, fres, a, c, why in fspecs:
            row = blank_test_row()
            row["row_type"] = "after_recomputed"
            row["source"] = "a12_rung%d" % rung
            row["ref_set"] = "human_anello"
            row["window_label"] = cfg["label"]
            row["utr_lo"] = fnum(cfg["lo"], 2) if cfg["on"] else "NA"
            row["utr_hi"] = fnum(cfg["hi"], 2) if cfg["on"] else "NA"
            row["rung"] = str(rung)
            row["rung_name"] = RUNG_NAMES[rung]
            row["metric"] = metric
            row["test"] = "fisher_exact_2x2"
            row["group1"] = g1
            row["n1"] = str(res["n1"])
            row["value_group1"] = "%d/%d" % (a, res["n1"])
            row["group2"] = g2
            row["n2"] = str(res["n2"])
            row["value_group2"] = "%d/%d" % (c, res["n2"])
            row["statistic"] = "fisher_exact_2x2_odds_ratio"
            row["statistic_value"] = odds_text(fres)
            note = "%s; filter at this rung: %s" % (why, filters[rung])
            if fres is None:
                row["p_value"] = "NA"
                note += ("; not tested: degenerate 2x2 table "
                         "[[%d,%d],[%d,%d]]" % (a, res["n1"] - a, c,
                                                res["n2"] - c))
            else:
                row["p_value"] = fp(fres["p"])
                note += ("; two-sided Fisher exact on [[%d,%d],[%d,%d]], "
                         "samples counted once each"
                         % (a, res["n1"] - a, c, res["n2"] - c))
            row["note"] = note
            rows.append([row[c] for c in TEST_COLUMNS])

    # ---- the chimpanzee residual, summarised, never tested ----------------- #
    chimp = ev["samples"].get("chimp_flagged", {})
    if chimp:
        for rung in (0, FINAL_RUNG):
            c1 = [v for v in chimp.values() if v["group"] == g1]
            c2 = [v for v in chimp.values() if v["group"] == g2]
            row = blank_test_row()
            row["row_type"] = "negative_control"
            row["source"] = "a12_rung%d" % rung
            row["ref_set"] = "chimp_flagged"
            row["window_label"] = cfg["label"]
            row["utr_lo"] = fnum(cfg["lo"], 2) if cfg["on"] else "NA"
            row["utr_hi"] = fnum(cfg["hi"], 2) if cfg["on"] else "NA"
            row["rung"] = str(rung)
            row["rung_name"] = RUNG_NAMES[rung]
            row["metric"] = "chimp_reads_surviving"
            row["test"] = "none_summary_only"
            row["group1"] = g1
            row["n1"] = str(len(c1))
            row["value_group1"] = str(sum(v["reads"][rung] for v in c1))
            row["group2"] = g2
            row["n2"] = str(len(c2))
            row["value_group2"] = str(sum(v["reads"][rung] for v in c2))
            row["note"] = ("no human sample carries chimpanzee TTV, so every "
                           "surviving read here is a false positive of this "
                           "filter; summarised, never tested")
            rows.append([row[c] for c in TEST_COLUMNS])
    return rows, per_rung


def write_group_test(path, rows, args, cfg, g1, g2):
    write_tsv(path, common_comments(args, cfg, [
        "the BEFORE / AFTER block for %s vs %s, one place to read the change"
        % (g1, g2),
        "row_type=before_published: a7's published numbers, QUOTED not "
        "recomputed, so the change is explicit and quotable",
        "row_type=before_a7_table: the same contrast recomputed from a7's own "
        "burden table when it is in --indir (a7's read threshold, not this "
        "module's)",
        "row_type=after_recomputed: this module's tests at every rung of the "
        "ladder, so the reader sees what each filter step costs",
        "row_type=negative_control: chimpanzee references, summarised only",
        "n1 and n2 sit beside every p value; both tests are the same "
        "standard-library implementations a7 and a10 use",
        "CANNOT support: a p value that survives is not proof of virus, and a "
        "p value that dies is proof only that the contrast is compatible with "
        "conserved-UTR cross-mapping",
    ]), TEST_COLUMNS, rows)


def sweep_rows(evs, args, g1, g2, hotspots):
    rows = []
    for ev in evs:
        cfg = ev["cfg"]
        total, inside = hotspots_in_window(hotspots, cfg)
        filters = rung_filters(cfg, args)
        human_items = scope_pairs(ev, "human_anello", "ALL")
        base_reads = sum(i["counts"][0] for i in human_items)
        chimp_items = scope_pairs(ev, "chimp_flagged", "ALL")
        chimp_recs = ev["samples"].get("chimp_flagged", {})
        human_recs = ev["samples"].get("human_anello", {})
        for rung in range(N_RUNGS):
            res = run_tests(ev, "human_anello", rung, args, g1, g2)
            n_reads = sum(i["counts"][rung] for i in human_items)
            row = dict((c, "NA") for c in SWEEP_COLUMNS)
            row["row_type"] = "window_rung"
            row["window_label"] = cfg["label"]
            row["utr_lo"] = fnum(cfg["lo"], 2) if cfg["on"] else "NA"
            row["utr_hi"] = fnum(cfg["hi"], 2) if cfg["on"] else "NA"
            row["utr_filter_applied"] = "1" if cfg["on"] else "0"
            row["is_primary_window"] = "1" if cfg["primary"] else "0"
            row["ref_set"] = "human_anello"
            row["rung"] = str(rung)
            row["rung_name"] = RUNG_NAMES[rung]
            row["n_samples_in_scope"] = str(len(human_recs))
            row["n_samples_with_surviving_read"] = str(
                sum(1 for v in human_recs.values() if v["reads"][rung] > 0))
            row["n_pairs_with_surviving_read"] = str(
                sum(1 for i in human_items if i["counts"][rung] > 0))
            row["n_pairs_at_min_reads"] = str(
                sum(1 for i in human_items
                    if i["counts"][rung] >= args.min_reads))
            row["n_reads"] = str(n_reads)
            row["frac_reads_of_rung0"] = fnum(
                (float(n_reads) / float(base_reads)) if base_reads else None, 4)
            row["richness_cut"] = str(args.richness_cut)
            row["min_reads_per_reference"] = str(args.min_reads)
            row["group1"] = g1
            row["n1"] = str(res["n1"])
            row["richness_ge_cut_group1"] = "%d/%d" % (res["rich_a"], res["n1"])
            row["group2"] = g2
            row["n2"] = str(res["n2"])
            row["richness_ge_cut_group2"] = "%d/%d" % (res["rich_c"], res["n2"])
            row["fisher_odds_ratio"] = odds_text(res["fisher_rich"])
            row["fisher_p"] = ("NA" if res["fisher_rich"] is None
                               else fp(res["fisher_rich"]["p"]))
            if res["mw_raw"] is not None:
                row["burden_median_group1"] = fnum(res["mw_raw"]["median1"], 3)
                row["burden_median_group2"] = fnum(res["mw_raw"]["median2"], 3)
                row["mann_whitney_U_group1"] = fnum(res["mw_raw"]["U1"], 1)
                row["mann_whitney_p"] = fp(res["mw_raw"]["p"])
            row["chimp_reads_surviving"] = str(
                sum(i["counts"][rung] for i in chimp_items))
            row["chimp_samples_with_read"] = str(
                sum(1 for v in chimp_recs.values() if v["reads"][rung] > 0))
            row["chimp_pairs_at_min_reads"] = str(
                sum(1 for i in chimp_items
                    if i["counts"][rung] >= args.min_reads))
            row["a10_shared_hotspots_total"] = ("NA" if total is None
                                                else str(total))
            row["a10_shared_hotspots_inside_window"] = ("NA" if inside is None
                                                        else str(inside))
            note = "filter at this rung: %s" % filters[rung]
            if rung == 0:
                note += ("; rung 0 does not depend on the window, so it is "
                         "identical in every window block, by construction")
            if not cfg["on"] and rung == 1:
                note += "; UTR rung disabled, so rung 1 equals rung 0"
            row["note"] = note
            rows.append([row[c] for c in SWEEP_COLUMNS])
    return rows


def write_sweep(path, rows, args, cfg, g1, g2):
    write_tsv(path, common_comments(args, cfg, [
        "the whole ladder repeated at every UTR window in --sweep-lo and once "
        "with the UTR rung disabled, so the conclusion cannot be an artefact "
        "of where the window boundary was drawn",
        "one row per (window, rung); the %s vs %s Fisher and Mann-Whitney "
        "columns are the same tests as in %s" % (
            g1, g2, out_name(args.prefix, "utr_exclusion_group_test.tsv")),
        "chimp_* columns are the negative control under the same window: "
        "surviving chimpanzee reads are the residual false-positive rate",
        "a10_shared_hotspots_inside_window counts how many of the shared "
        "hotspots a10 reported (ALL scope, >= 2 agreeing pairs) fall inside "
        "this window; a window that covers none is not testing a10's finding",
        "CANNOT support: the sweep shows robustness to the boundary, not that "
        "the excluded region really is the UTR. Check the coordinate against "
        "the reference annotation before naming it",
    ]), SWEEP_COLUMNS, rows)


def write_headers_only(paths, args, cfg, g1, g2):
    write_ladder(paths["ladder"], [], args, cfg)
    write_by_sample(paths["by_sample"], [], args, cfg)
    write_group_test(paths["test"], [], args, cfg, g1, g2)
    write_sweep(paths["sweep"], [], args, cfg, g1, g2)


# --------------------------------------------------------------------------- #
# stdout helpers
# --------------------------------------------------------------------------- #
def print_ladder_table(ev, args, ref_set, g1, g2, title, legend=True):
    items_all = scope_pairs(ev, ref_set, "ALL")
    if not items_all:
        return
    filters = rung_filters(ev["cfg"], args)
    recs1 = scope_samples(ev, ref_set, g1)
    recs2 = scope_samples(ev, ref_set, g2)
    base = sum(i["counts"][0] for i in items_all)
    print("")
    print("-- %s --" % title)
    print("   %-4s %-24s %9s %6s %8s %9s %9s %9s"
          % ("rung", "filter", "reads", "kept", "pairs>=%d" % args.min_reads,
             "samples", "%s r>=%d" % (g1, args.richness_cut),
             "%s r>=%d" % (g2, args.richness_cut)))
    for rung in range(N_RUNGS):
        n_reads = sum(i["counts"][rung] for i in items_all)
        n_pairs_min = sum(1 for i in items_all
                          if i["counts"][rung] >= args.min_reads)
        recs = scope_samples(ev, ref_set, "ALL")
        with_read = sum(1 for v in recs.values() if v["reads"][rung] > 0)
        a = sum(1 for v in recs1.values()
                if v["detected"][rung] >= args.richness_cut)
        c = sum(1 for v in recs2.values()
                if v["detected"][rung] >= args.richness_cut)
        kept = fnum((float(n_reads) / float(base)) if base else None, 3)
        print("   %-4d %-24s %9d %6s %8d %9s %9s %9s"
              % (rung, RUNG_NAMES[rung][:24], n_reads, kept, n_pairs_min,
                 "%d/%d" % (with_read, len(recs)),
                 "%d/%d" % (a, len(recs1)), "%d/%d" % (c, len(recs2))))
    if legend:
        for rung in range(N_RUNGS):
            print("   rung %d = %s" % (rung, filters[rung]))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Does the HIV vs HL anellovirus difference survive when "
                    "every read that could come from the conserved UTR is "
                    "removed? Re-counts the existing BAMs under a cumulative "
                    "evidence filter and re-runs a7's headline tests. Standard "
                    "library only, no realignment, no figures.",
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
                   help="directory holding the a7 and a10 outputs used for "
                        "context and for the BEFORE row")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR,
                   help="output directory (created if absent)")
    p.add_argument("--prefix", default=DEFAULT_PREFIX,
                   help="filename prefix for every output file")
    p.add_argument("--a7-burden", default=None,
                   help="explicit path to a7_virome_anellovirus_burden.tsv")
    p.add_argument("--a7-key", default=None,
                   help="explicit path to a7_virome_sample_key.tsv")
    p.add_argument("--a10-pooled", default=None,
                   help="explicit path to anello_read_audit_pooled_positions.tsv")
    p.add_argument("--anello-accessions", default=DEFAULT_ANELLO_ACC_FILE,
                   help="file of anellovirus accessions / reference ids")
    p.add_argument("--chimp-accessions", default=DEFAULT_CHIMP_ACC,
                   help="comma-separated chimpanzee-isolate accessions; these "
                        "are the negative control")
    p.add_argument("--include-chimp", action="store_true",
                   help="fold the chimpanzee references into the human metrics "
                        "instead of keeping them as the negative control "
                        "(destroys the control; not recommended)")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--samtools-threads", type=int, default=1,
                   help="samtools -@ value")
    p.add_argument("--bam-glob", default="bam/*.bam",
                   help="glob for BAMs relative to a run directory")
    p.add_argument("--exclude-flags", default="0x904",
                   help="samtools view -F value (unmapped, secondary, "
                        "supplementary); duplicates are kept, rung 2 handles them")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="samtools -q value; 0 on purpose, a MAPQ filter would "
                        "pre-empt the cross-mapping this module measures")
    p.add_argument("--utr-lo", type=float, default=0.75,
                   help="lower edge of the excluded UTR window, relative "
                        "position along each reference")
    p.add_argument("--utr-hi", type=float, default=1.00,
                   help="upper edge of the excluded UTR window")
    p.add_argument("--sweep-lo", default=DEFAULT_SWEEP_LO,
                   help="comma-separated --utr-lo values for the window sweep; "
                        "the primary --utr-lo is added if absent and the "
                        "filter-disabled case is always included")
    p.add_argument("--min-entropy", type=float, default=1.2,
                   help="k-mer Shannon entropy in bits below which a read is "
                        "called low complexity and dropped")
    p.add_argument("--entropy-k", type=int, default=3,
                   help="k for the entropy filter (the spec value is 3)")
    p.add_argument("--min-distinct-positions", type=int, default=2,
                   help="non-overlapping aligned segments required per "
                        "(sample, reference); 1 pile-up position is not evidence")
    p.add_argument("--min-reads", type=int, default=5,
                   help="surviving reads needed for a reference to count as "
                        "detected in a sample")
    p.add_argument("--richness-cut", type=int, default=3,
                   help="richness at or above which a sample is a Fisher "
                        "success (a7's headline cut)")
    p.add_argument("--norm-source",
                   choices=["idxstats_mapped", "idxstats_total",
                            "filtered_categories"],
                   default="idxstats_mapped",
                   help="denominator for the per-million normalisation "
                        "(a7's default is idxstats_mapped)")
    p.add_argument("--max-ref-len", type=int, default=100000,
                   help="skip any 'anellovirus' reference longer than this. "
                        "Anellovirus genomes are ~3.7 kb, so a longer hit is a "
                        "misclassified reference")
    p.add_argument("--a7-published-group1", default=DEFAULT_A7_PUB_G1,
                   help="a7's published successes/total for group 1, quoted in "
                        "the BEFORE row")
    p.add_argument("--a7-published-group2", default=DEFAULT_A7_PUB_G2,
                   help="a7's published successes/total for group 2")
    p.add_argument("--a7-published-p", default=DEFAULT_A7_PUB_P,
                   help="a7's published p value for that contrast")
    p.add_argument("--test-groups", default="HIV,HL",
                   help="the two group labels compared")
    p.add_argument("--limit", type=int, default=0,
                   help="stream only the first N BAMs per run (debugging)")
    p.add_argument("--no-run-name-group", action="store_true",
                   help="do not fall back to the run directory name for TCL")
    return p


def resolve_args(args):
    """Clamp the options that would otherwise produce nonsense silently."""
    if args.entropy_k < 1:
        print("WARN: --entropy-k must be >= 1, got %d; using 3" % args.entropy_k)
        args.entropy_k = 3
    if args.min_reads < 1:
        args.min_reads = 1
    if args.richness_cut < 1:
        args.richness_cut = 1
    if args.min_distinct_positions < 1:
        args.min_distinct_positions = 1
    if not (0.0 <= args.utr_lo < args.utr_hi <= 1.0):
        print("WARN: --utr-lo %.3f / --utr-hi %.3f is not a valid relative "
              "window; using 0.75 / 1.00" % (args.utr_lo, args.utr_hi))
        args.utr_lo, args.utr_hi = 0.75, 1.00
    return args


def build_configs(args):
    """The primary window first, then the sweep windows, then 'disabled'."""
    los = []
    for token in (args.sweep_lo or "").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            val = float(token)
        except ValueError:
            print("WARN: --sweep-lo entry %r is not a number, ignored"
                  % to_ascii(token))
            continue
        if not (0.0 <= val < args.utr_hi):
            print("WARN: --sweep-lo entry %.3f is not below --utr-hi %.3f, "
                  "ignored" % (val, args.utr_hi))
            continue
        los.append(round(val, 6))
    if round(args.utr_lo, 6) not in los:
        los.append(round(args.utr_lo, 6))
    los = sorted(set(los))
    cfgs = []
    for lo in los:
        cfgs.append({
            "label": "%.2f-%.2f" % (lo, args.utr_hi),
            "lo": lo, "hi": args.utr_hi, "on": True,
            "primary": abs(lo - args.utr_lo) < 1e-9,
        })
    cfgs.append({"label": "disabled", "lo": args.utr_lo, "hi": args.utr_hi,
                 "on": False, "primary": False})
    return cfgs


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = resolve_args(build_parser().parse_args(argv))

    test_groups = [g.strip() for g in args.test_groups.split(",") if g.strip()]
    if len(test_groups) != 2:
        print("WARN: --test-groups needs exactly two labels, got %r; using HIV,HL"
              % args.test_groups)
        test_groups = ["HIV", "HL"]
    g1, g2 = test_groups

    pub = {
        "g1": parse_fraction_pair(args.a7_published_group1, (18, 37)),
        "g2": parse_fraction_pair(args.a7_published_group2, (0, 23)),
        "p": to_ascii(args.a7_published_p) or DEFAULT_A7_PUB_P,
    }

    runs = args.run if args.run else list(DEFAULT_RUNS)
    runs = [r if os.path.isabs(r) else os.path.join(args.runs_root, r)
            for r in runs]

    cfgs = build_configs(args)
    primary_cfg = next((c for c in cfgs if c["primary"]), cfgs[0])

    print("")
    print("%s  %s" % (SCRIPT, TODAY))
    print("question does the %s vs %s anellovirus difference survive removing "
          "every UTR-compatible read?" % (g1, g2))
    print("runs     %s" % ", ".join(os.path.basename(r.rstrip("/\\")) or r
                                    for r in runs))
    print("filter   samtools view -F %s -q %d (MAPQ not filtered: "
          "multi-mapping is the artefact under test)"
          % (args.exclude_flags, args.min_mapq))
    print("ladder   UTR %.2f-%.2f excluded | POS+CIGAR dedup | %d-mer entropy "
          ">= %.2f | AS != XS | >= %d non-overlapping positions"
          % (args.utr_lo, args.utr_hi, args.entropy_k, args.min_entropy,
             args.min_distinct_positions))
    print("calls    reference detected at >= %d surviving reads; headline cut "
          "richness >= %d; norm %s"
          % (args.min_reads, args.richness_cut, args.norm_source))
    print("sweep    %s" % ", ".join(c["label"] for c in cfgs))
    print("")

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            warn_missing("writable output directory", args.outdir)
            return 0
    paths = {
        "ladder": os.path.join(args.outdir,
                               out_name(args.prefix, "utr_exclusion_ladder.tsv")),
        "by_sample": os.path.join(
            args.outdir, out_name(args.prefix, "utr_exclusion_by_sample.tsv")),
        "test": os.path.join(
            args.outdir, out_name(args.prefix, "utr_exclusion_group_test.tsv")),
        "sweep": os.path.join(args.outdir,
                              out_name(args.prefix, "utr_window_sweep.tsv")),
        "key": os.path.join(
            args.outdir, out_name(args.prefix, "utr_exclusion_sample_key.tsv")),
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

    a7_key = load_a7_key(args)
    a7_burden = load_a7_burden(args)
    hotspots = load_a10_hotspots(args)

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
    a7_key_of_anon = {}
    for name in sorted(real_names):
        key_rows.append([anon_of[name], name, group_of_sample.get(name, "NA"),
                         ";".join(runs_of.get(name, [])),
                         a7_key.get(name, "NA")])
        a7_key_of_anon[anon_of[name]] = a7_key.get(name, "NA")
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
        write_headers_only(paths, args, primary_cfg, g1, g2)
        print("wrote: %s, %s, %s, %s" % (paths["ladder"], paths["by_sample"],
                                         paths["test"], paths["sweep"]))
        return 0

    ver = check_samtools(args.samtools)
    if ver is None:
        print("Nothing can be re-counted without samtools; writing headed, "
              "empty tables.")
        write_headers_only(paths, args, primary_cfg, g1, g2)
        print("wrote: %s, %s, %s, %s" % (paths["ladder"], paths["by_sample"],
                                         paths["test"], paths["sweep"]))
        return 0
    print("samtools %s" % ver)
    print("")

    # ---- pass 2: stream the BAMs once, keeping per-read facts -------------- #
    pairs = []
    universe = {}
    oversize = set()
    n_no_index = 0
    n_failed = 0
    n_no_anello_ref = 0
    n_no_denom = 0
    counters = {"no_seq": 0, "no_as_tag": 0, "no_cigar": 0}
    for run_dir, base, bams in run_bams:
        cat_totals = None
        if args.norm_source == "filtered_categories":
            cat_path = find_filtered_category_counts(
                os.path.join(run_dir, "results"))
            if not cat_path:
                warn_missing("category count table for run " + base,
                             os.path.join(run_dir, "results",
                                          "*_category_counts.tsv"))
            else:
                cat_totals = parse_category_counts(cat_path)
                if cat_totals is None:
                    warn_missing("readable category count table for run " + base,
                                 cat_path)
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
            lens, total_mapped, total_unmapped = idxstats_table(run_dir, sample)
            from_idxstats = bool(lens)
            if not lens:
                lens = header_lengths(args.samtools, bam_path)
            if not lens:
                warn_missing("idxstats and BAM header for " + anon,
                             os.path.join(run_dir, "results",
                                          "<sample>.idxstats.tsv"))
                continue
            if args.norm_source == "idxstats_total":
                denom = total_mapped + total_unmapped
            elif args.norm_source == "filtered_categories":
                denom = (cat_totals or {}).get(sample)
            else:
                denom = total_mapped
            if not denom or denom <= 0:
                denom = None
                n_no_denom += 1

            wanted = {}
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
                wanted[rname] = (seqlen, mapped)
            if not wanted:
                # The sample still belongs to the cohort: register it with zero
                # counts so it stays in every denominator.
                n_no_anello_ref += 1
                universe[anon] = {"group": group, "run": base, "denom": denom}
                done += 1
                continue
            ref_lens = dict((r, v[0]) for r, v in wanted.items())
            reads_by_ref, err, cnt = stream_sample(
                args.samtools, bam_path, list(wanted), ref_lens, args,
                sample, anon)
            if reads_by_ref is None:
                n_failed += 1
                if n_failed == 1:
                    # err was already masked inside stream_sample
                    print("WARN: samtools view failed for %s (%s); that sample "
                          "is skipped" % (anon, err))
                continue
            for key in counters:
                counters[key] += cnt.get(key, 0)
            done += 1
            universe[anon] = {"group": group, "run": base, "denom": denom}
            for rname, (seqlen, mapped) in sorted(wanted.items()):
                chimp = is_chimp_ref(rname, refmap, chimp_norm, chimp_ids)
                ref_set = ("human_anello"
                           if (args.include_chimp or not chimp)
                           else "chimp_flagged")
                pairs.append({
                    "run": base,
                    "sample_anon": anon,
                    "group": group,
                    "reference_id": rname,
                    "ref_label": ref_label(rname, refmap),
                    "ref_len": seqlen,
                    "chimp": chimp,
                    "ref_set": ref_set,
                    "idxstats_mapped": mapped,
                    "reads": reads_by_ref.get(rname, []),
                })
        print("  %-46s %3d/%3d samples streamed" % (base[:46], done, len(bams)))

    if n_no_index:
        print("WARN: %d BAM(s) without an index were skipped" % n_no_index)
    if n_failed:
        print("WARN: %d sample(s) were skipped because samtools view failed"
              % n_failed)
    if n_no_anello_ref:
        print("NOTE: %d sample(s) had no anellovirus reference with a mapped "
              "read; they are kept with zero counts so they stay in the "
              "denominators" % n_no_anello_ref)
    if n_no_denom:
        print("NOTE: %d sample(s) have no usable %s denominator; they are "
              "dropped from the per-million rows only"
              % (n_no_denom, args.norm_source))
    if counters["no_seq"]:
        print("NOTE: %d read(s) carried no SEQ and could not be scored for "
              "entropy; they are KEPT by rung 3" % counters["no_seq"])
    if counters["no_as_tag"]:
        print("NOTE: %d read(s) carried no AS tag and could not be scored for "
              "ambiguity; they are KEPT by rung 4" % counters["no_as_tag"])
    if counters["no_cigar"]:
        print("NOTE: %d read(s) carried no usable CIGAR; their reference span "
              "was taken from the SEQ length" % counters["no_cigar"])
    if args.include_chimp:
        print("WARN: --include-chimp folds the chimpanzee references into the "
              "human metrics, so this run has NO negative control")

    if not universe:
        print("")
        print("No sample could be streamed; writing headed, empty tables.")
        write_headers_only(paths, args, primary_cfg, g1, g2)
        print("wrote: %s, %s, %s, %s" % (paths["ladder"], paths["by_sample"],
                                         paths["test"], paths["sweep"]))
        return 0

    # ---- pass 3: apply the ladder under every window ----------------------- #
    evs = [evaluate_config(pairs, cfg, universe, args) for cfg in cfgs]
    primary_ev = next((e for e in evs if e["cfg"]["primary"]), evs[0])

    scopes = ["ALL"]
    seen_groups = set(v["group"] for v in universe.values())
    scopes += [g for g in GROUP_ORDER if g in seen_groups]
    scopes += [g for g in sorted(seen_groups) if g not in scopes]

    write_ladder(paths["ladder"], ladder_rows(primary_ev, args, scopes), args,
                 primary_cfg)
    write_by_sample(paths["by_sample"],
                    by_sample_rows(primary_ev, args, a7_key_of_anon, a7_burden,
                                   universe), args, primary_cfg)
    test_rows, per_rung = build_test_rows(primary_ev, args, g1, g2, a7_burden,
                                          pub)
    write_group_test(paths["test"], test_rows, args, primary_cfg, g1, g2)
    write_sweep(paths["sweep"], sweep_rows(evs, args, g1, g2, hotspots), args,
                primary_cfg, g1, g2)

    # ---- stdout summary ---------------------------------------------------- #
    n_human_pairs = sum(1 for p in pairs if p["ref_set"] == "human_anello")
    n_chimp_pairs = sum(1 for p in pairs if p["ref_set"] == "chimp_flagged")
    print("")
    print("-- input --")
    print("   %d sample(s) streamed, %d (sample, reference) pair(s): %d human "
          "anellovirus + %d chimpanzee control"
          % (len(universe), len(pairs), n_human_pairs, n_chimp_pairs))
    group_counts = {}
    for v in universe.values():
        group_counts[v["group"]] = group_counts.get(v["group"], 0) + 1
    print("   groups: " + ", ".join("%s=%d" % (g, group_counts[g])
                                    for g in sorted(group_counts)))

    print_ladder_table(primary_ev, args, "human_anello", g1, g2,
                       "filter ladder, human anellovirus references (UTR %s)"
                       % primary_cfg["label"])
    print_ladder_table(primary_ev, args, "chimp_flagged", g1, g2,
                       "NEGATIVE CONTROL: the same ladder on chimpanzee "
                       "references (no human carries chimp TTV)",
                       legend=False)

    # ---- the BEFORE / AFTER block ------------------------------------------ #
    r0 = per_rung[0]
    rf = per_rung[FINAL_RUNG]
    print("")
    print("-- BEFORE / AFTER: richness >= %d, %s vs %s --"
          % (args.richness_cut, g1, g2))
    print("   %-42s %-9s %-9s %-9s %-11s %s"
          % ("source", g1, g2, "OR", "p", "n"))
    print("   %-42s %-9s %-9s %-9s %-11s %s"
          % ("a7 as published (quoted, not recomputed)",
             "%d/%d" % pub["g1"], "%d/%d" % pub["g2"], "NA", pub["p"],
             pub["g1"][1] + pub["g2"][1]))
    if a7_burden:
        b1 = [v for v in a7_burden.values() if v["group"] == g1]
        b2 = [v for v in a7_burden.values() if v["group"] == g2]
        a = sum(1 for v in b1 if v["richness_int"] is not None
                and v["richness_int"] >= args.richness_cut)
        c = sum(1 for v in b2 if v["richness_int"] is not None
                and v["richness_int"] >= args.richness_cut)
        fres = fisher_exact_2x2(a, len(b1) - a, c, len(b2) - c)
        print("   %-42s %-9s %-9s %-9s %-11s %s"
              % ("a7 burden table in --indir (a7 threshold)",
                 "%d/%d" % (a, len(b1)), "%d/%d" % (c, len(b2)),
                 odds_text(fres), "NA" if fres is None else fp(fres["p"]),
                 len(b1) + len(b2)))
    else:
        print("   %-42s (a7 burden table not found in --indir)"
              % "a7 burden table in --indir")
    print("   %-42s %-9s %-9s %-9s %-11s %s"
          % ("a12 rung 0  all reads on anello refs",
             "%d/%d" % (r0["rich_a"], r0["n1"]),
             "%d/%d" % (r0["rich_c"], r0["n2"]),
             odds_text(r0["fisher_rich"]),
             "NA" if r0["fisher_rich"] is None else fp(r0["fisher_rich"]["p"]),
             r0["n1"] + r0["n2"]))
    print("   %-42s %-9s %-9s %-9s %-11s %s"
          % ("a12 rung %d  full ladder, UTR %s"
             % (FINAL_RUNG, primary_cfg["label"]),
             "%d/%d" % (rf["rich_a"], rf["n1"]),
             "%d/%d" % (rf["rich_c"], rf["n2"]),
             odds_text(rf["fisher_rich"]),
             "NA" if rf["fisher_rich"] is None else fp(rf["fisher_rich"]["p"]),
             rf["n1"] + rf["n2"]))
    print("   a reference counts as detected at >= %d surviving reads; a7 used "
          "its own threshold," % args.min_reads)
    print("   so the a7 rows and the a12 rung 0 row are not required to agree "
          "exactly.")

    print("")
    print("-- BEFORE / AFTER: burden (reads per sample), %s vs %s --"
          % (g1, g2))
    print("   %-42s %-11s %-11s %-9s %-11s %s"
          % ("source", "%s median" % g1, "%s median" % g2, "U", "p", "n"))
    if a7_burden:
        x = [float(v["reads_int"]) for v in a7_burden.values()
             if v["group"] == g1 and v["reads_int"] is not None]
        y = [float(v["reads_int"]) for v in a7_burden.values()
             if v["group"] == g2 and v["reads_int"] is not None]
        mres = mann_whitney_u(x, y)
        if mres is not None:
            print("   %-42s %-11s %-11s %-9s %-11s %s"
                  % ("a7 burden table in --indir",
                     fnum(mres["median1"], 2), fnum(mres["median2"], 2),
                     fnum(mres["U1"], 1), fp(mres["p"]), len(x) + len(y)))
    for label, res in (("a12 rung 0  all reads on anello refs", r0),
                       ("a12 rung %d  full ladder, UTR %s"
                        % (FINAL_RUNG, primary_cfg["label"]), rf)):
        mres = res["mw_raw"]
        if mres is None:
            print("   %-42s not tested (n %s=%d, %s=%d)"
                  % (label, g1, res["n1"], g2, res["n2"]))
            continue
        print("   %-42s %-11s %-11s %-9s %-11s %s"
              % (label, fnum(mres["median1"], 2), fnum(mres["median2"], 2),
                 fnum(mres["U1"], 1), fp(mres["p"]), res["n1"] + res["n2"]))

    # ---- the sweep --------------------------------------------------------- #
    print("")
    print("-- UTR window sweep (full ladder at rung %d in each window) --"
          % FINAL_RUNG)
    print("   %-11s %8s %8s %-9s %-9s %-11s %5s %-11s %8s %s"
          % ("window", "reads", "pairs", "%s r>=%d" % (g1, args.richness_cut),
             "%s r>=%d" % (g2, args.richness_cut), "fisher p", "n",
             "burden p", "chimp", "a10 hs in win"))
    for ev in evs:
        cfg = ev["cfg"]
        res = run_tests(ev, "human_anello", FINAL_RUNG, args, g1, g2)
        items = scope_pairs(ev, "human_anello", "ALL")
        chimp_items = scope_pairs(ev, "chimp_flagged", "ALL")
        total, inside = hotspots_in_window(hotspots, cfg)
        print("   %-11s %8d %8d %-9s %-9s %-11s %5d %-11s %8d %s"
              % (cfg["label"],
                 sum(i["counts"][FINAL_RUNG] for i in items),
                 sum(1 for i in items
                     if i["counts"][FINAL_RUNG] >= args.min_reads),
                 "%d/%d" % (res["rich_a"], res["n1"]),
                 "%d/%d" % (res["rich_c"], res["n2"]),
                 ("NA" if res["fisher_rich"] is None
                  else fp(res["fisher_rich"]["p"])),
                 res["n1"] + res["n2"],
                 ("NA" if res["mw_raw"] is None else fp(res["mw_raw"]["p"])),
                 sum(i["counts"][FINAL_RUNG] for i in chimp_items),
                 ("NA" if total is None else "%s/%s" % (inside, total))))
    print("   pairs = (sample, reference) pairs with >= %d surviving reads; "
          "chimp = surviving reads on the" % args.min_reads)
    print("   chimpanzee control references, which is the residual "
          "false-positive rate of this filter.")
    if hotspots is None:
        print("   a10 hs in win = NA: a10's pooled-position table was not "
              "found in --indir, so the windows")
        print("   could not be checked against the hotspots a10 actually "
              "observed.")
    else:
        print("   a10 hs in win = how many of a10's %d shared hotspot(s) that "
              "window covers; a window that" % len(hotspots))
        print("   covers none is not testing a10's finding.")

    # ---- the caveat block -------------------------------------------------- #
    print("")
    print("READ THIS BEFORE QUOTING ANY OF IT:")
    print("   This tests ONE mechanism: cross-mapping onto the conserved "
          "terminal UTR (plus duplication,")
    print("   low complexity and alignment ties). SURVIVING THE FILTER IS "
          "NECESSARY BUT NOT SUFFICIENT")
    print("   to call real virus. Reads outside the UTR can still be "
          "cross-mapped from an anellovirus")
    print("   the panel does not carry, from another small circular DNA virus, "
          "or from an unmasked human")
    print("   repeat. The confirmatory test is a COMPETITIVE REALIGNMENT of the "
          "survivors against a fuller")
    print("   reference set; this module deliberately does not realign "
          "anything.")
    print("   Rung 1 is conservative against detection by construction: a true "
          "infection represented in")
    print("   the panel only through its conserved UTR is removed too. The "
          "sweep exists to show by how")
    print("   much, and the chimpanzee arm to show what the filter still lets "
          "through.")

    print("")
    print("wrote:")
    for key in ("ladder", "by_sample", "test", "sweep"):
        print("  %s" % paths[key])
    if key_rows:
        print("  %s" % paths["key"])
        print("REMINDER: %s contains real sample identifiers - do not commit or "
              "email it." % os.path.basename(paths["key"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
