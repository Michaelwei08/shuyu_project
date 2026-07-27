#!/usr/bin/env python3
"""a9_cd4_correlation.py -- anellovirus burden against CD4 count (dose-response).

WHAT IT COMPUTES
  a7 reports that anellovirus burden is higher in HIV+ than in HL. That is a
  two-group contrast. The TTV / immunocompetence literature is written in the
  form "burden rises as CD4 falls", so this module re-expresses the same burden
  numbers against a user-supplied CD4 count.

  It reads the per-sample anellovirus burden table a7 already wrote
  (a7_virome_anellovirus_burden.tsv, anonymous ids), joins a user-supplied CD4
  table to it through a7's sample key, and for the WGS cohort (HIV + HL pooled)
  and again for HIV only computes:

    1. SPEARMAN rank correlation between CD4 and each of
       anello_reads_human_total, anello_rpm_human, anello_richness_human and
       anello_shannon_human. Ranks are tie-corrected (average ranks); rho is the
       Pearson correlation of those ranks; the two-sided p comes from
       t = rho * sqrt((n - 2) / (1 - rho^2)) on n - 2 degrees of freedom, with
       the Student-t tail evaluated from the regularised incomplete beta
       function implemented here.
    2. KENDALL tau-b over the same pairs as a second, more tie-robust estimate,
       with the tie-corrected normal-approximation two-sided p.
    3. STRATA. The standard clinical CD4 bands (<200, 200-499, >=500 cells/uL,
       configurable with --strata-breaks): n, median and IQR of every burden
       metric per band, plus a tie-corrected Kruskal-Wallis H across bands with
       a chi-square p on k - 1 degrees of freedom.
    4. A TWO-GROUP CONTRAST on CD4 itself: CD4 in samples with
       anello_richness_human >= --richness-threshold (default 3, the cut used in
       the a7 headline) versus CD4 in samples below it, tested with the same
       standard-library Mann-Whitney U as a7 (tie- and continuity-corrected
       normal approximation, rank-biserial effect size).

  Every statistic is implemented in the standard library; scipy is not imported
  and no figure is drawn, so matplotlib is not imported either. Every reported
  statistic carries its n in the same row.

WHAT THIS CAN AND CANNOT SUPPORT
  CAN: it says whether, among these samples, a lower CD4 count goes with a
  higher measured anellovirus burden, and how big that association is.
  CANNOT, and the same three sentences are printed to stdout:
    - The outcome is mostly zero. Most samples carry no detected anellovirus at
      all, so a rank correlation is dominated by a large tied block of zeros; a
      handful of high-burden samples can move rho on their own. Treat any
      correlation over a mostly-zero outcome as weak evidence and read the
      per-stratum medians and the zero counts (n_zero_metric) beside it.
    - CD4 and HIV status are collinear in this cohort. HL samples are not
      CD4-depleted, so a pooled HIV+HL correlation largely re-discovers the
      group difference a7 already reported. THE HIV-ONLY ANALYSIS IS THE
      INFORMATIVE ONE; the pooled row is kept only for completeness.
    - Direction of causation is not identifiable here. Cross-sectional CD4 and
      cross-sectional viral read counts cannot separate "immunosuppression
      permits anellovirus expansion" from "anellovirus burden tracks something
      else that also tracks CD4" (ART exposure, time since diagnosis, sample
      handling, sequencing depth). Nothing in this module adjusts for depth
      beyond a7's RPM normalisation.
  The richness>=3 contrast is a re-expression of the same anellovirus numbers,
  not independent evidence, and in the pooled cohort it mostly recovers HIV
  status because all a7 richness>=3 samples were HIV+.

INPUT
  --cd4 <tsv/csv>  a user-supplied table with a sample identifier column
      (sample / sample_id / subject / id / ...) and a CD4 column
      (cd4 / cd4_count / cd4_abs / cd4_cells / ...). The delimiter is sniffed
      (tab, comma, semicolon or pipe) and the columns are matched case- and
      punctuation-insensitively. Extra clinical columns (viral_load, treatment,
      ...) are carried through to the joined table verbatim and are NEVER
      interpreted. Columns whose name looks like a direct identifier or a date
      (mrn, patient, name, dob, birth, date, address, phone, email, ssn,
      accession, barcode) are dropped rather than carried, and any carried value
      that is itself a known real sample name is replaced by its anonymous id.
  If --cd4 is absent or missing, the module writes cd4_input_template.tsv (the
  anonymised WGS samples with a blank cd4 column), prints how to fill it in, and
  exits 0 without computing anything.

WHAT IT WRITES (tab separated, pure ASCII, into --outdir)
  cd4_anello_correlation.tsv          Spearman and Kendall per cohort x metric.
  cd4_anello_strata.tsv               per-CD4-band n / median / IQR + Kruskal-Wallis.
  cd4_anello_richness_contrast.tsv    Mann-Whitney on CD4 by richness >= threshold.
  cd4_anello_joined.tsv               per-sample CD4 + burden + carried columns.
  cd4_input_template.tsv              only when no usable CD4 table was given.
  a9_cd4_sample_key.tsv               only when the CD4 file used real names.
                                      THE ONLY file that may contain a real
                                      sample identifier; its first line is
                                      "# CONTAINS IDENTIFIERS - DO NOT COMMIT OR
                                      EMAIL". Everything else is S01..Snn.
  With a non-empty --prefix every name above is prefixed (a9_cd4_correlation...).

EXAMPLE
  # 1. no CD4 data yet: get the template
  python3 a9_cd4_correlation.py --indir /path/to/suite_out --outdir /path/to/suite_out

  # 2. with a filled-in clinical table
  python3 a9_cd4_correlation.py \
      --indir /path/to/suite_out \
      --cd4 /path/to/clinical/hiv_cd4_counts.tsv \
      --outdir /path/to/suite_out \
      --richness-threshold 3 --strata-breaks 200,500

  SUITE_OUTDIR and CD4_TABLE are honoured as defaults for --indir and --cd4.
  A missing input is reported as "WARN: <what> missing at <path>, skipping" and
  the module still exits 0. No network access. Dates are YYYY-MM-DD.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import math
import os
import re
import sys

DEFAULT_INDIR = os.environ.get("SUITE_OUTDIR", "/path/to/suite_out")
DEFAULT_OUTDIR = os.environ.get("A9_OUTDIR", "") or DEFAULT_INDIR
DEFAULT_CD4 = os.environ.get("CD4_TABLE", "")
DEFAULT_A7_PREFIX = "a7_virome"

# The key file keeps this stem even when --prefix is empty, so that the one file
# that may carry identifiers is always obvious in a directory listing.
KEY_STEM = "a9_cd4"

# Burden metrics taken from the a7 table, with the scale label reported beside
# each statistic. "is_count" marks metrics that are integers in the a7 table.
METRIC_SPECS = [
    ("anello_reads_human_total", "raw_reads", True),
    ("anello_rpm_human", "normalised_rpm", False),
    ("anello_richness_human", "count_of_references", True),
    ("anello_shannon_human", "shannon_index_natural_log", False),
]
METRICS = [m[0] for m in METRIC_SPECS]
SCALE_OF = dict((m, s) for m, s, _i in METRIC_SPECS)

RICHNESS_METRIC = "anello_richness_human"

# Identifier-column detection, in priority order.
SAMPLE_COL_CANDIDATES = [
    "sample", "sample_id", "sampleid", "sample_name", "samplename",
    "subject", "subject_id", "subjectid", "id", "anon_sample", "sample_anon",
    "specimen", "specimen_id", "library", "library_id",
]
CD4_COL_CANDIDATES = [
    "cd4", "cd4_count", "cd4count", "cd4_abs", "cd4_absolute", "cd4_cells",
    "cd4_cells_ul", "cd4_cells_per_ul", "abs_cd4", "cd4_absolute_count",
    "cd4_t_cell_count", "cd4_nadir",
]

# Extra columns whose NAME suggests a direct identifier or a date. They are
# dropped instead of carried, because everything except the key file must be
# safe to share. HIPAA treats dates finer than a year as identifiers, hence the
# blunt "date" rule.
BLOCKED_EXTRA_TOKENS = [
    "mrn", "medical_record", "patient", "name", "initial", "dob", "birth",
    "date", "address", "zip", "postal", "phone", "email", "ssn", "nhs",
    "accession", "barcode", "record_number", "chart", "encounter",
]

DELIMITERS = ["\t", ",", ";", "|"]
DELIM_NAMES = {"\t": "tab", ",": "comma", ";": "semicolon", "|": "pipe"}

TODAY = datetime.date.today().isoformat()
SCRIPT = os.path.basename(__file__)


# --------------------------------------------------------------------------- #
# small utilities (same shapes as a5 / a7)
# --------------------------------------------------------------------------- #
def warn(what, path):
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


def out_name(prefix, base):
    return (prefix + "_" + base) if prefix else base


def write_tsv(path, comments, header, rows):
    with open(path, "w", encoding="ascii", errors="replace", newline="") as fh:
        for line in comments:
            fh.write("# " + to_ascii(line) + "\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(to_ascii(c) for c in row) + "\n")


def mean(values):
    if not values:
        return None
    return float(sum(values)) / float(len(values))


def median(values):
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return float(s[mid])
    return (float(s[mid - 1]) + float(s[mid])) / 2.0


def quantile(values, q):
    """Linear interpolation between order statistics (the numpy default)."""
    if not values:
        return None
    s = sorted(float(v) for v in values)
    if len(s) == 1:
        return s[0]
    pos = (len(s) - 1) * float(q)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return s[lo]
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


# --------------------------------------------------------------------------- #
# standard-library distribution tails
#
# Student t and chi-square tails are needed for the Spearman and Kruskal-Wallis
# p values. Both are evaluated from series / continued-fraction expansions of
# the incomplete beta and incomplete gamma functions (Numerical Recipes forms),
# so nothing outside math is required.
# --------------------------------------------------------------------------- #
def _betacf(a, b, x, itmax=300, eps=3.0e-12):
    tiny = 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betai(a, b, x):
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    try:
        lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                 + a * math.log(x) + b * math.log1p(-x))
        front = math.exp(lbeta)
    except (ValueError, OverflowError):
        return None
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t_stat, df):
    """Two-sided Student-t tail probability."""
    if df is None or df <= 0:
        return None
    if t_stat != t_stat:
        return None
    t2 = float(t_stat) * float(t_stat)
    if t2 == float("inf"):
        return 0.0
    p = betai(df / 2.0, 0.5, df / (df + t2))
    if p is None:
        return None
    return min(1.0, max(0.0, p))


def _gser(a, x, itmax=1000, eps=3.0e-14):
    ap = a
    total = 1.0 / a
    delta = total
    for _ in range(itmax):
        ap += 1.0
        delta *= x / ap
        total += delta
        if abs(delta) < abs(total) * eps:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a, x, itmax=1000, eps=3.0e-14):
    tiny = 1.0e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b if b != 0 else 1.0 / tiny
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def chi2_sf(x, df):
    """Upper tail of the chi-square distribution, Q(df/2, x/2)."""
    if df is None or df <= 0 or x is None:
        return None
    if x <= 0:
        return 1.0
    a = float(df) / 2.0
    z = float(x) / 2.0
    try:
        if z < a + 1.0:
            p = 1.0 - _gser(a, z)
        else:
            p = _gcf(a, z)
    except (ValueError, OverflowError):
        return None
    return min(1.0, max(0.0, p))


def normal_sf_two_sided(z):
    return min(1.0, max(0.0, math.erfc(abs(z) / math.sqrt(2.0))))


# --------------------------------------------------------------------------- #
# rank statistics (standard library only)
# --------------------------------------------------------------------------- #
def average_ranks(values):
    """Average (tie-corrected) ranks, plus the sizes of the tied groups."""
    n = len(values)
    order = sorted(range(n), key=lambda i: float(values[i]))
    ranks = [0.0] * n
    tie_sizes = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and float(values[order[j + 1]]) == float(values[order[i]]):
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        tie_sizes.append(j - i + 1)
        i = j + 1
    return ranks, tie_sizes


def pearson(x, y):
    n = len(x)
    if n < 2:
        return None
    mx, my = mean(x), mean(y)
    sxy = sxx = syy = 0.0
    for a, b in zip(x, y):
        da, db = float(a) - mx, float(b) - my
        sxy += da * db
        sxx += da * da
        syy += db * db
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    """Tie-corrected Spearman rho with a two-sided t-approximation p value.

    Returns a dict, or None when n < 3. rho is None (with a reason) when either
    variable has no variance, which happens whenever a burden metric is zero in
    every joined sample.
    """
    n = len(x)
    if n != len(y) or n < 3:
        return None
    rx, tx = average_ranks(x)
    ry, ty = average_ranks(y)
    rho = pearson(rx, ry)
    if rho is None:
        return {"n": n, "rho": None, "t": None, "df": n - 2, "p": None,
                "n_tied_groups_x": sum(1 for t in tx if t > 1),
                "n_tied_groups_y": sum(1 for t in ty if t > 1),
                "largest_tie_y": max(ty) if ty else 0,
                "note": "no variance in one variable; rho undefined"}
    rho = max(-1.0, min(1.0, rho))
    if abs(rho) >= 1.0:
        t_stat, p = None, 0.0
    else:
        t_stat = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
        p = t_sf_two_sided(t_stat, n - 2)
    return {"n": n, "rho": rho, "t": t_stat, "df": n - 2, "p": p,
            "n_tied_groups_x": sum(1 for t in tx if t > 1),
            "n_tied_groups_y": sum(1 for t in ty if t > 1),
            "largest_tie_y": max(ty) if ty else 0,
            "note": ""}


def kendall_tau_b(x, y):
    """Kendall tau-b with a tie-corrected normal-approximation two-sided p.

    Concordant / discordant pairs are counted directly (O(n^2)); the cohorts
    here are tens of samples, so the naive count is fast enough and is far
    easier to audit than a merge-sort count.
    """
    n = len(x)
    if n != len(y) or n < 3:
        return None
    con = dis = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            prod = (float(x[i]) - float(x[j])) * (float(y[i]) - float(y[j]))
            if prod > 0:
                con += 1
            elif prod < 0:
                dis += 1
    _rx, tx = average_ranks(x)
    _ry, ty = average_ranks(y)
    xtie = sum(t * (t - 1) / 2.0 for t in tx if t > 1)
    ytie = sum(t * (t - 1) / 2.0 for t in ty if t > 1)
    n0 = n * (n - 1) / 2.0
    denom = math.sqrt((n0 - xtie) * (n0 - ytie))
    tau = ((con - dis) / denom) if denom > 0 else None
    x0 = sum(t * (t - 1) for t in tx if t > 1)
    y0 = sum(t * (t - 1) for t in ty if t > 1)
    x1 = sum(t * (t - 1) * (2 * t + 5) for t in tx if t > 1)
    y1 = sum(t * (t - 1) * (2 * t + 5) for t in ty if t > 1)
    var = ((n * (n - 1) * (2.0 * n + 5) - x1 - y1) / 18.0
           + (2.0 * xtie * ytie) / (n * (n - 1))
           + (x0 * y0) / (9.0 * n * (n - 1) * (n - 2)))
    if var > 0:
        z = (con - dis) / math.sqrt(var)
        p = normal_sf_two_sided(z)
    else:
        z, p = None, None
    return {"n": n, "tau_b": tau, "concordant": con, "discordant": dis,
            "z": z, "p": p,
            "note": "" if tau is not None else "all pairs tied; tau undefined"}


def kruskal_wallis(groups):
    """Tie-corrected Kruskal-Wallis H over a list of value lists."""
    used = [g for g in groups if g]
    if len(used) < 2:
        return None
    pooled = []
    for gi, g in enumerate(used):
        for v in g:
            pooled.append((float(v), gi))
    n = len(pooled)
    if n < 3:
        return None
    ranks, tie_sizes = average_ranks([p[0] for p in pooled])
    rank_sum = [0.0] * len(used)
    counts = [0] * len(used)
    for (_val, gi), r in zip(pooled, ranks):
        rank_sum[gi] += r
        counts[gi] += 1
    h = 0.0
    for gi in range(len(used)):
        if counts[gi]:
            h += (rank_sum[gi] ** 2) / float(counts[gi])
    h = 12.0 / (n * (n + 1.0)) * h - 3.0 * (n + 1.0)
    tie_term = sum(t ** 3 - t for t in tie_sizes)
    correction = 1.0 - tie_term / float(n ** 3 - n) if n > 1 else 1.0
    if correction > 0:
        h = h / correction
    df = len(used) - 1
    return {"H": h, "df": df, "p": chi2_sf(h, df), "n_total": n,
            "n_groups_used": len(used), "counts": counts,
            "tie_correction": correction}


def mann_whitney_u(x, y):
    """Two-sided Mann-Whitney U with tie correction and normal approximation.

    Identical implementation to a7_virome_structure.py, so that the U / z / p
    columns of the two modules are directly comparable.
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
    var = (float(n1) * float(n2) / 12.0) * ((n + 1) - tie_term / float(n * (n - 1)))
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
        p = normal_sf_two_sided(z)
    return {
        "n1": n1, "n2": n2, "U1": u1, "U2": u2, "z": z, "p": p,
        "median1": median(x), "median2": median(y),
        "mean1": mean(x), "mean2": mean(y),
        "effect_r": 2.0 * u1 / (float(n1) * float(n2)) - 1.0,
        "n_ties_groups": sum(1 for t in tie_sizes if t > 1),
    }


# --------------------------------------------------------------------------- #
# sample naming (group_of copied from a7_virome_structure.py)
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


# --------------------------------------------------------------------------- #
# generic table reading: delimiter sniffing and flexible column detection
# --------------------------------------------------------------------------- #
def norm_header(name):
    """Header name -> lower-case, punctuation-collapsed comparison key."""
    text = to_ascii(name).strip().strip('"').strip("'").lower()
    return re.sub(r"[^0-9a-z]+", "_", text).strip("_")


def content_lines(path):
    """All non-blank, non-comment lines of a text file, or None if unreadable."""
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            raw = fh.read().splitlines()
    except OSError:
        return None
    return [ln for ln in raw if ln.strip() and not ln.lstrip().startswith("#")]


def sniff_delimiter(lines):
    """-> (delimiter, how_it_was_chosen) or (None, reason)."""
    if not lines:
        return None, "file has no data lines"
    header = lines[0]
    counts = dict((d, header.count(d)) for d in DELIMITERS)
    best = max(counts.values())
    if best > 0:
        # ties resolve in DELIMITERS order: tab, comma, semicolon, pipe
        for d in DELIMITERS:
            if counts[d] == best:
                return d, "counted %d %s separators in the header line" % (
                    best, DELIM_NAMES[d])
    try:
        dialect = csv.Sniffer().sniff("\n".join(lines[:20]),
                                      delimiters="".join(DELIMITERS))
        return dialect.delimiter, "csv.Sniffer"
    except (csv.Error, TypeError):
        return None, "no tab, comma, semicolon or pipe found in the header line"


def parse_table(lines, delim):
    """-> (raw header, normalised header, rows as lists of cells)."""
    rdr = csv.reader(lines, delimiter=delim)
    try:
        header = [to_ascii(h) for h in next(rdr)]
    except StopIteration:
        return [], [], []
    rows = []
    for parts in rdr:
        if not parts or all(not to_ascii(p) for p in parts):
            continue
        rows.append([to_ascii(p) for p in parts])
    return header, [norm_header(h) for h in header], rows


def pick_column(norm_head, candidates, tokens):
    """Index of the first header matching a candidate, then a prefix, then a
    substring. Returns (index, how) or (None, "")."""
    for cand in candidates:
        for i, h in enumerate(norm_head):
            if h == cand:
                return i, "exact header '%s'" % h
    for token in tokens:
        for i, h in enumerate(norm_head):
            if h.startswith(token):
                return i, "header '%s' starts with '%s'" % (h, token)
    for token in tokens:
        for i, h in enumerate(norm_head):
            if token in h:
                return i, "header '%s' contains '%s'" % (h, token)
    return None, ""


def as_float(text):
    """-> (value, flag). flag is '' , 'censored' or 'unparsed'.

    Tolerates thousands separators, a trailing unit ("cells/uL"), and a leading
    censoring operator ("<20"), which is kept as the bare number and counted.
    """
    raw = to_ascii(text).strip().strip('"').strip("'")
    if raw == "" or raw.upper() in ("NA", "N/A", "NAN", "NULL", "ND", "."):
        return None, "unparsed"
    cleaned = raw.replace(",", "").replace(" ", "")
    m = re.match(r"^([<>]=?)?(-?[0-9]*\.?[0-9]+)", cleaned)
    if not m:
        return None, "unparsed"
    try:
        value = float(m.group(2))
    except ValueError:
        return None, "unparsed"
    return value, ("censored" if m.group(1) else "")


def norm_id(text, loose=False):
    """Join key: trimmed, internal whitespace collapsed, lower-cased.

    loose=True additionally removes every non-alphanumeric character, which
    absorbs the usual "-" vs "_" vs "." drift between a clinical spreadsheet and
    a sequencing sample sheet.
    """
    s = to_ascii(text).strip().strip('"').strip("'")
    s = re.sub(r"\s+", " ", s).strip().lower()
    if loose:
        s = re.sub(r"[^0-9a-z]+", "", s)
    return s


# --------------------------------------------------------------------------- #
# a7 inputs
# --------------------------------------------------------------------------- #
def read_suite_tsv(path):
    """A tab-separated suite output: '#' comments stripped. -> (norm_head, rows)."""
    lines = content_lines(path)
    if lines is None:
        return None, None
    header, norm_head, rows = parse_table(lines, "\t")
    if not header:
        return [], []
    dict_rows = []
    for parts in rows:
        row = {}
        for i, key in enumerate(norm_head):
            row[key] = parts[i] if i < len(parts) else ""
        dict_rows.append(row)
    return norm_head, dict_rows


def load_sample_key(path):
    """-> (anon -> real, real -> anon, note). ({}, {}, reason) when unusable.

    Accepts both key layouts in the suite: a7's anon_sample / real_sample and
    a5's sample_anon / sample_real.
    """
    if not os.path.exists(path):
        warn("a7 sample key", path)
        return {}, {}, "sample key not found"
    norm_head, rows = read_suite_tsv(path)
    if not norm_head:
        warn("readable a7 sample key", path)
        return {}, {}, "sample key unreadable"
    anon_col = next((c for c in ("anon_sample", "sample_anon", "anon", "sample_id")
                     if c in norm_head), None)
    real_col = next((c for c in ("real_sample", "sample_real", "sample", "real_name")
                     if c in norm_head and c != anon_col), None)
    if anon_col is None or real_col is None:
        warn("anon / real columns in the a7 sample key", path)
        return {}, {}, "sample key has no anon/real column pair"
    anon_to_real, real_to_anon = {}, {}
    for row in rows:
        anon = (row.get(anon_col) or "").strip()
        real = (row.get(real_col) or "").strip()
        if not anon or not real:
            continue
        anon_to_real[anon] = real
        real_to_anon[real] = anon
    return anon_to_real, real_to_anon, ""


def load_burden(path):
    """-> list of per-sample dicts from a7_virome_anellovirus_burden.tsv."""
    if not os.path.exists(path):
        warn("a7 anellovirus burden table", path)
        return None
    norm_head, rows = read_suite_tsv(path)
    if not norm_head:
        warn("readable a7 anellovirus burden table", path)
        return None
    id_col = next((c for c in ("sample", "sample_anon", "anon_sample")
                   if c in norm_head), None)
    if id_col is None:
        warn("sample column in the a7 anellovirus burden table", path)
        return None
    missing = [m for m in METRICS if m not in norm_head]
    if len(missing) == len(METRICS):
        warn("anellovirus metric columns in the a7 burden table", path)
        return None
    if missing:
        print("NOTE: the a7 burden table has no %s column(s); those metrics are "
              "reported as NA" % ", ".join(missing))
    out = []
    for row in rows:
        anon = (row.get(id_col) or "").strip()
        if not anon:
            continue
        rec = {"anon": anon,
               "group": (row.get("group") or "NA").strip() or "NA",
               "run": (row.get("run") or "").strip()}
        for metric in METRICS:
            value, _flag = as_float(row.get(metric, ""))
            rec[metric] = value
        out.append(rec)
    return out


def build_lookup(real_to_anon):
    """Three join maps plus the set of keys that collide under normalisation."""
    exact, normed, loose = {}, {}, {}
    clashes = set()
    for real, anon in real_to_anon.items():
        exact[real] = anon
        for table, key in ((normed, norm_id(real)), (loose, norm_id(real, True))):
            if not key:
                continue
            if key in table and table[key] != anon:
                clashes.add(key)
            else:
                table[key] = anon
    for key in clashes:
        normed.pop(key, None)
        loose.pop(key, None)
    return exact, normed, loose, clashes


def resolve_id(raw, exact, normed, loose, anon_set):
    """-> (anon_id or None, how, is_real_identifier)."""
    text = to_ascii(raw).strip()
    if not text:
        return None, "blank", False
    if text in anon_set:
        return text, "anon_id_used_directly", False
    upper = text.upper()
    if upper in anon_set:
        return upper, "anon_id_case_normalised", False
    if text in exact:
        return exact[text], "exact", True
    key = norm_id(text)
    if key in normed:
        return normed[key], "whitespace_and_case_normalised", True
    key = norm_id(text, True)
    if key in loose:
        return loose[key], "alphanumeric_only_normalised", True
    return None, "unmatched", True


# --------------------------------------------------------------------------- #
# the CD4 table
# --------------------------------------------------------------------------- #
def load_cd4_table(path):
    """-> dict describing the user's CD4 table, or None when it is unusable."""
    if not path:
        return None
    if not os.path.exists(path):
        warn("CD4 table", path)
        return None
    lines = content_lines(path)
    if lines is None:
        warn("readable CD4 table", path)
        return None
    if len(lines) < 2:
        warn("data rows in the CD4 table", path)
        return None
    delim, how = sniff_delimiter(lines)
    if delim is None:
        warn("a usable delimiter in the CD4 table (%s)" % how, path)
        return None
    header, norm_head, rows = parse_table(lines, delim)
    if not header or not rows:
        warn("a header and at least one data row in the CD4 table", path)
        return None
    sample_idx, sample_how = pick_column(
        norm_head, SAMPLE_COL_CANDIDATES,
        ["sample", "subject", "specimen", "library", "patient_id"])
    if sample_idx is None:
        warn("a sample identifier column in the CD4 table (looked for %s)"
             % "/".join(SAMPLE_COL_CANDIDATES[:6]), path)
        return None
    cd4_idx, cd4_how = pick_column(norm_head, CD4_COL_CANDIDATES, ["cd4"])
    if cd4_idx is None:
        warn("a CD4 column in the CD4 table (looked for %s)"
             % "/".join(CD4_COL_CANDIDATES[:5]), path)
        return None

    extras, dropped = [], []
    used = set()
    for i, name in enumerate(norm_head):
        if i in (sample_idx, cd4_idx):
            continue
        if not name:
            continue
        if any(tok in name for tok in BLOCKED_EXTRA_TOKENS):
            dropped.append(name)
            continue
        col = "extra_" + name
        if col in used:                       # duplicate header in the source
            col = "%s_%d" % (col, i)
        used.add(col)
        extras.append((i, col, header[i]))
    return {
        "path": path, "delimiter": delim, "delimiter_how": how,
        "header": header, "norm_head": norm_head, "rows": rows,
        "sample_idx": sample_idx, "sample_how": sample_how,
        "sample_col": header[sample_idx],
        "cd4_idx": cd4_idx, "cd4_how": cd4_how, "cd4_col": header[cd4_idx],
        "extras": extras, "dropped": dropped,
        # a row with more fields than the header usually means an unquoted
        # separator inside a value ("1,094" in a comma-separated file), which
        # silently shifts every column to its right
        "ragged": sum(1 for r in rows if len(r) != len(header)),
    }


def build_redactor(real_to_anon):
    """-> redact(value) -> (safe_value, changed).

    Carried clinical columns are free text, so a real sample name can sit
    INSIDE a longer string ("prior=<name>", "aliquot <name> B"). Whole-cell
    equality is not enough, so every occurrence is substituted by the anonymous
    id, and a punctuation-drifted occurrence that survives that pass costs the
    whole cell. Names shorter than 4 characters are only matched whole-cell,
    because a two-character "id" would otherwise rewrite ordinary text.
    """
    long_names = sorted((n for n in real_to_anon if n and len(n) >= 4),
                        key=len, reverse=True)
    lower_map = dict((n.lower(), real_to_anon[n]) for n in long_names)
    loose_names = [(norm_id(n, True), real_to_anon[n]) for n in long_names]
    loose_names = [(k, v) for k, v in loose_names if len(k) >= 4]
    short_exact = dict((norm_id(n, True), a) for n, a in real_to_anon.items()
                       if n and len(n) < 4)
    pattern = (re.compile("|".join(re.escape(n) for n in long_names), re.IGNORECASE)
               if long_names else None)

    def redact(value):
        if not value:
            return value, False
        short_hit = short_exact.get(norm_id(value, True))
        if short_hit:
            return short_hit, True
        if pattern is None:
            return value, False
        new = pattern.sub(
            lambda m: lower_map.get(m.group(0).lower(), "REDACTED_IDENTIFIER"),
            value)
        loose_new = norm_id(new, True)
        for loose_name, _anon in loose_names:
            if loose_name in loose_new:
                return "REDACTED_CONTAINED_IDENTIFIER", True
        return new, new != value

    return redact


def join_cd4(table, burden_by_anon, exact, normed, loose, anon_set, real_to_anon):
    """Join CD4 rows onto the burden records. -> (records, report dict)."""
    records = []
    report = {
        "rows_read": len(table["rows"]), "matched": 0,
        "unmatched_known_anon": [], "unmatched_unknown_count": 0,
        "blank_id": 0, "bad_cd4": 0, "censored_cd4": 0, "implausible_cd4": 0,
        "negative_cd4": 0, "duplicate_anon": [], "redacted_extra_values": 0,
        "match_methods": {}, "real_identifier_rows": 0,
    }
    # a real sample name must never leave inside a carried-through column
    redact = build_redactor(real_to_anon)
    seen = set()
    for parts in table["rows"]:
        def cell(i):
            return parts[i] if i < len(parts) else ""

        raw_id = cell(table["sample_idx"])
        if not to_ascii(raw_id):
            report["blank_id"] += 1
            continue
        anon, how, is_real = resolve_id(raw_id, exact, normed, loose, anon_set)
        if anon is None:
            report["unmatched_unknown_count"] += 1
            continue
        if anon not in burden_by_anon:
            if anon not in report["unmatched_known_anon"]:
                report["unmatched_known_anon"].append(anon)
            continue
        if anon in seen:
            if anon not in report["duplicate_anon"]:
                report["duplicate_anon"].append(anon)
            continue
        cd4, flag = as_float(cell(table["cd4_idx"]))
        if cd4 is None:
            report["bad_cd4"] += 1
            continue
        if cd4 < 0:
            report["negative_cd4"] += 1
            continue
        if flag == "censored":
            report["censored_cd4"] += 1
        if cd4 > 5000:
            report["implausible_cd4"] += 1
        seen.add(anon)
        report["match_methods"][how] = report["match_methods"].get(how, 0) + 1
        if is_real:
            report["real_identifier_rows"] += 1
        rec = dict(burden_by_anon[anon])
        rec["cd4"] = cd4
        rec["cd4_censored"] = 1 if flag == "censored" else 0
        rec["id_match_method"] = how
        # True only when this row named the sample by its REAL name. The key
        # file must not re-export a real name the CD4 table never carried.
        rec["id_was_real"] = bool(is_real)
        rec["extras"] = {}
        for idx, col, _raw_name in table["extras"]:
            value, changed = redact(to_ascii(cell(idx)))
            if changed:
                report["redacted_extra_values"] += 1
            rec["extras"][col] = value
        records.append(rec)
        report["matched"] += 1
    records.sort(key=lambda r: r["anon"])
    return records, report


# --------------------------------------------------------------------------- #
# strata
# --------------------------------------------------------------------------- #
def fmt_break(value):
    return str(int(value)) if float(value).is_integer() else ("%g" % value)


def strata_labels(breaks):
    """-> list of (label, range_text) for the half-open CD4 bands."""
    out = [("CD4_lt_%s" % fmt_break(breaks[0]), "[0,%s)" % fmt_break(breaks[0]))]
    for i in range(len(breaks) - 1):
        lo, hi = breaks[i], breaks[i + 1]
        if float(lo).is_integer() and float(hi).is_integer():
            label = "CD4_%d_to_%d" % (int(lo), int(hi) - 1)
        else:
            label = "CD4_%s_to_%s" % (fmt_break(lo), fmt_break(hi))
        out.append((label, "[%s,%s)" % (fmt_break(lo), fmt_break(hi))))
    out.append(("CD4_ge_%s" % fmt_break(breaks[-1]),
                "[%s,inf)" % fmt_break(breaks[-1])))
    return out


def stratum_index(cd4, breaks):
    return sum(1 for b in breaks if cd4 >= b)


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def common_comments(args, key_basename, burden_basename, extra=()):
    lines = [
        "%s  generated %s" % (SCRIPT, TODAY),
        "no real sample identifiers in this file; ids are anonymous S01..Snn "
        "(a7 mapping in %s)" % burden_basename.replace(
            "anellovirus_burden.tsv", "sample_key.tsv"),
        "burden metrics are read verbatim from %s; CD4 is user-supplied and is "
        "assumed to be an absolute count in cells/uL" % burden_basename,
        "CAVEAT 1: the outcome is mostly zero, so a rank correlation over it is "
        "weak evidence; read n_zero_metric and the per-stratum medians beside it",
        "CAVEAT 2: CD4 and HIV status are collinear here, so the HIV_ONLY rows "
        "are the informative ones and the pooled WGS rows largely restate a7",
        "CAVEAT 3: cross-sectional association only; this cannot establish the "
        "direction of causation",
    ]
    if key_basename:
        lines.append("identifiers handled this run: see %s (do not commit or email)"
                     % key_basename)
    lines.extend(extra)
    return lines


def write_sample_key(path, records, anon_to_real):
    lines = [
        "CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
        "generated %s by %s" % (TODAY, SCRIPT),
        "only the samples whose real name appeared in the supplied CD4 table; "
        "rows joined by an anonymous id are deliberately not listed here",
    ]
    rows = []
    for rec in sorted(records, key=lambda r: r["anon"]):
        if not rec.get("id_was_real"):
            continue
        real = anon_to_real.get(rec["anon"], "")
        if not real:
            continue
        rows.append([rec["anon"], real, rec["group_final"], rec["id_match_method"]])
    write_tsv(path, lines,
              ["anon_sample", "real_sample", "group", "id_match_method"], rows)
    return len(rows)


def write_template(path, samples, args, burden_basename):
    header = ["sample_id", "group", "anello_richness_human", "anello_rpm_human",
              "cd4", "viral_load", "treatment", "notes"]
    rows = []
    for rec in samples:
        rows.append([rec["anon"], rec["group_final"],
                     fnum(rec.get(RICHNESS_METRIC), 0),
                     fnum(rec.get("anello_rpm_human"), 3),
                     "", "", "", ""])
    comments = [
        "%s  generated %s" % (SCRIPT, TODAY),
        "TEMPLATE - no CD4 table was supplied, so nothing was computed.",
        "sample_id is the anonymous a7 id. Resolve it to a real sample name with "
        "the a7 key (%s), which carries identifiers: keep it out of anything you "
        "commit or email." % burden_basename.replace(
            "anellovirus_burden.tsv", "sample_key.tsv"),
        "Fill in the cd4 column with the absolute CD4 count in cells/uL. A plain "
        "number is safest; '<20' is read as 20; a thousands separator ('1,234') "
        "is only safe here or in a quoted CSV field. Blank means unknown.",
        "viral_load / treatment / notes are optional free columns; they are "
        "carried through uninterpreted. Columns whose name looks like an "
        "identifier or a date are dropped on re-read.",
        "Then rerun: %s --cd4 <this file> --indir <suite outdir>" % SCRIPT,
        "The burden columns are shown only so the row is recognisable; do not "
        "let them influence what you enter.",
    ]
    write_tsv(path, comments, header, rows)


def write_joined(path, records, args, key_basename, burden_basename,
                 extra_cols, breaks, labels):
    header = (["sample", "group", "cd4_cells_ul", "cd4_censored", "cd4_stratum",
               "richness_class", "id_match_method"] + METRICS + extra_cols)
    rows = []
    for rec in records:
        idx = stratum_index(rec["cd4"], breaks)
        klass = ("richness_ge_%d" % args.richness_threshold
                 if (rec.get(RICHNESS_METRIC) or 0) >= args.richness_threshold
                 else "richness_lt_%d" % args.richness_threshold)
        row = [rec["anon"], rec["group_final"], fnum(rec["cd4"], 1),
               str(rec["cd4_censored"]), labels[idx][0], klass,
               rec["id_match_method"]]
        for metric, _scale, is_count in METRIC_SPECS:
            value = rec.get(metric)
            row.append(fnum(value, 0 if is_count else 4))
        for col in extra_cols:
            row.append(rec["extras"].get(col, ""))
        rows.append(row)
    write_tsv(path, common_comments(
        args, key_basename, burden_basename,
        ["one row per joined sample; extra_* columns are carried through from "
         "the CD4 table verbatim and are never interpreted",
         "id join normalisation, in order: exact, then trimmed + internal "
         "whitespace collapsed + lower-cased, then alphanumeric characters only; "
         "id_match_method records which one matched each row",
         "a carried value that matched a real sample name was replaced by its "
         "anonymous id",
         "samples whose group is outside the analysed cohorts are listed here "
         "for completeness and take part in no test"]), header, rows)


CORR_HEADER = [
    "cohort", "n_cohort", "metric", "metric_scale", "n_pairs",
    "n_zero_metric", "n_nonzero_metric", "frac_zero_metric",
    "cd4_median", "cd4_min", "cd4_max",
    "method", "statistic_name", "statistic",
    "test_statistic_name", "test_statistic", "df", "p_two_sided",
    "n_concordant_pairs", "n_discordant_pairs", "reading", "note",
]


def reading_of(stat):
    if stat is None:
        return "NA"
    if stat < 0:
        return "burden_higher_at_lower_cd4"
    if stat > 0:
        return "burden_higher_at_higher_cd4"
    return "no_monotone_trend"


def correlation_rows(cohorts, args):
    """-> (rows, {(cohort, metric): spearman_dict})."""
    rows = []
    keep = {}
    for cohort_name, recs in cohorts:
        for metric in METRICS:
            pairs = [(r["cd4"], r[metric]) for r in recs if r.get(metric) is not None]
            x = [p[0] for p in pairs]
            y = [p[1] for p in pairs]
            n_pairs = len(pairs)
            n_zero = sum(1 for v in y if float(v) == 0.0)
            frac_zero = (float(n_zero) / n_pairs) if n_pairs else None
            base_note = []
            if n_pairs < args.min_n:
                base_note.append("n=%d < %d; descriptive only" % (n_pairs, args.min_n))
            if frac_zero is not None and frac_zero >= 0.5:
                base_note.append("%d/%d samples are zero for this metric; the rank "
                                 "correlation is dominated by the tied zero block"
                                 % (n_zero, n_pairs))
            common = [cohort_name, str(len(recs)), metric, SCALE_OF[metric],
                      str(n_pairs), str(n_zero), str(n_pairs - n_zero),
                      fnum(frac_zero, 3), fnum(median(x), 1),
                      fnum(min(x) if x else None, 1), fnum(max(x) if x else None, 1)]

            sp = spearman(x, y)
            keep[(cohort_name, metric)] = sp
            if sp is None:
                rows.append(common + ["spearman_rank", "rho", "NA", "t", "NA", "NA",
                                      "NA", "NA", "NA", "NA",
                                      "; ".join(base_note + ["n < 3; test not run"])])
            else:
                note = base_note + ([sp["note"]] if sp["note"] else [])
                if sp["largest_tie_y"] > 1:
                    note.append("largest tied block in the metric: %d values"
                                % sp["largest_tie_y"])
                rows.append(common + [
                    "spearman_rank", "rho", fnum(sp["rho"], 4),
                    "t", fnum(sp["t"], 4), str(sp["df"]), fp(sp["p"]),
                    "NA", "NA", reading_of(sp["rho"]), "; ".join(note)])

            kd = kendall_tau_b(x, y)
            if kd is None:
                rows.append(common + ["kendall_tau_b", "tau_b", "NA", "z", "NA",
                                      "NA", "NA", "NA", "NA", "NA",
                                      "; ".join(base_note + ["n < 3; test not run"])])
            else:
                note = base_note + ([kd["note"]] if kd["note"] else [])
                note.append("tie-corrected normal approximation")
                rows.append(common + [
                    "kendall_tau_b", "tau_b", fnum(kd["tau_b"], 4),
                    "z", fnum(kd["z"], 4), "NA", fp(kd["p"]),
                    str(kd["concordant"]), str(kd["discordant"]),
                    reading_of(kd["tau_b"]), "; ".join(note)])
    return rows, keep


STRATA_HEADER = [
    "cohort", "metric", "metric_scale", "row_type", "stratum", "stratum_range",
    "n", "n_zero_metric", "median", "q1", "q3", "iqr", "mean", "min", "max",
    "cd4_median_in_stratum", "H", "df", "p_chi2_upper_tail", "strata_n", "note",
]


def strata_rows(cohorts, args, breaks, labels):
    rows = []
    for cohort_name, recs in cohorts:
        buckets = [[] for _ in labels]
        for rec in recs:
            buckets[stratum_index(rec["cd4"], breaks)].append(rec)
        for metric in METRICS:
            per_stratum = []
            for i, (label, rng) in enumerate(labels):
                vals = [r[metric] for r in buckets[i] if r.get(metric) is not None]
                per_stratum.append(vals)
                cd4s = [r["cd4"] for r in buckets[i]]
                note = "" if vals else "no sample in this CD4 band"
                if vals and len(vals) < 3:
                    note = "n=%d; median and IQR are unstable" % len(vals)
                rows.append([
                    cohort_name, metric, SCALE_OF[metric], "stratum_summary",
                    label, rng, str(len(vals)),
                    str(sum(1 for v in vals if float(v) == 0.0)),
                    fnum(median(vals), 4), fnum(quantile(vals, 0.25), 4),
                    fnum(quantile(vals, 0.75), 4),
                    fnum(None if not vals else quantile(vals, 0.75)
                         - quantile(vals, 0.25), 4),
                    fnum(mean(vals), 4),
                    fnum(min(vals) if vals else None, 4),
                    fnum(max(vals) if vals else None, 4),
                    fnum(median(cd4s), 1),
                    "NA", "NA", "NA", "NA", note,
                ])
            kw = kruskal_wallis(per_stratum)
            strata_n = "/".join(str(len(v)) for v in per_stratum)
            if kw is None:
                rows.append([cohort_name, metric, SCALE_OF[metric],
                             "kruskal_wallis", "ALL_STRATA",
                             "|".join(l for l, _r in labels)] + ["NA"] * 10
                            + ["NA", "NA", "NA", strata_n,
                               "fewer than two non-empty CD4 bands; test not run"])
                continue
            note = ["tie-corrected H over %d non-empty bands, n=%d"
                    % (kw["n_groups_used"], kw["n_total"])]
            if min(c for c in kw["counts"]) < 5:
                note.append("a band has n<5; the chi-square p is approximate")
            rows.append([
                cohort_name, metric, SCALE_OF[metric], "kruskal_wallis",
                "ALL_STRATA", "|".join(l for l, _r in labels),
                str(kw["n_total"]),
                str(sum(1 for v in per_stratum for x in v if float(x) == 0.0)),
                "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                fnum(kw["H"], 4), str(kw["df"]), fp(kw["p"]), strata_n,
                "; ".join(note),
            ])
    return rows


CONTRAST_HEADER = [
    "cohort", "outcome", "group1", "n1", "group2", "n2",
    "median_group1", "median_group2", "mean_group1", "mean_group2",
    "U_group1", "U_group2", "z", "p_two_sided_normal",
    "effect_rank_biserial_group1_vs_group2", "note",
]


def contrast_rows(cohorts, args):
    rows = []
    results = []
    thr = args.richness_threshold
    g1 = "richness_ge_%d" % thr
    g2 = "richness_lt_%d" % thr
    for cohort_name, recs in cohorts:
        hi = [r["cd4"] for r in recs
              if r.get(RICHNESS_METRIC) is not None and r[RICHNESS_METRIC] >= thr]
        lo = [r["cd4"] for r in recs
              if r.get(RICHNESS_METRIC) is not None and r[RICHNESS_METRIC] < thr]
        res = mann_whitney_u(hi, lo)
        if res is None:
            rows.append([cohort_name, "cd4_cells_ul", g1, str(len(hi)), g2,
                         str(len(lo))] + ["NA"] * 9
                        + ["one or both richness classes are empty; test not run"])
            results.append((cohort_name, None))
            continue
        note = ["two-sided normal approximation, tie- and continuity-corrected",
                "same implementation as a7; positive effect means %s has the "
                "higher CD4" % g1,
                "the richness class comes from the same anellovirus data, so this "
                "is a re-expression of the correlation, not independent evidence"]
        if min(res["n1"], res["n2"]) < 5:
            note.append("n<5 in one class, p is approximate")
        rows.append([
            cohort_name, "cd4_cells_ul", g1, str(res["n1"]), g2, str(res["n2"]),
            fnum(res["median1"], 1), fnum(res["median2"], 1),
            fnum(res["mean1"], 1), fnum(res["mean2"], 1),
            fnum(res["U1"], 1), fnum(res["U2"], 1), fnum(res["z"], 4),
            fp(res["p"]), fnum(res["effect_r"], 4), "; ".join(note),
        ])
        results.append((cohort_name, res))
    return rows, results


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    ap = argparse.ArgumentParser(
        description="Anellovirus burden against CD4 count: Spearman, Kendall "
                    "tau-b, clinical CD4 strata with Kruskal-Wallis, and a "
                    "Mann-Whitney contrast on the a7 richness cut. Standard "
                    "library only, no figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--indir", default=DEFAULT_INDIR,
                    help="directory holding the a7 outputs (env SUITE_OUTDIR)")
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR,
                    help="output directory (env A9_OUTDIR; created if absent)")
    ap.add_argument("--cd4", default=DEFAULT_CD4,
                    help="user-supplied CD4 table, tsv or csv (env CD4_TABLE). "
                         "Without it only cd4_input_template.tsv is written.")
    ap.add_argument("--burden", default="",
                    help="a7 burden table; default <indir>/<a7-prefix>"
                         "_anellovirus_burden.tsv")
    ap.add_argument("--key", default="",
                    help="a7 sample key; default <indir>/<a7-prefix>_sample_key.tsv")
    ap.add_argument("--a7-prefix", default=DEFAULT_A7_PREFIX,
                    help="filename prefix a7 was run with")
    ap.add_argument("--prefix", default="",
                    help="optional prefix for every file this module writes")
    ap.add_argument("--cohort-groups", default="HIV,HL",
                    help="group labels that make up the WGS cohort")
    ap.add_argument("--focus-group", default="HIV",
                    help="group analysed on its own (the informative analysis)")
    ap.add_argument("--strata-breaks", default="200,500",
                    help="ascending CD4 cut points in cells/uL")
    ap.add_argument("--richness-threshold", type=int, default=3,
                    help="anello_richness_human cut for the CD4 contrast")
    ap.add_argument("--min-n", type=int, default=10,
                    help="n below which a correlation is flagged descriptive only")
    ap.add_argument("--no-carry-extra", action="store_true",
                    help="do not carry extra clinical columns into the joined table")
    return ap


def parse_breaks(text):
    out = []
    for token in str(text).split(","):
        token = token.strip()
        if not token:
            continue
        try:
            out.append(float(token))
        except ValueError:
            return None
    if not out or out != sorted(out) or len(set(out)) != len(out):
        return None
    return out


CAVEATS = [
    ["The outcome is mostly zero, so a rank correlation over it is weak",
     "evidence: read n_zero_metric and the per-stratum medians beside",
     "every rho."],
    ["CD4 and HIV status are collinear in this cohort, so the pooled WGS",
     "rows largely restate the group difference a7 already reported. The",
     "HIV-only rows are the informative ones."],
    ["This is a cross-sectional association and cannot establish the",
     "direction of causation, nor rule out ART, time since diagnosis or",
     "sequencing depth as the driver."],
]


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = build_parser().parse_args(argv)

    breaks = parse_breaks(args.strata_breaks)
    if breaks is None:
        print("WARN: --strata-breaks %r is not an ascending list of distinct "
              "numbers; using 200,500" % args.strata_breaks)
        breaks = [200.0, 500.0]
    labels = strata_labels(breaks)
    if args.richness_threshold < 1:
        print("WARN: --richness-threshold must be >= 1, got %d; using 3"
              % args.richness_threshold)
        args.richness_threshold = 3

    burden_path = args.burden or os.path.join(
        args.indir, args.a7_prefix + "_anellovirus_burden.tsv")
    key_path = args.key or os.path.join(
        args.indir, args.a7_prefix + "_sample_key.tsv")

    if not os.path.exists(burden_path) and not os.path.isdir(args.indir):
        warn("a7 output directory", args.indir)
        return 0
    burden = load_burden(burden_path)
    if burden is None:
        return 0
    if not burden:
        print("WARN: %s has no sample rows; nothing to correlate"
              % os.path.basename(burden_path))
        return 0

    anon_to_real, real_to_anon, key_problem = load_sample_key(key_path)
    if key_problem:
        print("NOTE: %s; only CD4 rows that already carry anonymous ids "
              "(S01..Snn) can be joined" % key_problem)

    cohort_groups = [g.strip() for g in args.cohort_groups.split(",") if g.strip()]
    if not cohort_groups:
        cohort_groups = ["HIV", "HL"]
    for rec in burden:
        real = anon_to_real.get(rec["anon"], "")
        from_name = group_of(real, rec["run"]) if real else "NA"
        rec["group_final"] = from_name if from_name != "NA" else (rec["group"] or "NA")
    wgs = [r for r in burden if r["group_final"] in cohort_groups]
    if not wgs:
        print("WARN: no sample in %s carries a group label in {%s}; nothing to "
              "correlate" % (os.path.basename(burden_path), ",".join(cohort_groups)))
        return 0

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            warn("writable output directory", args.outdir)
            return 0

    def path_of(base):
        return os.path.join(args.outdir, out_name(args.prefix, base))

    template_path = path_of("cd4_input_template.tsv")
    corr_path = path_of("cd4_anello_correlation.tsv")
    strata_path = path_of("cd4_anello_strata.tsv")
    contrast_path = path_of("cd4_anello_richness_contrast.tsv")
    joined_path = path_of("cd4_anello_joined.tsv")
    key_out_path = os.path.join(
        args.outdir, out_name(args.prefix or KEY_STEM, "sample_key.tsv"))
    burden_base = os.path.basename(burden_path)

    def emit_template(reason):
        write_template(template_path, wgs, args, burden_base)
        print("")
        print("a9_cd4_correlation  %s" % TODAY)
        print("%s, so no correlation was computed." % reason)
        print("wrote %s (%d anonymised %s samples, blank cd4 column)"
              % (template_path, len(wgs), "/".join(cohort_groups)))
        print("")
        print("HOW TO FILL IT IN")
        print("  1. open %s" % os.path.basename(template_path))
        print("  2. map sample_id (S01..Snn) back to the real sample with %s"
              % os.path.basename(key_path))
        print("  3. put the absolute CD4 count in cells/uL in the cd4 column. A")
        print("     plain number is safest, '<20' is read as 20, blank means")
        print("     unknown; a thousands separator is only safe in this tab-")
        print("     separated file or in a quoted CSV field. viral_load /")
        print("     treatment / notes are optional and are carried through")
        print("     uninterpreted.")
        print("  4. rerun: python3 %s --indir %s --cd4 <filled file>"
              % (SCRIPT, args.indir))
        print("  You may keep real names in your own copy and pass that file with")
        print("  --cd4: the join tolerates them and nothing but %s will echo them."
              % os.path.basename(key_out_path))
        return 0

    if not args.cd4:
        return emit_template("No --cd4 table was given")
    table = load_cd4_table(args.cd4)
    if table is None:
        return emit_template("The CD4 table could not be used")

    burden_by_anon = dict((r["anon"], r) for r in burden)
    anon_set = set(burden_by_anon)
    exact, normed, loose, clashes = build_lookup(real_to_anon)
    records, report = join_cd4(table, burden_by_anon, exact, normed, loose,
                               anon_set, real_to_anon)
    if not records:
        print("WARN: none of the %d CD4 rows matched a sample in %s"
              % (report["rows_read"], burden_base))
        return emit_template("No CD4 row could be joined")

    if args.no_carry_extra:
        extra_cols = []
    else:
        extra_cols = [col for _i, col, _raw in table["extras"]]

    in_cohort = [r for r in records if r["group_final"] in cohort_groups]
    focus = [r for r in records if r["group_final"] == args.focus_group]
    cohorts = [("WGS_" + "_".join(cohort_groups), in_cohort),
               (args.focus_group + "_ONLY", focus)]
    cohorts = [(name, recs) for name, recs in cohorts if recs]
    if not cohorts:
        print("WARN: %d CD4 rows joined but none of them is in {%s}; nothing to "
              "correlate" % (len(records), ",".join(cohort_groups)))
        return 0

    key_rows = 0
    key_written = ""
    if report["real_identifier_rows"] > 0:
        key_rows = write_sample_key(key_out_path, records, anon_to_real)
        if key_rows:
            key_written = os.path.basename(key_out_path)

    corr, sp_by_key = correlation_rows(cohorts, args)
    write_tsv(corr_path, common_comments(
        args, key_written, burden_base,
        ["Spearman rho over tie-corrected average ranks, p from the t "
         "approximation on n-2 df; Kendall tau-b p from the tie-corrected "
         "normal approximation",
         "every statistic is reported with the n it was computed on"]),
        CORR_HEADER, corr)

    strata = strata_rows(cohorts, args, breaks, labels)
    write_tsv(strata_path, common_comments(
        args, key_written, burden_base,
        ["CD4 bands are half-open intervals in cells/uL: %s"
         % ", ".join("%s = %s" % (l, r) for l, r in labels),
         "q1 / q3 are linear-interpolation quantiles (the numpy default)",
         "Kruskal-Wallis H is tie-corrected; p is the chi-square upper tail on "
         "k-1 df"]), STRATA_HEADER, strata)

    contrast, contrast_res = contrast_rows(cohorts, args)
    write_tsv(contrast_path, common_comments(
        args, key_written, burden_base,
        ["outcome is CD4; the groups are the a7 richness classes",
         "Mann-Whitney U is the same standard-library implementation as a7"]),
        CONTRAST_HEADER, contrast)

    write_joined(joined_path, records, args, key_written, burden_base,
                 extra_cols, breaks, labels)

    # ---------------------------- stdout summary --------------------------- #
    print("")
    print("a9_cd4_correlation  %s" % TODAY)
    print("a7 burden      : %s (%d samples, %d labelled %s)"
          % (burden_base, len(burden), len(wgs), "/".join(cohort_groups)))
    print("cd4 table      : %s" % table["path"])
    print("                 delimiter %s (%s)"
          % (DELIM_NAMES.get(table["delimiter"], repr(table["delimiter"])),
             table["delimiter_how"]))
    print("                 id column   '%s' (%s)"
          % (to_ascii(table["sample_col"]), table["sample_how"]))
    print("                 cd4 column  '%s' (%s)"
          % (to_ascii(table["cd4_col"]), table["cd4_how"]))
    if table["ragged"]:
        print("                 WARNING: %d row(s) do not have the same number of "
              "fields as the header; an unquoted %s inside a value (a thousands "
              "separator?) shifts every column to its right"
              % (table["ragged"],
                 DELIM_NAMES.get(table["delimiter"], "separator")))
    print("join           : %d/%d CD4 rows joined to a burden sample"
          % (report["matched"], report["rows_read"]))
    print("                 id normalisation: trim, collapse internal whitespace, "
          "lower-case; then alphanumeric-only as a last resort")
    for how in sorted(report["match_methods"]):
        print("                 %-32s %d joined rows"
              % (how, report["match_methods"][how]))
    if report["unmatched_known_anon"]:
        print("                 %d CD4 row(s) name a sample that is in the a7 key "
              "but not in the burden table: %s"
              % (len(report["unmatched_known_anon"]),
                 ", ".join(sorted(report["unmatched_known_anon"]))))
    if report["unmatched_unknown_count"]:
        print("                 %d CD4 row(s) matched nothing in the a7 key; they "
              "are counted only, never printed" % report["unmatched_unknown_count"])
    for field, text in (("blank_id", "row(s) with a blank identifier"),
                        ("bad_cd4", "row(s) with an unparsable / empty CD4 value"),
                        ("negative_cd4", "row(s) with a negative CD4 value, dropped"),
                        ("censored_cd4", "censored CD4 value(s) ('<20'), kept as "
                                         "the bare number"),
                        ("implausible_cd4", "CD4 value(s) above 5000 cells/uL, kept"),
                        ("redacted_extra_values", "carried value(s) replaced by an "
                                                  "anonymous id")):
        if report[field]:
            print("                 %d %s" % (report[field], text))
    if report["duplicate_anon"]:
        print("                 duplicate CD4 rows for %s; the first was kept"
              % ", ".join(sorted(report["duplicate_anon"])))
    if clashes:
        print("                 %d real name(s) collide once normalised and were "
              "excluded from the loose join" % len(clashes))
    if table["dropped"]:
        print("dropped columns: %s (identifier- or date-like; never carried)"
              % ", ".join(sorted(set(table["dropped"]))))
    if extra_cols:
        print("carried columns: %s (uninterpreted)" % ", ".join(extra_cols))
    elif args.no_carry_extra:
        print("carried columns: none (--no-carry-extra)")

    print("")
    print("(1) rank correlation of CD4 with anellovirus burden")
    for cohort_name, recs in cohorts:
        print("    %s (n=%d)" % (cohort_name, len(recs)))
        for metric in METRICS:
            sp = sp_by_key.get((cohort_name, metric))
            if sp is None:
                print("      %-26s n<3, not tested" % metric)
                continue
            print("      %-26s rho=%7s  n=%2d  p=%-10s (%d/%d samples zero)"
                  % (metric, fnum(sp["rho"], 3), sp["n"], fp(sp["p"]),
                     sum(1 for r in recs
                         if r.get(metric) is not None and float(r[metric]) == 0.0),
                     sp["n"]))

    print("")
    print("(2) CD4 strata (%s), median %s per band"
          % (", ".join(l for l, _r in labels), RICHNESS_METRIC))
    for cohort_name, recs in cohorts:
        buckets = [[] for _ in labels]
        for rec in recs:
            buckets[stratum_index(rec["cd4"], breaks)].append(rec)
        cells = []
        for i, (label, _rng) in enumerate(labels):
            vals = [r[RICHNESS_METRIC] for r in buckets[i]
                    if r.get(RICHNESS_METRIC) is not None]
            cells.append("%s n=%d med=%s" % (label, len(vals), fnum(median(vals), 1)))
        print("    %-14s %s" % (cohort_name, " | ".join(cells)))
        kw = kruskal_wallis([[r[RICHNESS_METRIC] for r in b
                              if r.get(RICHNESS_METRIC) is not None]
                             for b in buckets])
        if kw is None:
            print("    %-14s Kruskal-Wallis not run (fewer than two non-empty bands)"
                  % "")
        else:
            print("    %-14s Kruskal-Wallis H=%s df=%d p=%s over n=%d (%s)"
                  % ("", fnum(kw["H"], 3), kw["df"], fp(kw["p"]), kw["n_total"],
                     "/".join(str(c) for c in kw["counts"])))

    print("")
    print("(3) CD4 by anellovirus richness class (threshold %d)"
          % args.richness_threshold)
    for cohort_name, res in contrast_res:
        if res is None:
            print("    %-14s not tested (a richness class is empty)" % cohort_name)
            continue
        print("    %-14s CD4 median %s (n=%d, richness>=%d) vs %s (n=%d, "
              "richness<%d) | U=%s p=%s"
              % (cohort_name, fnum(res["median1"], 1), res["n1"],
                 args.richness_threshold, fnum(res["median2"], 1), res["n2"],
                 args.richness_threshold, fnum(res["U1"], 1), fp(res["p"])))

    print("")
    print("READ THIS BEFORE QUOTING ANY NUMBER ABOVE")
    for block in CAVEATS:
        for i, line in enumerate(block):
            print(("  - " if i == 0 else "    ") + line)

    print("")
    print("wrote:")
    for path in (corr_path, strata_path, contrast_path, joined_path):
        print("  %s" % path)
    if key_rows:
        print("  %s" % key_out_path)
        print("REMINDER: %s contains real sample identifiers - do not commit or "
              "email it." % os.path.basename(key_out_path))
    else:
        print("  (no sample key written: no real sample identifier was handled)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
