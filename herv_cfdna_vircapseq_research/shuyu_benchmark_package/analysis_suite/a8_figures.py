#!/usr/bin/env python3
"""a8_figures.py -- figure renderer for the viral-sequencing analysis suite.

WHAT IT DOES
  Reads the .tsv outputs produced by the other suite modules (a1..a7) and draws
  one figure per input table. It recomputes nothing from BAM files and never
  calls samtools; it only reshapes numbers that already exist in the .tsv files.
  This module DOES make figures, so it imports matplotlib (Agg backend, no
  display needed). Everything else is Python standard library only. No network.

  It runs on the cluster or on a laptop -- point --indir at wherever the .tsv
  files were copied.

INPUTS LOOKED FOR IN --indir (each one is optional; a missing file prints
"WARN: ... missing at ..., skipping" and the module carries on, exit code 0)
  detection_threshold_sweep.tsv        -> fig_detection_threshold_sweep
  reference_comparison_by_category.tsv -> fig_reference_comparison_by_category
  depth_sensitivity_by_category.tsv    -> fig_depth_sensitivity_by_category
  refprofile_bins.tsv                  -> fig_refprofile_coverage_tracks
  anellovirus_burden.tsv               -> fig_anellovirus_burden
  coinfection_pairs.tsv                -> fig_coinfection_pairs

  A prefixed upstream filename is picked up too: if the exact name is absent,
  any file in --indir ending with it is used (so a7_virome_coinfection_pairs.tsv
  satisfies coinfection_pairs.tsv). Each path can also be given explicitly with
  --sweep-tsv / --refcmp-tsv / --depth-tsv / --profile-tsv / --anello-tsv /
  --coinf-tsv.

  Column names are resolved case-insensitively from lists of aliases, so the
  upstream modules can name a column "threshold" or "min_reads", "reads" or
  "n_reads", "sample" or "sample_anon", etc. The observed suite schemas are
  handled directly: reference_id / ref_label / sample_anon / bin_start /
  mean_depth (bin profiles), anello_rpm_human / anello_richness_human (virome),
  virus_group_a / virus_group_b / n_both / n_a / n_b (co-infection pairs),
  reads_5m / reads_full / n_pos_5m / n_pos_full (depth sensitivity). Both wide
  (viral_only_reads / hg38_reads) and long (reference + reads) shapes work for
  the reference and depth tables. When a table is stratified by cohort or group,
  the pooled "ALL" rows are used if present (stated on the figure), otherwise
  the stratum is folded into the row label.

WHAT IT WRITES (into --outdir)
  fig_<name>.png                       dpi 200 (--dpi)
  fig_<name>.svg                       same figure, vector
  <prefix>_figure_index.tsv            one row per figure: source, rows, status
  <prefix>_sample_key.tsv              anon_sample, real_sample, group.
                                       Written ONLY if a real name was seen in
                                       an input. This is the one file that
                                       contains identifiers; it carries a
                                       "DO NOT COMMIT OR EMAIL" header comment.

IDENTIFIER RULE
  No figure and no .tsv other than the sample key ever contains a real sample
  name. Any sample label that is not already of the form S<digits> is mapped to
  S01..Snn, ordered by the sorted real name (numbering continues past any
  already-anonymous S<n> labels found in the same inputs). Group labels
  (HIV / HL / TCL / NA) come from a group column when present, otherwise they
  are derived from the real sample name before that name is discarded:
  "_HIV" -> HIV, "_HL" -> HL, "targeted_htlv" or "TCL" -> TCL, else NA.

STYLE
  Light surface, ink/secondary/muted text, a 5-slot categorical palette and a
  4-stop sequential blue for magnitude. Bar charts are never drawn on a log
  axis -- log panels use dots / lollipops. Every panel carries direct value
  labels, so colour is never the only encoding.

EXAMPLE
  python3 a8_figures.py \
      --indir /path/to/runs/panel_report_20260725/suite_out \
      --outdir /path/to/runs/panel_report_20260725/suite_out/figures \
      --prefix a8 --default-threshold 100

Date: 2026-07-26
"""

from __future__ import annotations

import argparse
import os
import random
import re
import statistics
import sys

# ----------------------------------------------------------------- matplotlib
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    HAVE_MPL = True
    MPL_ERR = ""
except Exception as _exc:                                   # pragma: no cover
    HAVE_MPL = False
    MPL_ERR = str(_exc)

# ----------------------------------------------------------------- defaults
DEF_INDIR = "/path/to/runs/panel_report_20260725/suite_out"
DEF_OUTDIR = "/path/to/runs/panel_report_20260725/suite_out/figures"
DEF_PREFIX = "a8"
DEF_THRESHOLD = 100.0
TODAY = "2026-07-26"

# ----------------------------------------------------------------- palette
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"]
SEQ = ["#cde2fb", "#86b6ef", "#2a78d6", "#1c5cab"]

# values that mean "this row pools every stratum"
ALL_VALUES = ("ALL", "ALL_GROUPS", "ALL_SAMPLES", "ALL_CATEGORIES", "TOTAL",
              "OVERALL", "COMBINED")

GROUP_ORDER = ["HIV", "HL", "TCL", "NA"]
GROUP_COLOR = {"HIV": CAT[0], "HL": CAT[1], "TCL": CAT[2], "NA": MUTED}

# named panel references, for readable axis labels
NAMED_REFS = {
    "SHUYU_000096_NC_007605.1": "EBV type 1",
    "SHUYU_000101_NC_009334.1": "EBV type 2",
    "SHUYU_000001_NC_000898.1": "HHV-6B",
    "SHUYU_000080_NC_006273.2": "CMV",
    "SHUYU_000054_NC_002076.2": "TTV1",
}
ACC_NAMES = {}
for _rid, _nm in NAMED_REFS.items():
    ACC_NAMES[_rid.split("_", 2)[-1]] = _nm

# input table registry: key, canonical filename, human label, alternate
# filename suffixes (upstream modules prefix their outputs with --prefix; the
# canonical names here are what the default prefixes produce)
INPUTS = [
    ("sweep", "detection_threshold_sweep.tsv", "detection threshold sweep table",
     ["threshold_sweep.tsv"]),
    ("refcmp", "reference_comparison_by_category.tsv",
     "reference comparison table", []),
    ("depth", "depth_sensitivity_by_category.tsv", "depth sensitivity table",
     ["_by_category.tsv"]),
    ("profile", "refprofile_bins.tsv", "reference coverage bin table",
     ["_bins.tsv"]),
    ("anello", "anellovirus_burden.tsv", "anellovirus burden table", []),
    ("coinf", "coinfection_pairs.tsv", "co-infection pair table", []),
]
PRIMARY_NAMES = [row[1] for row in INPUTS]


def style_rcparams():
    plt.rcParams.update({
        "font.family": ["DejaVu Sans", "Segoe UI", "sans-serif"],
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS,
        "text.color": INK,
        "axes.labelcolor": INK2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "svg.fonttype": "none",
    })


# ===================================================================== io
def read_tsv(path, label):
    """Read a tab-separated file with a header row. Returns (header, rows) or
    (None, None) if the file is missing or unusable. Comment lines starting
    with '#' and blank lines are ignored."""
    if not os.path.isfile(path):
        print("WARN: %s missing at %s, skipping" % (label, path))
        return None, None
    header = None
    rows = []
    try:
        fh = open(path, "r", encoding="utf-8-sig", errors="replace")
    except OSError as exc:
        print("WARN: %s unreadable at %s (%s), skipping" % (label, path, exc))
        return None, None
    with fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                continue
            if len(parts) < len(header):
                parts = parts + [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, [p.strip() for p in parts])))
    if header is None:
        print("WARN: %s empty at %s, skipping" % (label, path))
        return None, None
    if not rows:
        print("WARN: %s has a header but no data rows at %s, skipping"
              % (label, path))
        return None, None
    return header, rows


def pick(header, cands, fuzzy=True):
    """Resolve a column name from a list of aliases (case-insensitive)."""
    if not header:
        return None
    low = {}
    for h in header:
        low.setdefault(h.strip().lower(), h)
    for c in cands:
        if c in low:
            return low[c]
    if fuzzy:
        for c in cands:
            for h in header:
                if c in h.strip().lower():
                    return h
    return None


def num(value):
    """Tolerant numeric parse. Returns float or None."""
    if value is None:
        return None
    s = str(value).strip().replace(",", "")
    if s.endswith("%"):
        s = s[:-1]
    if s in ("", "-", ".", "NA", "na", "N/A", "n/a", "None", "none",
             "nan", "NaN", "NULL", "null"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def fmt_int(v):
    if v is None:
        return "NA"
    return "{:,}".format(int(round(v)))


def fmt_val(v):
    """Integer-looking values as integers, otherwise 2-3 decimals."""
    if v is None:
        return "NA"
    if abs(v - round(v)) < 1e-9 and abs(v) < 1e15:
        return "{:,}".format(int(round(v)))
    if abs(v) >= 100:
        return "{:,.0f}".format(v)
    if abs(v) >= 1:
        return "%.2f" % v
    return "%.3f" % v


def shorten(text, maxlen=26):
    t = str(text).strip()
    return t if len(t) <= maxlen else t[:maxlen - 3] + "..."


def fold(a, b):
    """Fold difference between two non-negative values, floored at 1 read."""
    hi, lo = max(a, b), min(a, b)
    if hi <= 0:
        return 0.0
    return hi / max(lo, 1.0)


def fmt_pct(v):
    if v is None:
        return "NA"
    return "%.1f%%" % v


def thousands(v, _pos=None):
    try:
        return "{:,}".format(int(round(v)))
    except (ValueError, OverflowError):
        return str(v)


def ref_label(rid, maxlen=30):
    s = str(rid).strip()
    if s in NAMED_REFS:
        return NAMED_REFS[s]
    bare = re.sub(r"^SHUYU_\d+_", "", s)
    if bare in ACC_NAMES:
        return ACC_NAMES[bare]
    if len(bare) > maxlen:
        bare = bare[:maxlen - 3] + "..."
    return bare


def group_of(name):
    """Group label derived from the REAL sample name.

    Suite-wide rule, matched case-insensitively so it agrees with a1..a7:
    "_HIV" -> HIV, "_HL" -> HL, "TCL"/"targeted_htlv" -> TCL, else NA.
    """
    if not name:
        return "NA"
    up = str(name).upper()
    if "_HIV" in up:
        return "HIV"
    if "_HL" in up:
        return "HL"
    if "TARGETED_HTLV" in up or "TCL" in up:
        return "TCL"
    return "NA"


def norm_group(value, fallback_name=None):
    """Group from an explicit column value, falling back to the sample name."""
    s = "" if value is None else str(value).strip()
    up = s.upper()
    if up in ("HIV", "HIV1", "HIVPOS", "HIV+"):
        return "HIV"
    if up in ("HL", "HODGKIN", "HODGKIN_LYMPHOMA"):
        return "HL"
    if up in ("TCL", "HTLV", "TARGETED_HTLV", "T_CELL_LYMPHOMA"):
        return "TCL"
    if s:
        derived = group_of(s)
        if derived != "NA":
            return derived
    return group_of(fallback_name) if fallback_name else "NA"


# ============================================================ anonymisation
ANON_RE = re.compile(r"^S\d+$")


class Anon(object):
    """Maps real sample names to S01..Snn (ordered by the sorted real name).
    Labels that already look anonymous (S + digits) pass straight through and
    are never written to the key file."""

    def __init__(self):
        self._seen = set()
        self._map = {}
        self._passthrough = set()

    def observe(self, name):
        if name is None:
            return
        s = str(name).strip()
        if not s:
            return
        if ANON_RE.match(s):
            self._passthrough.add(s)
        else:
            self._seen.add(s)

    def finalize(self):
        start = 0
        for s in self._passthrough:
            try:
                start = max(start, int(s[1:]))
            except ValueError:
                pass
        total = start + len(self._seen)
        width = 2 if total <= 99 else len(str(total))
        for i, real in enumerate(sorted(self._seen), start=start + 1):
            self._map[real] = "S" + str(i).zfill(width)

    def get(self, name):
        if name is None:
            return "NA"
        s = str(name).strip()
        if not s:
            return "NA"
        if s in self._map:
            return self._map[s]
        if ANON_RE.match(s):
            return s
        return "S??"                 # never emit an unmapped real name

    @property
    def n_mapped(self):
        return len(self._map)

    def write_key(self, path):
        if not self._map:
            return None
        lines = [
            "# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
            "# written by a8_figures.py on %s" % TODAY,
            "# anon ids continue after any pre-anonymised S<n> labels already "
            "present in the inputs",
            "\t".join(["anon_sample", "real_sample", "group"]),
        ]
        for real in sorted(self._map):
            lines.append("\t".join([self._map[real], real, group_of(real)]))
        with open(path, "w", encoding="ascii", newline="") as fh:
            fh.write("\n".join(lines) + "\n")
        return path


SAMPLE_ALIASES = ["sample", "sample_id", "sample_name", "anon_sample",
                  "sample_anon", "sampleid", "library", "specimen"]
# deliberately strict: must not match aggregate columns such as "samples_full"
SAMPLE_RE = re.compile(r"^(anon[_-]?)?sample([_-]?(id|name|label))?$")


def sample_col(header):
    """The per-sample label column of a table, or None."""
    col = pick(header, SAMPLE_ALIASES, fuzzy=False)
    if col is not None:
        return col
    for h in header or []:
        if SAMPLE_RE.match(h.strip().lower()):
            return h
    return None


def prescan_samples(tables, anon):
    """Feed every sample label from every loaded table into the mapper, so the
    numbering is stable across figures."""
    for key in tables:
        got = tables.get(key)
        if not got:
            continue
        header, rows = got
        col = sample_col(header)
        if col is None:
            continue
        for r in rows:
            anon.observe(r.get(col))
    anon.finalize()


# ================================================================= plotting
def blue_for(value, vmax):
    if vmax is None or vmax <= 0 or value is None:
        return SEQ[0]
    r = float(value) / float(vmax)
    if r < 0.25:
        return SEQ[0]
    if r < 0.50:
        return SEQ[1]
    if r < 0.80:
        return SEQ[2]
    return SEQ[3]


def blue_cmap():
    return LinearSegmentedColormap.from_list("suite_blue", SEQ)


def clean_axes(ax, xgrid=True, ygrid=False, hide_left=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if hide_left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    if xgrid:
        ax.xaxis.grid(True, color=GRID, lw=1)
    if ygrid:
        ax.yaxis.grid(True, color=GRID, lw=1)
    ax.set_axisbelow(True)


def try_symlog(ax, axis="x", linthresh=1.0):
    """symlog if this matplotlib accepts it, otherwise leave the axis linear."""
    try:
        if axis == "x":
            ax.set_xscale("symlog", linthresh=linthresh)
        else:
            ax.set_yscale("symlog", linthresh=linthresh)
        return True
    except Exception:
        return False


def save_fig(fig, outdir, stem, dpi):
    png = os.path.join(outdir, stem + ".png")
    svg = os.path.join(outdir, stem + ".svg")
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png, svg


def sup(fig, sentence, sub=None):
    """One clear sentence as the title, plus a small method note under it."""
    fig.suptitle(sentence, x=0.02, y=1.075, ha="left", va="bottom",
                 fontsize=15.5, fontweight="bold", color=INK)
    if sub:
        fig.text(0.02, 1.055, sub, fontsize=9.2, color=MUTED, ha="left",
                 va="top")


# ---------------------------------------------------------------- figure 1
def fig_sweep(header, rows, args):
    tcol = pick(header, ["threshold", "min_reads", "read_threshold",
                         "reads_threshold", "cutoff", "thresh", "min_read"])
    pcol = pick(header, ["precision", "precision_pct", "ppv"])
    rcol = pick(header, ["recall", "recall_pct", "sensitivity", "tpr"])
    fcol = pick(header, ["f1", "f1_score", "fscore"])
    if tcol is None:
        print("WARN: detection_threshold_sweep.tsv has no threshold column, "
              "skipping")
        return None
    metrics = [(n, c, col) for n, c, col in
               [("Precision", CAT[0], pcol), ("Recall", CAT[1], rcol),
                ("F1", CAT[2], fcol)] if col is not None]
    if not metrics:
        print("WARN: detection_threshold_sweep.tsv has no precision / recall "
              "column, skipping")
        return None

    # a1 emits target x scope x threshold: keep the pooled scope if present
    scol_scope = pick(header, ["scope", "stratum", "subset"], fuzzy=False)
    scope_note = None
    if scol_scope:
        seen = sorted(set(str(r.get(scol_scope, "")).strip() for r in rows))
        seen = [v for v in seen if v]
        pooled = [v for v in seen
                  if v.upper() in ALL_VALUES + ("COMBINED", "ALL_RUNS")]
        if len(seen) > 1:
            keep = pooled[0] if pooled else seen[0]
            rows = [r for r in rows
                    if str(r.get(scol_scope, "")).strip() == keep]
            scope_note = "%s = %s rows only" % (scol_scope, keep)
    gcol = pick(header, ["category", "cohort", "reference", "target", "panel",
                         "group", "run"], fuzzy=False)
    series = {}
    for r in rows:
        t = num(r.get(tcol))
        if t is None:
            continue
        gname = (str(r.get(gcol)).strip() if gcol else "") or "all"
        series.setdefault(gname, []).append(
            (t, dict((n, num(r.get(col))) for n, _c, col in metrics)))
    if not series:
        print("WARN: detection_threshold_sweep.tsv has no numeric thresholds, "
              "skipping")
        return None

    # scale detection: values given as fractions become percent
    allv = [v for pts in series.values() for _t, d in pts
            for v in d.values() if v is not None]
    scale = 100.0 if (allv and max(allv) <= 1.0001) else 1.0

    names = sorted(series, key=lambda k: (-len(series[k]), k))
    dropped = names[5:]
    names = names[:5]
    n = len(names)
    fig, axarr = plt.subplots(1, n, figsize=(4.9 * n + 1.2, 5.2),
                              squeeze=False, sharey=True)
    axes = list(axarr[0])

    head_p = head_r = None
    for ai, gname in enumerate(names):
        ax = axes[ai]
        pts = sorted(series[gname], key=lambda x: x[0])
        xs = [p[0] for p in pts]
        pos = [x for x in xs if x > 0]
        uselog = bool(pos) and (max(pos) / min(pos) >= 50.0)
        for mname, mcol, _src in metrics:
            ys = [(p[1].get(mname) * scale)
                  if p[1].get(mname) is not None else None for p in pts]
            gx = [x for x, y in zip(xs, ys) if y is not None]
            gy = [y for y in ys if y is not None]
            if not gy:
                continue
            ax.plot(gx, gy, color=mcol, lw=2.0, marker="o", markersize=4.5,
                    markeredgecolor=SURFACE, markeredgewidth=0.8,
                    label=mname, zorder=3)
            ax.annotate(fmt_pct(gy[-1]), (gx[-1], gy[-1]),
                        textcoords="offset points", xytext=(6, 0),
                        fontsize=9, color=mcol, va="center")
        # default threshold marker plus its value labels
        thr = args.default_threshold
        ax.axvline(thr, color=INK2, lw=1.2, ls="--", zorder=2)
        near = min(pts, key=lambda p: abs(p[0] - thr)) if pts else None
        if near is not None and abs(near[0] - thr) <= max(1.0, 0.2 * thr):
            bits = []
            for mname, mcol, _src in metrics:
                v = near[1].get(mname)
                if v is None:
                    continue
                v = v * scale
                ax.scatter([near[0]], [v], s=95, color=mcol, zorder=5,
                           edgecolor=SURFACE, linewidth=1.4)
                bits.append("%s %s" % (mname[0], fmt_pct(v)))
                if ai == 0 and mname == "Precision":
                    head_p = v
                if ai == 0 and mname == "Recall":
                    head_r = v
            if bits:
                ax.annotate("default %s reads\n%s"
                            % (fmt_int(thr), "  ".join(bits)),
                            (near[0], 4), textcoords="offset points",
                            xytext=(8, 0), fontsize=9, color=INK2,
                            va="bottom", ha="left")
        if uselog:
            try:
                ax.set_xscale("log")
            except Exception:
                pass
            ax.set_xticks(list(pos))
            ax.xaxis.set_major_formatter(FuncFormatter(thousands))
            ax.tick_params(axis="x", labelrotation=45)
            for lab in ax.get_xticklabels():
                lab.set_horizontalalignment("right")
        ax.set_ylim(0, 105)
        ax.set_xlabel("Detection threshold (reads)", fontsize=10, color=INK2)
        if ai == 0:
            ax.set_ylabel("Percent", fontsize=10, color=INK2)
        ax.set_title(gname if gname != "all" else "All categories",
                     loc="left", fontsize=11.5, fontweight="bold", color=INK)
        clean_axes(ax, xgrid=True, ygrid=True, hide_left=False)
    axes[0].legend(frameon=False, fontsize=10, loc="lower left")

    if head_p is not None and head_r is not None:
        sentence = ("At the default threshold of %s reads, precision is %s "
                    "and recall is %s" % (fmt_int(args.default_threshold),
                                          fmt_pct(head_p), fmt_pct(head_r)))
    else:
        sentence = ("Precision and recall across detection thresholds, with "
                    "the default of %s reads marked"
                    % fmt_int(args.default_threshold))
    subnote = "Dashed line = default threshold."
    if scope_note:
        subnote += "  " + scope_note + "."
    if dropped:
        subnote += ("  Showing the %d largest series; %d not shown."
                    % (len(names), len(dropped)))
    sup(fig, sentence, subnote)
    fig.tight_layout(rect=[0.015, 0.02, 1, 0.94])
    return save_fig(fig, args.outdir, "fig_detection_threshold_sweep",
                    args.dpi)


# ---------------------------------------------------------------- pair logic
VIRAL_WORDS = ("viral_only", "viralonly", "viral-only", "competitive",
               "viral_ref", "full_competitive")
HG38_WORDS = ("hg38", "inclusive", "with_human", "human_inclusive", "refseq",
              "panel", "masked")
FULL_WORDS = ("full", "all_reads", "total_depth", "fulldepth")
SUB_WORDS = ("5m", "5_m", "sub", "subsample", "downsample")


def side_of(text, a_words, b_words):
    low = str(text).strip().lower()
    a = any(w in low for w in a_words)
    b = any(w in low for w in b_words)
    if a and not b:
        return "a"
    if b and not a:
        return "b"
    if a and b:
        ai = min(low.find(w) for w in a_words if w in low)
        bi = min(low.find(w) for w in b_words if w in low)
        return "a" if ai <= bi else "b"
    return None


def metric_of(text):
    low = str(text).strip().lower()
    if "sample" in low or "detect" in low or low.startswith("n_"):
        return "samples"
    return "reads"


SIDE_DESC_STEMS = ("arm", "variant", "reference", "ref", "run", "mode",
                   "build", "label", "name")


def suffix_pairs(header, rows, a_words, b_words, catcol, cohcol):
    """Handle the "<stem>_a / <stem>_b" layout, where which reference is side a
    lives in a descriptor column such as arm_a / arm_b (a2's by_category
    output). Returns (pairs, note, legend) or None."""
    stems = {}
    for h in header:
        low = h.strip().lower()
        if len(low) > 2 and low[-2] == "_" and low[-1] in ("a", "b"):
            stems.setdefault(low[:-2], {})[low[-1]] = h
    stems = dict((k, v) for k, v in stems.items() if "a" in v and "b" in v)
    descs = [k for k in SIDE_DESC_STEMS if k in stems]
    if not stems or not descs:
        return None

    groupcol = pick(header, ["comparison", "contrast", "pair"], fuzzy=False)
    groups = {}
    for r in rows:
        g = str(r.get(groupcol, "")).strip() if groupcol else ""
        groups.setdefault(g, []).append(r)

    best = None                      # (matched, n_rows, group, orient, desc)
    for g in sorted(groups):
        grows = groups[g]
        for d in descs:
            ca, cb = stems[d]["a"], stems[d]["b"]
            va = set(str(r.get(ca, "")).strip() for r in grows)
            vb = set(str(r.get(cb, "")).strip() for r in grows)
            sa = set(side_of(v, a_words, b_words) for v in va if v)
            sb = set(side_of(v, a_words, b_words) for v in vb if v)
            if sa == set(["a"]) and sb == set(["b"]):
                cand = (1, len(grows), g, ("a", "b"), d)
            elif sa == set(["b"]) and sb == set(["a"]):
                cand = (1, len(grows), g, ("b", "a"), d)
            else:
                cand = (0, len(grows), g, ("a", "b"), d)
            if best is None or cand[:2] > best[:2]:
                best = cand
            break                    # first descriptor stem decides
    if best is None:
        return None
    matched, _n, gname, orient, desc = best
    grows = groups[gname]
    bits = []
    if gname:
        bits.append("comparison = " + gname)
    arm = {}
    for side in ("a", "b"):
        vals = sorted(set(str(r.get(stems[desc][side], "")).strip()
                          for r in grows))
        vals = [v for v in vals if v]
        arm[side] = vals[0] if len(vals) == 1 else (vals[0] if vals else side)
    legend = (arm[orient[0]], arm[orient[1]])
    if not matched:
        bits.append("arms taken as given (%s vs %s)" % legend)

    out = {}
    for stem, cols in stems.items():
        if stem in SIDE_DESC_STEMS:
            continue
        ca, cb = cols[orient[0]], cols[orient[1]]
        if not any(num(r.get(ca)) is not None or num(r.get(cb)) is not None
                   for r in grows):
            continue
        metric = metric_of(stem)
        if metric in out:
            continue
        per = {}
        for r in grows:
            lab = str(r.get(catcol, "")).strip()
            if not lab:
                continue
            if cohcol:
                coh = str(r.get(cohcol, "")).strip()
                if coh:
                    lab = "%s  (%s)" % (lab, coh)
            va, vb = num(r.get(ca)), num(r.get(cb))
            if va is None and vb is None:
                continue
            per[lab] = (va or 0.0, vb or 0.0)
        if per:
            out[metric] = per
    if not out:
        return None
    return out, ("; ".join(bits) if bits else None), legend


def extract_pairs(header, rows, a_words, b_words, cat_aliases):
    """Return (pairs, catcol, note) where pairs is
    {metric: {row_label: (a_value, b_value)}}. Wide and long shapes both work.
    If the table is stratified (a cohort or group column) the pooled "ALL" rows
    are used when they exist, otherwise the stratum is folded into the label."""
    catcol = pick(header, cat_aliases)
    if catcol is None:
        return None, None, None, None
    cohcol = pick(header, ["cohort", "group", "run_group", "assay", "arm",
                           "panel_cohort"], fuzzy=False)
    note = None
    if cohcol:
        vals = sorted(set(str(r.get(cohcol, "")).strip() for r in rows))
        vals = [v for v in vals if v]
        pooled = [v for v in vals if v.upper() in ALL_VALUES]
        if len(vals) > 1 and pooled:
            keep = pooled[0]
            rows = [r for r in rows
                    if str(r.get(cohcol, "")).strip() == keep]
            note = "%s = %s rows only" % (cohcol, keep)
            cohcol = None                      # no need to label the stratum
        elif len(vals) <= 1:
            cohcol = None

    def row_label(r):
        lab = str(r.get(catcol, "")).strip()
        if not lab:
            return None
        if cohcol:
            coh = str(r.get(cohcol, "")).strip()
            if coh:
                lab = "%s  (%s)" % (lab, coh)
        return lab

    # "<stem>_a / <stem>_b" shape first: it is the most specific
    got = suffix_pairs(header, rows, a_words, b_words, catcol, cohcol)
    if got is not None:
        out, snote, legend = got
        both = [x for x in (note, snote) if x]
        return out, catcol, ("; ".join(both) if both else None), legend

    # wide shape: both sides live in columns of the same row
    wide = {}
    for h in header:
        if h in (catcol, cohcol):
            continue
        s = side_of(h, a_words, b_words)
        if s is None:
            continue
        wide.setdefault(metric_of(h), {})[s] = h
    wide = dict((m, d) for m, d in wide.items() if "a" in d and "b" in d)
    if wide:
        out = {}
        for metric, cols in wide.items():
            per = {}
            for r in rows:
                lab = row_label(r)
                if lab is None:
                    continue
                va, vb = num(r.get(cols["a"])), num(r.get(cols["b"]))
                if va is None and vb is None:
                    continue
                per[lab] = (va or 0.0, vb or 0.0)
            if per:
                out[metric] = per
        if out:
            return out, catcol, note, None

    # long shape: a side column plus one or more value columns
    sidecol = pick(header, ["reference", "reference_set", "ref", "run", "mode",
                            "build", "depth", "subset", "arm"], fuzzy=False)
    if sidecol is None:
        for h in header:
            vals = set(str(r.get(h, "")).strip() for r in rows)
            vals = set(v for v in vals if v)
            if 2 <= len(vals) <= 6 and all(side_of(v, a_words, b_words)
                                           for v in vals):
                sidecol = h
                break
    if sidecol is None:
        return None, catcol, note, None
    valcols = {}
    for h in header:
        if h in (catcol, cohcol, sidecol):
            continue
        if any(num(r.get(h)) is not None for r in rows):
            valcols.setdefault(metric_of(h), h)
    out = {}
    for metric, vcol in valcols.items():
        per = {}
        for r in rows:
            lab = row_label(r)
            if lab is None:
                continue
            s = side_of(r.get(sidecol), a_words, b_words)
            v = num(r.get(vcol))
            if s is None or v is None:
                continue
            cur = list(per.get(lab, (None, None)))
            cur[0 if s == "a" else 1] = v
            per[lab] = tuple(cur)
        per = dict((k, (v[0] or 0.0, v[1] or 0.0)) for k, v in per.items())
        if per:
            out[metric] = per
    return (out or None), catcol, note, None


def dumbbell_panel(ax, labels, a_vals, b_vals, a_color, b_color, xlabel,
                   title, log=True):
    """Dots joined by a line -- used instead of bars wherever the axis is log
    or symlog."""
    ypos = list(range(len(labels)))[::-1]
    for y, av, bv in zip(ypos, a_vals, b_vals):
        ax.plot([av, bv], [y, y], color=AXIS, lw=2.0, zorder=1,
                solid_capstyle="round")
        ax.scatter([av], [y], s=150, color=a_color, zorder=3,
                   edgecolor=SURFACE, linewidth=1.5)
        ax.scatter([bv], [y], s=78, color=b_color, zorder=4,
                   edgecolor=SURFACE, linewidth=1.2)
        ax.annotate("%s -> %s" % (fmt_val(av), fmt_val(bv)),
                    (max(av, bv), y), textcoords="offset points",
                    xytext=(13, 0), va="center", fontsize=9, color=INK2)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=10, color=INK2)
    ax.set_ylim(-0.8, len(labels) - 0.2)
    hi = max([v for v in list(a_vals) + list(b_vals)] + [1.0])
    # log-shaped headroom and the "symlog" wording are only correct if the
    # scale change actually took; otherwise fall back to a linear axis and say so
    went_log = try_symlog(ax, "x", 1.0) if log else False
    if went_log:
        ax.set_xlim(-0.35, hi * 40)
        ax.xaxis.set_major_formatter(FuncFormatter(thousands))
        ax.tick_params(axis="x", labelrotation=45)
        for lab in ax.get_xticklabels():
            lab.set_horizontalalignment("right")
    else:
        ax.set_xlim(0, hi * 1.45)
        xlabel = xlabel.replace("(symlog scale)", "(linear scale)")
    ax.set_xlabel(xlabel, fontsize=10, color=INK2)
    ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold",
                 color=INK, pad=10)
    clean_axes(ax)


# ---------------------------------------------------------------- figure 2
def headline_excluded(header, rows):
    """Categories whose A-vs-B fold change must not become the headline.

    Two kinds, both still plotted but never quoted as a finding:
      * structural=yes -- a2 marks a category that exists in only one arm's
        reference (HUMAN / HERV / LINE1). "0 -> 24,000,000 reads" there is
        reference content, not a call the stricter arm changed.
      * TOTAL_* -- aggregate rows that are sums of the other rows.
    """
    out = set()
    scol = pick(header, ["structural"], fuzzy=False)
    ccol = pick(header, ["category", "class", "cat", "categories"])
    for r in rows:
        cat = str(r.get(ccol, "")).strip() if ccol else ""
        if not cat:
            continue
        if cat.upper().startswith("TOTAL"):
            out.add(cat)
        if scol and str(r.get(scol, "")).strip().lower() in ("yes", "y", "true", "1"):
            out.add(cat)
    return out


def fig_refcmp(header, rows, args):
    skip_headline = headline_excluded(header, rows)
    pairs, _catcol, note, legend = extract_pairs(
        header, rows, VIRAL_WORDS, HG38_WORDS,
        ["category", "class", "cat", "categories"])
    lab_a, lab_b = legend if legend else ("viral-only reference (competitive)",
                                          "hg38-inclusive reference")
    if not pairs:
        print("WARN: reference_comparison_by_category.tsv has no "
              "viral-only / hg38-inclusive column pair, skipping")
        return None
    metrics = [m for m in ("reads", "samples") if m in pairs]
    n = len(metrics)
    fig, axarr = plt.subplots(1, n, figsize=(7.6 * n, 5.4), squeeze=False)
    # one row order for every panel, taken from the primary metric
    primary = pairs[metrics[0]]
    order = [k for k, _v in sorted(primary.items(), key=lambda kv: -max(kv[1]))]
    for m in metrics[1:]:
        for k in sorted(pairs[m], key=lambda k2: -max(pairs[m][k2])):
            if k not in order:
                order.append(k)
    biggest = None
    for ai, metric in enumerate(metrics):
        per = pairs[metric]
        items = [(k, per[k]) for k in order if k in per]
        labels = [k for k, _v in items]
        av = [v[0] for _k, v in items]
        bv = [v[1] for _k, v in items]
        dumbbell_panel(
            axarr[0][ai], labels, av, bv, SEQ[1], SEQ[3],
            "Reads (symlog scale)" if metric == "reads"
            else "Samples with a call (symlog scale)",
            "%s. %s per category" % ("AB"[ai],
                                     "Reads" if metric == "reads"
                                     else "Samples detected"),
            log=True)
        if metric == "reads" or biggest is None:
            for k, v in items:
                if max(v) < 10:                      # ignore near-zero noise
                    continue
                cat = k.split("  (")[0]
                if cat in skip_headline:             # structural / aggregate row
                    continue
                score = (fold(v[0], v[1]), abs(v[1] - v[0]))
                if biggest is None or score > biggest[0]:
                    biggest = (score, cat, v[0], v[1], metric)
    if biggest:
        _sc, blab, bva, bvb, bmetric = biggest
        unit = "reads" if bmetric == "reads" else "samples"
        sentence = ("%s differs most between the two references: %s %s (%s) "
                    "versus %s (%s), %sx"
                    % (blab, fmt_val(bva), unit, shorten(lab_a),
                       fmt_val(bvb), shorten(lab_b),
                       fmt_val(fold(bva, bvb))))
    else:
        sentence = ("%s versus %s, per category" % (lab_a, lab_b))
    h = [plt.Line2D([], [], marker="o", ls="", markersize=10, color=SEQ[1],
                    markeredgecolor=SURFACE, markeredgewidth=2,
                    label=lab_a),
         plt.Line2D([], [], marker="o", ls="", markersize=9, color=SEQ[3],
                    markeredgecolor=SURFACE, markeredgewidth=2,
                    label=lab_b)]
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.02, -0.055),
               ncol=2, frameon=False, fontsize=10)
    sub = "Dots, not bars, because the axis is symlog; 0 sits on the origin."
    if note:
        sub = sub + "  " + note + "."
    if skip_headline:
        sub = sub + ("  Reference-only and TOTAL_* rows are plotted but "
                     "excluded from the headline.")
    sup(fig, sentence, sub)
    fig.tight_layout(rect=[0.015, 0.04, 1, 0.94])
    return save_fig(fig, args.outdir, "fig_reference_comparison_by_category",
                    args.dpi)


# ---------------------------------------------------------------- figure 3
def fig_depth(header, rows, args):
    pairs, catcol, note, _legend = extract_pairs(
        header, rows, SUB_WORDS, FULL_WORDS,
        ["category", "class", "cat", "categories"])
    retcol = pick(header, ["retention", "retention_pct", "retained_pct",
                           "pct_retained", "recovery", "recovery_pct"])
    if not pairs and retcol is None:
        print("WARN: depth_sensitivity_by_category.tsv has no 5M / full column "
              "pair or retention column, skipping")
        return None

    # headline retention metric: detections first, then reads
    ret = {}
    counts = {}
    src = None
    for metric in ("samples", "reads"):
        if pairs and metric in pairs:
            src = metric
            for lab, (sub, full) in pairs[metric].items():
                if full and full > 0:
                    ret[lab] = 100.0 * sub / full
                    counts[lab] = (sub, full)
                elif not sub:
                    ret[lab] = 0.0
                    counts[lab] = (0.0, 0.0)
            break
    if not ret and retcol is not None and catcol is not None:
        for r in rows:
            lab = str(r.get(catcol, "")).strip()
            v = num(r.get(retcol))
            if lab and v is not None:
                ret[lab] = v * 100.0 if v <= 1.0001 else v
        src = "reported"
    if not ret:
        print("WARN: depth_sensitivity_by_category.tsv has no usable "
              "retention values, skipping")
        return None

    have_reads = bool(pairs and "reads" in pairs)
    ncol = 2 if have_reads else 1
    fig, axarr = plt.subplots(1, ncol, figsize=(7.8 * ncol, 5.4),
                              squeeze=False)
    ax = axarr[0][0]
    items = sorted(ret.items(), key=lambda kv: -kv[1])
    ypos = list(range(len(items)))[::-1]
    vmax = max([v for _k, v in items] + [1.0])
    for y, (lab, v) in zip(ypos, items):
        ax.plot([0, v], [y, y], color=GRID, lw=2.4, zorder=1,
                solid_capstyle="round")
        ax.scatter([v], [y], s=135, color=blue_for(v, 100.0), zorder=3,
                   edgecolor=SURFACE, linewidth=1.4)
        txt = fmt_pct(v)
        if lab in counts and counts[lab][1]:
            txt += "  (%s / %s)" % (fmt_val(counts[lab][0]),
                                    fmt_val(counts[lab][1]))
        ax.annotate(txt, (v, y), textcoords="offset points", xytext=(12, 0),
                    va="center", fontsize=9.5, color=INK2)
    ax.set_yticks(ypos)
    ax.set_yticklabels([k for k, _v in items], fontsize=10, color=INK2)
    ax.set_ylim(-0.8, len(items) - 0.2)
    ax.set_xlim(0, max(105.0, vmax * 1.35))
    ax.set_xlabel("Percent retained at 5M reads", fontsize=10, color=INK2)
    ax.axvline(100.0, color=AXIS, lw=1.0, ls=":")
    ax.set_title("A. Retention of the full-depth signal at 5M reads",
                 loc="left", fontsize=11.5, fontweight="bold", color=INK,
                 pad=10)
    clean_axes(ax)

    if have_reads:
        per = pairs["reads"]
        it2 = sorted(per.items(), key=lambda kv: -max(kv[1]))
        dumbbell_panel(axarr[0][1], [k for k, _v in it2],
                       [v[1] for _k, v in it2], [v[0] for _k, v in it2],
                       SEQ[3], SEQ[1], "Reads (symlog scale)",
                       "B. Reads at full depth versus 5M subsample", log=True)
        h = [plt.Line2D([], [], marker="o", ls="", markersize=10, color=SEQ[3],
                        markeredgecolor=SURFACE, markeredgewidth=2,
                        label="Full depth"),
             plt.Line2D([], [], marker="o", ls="", markersize=9, color=SEQ[1],
                        markeredgecolor=SURFACE, markeredgewidth=2,
                        label="5M-read subsample")]
        fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.02, -0.055),
                   ncol=2, frameon=False, fontsize=10)

    best, worst = items[0], items[-1]
    if src == "samples":
        sentence = ("Subsampling to 5M reads keeps %s of the %s detections but "
                    "only %s of the %s detections"
                    % (fmt_pct(best[1]), best[0], fmt_pct(worst[1]), worst[0]))
    else:
        sentence = ("Subsampling to 5M reads retains %s of the %s signal and "
                    "%s of the %s signal"
                    % (fmt_pct(best[1]), best[0], fmt_pct(worst[1]), worst[0]))
    sub = ("Retention = 5M value / full-depth value, per category "
           "(source metric: %s)." % src)
    if note:
        sub = sub + "  " + note + "."
    sup(fig, sentence, sub)
    fig.tight_layout(rect=[0.015, 0.04, 1, 0.94])
    return save_fig(fig, args.outdir, "fig_depth_sensitivity_by_category",
                    args.dpi)


# ---------------------------------------------------------------- figure 4
def fig_profile(header, rows, args, anon):
    rcol = pick(header, ["reference_id", "reference", "refname", "ref",
                         "rname"])
    scol = sample_col(header)
    bcol = pick(header, ["bin_start", "start", "bin_index", "bin", "pos",
                         "position", "window_start"])
    ecol = pick(header, ["bin_end", "end", "window_end"], fuzzy=False)
    vcol = pick(header, ["mean_depth", "depth", "depth_sum", "n_reads",
                         "reads", "count", "coverage", "cov", "value"])
    lcol = pick(header, ["ref_label", "reference_label", "ref_name"],
                fuzzy=False)
    runcol = pick(header, ["run", "run_dir", "run_name"], fuzzy=False)
    if rcol is None or scol is None or bcol is None or vcol is None:
        print("WARN: refprofile_bins.tsv is missing reference / sample / bin / "
              "value columns, skipping")
        return None
    tracks = {}
    disp = {}
    runs = set()
    for r in rows:
        rid = str(r.get(rcol, "")).strip()
        smp = anon.get(r.get(scol))
        b = num(r.get(bcol))
        v = num(r.get(vcol))
        if not rid or b is None:
            continue
        run = str(r.get(runcol, "")).strip() if runcol else ""
        runs.add(run)
        if lcol:
            lab = str(r.get(lcol, "")).strip()
            if lab:
                disp.setdefault(rid, lab)
        e = num(r.get(ecol)) if ecol else None
        tracks.setdefault((rid, smp, run), []).append((b, e, v or 0.0))
    if not tracks:
        print("WARN: refprofile_bins.tsv has no usable bin rows, skipping")
        return None
    multirun = len(runs) > 1

    def rdisp(rid):
        lab = disp.get(rid)
        return lab if lab else ref_label(rid, 20)

    ref_tot = {}
    for (rid, _s, _run), pts in tracks.items():
        ref_tot[rid] = ref_tot.get(rid, 0.0) + sum(p[2] for p in pts)
    chosen = []
    for rid in sorted(ref_tot, key=lambda k: -ref_tot[k]):
        mine = [(k, v) for k, v in tracks.items() if k[0] == rid]
        mine.sort(key=lambda kv: (-sum(p[2] for p in kv[1]), kv[0][1]))
        for rank, (k, v) in enumerate(mine[:max(1, args.tracks_per_ref)]):
            chosen.append((k, v, rank))
        if len(chosen) >= args.max_tracks:
            break
    chosen = chosen[:max(1, args.max_tracks)]

    nrow = len(chosen)
    fig, axarr = plt.subplots(nrow, 1, figsize=(11.6, 1.28 * nrow + 1.5),
                              squeeze=False)
    unit = "bp" if any(p[1] is not None for _k, v, _r in chosen for p in v) \
        else "bin index"
    for i, ((rid, smp, run), pts, rank) in enumerate(chosen):
        ax = axarr[i][0]
        pts.sort(key=lambda p: p[0])
        xs = [p[0] for p in pts]
        ys = [p[2] for p in pts]
        last = pts[-1]
        step = (xs[-1] - xs[-2]) if len(xs) > 1 else 1.0
        xs_step = xs + [last[1] if last[1] is not None else xs[-1] + step]
        ys_step = ys + [ys[-1]]
        col = SEQ[min(rank + 1, len(SEQ) - 1)]
        ax.fill_between(xs_step, ys_step, step="post", color=col, linewidth=0)
        ax.plot(xs_step, ys_step, drawstyle="steps-post", color=SEQ[3],
                lw=0.8, alpha=0.75)
        nz = sum(1 for y in ys if y > 0)
        side = "%s\n%s" % (smp, rdisp(rid))
        if multirun and run:
            side = side + "\n" + run[:22]
        ax.set_ylabel(side, fontsize=9, color=INK2, rotation=0,
                      ha="right", va="center", labelpad=10)
        ax.annotate("max %s  |  %d / %d bins covered (%s)"
                    % (fmt_val(max(ys) if ys else 0.0), nz, len(ys),
                       fmt_pct(100.0 * nz / len(ys) if ys else 0.0)),
                    (1.0, 1.0), xycoords="axes fraction",
                    textcoords="offset points", xytext=(-2, -2),
                    ha="right", va="top", fontsize=8.6, color=INK2)
        top = max(ys) if ys else 1.0
        ax.set_ylim(0, (top * 1.35) if top > 0 else 1.0)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
        ax.tick_params(axis="y", labelsize=8, length=2)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.yaxis.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        # each panel keeps its own x tick labels: references differ in length
        ax.xaxis.set_major_formatter(FuncFormatter(thousands))
        ax.tick_params(axis="x", labelsize=8, length=2)
        if i == nrow - 1:
            ax.set_xlabel("Position on reference (%s)" % unit, fontsize=10,
                          color=INK2)
    sup(fig, "Per-bin coverage for the %d highest-signal sample tracks across "
             "%d reference(s)"
             % (nrow, len(set(k[0] for k, _v, _r in chosen))),
        "One panel per sample and reference; anonymous sample IDs; each panel "
        "is scaled independently.")
    fig.tight_layout(rect=[0.015, 0.01, 1, 0.955])
    return save_fig(fig, args.outdir, "fig_refprofile_coverage_tracks",
                    args.dpi)


# ---------------------------------------------------------------- figure 5
def fig_anello(header, rows, args, anon):
    scol = sample_col(header)
    gcol = pick(header, ["group", "cohort", "class", "arm"], fuzzy=False)
    # normalised burden first (comparable across samples), raw reads after
    bcol = pick(header, ["anello_rpm_human", "anello_rpm", "burden_rpm",
                         "reads_per_million", "rpm", "burden",
                         "anello_reads_human_total", "anello_reads_human",
                         "anello_reads", "total_reads", "reads", "n_reads",
                         "count"])
    rcol = pick(header, ["anello_richness_human", "anello_richness",
                         "richness", "n_references", "n_refs",
                         "distinct_references", "n_distinct", "n_species"])
    if scol is None or (bcol is None and rcol is None):
        print("WARN: anellovirus_burden.tsv is missing sample or burden / "
              "richness columns, skipping")
        return None
    data = []
    for r in rows:
        real = r.get(scol)
        grp = norm_group(r.get(gcol) if gcol else None, real)
        data.append((anon.get(real), grp,
                     num(r.get(bcol)) if bcol else None,
                     num(r.get(rcol)) if rcol else None))
    if not data:
        print("WARN: anellovirus_burden.tsv has no usable rows, skipping")
        return None
    data.sort(key=lambda d: d[0])
    groups = [g for g in GROUP_ORDER if any(d[1] == g for d in data)] or ["NA"]

    panels = []
    if bcol is not None:
        panels.append(("burden", "Anellovirus burden", 2))
    if rcol is not None:
        panels.append(("richness", "Distinct anellovirus references", 3))
    fig, axarr = plt.subplots(1, len(panels), figsize=(6.4 * len(panels), 5.4),
                              squeeze=False)
    rng = random.Random(20260726)
    jitter = dict((d[0], rng.uniform(-0.17, 0.17)) for d in data)
    med_txt = {}
    for pi, (kind, ylab, idx) in enumerate(panels):
        ax = axarr[0][pi]
        meds = {}
        for gi, g in enumerate(groups):
            vals = [d[idx] for d in data if d[1] == g and d[idx] is not None]
            xs = [gi + jitter[d[0]] for d in data
                  if d[1] == g and d[idx] is not None]
            if not vals:
                continue
            ax.scatter(xs, vals, s=46, color=GROUP_COLOR.get(g, MUTED),
                       alpha=0.85, edgecolor=SURFACE, linewidth=0.7, zorder=3)
            m = statistics.median(vals)
            meds[g] = m
            ax.plot([gi - 0.29, gi + 0.29], [m, m], color=INK, lw=2.0,
                    zorder=4, solid_capstyle="round")
            ax.annotate("median %s" % fmt_val(m), (gi + 0.31, m),
                        textcoords="offset points", xytext=(2, 0),
                        fontsize=9, color=INK, va="center")
            ax.annotate("n = %d" % len(vals), (gi, 0.0),
                        xycoords=("data", "axes fraction"),
                        textcoords="offset points", xytext=(0, -32),
                        ha="center", fontsize=9, color=MUTED)
        allv = [d[idx] for d in data if d[idx] is not None]
        pos = [v for v in allv if v > 0]
        if kind == "burden" and pos and (max(pos) / min(pos)) >= 100.0:
            if min(allv) <= 0:
                if try_symlog(ax, "y", 1.0):
                    ylab = ylab + " (symlog scale)"
            else:
                try:
                    ax.set_yscale("log")
                    ylab = ylab + " (log scale)"
                except Exception:
                    pass
            ax.yaxis.set_major_formatter(FuncFormatter(thousands))
        else:
            top = max(allv) if allv else 1.0
            ax.set_ylim(0, (top * 1.25) if top > 0 else 1.0)
        ax.set_xticks(list(range(len(groups))))
        ax.set_xticklabels(groups, fontsize=11, color=INK2)
        ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.set_ylabel(ylab, fontsize=10, color=INK2)
        ax.set_title("%s. %s by group"
                     % ("AB"[pi], "Burden" if kind == "burden" else "Richness"),
                     loc="left", fontsize=11.5, fontweight="bold", color=INK,
                     pad=10)
        clean_axes(ax, xgrid=False, ygrid=True, hide_left=False)
        ax.spines["bottom"].set_color(AXIS)
        med_txt[kind] = meds

    bm = med_txt.get("burden") or med_txt.get("richness") or {}
    if len(bm) >= 2:
        top_g = max(bm, key=lambda k: bm[k])
        low_g = min(bm, key=lambda k: bm[k])
        sentence = ("Anellovirus burden is highest in the %s group "
                    "(median %s versus %s in %s)"
                    % (top_g, fmt_val(bm[top_g]), fmt_val(bm[low_g]), low_g))
    elif bm:
        only = list(bm)[0]
        sentence = ("Anellovirus burden in the %s group has a median of %s"
                    % (only, fmt_val(bm[only])))
    else:
        sentence = "Anellovirus burden and richness by group"
    sub = ("One dot per sample (anonymous IDs), horizontal jitter only; "
           "black bar = group median.")
    cols = []
    if bcol is not None:
        cols.append("burden = " + str(bcol))
    if rcol is not None:
        cols.append("richness = " + str(rcol))
    if cols:
        sub = sub + "  Source columns: " + ", ".join(cols) + "."
    sup(fig, sentence, sub)
    fig.tight_layout(rect=[0.015, 0.05, 1, 0.94])
    return save_fig(fig, args.outdir, "fig_anellovirus_burden", args.dpi)


# ---------------------------------------------------------------- figure 6
def fig_coinf(header, rows, args):
    # exact matches only: single-letter aliases would otherwise capture
    # unrelated columns such as "n_both" for side "b"
    acol = pick(header, ["virus_group_a", "ref_a", "reference_a", "virus_a",
                         "group_a", "ref1", "ref_1", "reference_1", "virus_1",
                         "a", "pair_a", "left"], fuzzy=False)
    bcol = pick(header, ["virus_group_b", "ref_b", "reference_b", "virus_b",
                         "group_b", "ref2", "ref_2", "reference_2", "virus_2",
                         "b", "pair_b", "right"], fuzzy=False)
    vcol = pick(header, ["n_both", "samples_both", "both", "n_co",
                         "co_occurrence", "cooccurrence", "jaccard", "count",
                         "n_shared", "n"], fuzzy=False)
    if len(set([acol, bcol, vcol]) - set([None])) < 3:
        print("WARN: coinfection_pairs.tsv is missing distinct pair or count "
              "columns, skipping")
        return None
    # a stratified table: use the pooled rows if they exist, else one stratum
    cohcol = pick(header, ["cohort", "group", "run", "arm"], fuzzy=False)
    note = None
    if cohcol:
        seen = sorted(set(str(r.get(cohcol, "")).strip() for r in rows))
        seen = [v for v in seen if v]
        pooled = [v for v in seen if v.upper() in ALL_VALUES]
        if len(seen) > 1:
            keep = pooled[0] if pooled else seen[0]
            rows = [r for r in rows
                    if str(r.get(cohcol, "")).strip() == keep]
            note = "%s = %s rows only" % (cohcol, keep)
    ncol_a = pick(header, ["n_a", "n_samples_a", "count_a"], fuzzy=False)
    ncol_b = pick(header, ["n_b", "n_samples_b", "count_b"], fuzzy=False)
    vals = {}
    marg = {}
    for r in rows:
        a = str(r.get(acol, "")).strip()
        b = str(r.get(bcol, "")).strip()
        v = num(r.get(vcol))
        if not a or not b or v is None:
            continue
        vals[(a, b)] = v
        vals[(b, a)] = v
        if ncol_a:                       # single-reference counts, if given
            na = num(r.get(ncol_a))
            if na is not None:
                vals[(a, a)] = na
        if ncol_b:
            nb = num(r.get(ncol_b))
            if nb is not None:
                vals[(b, b)] = nb
        if a != b:
            marg[a] = marg.get(a, 0.0) + v
            marg[b] = marg.get(b, 0.0) + v
        else:
            marg.setdefault(a, 0.0)
    if not vals:
        print("WARN: coinfection_pairs.tsv has no usable pair rows, skipping")
        return None
    labs = sorted(marg, key=lambda k: (-marg[k], k))[:max(2, args.max_refs)]
    n = len(labs)
    mat = []
    for i in range(n):
        mat.append([vals.get((labs[i], labs[j])) if j <= i else None
                    for j in range(n)])
    flat = [v for row in mat for v in row if v is not None]
    vmax = max(flat) if flat else 1.0
    integral = all(abs(v - round(v)) < 1e-9 for v in flat)

    size = max(6.0, 0.72 * n + 3.4)
    fig, ax = plt.subplots(figsize=(size + 1.8, size))
    cmap = blue_cmap()
    try:
        import numpy                                   # ships with matplotlib
        arr = numpy.ma.masked_invalid(numpy.array(
            [[numpy.nan if v is None else float(v) for v in row]
             for row in mat], dtype=float))
        cmap.set_bad(SURFACE)
    except Exception:
        arr = [[(0.0 if v is None else v) for v in row] for row in mat]
    im = ax.imshow(arr, cmap=cmap, vmin=0, vmax=vmax, aspect="equal")
    for i in range(n):
        for j in range(n):
            v = mat[i][j]
            if v is None:
                continue
            frac = (v / vmax) if vmax else 0.0
            ax.text(j, i, fmt_int(v) if integral else "%.2f" % v,
                    ha="center", va="center", fontsize=9,
                    color="#ffffff" if frac > 0.55 else INK,
                    fontweight="bold" if frac > 0.55 else "normal")
    names = [ref_label(x, 22) for x in labs]
    ax.set_xticks(list(range(n)))
    ax.set_yticks(list(range(n)))
    ax.set_xticklabels(names, fontsize=9, color=INK2, rotation=45, ha="right")
    ax.set_yticklabels(names, fontsize=9, color=INK2)
    ax.set_xticks([x - 0.5 for x in range(1, n)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n)], minor=True)
    ax.grid(which="minor", color=SURFACE, lw=1.6)
    ax.tick_params(which="both", length=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
    cb.outline.set_visible(False)
    cb.set_label("Samples carrying both references", fontsize=9.5, color=INK2)
    cb.ax.tick_params(labelsize=8.5)

    off = [(k, v) for k, v in vals.items()
           if k[0] != k[1] and k[0] in labs and k[1] in labs]
    if off:
        (pa, pb), pv = max(off, key=lambda kv: kv[1])
        sentence = ("Co-infection is dominated by the %s + %s pair (%s samples "
                    "carry both)" % (ref_label(pa, 22), ref_label(pb, 22),
                                     fmt_val(pv)))
    else:
        sentence = "Co-occurrence of viral references across samples"
    sub = ("Lower triangle only; the diagonal is the single-reference count. "
           "Every cell is labelled, so colour is redundant.")
    if note:
        sub = sub + "  " + note + "."
    sup(fig, sentence, sub)
    fig.tight_layout(rect=[0.015, 0.01, 1, 0.94])
    return save_fig(fig, args.outdir, "fig_coinfection_pairs", args.dpi)


# ===================================================================== main
def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Render the suite figures from the .tsv outputs of the "
                    "other modules (no BAM access, no network).")
    p.add_argument("--indir", default=DEF_INDIR,
                   help="directory holding the suite .tsv outputs "
                        "(default: %(default)s; falls back to the current "
                        "directory if the default does not exist)")
    p.add_argument("--outdir", default=DEF_OUTDIR,
                   help="directory for fig_*.png / fig_*.svg and the index "
                        "tsv (default: %(default)s)")
    p.add_argument("--prefix", default=DEF_PREFIX,
                   help="prefix for this module's .tsv outputs "
                        "(default: %(default)s)")
    p.add_argument("--sweep-tsv", default=None,
                   help="override path to detection_threshold_sweep.tsv")
    p.add_argument("--refcmp-tsv", default=None,
                   help="override path to "
                        "reference_comparison_by_category.tsv")
    p.add_argument("--depth-tsv", default=None,
                   help="override path to depth_sensitivity_by_category.tsv")
    p.add_argument("--profile-tsv", default=None,
                   help="override path to refprofile_bins.tsv")
    p.add_argument("--anello-tsv", default=None,
                   help="override path to anellovirus_burden.tsv")
    p.add_argument("--coinf-tsv", default=None,
                   help="override path to coinfection_pairs.tsv")
    p.add_argument("--default-threshold", type=float, default=DEF_THRESHOLD,
                   help="detection threshold to mark on the sweep figure "
                        "(default: %(default)s reads)")
    p.add_argument("--max-tracks", type=int, default=9,
                   help="max coverage tracks in the profile figure "
                        "(default: %(default)s)")
    p.add_argument("--tracks-per-ref", type=int, default=3,
                   help="max sample tracks per reference "
                        "(default: %(default)s)")
    p.add_argument("--max-refs", type=int, default=12,
                   help="max references in the co-infection heatmap "
                        "(default: %(default)s)")
    p.add_argument("--dpi", type=int, default=200,
                   help="png resolution (default: %(default)s)")
    p.add_argument("--only", default=None,
                   help="comma-separated subset of figures to draw: "
                        "sweep,refcmp,depth,profile,anello,coinf")
    return p.parse_args(argv)


def resolve_input_path(indir, fname, alts, override):
    """Exact filename first, then any "*<fname>" in indir, then the alternate
    suffixes, so prefixed upstream outputs such as
    a7_virome_anellovirus_burden.tsv or a1_threshold_sweep.tsv are found.
    A file that matches another input's canonical name is never borrowed."""
    if override:
        return override
    exact = os.path.join(indir, fname)
    if os.path.isfile(exact):
        return exact
    try:
        here = sorted(os.listdir(indir))
    except OSError:
        here = []
    others = [n for n in PRIMARY_NAMES if n != fname]
    for suffix in [fname] + list(alts):
        for g in here:
            if g == fname or not g.endswith(suffix):
                continue
            if any(g == o or g.endswith(o) for o in others):
                continue
            print("NOTE: using %s in place of %s (prefixed upstream output)"
                  % (g, fname))
            return os.path.join(indir, g)
    return exact


def resolve_indir(args, argv):
    """If --indir was not given and the cluster default is absent, fall back to
    the current directory so the module also works on a laptop."""
    given = any(a == "--indir" or a.startswith("--indir=") for a in argv)
    if given or os.path.isdir(args.indir):
        return args.indir
    cwd = os.getcwd()
    try:
        here = os.listdir(cwd)
    except OSError:
        here = []
    hits = [fn for _k, fn, _l, _alt in INPUTS
            if any(g == fn or g.endswith(fn) for g in here)]
    if hits:
        print("NOTE: default --indir %s not present, using the current "
              "directory (%d input file(s) found there)"
              % (args.indir, len(hits)))
        return cwd
    return args.indir


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parse_args(argv)
    if not HAVE_MPL:
        print("WARN: matplotlib unavailable (%s), no figures drawn, skipping"
              % MPL_ERR)
        return 0
    style_rcparams()

    args.indir = resolve_indir(args, argv)
    if not os.path.isdir(args.indir):
        print("WARN: input directory missing at %s, skipping" % args.indir)
    try:
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir)
    except OSError as exc:
        print("WARN: outdir %s could not be created (%s), skipping"
              % (args.outdir, exc))
        return 0

    overrides = {"sweep": args.sweep_tsv, "refcmp": args.refcmp_tsv,
                 "depth": args.depth_tsv, "profile": args.profile_tsv,
                 "anello": args.anello_tsv, "coinf": args.coinf_tsv}
    only = None
    if args.only:
        only = set(w.strip().lower() for w in args.only.split(",") if w.strip())

    paths = {}
    tables = {}
    for key, fname, label, alts in INPUTS:
        path = resolve_input_path(args.indir, fname, alts, overrides.get(key))
        paths[key] = path
        if only is not None and key not in only:
            tables[key] = None
            continue
        header, rows = read_tsv(path, label)
        tables[key] = (header, rows) if header else None

    anon = Anon()
    prescan_samples(tables, anon)
    key_path = None
    if anon.n_mapped:
        try:
            key_path = anon.write_key(os.path.join(
                args.outdir, "%s_sample_key.tsv" % args.prefix))
        except OSError as exc:
            print("WARN: sample key not written (%s), skipping" % exc)

    plans = [
        ("detection_threshold_sweep", "sweep",
         lambda h, r: fig_sweep(h, r, args)),
        ("reference_comparison_by_category", "refcmp",
         lambda h, r: fig_refcmp(h, r, args)),
        ("depth_sensitivity_by_category", "depth",
         lambda h, r: fig_depth(h, r, args)),
        ("refprofile_coverage_tracks", "profile",
         lambda h, r: fig_profile(h, r, args, anon)),
        ("anellovirus_burden", "anello",
         lambda h, r: fig_anello(h, r, args, anon)),
        ("coinfection_pairs", "coinf",
         lambda h, r: fig_coinf(h, r, args)),
    ]

    index = []
    made = 0
    for name, key, fn in plans:
        got = tables.get(key)
        nrows = len(got[1]) if got else 0
        if only is not None and key not in only:
            index.append((name, paths[key], nrows, "not_selected", "", ""))
            continue
        if not got:
            index.append((name, paths[key], nrows, "skipped_missing_input",
                          "", ""))
            continue
        try:
            out = fn(got[0], got[1])
        except Exception as exc:
            print("WARN: figure %s failed (%s: %s), skipping"
                  % (name, type(exc).__name__, exc))
            out = None
        if out:
            made += 1
            index.append((name, paths[key], nrows, "written",
                          os.path.basename(out[0]), os.path.basename(out[1])))
        else:
            index.append((name, paths[key], nrows, "skipped_unusable_input",
                          "", ""))

    idx_path = os.path.join(args.outdir, "%s_figure_index.tsv" % args.prefix)
    lines = ["# a8_figures.py figure index, written %s" % TODAY,
             "# anonymous sample IDs only; no real sample names in this file",
             "\t".join(["figure", "source_tsv", "source_rows", "status",
                        "png", "svg"])]
    for row in index:
        lines.append("\t".join([str(x) for x in row]))
    try:
        with open(idx_path, "w", encoding="ascii", newline="") as fh:
            fh.write("\n".join(lines) + "\n")
    except (OSError, UnicodeEncodeError) as exc:
        print("WARN: figure index not written to %s (%s), skipping"
              % (idx_path, exc))
        idx_path = None

    print("")
    print("a8_figures.py  %s" % TODAY)
    print("  indir : %s" % args.indir)
    print("  outdir: %s" % args.outdir)
    print("  figures written: %d of %d" % (made, len(plans)))
    for name, src, nrows, status, png, _svg in index:
        print("    %-34s %-40s rows=%-6d %s"
              % (name, os.path.basename(src), nrows,
                 status if not png else png + " (+ .svg)"))
    if key_path:
        print("  sample key: %s" % key_path)
        print("              %d real name(s) mapped to S-IDs -- CONTAINS "
              "IDENTIFIERS, do not commit or email" % anon.n_mapped)
    else:
        print("  sample key: not written (no real sample names in the inputs)")
    if idx_path:
        print("  index     : %s" % idx_path)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("WARN: interrupted, stopping")
        sys.exit(0)
