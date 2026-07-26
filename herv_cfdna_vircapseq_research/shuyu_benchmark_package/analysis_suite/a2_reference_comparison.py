#!/usr/bin/env python3
"""a2_reference_comparison.py -- what does ADDING THE HUMAN GENOME as a competitor do?

WHAT IT COMPUTES
  Two families of paired comparisons, both driven by the headline deduplicated
  category tables (results/[<prefix>_]filtered_category_counts.tsv):

  1. "viral_only_vs_with_human"  (the main question)
     For each cohort, the VIRAL-ONLY reference run (*_full_competitive, arm A =
     permissive) is compared against the hg38-inclusive run of the same cohort
     (*_hg38_refseq_mapq_human60_viral40_coord, arm B = stringent), on the
     intersection of samples present in both runs. Per category we report reads
     and nonzero-sample counts in each arm, the absolute and percentage change,
     and the number of samples that FLIP positive -> negative, i.e. calls that
     the human competitor removes (candidate false positives).

  2. "all_vs_primary_only"  (multi-mapping inflation)
     For any run that ships BOTH filename prefixes of the headline table
     (filtered_category_counts.tsv and primary_only_filtered_category_counts.tsv)
     the same comparison is run inside that single run, arm A = all alignments,
     arm B = primary-only. The delta is the secondary/supplementary-alignment
     inflation.

  In both families arm A is the more permissive setting and arm B the stricter
  one, so delta_reads = reads_b - reads_a is normally <= 0 and a pos -> neg flip
  is always "a call the stricter setting removes". Two aggregate rows are added
  per comparison: TOTAL_VIRAL (categories present in BOTH references and not in
  --nonviral-categories) and TOTAL_ALL (every category in either reference, so
  it is dominated by HUMAN in the hg38-inclusive arm). Categories that exist in
  only one arm's reference are flagged structural=yes and their zero ->
  positive transitions are excluded from the flip counters and the flips file,
  because they reflect the reference content rather than a changed call.

  Optional (--verify-flips N, off by default): for up to N pos -> neg flips of
  family 1, re-count unique-best reads (AS > XS, MAPQ >= --mapq, -F 0x904)
  directly on the arm-B BAM, restricted to the reference_ids that carried the
  signal in arm A. This calls samtools via subprocess; nothing else here does.

WHAT IT WRITES (tab-separated, into --outdir, anon sample IDs only)
  <prefix>_reference_comparison_by_category.tsv   per comparison x category
  <prefix>_reference_comparison_by_sample.tsv     per comparison x sample x category
  <prefix>_reference_comparison_flips.tsv         one row per flip (both directions)
  <prefix>_sample_key.tsv                         real -> anon map; CONTAINS IDENTIFIERS

  Every sample is anonymised to S01..Snn (sorted by the real sample name) in all
  outputs except the sample key. No figures are produced, so matplotlib is not
  imported; standard library only.

EXAMPLE
  python3 a2_reference_comparison.py \
      --runs-root /path/to/runs \
      --outdir /path/to/runs/reply_to_shuyu_primary_only/a2_out \
      --prefix a2 --verify-flips 12

Date: 2026-07-26
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import subprocess
import sys

# ----------------------------------------------------------------------------
# defaults
# ----------------------------------------------------------------------------
DEF_RUNS_ROOT = "/path/to/runs"
DEF_OUTDIR = "./a2_reference_comparison_out"
DEF_PREFIX = "a2"

DEF_PANEL_REFMAP = ("/path/to/runs/shuyu_masked_panel_hg38_herv_line1_refixed/"
                    "ref/hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv")
DEF_BASE_REFMAP = ("/path/to/runs/retro_reference_hg38_refseq/"
                   "ref/hg38_plus_retro.refseq.reference_map.csv")

# cohort, viral-only run (arm A, permissive), hg38-inclusive run (arm B, stringent)
DEF_PAIRS = [
    ("targeted", "targeted_htlv_full_competitive",
     "targeted_htlv_hg38_refseq_mapq_human60_viral40_coord"),
    ("wgs", "wgs_hiv_hl_full_competitive",
     "wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord"),
]

HEADLINE_SUFFIX = "filtered_category_counts.tsv"
VARIANT_ALL = "all_alignments"
VARIANT_PRIMARY = "primary_only"
VARIANT_PREFERENCE = [VARIANT_PRIMARY, VARIANT_ALL]

DEF_NONVIRAL = "HUMAN,HERV,LINE1"

CMP_HUMAN = "viral_only_vs_with_human"
CMP_PRIMARY = "all_vs_primary_only"

KEY_WARNING = "# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL"


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
def warn(what, path):
    print("WARN: %s missing at %s, skipping" % (what, path))


def ascii_safe(text):
    """Force any string into pure ASCII before it reaches an output file."""
    if text is None:
        return ""
    return str(text).encode("ascii", "replace").decode("ascii")


def parse_count(raw):
    raw = (raw or "").strip()
    if not raw or raw in (".", "NA", "NaN", "nan"):
        return 0
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return int(round(float(raw)))
    except ValueError:
        return 0


def pct_change(a, b):
    if a == 0:
        return "NA"
    return "%.2f" % (100.0 * (b - a) / float(a))


def group_from_name(sample, context=""):
    """Group label from the real sample name; run name used only as a fallback.

    Suite-wide rule, matched case-insensitively so it behaves the same in every
    module: "_HIV" -> HIV, "_HL" -> HL, "TCL"/"targeted_htlv" -> TCL, else NA.
    """
    up = (sample or "").upper()
    if "_HIV" in up:
        return "HIV"
    if "_HL" in up:
        return "HL"
    if "TARGETED_HTLV" in up or "TCL" in up:
        return "TCL"
    ctx = (context or "").upper()
    if "TARGETED_HTLV" in ctx or "TCL" in ctx:
        return "TCL"
    return "NA"


def norm_sample(raw):
    """Sample names are assumed identical across runs; only trivial noise is stripped."""
    s = (raw or "").strip()
    if s.endswith(".bam"):
        s = s[:-4]
    return s


def cohort_of_run(run_name):
    low = os.path.basename(run_name).lower()
    if "targeted" in low:
        return "targeted"
    if "wgs" in low:
        return "wgs"
    return "other"


# ----------------------------------------------------------------------------
# input discovery / loading
# ----------------------------------------------------------------------------
def resolve_run(runs_root, name):
    if os.path.isabs(name):
        return name
    return os.path.join(runs_root, name)


def headline_variants(run_dir):
    """Map variant label -> path for every *filtered_category_counts.tsv in a run."""
    out = {}
    res = os.path.join(run_dir, "results")
    if not os.path.isdir(res):
        return out
    for path in sorted(glob.glob(os.path.join(res, "*" + HEADLINE_SUFFIX))):
        base = os.path.basename(path)
        stem = base[:-len(HEADLINE_SUFFIX)].strip("_")
        if not stem:
            out[VARIANT_ALL] = path
        elif stem == "primary_only":
            out[VARIANT_PRIMARY] = path
        else:
            out[stem] = path
    return out


def load_category_table(path):
    """-> (categories, {sample: {category: reads}}) or (None, None) on failure."""
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except (IOError, OSError):
        return None, None
    rows = [r for r in rows if r and any((c or "").strip() for c in r)]
    if len(rows) < 1:
        return None, None
    header = [(c or "").strip() for c in rows[0]]
    cats = [c for c in header[1:] if c]
    table = {}
    for row in rows[1:]:
        if not row:
            continue
        sample = norm_sample(row[0])
        if not sample or sample.startswith("#"):
            continue
        vals = {}
        for i, cat in enumerate(header[1:], start=1):
            if not cat:
                continue
            vals[cat] = parse_count(row[i] if i < len(row) else "")
        prev = table.get(sample)
        if prev is None:
            table[sample] = vals
        else:                                    # duplicate sample row: sum it
            for cat, v in vals.items():
                prev[cat] = prev.get(cat, 0) + v
    return cats, table


def refmap_from_manifest(run_dir):
    man = os.path.join(run_dir, "results", "run_manifest.tsv")
    if not os.path.exists(man):
        return None
    try:
        with open(man, newline="", encoding="utf-8", errors="replace") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except (IOError, OSError):
        return None
    if not rows:
        return None
    for key in rows[0]:
        if key and "reference_map" in key.lower():
            val = (rows[0][key] or "").strip()
            if val and os.path.exists(val):
                return val
    for val in rows[0].values():
        val = (val or "").strip()
        if val.endswith(".csv") and "reference_map" in val and os.path.exists(val):
            return val
    return None


def load_refmap(path):
    """-> {reference_id: category}"""
    out = {}
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh):
                rid = (row.get("reference_id") or "").strip()
                cat = (row.get("category") or "").strip()
                if rid:
                    out[rid] = cat
    except (IOError, OSError):
        return {}
    return out


def refmap_for_run(run_dir, args, cache):
    if run_dir in cache:
        return cache[run_dir]
    path = refmap_from_manifest(run_dir)
    if not path:
        guess = args.panel_refmap if "panel" in os.path.basename(run_dir).lower() \
            else args.base_refmap
        path = guess if guess and os.path.exists(guess) else None
    info = load_refmap(path) if path else {}
    if not info:
        warn("reference map for run %s" % os.path.basename(run_dir),
             path or "<not resolved>")
    cache[run_dir] = (path, info)
    return cache[run_dir]


def refs_with_signal(run_dir, sample, category, refmap):
    """reference_ids of one category that carried mapped reads for this sample."""
    idx = os.path.join(run_dir, "results", sample + ".idxstats.tsv")
    if not os.path.exists(idx):
        return None
    hits = []
    try:
        with open(idx, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3 or parts[0] == "*":
                    continue
                if parse_count(parts[2]) <= 0:
                    continue
                if refmap.get(parts[0], "") == category:
                    hits.append(parts[0])
    except (IOError, OSError):
        return None
    return hits


def count_unique_best(samtools, bam, refs, mapq):
    """Unique-best (AS>XS or no XS) primary reads at MAPQ>=mapq on given refs."""
    cmd = [samtools, "view", "-F", "0x904", "-q", str(mapq), bam] + list(refs)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL,
                                universal_newlines=True)
    except (OSError, ValueError):
        return None
    total = 0
    try:
        for line in proc.stdout:
            if not line or line[0] == "@":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            a_score = None
            x_score = None
            for tag in fields[11:]:
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
            if a_score is not None and (x_score is None or a_score > x_score):
                total += 1
    finally:
        try:
            proc.stdout.close()
        except (IOError, OSError):
            pass
        rc = proc.wait()
    if rc != 0:
        return None
    return total


# ----------------------------------------------------------------------------
# comparison engine
# ----------------------------------------------------------------------------
class Arm(object):
    def __init__(self, label, run_dir, variant, path, cats, table):
        self.label = label
        self.run_dir = run_dir
        self.run = os.path.basename(run_dir.rstrip("/\\"))
        self.variant = variant
        self.path = path
        self.cats = cats
        self.table = table


def load_arm(label, run_dir, variant, variants):
    path = variants.get(variant)
    if not path:
        return None
    cats, table = load_category_table(path)
    if table is None:
        warn("headline category table for %s" % os.path.basename(run_dir), path)
        return None
    if not table:
        warn("sample rows in headline category table for %s"
             % os.path.basename(run_dir), path)
        return None
    return Arm(label, run_dir, variant, path, cats, table)


def pick_common_variant(vars_a, vars_b):
    for v in VARIANT_PREFERENCE:
        if v in vars_a and v in vars_b:
            return v, v
    va = sorted(vars_a)[0] if vars_a else None
    vb = sorted(vars_b)[0] if vars_b else None
    return va, vb


def compare(cmp_name, cohort, arm_a, arm_b, nonviral, keep_zero):
    """Build (category_rows, sample_rows, flip_rows, ...) for one paired comparison.

    A category that exists in only one arm's reference (HUMAN/HERV/LINE1 in the
    viral-only vs with-human family) is marked structural=yes; its per-sample
    zero -> positive transitions are a property of the reference, not a call
    change, so they are excluded from the flip counters and from the flips file.
    """
    samples_a = set(arm_a.table)
    samples_b = set(arm_b.table)
    shared = sorted(samples_a & samples_b)
    only_a = len(samples_a - samples_b)
    only_b = len(samples_b - samples_a)

    cats = list(arm_a.cats)
    for c in arm_b.cats:
        if c not in cats:
            cats.append(c)
    both = [c for c in cats if c in arm_a.cats and c in arm_b.cats]
    viral_cats = [c for c in both if c.upper() not in nonviral]

    cat_rows, sample_rows, flip_rows = [], [], []

    def agg(label, use_cats, in_a, in_b, structural):
        reads_a = reads_b = nz_a = nz_b = 0
        f_pn = f_np = 0
        for s in shared:
            va = sum(arm_a.table[s].get(c, 0) for c in use_cats)
            vb = sum(arm_b.table[s].get(c, 0) for c in use_cats)
            reads_a += va
            reads_b += vb
            nz_a += 1 if va > 0 else 0
            nz_b += 1 if vb > 0 else 0
            if structural:
                continue
            if va > 0 and vb == 0:
                f_pn += 1
            elif va == 0 and vb > 0:
                f_np += 1
        cat_rows.append([cmp_name, cohort, arm_a.label, arm_b.label,
                         arm_a.run, arm_b.run, arm_a.variant, arm_b.variant,
                         label, in_a, in_b, "yes" if structural else "no",
                         len(shared),
                         reads_a, reads_b, reads_b - reads_a,
                         pct_change(reads_a, reads_b),
                         nz_a, nz_b, nz_b - nz_a,
                         "NA" if structural else f_pn,
                         "NA" if structural else f_np,
                         only_a, only_b])

    for cat in cats:
        in_a = "yes" if cat in arm_a.cats else "no"
        in_b = "yes" if cat in arm_b.cats else "no"
        structural = (in_a == "no" or in_b == "no")
        agg(cat, [cat], in_a, in_b, structural)
        for s in shared:
            va = arm_a.table[s].get(cat, 0)
            vb = arm_b.table[s].get(cat, 0)
            if va == 0 and vb == 0 and not keep_zero:
                continue
            if structural:
                flip = "structural_only_in_b" if in_a == "no" else "structural_only_in_a"
            elif va > 0 and vb == 0:
                flip = "pos_to_neg"
            elif va == 0 and vb > 0:
                flip = "neg_to_pos"
            else:
                flip = "none"
            sample_rows.append([cmp_name, cohort, s, arm_a.run, arm_b.run,
                                arm_a.variant, arm_b.variant, cat,
                                va, vb, vb - va, pct_change(va, vb), flip])
            if flip in ("pos_to_neg", "neg_to_pos"):
                flip_rows.append([cmp_name, cohort, s, cat, flip, va, vb,
                                  arm_a.run, arm_b.run,
                                  arm_a.variant, arm_b.variant])
    if viral_cats:
        agg("TOTAL_VIRAL", viral_cats, "yes", "yes", False)
    agg("TOTAL_ALL", cats, "yes", "yes", False)
    return cat_rows, sample_rows, flip_rows, shared, only_a, only_b


# ----------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------
CAT_HEADER = ["comparison", "cohort", "arm_a", "arm_b", "run_a", "run_b",
              "variant_a", "variant_b", "category", "category_in_a",
              "category_in_b", "structural", "n_samples_compared",
              "reads_a", "reads_b", "delta_reads", "pct_change_reads",
              "nonzero_samples_a", "nonzero_samples_b",
              "delta_nonzero_samples", "flips_pos_to_neg", "flips_neg_to_pos",
              "samples_only_in_a", "samples_only_in_b"]

SAMPLE_HEADER = ["comparison", "cohort", "sample", "group", "run_a", "run_b",
                 "variant_a", "variant_b", "category", "reads_a", "reads_b",
                 "delta_reads", "pct_change_reads", "flip"]

FLIP_HEADER = ["comparison", "cohort", "sample", "group", "category",
               "direction", "reads_a", "reads_b", "run_a", "run_b",
               "variant_a", "variant_b", "verify_refs_checked",
               "verify_unique_best_reads_in_b", "verify_mapq"]


def write_tsv(path, header, rows):
    with open(path, "w", encoding="ascii", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([ascii_safe(c) for c in header])
        for row in rows:
            w.writerow([ascii_safe(c) for c in row])


def write_sample_key(path, anon, groups, seen_runs):
    with open(path, "w", encoding="ascii", newline="") as fh:
        fh.write(KEY_WARNING + "\n")
        fh.write("# generated 2026-07-26 by a2_reference_comparison.py\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["anon_id", "real_sample", "group", "runs"])
        for real in sorted(anon, key=lambda s: anon[s]):
            w.writerow([anon[real], ascii_safe(real), groups[real],
                        ascii_safe(";".join(sorted(seen_runs.get(real, []))))])


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def build_args():
    ap = argparse.ArgumentParser(
        description="Quantify the effect of adding hg38 as a competitor "
                    "(viral-only vs hg38-inclusive runs), plus all-alignment "
                    "vs primary-only multi-mapping inflation.")
    ap.add_argument("--runs-root", default=DEF_RUNS_ROOT,
                    help="root holding the run directories (default: %s)" % DEF_RUNS_ROOT)
    ap.add_argument("--outdir", default=DEF_OUTDIR,
                    help="output directory (default: %s)" % DEF_OUTDIR)
    ap.add_argument("--prefix", default=DEF_PREFIX,
                    help="output filename prefix (default: %s)" % DEF_PREFIX)
    ap.add_argument("--pair", action="append", default=None, metavar="COHORT,RUN_A,RUN_B",
                    help="comma-separated cohort,viral_only_run,with_human_run; "
                         "repeatable; run names may be basenames under --runs-root "
                         "or absolute paths. Overrides the built-in pairs.")
    ap.add_argument("--primary-runs", action="append", default=None, metavar="RUN",
                    help="run to use for the all-vs-primary_only comparison; "
                         "repeatable. Default: every run under --runs-root that "
                         "ships both headline filename prefixes.")
    ap.add_argument("--no-primary-comparison", action="store_true",
                    help="skip the all-alignment vs primary-only family")
    ap.add_argument("--base-refmap", default=DEF_BASE_REFMAP,
                    help="fallback reference map for non-panel runs")
    ap.add_argument("--panel-refmap", default=DEF_PANEL_REFMAP,
                    help="fallback reference map for panel runs")
    ap.add_argument("--nonviral-categories", default=DEF_NONVIRAL,
                    help="comma-separated categories excluded from TOTAL_VIRAL "
                         "(default: %s)" % DEF_NONVIRAL)
    ap.add_argument("--keep-zero-rows", action="store_true",
                    help="keep per-sample rows where both arms are zero")
    ap.add_argument("--verify-flips", type=int, default=0, metavar="N",
                    help="re-count unique-best reads on the arm-B BAM for up to N "
                         "pos->neg flips of the viral-only vs with-human family "
                         "(0 = off, the default; needs samtools + .bai)")
    ap.add_argument("--verify-max-refs", type=int, default=25,
                    help="cap on reference_ids passed to samtools per flip (default 25)")
    ap.add_argument("--mapq", type=int, default=40,
                    help="MAPQ floor for the optional samtools verification (default 40)")
    ap.add_argument("--samtools", default="samtools",
                    help="samtools executable (default: samtools)")
    return ap.parse_args()


def parse_pairs(args):
    if not args.pair:
        return list(DEF_PAIRS)
    pairs = []
    for spec in args.pair:
        parts = [p.strip() for p in spec.split(",")]
        if len(parts) != 3 or not all(parts):
            print("WARN: cannot parse --pair '%s', expected COHORT,RUN_A,RUN_B, skipping"
                  % ascii_safe(spec))
            continue
        pairs.append((parts[0], parts[1], parts[2]))
    return pairs


def discover_primary_runs(runs_root):
    out = []
    if not os.path.isdir(runs_root):
        return out
    for d in sorted(glob.glob(os.path.join(runs_root, "*"))):
        if not os.path.isdir(d):
            continue
        v = headline_variants(d)
        if VARIANT_ALL in v and VARIANT_PRIMARY in v:
            out.append(d)
    return out


def print_block(cmp_name, cohort, arm_a, arm_b, cat_rows, n_shared, only_a, only_b):
    print("")
    print("=" * 78)
    print("%s | cohort=%s" % (cmp_name, cohort))
    print("  arm A (%s): %s [%s]" % (arm_a.label, arm_a.run, arm_a.variant))
    print("  arm B (%s): %s [%s]" % (arm_b.label, arm_b.run, arm_b.variant))
    print("  samples compared: %d (only in A: %d, only in B: %d)"
          % (n_shared, only_a, only_b))
    print("  %-14s %14s %14s %14s %9s %13s %9s %9s"
          % ("category", "reads_A", "reads_B", "delta", "pct_chg",
             "nonzero A->B", "flip_p2n", "flip_n2p"))
    ordered = sorted(cat_rows, key=lambda r: -abs(r[15]))
    for r in ordered:
        pct = r[16] if r[16] == "NA" else r[16] + "%"
        struct = "  (reference-only category)" if r[11] == "yes" else ""
        print("  %-14s %14d %14d %14d %9s %6d -> %-4d %9s %9s%s"
              % (r[8][:14], r[13], r[14], r[15], pct, r[17], r[18],
                 r[20], r[21], struct))


def main():
    args = build_args()
    nonviral = set(c.strip().upper() for c in args.nonviral_categories.split(",")
                   if c.strip())

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            warn("output directory (cannot create)", args.outdir)
            return 0

    jobs = []            # (cmp_name, cohort, arm_a, arm_b)

    # ---- family 1: viral-only vs hg38-inclusive ----
    for cohort, name_a, name_b in parse_pairs(args):
        run_a = resolve_run(args.runs_root, name_a)
        run_b = resolve_run(args.runs_root, name_b)
        ok = True
        for run in (run_a, run_b):
            if not os.path.isdir(run):
                warn("run directory", run)
                ok = False
            elif not os.path.isdir(os.path.join(run, "results")):
                warn("results directory", os.path.join(run, "results"))
                ok = False
        if not ok:
            continue
        vars_a = headline_variants(run_a)
        vars_b = headline_variants(run_b)
        if not vars_a:
            warn("headline *%s table" % HEADLINE_SUFFIX, os.path.join(run_a, "results"))
            continue
        if not vars_b:
            warn("headline *%s table" % HEADLINE_SUFFIX, os.path.join(run_b, "results"))
            continue
        va, vb = pick_common_variant(vars_a, vars_b)
        if va != vb:
            print("NOTE: %s uses variant '%s' but %s uses '%s'; the two arms are "
                  "not filename-matched." % (os.path.basename(run_a), va,
                                             os.path.basename(run_b), vb))
        arm_a = load_arm("viral_only", run_a, va, vars_a)
        arm_b = load_arm("with_human", run_b, vb, vars_b)
        if arm_a is None or arm_b is None:
            continue
        jobs.append((CMP_HUMAN, cohort, arm_a, arm_b))

    # ---- family 2: all alignments vs primary-only, inside one run ----
    if not args.no_primary_comparison:
        if args.primary_runs:
            cand = [resolve_run(args.runs_root, r) for r in args.primary_runs]
        else:
            cand = discover_primary_runs(args.runs_root)
            if not cand:
                print("NOTE: no run under %s ships both headline filename prefixes; "
                      "the all-vs-primary_only family is empty."
                      % ascii_safe(args.runs_root))
        for run in cand:
            if not os.path.isdir(run):
                warn("run directory", run)
                continue
            v = headline_variants(run)
            if VARIANT_ALL not in v or VARIANT_PRIMARY not in v:
                warn("both headline filename prefixes in run %s"
                     % os.path.basename(run), os.path.join(run, "results"))
                continue
            arm_a = load_arm("all_alignments", run, VARIANT_ALL, v)
            arm_b = load_arm("primary_only", run, VARIANT_PRIMARY, v)
            if arm_a is None or arm_b is None:
                continue
            jobs.append((CMP_PRIMARY, cohort_of_run(run), arm_a, arm_b))

    if not jobs:
        print("Nothing to compare -- no usable run pair was found. Exiting cleanly.")
        return 0

    # ---- anonymisation over the union of all samples actually used ----
    seen_runs = {}
    for _cmp, _coh, arm_a, arm_b in jobs:
        for arm in (arm_a, arm_b):
            for s in arm.table:
                seen_runs.setdefault(s, set()).add(arm.run)
    reals = sorted(seen_runs)
    width = max(2, len(str(len(reals))))
    anon = {}
    groups = {}
    for i, real in enumerate(reals, start=1):
        anon[real] = "S" + str(i).zfill(width)
        groups[real] = group_from_name(real, " ".join(sorted(seen_runs[real])))

    # ---- run the comparisons ----
    cat_rows_all, sample_rows_all, flip_rows_all = [], [], []
    for cmp_name, cohort, arm_a, arm_b in jobs:
        cat_rows, sample_rows, flip_rows, shared, only_a, only_b = compare(
            cmp_name, cohort, arm_a, arm_b, nonviral, args.keep_zero_rows)
        if not shared:
            print("NOTE: %s / %s vs %s share no sample names; nothing comparable."
                  % (cmp_name, arm_a.run, arm_b.run))
        print_block(cmp_name, cohort, arm_a, arm_b, cat_rows,
                    len(shared), only_a, only_b)

        # anonymise sample columns before anything is stored for output
        for row in sample_rows:
            real = row[2]
            row[2] = anon.get(real, "S??")
            row.insert(3, groups.get(real, "NA"))
        for row in flip_rows:
            real = row[2]
            row[2] = anon.get(real, "S??")
            row.insert(3, groups.get(real, "NA"))
            row.extend(["NA", "NA", "NA"])        # verify_* placeholders
            row.append(real)                      # temp: real name, popped later
        cat_rows_all.extend(cat_rows)
        sample_rows_all.extend(sample_rows)
        flip_rows_all.extend((row, arm_a, arm_b) for row in flip_rows)

    # ---- optional samtools verification of pos->neg flips (family 1 only) ----
    verified = 0
    if args.verify_flips > 0:
        cache = {}
        todo = [(row, arm_a, arm_b) for row, arm_a, arm_b in flip_rows_all
                if row[0] == CMP_HUMAN and row[5] == "pos_to_neg"]
        todo.sort(key=lambda t: -t[0][6])         # largest arm-A signal first
        for row, arm_a, arm_b in todo[:args.verify_flips]:
            real = row[-1]
            _path, refmap = refmap_for_run(arm_a.run_dir, args, cache)
            if not refmap:
                continue
            refs = refs_with_signal(arm_a.run_dir, real, row[4], refmap)
            if refs is None:
                warn("idxstats for %s in run %s" % (row[2], arm_a.run),
                     os.path.join(arm_a.run_dir, "results"))
                continue
            if not refs:
                continue
            refs = refs[:args.verify_max_refs]
            bam = os.path.join(arm_b.run_dir, "bam", real + ".bam")
            if not os.path.exists(bam):
                warn("BAM for %s in run %s" % (row[2], arm_b.run),
                     os.path.join(arm_b.run_dir, "bam"))
                continue
            if not (os.path.exists(bam + ".bai") or
                    os.path.exists(bam[:-4] + ".bai")):
                warn("BAM index for %s in run %s" % (row[2], arm_b.run),
                     os.path.join(arm_b.run_dir, "bam"))
                continue
            n = count_unique_best(args.samtools, bam, refs, args.mapq)
            if n is None:
                # never print the real BAM filename: row[2] is the anon id
                warn("samtools output for %s (executable '%s')"
                     % (row[2], args.samtools),
                     os.path.join(arm_b.run_dir, "bam", row[2] + ".bam"))
                continue
            row[12] = len(refs)
            row[13] = n
            row[14] = args.mapq
            verified += 1

    flip_rows_final = []
    for row, _arm_a, _arm_b in flip_rows_all:
        row = list(row)
        row.pop()                                  # drop the real sample name
        flip_rows_final.append(row)

    # ---- write ----
    stem = os.path.join(args.outdir, args.prefix)
    p_cat = stem + "_reference_comparison_by_category.tsv"
    p_samp = stem + "_reference_comparison_by_sample.tsv"
    p_flip = stem + "_reference_comparison_flips.tsv"
    p_key = stem + "_sample_key.tsv"
    write_tsv(p_cat, CAT_HEADER, cat_rows_all)
    write_tsv(p_samp, SAMPLE_HEADER, sample_rows_all)
    write_tsv(p_flip, FLIP_HEADER, flip_rows_final)
    write_sample_key(p_key, anon, groups, seen_runs)

    # ---- stdout summary ----
    n_pn = sum(1 for r in flip_rows_final if r[5] == "pos_to_neg")
    n_np = sum(1 for r in flip_rows_final if r[5] == "neg_to_pos")
    print("")
    print("=" * 78)
    print("SUMMARY  (arm A = permissive, arm B = stricter; delta = B - A)")
    print("  comparisons run      : %d" % len(jobs))
    print("  samples anonymised   : %d (S%s..%s)"
          % (len(reals), "1".zfill(width), anon[reals[-1]] if reals else "-"))
    print("  flip rows pos->neg   : %d   (calls the stricter arm removes)" % n_pn)
    print("  flip rows neg->pos   : %d" % n_np)
    if args.verify_flips > 0:
        print("  flips re-counted     : %d of %d requested (samtools, MAPQ>=%d, AS>XS)"
              % (verified, args.verify_flips, args.mapq))
    print("  wrote %s" % p_cat)
    print("  wrote %s" % p_samp)
    print("  wrote %s" % p_flip)
    print("  wrote %s   <- %s" % (p_key, KEY_WARNING.lstrip("# ")))
    print("")
    print("Caveat: viral-only and hg38-inclusive runs are different experiments, not")
    print("two biological measurements; read the delta as the human-competitor effect")
    print("(and the all-vs-primary_only delta as multi-mapping inflation), not as a")
    print("change in viral load.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
