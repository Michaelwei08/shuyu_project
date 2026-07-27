#!/usr/bin/env python3
"""a1_detection_performance.py - HTLV-1 / HIV-1 detection performance.

WHAT IT COMPUTES
  Reads the headline deduplicated per-sample category counts
  (results/*filtered_category_counts.tsv) from one or more run directories and
  turns the target-virus read counts into detection performance:

    * a 2x2 confusion matrix (TP / FN / FP / TN) at a single calling threshold
      (default: >= 100 reads called positive), with sensitivity, specificity,
      PPV, NPV, F1, accuracy and Youden's J;
    * a rank-based (Mann-Whitney) ROC AUC over the raw read counts, which uses
      every threshold rather than only the ladder points;
    * a threshold sweep over a read ladder (0,1,2,5,10,25,50,100,250,500,1000)
      giving sensitivity/specificity/PPV/NPV/F1 plus TPR/FPR columns so the
      sweep doubles as a coarse ROC table. The currently used default (100) is
      flagged in an is_default column.

  Metrics are reported per run (scope = <run name>) and for the pooled
  cohort (scope = COMBINED). Per-run rows usually contain only one truth class,
  so half their metrics are legitimately NA; COMBINED is the headline row.

TRUTH LABELS
  --labels <tsv> with columns sample_id, expected_status (and an optional
  target column). Values are case-insensitive: positive/pos/1/yes/+ and
  negative/neg/0/no/- ; blank, NA or unknown leaves a sample unlabeled and it
  is excluded from the metrics. sample_id may be either a real sample name or
  the anonymised S01..Snn id from the sample key.

  With no --labels the script falls back to COHORT-AS-LABEL: the targeted HTLV
  cohort (group TCL) is treated as expected HTLV1-positive and the WGS HL
  controls (group HL) as expected HTLV1-negative; for HIV1 the HIV+ WGS group
  is positive and the HL group negative. This is a cohort-level proxy, NOT
  per-sample ground truth, and the script says so loudly on stdout and in the
  label_source column of every output row.

WHAT IT WRITES  (tab separated, into --outdir, pure ASCII, anonymised ids)
  <prefix>_confusion.tsv         2x2 counts + metrics per target and scope
  <prefix>_threshold_sweep.tsv   the read ladder, one row per target/scope/threshold
  <prefix>_label_template.tsv    one row per sample x target with the observed
                                 read count, the cohort proxy call and a blank
                                 expected_status column to fill in
  <prefix>_sample_key.tsv        real -> S01..Snn map. CONTAINS IDENTIFIERS.
                                 Never commit or email this file; it is the
                                 only output that carries real sample names.

NOTES
  Standard library only; no figures, so matplotlib is not imported. No network
  access. Missing runs, missing counts files and missing categories produce a
  "WARN: ... skipping" line and the script continues, exiting 0.

EXAMPLE
  python3 a1_detection_performance.py \
      --runs-root /path/to/runs \
      --run targeted_htlv_hg38_refseq_mapq_human60_viral40_coord \
      --run wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord \
      --targets HTLV1,HIV1 --threshold 100 \
      --outdir /path/to/runs/reply_to_shuyu_primary_only/a1_out
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import re
import os
import sys

RUNS_ROOT = "/path/to/runs"

# Default cohort: the current-reference targeted run (HTLV-positive cohort)
# plus the current-reference WGS run (HIV+ and HL controls). Sample sets are
# disjoint, so they can be pooled into one COMBINED confusion matrix.
DEFAULT_RUNS = [
    "targeted_htlv_hg38_refseq_mapq_human60_viral40_coord",
    "wgs_hiv_hl_hg38_refseq_mapq_human60_viral40_coord",
]

DEFAULT_TARGETS = ["HTLV1", "HIV1"]
DEFAULT_LADDER = [0, 1, 2, 5, 10, 25, 50, 100, 250, 500, 1000]
DEFAULT_THRESHOLD = 100

# Cohort-as-label proxy: target category -> {group: expected status}.
# Groups not listed for a target stay unlabeled and are dropped from metrics.
COHORT_PROXY = {
    "HTLV1": {"TCL": "POSITIVE", "HL": "NEGATIVE"},
    "HTLV2": {"TCL": "POSITIVE", "HL": "NEGATIVE"},
    "HIV1": {"HIV": "POSITIVE", "HL": "NEGATIVE"},
    "HIV2": {"HIV": "POSITIVE", "HL": "NEGATIVE"},
}

POS_WORDS = {"positive", "pos", "p", "1", "yes", "y", "true", "t", "+", "case"}
NEG_WORDS = {"negative", "neg", "n", "0", "no", "false", "f", "-", "control"}
BLANK_WORDS = {"", "na", "n/a", "nan", "none", "null", "unknown", "unk", "?", "."}

TODAY = datetime.date.today().strftime("%Y-%m-%d")


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def warn(what, path):
    print("WARN: %s missing at %s, skipping" % (what, path))


def group_label(sample_name, run_name=""):
    """Group from the real sample name; run name only as a TCL fallback."""
    # The sample-specific label is HIV/HL immediately followed by a digit
    # (HIV<ID>, HL<ID>). Match that first and case-sensitively: the WGS cohort
    # prefix "wgs_60samples_hiv_hl_" is lowercase, and an upper()/lower()
    # test for "_HIV" would otherwise label every WGS sample HIV and leave
    # the HL group empty.
    _m = re.search(r"(?:^|_)(HIV|HL)[0-9]", sample_name or "")
    if _m:
        return _m.group(1)
    up = (sample_name or "").upper()
    if "_HIV" in up:
        return "HIV"
    if "_HL" in up:
        return "HL"
    if "TCL" in up or "TARGETED_HTLV" in up:
        return "TCL"
    run_up = (run_name or "").upper()
    if "TARGETED_HTLV" in run_up or "TCL" in run_up:
        return "TCL"
    return "NA"


def anon_ids(real_names):
    """S01..Snn, numbered by sorted real name. Widened past 99 samples."""
    names = sorted(set(real_names))
    width = max(2, len(str(len(names))))
    return {name: "S" + str(i + 1).zfill(width) for i, name in enumerate(names)}


def fmt(x, nd=4):
    return "NA" if x is None else ("%.*f" % (nd, x))


def as_int(text):
    text = (text or "").strip()
    if text in ("", "NA", "na", "."):
        return 0
    try:
        return int(text)
    except ValueError:
        try:
            return int(round(float(text)))
        except ValueError:
            return 0


def read_rows(path, delimiter="\t"):
    """DictReader with '#' comment lines stripped and header keys normalised."""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    if not lines:
        return [], []
    rdr = csv.reader(lines, delimiter=delimiter)
    header = [h.strip() for h in next(rdr)]
    rows = []
    for parts in rdr:
        if not parts or all(not p.strip() for p in parts):
            continue
        row = {}
        for i, key in enumerate(header):
            row[key] = parts[i] if i < len(parts) else ""
        rows.append(row)
    return header, rows


def write_tsv(path, header, rows, comments=()):
    with open(path, "w", encoding="ascii", newline="") as fh:
        for c in comments:
            fh.write("# " + c + "\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(h, "")) for h in header) + "\n")


# --------------------------------------------------------------------------- #
# input discovery
# --------------------------------------------------------------------------- #
def find_counts_file(run_dir):
    """results/*filtered_category_counts.tsv (handles the primary_only_ prefix).

    The glob deliberately does not match filtered_record_category_counts.tsv,
    raw_idxstats_category_counts.tsv or dedup_removed_category_counts.tsv.
    """
    res = os.path.join(run_dir, "results")
    if not os.path.isdir(res):
        warn("results/ directory", res)
        return None
    hits = sorted(glob.glob(os.path.join(res, "*filtered_category_counts.tsv")))
    hits = [h for h in hits if not os.path.basename(h).startswith("filtered_record")]
    if not hits:
        warn("filtered_category_counts.tsv", os.path.join(res, "*filtered_category_counts.tsv"))
        return None
    exact = [h for h in hits if os.path.basename(h) == "filtered_category_counts.tsv"]
    chosen = exact[0] if exact else hits[0]
    if len(hits) > 1:
        print("NOTE: %d filtered count tables in %s; using %s"
              % (len(hits), res, os.path.basename(chosen)))
    return chosen


def load_counts(path):
    """-> (categories, {real_sample: {category: reads}})."""
    header, rows = read_rows(path)
    if not header:
        warn("header in category counts table", path)
        return [], {}
    key = header[0]
    if key.strip().lower() != "sample":
        print("NOTE: first column of %s is '%s', treating it as the sample column"
              % (path, key))
    cats = [h for h in header[1:] if h.strip()]
    out = {}
    for row in rows:
        name = (row.get(key) or "").strip()
        if not name:
            continue
        out[name] = dict((c, as_int(row.get(c))) for c in cats)
    return cats, out


def load_labels(path, targets):
    """-> ({(target_or_'', id_lower): status}, n_rows, n_bad)."""
    header, rows = read_rows(path)
    if not header:
        warn("header in labels table", path)
        return {}, 0, 0
    lower = dict((h.strip().lower(), h) for h in header)

    def pick(*names):
        for n in names:
            if n in lower:
                return lower[n]
        return None

    c_id = pick("sample_id", "sample", "id", "sample_name")
    c_st = pick("expected_status", "status", "expected", "truth", "label")
    c_tg = pick("target", "virus", "category")
    if not c_id or not c_st:
        print("WARN: labels table %s lacks sample_id/expected_status columns, "
              "falling back to cohort-as-label" % path)
        return {}, 0, 0
    labels, bad = {}, 0
    for row in rows:
        sid = (row.get(c_id) or "").strip()
        if not sid:
            continue
        raw = (row.get(c_st) or "").strip().lower()
        if raw in BLANK_WORDS:
            continue
        if raw in POS_WORDS:
            status = "POSITIVE"
        elif raw in NEG_WORDS:
            status = "NEGATIVE"
        else:
            bad += 1
            continue
        tgt = (row.get(c_tg) or "").strip().upper() if c_tg else ""
        if tgt and tgt not in targets:
            continue
        labels[(tgt, sid.lower())] = status
    return labels, len(rows), bad


def resolve_label(labels, target, real_name, anon_id):
    """Target-specific label first, then a target-agnostic one; real id first."""
    for tgt in (target.upper(), ""):
        for key in (real_name.lower(), anon_id.lower()):
            if (tgt, key) in labels:
                return labels[(tgt, key)]
    return None


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def ratio(num, den):
    return (float(num) / den) if den else None


def confusion(pairs, threshold):
    """pairs = [(reads, status)] -> (tp, fn, fp, tn). reads >= threshold = call."""
    tp = fn = fp = tn = 0
    for reads, status in pairs:
        called = reads >= threshold
        if status == "POSITIVE":
            if called:
                tp += 1
            else:
                fn += 1
        else:
            if called:
                fp += 1
            else:
                tn += 1
    return tp, fn, fp, tn


def metrics(tp, fn, fp, tn):
    sens = ratio(tp, tp + fn)
    spec = ratio(tn, tn + fp)
    ppv = ratio(tp, tp + fp)
    npv = ratio(tn, tn + fn)
    acc = ratio(tp + tn, tp + fn + fp + tn)
    f1 = ratio(2 * tp, 2 * tp + fp + fn)
    j = (sens + spec - 1.0) if (sens is not None and spec is not None) else None
    fpr = (1.0 - spec) if spec is not None else None
    return {"sensitivity": sens, "specificity": spec, "ppv": ppv, "npv": npv,
            "accuracy": acc, "f1": f1, "youden_j": j, "fpr": fpr}


def auc_rank(pairs):
    """Mann-Whitney ROC AUC over the raw counts, ties counted as 0.5."""
    pos = [r for r, s in pairs if s == "POSITIVE"]
    neg = [r for r, s in pairs if s == "NEGATIVE"]
    if not pos or not neg:
        return None
    wins = 0.0
    for p in pos:
        for n in neg:
            if p > n:
                wins += 1.0
            elif p == n:
                wins += 0.5
    return wins / (len(pos) * len(neg))


def auc_trapezoid(points):
    """Trapezoid over swept (fpr, tpr) points, with (0,0) and (1,1) added."""
    pts = [p for p in points if p[0] is not None and p[1] is not None]
    if not pts:
        return None
    pts = sorted(set(pts + [(0.0, 0.0), (1.0, 1.0)]))
    area = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        area += (x1 - x0) * (y0 + y1) / 2.0
    return area


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def build_parser():
    ap = argparse.ArgumentParser(
        description="HTLV-1 / HIV-1 detection performance from filtered category counts.")
    ap.add_argument("--runs-root", default=RUNS_ROOT,
                    help="root holding the run directories (default: %s)" % RUNS_ROOT)
    ap.add_argument("--run", action="append", default=None, metavar="RUN",
                    help="run directory name under --runs-root, or an absolute "
                         "path. Repeatable. Default: %s" % ", ".join(DEFAULT_RUNS))
    ap.add_argument("--targets", default=",".join(DEFAULT_TARGETS),
                    help="comma-separated target categories (default: %s)"
                         % ",".join(DEFAULT_TARGETS))
    ap.add_argument("--labels", default=None,
                    help="truth table tsv: sample_id, expected_status[, target]. "
                         "Omit to use the cohort-as-label proxy.")
    ap.add_argument("--threshold", type=int, default=DEFAULT_THRESHOLD,
                    help="reads >= threshold is a positive call (default: %d)"
                         % DEFAULT_THRESHOLD)
    ap.add_argument("--thresholds",
                    default=",".join(str(t) for t in DEFAULT_LADDER),
                    help="sweep ladder (default: %s)"
                         % ",".join(str(t) for t in DEFAULT_LADDER))
    ap.add_argument("--outdir", default="suite_out",
                    help="output directory (default: suite_out)")
    ap.add_argument("--prefix", default="detection",
                    help="output filename prefix (default: detection)")
    ap.add_argument("--samtools", default="samtools",
                    help="samtools executable, accepted for suite-wide CLI "
                         "consistency; this module reads only text tables and "
                         "never shells out")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    targets = [t.strip().upper() for t in args.targets.split(",") if t.strip()]
    if not targets:
        print("WARN: no targets requested, nothing to do")
        return 0
    ladder = sorted(set(int(t) for t in args.thresholds.replace(" ", "").split(",")
                        if t not in ("",)))
    if args.threshold not in ladder:
        ladder = sorted(set(ladder + [args.threshold]))

    run_names = args.run if args.run else list(DEFAULT_RUNS)
    outdir = args.outdir
    try:
        os.makedirs(outdir)
    except OSError:
        if not os.path.isdir(outdir):
            print("WARN: could not create outdir at %s, skipping all output" % outdir)
            return 0

    # ---- collect per-run, per-sample counts -------------------------------- #
    # records[run_label] = {real_sample: {category: reads}}
    records, run_order = {}, []
    for name in run_names:
        run_dir = name if os.path.isabs(name) else os.path.join(args.runs_root, name)
        label = os.path.basename(run_dir.rstrip("/\\")) or run_dir
        if not os.path.isdir(run_dir):
            warn("run directory", run_dir)
            continue
        path = find_counts_file(run_dir)
        if not path:
            continue
        cats, per_sample = load_counts(path)
        if not per_sample:
            warn("per-sample rows in category counts table", path)
            continue
        missing = [t for t in targets if t not in cats]
        for t in missing:
            print("WARN: category %s missing at %s, skipping (this run cannot "
                  "score %s)" % (t, path, t))
        records[label] = per_sample
        run_order.append(label)
        print("OK  : %-58s %3d samples, categories: %s"
              % (label, len(per_sample), ",".join(cats) if cats else "none"))

    if not records:
        print("WARN: usable run data missing at %s, skipping (no outputs written)"
              % args.runs_root)
        return 0

    # ---- anonymise --------------------------------------------------------- #
    all_names = set()
    for per_sample in records.values():
        all_names.update(per_sample)
    anon = anon_ids(all_names)

    groups, sample_runs = {}, {}
    for label in run_order:
        for name in records[label]:
            sample_runs.setdefault(name, []).append(label)
            g = group_label(name, label)
            if groups.get(name, "NA") == "NA":
                groups[name] = g

    dupes = [n for n, rs in sample_runs.items() if len(rs) > 1]
    if dupes:
        print("NOTE: %d sample name(s) occur in more than one run; the COMBINED "
              "scope keeps the first run in --run order to avoid double counting."
              % len(dupes))

    key_path = os.path.join(outdir, "%s_sample_key.tsv" % args.prefix)
    write_tsv(key_path,
              ["anon_id", "real_sample_id", "group", "runs"],
              [{"anon_id": anon[n], "real_sample_id": n, "group": groups[n],
                "runs": ";".join(sample_runs[n])} for n in sorted(all_names)],
              comments=["CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
                        "a1_detection_performance.py  generated %s" % TODAY,
                        "maps the anonymised ids used in every other output "
                        "back to real sample names"])

    # ---- truth labels ------------------------------------------------------ #
    labels, n_label_rows, n_bad = {}, 0, 0
    if args.labels:
        if os.path.exists(args.labels):
            labels, n_label_rows, n_bad = load_labels(args.labels, targets)
            if n_bad:
                print("NOTE: %d label row(s) had an unrecognised expected_status "
                      "and were ignored" % n_bad)
        else:
            warn("labels table", args.labels)
    label_source = "user_labels" if labels else "cohort_proxy"
    proxy = not labels

    # ---- per-sample table (drives every output) ---------------------------- #
    # per_scope[(target, scope)] = [(reads, status)]
    per_scope, template_rows = {}, []
    combined_seen = set()
    for target in targets:
        for label in run_order:
            for name in sorted(records[label]):
                if target not in records[label][name]:
                    continue
                reads = records[label][name][target]
                aid, grp = anon[name], groups[name]
                proxy_status = COHORT_PROXY.get(target, {}).get(grp, "NA")
                truth = (resolve_label(labels, target, name, aid) if labels
                         else (proxy_status if proxy_status in ("POSITIVE", "NEGATIVE")
                               else None))
                if truth:
                    per_scope.setdefault((target, label), []).append((reads, truth))
                    if (target, name) not in combined_seen:
                        combined_seen.add((target, name))
                        per_scope.setdefault((target, "COMBINED"), []).append((reads, truth))
                template_rows.append({
                    "sample_id": aid, "group": grp, "run": label,
                    "target": target, "reads": reads,
                    "cohort_proxy_status": proxy_status,
                    "resolved_expected_status": truth or "",
                    "expected_status": "",
                })

    write_tsv(os.path.join(outdir, "%s_label_template.tsv" % args.prefix),
              ["sample_id", "group", "run", "target", "reads",
               "cohort_proxy_status", "resolved_expected_status",
               "expected_status"],
              template_rows,
              comments=["a1_detection_performance.py  generated %s" % TODAY,
                        "anonymised ids only - resolve them with %s_sample_key.tsv"
                        % args.prefix,
                        "fill in expected_status with positive / negative and pass "
                        "this file back via --labels; leave it blank for unknown",
                        "the target column is honoured on re-read, so one sample "
                        "can carry a different status per virus",
                        "cohort_proxy_status is the cohort-level guess, NOT ground "
                        "truth; resolved_expected_status is what this run scored"])

    # ---- confusion + sweep ------------------------------------------------- #
    scopes = ["COMBINED"] + run_order
    conf_rows, sweep_rows, summary = [], [], []
    for target in targets:
        for scope in scopes:
            pairs = per_scope.get((target, scope), [])
            if not pairs:
                continue
            n_pos = sum(1 for _r, s in pairs if s == "POSITIVE")
            n_neg = len(pairs) - n_pos
            n_total_scope = len(records[scope]) if scope in records else len(all_names)
            roc_pts = []
            for thr in ladder:
                tp, fn, fp, tn = confusion(pairs, thr)
                m = metrics(tp, fn, fp, tn)
                roc_pts.append((m["fpr"], m["sensitivity"]))
                sweep_rows.append({
                    "target": target, "scope": scope, "label_source": label_source,
                    "threshold": thr,
                    "is_default": "yes" if thr == args.threshold else "",
                    "TP": tp, "FN": fn, "FP": fp, "TN": tn,
                    "tpr_sensitivity": fmt(m["sensitivity"]),
                    "fpr_1_minus_specificity": fmt(m["fpr"]),
                    "specificity": fmt(m["specificity"]),
                    "ppv": fmt(m["ppv"]), "npv": fmt(m["npv"]),
                    "f1": fmt(m["f1"]), "youden_j": fmt(m["youden_j"]),
                })
            tp, fn, fp, tn = confusion(pairs, args.threshold)
            m = metrics(tp, fn, fp, tn)
            a_rank = auc_rank(pairs)
            conf_rows.append({
                "target": target, "scope": scope, "label_source": label_source,
                "threshold": args.threshold,
                "n_samples_in_scope": n_total_scope,
                "n_scored": len(pairs),
                "n_expected_pos": n_pos, "n_expected_neg": n_neg,
                "n_unlabeled": max(0, n_total_scope - len(pairs)),
                "TP": tp, "FN": fn, "FP": fp, "TN": tn,
                "sensitivity": fmt(m["sensitivity"]),
                "specificity": fmt(m["specificity"]),
                "ppv": fmt(m["ppv"]), "npv": fmt(m["npv"]),
                "f1": fmt(m["f1"]), "accuracy": fmt(m["accuracy"]),
                "youden_j": fmt(m["youden_j"]),
                "auc_rank_all_thresholds": fmt(a_rank),
                "auc_trapezoid_over_ladder": fmt(auc_trapezoid(roc_pts)),
            })
            if scope == "COMBINED":
                best = None
                for row in sweep_rows:
                    if row["target"] != target or row["scope"] != "COMBINED":
                        continue
                    if row["f1"] == "NA":
                        continue
                    val = float(row["f1"])
                    if best is None or val > best[0]:
                        best = (val, row["threshold"])
                summary.append((target, n_pos, n_neg, tp, fn, fp, tn, m, a_rank, best))

    conf_header = ["target", "scope", "label_source", "threshold",
                   "n_samples_in_scope", "n_scored", "n_expected_pos",
                   "n_expected_neg", "n_unlabeled", "TP", "FN", "FP", "TN",
                   "sensitivity", "specificity", "ppv", "npv", "f1", "accuracy",
                   "youden_j", "auc_rank_all_thresholds",
                   "auc_trapezoid_over_ladder"]
    sweep_header = ["target", "scope", "label_source", "threshold", "is_default",
                    "TP", "FN", "FP", "TN", "tpr_sensitivity",
                    "fpr_1_minus_specificity", "specificity", "ppv", "npv",
                    "f1", "youden_j"]
    common = ["a1_detection_performance.py  generated %s" % TODAY,
              "positive call = target reads >= threshold, from the deduplicated "
              "filtered_category_counts.tsv of each run",
              "scope COMBINED pools all runs (one row per sample); per-run scopes "
              "often hold a single truth class, so half their metrics are NA",
              "anonymised sample ids only; no identifiers in this file"]
    if proxy:
        common.append("label_source=cohort_proxy - CAVEAT: cohort-level proxy "
                      "labels, NOT per-sample ground truth")
    write_tsv(os.path.join(outdir, "%s_confusion.tsv" % args.prefix),
              conf_header, conf_rows, comments=common)
    write_tsv(os.path.join(outdir, "%s_threshold_sweep.tsv" % args.prefix),
              sweep_header, sweep_rows,
              comments=common + ["is_default marks the threshold currently used "
                                 "elsewhere in the suite (%d reads)" % args.threshold])

    # ---- stdout ------------------------------------------------------------ #
    print("")
    print("=" * 78)
    print("DETECTION PERFORMANCE  (%s)" % TODAY)
    print("runs scored : %d   samples seen: %d   thresholds: %s"
          % (len(run_order), len(all_names),
             ",".join(str(t) for t in ladder)))
    print("call rule   : target reads >= %d" % args.threshold)
    print("label source: %s%s" % (label_source,
                                  ("  (%d rows read)" % n_label_rows) if labels else ""))
    if proxy:
        print("")
        print("!" * 78)
        print("!! CAVEAT: no --labels supplied, so labels are a COHORT-LEVEL PROXY.")
        print("!! The targeted HTLV cohort is assumed HTLV1-positive and the WGS HL")
        print("!! controls HTLV1-negative (HIV+ group positive for HIV1). This is")
        print("!! NOT per-sample ground truth: sensitivity/PPV computed this way")
        print("!! measure agreement with a cohort assignment, not with a clinical")
        print("!! assay. Fill in %s_label_template.tsv and rerun with --labels."
              % args.prefix)
        print("!" * 78)
    print("")
    if not summary:
        print("No target/label overlap, so no metrics could be computed.")
    for target, n_pos, n_neg, tp, fn, fp, tn, m, a_rank, best in summary:
        print("-- %s  (COMBINED scope) --" % target)
        print("   expected positive: %d    expected negative: %d" % (n_pos, n_neg))
        print("   confusion at >=%d reads:  TP=%d  FN=%d  FP=%d  TN=%d"
              % (args.threshold, tp, fn, fp, tn))
        print("   sens=%s  spec=%s  PPV=%s  NPV=%s  F1=%s  J=%s  AUC=%s"
              % (fmt(m["sensitivity"], 3), fmt(m["specificity"], 3),
                 fmt(m["ppv"], 3), fmt(m["npv"], 3), fmt(m["f1"], 3),
                 fmt(m["youden_j"], 3), fmt(a_rank, 3)))
        if best:
            print("   best F1 on the ladder: %s at >=%s reads"
                  % (fmt(best[0], 3), best[1]))
        print("")
    print("wrote %s_confusion.tsv, %s_threshold_sweep.tsv, %s_label_template.tsv"
          % (args.prefix, args.prefix, args.prefix))
    print("wrote %s  <- CONTAINS IDENTIFIERS, do not commit or email" % key_path)
    print("outdir: %s" % os.path.abspath(outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
