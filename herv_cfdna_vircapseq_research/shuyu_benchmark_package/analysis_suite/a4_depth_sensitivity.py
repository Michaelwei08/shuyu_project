#!/usr/bin/env python3
"""A4 -- depth sensitivity of viral detection: 5M-read subsample vs full depth.

WHAT IT COMPUTES
    Compares a read-subsampled run (default wgs_hiv_hl_5m_competitive, 5M reads
    per sample) against the matching full-depth run (default
    wgs_hiv_hl_full_competitive) over the same samples and the same categories.
    Both defaults are the viral-only competitive reference, so the category sets
    are directly comparable.

    Counts are read from the headline deduplicated per-sample table
    results/*filtered_category_counts.tsv (this glob also picks up runs that use
    the "primary_only_" filename prefix). No BAM is touched, so no samtools and
    no figures: standard library only.

    Per category (and per cohort group) it reports
      - total reads at 5M vs at full depth, and the summed read ratio;
      - number of samples with a nonzero call (>= --min-reads) at each depth;
      - detection retention  = (samples positive at BOTH) / (samples positive at
        full depth): the fraction of full-depth detections that survive at 5M;
      - per-sample read-ratio statistics (mean / median / min / max);
      - ratio_vs_expected_median, the per-sample read ratio divided by that
        sample's overall 5M/full total-read ratio. ~1.0 means the category
        simply scales with depth; << 1.0 means it is lost faster than depth.
      - a verdict: SATURATED (retention >= --saturated-min), DEPTH_LIMITED
        (retention <= --depth-limited-max), PARTIAL in between, or NO_SIGNAL.

    If the two runs do not cover the same samples, the intersection is analysed
    and the difference is reported on stdout and in the sample key.

WHAT IT WRITES (tab-separated, into --outdir)
    <prefix>_by_category.tsv   one row per (category, group); group ALL first.
                               Also carries a row for the pseudo-category
                               ALL_CATEGORIES (the per-sample totals).
    <prefix>_by_sample.tsv     long format: one row per (sample, category) with
                               reads at each depth, the ratio, and a status of
                               RETAINED / LOST / GAINED / ABSENT.
    <prefix>_sample_key.tsv    real -> anonymous ID mapping. CONTAINS
                               IDENTIFIERS; every other output uses S01..Snn
                               only (assigned by sorted real sample name).

    Missing inputs never crash the module: it prints
    "WARN: <what> missing at <path>, skipping" and exits 0.

EXAMPLE
    python3 a4_depth_sensitivity.py \
        --runs-root /path/to/runs \
        --sub-run wgs_hiv_hl_5m_competitive \
        --full-run wgs_hiv_hl_full_competitive \
        --outdir ./a4_depth_sensitivity_out
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import statistics
import sys
from datetime import date

CAT_TOTAL = "ALL_CATEGORIES"
GROUP_ALL = "ALL"
GROUP_RANK = {"ALL": 0, "HIV": 1, "HL": 2, "TCL": 3, "NA": 5}

# trailing depth markers a subsampled run may append to a sample name,
# e.g. "..._5M", "...-sub5m", "..._downsampled_5Mreads"
DEPTH_SUFFIX_RE = re.compile(
    r"[._-]+(?:sub|subsample|subsampled|down|downsample|downsampled|ds)?"
    r"[._-]*\d+(?:\.\d+)?[mkg](?:_?reads?)?$")


# ----------------------------------------------------------------- utilities

def warn(what, path):
    """The one warning shape this suite uses for absent inputs."""
    print("WARN: %s missing at %s, skipping" % (what, path))


def ascii_safe(text):
    """Everything written out is forced to pure ASCII."""
    if text is None:
        return ""
    return str(text).encode("ascii", "replace").decode("ascii")


def f4(value):
    return "NA" if value is None else "%.4f" % value


def to_int(cell):
    cell = (cell or "").strip()
    if not cell or cell.upper() in ("NA", "NAN", "NONE", "."):
        return 0
    try:
        return int(cell)
    except ValueError:
        pass
    try:
        return int(round(float(cell)))
    except ValueError:
        return 0


def group_of(sample_name):
    """Cohort label derived from the real sample name (suite-wide rule).

    Matched case-insensitively so every module in the suite agrees:
    "_HIV" -> HIV, "_HL" -> HL, "TCL"/"targeted_htlv" -> TCL, else NA.
    """
    up = (sample_name or "").upper()
    if "_HIV" in up:
        return "HIV"
    if "_HL" in up:
        return "HL"
    if "TARGETED_HTLV" in up or "TCL" in up:
        return "TCL"
    return "NA"


def norm_name(name):
    """Lowercase name with BAM-ish and depth suffixes stripped (fuzzy match)."""
    text = name.strip().lower()
    for suffix in (".bam", ".cram", ".sorted", ".dedup", ".markdup"):
        if text.endswith(suffix):
            text = text[:-len(suffix)]
    for _ in range(3):
        match = DEPTH_SUFFIX_RE.search(text)
        if not match or match.start() == 0:
            break
        text = text[:match.start()]
    return text.strip("._-")


def resolve_run(runs_root, run):
    if os.path.isabs(run) or os.sep in run or (os.altsep and os.altsep in run):
        return run
    return os.path.join(runs_root, run)


# ------------------------------------------------------------- input loading

def load_category_counts(run_dir, pattern, drop_cols):
    """-> (path, categories, {sample: {category: reads}}) or None if absent."""
    if not os.path.isdir(run_dir):
        warn("run directory", run_dir)
        return None
    results = os.path.join(run_dir, "results")
    if not os.path.isdir(results):
        warn("results/ subdirectory", results)
        return None
    candidates = sorted(glob.glob(os.path.join(results, pattern)))
    if not candidates:
        warn("category counts table (%s)" % pattern, results)
        return None
    # prefer the unprefixed filename when both it and primary_only_* exist
    path = sorted(candidates, key=lambda p: (len(os.path.basename(p)), p))[0]
    if len(candidates) > 1:
        print("NOTE: %d category tables matched in %s; using %s"
              % (len(candidates), results, os.path.basename(path)))

    header, body = None, []
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if not row or all(not cell.strip() for cell in row):
                    continue
                if row[0].lstrip().startswith("#"):
                    continue
                if header is None:
                    header = [ascii_safe(cell).strip() for cell in row]
                    continue
                body.append(row)
    except OSError:
        warn("readable category counts table", path)
        return None
    if not header:
        warn("header row in category counts table", path)
        return None

    lowered = [cell.lower() for cell in header]
    sample_idx = lowered.index("sample") if "sample" in lowered else 0
    categories, col_of = [], {}
    for idx, name in enumerate(header):
        if idx == sample_idx or not name:
            continue
        if name.lower() in drop_cols:
            continue
        categories.append(name)
        col_of[name] = idx

    counts, duplicates = {}, 0
    for row in body:
        if sample_idx >= len(row):
            continue
        sample = ascii_safe(row[sample_idx]).strip()
        if not sample or sample == "*":
            continue
        if sample in counts:
            duplicates += 1
        target = counts.setdefault(sample, dict((c, 0) for c in categories))
        for cat in categories:
            idx = col_of[cat]
            if idx < len(row):
                target[cat] += to_int(row[idx])
    if duplicates:
        print("NOTE: %d duplicate sample rows in %s were summed"
              % (duplicates, os.path.basename(path)))
    if not counts:
        warn("sample rows in category counts table", path)
        return None
    return path, categories, counts


def match_samples(names_sub, names_full, mode):
    """-> (pairs[(full_name, sub_name)], mode_used, only_full, only_sub, notes)"""
    set_sub, set_full = set(names_sub), set(names_full)
    exact = sorted(set_sub & set_full)
    if mode == "exact" or (mode == "auto" and exact):
        return ([(n, n) for n in exact], "exact",
                sorted(set_full - set_sub), sorted(set_sub - set_full), [])

    def index(names):
        idx = {}
        for name in names:
            idx.setdefault(norm_name(name), []).append(name)
        return idx

    idx_sub, idx_full = index(set_sub), index(set_full)
    pairs, used_sub, used_full, notes = [], set(), set(), []
    for key in sorted(set(idx_sub) & set(idx_full)):
        if len(idx_sub[key]) > 1 or len(idx_full[key]) > 1:
            notes.append("ambiguous normalised sample key, not matched")
            continue
        pairs.append((idx_full[key][0], idx_sub[key][0]))
        used_full.add(idx_full[key][0])
        used_sub.add(idx_sub[key][0])
    pairs.sort(key=lambda pair: pair[0])
    return (pairs, "normalised", sorted(set_full - used_full),
            sorted(set_sub - used_sub), notes)


# ---------------------------------------------------------------- statistics

def summarise(records, min_reads, saturated_min, depth_limited_max, few_n):
    """Aggregate per-(sample, category) records for one category/group cell."""
    reads_sub = sum(r["reads_sub"] for r in records)
    reads_full = sum(r["reads_full"] for r in records)
    n_pos_full = sum(1 for r in records if r["reads_full"] >= min_reads)
    n_pos_sub = sum(1 for r in records if r["reads_sub"] >= min_reads)
    n_retained = sum(1 for r in records if r["status"] == "RETAINED")
    n_lost = sum(1 for r in records if r["status"] == "LOST")
    n_gained = sum(1 for r in records if r["status"] == "GAINED")

    ratios = [r["ratio"] for r in records if r["ratio"] is not None]
    expected = [r["ratio_vs_expected"] for r in records
                if r["ratio_vs_expected"] is not None]

    retention = float(n_retained) / n_pos_full if n_pos_full else None
    if n_pos_full == 0 and n_pos_sub == 0:
        verdict = "NO_SIGNAL"
    elif n_pos_full == 0:
        verdict = "SUB_ONLY_SIGNAL"
    elif retention >= saturated_min:
        verdict = "SATURATED"
    elif retention <= depth_limited_max:
        verdict = "DEPTH_LIMITED"
    else:
        verdict = "PARTIAL"

    notes = []
    if 0 < n_pos_full < few_n:
        notes.append("FEW_POSITIVES")
    if n_gained:
        notes.append("GAINED_AT_SUBSAMPLE")
    exp_median = statistics.median(expected) if expected else None
    if exp_median is not None:
        if exp_median < 0.75:
            notes.append("SUBLINEAR_YIELD")
        elif exp_median > 1.25:
            notes.append("SUPERLINEAR_YIELD")

    return {
        "n_samples": len(records),
        "reads_sub": reads_sub,
        "reads_full": reads_full,
        "reads_ratio_total": (float(reads_sub) / reads_full) if reads_full else None,
        "n_pos_full": n_pos_full,
        "n_pos_sub": n_pos_sub,
        "n_retained": n_retained,
        "n_lost": n_lost,
        "n_gained": n_gained,
        "detection_retention": retention,
        "n_ratio_samples": len(ratios),
        "ratio_mean": statistics.mean(ratios) if ratios else None,
        "ratio_median": statistics.median(ratios) if ratios else None,
        "ratio_min": min(ratios) if ratios else None,
        "ratio_max": max(ratios) if ratios else None,
        "ratio_vs_expected_median": exp_median,
        "verdict": verdict,
        "note": ";".join(notes) if notes else "-",
    }


# ------------------------------------------------------------------- writing

def write_tsv(path, header, rows):
    with open(path, "w", newline="", encoding="ascii", errors="replace") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([ascii_safe(cell) for cell in row])


def write_sample_key(path, key_rows, today):
    with open(path, "w", newline="", encoding="ascii", errors="replace") as handle:
        handle.write("# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL\n")
        handle.write("# generated %s by a4_depth_sensitivity.py\n" % today)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["anon_sample", "real_sample_full_run",
                         "real_sample_sub_run", "group", "in_full_run",
                         "in_sub_run", "analysed"])
        for row in key_rows:
            writer.writerow([ascii_safe(cell) for cell in row])


# ---------------------------------------------------------------------- main

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Depth sensitivity: subsampled run vs full-depth run.")
    parser.add_argument("--runs-root", default="/path/to/runs",
                        help="directory holding the run directories")
    parser.add_argument("--sub-run", default="wgs_hiv_hl_5m_competitive",
                        help="subsampled run (name under --runs-root, or a path)")
    parser.add_argument("--full-run", default="wgs_hiv_hl_full_competitive",
                        help="full-depth run (name under --runs-root, or a path)")
    parser.add_argument("--counts-glob", default="*filtered_category_counts.tsv",
                        help="counts table to read inside <run>/results/")
    parser.add_argument("--outdir", default="./a4_depth_sensitivity_out",
                        help="directory for the .tsv outputs")
    parser.add_argument("--prefix", default="depth_sensitivity",
                        help="output filename prefix")
    parser.add_argument("--min-reads", type=int, default=1,
                        help="reads needed to call a category present")
    parser.add_argument("--saturated-min", type=float, default=0.95,
                        help="detection retention at or above this = SATURATED")
    parser.add_argument("--depth-limited-max", type=float, default=0.60,
                        help="detection retention at or below this = DEPTH_LIMITED")
    parser.add_argument("--few-positives", type=int, default=3,
                        help="flag FEW_POSITIVES below this many full-depth calls")
    parser.add_argument("--match-mode", choices=("auto", "exact", "normalised"),
                        default="auto",
                        help="how sample names are matched across the two runs")
    parser.add_argument("--drop-columns", default="total,total_reads,sum",
                        help="comma-separated counts columns that are not categories")
    parser.add_argument("--samtools", default="samtools",
                        help="samtools executable (accepted for suite-wide CLI "
                             "consistency; this module reads no BAMs)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    today = date.today().isoformat()
    drop_cols = set(c.strip().lower() for c in args.drop_columns.split(",")
                    if c.strip())

    sub_dir = resolve_run(args.runs_root, args.sub_run)
    full_dir = resolve_run(args.runs_root, args.full_run)
    if not os.path.isdir(args.runs_root) and not (
            os.path.isdir(sub_dir) or os.path.isdir(full_dir)):
        warn("runs root", args.runs_root)
        return 0

    loaded_sub = load_category_counts(sub_dir, args.counts_glob, drop_cols)
    loaded_full = load_category_counts(full_dir, args.counts_glob, drop_cols)
    if loaded_sub is None or loaded_full is None:
        print("Nothing to compare; a4_depth_sensitivity exits without output.")
        return 0
    sub_path, cats_sub, counts_sub = loaded_sub
    full_path, cats_full, counts_full = loaded_full

    try:
        os.makedirs(args.outdir, exist_ok=True)
    except OSError:
        warn("writable output directory", args.outdir)
        return 0

    pairs, mode_used, only_full, only_sub, match_notes = match_samples(
        list(counts_sub), list(counts_full), args.match_mode)

    # ---- anonymise the union of both runs, sorted by real sample name -------
    canonical = sorted(set([p[0] for p in pairs]) | set(only_full) | set(only_sub))
    width = max(2, len(str(len(canonical))))
    anon_of = {}
    for idx, name in enumerate(canonical, start=1):
        anon_of[name] = "S%s" % str(idx).zfill(width)
    sub_real_of = dict((full_name, sub_name) for full_name, sub_name in pairs)

    key_rows = []
    for name in canonical:
        matched = name in sub_real_of
        in_full = name in counts_full
        in_sub = matched or (name in counts_sub)
        key_rows.append([anon_of[name],
                         name if in_full else "-",
                         sub_real_of.get(name, name if name in counts_sub else "-"),
                         group_of(name),
                         "yes" if in_full else "no",
                         "yes" if in_sub else "no",
                         "yes" if matched else "no"])
    key_path = os.path.join(args.outdir, "%s_sample_key.tsv" % args.prefix)
    write_sample_key(key_path, key_rows, today)

    # ---- category set -------------------------------------------------------
    categories = list(cats_full) + [c for c in cats_sub if c not in cats_full]
    cat_note = {}
    for cat in categories:
        flags = []
        if cat not in cats_sub:
            flags.append("CATEGORY_ABSENT_IN_SUB_RUN")
        if cat not in cats_full:
            flags.append("CATEGORY_ABSENT_IN_FULL_RUN")
        cat_note[cat] = flags

    # ---- per-(sample, category) records ------------------------------------
    records = []
    for full_name, sub_name in pairs:
        anon = anon_of[full_name]
        grp = group_of(full_name)
        row_full = counts_full.get(full_name, {})
        row_sub = counts_sub.get(sub_name, {})
        total_full = sum(row_full.get(c, 0) for c in categories)
        total_sub = sum(row_sub.get(c, 0) for c in categories)
        expected = (float(total_sub) / total_full) if total_full else None

        for cat in [CAT_TOTAL] + categories:
            if cat == CAT_TOTAL:
                reads_full, reads_sub = total_full, total_sub
            else:
                reads_full = row_full.get(cat, 0)
                reads_sub = row_sub.get(cat, 0)
            pos_full = reads_full >= args.min_reads
            pos_sub = reads_sub >= args.min_reads
            if pos_full and pos_sub:
                status = "RETAINED"
            elif pos_full:
                status = "LOST"
            elif pos_sub:
                status = "GAINED"
            else:
                status = "ABSENT"
            ratio = (float(reads_sub) / reads_full) if reads_full > 0 else None
            rve = None
            if ratio is not None and expected not in (None, 0.0):
                rve = ratio / expected
            records.append({
                "anon": anon, "group": grp, "category": cat,
                "reads_sub": reads_sub, "reads_full": reads_full,
                "ratio": ratio, "ratio_vs_expected": rve,
                "pos_full": pos_full, "pos_sub": pos_sub, "status": status,
            })

    # ---- by-sample table ----------------------------------------------------
    sample_header = ["anon_sample", "group", "category", "reads_5m",
                     "reads_full", "read_ratio", "ratio_vs_expected",
                     "detected_5m", "detected_full", "status"]
    sample_rows = []
    for rec in sorted(records, key=lambda r: (r["anon"],
                                             0 if r["category"] == CAT_TOTAL else 1,
                                             r["category"])):
        sample_rows.append([
            rec["anon"], rec["group"], rec["category"], rec["reads_sub"],
            rec["reads_full"], f4(rec["ratio"]), f4(rec["ratio_vs_expected"]),
            "1" if rec["pos_sub"] else "0", "1" if rec["pos_full"] else "0",
            rec["status"]])
    sample_path = os.path.join(args.outdir, "%s_by_sample.tsv" % args.prefix)
    write_tsv(sample_path, sample_header, sample_rows)

    # ---- by-category table --------------------------------------------------
    groups_present = sorted(set(r["group"] for r in records),
                            key=lambda g: (GROUP_RANK.get(g, 4), g))
    cat_header = ["category", "group", "n_samples", "reads_5m", "reads_full",
                  "reads_ratio_total", "n_pos_full", "n_pos_5m", "n_retained",
                  "n_lost", "n_gained", "detection_retention",
                  "n_ratio_samples", "ratio_mean", "ratio_median", "ratio_min",
                  "ratio_max", "ratio_vs_expected_median", "verdict", "note"]
    cat_rows, summary_all = [], {}
    for cat in [CAT_TOTAL] + categories:
        cat_records = [r for r in records if r["category"] == cat]
        for grp in [GROUP_ALL] + groups_present:
            subset = (cat_records if grp == GROUP_ALL
                      else [r for r in cat_records if r["group"] == grp])
            if not subset:
                continue
            stats = summarise(subset, args.min_reads, args.saturated_min,
                              args.depth_limited_max, args.few_positives)
            note_bits = [n for n in cat_note.get(cat, [])]
            if stats["note"] != "-":
                note_bits.append(stats["note"])
            note = ";".join(note_bits) if note_bits else "-"
            if grp == GROUP_ALL:
                summary_all[cat] = stats
            cat_rows.append([
                cat, grp, stats["n_samples"], stats["reads_sub"],
                stats["reads_full"], f4(stats["reads_ratio_total"]),
                stats["n_pos_full"], stats["n_pos_sub"], stats["n_retained"],
                stats["n_lost"], stats["n_gained"],
                f4(stats["detection_retention"]), stats["n_ratio_samples"],
                f4(stats["ratio_mean"]), f4(stats["ratio_median"]),
                f4(stats["ratio_min"]), f4(stats["ratio_max"]),
                f4(stats["ratio_vs_expected_median"]), stats["verdict"], note])
    cat_rows.sort(key=lambda row: (
        0 if row[0] == CAT_TOTAL else 1,
        GROUP_RANK.get(row[1], 4),
        -int(row[4]), row[0]))
    cat_path = os.path.join(args.outdir, "%s_by_category.tsv" % args.prefix)
    write_tsv(cat_path, cat_header, cat_rows)

    # ---- stdout summary -----------------------------------------------------
    print("")
    print("A4 depth sensitivity (subsample vs full depth)   %s" % today)
    print("  subsampled run : %s" % sub_dir)
    print("                   %s (%d samples)" % (sub_path, len(counts_sub)))
    print("  full-depth run : %s" % full_dir)
    print("                   %s (%d samples)" % (full_path, len(counts_full)))
    print("  categories     : %d (%s)" % (len(categories), ", ".join(categories)))
    print("  sample matching: %s; %d samples analysed" % (mode_used, len(pairs)))
    for note in match_notes:
        print("  NOTE: %s" % note)
    if only_full or only_sub:
        print("  SAMPLE SETS DIFFER: intersected %d samples; %d only in the "
              "full-depth run, %d only in the subsampled run."
              % (len(pairs), len(only_full), len(only_sub)))
        dropped = [anon_of[n] for n in sorted(set(only_full) | set(only_sub))]
        print("  dropped (anon): %s" % ", ".join(dropped))
    else:
        print("  sample sets are identical across the two runs.")
    if not pairs:
        print("  WARN: no shared samples, tables written with headers only.")

    total = summary_all.get(CAT_TOTAL)
    if total and total["reads_ratio_total"] is not None:
        print("  observed depth fraction (all categories): %.4f "
              "(%d / %d reads)" % (total["reads_ratio_total"],
                                   total["reads_sub"], total["reads_full"]))
    print("")
    print("  %-14s %14s %14s %7s %8s %9s  %s"
          % ("category", "reads_full", "reads_5m", "ratio", "pos f>5m",
             "retention", "verdict"))
    for cat in sorted(summary_all,
                      key=lambda c: (0 if c == CAT_TOTAL else 1,
                                     -summary_all[c]["reads_full"], c)):
        st = summary_all[cat]
        print("  %-14s %14d %14d %7s %4d>%-3d %9s  %s"
              % (cat[:14], st["reads_full"], st["reads_sub"],
                 f4(st["reads_ratio_total"]), st["n_pos_full"],
                 st["n_pos_sub"], f4(st["detection_retention"]),
                 st["verdict"]))

    limited = [c for c in summary_all
               if summary_all[c]["verdict"] == "DEPTH_LIMITED" and c != CAT_TOTAL]
    partial = [c for c in summary_all
               if summary_all[c]["verdict"] == "PARTIAL" and c != CAT_TOTAL]
    saturated = [c for c in summary_all
                 if summary_all[c]["verdict"] == "SATURATED" and c != CAT_TOTAL]
    nosig = [c for c in summary_all
             if summary_all[c]["verdict"] == "NO_SIGNAL" and c != CAT_TOTAL]
    print("")
    print("  DEPTH-LIMITED (detections lost at 5M): %s"
          % (", ".join(sorted(limited)) if limited else "none"))
    print("  PARTIAL                              : %s"
          % (", ".join(sorted(partial)) if partial else "none"))
    print("  SATURATED (detection holds at 5M)    : %s"
          % (", ".join(sorted(saturated)) if saturated else "none"))
    if nosig:
        print("  no signal at either depth            : %s"
              % ", ".join(sorted(nosig)))
    print("")
    print("  wrote %s" % cat_path)
    print("  wrote %s" % sample_path)
    print("  wrote %s  <- CONTAINS IDENTIFIERS, do not commit or email"
          % key_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
