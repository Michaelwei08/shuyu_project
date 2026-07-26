#!/usr/bin/env python3
"""Per-sample coverage profiles for named panel references (real virus vs pile-up).

WHAT IT COMPUTES
  For every requested reference (--ref, repeatable; default EBV type 1, EBV type 2,
  HHV-6B) and every BAM in every requested run directory (--run, repeatable):

    1. streams  samtools view -F 0x904 -q <min-mapq> <bam> <ref>
       and keeps only UNIQUE-BEST reads, i.e. records that carry AS:i: and have
       either no XS:i: or AS > XS (the same test as the awk idiom used elsewhere
       in this suite; it is done in python here so no awk is required),
    2. rebuilds per-base depth on that reference from POS + CIGAR
       (M/=/X/D consume reference and count as covered, N consumes reference
       without coverage, I/S/H/P do not consume reference),
    3. reports breadth (fraction of the reference at >=1x and >=5x), mean depth,
       median depth (over all reference positions, zeros included), mean depth over
       covered positions only, the coefficient of variation of depth across
       --bin-size (default 1000 bp) bins, and the fraction of all depth that sits
       in the single highest bin -- the pile-up detector,
    4. emits the binned coverage profile itself so a coverage track can be plotted
       downstream (this module draws no figures, so it is standard library only:
       no matplotlib, no numpy, no pysam),
    5. adds a ciHHV-6 heuristic call per (sample, reference):
         no_signal   fewer than --min-reads unique-best reads
         pileup_like max-bin depth fraction >= --pileup-max-bin-frac, or breadth
                     below --pileup-max-breadth: one hot spot, not a genome
         ciHHV6_candidate
                     breadth >= --cihhv6-min-breadth AND depth CV <= --cihhv6-max-cv
                     AND mean depth inside [--cihhv6-min-mean-depth,
                     --cihhv6-max-mean-depth] -- flat, genome-wide, about 1x, which
                     is what inherited chromosomally integrated HHV-6 looks like
         active_like breadth >= --active-min-breadth with mean depth above the
                     ciHHV-6 window: real, replicating, usually uneven
         low_level   anything else (a handful of scattered reads)
       Column ref_is_hhv6 marks the references where the ciHHV-6 reading is the
       biologically meaningful one; the same signature on EBV means "flat 1x", not
       inherited integration. An estimated host depth (from the run's headline
       deduplicated HUMAN counts, --read-length and --genome-size) and the
       virus/host depth ratio are reported as context: germline ciHHV-6 is expected
       near 0.5 x autosomal depth, one copy on one chromosome.

WHAT IT WRITES (all into --outdir, tab separated, pure ASCII, anon IDs only)
  <prefix>_summary.tsv     one row per run x reference x sample, all metrics + call
  <prefix>_bins.tsv        one row per run x reference x sample x bin, the profile
  <prefix>_sample_key.tsv  real -> anonymous mapping; CARRIES IDENTIFIERS, never
                           commit or email this file
  With the default --prefix that is refprofile_summary.tsv, refprofile_bins.tsv and
  refprofile_sample_key.tsv. Samples are anonymised to S01..Snn sorted by real name
  across all runs processed in one invocation. A short summary goes to stdout.

EXAMPLE
  python3 a5_reference_depth_profiles.py \
      --run /path/to/runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only \
      --run /path/to/runs/targeted_htlv_hg38_shuyu_masked_panel_refixed_primary_only \
      --ref hhv6b --ref ebv1 --ref ebv2 \
      --outdir /path/to/tmp/panel_report_20260725 --prefix refprofile

  Missing runs, missing BAMs, missing indexes, missing references and a missing
  samtools are all reported as "WARN: ... missing at ..., skipping" and the script
  still exits 0 with header-only tables if nothing could be computed.
  No network access. Date format YYYY-MM-DD.
"""
from __future__ import annotations

import argparse
import array
import csv
import datetime
import glob
import os
import re
import subprocess
import sys
import tempfile

DEFAULT_ROOT = "/path/to/runs"

DEFAULT_RUNS = [
    DEFAULT_ROOT + "/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only",
    DEFAULT_ROOT + "/targeted_htlv_hg38_shuyu_masked_panel_refixed_primary_only",
]

DEFAULT_REFMAP = (DEFAULT_ROOT + "/shuyu_masked_panel_hg38_herv_line1_refixed/ref/"
                  "hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv")

# Named panel references (reference_id -> short ASCII label).
KNOWN_REFS = [
    ("SHUYU_000096_NC_007605.1", "EBV_type1"),
    ("SHUYU_000101_NC_009334.1", "EBV_type2"),
    ("SHUYU_000001_NC_000898.1", "HHV6B"),
    ("SHUYU_000080_NC_006273.2", "CMV"),
    ("SHUYU_000054_NC_002076.2", "TTV1"),
]
LABEL_BY_ID = dict(KNOWN_REFS)

# Convenience aliases accepted by --ref.
REF_ALIASES = {
    "ebv1": "SHUYU_000096_NC_007605.1",
    "ebv_type1": "SHUYU_000096_NC_007605.1",
    "ebvtype1": "SHUYU_000096_NC_007605.1",
    "nc_007605.1": "SHUYU_000096_NC_007605.1",
    "ebv2": "SHUYU_000101_NC_009334.1",
    "ebv_type2": "SHUYU_000101_NC_009334.1",
    "ebvtype2": "SHUYU_000101_NC_009334.1",
    "nc_009334.1": "SHUYU_000101_NC_009334.1",
    "hhv6b": "SHUYU_000001_NC_000898.1",
    "hhv_6b": "SHUYU_000001_NC_000898.1",
    "hhv6": "SHUYU_000001_NC_000898.1",
    "nc_000898.1": "SHUYU_000001_NC_000898.1",
    "cmv": "SHUYU_000080_NC_006273.2",
    "nc_006273.2": "SHUYU_000080_NC_006273.2",
    "ttv1": "SHUYU_000054_NC_002076.2",
    "ttv": "SHUYU_000054_NC_002076.2",
    "nc_002076.2": "SHUYU_000054_NC_002076.2",
}

DEFAULT_REFS = ["SHUYU_000096_NC_007605.1",
                "SHUYU_000101_NC_009334.1",
                "SHUYU_000001_NC_000898.1"]

HHV6_HINTS = ("nc_000898", "nc_001664", "hhv-6", "hhv6",
              "herpesvirus 6", "betaherpesvirus 6", "roseolovirus")

SUMMARY_COLUMNS = [
    "run", "reference_id", "ref_label", "ref_len", "sample_anon", "group",
    "idxstats_mapped", "reads_pass_filters", "reads_uniqbest", "reads_ambiguous",
    "reads_no_as_tag", "aligned_bases", "covered_bases_1x", "breadth_1x",
    "breadth_5x", "mean_depth", "median_depth", "mean_depth_covered",
    "bin_size", "n_bins", "n_bins_nonzero", "depth_cv_bins",
    "max_bin_depth_frac", "max_bin_start", "max_bin_end",
    "human_reads_dedup", "est_host_depth", "depth_vs_host_ratio",
    "ref_is_hhv6", "call", "call_metrics",
]

BIN_COLUMNS = [
    "run", "reference_id", "ref_label", "sample_anon", "group",
    "bin_index", "bin_start", "bin_end", "bin_len",
    "depth_sum", "mean_depth", "covered_bases", "breadth_bin",
    "frac_total_depth",
]


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
def warn_missing(what, path):
    print("WARN: %s missing at %s, skipping" % (what, path))


def redact(path, sample, anon):
    """Replace a real sample name inside a path with its anonymous ID."""
    if not sample:
        return path
    return str(path).replace(sample, anon)


def fnum(value, nd=4):
    if value is None:
        return "NA"
    return "%.*f" % (nd, value)


def sanitize_label(text):
    out = re.sub(r"[^0-9A-Za-z]+", "_", str(text)).strip("_")
    return out[:48] if out else "NA"


def group_from_text(text):
    """Group label from a real sample name (falls back to the run name).

    Suite-wide rule, matched case-insensitively: "_HIV" -> HIV, "_HL" -> HL,
    "TCL"/"targeted_htlv" -> TCL. Returns None (not "NA") so the caller can
    chain a fallback.
    """
    if not text:
        return None
    up = str(text).upper()
    if "_HIV" in up:
        return "HIV"
    if "_HL" in up:
        return "HL"
    if "TARGETED_HTLV" in up or "TCL" in up:
        return "TCL"
    return None


def median_from_hist(hist, total_positions):
    """Median of a multiset given as {value: count}."""
    if total_positions <= 0:
        return None
    keys = sorted(hist)
    lo_target = (total_positions - 1) // 2
    hi_target = total_positions // 2
    seen = 0
    lo_val = hi_val = None
    for k in keys:
        seen += hist[k]
        if lo_val is None and seen > lo_target:
            lo_val = k
        if hi_val is None and seen > hi_target:
            hi_val = k
            break
    if lo_val is None or hi_val is None:
        return None
    return (lo_val + hi_val) / 2.0


def pop_cv(values):
    """Population coefficient of variation; None when the mean is 0."""
    n = len(values)
    if n == 0:
        return None
    mean = sum(values) / float(n)
    if mean <= 0:
        return None
    var = sum((v - mean) ** 2 for v in values) / float(n)
    return (var ** 0.5) / mean


# ----------------------------------------------------------------------------
# inputs
# ----------------------------------------------------------------------------
def load_refmap(path):
    """reference_id -> (category, description) from the reference map CSV."""
    info = {}
    if not path:
        return info
    if not os.path.exists(path):
        warn_missing("reference map", path)
        return info
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh):
                rid = (row.get("reference_id") or "").strip()
                if not rid:
                    continue
                desc = (row.get("description") or "").strip()
                if not desc:
                    desc = (row.get("original_header")
                            or row.get("original_reference_id") or "").strip()
                info[rid] = ((row.get("category") or "").strip(), desc)
    except Exception as exc:                      # unreadable map is not fatal
        print("WARN: could not parse reference map at %s (%s), continuing without it"
              % (path, exc.__class__.__name__))
    return info


def ref_label(reference_id, refmap):
    if reference_id in LABEL_BY_ID:
        return LABEL_BY_ID[reference_id]
    desc = refmap.get(reference_id, ("", ""))[1]
    return sanitize_label(desc) if desc else sanitize_label(reference_id)


def ref_is_hhv6(reference_id, refmap):
    blob = (reference_id + " " + refmap.get(reference_id, ("", ""))[1]
            + " " + ref_label(reference_id, refmap)).lower()
    return any(h in blob for h in HHV6_HINTS)


def resolve_ref(token):
    key = token.strip().lower().replace(" ", "_").replace("-", "_")
    return REF_ALIASES.get(key, token.strip())


def list_bams(run_dir, bam_glob):
    pattern = os.path.join(run_dir, bam_glob)
    return sorted(p for p in glob.glob(pattern) if os.path.isfile(p))


def sample_of_bam(bam_path):
    base = os.path.basename(bam_path)
    return base[:-4] if base.endswith(".bam") else base


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
                if len(p) < 3 or p[0] == "*":
                    continue
                try:
                    out[p[0]] = (int(p[1]), int(p[2]))
                except ValueError:
                    continue
    except Exception:
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


def human_counts(run_dir):
    """{sample: HUMAN reads} from the headline deduplicated category table."""
    res = os.path.join(run_dir, "results")
    hits = sorted(glob.glob(os.path.join(res, "*filtered_category_counts.tsv")))
    hits = [h for h in hits if "record_category_counts" not in os.path.basename(h)]
    if not hits:
        warn_missing("filtered_category_counts.tsv",
                     os.path.join(res, "*filtered_category_counts.tsv"))
        return {}
    path = hits[0]
    out = {}
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sample = (row.get("sample") or "").strip()
                if not sample:
                    continue
                raw = (row.get("HUMAN") or "").strip()
                try:
                    out[sample] = int(float(raw))
                except (TypeError, ValueError):
                    continue
    except Exception:
        warn_missing("readable filtered_category_counts.tsv", path)
        return {}
    return out


# ----------------------------------------------------------------------------
# core: unique-best depth on one reference of one BAM
# ----------------------------------------------------------------------------
def profile_bam_ref(samtools, bam_path, reference_id, ref_len, args, anon):
    """Return (stats_dict, bins_list) or None when samtools failed.

    anon is used only so that failures can be reported without printing a real
    sample name to stdout.
    """
    cmd = [samtools, "view", "-F", str(args.exclude_flags), "-q", str(args.min_mapq)]
    if args.samtools_threads > 1:
        cmd += ["-@", str(args.samtools_threads)]
    cmd += [bam_path, reference_id]

    delta = array.array("i", [0]) * (ref_len + 2)
    n_pass = n_uniq = n_amb = n_no_as = 0
    aligned_bases = 0

    # stderr goes to a temp file, not a pipe: stderr is only drained after the
    # stdout loop finishes, and a pipe that fills up first would deadlock.
    err_fh = tempfile.TemporaryFile(mode="w+", errors="replace")
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err_fh,
                                universal_newlines=True, errors="replace")
    except OSError:
        err_fh.close()
        warn_missing("samtools executable", samtools)
        return None

    for line in proc.stdout:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        n_pass += 1
        as_v = xs_v = None
        for tag in f[11:]:
            if tag.startswith("AS:i:"):
                as_v = tag[5:]
            elif tag.startswith("XS:i:"):
                xs_v = tag[5:]
        if as_v is None:
            n_no_as += 1
            continue
        try:
            as_i = int(as_v)
        except ValueError:
            n_no_as += 1
            continue
        if xs_v is not None:
            try:
                if as_i <= int(xs_v):
                    n_amb += 1
                    continue
            except ValueError:
                pass
        n_uniq += 1

        cigar = f[5]
        if cigar == "*":
            continue
        try:
            pos = int(f[3]) - 1
        except ValueError:
            continue
        num = 0
        for ch in cigar:
            if "0" <= ch <= "9":
                num = num * 10 + (ord(ch) - 48)
                continue
            if ch in "M=XD":
                start = pos if pos > 0 else 0
                end = pos + num
                if end > ref_len:
                    end = ref_len
                if end > start:
                    delta[start] += 1
                    delta[end] -= 1
                    aligned_bases += end - start
                pos += num
            elif ch == "N":
                pos += num
            # I, S, H, P consume query or nothing on the reference
            num = 0

    proc.stdout.close()
    rc = proc.wait()
    try:
        err_fh.seek(0)
        stderr = err_fh.read()
    except (IOError, OSError, ValueError):
        stderr = ""
    err_fh.close()
    if rc != 0:
        # samtools echoes the BAM path, which carries the real sample name
        stderr = redact(stderr, sample_of_bam(bam_path), anon)
        print("WARN: samtools view failed (rc=%d) for %s on %s: %s"
              % (rc, anon, reference_id,
                 " ".join(stderr.split())[:160] or "no stderr"))
        return None

    bin_size = args.bin_size
    n_bins = (ref_len + bin_size - 1) // bin_size
    bin_sum = [0] * n_bins
    bin_cov = [0] * n_bins

    hist = {}
    total_depth = 0
    covered = 0
    covered5 = 0
    if n_uniq > 0:
        depth = 0
        for i in range(ref_len):
            depth += delta[i]
            if depth:
                total_depth += depth
                covered += 1
                if depth >= args.breadth_depth2:
                    covered5 += 1
                b = i // bin_size
                bin_sum[b] += depth
                bin_cov[b] += 1
            hist[depth] = hist.get(depth, 0) + 1
    else:
        hist[0] = ref_len

    bin_rows = []
    bin_means = []
    max_bin_sum = 0
    max_bin_idx = 0
    for b in range(n_bins):
        start = b * bin_size
        blen = min(bin_size, ref_len - start)
        mean_b = bin_sum[b] / float(blen) if blen else 0.0
        bin_means.append(mean_b)
        if bin_sum[b] > max_bin_sum:
            max_bin_sum = bin_sum[b]
            max_bin_idx = b
        bin_rows.append({
            "bin_index": b,
            "bin_start": start + 1,
            "bin_end": start + blen,
            "bin_len": blen,
            "depth_sum": bin_sum[b],
            "mean_depth": mean_b,
            "covered_bases": bin_cov[b],
            "breadth_bin": (bin_cov[b] / float(blen)) if blen else None,
        })
    for row in bin_rows:
        row["frac_total_depth"] = (row["depth_sum"] / float(total_depth)
                                  if total_depth else None)

    stats = {
        "ref_len": ref_len,
        "reads_pass_filters": n_pass,
        "reads_uniqbest": n_uniq,
        "reads_ambiguous": n_amb,
        "reads_no_as_tag": n_no_as,
        "aligned_bases": aligned_bases,
        "covered_bases_1x": covered,
        "breadth_1x": covered / float(ref_len) if ref_len else None,
        "breadth_5x": covered5 / float(ref_len) if ref_len else None,
        "mean_depth": total_depth / float(ref_len) if ref_len else None,
        "median_depth": median_from_hist(hist, ref_len),
        "mean_depth_covered": (total_depth / float(covered)) if covered else None,
        "n_bins": n_bins,
        "n_bins_nonzero": sum(1 for v in bin_sum if v > 0),
        "depth_cv_bins": pop_cv(bin_means),
        "max_bin_depth_frac": (max_bin_sum / float(total_depth)
                               if total_depth else None),
        "max_bin_start": max_bin_idx * bin_size + 1 if total_depth else None,
        "max_bin_end": (min((max_bin_idx + 1) * bin_size, ref_len)
                        if total_depth else None),
    }
    return stats, bin_rows


# ----------------------------------------------------------------------------
# ciHHV-6 / pile-up heuristic
# ----------------------------------------------------------------------------
def classify(stats, args):
    reads = stats["reads_uniqbest"]
    breadth = stats["breadth_1x"]
    cv = stats["depth_cv_bins"]
    mean_depth = stats["mean_depth"]
    frac = stats["max_bin_depth_frac"]

    metrics = "reads=%d;breadth=%s;cv=%s;mean_depth=%s;max_bin_frac=%s" % (
        reads, fnum(breadth, 4), fnum(cv, 3), fnum(mean_depth, 4), fnum(frac, 3))

    if reads < args.min_reads:
        return "no_signal", metrics
    if (frac is not None and frac >= args.pileup_max_bin_frac) or \
       (breadth is not None and breadth < args.pileup_max_breadth):
        return "pileup_like", metrics
    if (breadth is not None and breadth >= args.cihhv6_min_breadth
            and cv is not None and cv <= args.cihhv6_max_cv
            and mean_depth is not None
            and args.cihhv6_min_mean_depth <= mean_depth <= args.cihhv6_max_mean_depth):
        return "ciHHV6_candidate", metrics
    if (breadth is not None and breadth >= args.active_min_breadth
            and mean_depth is not None and mean_depth > args.cihhv6_max_mean_depth):
        return "active_like", metrics
    return "low_level", metrics


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="Unique-best coverage profiles on named panel references; "
                    "pile-up detection and a ciHHV-6 heuristic. Standard library "
                    "only, no figures.")
    p.add_argument("--run", action="append", default=None, metavar="DIR",
                   help="run directory containing bam/ and results/ (repeatable). "
                        "Default: the two panel runs (WGS then targeted).")
    p.add_argument("--ref", action="append", default=None, metavar="REF",
                   help="reference_id or alias (ebv1, ebv2, hhv6b, cmv, ttv1); "
                        "repeatable. Default: ebv1, ebv2, hhv6b.")
    p.add_argument("--refmap", default=DEFAULT_REFMAP,
                   help="panel reference map CSV (labels/descriptions only)")
    p.add_argument("--outdir", default="/path/to/tmp/panel_report_20260725",
                   help="output directory (created if absent)")
    p.add_argument("--prefix", default="refprofile", help="output filename prefix")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--samtools-threads", type=int, default=1,
                   help="samtools -@ value")
    p.add_argument("--bam-glob", default="bam/*.bam",
                   help="glob for BAMs relative to a run directory")
    p.add_argument("--bin-size", type=int, default=1000, help="bin size in bp")
    p.add_argument("--min-mapq", type=int, default=40, help="samtools -q value")
    p.add_argument("--exclude-flags", default="0x904", help="samtools -F value")
    p.add_argument("--breadth-depth2", type=int, default=5,
                   help="second breadth threshold reported as breadth_5x")
    p.add_argument("--min-reads", type=int, default=3,
                   help="unique-best reads below which the call is no_signal")
    p.add_argument("--bins-min-reads", type=int, default=1,
                   help="emit the binned profile only for samples with at least "
                        "this many unique-best reads (0 = every sample)")
    p.add_argument("--limit", type=int, default=0,
                   help="process only the first N BAMs per run (debugging)")
    p.add_argument("--read-length", type=int, default=150,
                   help="assumed read length for the host depth estimate")
    p.add_argument("--genome-size", type=float, default=3.1e9,
                   help="assumed haploid genome size for the host depth estimate")
    p.add_argument("--cihhv6-min-breadth", type=float, default=0.90)
    p.add_argument("--cihhv6-max-cv", type=float, default=0.50)
    p.add_argument("--cihhv6-min-mean-depth", type=float, default=0.30)
    p.add_argument("--cihhv6-max-mean-depth", type=float, default=3.00)
    p.add_argument("--pileup-max-bin-frac", type=float, default=0.50,
                   help="depth fraction in one bin at or above which a profile is "
                        "called pileup_like")
    p.add_argument("--pileup-max-breadth", type=float, default=0.05,
                   help="breadth below which a profile is called pileup_like")
    p.add_argument("--active-min-breadth", type=float, default=0.30,
                   help="breadth needed for an active_like call")
    return p


def check_samtools(samtools):
    try:
        proc = subprocess.run([samtools, "--version"], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        warn_missing("samtools executable", samtools)
        return None
    first = (proc.stdout or proc.stderr or "").splitlines()
    return first[0].strip() if first else "samtools (version unknown)"


def write_headers_only(paths):
    with open(paths["summary"], "w", encoding="ascii", newline="") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS) + "\n")
    with open(paths["bins"], "w", encoding="ascii", newline="") as fh:
        fh.write("\t".join(BIN_COLUMNS) + "\n")


def main(argv=None):
    args = build_parser().parse_args(argv)
    today = datetime.date.today().strftime("%Y-%m-%d")

    runs = args.run if args.run else list(DEFAULT_RUNS)
    refs = [resolve_ref(r) for r in (args.ref if args.ref else DEFAULT_REFS)]
    seen = set()
    refs = [r for r in refs if not (r in seen or seen.add(r))]
    if args.bin_size < 1:
        args.bin_size = 1000

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            print("WARN: output directory %s is not usable, nothing written"
                  % args.outdir)
            return 0
    paths = {
        "summary": os.path.join(args.outdir, args.prefix + "_summary.tsv"),
        "bins": os.path.join(args.outdir, args.prefix + "_bins.tsv"),
        "key": os.path.join(args.outdir, args.prefix + "_sample_key.tsv"),
    }

    print("a5_reference_depth_profiles.py  %s" % today)
    print("outdir   %s" % args.outdir)
    print("refs     %s" % ", ".join(refs))
    print("filter   samtools view -F %s -q %d <bam> <ref>, then unique-best "
          "(AS present and (no XS or AS>XS))" % (args.exclude_flags, args.min_mapq))
    print("bins     %d bp" % args.bin_size)
    print()

    ver = check_samtools(args.samtools)
    if ver is None:
        print("Nothing can be computed without samtools; writing empty tables.")
        write_headers_only(paths)
        with open(paths["key"], "w", encoding="ascii", newline="") as fh:
            fh.write("# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL\n")
            fh.write("# generated %s by a5_reference_depth_profiles.py\n" % today)
            fh.write("sample_anon\tsample_real\tgroup\truns\n")
        return 0
    print("samtools %s" % ver)

    # ---- pass 1: discover runs / BAMs, build the anonymisation map ----------
    run_bams = []          # [(run_dir, run_base, [bam_path, ...])]
    for run_dir in runs:
        run_dir = run_dir.rstrip("/\\")
        base = os.path.basename(run_dir) or run_dir
        if not os.path.isdir(run_dir):
            warn_missing("run directory", run_dir)
            continue
        bams = list_bams(run_dir, args.bam_glob)
        if not bams:
            warn_missing("BAM files", os.path.join(run_dir, args.bam_glob))
            continue
        if args.limit > 0:
            bams = bams[:args.limit]
        run_bams.append((run_dir, base, bams))

    real_names = set()
    for _run_dir, _base, bams in run_bams:
        for b in bams:
            real_names.add(sample_of_bam(b))
    ordered = sorted(real_names)
    width = max(2, len(str(len(ordered))))
    anon_of = {}
    for i, name in enumerate(ordered, start=1):
        anon_of[name] = "S" + str(i).zfill(width)

    group_of = {}
    runs_of = {}
    for _run_dir, base, bams in run_bams:
        for b in bams:
            s = sample_of_bam(b)
            runs_of.setdefault(s, [])
            if base not in runs_of[s]:
                runs_of[s].append(base)
            if group_of.get(s) in (None, "NA"):
                group_of[s] = group_from_text(s) or group_from_text(base) or "NA"

    with open(paths["key"], "w", encoding="ascii", newline="") as fh:
        fh.write("# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL\n")
        fh.write("# generated %s by a5_reference_depth_profiles.py\n" % today)
        fh.write("sample_anon\tsample_real\tgroup\truns\n")
        for name in ordered:
            fh.write("%s\t%s\t%s\t%s\n" % (anon_of[name], name,
                                           group_of.get(name, "NA"),
                                           ";".join(runs_of.get(name, []))))
    print("samples  %d unique real names -> %s..%s"
          % (len(ordered), anon_of[ordered[0]] if ordered else "-",
             anon_of[ordered[-1]] if ordered else "-"))
    print("key      %s (identifiers; do not commit or email)" % paths["key"])
    print()

    # Drop BAMs without an index once, so the warning is not repeated per
    # reference. Paths are redacted to the anonymous ID before printing.
    indexed = []
    for run_dir, base, bams in run_bams:
        keep = []
        for bam_path in bams:
            sample = sample_of_bam(bam_path)
            if os.path.exists(bam_path + ".bai") or \
                    os.path.exists(bam_path[:-4] + ".bai") or \
                    os.path.exists(bam_path + ".csi"):
                keep.append(bam_path)
            else:
                warn_missing("BAM index for " + anon_of[sample],
                             redact(bam_path + ".bai", sample, anon_of[sample]))
        if keep:
            indexed.append((run_dir, base, keep))
        else:
            warn_missing("indexed BAMs in run " + base,
                         os.path.join(run_dir, args.bam_glob))
    run_bams = indexed

    if not run_bams:
        print("No usable run directory; writing empty tables.")
        write_headers_only(paths)
        return 0

    refmap = load_refmap(args.refmap)

    # ---- pass 2: compute --------------------------------------------------
    summary_rows = []
    bins_fh = open(paths["bins"], "w", encoding="ascii", newline="")
    bins_fh.write("\t".join(BIN_COLUMNS) + "\n")

    for run_dir, base, bams in run_bams:
        humans = human_counts(run_dir)
        len_cache = {}
        for reference_id in refs:
            label = ref_label(reference_id, refmap)
            is_hhv6 = "1" if ref_is_hhv6(reference_id, refmap) else "0"
            done = 0
            warned_absent = False
            for bam_path in bams:
                sample = sample_of_bam(bam_path)
                anon = anon_of[sample]
                group = group_from_text(sample) or group_from_text(base) or "NA"
                if sample not in len_cache:
                    lens = idxstats_lengths(run_dir, sample)
                    if not lens:
                        lens = header_lengths(args.samtools, bam_path)
                    len_cache[sample] = lens
                lens = len_cache[sample]
                ref_len, idx_mapped = lens.get(reference_id, (None, None))
                if not ref_len or ref_len < 1:
                    if not warned_absent:
                        warned_absent = True
                        warn_missing("reference %s (%s) in run %s"
                                     % (reference_id, label, base),
                                     os.path.join(run_dir, "bam", "*.bam"))
                    continue

                result = profile_bam_ref(args.samtools, bam_path, reference_id,
                                         ref_len, args, anon)
                if result is None:
                    continue
                stats, bin_rows = result
                done += 1

                hr = humans.get(sample)
                est_host = None
                ratio = None
                if hr is not None and args.genome_size > 0:
                    est_host = hr * float(args.read_length) / float(args.genome_size)
                    if est_host > 0 and stats["mean_depth"] is not None:
                        ratio = stats["mean_depth"] / est_host
                call, metrics = classify(stats, args)

                row = {
                    "run": base,
                    "reference_id": reference_id,
                    "ref_label": label,
                    "ref_len": ref_len,
                    "sample_anon": anon,
                    "group": group,
                    "idxstats_mapped": "NA" if idx_mapped is None else idx_mapped,
                    "reads_pass_filters": stats["reads_pass_filters"],
                    "reads_uniqbest": stats["reads_uniqbest"],
                    "reads_ambiguous": stats["reads_ambiguous"],
                    "reads_no_as_tag": stats["reads_no_as_tag"],
                    "aligned_bases": stats["aligned_bases"],
                    "covered_bases_1x": stats["covered_bases_1x"],
                    "breadth_1x": fnum(stats["breadth_1x"], 6),
                    "breadth_5x": fnum(stats["breadth_5x"], 6),
                    "mean_depth": fnum(stats["mean_depth"], 6),
                    "median_depth": fnum(stats["median_depth"], 2),
                    "mean_depth_covered": fnum(stats["mean_depth_covered"], 4),
                    "bin_size": args.bin_size,
                    "n_bins": stats["n_bins"],
                    "n_bins_nonzero": stats["n_bins_nonzero"],
                    "depth_cv_bins": fnum(stats["depth_cv_bins"], 4),
                    "max_bin_depth_frac": fnum(stats["max_bin_depth_frac"], 4),
                    "max_bin_start": ("NA" if stats["max_bin_start"] is None
                                      else stats["max_bin_start"]),
                    "max_bin_end": ("NA" if stats["max_bin_end"] is None
                                    else stats["max_bin_end"]),
                    "human_reads_dedup": "NA" if hr is None else hr,
                    "est_host_depth": fnum(est_host, 4),
                    "depth_vs_host_ratio": fnum(ratio, 4),
                    "ref_is_hhv6": is_hhv6,
                    "call": call,
                    "call_metrics": metrics,
                }
                summary_rows.append(row)

                if stats["reads_uniqbest"] >= args.bins_min_reads:
                    for br in bin_rows:
                        bins_fh.write("\t".join([
                            base, reference_id, label, anon, group,
                            str(br["bin_index"]), str(br["bin_start"]),
                            str(br["bin_end"]), str(br["bin_len"]),
                            str(br["depth_sum"]), fnum(br["mean_depth"], 4),
                            str(br["covered_bases"]), fnum(br["breadth_bin"], 6),
                            fnum(br["frac_total_depth"], 6),
                        ]) + "\n")
            print("  %-46s %-11s %3d/%3d samples profiled"
                  % (base[:46], label, done, len(bams)))
    bins_fh.close()

    with open(paths["summary"], "w", encoding="ascii", newline="") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for row in summary_rows:
            fh.write("\t".join(str(row[c]) for c in SUMMARY_COLUMNS) + "\n")

    # ---- stdout summary ---------------------------------------------------
    print()
    print("wrote %s (%d rows)" % (paths["summary"], len(summary_rows)))
    print("wrote %s" % paths["bins"])
    if not summary_rows:
        print("No (sample, reference) pair could be profiled.")
        return 0

    calls = ["no_signal", "low_level", "pileup_like", "ciHHV6_candidate",
             "active_like"]
    print()
    print("-- per run x reference (medians over samples with >= %d unique-best reads) --"
          % args.min_reads)
    hdr = "%-30s %-10s %5s %6s %8s %8s %7s  " % (
        "run", "ref", "n", "signal", "med_brd", "med_dpth", "med_cv")
    print(hdr + " ".join("%*s" % (max(5, len(c)), c) for c in calls))
    keys = []
    for row in summary_rows:
        k = (row["run"], row["ref_label"])
        if k not in keys:
            keys.append(k)
    for run_name, label in keys:
        grp = [r for r in summary_rows
               if r["run"] == run_name and r["ref_label"] == label]
        sig = [r for r in grp if int(r["reads_uniqbest"]) >= args.min_reads]

        def med(field):
            vals = sorted(float(r[field]) for r in sig if r[field] != "NA")
            if not vals:
                return None
            n = len(vals)
            return (vals[(n - 1) // 2] + vals[n // 2]) / 2.0

        counts = " ".join("%*d" % (max(5, len(c)),
                                   sum(1 for r in grp if r["call"] == c))
                          for c in calls)
        print("%-30s %-10s %5d %6d %8s %8s %7s  %s"
              % (run_name[:30], label[:10], len(grp), len(sig),
                 fnum(med("breadth_1x"), 3), fnum(med("mean_depth"), 3),
                 fnum(med("depth_cv_bins"), 2), counts))

    cand = [r for r in summary_rows if r["call"] == "ciHHV6_candidate"]
    print()
    print("-- ciHHV-6 heuristic candidates (breadth >= %.2f, CV <= %.2f, mean depth "
          "%.2f-%.2f x) --" % (args.cihhv6_min_breadth, args.cihhv6_max_cv,
                               args.cihhv6_min_mean_depth, args.cihhv6_max_mean_depth))
    if not cand:
        print("   none")
    else:
        print("   %-6s %-6s %-10s %-26s %7s %7s %6s %7s %8s"
              % ("sample", "group", "ref", "run", "breadth", "depth", "cv",
                 "maxbin", "vs_host"))
        for r in sorted(cand, key=lambda x: (x["ref_label"], x["sample_anon"])):
            flag = "" if r["ref_is_hhv6"] == "1" else "   (not an HHV-6 reference)"
            print("   %-6s %-6s %-10s %-26s %7s %7s %6s %7s %8s%s"
                  % (r["sample_anon"], r["group"], r["ref_label"][:10],
                     r["run"][:26], r["breadth_1x"][:7], r["mean_depth"][:7],
                     r["depth_cv_bins"][:6], r["max_bin_depth_frac"][:7],
                     r["depth_vs_host_ratio"][:8], flag))
    pile = [r for r in summary_rows if r["call"] == "pileup_like"]
    act = [r for r in summary_rows if r["call"] == "active_like"]
    print()
    print("totals: %d rows, %d pileup_like, %d active_like, %d ciHHV6_candidate"
          % (len(summary_rows), len(pile), len(act), len(cand)))
    print("Reminder: only %s carries real sample names."
          % os.path.basename(paths["key"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
