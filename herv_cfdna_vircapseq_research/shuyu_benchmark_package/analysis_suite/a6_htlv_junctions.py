#!/usr/bin/env python3
"""a6_htlv_junctions.py - CANDIDATE HTLV-1 host-virus integration junctions.

WHAT THIS COMPUTES
------------------
Integration junctions are the strongest read-level evidence that an HTLV-1 call
is a real infection (a provirus physically joined to host DNA) rather than
index hopping, library carry-over, or reference contamination. The number of
distinct junction clusters in a sample is also a crude clonality proxy: one
dominant cluster looks monoclonal, many low-support clusters look polyclonal.

Two independent signals are read from the BAMs of a run, on every reference
whose reference-map category is the virus category (default HTLV1):

  (1) DISCORDANT PAIRS. A read aligned to the HTLV-1 reference whose mate
      (RNEXT/PNEXT) is on a different reference whose category is HUMAN. The
      junction anchor on the virus is taken as the read's alignment END when
      the read is forward (the mate lies downstream) and as the read's POS when
      the read is reverse, i.e. the breakpoint lies just beyond the read in the
      direction it points. Mates on non-human, non-viral-of-interest references
      are counted separately (n_discordant_nonhuman) as an artefact gauge.

  (2) SOFT-CLIPPED READS. A read on the HTLV-1 reference with a leading or
      trailing S operation of >= --min-clip bp. The clip position (POS for a
      left clip, alignment end for a right clip) is the putative provirus/host
      boundary. Clips whose sequence is >= --max-homopolymer-frac one base are
      dropped as adapter/poly-A artefacts (SEQ "*" records are kept).

Reads are primary, non-duplicate, MAPQ >= --mapq, and (unless
--allow-ambiguous) unique-best: AS present and (XS absent or AS > XS), the same
filter used elsewhere in this suite. Positions from both signals are pooled per
(sample, virus reference) and clustered greedily: a new cluster starts when the
next position is more than --cluster-bp away from the previous one. Human mate
coordinates inside a cluster are themselves clustered (--mate-cluster-bp) to
report the most-supported host site.

THESE ARE CANDIDATES, NOT CALLS. Every cluster in the output needs manual /
IGV review before it is believed. Targeted capture makes chimeric-artefact
junctions genuinely likely: capture probes co-hybridise virus and host
fragments, PCR chimeras form on shared/repeat sequence, and soft clips pile up
at reference ends, low-complexity tracts, and provirus LTR boundaries. Clusters
whose host mate lands in a repeat, whose support is one read, or that recur at
the same virus coordinate across many unrelated samples are suspect.

WHAT IT WRITES (tab-separated, into --outdir)
--------------------------------------------
  <prefix>_candidates.tsv   one row per junction cluster (anon sample IDs only)
  <prefix>_per_sample.tsv   one row per sample: junction counts, clonality proxy
  <prefix>_sample_key.tsv   real -> anon mapping; CONTAINS IDENTIFIERS

With the default --prefix htlv_junction these are htlv_junction_candidates.tsv,
htlv_junction_per_sample.tsv and htlv_junction_sample_key.tsv. No real sample
name is written to any file other than the key, and none is printed to stdout.

Standard library only; no figures; no network. samtools is called through
subprocess (--samtools). Missing runs, BAMs, indexes or reference maps produce a
WARN line and are skipped; the script still writes headed tables and exits 0.

EXAMPLE
-------
  python3 a6_htlv_junctions.py \
      --run-dir /path/to/runs/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord \
      --outdir  /path/to/runs/reports/a6_htlv_junctions_2026-07-26 \
      --min-clip 20 --cluster-bp 50 --min-support 2

Written 2026-07-26.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile

SHUYU_ROOT = "/path/to/runs"

# Cluster paths are POSIX; keep them literal so --help reads the same anywhere.
DEFAULT_RUNS = [
    SHUYU_ROOT + "/targeted_htlv_hg38_refseq_mapq_human60_viral40_coord",
]
DEFAULT_BASE_REFMAP = (
    SHUYU_ROOT + "/retro_reference_hg38_refseq/ref/"
    "hg38_plus_retro.refseq.reference_map.csv")
DEFAULT_PANEL_REFMAP = (
    SHUYU_ROOT + "/shuyu_masked_panel_hg38_herv_line1_refixed/ref/"
    "hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv")
DEFAULT_OUTDIR = SHUYU_ROOT + "/reports/a6_htlv_junctions"

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = frozenset("MDN=X")
HUMAN_FALLBACK_RE = re.compile(r"^(chr)?([0-9]{1,2}|X|Y|M|MT)$", re.IGNORECASE)

KEY_WARNING = "# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL"


# ---------------------------------------------------------------- utilities

# Suffixes some pipelines add to the BAM basename but not to the sibling
# idxstats / count-table sample name. Stripped when deriving a sample id.
BAM_NAME_SUFFIXES = (".retrovirus", ".retro", ".markdup", ".dedup",
                     ".sorted", ".filtered")


def sample_from_bam(bam_path):
    """Sample id from a BAM path, minus any pipeline-added suffix."""
    base = os.path.basename(bam_path)
    if base.endswith(".bam"):
        base = base[:-4]
    for _suf in BAM_NAME_SUFFIXES:
        if base.endswith(_suf):
            base = base[: -len(_suf)]
            break
    return base


def warn(what, path):
    """The suite's standard missing-input line."""
    print("WARN: %s missing at %s, skipping" % (what, path))


def group_label(sample_name, run_basename=""):
    """HIV / HL / TCL / NA from the real sample name (run name as fallback).

    Suite-wide rule, matched case-insensitively so every module agrees:
    "_HIV" -> HIV, "_HL" -> HL, "TCL"/"targeted_htlv" -> TCL, else NA.
    """
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
    hay = up + " " + (run_basename or "").upper()
    if "TARGETED_HTLV" in hay or "TCL" in hay:
        return "TCL"
    return "NA"


def anon_ids(real_names):
    """Sorted real names -> S01..Snn (width grows past 99 samples)."""
    uniq = sorted(set(real_names))
    width = max(2, len(str(len(uniq))))
    return dict((n, "S" + str(i + 1).zfill(width)) for i, n in enumerate(uniq))


def masked(text, sample, anon_id):
    """Text with the real sample name swapped for its anonymous id.

    Used on every path and samtools message that reaches stdout, so a run log
    can be pasted into an email without carrying an identifier.
    """
    out = str(text)
    if sample and anon_id:
        out = out.replace(sample, anon_id)
    return out


def resolve_refmap(run_dir, explicit):
    """Explicit path, then the map the run recorded, then a name heuristic."""
    if explicit:
        return explicit if os.path.exists(explicit) else None
    man = os.path.join(run_dir, "results", "run_manifest.tsv")
    if os.path.exists(man):
        try:
            with open(man, newline="", encoding="utf-8", errors="replace") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        except Exception:
            rows = []
        for row in rows[:1]:
            for key, val in row.items():
                val = (val or "").strip()
                if not val or not val.endswith(".csv"):
                    continue
                if "reference_map" in (key or "").lower() or "reference_map" in val:
                    if os.path.exists(val):
                        return val
    base = os.path.basename(os.path.normpath(run_dir)).lower()
    guess = DEFAULT_PANEL_REFMAP if "panel" in base else DEFAULT_BASE_REFMAP
    return guess if os.path.exists(guess) else None


def load_refmap(path):
    """reference_id -> category."""
    cats = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            rid = (row.get("reference_id") or "").strip()
            if rid:
                cats[rid] = (row.get("category") or "").strip()
    return cats


def bam_index_present(bam):
    for cand in (bam + ".bai", bam + ".csi",
                 os.path.splitext(bam)[0] + ".bai",
                 os.path.splitext(bam)[0] + ".csi"):
        if os.path.exists(cand):
            return True
    return False


def samtools_ok(samtools):
    try:
        proc = subprocess.run([samtools, "--version"], stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    except OSError:
        return False
    return proc.returncode == 0


def bam_references(samtools, bam):
    """@SQ names present in the BAM header (order preserved)."""
    try:
        proc = subprocess.run([samtools, "view", "-H", bam],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        return None
    if proc.returncode != 0:
        return None
    names = []
    for line in proc.stdout.decode("utf-8", "replace").splitlines():
        if not line.startswith("@SQ"):
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                names.append(field[3:])
                break
    return names


def idxstats_virus_reads(run_dir, sample, virus_refs):
    """Raw mapped reads on the virus refs from results/<sample>.idxstats.tsv.

    Returns None when the file is absent. The raw count is a superset of what
    the filters below keep, so 0 is a safe reason to skip a BAM entirely.
    """
    path = os.path.join(run_dir, "results", sample + ".idxstats.tsv")
    if not os.path.exists(path):
        return None
    total = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3 or parts[0] == "*":
                    continue
                if parts[0] in virus_refs:
                    try:
                        total += int(parts[2])
                    except ValueError:
                        continue
    except OSError:
        return None
    return total


# ---------------------------------------------------------------- SAM parsing

def cigar_ref_span(cigar):
    """(reference bases consumed, leading soft clip, trailing soft clip)."""
    ops = CIGAR_RE.findall(cigar)
    if not ops:
        return 0, 0, 0
    span = 0
    for length, op in ops:
        if op in REF_CONSUMING:
            span += int(length)
    lead = 0
    for length, op in ops:
        if op == "H":
            continue
        if op == "S":
            lead = int(length)
        break
    trail = 0
    for length, op in reversed(ops):
        if op == "H":
            continue
        if op == "S":
            trail = int(length)
        break
    return span, lead, trail


def homopolymer_frac(seq):
    if not seq:
        return 0.0
    counts = {}
    for base in seq.upper():
        counts[base] = counts.get(base, 0) + 1
    return max(counts.values()) / float(len(seq))


def scan_virus_reference(samtools, bam, ref, args, human_refs, virus_refs):
    """Collect junction events for one BAM x one virus reference.

    Returns (events, n_reads_kept, n_mates_on_other_refs, error_or_None). Each
    event is a dict with pos, kind ("discordant" | "clip"), side, mate_ref,
    mate_pos, clip_len.
    """
    cmd = [samtools, "view", "-F", str(args.exclude_flags),
           "-q", str(args.mapq)]
    if args.threads > 0:
        cmd += ["-@", str(args.threads)]
    cmd += [bam, ref]
    # stderr goes to a temp file, not a pipe: it is only drained after the
    # stdout loop finishes, and a pipe that filled up first would deadlock.
    err_fh = tempfile.TemporaryFile(mode="w+b")
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err_fh)
    except OSError as exc:
        err_fh.close()
        return [], 0, 0, str(exc)

    events = []
    kept = 0
    nonhuman = 0
    assert proc.stdout is not None
    for raw in proc.stdout:
        line = raw.decode("utf-8", "replace")
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        try:
            flag = int(f[1])
            pos = int(f[3])
            pnext = int(f[7])
        except ValueError:
            continue
        rname, cigar, rnext, seq = f[2], f[5], f[6], f[9]
        if rname != ref or cigar == "*":
            continue

        # unique-best filter: AS present and better than any XS
        if not args.allow_ambiguous:
            a_score = x_score = None
            for tag in f[11:]:
                if tag.startswith("AS:i:"):
                    a_score = tag[5:].strip()
                elif tag.startswith("XS:i:"):
                    x_score = tag[5:].strip()
            if a_score is None:
                continue
            try:
                if x_score is not None and float(a_score) <= float(x_score):
                    continue
            except ValueError:
                continue

        kept += 1
        span, lead_clip, trail_clip = cigar_ref_span(cigar)
        aln_end = pos + max(span, 1) - 1
        reverse = bool(flag & 0x10)

        # (1) discordant pair with a human mate
        if (flag & 0x1) and not (flag & 0x8) and rnext not in ("=", "*", ""):
            mate_mapq = None
            for tag in f[11:]:
                if tag.startswith("MQ:i:"):
                    try:
                        mate_mapq = int(tag[5:].strip())
                    except ValueError:
                        mate_mapq = None
                    break
            mate_ok = mate_mapq is None or mate_mapq >= args.min_mate_mapq
            if rnext in human_refs:
                if mate_ok:
                    events.append({
                        "pos": pos if reverse else aln_end,
                        "kind": "discordant", "side": "R" if not reverse else "L",
                        "mate_ref": rnext, "mate_pos": pnext, "clip_len": 0,
                    })
            elif rnext not in virus_refs:
                nonhuman += 1

        # (2) soft clip marking a provirus/host boundary
        for clip_len, side, clip_pos, clip_seq in (
                (lead_clip, "L", pos, seq[:lead_clip] if seq != "*" else ""),
                (trail_clip, "R", aln_end,
                 seq[len(seq) - trail_clip:] if seq != "*" and trail_clip else "")):
            if clip_len < args.min_clip:
                continue
            if clip_seq and homopolymer_frac(clip_seq) >= args.max_homopolymer_frac:
                continue
            events.append({"pos": clip_pos, "kind": "clip", "side": side,
                           "mate_ref": "", "mate_pos": 0, "clip_len": clip_len})

    proc.stdout.close()
    rc = proc.wait()
    try:
        err_fh.seek(0)
        stderr = err_fh.read().decode("utf-8", "replace")
    except (IOError, OSError, ValueError):
        stderr = ""
    err_fh.close()
    if rc != 0:
        msg = (stderr.strip().splitlines() or ["samtools rc=%d" % rc])[-1]
        return events, kept, nonhuman, msg
    return events, kept, nonhuman, None


# ---------------------------------------------------------------- clustering

def cluster_positions(events, window):
    """Greedy single-linkage clustering of events on one reference."""
    clusters = []
    for ev in sorted(events, key=lambda e: e["pos"]):
        if clusters and ev["pos"] - clusters[-1][-1]["pos"] <= window:
            clusters[-1].append(ev)
        else:
            clusters.append([ev])
    return clusters


def cluster_host_sites(events, window):
    """Group human mate coordinates within a junction cluster."""
    mates = [(e["mate_ref"], e["mate_pos"]) for e in events
             if e["kind"] == "discordant" and e["mate_ref"]]
    sites = []
    for chrom in sorted(set(m[0] for m in mates)):
        run = []
        for _c, pos in sorted((m for m in mates if m[0] == chrom),
                              key=lambda m: m[1]):
            if run and pos - run[-1] <= window:
                run.append(pos)
            else:
                if run:
                    sites.append((chrom, run[0], run[-1], len(run)))
                run = [pos]
        if run:
            sites.append((chrom, run[0], run[-1], len(run)))
    sites.sort(key=lambda s: (-s[3], s[0], s[1]))
    return sites


def median_int(values):
    vals = sorted(values)
    n = len(vals)
    if not n:
        return 0
    if n % 2:
        return vals[n // 2]
    return int((vals[n // 2 - 1] + vals[n // 2]) / 2)


def summarise_cluster(events, mate_window):
    n_disc = sum(1 for e in events if e["kind"] == "discordant")
    n_clip = sum(1 for e in events if e["kind"] == "clip")
    n_left = sum(1 for e in events if e["kind"] == "clip" and e["side"] == "L")
    n_right = sum(1 for e in events if e["kind"] == "clip" and e["side"] == "R")
    clip_lens = [e["clip_len"] for e in events if e["kind"] == "clip"]
    sites = cluster_host_sites(events, mate_window)
    top = sites[0] if sites else ("NA", 0, 0, 0)
    site_str = ";".join("%s:%d-%d(n=%d)" % s for s in sites[:3]) or "NA"
    if n_disc and n_clip:
        support_class = "both"
    elif n_disc:
        support_class = "discordant_only"
    else:
        support_class = "clip_only"
    return {
        "junction_pos": median_int([e["pos"] for e in events]),
        "cluster_start": min(e["pos"] for e in events),
        "cluster_end": max(e["pos"] for e in events),
        "n_support_total": len(events),
        "n_discordant": n_disc,
        "n_clipped": n_clip,
        "n_clip_left": n_left,
        "n_clip_right": n_right,
        "max_clip_len": max(clip_lens) if clip_lens else 0,
        "host_sites_top3": site_str,
        "n_host_sites": len(sites),
        "best_host_chrom": top[0],
        "best_host_pos": median_int([top[1], top[2]]) if sites else 0,
        "best_host_support": top[3],
        "support_class": support_class,
    }


# ---------------------------------------------------------------- per run

def list_run_bams(run_dir, args):
    """(bam_dir, bam paths) for one run, in the order process_run will use them.

    Called once up front by main() so the anonymisation map exists before any
    per-sample warning is printed, and once inside process_run.
    """
    bam_dir = os.path.join(run_dir, args.bam_subdir)
    if not os.path.isdir(bam_dir):
        return bam_dir, []
    bams = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))
    if args.limit > 0:
        bams = bams[:args.limit]
    return bam_dir, bams


def process_run(run_dir, args, samtools, anon):
    """Returns (sample_rows, real_names) for one run directory.

    Each sample row carries its real name plus a "clusters" list; the caller
    anonymises before anything is written. anon is the real -> S01..Snn map
    built before this call, so nothing printed here is an identifier.
    """
    run_name = os.path.basename(os.path.normpath(run_dir))
    if not os.path.isdir(run_dir):
        warn("run directory", run_dir)
        return [], []

    refmap_path = resolve_refmap(run_dir, args.refmap)
    if not refmap_path:
        warn("reference map", os.path.join(run_dir, "results", "run_manifest.tsv"))
        cats = {}
    else:
        try:
            cats = load_refmap(refmap_path)
        except OSError:
            warn("reference map", refmap_path)
            cats = {}

    virus_cat = args.virus_category.upper()
    virus_refs = set(r for r, c in cats.items() if c.upper() == virus_cat)
    human_refs = set(r for r, c in cats.items() if c.upper() == args.human_category.upper())
    fallback = ""
    if not virus_refs or not human_refs:
        fallback = " (name-based fallback in use)"

    bam_dir, bams = list_run_bams(run_dir, args)
    if not os.path.isdir(bam_dir):
        warn("bam directory", bam_dir)
        return [], []
    if not bams:
        warn("BAM files", os.path.join(bam_dir, "*.bam"))
        return [], []

    print("RUN     %s" % run_dir)
    print("REFMAP  %s%s" % (refmap_path or "NONE", fallback))
    print("SAMPLES %d BAM files | %s refs in map: %d | HUMAN refs: %d"
          % (len(bams), virus_cat, len(virus_refs), len(human_refs)))

    sample_rows, real_names = [], []
    n_prefiltered = 0
    for bam in bams:
        sample = sample_from_bam(bam)
        sid = anon.get(sample, "S??")
        real_names.append(sample)
        row = {
            "real": sample, "run": run_name, "group": group_label(sample, run_name),
            "bam": bam, "bam_present": "yes", "prefiltered_zero": "no",
            "virus_refs_queried": 0, "reads_passing_filter": 0,
            "n_discordant_human": 0, "n_discordant_nonhuman": 0, "n_clipped_reads": 0,
            "clusters": [], "note": "",
        }

        if not bam_index_present(bam):
            warn("BAM index for " + sid, masked(bam + ".bai", sample, sid))
            row["bam_present"] = "no_index"
            row["note"] = "no_index"
            sample_rows.append(row)
            continue

        present = bam_references(samtools, bam)
        if present is None:
            warn("BAM header for %s (samtools view -H failed)" % sid,
                 masked(bam, sample, sid))
            row["bam_present"] = "header_unreadable"
            row["note"] = "header_unreadable"
            sample_rows.append(row)
            continue

        if virus_refs:
            targets = [r for r in present if r in virus_refs]
        else:
            token = virus_cat.replace("-", "").replace("_", "")
            targets = [r for r in present
                       if token in r.upper().replace("-", "").replace("_", "")]
        if human_refs:
            human_now = human_refs
        else:
            human_now = set(r for r in present if HUMAN_FALLBACK_RE.match(r))
        if not targets:
            row["note"] = "no_%s_reference_in_bam" % virus_cat.lower()
            sample_rows.append(row)
            continue
        row["virus_refs_queried"] = len(targets)

        if not args.no_prefilter:
            raw = idxstats_virus_reads(run_dir, sample, set(targets))
            if raw == 0:
                row["prefiltered_zero"] = "yes"
                row["note"] = "zero_raw_virus_reads"
                n_prefiltered += 1
                sample_rows.append(row)
                continue

        for ref in targets:
            events, kept, nonhuman, err = scan_virus_reference(
                samtools, bam, ref, args, human_now, virus_refs or set(targets))
            row["reads_passing_filter"] += kept
            row["n_discordant_nonhuman"] += nonhuman
            if err:
                # samtools echoes the BAM path in its own errors: mask both
                warn("samtools view on %s for %s (%s)"
                     % (ref, sid, masked(err, sample, sid)),
                     masked(bam, sample, sid))
                row["note"] = (row["note"] + ";samtools_error").strip(";")
                continue
            row["n_discordant_human"] += sum(1 for e in events
                                             if e["kind"] == "discordant")
            row["n_clipped_reads"] += sum(1 for e in events if e["kind"] == "clip")
            for events_in_cluster in cluster_positions(events, args.cluster_bp):
                summary = summarise_cluster(events_in_cluster, args.mate_cluster_bp)
                summary["virus_ref"] = ref
                row["clusters"].append(summary)
        sample_rows.append(row)

    if n_prefiltered:
        print("PREFILTER %d sample(s) had zero raw %s reads in idxstats and were "
              "not opened" % (n_prefiltered, virus_cat))
    print("")
    return sample_rows, real_names


# ---------------------------------------------------------------- output

def write_tsv(path, header, rows):
    with open(path, "w", encoding="ascii", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(h, "")) for h in header) + "\n")


CAND_HEADER = [
    "sample", "group", "run", "virus_ref", "cluster_id", "junction_pos",
    "cluster_start", "cluster_end", "n_support_total", "n_discordant",
    "n_clipped", "n_clip_left", "n_clip_right", "max_clip_len",
    "support_class", "n_host_sites", "best_host_chrom", "best_host_pos",
    "best_host_support", "host_sites_top3", "pass_min_support", "review_status",
]

SAMPLE_HEADER = [
    "sample", "group", "run", "bam_present", "prefiltered_zero",
    "virus_refs_queried", "reads_passing_filter", "n_discordant_human",
    "n_discordant_nonhuman", "n_clipped_reads", "n_clusters_all",
    "n_clusters_pass", "n_clusters_both_signals", "top_cluster_virus_ref",
    "top_cluster_pos", "top_cluster_support", "top_cluster_support_frac",
    "distinct_host_chroms", "junction_call", "note",
]


def main():
    ap = argparse.ArgumentParser(
        description="Candidate HTLV-1 host-virus integration junctions "
                    "(discordant pairs + soft clips). Candidates only - "
                    "IGV review required.")
    ap.add_argument("--run-dir", action="append", default=None,
                    help="run directory; repeatable. Default: the targeted "
                         "current-reference run under " + SHUYU_ROOT)
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR,
                    help="output directory (created if absent)")
    ap.add_argument("--prefix", default="htlv_junction",
                    help="output filename prefix (default: htlv_junction)")
    ap.add_argument("--refmap", default=None,
                    help="reference map CSV; default is the map the run "
                         "recorded in results/run_manifest.tsv, else "
                         + DEFAULT_BASE_REFMAP)
    ap.add_argument("--bam-subdir", default="bam",
                    help="BAM subdirectory inside a run (default: bam)")
    ap.add_argument("--samtools", default="samtools",
                    help="samtools executable (default: samtools)")
    ap.add_argument("--threads", type=int, default=0,
                    help="samtools decompression threads, 0 = do not pass -@")
    ap.add_argument("--virus-category", default="HTLV1",
                    help="reference-map category of the virus (default: HTLV1)")
    ap.add_argument("--human-category", default="HUMAN",
                    help="reference-map category of host refs (default: HUMAN)")
    ap.add_argument("--mapq", type=int, default=40,
                    help="minimum MAPQ on the virus reference (default: 40)")
    ap.add_argument("--min-mate-mapq", type=int, default=20,
                    help="minimum mate MAPQ when the MQ tag is present "
                         "(default: 20)")
    ap.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0xD04,
                    help="samtools -F value; default 0xD04 = unmapped, "
                         "secondary, duplicate, supplementary")
    ap.add_argument("--min-clip", type=int, default=20,
                    help="minimum soft-clip length in bp (default: 20)")
    ap.add_argument("--max-homopolymer-frac", type=float, default=0.9,
                    help="drop clips whose sequence is at least this fraction "
                         "one base (default: 0.9)")
    ap.add_argument("--cluster-bp", type=int, default=50,
                    help="junction clustering window on the virus (default: 50)")
    ap.add_argument("--mate-cluster-bp", type=int, default=500,
                    help="host mate clustering window in bp (default: 500)")
    ap.add_argument("--min-support", type=int, default=2,
                    help="supporting reads for a cluster to be reportable "
                         "(default: 2)")
    ap.add_argument("--strong-support", type=int, default=5,
                    help="support needed, with both signals, for a 'strong' "
                         "per-sample junction call (default: 5)")
    ap.add_argument("--allow-ambiguous", action="store_true",
                    help="keep reads with AS<=XS (default: unique-best only)")
    ap.add_argument("--no-prefilter", action="store_true",
                    help="do not use results/<sample>.idxstats.tsv to skip "
                         "samples with zero raw virus reads")
    ap.add_argument("--limit", type=int, default=0,
                    help="process only the first N BAMs per run (debugging)")
    args = ap.parse_args()

    runs = args.run_dir if args.run_dir else list(DEFAULT_RUNS)

    try:
        os.makedirs(args.outdir, exist_ok=True)
    except OSError:
        warn("output directory (cannot create)", args.outdir)
        return 0

    cand_path = os.path.join(args.outdir, args.prefix + "_candidates.tsv")
    samp_path = os.path.join(args.outdir, args.prefix + "_per_sample.tsv")
    key_path = os.path.join(args.outdir, args.prefix + "_sample_key.tsv")

    print("a6_htlv_junctions - CANDIDATE HTLV-1 integration junctions")
    print("date 2026-07-26 | virus category %s | MAPQ>=%d | -F 0x%X | "
          "min clip %d bp | cluster %d bp | %s"
          % (args.virus_category.upper(), args.mapq, args.exclude_flags,
             args.min_clip, args.cluster_bp,
             "any alignment" if args.allow_ambiguous else "unique-best AS>XS"))
    print("")

    if not samtools_ok(args.samtools):
        warn("samtools executable", args.samtools)
        write_tsv(cand_path, CAND_HEADER, [])
        write_tsv(samp_path, SAMPLE_HEADER, [])
        print("Wrote headers only: %s" % cand_path)
        print("Wrote headers only: %s" % samp_path)
        return 0

    # Enumerate every BAM first so the real -> S01..Snn map exists before any
    # per-sample warning is printed; nothing on stdout is then an identifier.
    pre_names = []
    for run_dir in runs:
        if not os.path.isdir(run_dir):
            continue
        for bam in list_run_bams(run_dir, args)[1]:
            pre_names.append(sample_from_bam(bam))
    anon = anon_ids(pre_names)

    all_rows, all_names = [], []
    for run_dir in runs:
        rows, names = process_run(run_dir, args, args.samtools, anon)
        all_rows.extend(rows)
        all_names.extend(names)

    for name in all_names:                 # belt and braces: nothing unmapped
        if name not in anon:
            anon = anon_ids(pre_names + all_names)
            break

    cand_rows, samp_rows = [], []
    for row in sorted(all_rows, key=lambda r: (anon.get(r["real"], "ZZZ"), r["run"])):
        sid = anon.get(row["real"], "S00")
        clusters = sorted(row["clusters"],
                          key=lambda c: (-c["n_support_total"], c["virus_ref"],
                                         c["junction_pos"]))
        passing = [c for c in clusters if c["n_support_total"] >= args.min_support]
        both = [c for c in passing if c["support_class"] == "both"]
        strong = [c for c in both if c["n_support_total"] >= args.strong_support]
        for i, clu in enumerate(clusters, start=1):
            ok = clu["n_support_total"] >= args.min_support
            out = dict(clu)
            out.update({
                "sample": sid, "group": row["group"], "run": row["run"],
                "cluster_id": "%s_%s_%d" % (sid, clu["virus_ref"], i),
                "pass_min_support": "yes" if ok else "no",
                "review_status": "CANDIDATE_IGV_REVIEW_REQUIRED" if ok
                                 else "LOW_SUPPORT",
            })
            cand_rows.append(out)

        total_support = sum(c["n_support_total"] for c in passing)
        top = passing[0] if passing else None
        chroms = set(c["best_host_chrom"] for c in passing
                     if c["best_host_chrom"] not in ("NA", ""))
        if strong:
            call = "strong"
        elif passing:
            call = "candidate"
        else:
            call = "none"
        samp_rows.append({
            "sample": sid, "group": row["group"], "run": row["run"],
            "bam_present": row["bam_present"],
            "prefiltered_zero": row["prefiltered_zero"],
            "virus_refs_queried": row["virus_refs_queried"],
            "reads_passing_filter": row["reads_passing_filter"],
            "n_discordant_human": row["n_discordant_human"],
            "n_discordant_nonhuman": row["n_discordant_nonhuman"],
            "n_clipped_reads": row["n_clipped_reads"],
            "n_clusters_all": len(clusters),
            "n_clusters_pass": len(passing),
            "n_clusters_both_signals": len(both),
            "top_cluster_virus_ref": top["virus_ref"] if top else "NA",
            "top_cluster_pos": top["junction_pos"] if top else 0,
            "top_cluster_support": top["n_support_total"] if top else 0,
            "top_cluster_support_frac": ("%.3f" % (top["n_support_total"]
                                                  / float(total_support)))
                                        if top and total_support else "NA",
            "distinct_host_chroms": len(chroms),
            "junction_call": call,
            "note": row["note"] or "ok",
        })

    write_tsv(cand_path, CAND_HEADER, cand_rows)
    write_tsv(samp_path, SAMPLE_HEADER, samp_rows)

    with open(key_path, "w", encoding="ascii", newline="") as fh:
        fh.write(KEY_WARNING + "\n")
        fh.write("# generated 2026-07-26 by a6_htlv_junctions.py\n")
        fh.write("anon_sample\treal_sample\tgroup\trun\tbam_path\n")
        seen = set()
        for row in sorted(all_rows, key=lambda r: r["real"]):
            tag = (row["real"], row["run"])
            if tag in seen:
                continue
            seen.add(tag)
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (anon.get(row["real"], "S00"),
                                               row["real"], row["group"],
                                               row["run"], row["bam"]))

    # ------------------------------------------------------------ stdout
    n_samples = len(samp_rows)
    with_cand = [r for r in samp_rows if r["junction_call"] != "none"]
    with_strong = [r for r in samp_rows if r["junction_call"] == "strong"]
    print("-- summary (anonymised) --")
    print("samples processed            %d" % n_samples)
    print("samples with a passing junction cluster (>=%d reads)  %d"
          % (args.min_support, len(with_cand)))
    print("samples with both signals and >=%d reads              %d"
          % (args.strong_support, len(with_strong)))
    print("junction clusters written    %d (%d passing --min-support)"
          % (len(cand_rows), sum(1 for c in cand_rows
                                 if c["pass_min_support"] == "yes")))
    by_group = {}
    for row in samp_rows:
        key = row["group"]
        hit = 1 if row["junction_call"] != "none" else 0
        tot, pos = by_group.get(key, (0, 0))
        by_group[key] = (tot + 1, pos + hit)
    for key in sorted(by_group):
        tot, pos = by_group[key]
        print("  group %-4s %3d samples, %3d with candidates" % (key, tot, pos))
    if with_cand:
        print("")
        print("-- top samples by cluster count (clonality proxy) --")
        print("%-6s %-4s %8s %8s %10s %12s %9s"
              % ("sample", "grp", "clusters", "both", "top_supp", "top_pos", "hostchr"))
        ranked = sorted(with_cand, key=lambda r: (-r["n_clusters_pass"],
                                                  -r["top_cluster_support"]))
        for row in ranked[:15]:
            print("%-6s %-4s %8d %8d %10d %12s %9d"
                  % (row["sample"], row["group"], row["n_clusters_pass"],
                     row["n_clusters_both_signals"], row["top_cluster_support"],
                     row["top_cluster_pos"], row["distinct_host_chroms"]))
    print("")
    print("wrote %s" % cand_path)
    print("wrote %s" % samp_path)
    print("wrote %s  <- %s" % (key_path, KEY_WARNING.lstrip("# ")))
    print("")
    print("NOTE: every row is a CANDIDATE. Targeted capture produces chimeric "
          "artefacts, so confirm each cluster in IGV (clip consensus, mate "
          "mapping quality, host repeat context) before calling integration.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
