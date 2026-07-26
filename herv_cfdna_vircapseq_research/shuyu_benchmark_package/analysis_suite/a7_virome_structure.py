#!/usr/bin/env python3
"""a7_virome_structure.py -- anellovirus burden and coinfection structure.

WHAT IT COMPUTES
  Two analyses over the per-sample x per-reference counts of the masked-panel runs.

  (1) ANELLOVIRUS (TTV / TTMV / TTMDV) BURDEN, used as an immunocompetence proxy.
      Per sample: total anellovirus reads, number of distinct anellovirus
      references detected (richness), and Shannon diversity (natural log) over the
      detected references. Raw counts and counts normalised per million usable
      reads (RPM) are both reported; the denominator comes from the sample's own
      samtools idxstats totals (or, with --norm-source filtered_categories, from
      the run's headline deduplicated category table). HIV vs HL is then compared
      with a Mann-Whitney U test implemented here in the standard library
      (U for each group, tie-corrected normal-approximation two-sided p with
      continuity correction, group medians and means, rank-biserial effect size).

  (2) COINFECTION STRUCTURE. A sample x virus-group presence/absence matrix over
      EBV, OTHER_HERPES, ANELLO, ADENO, POLYOMA, HBV; pairwise co-occurrence
      counts with the Jaccard index (whole cohort and within each group label);
      and the per-sample number of distinct virus groups detected, again compared
      between HIV and HL with the same Mann-Whitney U test.

  Anelloviruses are identified from the reference map (description / accession
  keywords such as "Torque teno", "anello", TTV/TTMV/TTMDV); extra accessions or
  reference ids can be supplied with --anello-accessions. The chimpanzee-isolate
  references (default NC_014069.1, NC_014077.1, NC_014480.2) are reported in
  separate flagged columns and are excluded from every human-diversity metric and
  from the ANELLO coinfection column.

  Per-reference counts come from results/<sample>.idxstats.tsv (mapped column) by
  default. With --counts unique_best the module instead calls samtools on the BAMs
  and counts primary, non-supplementary, MAPQ>=--mapq reads whose alignment is
  unambiguous (AS>XS, or no XS tag); the normalisation denominators still come
  from idxstats. No figures are produced, so matplotlib is not imported.

WHAT IT WRITES (tab-separated, into --outdir, names prefixed with --prefix)
  <prefix>_sample_key.tsv           real -> anonymous sample mapping. THE ONLY
                                    file that contains real sample identifiers.
  <prefix>_anellovirus_burden.tsv   per-sample anellovirus burden, raw + RPM.
  <prefix>_anellovirus_group_test.tsv  Mann-Whitney results for the burden
                                    metrics and for n_virus_groups_detected.
  <prefix>_coinfection_matrix.tsv   per-sample presence/absence + read counts.
  <prefix>_coinfection_pairs.tsv    pairwise co-occurrence counts and Jaccard.
  <prefix>_virus_group_refs.tsv     audit trail: which reference went to which
                                    virus group and why (no sample identifiers).
  Every file except the sample key uses anonymous ids S01..Snn only.

EXAMPLE
  python3 a7_virome_structure.py \
      --runs /path/to/runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only \
             /path/to/runs/targeted_htlv_hg38_shuyu_masked_panel_refixed_primary_only \
      --outdir ~/shuyu_project/local_work/panel_report_20260725/a7_virome_out \
      --min-reads 10

  # unique-best recount from the BAMs instead of raw idxstats
  python3 a7_virome_structure.py --counts unique_best --mapq 40 --samtools samtools
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import itertools
import math
import os
import re
import subprocess
import sys
import tempfile

SHUYU_RUNS = "/path/to/runs"

DEFAULT_RUNS = [
    os.path.join(SHUYU_RUNS, "wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only"),
    os.path.join(SHUYU_RUNS, "targeted_htlv_hg38_shuyu_masked_panel_refixed_primary_only"),
]
DEFAULT_REFMAP = os.path.join(
    SHUYU_RUNS,
    "shuyu_masked_panel_hg38_herv_line1_refixed",
    "ref",
    "hg38_herv_line1_plus_shuyu_masked_panel.reference_map.csv",
)
DEFAULT_BASE_REFMAP = os.path.join(
    SHUYU_RUNS,
    "retro_reference_hg38_refseq",
    "ref",
    "hg38_plus_retro.refseq.reference_map.csv",
)
DEFAULT_OUTDIR = os.path.expanduser(
    "~/shuyu_project/local_work/panel_report_20260725/a7_virome_out")
DEFAULT_ANELLO_ACC_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "anello_accessions.txt")
DEFAULT_CHIMP_ACC = "NC_014069.1,NC_014077.1,NC_014480.2"

# Categories in the reference map that can never be a virus group of interest.
NON_VIRAL_CATEGORIES = ("HUMAN", "HERV", "LINE1")

# Virus groups, in priority order: the first pattern set that matches wins, so
# EBV must be tested before the generic herpesvirus patterns.
GROUP_ORDER = ["EBV", "OTHER_HERPES", "ANELLO", "ADENO", "POLYOMA", "HBV"]
CHIMP_COL = "ANELLO_CHIMP_FLAGGED"

GROUP_KEYWORDS = {
    "ANELLO": [
        "anello", "torque teno", "torque-teno", "torquetenovirus",
        "transfusion transmitted virus", "transfusion-transmitted virus",
        "small anellovirus", "tt virus", "ttv", "ttmv", "ttmdv", "sen virus",
    ],
    "EBV": [
        "epstein", "epstein-barr", "human herpesvirus 4",
        "human gammaherpesvirus 4", "lymphocryptovirus", "ebv",
    ],
    "OTHER_HERPES": [
        "herpesvirus", "herpes simplex", "cytomegalovirus", "roseolovirus",
        "varicella", "simplexvirus", "rhadinovirus", "muromegalovirus",
        "hhv", "cmv", "hsv", "vzv", "kshv",
    ],
    "ADENO": ["adenovirus", "mastadenovirus"],
    "POLYOMA": [
        "polyomavirus", "polyoma virus", "merkel cell", "jc virus", "bk virus",
        "simian virus 40", "sv40", "mcpyv",
    ],
    "HBV": ["hepatitis b virus", "hepadnavirus", "hbv"],
}

# Short keywords that must not match inside a longer word.
SHORT_KEYWORDS = set([
    "ttv", "ttmv", "ttmdv", "ebv", "hhv", "cmv", "hsv", "vzv", "kshv", "hbv",
    "sv40", "mcpyv",
])

# Accessions used as a backstop when a panel header carries no description.
GROUP_ACCESSIONS = {
    "EBV": ["NC_007605", "NC_009334"],
    "OTHER_HERPES": [
        "NC_001664", "NC_000898", "NC_001716", "NC_006273", "NC_001806",
        "NC_001798", "NC_001348", "NC_009333",
    ],
    "ANELLO": ["NC_002076"],
}

ACC_RE = re.compile(
    r"(?<![A-Za-z0-9])([A-Z]{2}_\d{5,8}|[A-Z]{1,2}\d{5,6})(?:\.\d+)?(?![A-Za-z0-9])")

TODAY = datetime.date.today().isoformat()
SCRIPT = os.path.basename(__file__)


# --------------------------------------------------------------------------- #
# small utilities
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


def shannon(counts):
    """Shannon diversity (natural log) over a list of positive counts."""
    vals = [float(c) for c in counts if c > 0]
    total = sum(vals)
    if total <= 0 or len(vals) < 2:
        return 0.0
    h = 0.0
    for v in vals:
        p = v / total
        h -= p * math.log(p)
    return h


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
        p = math.erfc(abs(z) / math.sqrt(2.0))
        p = min(1.0, max(0.0, p))
    return {
        "n1": n1, "n2": n2, "U1": u1, "U2": u2, "z": z, "p": p,
        "median1": median(x), "median2": median(y),
        "mean1": mean(x), "mean2": mean(y),
        "effect_r": 2.0 * u1 / (float(n1) * float(n2)) - 1.0,
        "n_ties_groups": sum(1 for t in tie_sizes if t > 1),
    }


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
# sample naming
# --------------------------------------------------------------------------- #
def group_of(sample, run_name, use_run_name=True):
    """HIV / HL / TCL / NA from the real sample name (spec-defined rules)."""
    low = sample.lower()
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
    """real sample name -> S01..Snn, ordered by sorted real name."""
    uniq = sorted(set(real_names))
    width = max(2, len(str(len(uniq))))
    return dict((name, "S" + str(i + 1).zfill(width)) for i, name in enumerate(uniq))


# --------------------------------------------------------------------------- #
# reference map + virus-group classification
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


def refmap_from_manifest(run_dir):
    man = os.path.join(run_dir, "results", "run_manifest.tsv")
    if not os.path.exists(man):
        return None
    try:
        with open(man, newline="", encoding="utf-8", errors="replace") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except (OSError, csv.Error):
        return None
    if not rows:
        return None
    for key in rows[0]:
        if "reference_map" in (key or "").lower():
            val = (rows[0][key] or "").strip()
            if val and os.path.exists(val):
                return val
    for val in rows[0].values():
        val = (val or "").strip()
        if val.endswith(".csv") and "reference_map" in val and os.path.exists(val):
            return val
    return None


def load_refmap(path):
    """reference_id -> dict(category, text, accs). Empty dict if unreadable."""
    refs = {}
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
                text = " ".join(p for p in parts if p)
                refs[rid] = {
                    "category": (row.get("category") or "").strip().upper(),
                    "description": to_ascii(row.get("description") or ""),
                    "text": to_ascii(text),
                    "accs": set(norm_acc(a) for a in accessions_in(to_ascii(text))),
                }
    except OSError:
        return {}
    return refs


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


def classify_references(refs, extra_anello, chimp_acc):
    """reference_id -> (group_or_None, chimp_flag, reason)."""
    extra_norm = set(norm_acc(a) for a in extra_anello)
    extra_ids = set(a.strip() for a in extra_anello)
    chimp_norm = set(norm_acc(a) for a in chimp_acc)
    assign = {}
    for rid, info in refs.items():
        if info["category"] in NON_VIRAL_CATEGORIES:
            assign[rid] = (None, False, "non_viral_category")
            continue
        low = info["text"].lower()
        group, reason = None, ""
        is_chimp = bool(info["accs"] & chimp_norm) or rid in chimp_acc
        for cand in GROUP_ORDER:
            for kw in GROUP_KEYWORDS.get(cand, []):
                if keyword_hit(low, kw):
                    group, reason = cand, "keyword:" + kw
                    break
            if group:
                break
            for acc in GROUP_ACCESSIONS.get(cand, []):
                if norm_acc(acc) in info["accs"]:
                    group, reason = cand, "accession:" + acc
                    break
            if group:
                break
        if group is None and (info["accs"] & extra_norm or rid in extra_ids):
            group, reason = "ANELLO", "anello_accession_file"
        if group is None and is_chimp:
            group, reason = "ANELLO", "chimp_accession_list"
        if group != "ANELLO":
            is_chimp = False
        elif not is_chimp:
            if "chimpanzee" in low or "pan troglodytes" in low or "chimp" in low:
                is_chimp = True
        if group is None:
            assign[rid] = (None, False, "unclassified")
        else:
            assign[rid] = (group, is_chimp, reason or "matched")
    return assign


# --------------------------------------------------------------------------- #
# per-run data loading
# --------------------------------------------------------------------------- #
def parse_idxstats(path, wanted_refs):
    """-> (per_ref counts for wanted_refs, total_mapped, star_unmapped)."""
    per_ref = {}
    total_mapped = 0
    star_unmapped = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                name = parts[0]
                try:
                    mapped = int(parts[2])
                    unmapped = int(parts[3])
                except ValueError:
                    continue
                if name == "*":
                    star_unmapped += unmapped
                    continue
                total_mapped += mapped
                if mapped and name in wanted_refs:
                    per_ref[name] = per_ref.get(name, 0) + mapped
    except OSError:
        return None, 0, 0
    return per_ref, total_mapped, star_unmapped


def parse_category_counts(path):
    """-> dict sample -> dict category -> int. None if unreadable."""
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
                sample = row[0].strip()
                cats = {}
                for idx in range(1, min(len(row), len(header))):
                    try:
                        cats[header[idx]] = int(float(row[idx]))
                    except (ValueError, TypeError):
                        continue
                table[sample] = cats
    except OSError:
        return None
    return table


def find_filtered_category_counts(results_dir):
    """Headline deduplicated category table, with or without a filename prefix."""
    hits = sorted(glob.glob(os.path.join(results_dir, "*filtered_category_counts.tsv")))
    hits = [h for h in hits if "record_category_counts" not in os.path.basename(h)]
    return hits[0] if hits else None


def bam_for_sample(run_dir, sample):
    cand = os.path.join(run_dir, "bam", sample + ".bam")
    if os.path.exists(cand):
        return cand
    hits = sorted(glob.glob(os.path.join(run_dir, "bam", sample + "*.bam")))
    return hits[0] if hits else None


def bam_is_indexed(bam):
    return os.path.exists(bam + ".bai") or os.path.exists(
        os.path.splitext(bam)[0] + ".bai") or os.path.exists(bam + ".csi")


def unique_best_counts(bam, refs, samtools, mapq, exclude_flags, chunk=200):
    """Per-reference unique-best (AS>XS or no XS) primary read counts.

    Returns (counts, error_string). counts is None when samtools could not run.
    """
    counts = dict((r, 0) for r in refs)
    ref_list = sorted(refs)
    for start in range(0, len(ref_list), chunk):
        regions = ref_list[start:start + chunk]
        cmd = [samtools, "view", "-F", str(exclude_flags), "-q", str(mapq), bam] + regions
        # stderr goes to a temp file, not a pipe: it is only drained after the
        # stdout loop finishes, and a pipe that filled up first would deadlock.
        err_fh = tempfile.TemporaryFile(mode="w+", errors="replace")
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=err_fh, universal_newlines=True)
        except OSError as exc:
            err_fh.close()
            return None, "samtools not runnable (%s)" % exc
        assert proc.stdout is not None
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            rname = fields[2]
            if rname not in counts:
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
            if a_score is None:
                continue
            if x_score is None or a_score > x_score:
                counts[rname] = counts[rname] + 1
        proc.stdout.close()
        rc = proc.wait()
        try:
            err_fh.seek(0)
            err = err_fh.read()
        except (IOError, OSError, ValueError):
            err = ""
        err_fh.close()
        if rc != 0:
            return None, "samtools view rc=%d %s" % (rc, to_ascii(err)[:160])
    return counts, ""


# --------------------------------------------------------------------------- #
# main assembly
# --------------------------------------------------------------------------- #
def collect_run(run_dir, args, refs_cache, samtools_state):
    """-> (list of per-sample records, assignment dict, refmap path) for one run."""
    run_name = os.path.basename(os.path.normpath(run_dir))
    if not os.path.isdir(run_dir):
        warn("run directory", run_dir)
        return [], None, None
    results = os.path.join(run_dir, "results")
    if not os.path.isdir(results):
        warn("results directory for run " + run_name, results)
        return [], None, None
    idx_files = sorted(glob.glob(os.path.join(results, "*.idxstats.tsv")))
    if not idx_files:
        warn("per-sample idxstats for run " + run_name,
             os.path.join(results, "*.idxstats.tsv"))
        return [], None, None

    refmap_path = None
    if args.use_manifest_refmap:
        refmap_path = refmap_from_manifest(run_dir)
    if not refmap_path:
        refmap_path = args.refmap
    if not os.path.exists(refmap_path):
        warn("reference map for run " + run_name, refmap_path)
        if os.path.exists(args.base_refmap):
            print("NOTE: falling back to the base reference map %s" % args.base_refmap)
            refmap_path = args.base_refmap
        else:
            warn("fallback base reference map", args.base_refmap)
            return [], None, None

    if refmap_path not in refs_cache:
        refs = load_refmap(refmap_path)
        if not refs:
            warn("readable reference map for run " + run_name, refmap_path)
            return [], None, None
        assign = classify_references(refs, args._anello_extra, args._chimp_acc)
        refs_cache[refmap_path] = (refs, assign)
    refs, assign = refs_cache[refmap_path]

    group_refs = {}
    for name in GROUP_ORDER:
        group_refs[name] = []
    chimp_refs = []
    for rid, (grp, chimp, _reason) in assign.items():
        if grp is None:
            continue
        if grp == "ANELLO" and chimp:
            chimp_refs.append(rid)
        else:
            group_refs[grp].append(rid)
    wanted = set(chimp_refs)
    for name in GROUP_ORDER:
        wanted.update(group_refs[name])

    cat_path = find_filtered_category_counts(results)
    cat_table = None
    if cat_path:
        cat_table = parse_category_counts(cat_path)
        if cat_table is None:
            warn("readable filtered category counts for run " + run_name, cat_path)
    else:
        warn("filtered_category_counts.tsv for run " + run_name,
             os.path.join(results, "*filtered_category_counts.tsv"))

    records = []
    no_bam = 0
    no_index = 0
    for idx_path in idx_files:
        sample = os.path.basename(idx_path)
        for suffix in (".idxstats.tsv",):
            if sample.endswith(suffix):
                sample = sample[: -len(suffix)]
        per_ref, total_mapped, star_unmapped = parse_idxstats(idx_path, wanted)
        if per_ref is None:
            # idx_path embeds the real sample name; print a masked path instead
            warn("readable idxstats for a sample in run " + run_name,
                 os.path.join(results, "<sample>.idxstats.tsv"))
            continue
        counts_source = "idxstats_mapped"
        if args.counts == "unique_best":
            bam = bam_for_sample(run_dir, sample)
            if not bam:
                no_bam += 1
                if no_bam == 1:                      # warn once per run, not per sample
                    warn("bam for a sample in run " + run_name,
                         os.path.join(run_dir, "bam", "<sample>.bam"))
            elif not bam_is_indexed(bam):
                no_index += 1
                if no_index == 1:
                    # bam embeds the real sample name; print a masked path
                    warn("bam index for a sample in run " + run_name,
                         os.path.join(run_dir, "bam", "<sample>.bam.bai"))
            elif samtools_state["broken"]:
                pass
            else:
                ub, err = unique_best_counts(bam, wanted, args.samtools, args.mapq,
                                             args.exclude_flags)
                if ub is None:
                    if not samtools_state["broken"]:
                        # samtools echoes the BAM path in its errors: mask it
                        print("WARN: unique-best counting unavailable (%s); "
                              "using idxstats mapped counts instead"
                              % str(err).replace(sample, "<sample>"))
                    samtools_state["broken"] = True
                else:
                    per_ref = dict((k, v) for k, v in ub.items() if v)
                    counts_source = "unique_best_mapq%d" % args.mapq

        filtered_total = None
        if cat_table and sample in cat_table:
            filtered_total = sum(cat_table[sample].values())

        records.append({
            "sample": sample,
            "run": run_name,
            "group": group_of(sample, run_name, not args.no_run_name_group),
            "counts": per_ref,
            "counts_source": counts_source,
            "total_mapped": total_mapped,
            "star_unmapped": star_unmapped,
            "total_reads": total_mapped + star_unmapped,
            "filtered_total": filtered_total,
            "group_refs": group_refs,
            "chimp_refs": chimp_refs,
        })
    if no_bam or no_index:
        print("WARN: run %s -- %d samples without a bam and %d without a bam index; "
              "those samples fall back to idxstats mapped counts"
              % (run_name, no_bam, no_index))
    return records, assign, refmap_path


def denominator(rec, norm_source):
    if norm_source == "idxstats_total":
        val = rec["total_reads"]
    elif norm_source == "filtered_categories":
        val = rec["filtered_total"]
    else:
        val = rec["total_mapped"]
    if val is None or val <= 0:
        return None
    return float(val)


def summarise_sample(rec, min_reads, norm_source):
    counts = rec["counts"]
    group_refs = rec["group_refs"]
    chimp_refs = rec["chimp_refs"]

    anello_all = dict((r, counts.get(r, 0)) for r in group_refs["ANELLO"])
    anello_pos = dict((r, c) for r, c in anello_all.items() if c >= min_reads)
    chimp_all = dict((r, counts.get(r, 0)) for r in chimp_refs)
    chimp_pos = dict((r, c) for r, c in chimp_all.items() if c >= min_reads)

    total_anello = sum(anello_all.values())
    detected_anello = sum(anello_pos.values())
    top_ref, top_reads = "NA", 0
    if anello_pos:
        top_ref, top_reads = sorted(anello_pos.items(), key=lambda t: (-t[1], t[0]))[0]
    top_share = (float(top_reads) / float(detected_anello)) if detected_anello else None

    denom = denominator(rec, norm_source)
    rpm_h = (1e6 * total_anello / denom) if denom else None
    rpm_c = (1e6 * sum(chimp_all.values()) / denom) if denom else None

    group_reads = {}
    present = {}
    for name in GROUP_ORDER:
        reads = sum(counts.get(r, 0) for r in group_refs[name])
        group_reads[name] = reads
        present[name] = 1 if reads >= min_reads else 0
    chimp_reads = sum(chimp_all.values())
    chimp_present = 1 if chimp_reads >= min_reads else 0

    out = dict(rec)
    out.update({
        "anello_total": total_anello,
        "anello_detected_reads": detected_anello,
        "anello_richness": len(anello_pos),
        "anello_shannon": shannon(list(anello_pos.values())),
        "anello_top_ref": top_ref,
        "anello_top_reads": top_reads,
        "anello_top_share": top_share,
        "chimp_reads": chimp_reads,
        "chimp_richness": len(chimp_pos),
        "chimp_present": chimp_present,
        "anello_rpm": rpm_h,
        "chimp_rpm": rpm_c,
        "denom": denom,
        "group_reads": group_reads,
        "present": present,
        "n_groups": sum(present.values()),
    })
    return out


# --------------------------------------------------------------------------- #
# writers
# --------------------------------------------------------------------------- #
def write_sample_key(path, samples, anon):
    lines = [
        "CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL",
        "generated %s by %s" % (TODAY, SCRIPT),
    ]
    rows = []
    for rec in sorted(samples, key=lambda r: (r["sample"], r["run"])):
        rows.append([anon[rec["sample"]], rec["sample"], rec["group"], rec["run"]])
    write_tsv(path, lines, ["anon_sample", "real_sample", "group", "run"], rows)


def common_comments(args, extra=()):
    lines = [
        "%s  generated %s" % (SCRIPT, TODAY),
        "no real sample identifiers in this file; ids are anonymous S01..Snn "
        "(mapping in %s)" % out_name(args.prefix, "sample_key.tsv"),
        "presence / richness threshold: reads >= %d per reference or virus group"
        % args.min_reads,
        "normalisation denominator: %s (RPM = 1e6 * reads / denominator)"
        % args.norm_source,
    ]
    lines.extend(extra)
    return lines


def write_burden(path, samples, anon, args):
    header = [
        "sample", "group", "run", "counts_source", "min_reads",
        "anello_reads_human_total", "anello_reads_human_detected",
        "anello_richness_human", "anello_shannon_human",
        "anello_top_ref", "anello_top_ref_reads", "anello_top_ref_share",
        "anello_reads_chimp_flagged", "anello_richness_chimp_flagged",
        "total_mapped_idxstats", "total_unmapped_idxstats", "total_reads_idxstats",
        "total_filtered_categories", "norm_source", "norm_denominator",
        "anello_rpm_human", "anello_rpm_chimp_flagged",
    ]
    rows = []
    for rec in samples:
        rows.append([
            anon[rec["sample"]], rec["group"], rec["run"], rec["counts_source"],
            str(args.min_reads),
            str(rec["anello_total"]), str(rec["anello_detected_reads"]),
            str(rec["anello_richness"]), fnum(rec["anello_shannon"], 4),
            rec["anello_top_ref"], str(rec["anello_top_reads"]),
            fnum(rec["anello_top_share"], 4),
            str(rec["chimp_reads"]), str(rec["chimp_richness"]),
            str(rec["total_mapped"]), str(rec["star_unmapped"]),
            str(rec["total_reads"]),
            "NA" if rec["filtered_total"] is None else str(rec["filtered_total"]),
            args.norm_source,
            "NA" if rec["denom"] is None else str(int(rec["denom"])),
            fnum(rec["anello_rpm"], 3), fnum(rec["chimp_rpm"], 3),
        ])
    write_tsv(path, common_comments(
        args, ["chimpanzee-isolate anellovirus references (%s) are reported in the "
               "*_chimp_flagged columns only and excluded from the human metrics"
               % ",".join(args._chimp_acc)]), header, rows)


METRIC_SPECS = [
    ("anello_reads_human_total", "raw_reads", "anello_total"),
    ("anello_rpm_human", "normalised_rpm", "anello_rpm"),
    ("anello_richness_human", "count_of_references", "anello_richness"),
    ("anello_shannon_human", "shannon_index_natural_log", "anello_shannon"),
    ("n_virus_groups_detected", "count_of_virus_groups", "n_groups"),
]


def group_test_rows(samples, g1, g2):
    rows = []
    results = []
    for metric, scale, key in METRIC_SPECS:
        x = [r[key] for r in samples if r["group"] == g1 and r[key] is not None]
        y = [r[key] for r in samples if r["group"] == g2 and r[key] is not None]
        res = mann_whitney_u(x, y)
        if res is None:
            rows.append([metric, scale, g1, str(len(x)), g2, str(len(y))]
                        + ["NA"] * 9
                        + ["one or both groups empty; test not run"])
            results.append((metric, None))
            continue
        note = "two-sided normal approximation, tie- and continuity-corrected"
        if min(res["n1"], res["n2"]) < 5:
            note += "; n<5 in one group, p is approximate"
        rows.append([
            metric, scale, g1, str(res["n1"]), g2, str(res["n2"]),
            fnum(res["median1"], 4), fnum(res["median2"], 4),
            fnum(res["mean1"], 4), fnum(res["mean2"], 4),
            fnum(res["U1"], 1), fnum(res["U2"], 1), fnum(res["z"], 4),
            fp(res["p"]), fnum(res["effect_r"], 4), note,
        ])
        results.append((metric, res))
    return rows, results


def write_group_test(path, samples, args, g1, g2):
    header = [
        "metric", "scale", "group1", "n1", "group2", "n2",
        "median_group1", "median_group2", "mean_group1", "mean_group2",
        "U_group1", "U_group2", "z", "p_two_sided_normal",
        "effect_rank_biserial_group1_vs_group2", "note",
    ]
    rows, results = group_test_rows(samples, g1, g2)
    write_tsv(path, common_comments(
        args, ["Mann-Whitney U implemented in the standard library; U_group1 is the "
               "U statistic for %s" % g1,
               "positive rank-biserial effect means %s tends to be higher than %s"
               % (g1, g2)]), header, rows)
    return results


def write_matrix(path, samples, anon, args):
    header = (["sample", "group", "run", "counts_source"]
              + GROUP_ORDER + ["n_virus_groups_detected", CHIMP_COL]
              + [g + "_reads" for g in GROUP_ORDER]
              + [CHIMP_COL + "_reads"])
    rows = []
    for rec in samples:
        row = [anon[rec["sample"]], rec["group"], rec["run"], rec["counts_source"]]
        row += [str(rec["present"][g]) for g in GROUP_ORDER]
        row += [str(rec["n_groups"]), str(rec["chimp_present"])]
        row += [str(rec["group_reads"][g]) for g in GROUP_ORDER]
        row += [str(rec["chimp_reads"])]
        rows.append(row)
    write_tsv(path, common_comments(
        args, ["1 = present, 0 = absent; %s is a flagged column and is NOT counted "
               "in n_virus_groups_detected" % CHIMP_COL,
               "ANELLO here means human anelloviruses only"]), header, rows)


def write_pairs(path, samples, args):
    header = [
        "cohort", "n_samples", "virus_group_a", "virus_group_b",
        "n_a", "n_b", "n_both", "n_either", "n_neither",
        "jaccard", "expected_both_if_independent", "obs_over_exp",
    ]
    cohorts = [("ALL", samples)]
    for label in ["HIV", "HL", "TCL", "NA"]:
        sub = [r for r in samples if r["group"] == label]
        if len(sub) >= 2:
            cohorts.append((label, sub))
    rows = []
    for cohort, sub in cohorts:
        n = len(sub)
        for a, b in itertools.combinations(GROUP_ORDER, 2):
            na = sum(r["present"][a] for r in sub)
            nb = sum(r["present"][b] for r in sub)
            both = sum(1 for r in sub if r["present"][a] and r["present"][b])
            either = sum(1 for r in sub if r["present"][a] or r["present"][b])
            neither = n - either
            jac = (float(both) / float(either)) if either else None
            exp = (float(na) * float(nb) / float(n)) if n else None
            ratio = None
            if exp:
                ratio = float(both) / exp
            rows.append([
                cohort, str(n), a, b, str(na), str(nb), str(both), str(either),
                str(neither), fnum(jac, 4), fnum(exp, 3), fnum(ratio, 3),
            ])
    write_tsv(path, common_comments(
        args, ["jaccard = n_both / n_either; NA when neither group is present",
               "expected_both_if_independent = n_a * n_b / n_samples"]),
        header, rows)


def write_ref_audit(path, assign, refs, args, refmap_paths):
    header = ["reference_id", "map_category", "virus_group", "chimp_flagged",
              "reason", "description"]
    rows = []
    for rid in sorted(assign):
        grp, chimp, reason = assign[rid]
        if grp is None:
            continue
        info = refs.get(rid, {})
        rows.append([rid, info.get("category", "NA"), grp,
                     "1" if chimp else "0", reason,
                     info.get("description", "")[:180]])
    write_tsv(path, common_comments(
        args, ["reference maps used: %s" % ("; ".join(refmap_paths) or "none"),
               "audit of the keyword / accession classification; no sample data"]),
        header, rows)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    ap = argparse.ArgumentParser(
        description="Anellovirus burden and coinfection structure from panel runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--runs", nargs="+", default=DEFAULT_RUNS,
                    help="run directories to pool (each needs results/*.idxstats.tsv)")
    ap.add_argument("--refmap", default=DEFAULT_REFMAP,
                    help="panel reference map CSV used when the run manifest has none")
    ap.add_argument("--base-refmap", default=DEFAULT_BASE_REFMAP,
                    help="fallback reference map CSV")
    ap.add_argument("--no-manifest-refmap", dest="use_manifest_refmap",
                    action="store_false", default=True,
                    help="do not read the reference map path from run_manifest.tsv")
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR, help="output directory")
    ap.add_argument("--prefix", default="a7_virome",
                    help="filename prefix for every output file")
    ap.add_argument("--anello-accessions", default=DEFAULT_ANELLO_ACC_FILE,
                    help="optional file of extra anellovirus accessions / reference ids")
    ap.add_argument("--chimp-accessions", default=DEFAULT_CHIMP_ACC,
                    help="comma-separated chimpanzee-isolate accessions to flag apart")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="reads needed to call a reference or virus group present")
    ap.add_argument("--counts", choices=["idxstats", "unique_best"], default="idxstats",
                    help="per-reference counts from idxstats, or recounted from BAMs")
    ap.add_argument("--samtools", default="samtools", help="samtools executable")
    ap.add_argument("--mapq", type=int, default=40,
                    help="MAPQ floor for --counts unique_best")
    ap.add_argument("--exclude-flags", default="0x904",
                    help="samtools view -F value for --counts unique_best")
    ap.add_argument("--norm-source",
                    choices=["idxstats_mapped", "idxstats_total", "filtered_categories"],
                    default="idxstats_mapped",
                    help="denominator for the per-million normalisation")
    ap.add_argument("--test-groups", default="HIV,HL",
                    help="the two group labels compared with Mann-Whitney U")
    ap.add_argument("--no-run-name-group", action="store_true",
                    help="do not fall back to the run directory name for TCL labelling")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    args._chimp_acc = [a.strip() for a in args.chimp_accessions.split(",") if a.strip()]
    extra = read_accession_list(args.anello_accessions)
    if extra is None:
        warn("anellovirus accession list", args.anello_accessions)
        extra = []
    args._anello_extra = extra

    test_groups = [g.strip() for g in args.test_groups.split(",") if g.strip()]
    if len(test_groups) != 2:
        print("WARN: --test-groups needs exactly two labels, got %r; using HIV,HL"
              % args.test_groups)
        test_groups = ["HIV", "HL"]

    refs_cache = {}
    samtools_state = {"broken": False}
    raw_records = []
    assign_all = {}
    refs_all = {}
    refmap_paths = []
    for run in args.runs:
        recs, assign, refmap_path = collect_run(run, args, refs_cache, samtools_state)
        raw_records.extend(recs)
        if assign:
            assign_all.update(assign)
        if refmap_path:
            if refmap_path not in refmap_paths:
                refmap_paths.append(refmap_path)
            refs_all.update(refs_cache[refmap_path][0])

    if not raw_records:
        print("WARN: no usable samples were found in any run; nothing to do")
        return 0

    try:
        os.makedirs(args.outdir)
    except OSError:
        if not os.path.isdir(args.outdir):
            warn("writable output directory", args.outdir)
            return 0

    samples = [summarise_sample(r, args.min_reads, args.norm_source)
               for r in raw_records]
    samples.sort(key=lambda r: (r["sample"], r["run"]))
    anon = anonymise([r["sample"] for r in samples])

    n_anello_refs = sum(1 for v in assign_all.values() if v[0] == "ANELLO" and not v[1])
    n_chimp_refs = sum(1 for v in assign_all.values() if v[0] == "ANELLO" and v[1])
    if n_anello_refs == 0:
        print("WARN: no human anellovirus reference was identified in the reference "
              "map; supply --anello-accessions to add them explicitly")

    key_path = os.path.join(args.outdir, out_name(args.prefix, "sample_key.tsv"))
    burden_path = os.path.join(args.outdir,
                               out_name(args.prefix, "anellovirus_burden.tsv"))
    test_path = os.path.join(args.outdir,
                             out_name(args.prefix, "anellovirus_group_test.tsv"))
    matrix_path = os.path.join(args.outdir,
                               out_name(args.prefix, "coinfection_matrix.tsv"))
    pairs_path = os.path.join(args.outdir,
                              out_name(args.prefix, "coinfection_pairs.tsv"))
    audit_path = os.path.join(args.outdir,
                              out_name(args.prefix, "virus_group_refs.tsv"))

    write_sample_key(key_path, samples, anon)
    write_burden(burden_path, samples, anon, args)
    results = write_group_test(test_path, samples, args, test_groups[0], test_groups[1])
    write_matrix(matrix_path, samples, anon, args)
    write_pairs(pairs_path, samples, args)
    write_ref_audit(audit_path, assign_all, refs_all, args, refmap_paths)

    # ---------------------------- stdout summary --------------------------- #
    print("")
    print("a7_virome_structure  %s" % TODAY)
    print("runs           : %d requested, %d contributed samples"
          % (len(args.runs), len(set(r["run"] for r in samples))))
    print("samples        : %d (anonymised %s..%s)"
          % (len(samples), anon[samples[0]["sample"]], anon[samples[-1]["sample"]]))
    counts_by_group = {}
    for rec in samples:
        counts_by_group[rec["group"]] = counts_by_group.get(rec["group"], 0) + 1
    print("groups         : %s"
          % ", ".join("%s=%d" % (g, counts_by_group[g]) for g in sorted(counts_by_group)))
    print("counts source  : %s | presence threshold %d reads | norm %s"
          % (samples[0]["counts_source"], args.min_reads, args.norm_source))
    print("anello refs    : %d human, %d chimpanzee-flagged" % (n_anello_refs, n_chimp_refs))
    zero_mapped = sum(1 for r in samples if not r["total_mapped"])
    if zero_mapped:
        print("NOTE           : %d sample(s) have 0 mapped reads in idxstats; they are "
              "kept with zero counts and get RPM = NA" % zero_mapped)
    no_denom = sum(1 for r in samples if r["denom"] is None)
    if no_denom:
        print("NOTE           : %d sample(s) have no usable %s denominator and are "
              "dropped from the RPM test only" % (no_denom, args.norm_source))

    det = sum(1 for r in samples if r["anello_richness"] > 0)
    print("")
    print("(1) anellovirus burden")
    print("    detected in %d/%d samples; max richness %d, max RPM %s"
          % (det, len(samples),
             max([r["anello_richness"] for r in samples] or [0]),
             fnum(max([r["anello_rpm"] for r in samples if r["anello_rpm"] is not None]
                      or [0.0]), 2)))
    for metric, res in results:
        if res is None:
            print("    %-26s not tested (a group was empty)" % metric)
            continue
        print("    %-26s %s median %s vs %s median %s | U=%s p=%s"
              % (metric, test_groups[0], fnum(res["median1"], 3),
                 test_groups[1], fnum(res["median2"], 3),
                 fnum(res["U1"], 1), fp(res["p"])))

    print("")
    print("(2) coinfection structure")
    for grp in GROUP_ORDER:
        npos = sum(r["present"][grp] for r in samples)
        print("    %-13s present in %3d/%d samples" % (grp, npos, len(samples)))
    nchimp = sum(r["chimp_present"] for r in samples)
    print("    %-13s present in %3d/%d samples (flagged, not counted as human)"
          % (CHIMP_COL, nchimp, len(samples)))
    top = []
    for a, b in itertools.combinations(GROUP_ORDER, 2):
        both = sum(1 for r in samples if r["present"][a] and r["present"][b])
        either = sum(1 for r in samples if r["present"][a] or r["present"][b])
        if both:
            top.append((both, float(both) / float(either), a, b))
    top.sort(key=lambda t: (-t[0], -t[1]))
    if top:
        print("    strongest co-occurrences (whole cohort):")
        for both, jac, a, b in top[:5]:
            print("      %-13s + %-13s n_both=%3d jaccard=%.3f" % (a, b, both, jac))
    else:
        print("    no virus-group pair co-occurs above the read threshold")

    print("")
    print("wrote:")
    for path in (key_path, burden_path, test_path, matrix_path, pairs_path, audit_path):
        print("  %s" % path)
    print("REMINDER: %s contains real sample identifiers - do not commit or email it."
          % os.path.basename(key_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
