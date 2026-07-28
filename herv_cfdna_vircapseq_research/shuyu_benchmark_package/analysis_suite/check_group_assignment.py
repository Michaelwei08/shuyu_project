#!/usr/bin/env python3
"""check_group_assignment.py -- is the HIV/HL split the suite reports actually
the split the manifest defines?

THE QUESTION
  Every conclusion that compares HIV+ against HL controls -- the anellovirus
  presence test (21/37 vs 3/23), the EBV type 2 concentration, the WGS HIV
  depth-limit result -- rests on a group label that no module ever reads from
  the manifest. All twelve modules derive it by parsing the sample FILENAME.

  A differential test of those twelve derivations (2026-07-27) found they agree
  with each other on every case. That rules out divergence between modules. It
  does NOT rule out all twelve being wrong the same way, because they share one
  algorithm:

      m = re.search(r"(?:^|_)(HIV|HL)[0-9]", sample)   # explicit label
      if m: return m.group(1)
      low = sample.lower()
      if "_hiv" in low: return "HIV"                   # FALLBACK
      if "_hl"  in low: return "HL"

  The fallback is the hazard. WGS sample names carry the run prefix
  "wgs_60samples_hiv_hl_", which contains BOTH "_hiv" and "_hl", and "_hiv" is
  tested first. So any WGS sample whose name lacks an explicit HIV<digit> or
  HL<digit> label is silently labelled HIV, and can never be labelled HL. One
  such name inflates the HIV arm and deflates nothing -- it just disappears
  from the control group.

WHAT THIS CHECKS
  1. FALLBACK    for each sample, did the explicit case-sensitive regex match,
                 or did the label come from the prefix fallback? Any fallback
                 hit on a WGS sample is a finding, not a warning.
  2. MANIFEST    derived label vs the manifest's own group column, per sample.
  3. COUNTS      derived totals against --expect (default 37 HIV / 23 HL).

  Exit status is 0 only if there are no fallback hits, no manifest mismatches,
  and the counts match --expect. Anything else exits 1, so it can gate a run.

INPUTS
  --names-from   one sample name per line, OR a *_sample_key.tsv written by
                 a7/a10/a11/a12 (the real->anonymous map; the real name is
                 taken from the column named by --name-column, default
                 "sample"). Alternatively use --bam-dir.
  --bam-dir      derive names from *.bam in this directory, stripping the same
                 pipeline suffixes as the suite (commit a420e72).
  --manifest     CSV or TSV with a sample-name column (--manifest-name-column,
                 default "sample_id") and a group column
                 (--manifest-group-column, default "expected_group").
                 Optional: without it, checks 1 and 3 still run.

OUTPUT
  A TSV to --out (default stdout) with one row per sample, plus a summary to
  stderr. The TSV carries REAL sample names -- treat it like the sample_key
  files: never commit it, never email it.

EXAMPLE (cluster)
  python3 check_group_assignment.py \
      --names-from /drive3/cpwei/tmp/suite_out/a12_anello_utr_sample_key.tsv \
      --manifest   /path/to/wgs_manifest.csv \
      --out        /drive3/cpwei/tmp/suite_out/group_assignment_check.tsv
"""
import argparse
import csv
import os
import re
import sys

EXPLICIT = re.compile(r"(?:^|_)(HIV|HL)[0-9]")

# Same suffix list the suite strips when deriving a sample id from a BAM path.
BAM_SUFFIXES = (
    ".retrovirus", ".retro", ".markdup", ".dedup", ".sorted", ".filtered",
)


def strip_bam_suffixes(stem):
    changed = True
    while changed:
        changed = False
        for suf in BAM_SUFFIXES:
            if stem.endswith(suf):
                stem = stem[: -len(suf)]
                changed = True
    return stem


def derive(sample, run_name="", use_run_name=True):
    """The suite's group_of, verbatim, plus how the label was reached."""
    m = EXPLICIT.search(sample or "")
    if m:
        return m.group(1), "explicit"
    low = (sample or "").lower()
    if "_hiv" in low:
        return "HIV", "fallback_sample_hiv"
    if "_hl" in low:
        return "HL", "fallback_sample_hl"
    if "targeted_htlv" in low or "tcl" in low:
        return "TCL", "fallback_sample_tcl"
    if use_run_name:
        rlow = (run_name or "").lower()
        if "targeted_htlv" in rlow or "tcl" in rlow:
            return "TCL", "fallback_run_name"
    return "NA", "unresolved"


def sniff_read(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        head = fh.read(8192)
        fh.seek(0)
        delim = "\t" if head.count("\t") >= head.count(",") else ","
        return list(csv.DictReader(fh, delimiter=delim))


def load_names(args):
    if args.bam_dir:
        names = []
        for entry in sorted(os.listdir(args.bam_dir)):
            if entry.endswith(".bam"):
                names.append(strip_bam_suffixes(entry[: -len(".bam")]))
        return names
    with open(args.names_from, newline="", encoding="utf-8-sig") as fh:
        first = fh.readline()
    if "\t" in first or "," in first:
        rows = sniff_read(args.names_from)
        if rows and args.name_column in rows[0]:
            return [r[args.name_column] for r in rows if r.get(args.name_column)]
        if rows:
            key = list(rows[0])[0]
            sys.stderr.write(
                "note: column %r not found in %s, using first column %r\n"
                % (args.name_column, args.names_from, key))
            return [r[key] for r in rows if r.get(key)]
        return []
    with open(args.names_from, encoding="utf-8-sig") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--names-from")
    src.add_argument("--bam-dir")
    ap.add_argument("--name-column", default="sample")
    ap.add_argument("--manifest")
    ap.add_argument("--manifest-name-column", default="sample_id")
    ap.add_argument("--manifest-group-column", default="expected_group")
    ap.add_argument("--run-name", default="wgs_60samples_hiv_hl")
    ap.add_argument("--expect", default="HIV=37,HL=23",
                    help='expected counts, e.g. "HIV=37,HL=23" or "" to skip')
    ap.add_argument("--out", default="-")
    args = ap.parse_args()

    names = load_names(args)
    if not names:
        sys.stderr.write("ERROR: no sample names loaded\n")
        return 1

    manifest = {}
    if args.manifest:
        for row in sniff_read(args.manifest):
            key = (row.get(args.manifest_name_column) or "").strip()
            if key:
                manifest[key] = (row.get(args.manifest_group_column) or "").strip()
        if not manifest:
            sys.stderr.write(
                "ERROR: manifest %s yielded no rows under column %r\n"
                % (args.manifest, args.manifest_name_column))
            return 1

    rows, counts = [], {}
    fallbacks, mismatches, unmatched = [], [], []
    for name in names:
        group, how = derive(name, args.run_name)
        counts[group] = counts.get(group, 0) + 1
        expected, verdict = "", "no_manifest"
        if manifest:
            if name in manifest:
                expected = manifest[name]
                verdict = "match" if expected == group else "MISMATCH"
            else:
                # tolerate manifest ids that are a prefix/suffix of the run name
                hits = [k for k in manifest if k and (k in name or name in k)]
                if len(hits) == 1:
                    expected = manifest[hits[0]]
                    verdict = "match" if expected == group else "MISMATCH"
                else:
                    verdict = "NOT_IN_MANIFEST"
        if how != "explicit":
            fallbacks.append((name, group, how))
        if verdict == "MISMATCH":
            mismatches.append((name, group, expected))
        if verdict == "NOT_IN_MANIFEST":
            unmatched.append(name)
        rows.append({
            "sample": name,
            "derived_group": group,
            "how_derived": how,
            "manifest_group": expected,
            "verdict": verdict,
        })

    fh = sys.stdout if args.out == "-" else open(args.out, "w", newline="",
                                                 encoding="utf-8")
    try:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    finally:
        if fh is not sys.stdout:
            fh.close()

    err = sys.stderr.write
    err("\nsamples: %d\n" % len(names))
    err("derived counts: %s\n"
        % ", ".join("%s=%d" % kv for kv in sorted(counts.items())))

    ok = True
    if fallbacks:
        ok = False
        err("\nFALLBACK HITS (%d) - label came from the cohort prefix, not the\n"
            "sample label. On a WGS name this silently forces HIV:\n" % len(fallbacks))
        for name, group, how in fallbacks:
            err("  %-45s -> %-4s via %s\n" % (name, group, how))
    else:
        err("fallback hits: 0 (every label came from an explicit HIV<n>/HL<n>)\n")

    if manifest:
        if mismatches:
            ok = False
            err("\nMANIFEST MISMATCHES (%d):\n" % len(mismatches))
            for name, got, want in mismatches:
                err("  %-45s derived=%-4s manifest=%s\n" % (name, got, want))
        else:
            err("manifest mismatches: 0\n")
        if unmatched:
            ok = False
            err("\nNOT FOUND IN MANIFEST (%d):\n" % len(unmatched))
            for name in unmatched:
                err("  %s\n" % name)

    if args.expect:
        want = {}
        for part in args.expect.split(","):
            if "=" in part:
                k, v = part.split("=", 1)
                want[k.strip()] = int(v)
        bad = [(k, want[k], counts.get(k, 0)) for k in want
               if counts.get(k, 0) != want[k]]
        if bad:
            ok = False
            err("\nCOUNT MISMATCH:\n")
            for k, w_, g in bad:
                err("  %-4s expected %d, derived %d\n" % (k, w_, g))
        else:
            err("expected counts (%s): OK\n" % args.expect)

    err("\n%s\n" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
