#!/usr/bin/env python3
"""a13_anello_chimp_control.py -- is the chimpanzee anellovirus signal derivative
of a real human one, or evidence of a shared artefact?

THE QUESTION
  a10/a12 used the three chimpanzee anellovirus isolates (NC_014069.1,
  NC_014077.1, NC_014480.2) as a negative control: no human carries chimp TTV,
  so reads landing there must be artefact. The observation was that 729 chimp
  reads survive a12's full filter ladder in the HIV+ group against 0 in HL, and
  that per reference the chimp set accumulates about as many qualifying pairs as
  the human set (3.3 vs 3.6). The report reads that as a reason not to trust
  species-level statements.

  That reasoning has a gap. Chimp anelloviruses are congeners of the human ones
  and share the conserved UTR, so if a patient genuinely carries human
  anellovirus, reads cross-map onto the chimp references BY CONSTRUCTION. Chimp
  reads present in HIV+ and absent in HL is exactly what you would predict IF
  the human signal is real. The chimp signal may be a downstream consequence of
  true positives rather than evidence against them.

  This module runs the tests that separate those two readings. It does NOT
  re-align anything; it re-reads the same counts a7 reads.

THE THREE TESTS

  1 SUBSET TEST -- necessary, but NOT sufficient.
      Are the chimp-positive samples a SUBSET of the human-positive samples?
        chimp+ / human+   chimp reads only where human signal exists
        chimp+ / human-   chimp reads where there is NO human signal
                          -> a chimp-preferential source that is neither human
                             anellovirus nor cross-mapping from it: a problem
      Reported as a 2x2 with a Fisher exact p. The cell that matters is
      chimp+/human-; a p value alone does not answer the question.

      READ THIS TEST NARROWLY. An empty chimp+/human- cell rules out a source
      that hits the chimp references PREFERENTIALLY. It does NOT distinguish
      "derivative cross-mapping from real human virus" from "a shared source
      spraying human and chimp references together", because the shared-source
      scenario also produces no orphan samples. Verified on synthetic data:
      both scenarios pass test 1. Test 2 is what separates them.

  2 BEST-REFERENCE TEST.
      "3.3 vs 3.6 qualifying pairs per reference" compares the MEAN across
      references, and the mean is the wrong statistic: most of the 20 human
      references are not the patient's strain either, so they are cross-mapping
      sinks too and they drag the human mean down toward the chimp rate. A real
      infection CONCENTRATES on its best-matching reference; artefact spreads
      evenly. So compare the best reference, not the average one.

      The naive version of that comparison is biased: the maximum of 20 human
      references beats the maximum of 3 chimp references even under a pure
      artefact null, purely because 20 > 3. This module corrects for it exactly,
      comparing best_chimp against the EXPECTED MAXIMUM OF A RANDOM 3-SUBSET of
      the human references (closed form over all C(n,3) subsets, no sampling).
      Paired per sample, then an exact binomial sign test.

  3 READ QUALITY VS THE NEGATIVE CONTROL (optional, --a11-groups).
      Surfaces a11's own refset_contrast rows: human anellovirus references
      against the chimpanzee control, PAIRED Wilcoxon signed-rank within each
      sample carrying reads on both sets.

      This is diagnostic, in one direction. Under DERIVATIVE cross-mapping a
      genuine read maps best to its human reference and reaches a chimp
      reference only as a worse, more ambiguous alignment, so the human set
      should win on MAPQ, clipping and alignment entropy. Under a SHARED
      artefact the two sets should look alike. On the real data (2026-07-28)
      a11 gives median MAPQ 57 vs 40 (p = 4.4e-04), reads clipped 0.78 vs 1.00
      (p = 1.0e-04) and alignment entropy 5.13 vs 4.50 (p = 2.0e-05) -- all
      three favour the human references.

      MATE PAIRING IS NOT DIAGNOSTIC and is not scored: when a real pair
      cross-maps onto a chimp reference both mates move together, so chimp
      reads pair as tidily as human ones either way (a11: 0.753 vs 0.792,
      p = 0.71, as predicted).

      HISTORY, because the mistake is instructive. a13 originally recomputed
      this from a11_forensics_by_pair.tsv with an UNPAIRED Mann-Whitney over
      (sample, reference) rows. With 20 human references against 3 chimpanzee
      ones, the human distribution is dominated by marginal references drawn
      from different samples than the chimp rows, and that manufactured a null
      (MAPQ 42.25 vs 43.0, p = 0.26) out of a real 17-point paired gap. Do not
      recompute this contrast here; read a11's paired version.

WHAT THIS CAN AND CANNOT SHOW
  CAN: show whether chimp signal ever occurs without human signal; whether the
       best human reference outperforms the best chimp reference once the
       set-size advantage is removed; whether human-reference reads pair better.
  CANNOT: prove the reads are anellovirus. A shared non-anellovirus source that
       resembles the conserved UTR would inflate human and chimp references
       alike and would pass tests 1 and 2 if it tracked HIV status. Test 3 is
       the one that speaks to that, and competitive realignment against
       hg38 + CHM13 + microbial + adapter/probe references is still the
       confirmatory experiment.

WHAT IT WRITES (tab separated, pure ASCII, into --outdir)
  <prefix>_by_sample.tsv    one row per sample: human and chimp read totals,
                            positivity flags, best human / best chimp reference
                            counts, and the set-size-corrected human value
  <prefix>_tests.tsv        one row per test with the statistic and p
  <prefix>_sample_key.tsv   real -> anonymous mapping. THE ONLY file that
                            carries real sample identifiers. Never commit or
                            email it.
  A verdict for each test goes to stdout.

EXAMPLE
  python3 a13_anello_chimp_control.py \
      --runs /path/to/runs/wgs_hiv_hl_hg38_shuyu_masked_panel_refixed_primary_only \
      --counts unique_best --mapq 40 \
      --a11-pairs /path/to/suite_out/a11_forensics_by_pair.tsv \
      --outdir /path/to/suite_out --prefix a13_chimp
"""
from __future__ import annotations

import csv
import math
import os
import sys

import a7_virome_structure as a7
from a12_anello_utr_exclusion import fisher_exact_2x2, mann_whitney_u


# --------------------------------------------------------------------------- #
# statistics not already in a7 / a12
# --------------------------------------------------------------------------- #
def expected_max_of_k(values, k):
    """Exact E[max] of a uniformly random k-subset of `values`.

    Removes the set-size advantage when comparing the best of n human
    references against the best of k chimp references. Sorting descending, the
    number of k-subsets whose maximum sits at position i is C(n-1-i, k-1), so
    E[max] = sum_i v[i] * C(n-1-i, k-1) / C(n, k). Ties are broken by position,
    which does not change the value of the maximum.
    """
    vals = sorted(values, reverse=True)
    n = len(vals)
    if n == 0 or k <= 0:
        return None
    if k >= n:
        return float(vals[0])
    total = math.comb(n, k)
    acc = 0.0
    for i, v in enumerate(vals):
        ways = math.comb(n - 1 - i, k - 1) if (n - 1 - i) >= (k - 1) else 0
        if ways:
            acc += float(v) * ways
    return acc / total


def sign_test(wins, losses):
    """Exact two-sided binomial sign test at p = 0.5. Ties are excluded."""
    n = wins + losses
    if n == 0:
        return None
    obs = min(wins, losses)
    tail = sum(math.comb(n, i) for i in range(0, obs + 1)) / float(2 ** n)
    return min(1.0, 2.0 * tail)


def fnum(value, digits=4):
    if value is None:
        return "NA"
    return ("%." + str(digits) + "f") % value


# --------------------------------------------------------------------------- #
# test 3: reuse a11's per-pair forensics
# --------------------------------------------------------------------------- #
# Metrics where a difference is interpretable, and the direction that favours the
# human references being a genuine source rather than a cross-mapping sink.
# +1 = higher is better for human, -1 = lower is better for human.
A11_DIRECTION = {
    "median_mapq": +1,
    "aln_median_entropy3": +1,
    "frac_reads_clipped": -1,
    "frac_mapq0": -1,
    "frac_as_eq_xs": -1,
}


def read_a11_contrast(path):
    """-> list of (metric, human, chimp, p) from a11_forensics_by_group.tsv.

    a11 already computes this contrast correctly: a PAIRED Wilcoxon signed-rank
    within each sample carrying reads on both reference sets. a13 must not
    recompute it from the per-pair table -- an unpaired comparison across a
    design with 20 human references against 3 chimpanzee ones is dominated by
    marginal human references drawn from different samples than the chimp rows,
    which manufactures a null. (That bug was in a13 until 2026-07-28.)
    """
    if not path or not os.path.isfile(path):
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rows = [r for r in csv.DictReader(
            (ln for ln in fh if not ln.startswith("#")), delimiter="\t")]
    if not rows or "row_type" not in rows[0]:
        sys.stderr.write("WARN: %s is not an a11 by_group table; skipping test 3\n"
                         % os.path.basename(path))
        return []
    out = []
    for r in rows:
        if (r.get("row_type") or "").strip() != "refset_contrast":
            continue
        metric = (r.get("label") or "").strip()
        try:
            h = float(r["value_group1"])
            c = float(r["value_group2"])
            p = float(r["p_value"])
        except (KeyError, TypeError, ValueError):
            continue
        out.append((metric, h, c, p))
    return out


# --------------------------------------------------------------------------- #
def build_parser():
    ap = a7.build_parser()
    ap.set_defaults(prefix="a13_chimp")
    ap.add_argument("--a11-groups", "--a11-pairs", dest="a11_groups", default=None,
                    help="a11_forensics_by_group.tsv, enabling the read-quality "
                         "contrast (test 3). Optional. Pass the BY_GROUP file: a11 "
                         "computes the human-vs-chimp contrast paired within sample, "
                         "which is the correct test.")
    ap.description = "a13: chimpanzee-control tests for the anellovirus signal"
    return ap


def main():
    args = build_parser().parse_args()

    # a7.collect_run expects these two derived fields, which a7.main() normally
    # sets up before calling it.
    args._chimp_acc = [a.strip() for a in args.chimp_accessions.split(",") if a.strip()]
    extra = a7.read_accession_list(args.anello_accessions)
    if extra is None:
        a7.warn("anellovirus accession list", args.anello_accessions)
        extra = []
    args._anello_extra = extra

    groups = [g.strip() for g in args.test_groups.split(",") if g.strip()]
    if len(groups) != 2:
        print("WARN: --test-groups needs exactly two labels, got %r; using HIV,HL"
              % args.test_groups)
        groups = ["HIV", "HL"]
    g1, g2 = groups

    refs_cache, samtools_state = {}, {"broken": False}
    records = []
    for run in args.runs:
        recs, _assign, _refmap = a7.collect_run(run, args, refs_cache, samtools_state)
        records.extend(recs)
    if not records:
        sys.stderr.write("ERROR: no samples collected from --runs\n")
        return 1

    records = [r for r in records if r["group"] in (g1, g2)]
    if not records:
        sys.stderr.write("ERROR: no samples in groups %s / %s\n" % (g1, g2))
        return 1

    anon = a7.anonymise([r["sample"] for r in records])
    minr = args.min_reads

    rows = []
    for rec in sorted(records, key=lambda r: anon[r["sample"]]):
        counts = rec["counts"]
        human_refs = list(rec["group_refs"]["ANELLO"])
        chimp_refs = list(rec["chimp_refs"])
        hv = [counts.get(r, 0) for r in human_refs]
        cv = [counts.get(r, 0) for r in chimp_refs]

        best_h = max(hv) if hv else 0
        best_c = max(cv) if cv else 0
        # like-for-like: best of a random len(chimp_refs)-subset of the human set
        exp_h_k = expected_max_of_k(hv, len(chimp_refs)) if hv and chimp_refs else None

        rows.append({
            "sample_anon": anon[rec["sample"]],
            "run": rec["run"],
            "group": rec["group"],
            "counts_source": rec["counts_source"],
            "n_human_refs": len(human_refs),
            "n_chimp_refs": len(chimp_refs),
            "human_reads": sum(hv),
            "chimp_reads": sum(cv),
            "human_pos": 1 if any(v >= minr for v in hv) else 0,
            "chimp_pos": 1 if any(v >= minr for v in cv) else 0,
            "best_human_ref_reads": best_h,
            "best_chimp_ref_reads": best_c,
            "expected_best_human_k_subset": fnum(exp_h_k, 3),
        })

    # ---------------------------------------------------------------- test 1
    a = sum(1 for r in rows if r["chimp_pos"] and r["human_pos"])
    b = sum(1 for r in rows if r["chimp_pos"] and not r["human_pos"])
    c = sum(1 for r in rows if not r["chimp_pos"] and r["human_pos"])
    d = sum(1 for r in rows if not r["chimp_pos"] and not r["human_pos"])
    fisher = fisher_exact_2x2(a, b, c, d)
    p_subset = fisher["p"] if fisher else None
    orphan_reads = sum(r["chimp_reads"] for r in rows
                       if r["chimp_pos"] and not r["human_pos"])
    subset_holds = (b == 0)

    # ---------------------------------------------------------------- test 2
    wins = losses = ties = 0
    ratios = []
    for r in rows:
        exp_h = r["expected_best_human_k_subset"]
        if exp_h == "NA" or (r["best_chimp_ref_reads"] == 0 and r["best_human_ref_reads"] == 0):
            continue
        eh = float(exp_h)
        bc = float(r["best_chimp_ref_reads"])
        if eh > bc:
            wins += 1
        elif bc > eh:
            losses += 1
        else:
            ties += 1
        if bc > 0:
            ratios.append(eh / bc)
    p_sign = sign_test(wins, losses)
    med_ratio = a7.median(ratios) if ratios else None

    # per-group best-reference contrast, uncorrected but informative
    bh = {g: [r["best_human_ref_reads"] for r in rows if r["group"] == g] for g in (g1, g2)}
    mw_best = mann_whitney_u(bh[g1], bh[g2]) if bh[g1] and bh[g2] else None

    # ---------------------------------------------------------------- test 3
    contrast = read_a11_contrast(args.a11_groups)
    favours_human = [(m, h, c, p) for (m, h, c, p) in contrast
                     if p < 0.05 and m in A11_DIRECTION
                     and (h - c) * A11_DIRECTION[m] > 0]
    favours_chimp = [(m, h, c, p) for (m, h, c, p) in contrast
                     if p < 0.05 and m in A11_DIRECTION
                     and (h - c) * A11_DIRECTION[m] < 0]

    # ---------------------------------------------------------------- output
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    prefix = args.prefix

    def out(name):
        return os.path.join(args.outdir, "%s_%s" % (prefix, name))

    comments = [
        "a13_anello_chimp_control.py",
        "no real sample identifiers in this file; ids are anonymous "
        "(mapping in %s_sample_key.tsv)" % prefix,
        "counts source: %s ; presence threshold: reads >= %d per reference"
        % (rows[0]["counts_source"] if rows else "NA", minr),
    ]
    # a7.write_tsv takes rows as sequences of VALUES, not dicts
    sample_header = list(rows[0])
    a7.write_tsv(out("by_sample.tsv"), comments, sample_header,
                 [[r[k] for k in sample_header] for r in rows])

    trows = [
        {"test": "1_subset", "statistic": "chimp+/human- samples", "value": b,
         "p": fp(p_subset),
         "detail": "2x2 chimp+/human+=%d chimp+/human-=%d chimp-/human+=%d "
                   "chimp-/human-=%d ; %d chimp reads in human-negative samples"
                   % (a, b, c, d, orphan_reads)},
        {"test": "2_best_reference", "statistic": "samples where E[best human 3-subset] > best chimp",
         "value": "%d win / %d loss / %d tie" % (wins, losses, ties),
         "p": fp(p_sign),
         "detail": "median ratio E[best human k-subset] / best chimp = %s ; "
                   "set-size corrected, exact over all C(n,k) subsets"
                   % fnum(med_ratio, 3)},
    ]
    if mw_best:
        trows.append({
            "test": "2b_best_human_by_group",
            "statistic": "median best human ref reads %s vs %s" % (g1, g2),
            "value": "%s vs %s" % (fnum(mw_best["median1"], 1), fnum(mw_best["median2"], 1)),
            "p": fp(mw_best["p"]),
            "detail": "Mann-Whitney, rank-biserial %s" % fnum(mw_best["effect_r"], 3)})
    for m, h, c, p in contrast:
        direction = A11_DIRECTION.get(m)
        if direction is None:
            verdict = "not directional"
        elif p >= 0.05:
            verdict = "indistinguishable from the negative control"
        else:
            verdict = "favours human" if (h - c) * direction > 0 else "favours chimp"
        trows.append({
            "test": "3_read_quality_a11",
            "statistic": "%s human vs chimp references" % m,
            "value": "%s vs %s" % (fnum(h, 4), fnum(c, 4)),
            "p": fp(p),
            "detail": "a11 paired Wilcoxon signed-rank within sample; %s" % verdict})
    test_header = ["test", "statistic", "value", "p", "detail"]
    a7.write_tsv(out("tests.tsv"), comments, test_header,
                 [[str(t[k]) for k in test_header] for t in trows])

    key_path = out("sample_key.tsv")
    a7.write_sample_key(key_path, records, anon)

    # ---------------------------------------------------------------- stdout
    print("")
    print("a13 chimpanzee-control tests   (%d samples: %d %s, %d %s)"
          % (len(rows),
             sum(1 for r in rows if r["group"] == g1), g1,
             sum(1 for r in rows if r["group"] == g2), g2))
    print("counts source: %s" % (rows[0]["counts_source"] if rows else "NA"))
    print("")
    print("TEST 1  subset test -- are chimp+ samples a subset of human+ samples?")
    print("          chimp+ human+ : %3d" % a)
    print("          chimp+ human- : %3d   <-- the cell that matters" % b)
    print("          chimp- human+ : %3d" % c)
    print("          chimp- human- : %3d" % d)
    print("          Fisher p = %s ; %d chimp reads sit in human-negative samples"
          % (fp(p_subset), orphan_reads))
    if subset_holds:
        print("          VERDICT: no chimp signal without human signal. This rules out a")
        print("                   source that hits the chimp references preferentially.")
        print("                   It does NOT rule out a shared source hitting both -")
        print("                   that scenario passes this test too. See test 2.")
    else:
        print("          VERDICT: %d sample(s) carry chimp reads with NO human signal." % b)
        print("                   A chimp-preferential source is contributing. Something")
        print("                   other than human anellovirus is generating these reads.")
    print("")
    print("TEST 2  best-reference test -- corrected for the 20-vs-3 set-size advantage")
    print("          E[best of a random %d-subset of human refs] vs best chimp ref"
          % (rows[0]["n_chimp_refs"] if rows else 3))
    print("          human higher in %d samples, chimp higher in %d, tied in %d"
          % (wins, losses, ties))
    print("          exact sign test p = %s ; median ratio = %s"
          % (fp(p_sign), fnum(med_ratio, 3)))
    if p_sign is not None and p_sign < 0.05 and wins > losses:
        print("          VERDICT: the human references still win once set size is")
        print("                   equalised. The '3.3 vs 3.6 per reference' equivalence")
        print("                   was an artefact of averaging over 20 references.")
    else:
        print("          VERDICT: no advantage for the human references once set size is")
        print("                   equalised. This supports the cautious reading.")
    print("")
    if contrast:
        print("TEST 3  read quality vs the negative control (a11, paired within sample)")
        for m, h, c, p in contrast:
            mark = ""
            if m in A11_DIRECTION and p < 0.05:
                mark = ("  <-- favours human" if (h - c) * A11_DIRECTION[m] > 0
                        else "  <-- favours chimp")
            print("          %-26s human %-9s chimp %-9s p = %-11s%s"
                  % (m, fnum(h, 4), fnum(c, 4), fp(p), mark))
        print("")
        print("          %d directional metric(s) favour the human references, %d favour"
              % (len(favours_human), len(favours_chimp)))
        print("          the chimpanzee control. Under DERIVATIVE cross-mapping a genuine")
        print("          read maps best to its human reference and reaches a chimp")
        print("          reference only as a worse, more ambiguous alignment, so a human")
        print("          advantage on MAPQ, clipping and entropy is predicted. Under a")
        print("          SHARED artefact the two sets should look alike. Mate-pairing is")
        print("          NOT diagnostic either way: a cross-mapping pair moves both mates.")
    else:
        print("TEST 3  skipped -- pass --a11-groups a11_forensics_by_group.tsv to run it.")
    print("")

    t2_ok = (p_sign is not None and p_sign < 0.05 and wins > losses)
    print("OVERALL")
    if not subset_holds:
        print("  A chimp-preferential source is present (test 1). Neither the species")
        print("  nor the presence claim is safe until that source is identified.")
    elif t2_ok:
        print("  Test 1 excludes a chimp-preferential source and test 2 shows the human")
        print("  references still win once the 20-vs-3 set-size advantage is removed.")
        print("  Together these are consistent with the chimp signal being DERIVATIVE")
        print("  cross-mapping from real human anellovirus, not evidence against it.")
        print("  On this reading, presence-only is more conservative than the data")
        print("  require and burden is also defensible. Species assignment stays out:")
        print("  nothing here identifies WHICH anellovirus, only that one is present.")
    else:
        print("  Test 1 excludes a chimp-preferential source, but test 2 finds no")
        print("  advantage for the human references once set size is equalised. A")
        print("  shared source spraying both reference sets is not excluded. Keep the")
        print("  cautious position: presence only, no burden, no species.")
    print("  Neither test can prove the reads are anellovirus. Competitive realignment")
    print("  against hg38 + CHM13 + microbial + adapter/probe references remains the")
    print("  confirmatory experiment.")
    print("")
    print("wrote:")
    for nm in ("by_sample.tsv", "tests.tsv", "sample_key.tsv"):
        print("  %s" % out(nm))
    print("REMINDER: %s contains real sample identifiers - do not commit or email it."
          % os.path.basename(key_path))
    return 0


def fp(value):
    return a7.fp(value)


if __name__ == "__main__":
    sys.exit(main())
