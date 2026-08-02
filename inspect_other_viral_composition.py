#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from pathlib import Path

import pandas as pd


CORE = {"HIV1", "HIV2", "HTLV1", "HTLV2", "HERV", "LINE1", "HUMAN", "HG38"}


REF_COLS = [
    "reference_id", "ref_id", "reference", "ref", "contig", "rname",
    "target", "target_name", "sequence_id", "seq_id", "accession", "chrom"
]

SAMPLE_COLS = [
    "sample_id", "sample", "sample_name", "sample_key", "library", "bam", "run"
]

COUNT_COLS = [
    "filtered_count", "primary_count", "count", "counts", "read_count",
    "n_reads", "num_reads", "reads", "total_reads", "n"
]

CATEGORY_COLS = [
    "species", "category", "species_group", "reference_category",
    "group", "class", "type", "label"
]

NAME_COLS = [
    "name", "description", "reference_name", "organism", "virus", "taxon"
]


def norm_col(x):
    return str(x).strip().lower().replace("-", "_").replace(" ", "_")


def find_col(cols, candidates):
    norm_to_orig = {norm_col(c): c for c in cols}
    for cand in candidates:
        if cand in norm_to_orig:
            return norm_to_orig[cand]
    return None


def sniff_sep(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        line = f.readline()
    if line.count("\t") >= line.count(","):
        return "\t"
    return ","


def read_table(path, nrows=None):
    sep = sniff_sep(path)
    return pd.read_csv(path, sep=sep, nrows=nrows, dtype=str)


def load_refmap(path):
    if not path or not Path(path).exists():
        return {}

    df = read_table(Path(path))
    ref_col = find_col(df.columns, REF_COLS)
    if ref_col is None:
        print(f"[WARN] Could not find reference column in refmap: {path}")
        return {}

    cat_col = find_col(df.columns, CATEGORY_COLS)
    name_col = find_col(df.columns, NAME_COLS)

    out = {}
    for _, row in df.iterrows():
        rid = str(row.get(ref_col, "")).strip()
        if not rid or rid.lower() == "nan":
            continue
        out[rid] = {
            "refmap_category": str(row.get(cat_col, "")).strip() if cat_col else "",
            "refmap_name": str(row.get(name_col, "")).strip() if name_col else "",
        }
        # also map without version, e.g. NC_001802.1 -> NC_001802
        if "." in rid:
            out.setdefault(rid.split(".")[0], out[rid])
    return out


def infer_sample_from_path(path):
    parts = list(path.parts)
    # prefer parent directory, usually sample-specific
    for p in reversed(parts[:-1]):
        if re.search(r"(S\d+|P1|P2|HIV|HTLV|SRR|ERR|wgs|target)", p, re.I):
            return p
    return path.stem


def is_count_table(path):
    name = path.name.lower()

    if any(x in name for x in ["bam", "bai", "fasta", ".fa", ".fai", ".dict", ".log"]):
        return False
    if not any(name.endswith(x) for x in [".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz"]):
        return False
    if any(x in name for x in ["comparison", "reply_metrics", "cohort_summary"]):
        return False
    if not any(x in name for x in ["count", "summary", "reference", "viral", "retro", "filtered"]):
        return False

    try:
        df0 = read_table(path, nrows=3)
    except Exception:
        return False

    ref_col = find_col(df0.columns, REF_COLS)
    count_col = find_col(df0.columns, COUNT_COLS)

    return ref_col is not None and count_col is not None


def classify_other_viral(row, refmap):
    vals = []
    for c in CATEGORY_COLS:
        if c in row and pd.notna(row[c]):
            vals.append(str(row[c]).strip().upper())

    rid = str(row.get("reference_id", "")).strip()
    info = refmap.get(rid) or refmap.get(rid.split(".")[0]) or {}

    if info.get("refmap_category"):
        vals.append(info["refmap_category"].strip().upper())

    joined = " ".join(vals)

    if "OTHER_VIRAL" in joined:
        return True

    # Exclude known core buckets.
    for core in CORE:
        if core in joined:
            return False
        if core in rid.upper():
            return False

    # If it exists in the added viral panel / refmap but is not core, treat as OTHER_VIRAL.
    if info:
        return True

    return False


def normalize_table(path, cohort, refmap):
    df = read_table(path)
    df.columns = [norm_col(c) for c in df.columns]

    ref_col = find_col(df.columns, REF_COLS)
    sample_col = find_col(df.columns, SAMPLE_COLS)
    count_col = find_col(df.columns, COUNT_COLS)

    if ref_col is None or count_col is None:
        return pd.DataFrame()

    out = pd.DataFrame()
    out["cohort"] = cohort
    out["source_file"] = str(path)
    out["sample_id"] = df[sample_col].astype(str) if sample_col else infer_sample_from_path(path)
    out["reference_id"] = df[ref_col].astype(str).str.strip()

    counts = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    out["reads"] = counts.astype(int)

    for c in CATEGORY_COLS:
        if c in df.columns:
            out[c] = df[c].astype(str)
        else:
            out[c] = ""

    out["is_other_viral"] = out.apply(lambda r: classify_other_viral(r, refmap), axis=1)

    def get_refmap_cat(rid):
        info = refmap.get(rid) or refmap.get(str(rid).split(".")[0]) or {}
        return info.get("refmap_category", "")

    def get_refmap_name(rid):
        info = refmap.get(rid) or refmap.get(str(rid).split(".")[0]) or {}
        return info.get("refmap_name", "")

    out["refmap_category"] = out["reference_id"].map(get_refmap_cat)
    out["refmap_name"] = out["reference_id"].map(get_refmap_name)

    out = out[(out["reads"] > 0) & (out["is_other_viral"])]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-refmap", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--cohort",
        action="append",
        required=True,
        help="Format: cohort_name=/path/to/panel_run_dir"
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    refmap = load_refmap(args.panel_refmap)
    print(f"[INFO] Loaded refmap entries: {len(refmap)}")

    all_rows = []
    scanned = []

    for item in args.cohort:
        cohort, run_dir = item.split("=", 1)
        run_dir = Path(run_dir)

        if not run_dir.exists():
            print(f"[WARN] Missing run dir: {run_dir}")
            continue

        files = [p for p in run_dir.rglob("*") if p.is_file() and is_count_table(p)]

        # Prefer detailed filtered count tables if present.
        preferred = [
            p for p in files
            if re.search(r"(filtered.*count|count.*filtered|primary.*count|per.*reference|reference.*count)", p.name, re.I)
        ]
        use_files = preferred if preferred else files

        print(f"\n[INFO] {cohort}")
        print(f"  candidate count tables: {len(files)}")
        print(f"  using tables: {len(use_files)}")

        for p in use_files:
            scanned.append({"cohort": cohort, "file": str(p)})
            try:
                sub = normalize_table(p, cohort, refmap)
                if not sub.empty:
                    all_rows.append(sub)
                    print(f"  [USED] {p}  other_viral_rows={len(sub)}  reads={sub['reads'].sum()}")
            except Exception as e:
                print(f"  [SKIP] {p}  error={e}")

    pd.DataFrame(scanned).to_csv(outdir / "files_scanned.tsv", sep="\t", index=False)

    if not all_rows:
        print("\n[ERROR] No OTHER_VIRAL rows found.")
        print("Check files_scanned.tsv and make sure the script is reading per-reference count tables, not only species-level summaries.")
        return

    rows = pd.concat(all_rows, ignore_index=True)

    # Prevent exact duplicate rows from double-counting if the same table was discovered twice.
    rows = rows.drop_duplicates(
        subset=["cohort", "source_file", "sample_id", "reference_id", "reads"]
    )

    by_ref = (
        rows.groupby(["cohort", "reference_id", "refmap_category", "refmap_name"], dropna=False)
        .agg(
            samples_nonzero=("sample_id", "nunique"),
            total_reads=("reads", "sum"),
            max_reads_in_one_sample=("reads", "max"),
            files=("source_file", "nunique"),
        )
        .reset_index()
        .sort_values(["cohort", "total_reads"], ascending=[True, False])
    )

    by_sample_ref = (
        rows.groupby(["cohort", "sample_id", "reference_id", "refmap_category", "refmap_name"], dropna=False)
        .agg(reads=("reads", "sum"))
        .reset_index()
        .sort_values(["cohort", "sample_id", "reads"], ascending=[True, True, False])
    )

    by_sample = (
        rows.groupby(["cohort", "sample_id"], dropna=False)
        .agg(
            other_viral_references=("reference_id", "nunique"),
            other_viral_reads=("reads", "sum"),
        )
        .reset_index()
        .sort_values(["cohort", "other_viral_reads"], ascending=[True, False])
    )

    by_ref.to_csv(outdir / "other_viral_by_reference.tsv", sep="\t", index=False)
    by_sample_ref.to_csv(outdir / "other_viral_by_sample_reference.tsv", sep="\t", index=False)
    by_sample.to_csv(outdir / "other_viral_by_sample.tsv", sep="\t", index=False)

    print("\n[DONE] Wrote:")
    print(f"  {outdir / 'other_viral_by_reference.tsv'}")
    print(f"  {outdir / 'other_viral_by_sample_reference.tsv'}")
    print(f"  {outdir / 'other_viral_by_sample.tsv'}")
    print(f"  {outdir / 'files_scanned.tsv'}")

    print("\nTop OTHER_VIRAL references by cohort:")
    for cohort, g in by_ref.groupby("cohort"):
        print(f"\n### {cohort}")
        print(
            g[["reference_id", "samples_nonzero", "total_reads", "max_reads_in_one_sample", "refmap_name"]]
            .head(25)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
