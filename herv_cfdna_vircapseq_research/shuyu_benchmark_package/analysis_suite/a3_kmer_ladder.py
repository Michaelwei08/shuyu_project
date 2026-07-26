#!/usr/bin/env python3
"""a3_kmer_ladder.py -- parse and summarise the ALREADY-COMPUTED k-mer masking ladder.

WHAT THIS COMPUTES
------------------
Nothing is realigned or re-masked here. This module only reads tables that a
previous step already wrote, normalises them, and reports the trade-off across
the k-mer size k:

  1. How much HIV1 / HTLV1 sequence is masked at each k (the cost: real
     retroviral sequence removed from the reference, i.e. lost sensitivity).
  2. What happens to HERV cross-mapping at each k (the benefit: residual
     HIV1/HTLV1-vs-HERV exact-kmer similarity that can still mis-assign reads).

The input schema is treated as UNKNOWN. Every table is sniffed: delimiter,
whether the first row is a header, which column carries k, which columns are
identifiers, and which are numeric metrics. The detected columns are printed to
stdout and recorded in the report. Rows are melted into a normalised long form
(k, scope, metric, value). Nothing is guessed: if a column this module would
have liked is absent, the report lists what IS present instead (see the
"COLUMN AUDIT" section) rather than inventing a value.

The column vocabulary this module recognises comes from the lineage of
mask_shared_retro_regions.py (masked_pct / retained_pct / query_similar_pct /
post_mask_query_similar_pct / shared_query_kmers ...), but no such column is
required.

INPUTS (all optional; a missing input is WARNed about and skipped)
-----------------------------------------------------------------
  --ladder-summary    the ladder sweep table
                      (default .../reply_to_shuyu_primary_only/kmer_ladder_summary.tsv)
  --ladder-dir        per-k directory beside it (default .../kmer_ladder)
  --mask-metrics-dir  default .../retro_reference_hg38_refseq_mask_metrics_k40
  --masked-build-dir  optional masked reference build dir (its ref/ mask tables)

WHAT THIS WRITES (into --outdir, tab-separated, pure ASCII)
-----------------------------------------------------------
  kmer_ladder_long.tsv            normalised long form:
                                  k, scope, metric, value, value_num, source, row
  kmer_ladder_summary_report.txt  human-readable report: detected schema per
                                  file, k values found, per-k pivot of the
                                  mask-cost and HERV-cross series, direction and
                                  delta of each series, column audit, caveats.
  kmer_ladder_sample_key.tsv      ONLY if a sample-level field was found.
                                  Contains real identifiers -- never commit it.

IDENTIFIERS
-----------
No real sample identifier is ever written to kmer_ladder_long.tsv or to the
report. Any sample-like field is anonymised to S01..Snn (sorted by the real
name) and the real->anon mapping goes to kmer_ladder_sample_key.tsv alone. A
file basename is also treated as a sample name when it yields a cohort group
label and does not look generic (see SAMPLE_STEM_GENERIC). Every string written
out is additionally scrubbed against the collected real names.

Group labels are derived from the real sample name: "_HIV" -> HIV, "_HL" -> HL,
"targeted_htlv" or "TCL" -> TCL, otherwise NA.

OTHER NOTES
-----------
Standard library only. This module makes NO figures, so matplotlib is not
imported. It reads no BAMs, so samtools is not invoked. No network access.
Exits 0 even when every input is missing.

EXAMPLE
-------
  python3 a3_kmer_ladder.py --outdir ./a3_out

  python3 a3_kmer_ladder.py \
      --ladder-summary /path/to/runs/reply_to_shuyu_primary_only/kmer_ladder_summary.tsv \
      --ladder-dir     /path/to/runs/reply_to_shuyu_primary_only/kmer_ladder \
      --mask-metrics-dir /path/to/runs/retro_reference_hg38_refseq_mask_metrics_k40 \
      --include-md --outdir ./a3_out
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import time

# --------------------------------------------------------------------------
# defaults
# --------------------------------------------------------------------------

RUNS_ROOT = "/path/to/runs"
DEF_LADDER_SUMMARY = os.path.join(RUNS_ROOT, "reply_to_shuyu_primary_only",
                                  "kmer_ladder_summary.tsv")
DEF_LADDER_DIR = os.path.join(RUNS_ROOT, "reply_to_shuyu_primary_only", "kmer_ladder")
DEF_MASK_METRICS_DIR = os.path.join(RUNS_ROOT, "retro_reference_hg38_refseq_mask_metrics_k40")
DEF_MASKED_BUILD_DIR = os.path.join(RUNS_ROOT,
                                    "retro_reference_hg38_refseq_masked_hiv1_htlv1_vs_herv_k40")

TABLE_EXT = (".tsv", ".csv", ".txt", ".tab")
DELIMITERS = ["\t", ",", ";", "|"]

# Column names this module recognises, from the lineage of
# mask_shared_retro_regions.py.  None of them is required.
EXPECTED_COLUMNS = [
    "kmer_size", "reference_id", "category", "length",
    "masked_bases", "masked_pct", "retained_bases", "retained_pct",
    "query_reference_id", "query_category", "against_reference_id", "against_category",
    "shared_query_kmers", "query_length", "query_similar_bases", "query_similar_pct",
    "against_length", "against_similar_bases", "against_similar_pct",
    "post_mask_shared_query_kmers", "post_mask_query_similar_bases",
    "post_mask_query_similar_pct",
]

K_COLUMN_NAMES = {
    "k", "kmer", "kmersize", "ksize", "kmerlen", "kmerlength", "klen",
    "kvalue", "kmervalue", "kmerk", "kk", "sizek",
}

# Non-numeric columns that describe WHICH thing a row is about.
ID_COLUMN_NAMES = {
    "reference", "referenceid", "refid", "ref", "refname", "referencename",
    "queryreferenceid", "againstreferenceid", "queryref", "againstref",
    "category", "querycategory", "againstcategory", "cat", "class",
    "sample", "sampleid", "samplename", "samples",
    "group", "cohort", "run", "rundir", "rundirectory", "runname",
    "dataset", "panel", "build", "reference2", "label", "name", "type",
    "status", "mode", "filter", "stage", "step", "comparison", "pair",
}

# Identifier-ish columns kept out of the scope string: long, or path-like and so
# a possible identifier carrier.
SCOPE_EXCLUDE = {
    "description", "source", "path", "filepath", "file", "filename", "dir",
    "directory", "bam", "bampath", "fasta", "fastapath", "cmd", "command",
    "note", "notes", "comment", "comments", "againstreferenceids",
    "referenceids", "url",
}

SAMPLE_COLUMN_NAMES = {"sample", "sampleid", "samplename", "samples", "subjectid", "subject"}

# A file stem that yields a group label but matches one of these is a summary
# artefact, not a sample. group_label() matches case-insensitively, so this
# list also has to exclude the reference-space vocabulary of this module
# (kmer / herv / shared / region ...) or a table name such as
# "post_mask_hiv1_vs_herv" would be mistaken for a patient identifier.
SAMPLE_STEM_GENERIC = (
    "ladder", "summary", "report", "metric", "mask", "pair", "similar",
    "count", "manifest", "merged", "combined", "index", "readme", "all",
    "aggregate", "overall", "total", "stats", "final", "log", "config",
    "reference", "refmap", "categor",
    "kmer", "mer", "herv", "line1", "shared", "region", "retain", "audit",
    "profile", "table", "matrix", "bed", "fasta", "build",
)

MASK_COST = "MASK_COST"
HERV_CROSS = "HERV_CROSS"
OTHER = "OTHER"

BUCKET_TITLE = {
    MASK_COST: "MASK COST -- HIV1 / HTLV1 sequence removed or flagged as HERV-shared",
    HERV_CROSS: "HERV CROSS-MAPPING -- residual HERV-side similarity or HERV signal",
    OTHER: "OTHER SERIES THAT VARY WITH k",
}


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def norm(name):
    """Normalise a column name for matching: lowercase, alphanumerics only."""
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


def ascii_only(text):
    """Force pure ASCII; every output file must be ASCII."""
    if text is None:
        return ""
    return str(text).encode("ascii", "replace").decode("ascii")


def clean_cell(text):
    """One-line, tab-free, ASCII cell value."""
    return ascii_only(re.sub(r"[\r\n\t]+", " ", str(text if text is not None else ""))).strip()


def to_num(text):
    """Parse a cell as a number, tolerating thousands separators and '%'. None if not numeric."""
    if text is None:
        return None
    t = str(text).strip().replace(",", "")
    if t.endswith("%"):
        t = t[:-1].strip()
    if t in ("", "NA", "na", "N/A", "n/a", "-", ".", "None", "null", "NaN", "nan"):
        return None
    try:
        return float(t)
    except ValueError:
        return None


def fmt_num(value):
    """Format a float for a TSV cell without losing integrality."""
    if value is None:
        return ""
    if abs(value) < 1e15 and float(value).is_integer():
        return str(int(value))
    return repr(round(value, 6))


def fmt_report_num(value):
    """Format a float for the human report."""
    if value is None:
        return "NA"
    if abs(value) < 1e15 and float(value).is_integer():
        return "{0:,}".format(int(value))
    return "{0:.4f}".format(value)


def group_label(sample_name):
    """Cohort group from a real sample name. Order is significant.

    Suite-wide rule, matched case-insensitively so every module agrees:
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


def warn(what, path):
    print('WARN: {0} missing at {1}, skipping'.format(ascii_only(what), ascii_only(path)))


def looks_like_sample_stem(stem):
    """True when a bare file stem should be treated as a real sample identifier."""
    if group_label(stem) == "NA":
        return False
    low = stem.lower()
    for token in SAMPLE_STEM_GENERIC:
        if token in low:
            return False
    return True


def derive_k_from_text(text):
    """Pull a k value out of a filename, directory name, or cell (k40, k_40, k=40, kmer40)."""
    if not text:
        return None
    for pattern in (r"(?:^|[^0-9A-Za-z])k(?:mer)?[ _\-=:]?(\d{1,3})(?:$|[^0-9])",
                    r"(?:^|[^0-9A-Za-z])(\d{1,3})[ _\-]?mer(?:$|[^0-9A-Za-z])"):
        m = re.search(pattern, str(text))
        if m:
            try:
                value = int(m.group(1))
            except ValueError:
                continue
            if 5 <= value <= 500:
                return str(value)
    return None


def derive_k_from_path(path):
    """Look for k in the basename first, then walk up the directory components."""
    parts = os.path.normpath(path).replace("\\", "/").split("/")
    for part in reversed([p for p in parts if p]):
        found = derive_k_from_text(part)
        if found:
            return found
    return None


def k_sort_key(k):
    """Numeric k first in ascending order, non-numeric k last."""
    try:
        return (0, float(k), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(k))


# --------------------------------------------------------------------------
# delimiter / header sniffing
# --------------------------------------------------------------------------

def read_lines(path):
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return handle.read().splitlines()


def split_comments(lines):
    """Separate leading '#' comment lines and blanks from the data body."""
    comments, body = [], []
    for line in lines:
        if not body and (not line.strip() or line.lstrip().startswith("#")):
            if line.strip():
                comments.append(line.strip())
            continue
        body.append(line)
    while body and not body[-1].strip():
        body.pop()
    return comments, body


def sniff_delimiter(body):
    """Pick the delimiter that splits the sample lines most consistently."""
    sample_lines = [ln for ln in body[:40] if ln.strip()]
    if not sample_lines:
        return "\t"
    try:
        dialect = csv.Sniffer().sniff("\n".join(sample_lines[:20]), delimiters="\t,;|")
        if dialect.delimiter in DELIMITERS:
            counts = [ln.count(dialect.delimiter) for ln in sample_lines]
            if counts and min(counts) >= 1:
                return dialect.delimiter
    except (csv.Error, TypeError):
        pass
    best, best_score = "\t", (-1, -1)
    for delim in DELIMITERS:
        counts = [ln.count(delim) for ln in sample_lines]
        if not counts or min(counts) < 1:
            continue
        consistent = 1 if len(set(counts)) == 1 else 0
        score = (consistent, min(counts))
        if score > best_score:
            best, best_score = delim, score
    return best


def is_markdown_table(body):
    rows = [ln.strip() for ln in body if ln.strip()]
    if len(rows) < 3:
        return False
    if not rows[0].startswith("|"):
        return False
    return bool(re.match(r"^\|[\s:\-\|]+\|$", rows[1]))


def parse_markdown_table(body):
    """Rows of a leading GitHub-style markdown table; the dashes row is dropped."""
    rows = []
    for line in body:
        line = line.strip()
        if not line.startswith("|"):
            if rows:
                break
            continue
        if re.match(r"^\|[\s:\-\|]+\|$", line):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        rows.append(cells)
    return rows


def looks_like_header(cells):
    """A row is a header when its cells are mostly non-numeric text."""
    filled = [c for c in cells if str(c).strip()]
    if not filled:
        return False
    for cell in filled:
        if norm(cell) in K_COLUMN_NAMES or norm(cell) in ID_COLUMN_NAMES:
            return True
        if norm(cell) in (norm(e) for e in EXPECTED_COLUMNS):
            return True
    non_numeric = sum(1 for c in filled if to_num(c) is None)
    return non_numeric >= max(1, int(0.6 * len(filled)))


def sniff_table(path, include_md=False):
    """Sniff and read one tabular file.

    Returns a schema dict, or None when the file cannot be read as a table.
    """
    try:
        lines = read_lines(path)
    except OSError as exc:
        print("WARN: cannot read {0} ({1}), skipping".format(ascii_only(path), ascii_only(exc)))
        return None

    comments, body = split_comments(lines)
    if not body:
        return None

    md = is_markdown_table(body)
    if md and not include_md:
        return None
    if md:
        rows = parse_markdown_table(body)
        delimiter = "|"
    else:
        delimiter = sniff_delimiter(body)
        rows = [r for r in csv.reader(body, delimiter=delimiter) if any(str(c).strip() for c in r)]
    if not rows:
        return None

    if looks_like_header(rows[0]):
        raw_columns = [clean_cell(c) for c in rows[0]]
        data_rows = rows[1:]
        header_detected = True
    else:
        raw_columns = ["col{0}".format(i + 1) for i in range(len(rows[0]))]
        data_rows = rows
        header_detected = False

    columns = []
    for index, name in enumerate(raw_columns):
        name = name or "col{0}".format(index + 1)
        candidate, bump = name, 2
        while candidate in columns:
            candidate = "{0}_{1}".format(name, bump)
            bump += 1
        columns.append(candidate)

    records = []
    for row in data_rows:
        record = {}
        for index, column in enumerate(columns):
            record[column] = clean_cell(row[index]) if index < len(row) else ""
        records.append(record)

    return {
        "path": path,
        "delimiter": delimiter,
        "markdown": md,
        "header_detected": header_detected,
        "comments": comments,
        "columns": columns,
        "records": records,
        "preview": body[:6],
    }


# --------------------------------------------------------------------------
# column classification
# --------------------------------------------------------------------------

def classify_columns(schema):
    """Split columns into: the k column, identifier columns, and metric columns."""
    columns, records = schema["columns"], schema["records"]

    k_col = None
    for column in columns:
        n = norm(column)
        if n in K_COLUMN_NAMES or ("kmer" in n and ("size" in n or "len" in n)):
            values = [r.get(column, "") for r in records]
            if any(to_num(v) is not None or derive_k_from_text(v) for v in values):
                k_col = column
                break

    metric_col = value_col = None
    for column in columns:
        if norm(column) in ("metric", "measure", "statistic", "stat", "variable", "key"):
            metric_col = column
            break
    if metric_col:
        for column in columns:
            if norm(column) in ("value", "val", "number", "result", "amount"):
                value_col = column
                break
        if not value_col:
            metric_col = None

    id_cols, metric_cols = [], []
    for column in columns:
        if column in (k_col, metric_col, value_col):
            continue
        n = norm(column)
        values = [r.get(column, "") for r in records if str(r.get(column, "")).strip()]
        numeric = sum(1 for v in values if to_num(v) is not None)
        mostly_numeric = bool(values) and numeric >= 0.8 * len(values)
        if n in ID_COLUMN_NAMES or (not mostly_numeric and n not in SCOPE_EXCLUDE):
            id_cols.append(column)
        elif n in SCOPE_EXCLUDE and not mostly_numeric:
            metric_cols.append(column)      # kept as a text metric, out of scope
        else:
            metric_cols.append(column)

    schema["k_col"] = k_col
    schema["metric_col"] = metric_col
    schema["value_col"] = value_col
    schema["id_cols"] = id_cols
    schema["metric_cols"] = metric_cols
    schema["sample_cols"] = [c for c in columns if norm(c) in SAMPLE_COLUMN_NAMES
                             or ("sample" in norm(c) and norm(c) not in ("samplecount", "nsamples"))]
    return schema


def scope_pairs_for(schema, record):
    """Ordered key=value pairs describing which entity a row is about."""
    pairs = []
    for column in schema["id_cols"]:
        if norm(column) in SCOPE_EXCLUDE:
            continue
        value = str(record.get(column, "")).strip()
        if value:
            pairs.append((column, value))
    return pairs


def scope_string(pairs):
    if not pairs:
        return "ALL"
    return ";".join("{0}={1}".format(k, v) for k, v in pairs)


# --------------------------------------------------------------------------
# melting to long form
# --------------------------------------------------------------------------

def melt(schema, source_label):
    """Turn one sniffed table into long-form records."""
    out = []
    k_col = schema["k_col"]
    path_k = derive_k_from_path(schema["path"])

    for row_index, record in enumerate(schema["records"], start=1):
        k = None
        if k_col:
            raw = str(record.get(k_col, "")).strip()
            num = to_num(raw)
            if num is not None:
                k = fmt_num(num)
            else:
                k = derive_k_from_text(raw)
        if k is None:
            for column in schema["id_cols"]:
                k = derive_k_from_text(record.get(column, ""))
                if k:
                    break
        if k is None:
            k = path_k
        if k is None:
            k = "NA"

        pairs = scope_pairs_for(schema, record)

        if schema["metric_col"] and schema["value_col"]:
            metric = str(record.get(schema["metric_col"], "")).strip() or "value"
            value = str(record.get(schema["value_col"], "")).strip()
            out.append({"k": k, "scope_pairs": pairs, "metric": metric, "value": value,
                        "source": source_label, "row": row_index})
            continue

        for column in schema["metric_cols"]:
            value = str(record.get(column, "")).strip()
            if value == "":
                continue
            out.append({"k": k, "scope_pairs": pairs, "metric": column, "value": value,
                        "source": source_label, "row": row_index})
    return out


def bed_records(path, source_label):
    """Masked bp and interval count per reference from a 4-column mask BED."""
    per_ref = {}
    try:
        with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start, end = int(parts[1]), int(parts[2])
                except ValueError:
                    continue
                bucket = per_ref.setdefault(parts[0], [0, 0])
                bucket[0] += 1
                bucket[1] += max(0, end - start)
    except OSError as exc:
        print("WARN: cannot read {0} ({1}), skipping".format(ascii_only(path), ascii_only(exc)))
        return []

    k = derive_k_from_path(path) or "NA"
    out = []
    for row_index, ref in enumerate(sorted(per_ref), start=1):
        intervals, bases = per_ref[ref]
        pairs = [("reference_id", ref)]
        out.append({"k": k, "scope_pairs": pairs, "metric": "mask_bed_intervals",
                    "value": str(intervals), "source": source_label, "row": row_index})
        out.append({"k": k, "scope_pairs": pairs, "metric": "mask_bed_masked_bp",
                    "value": str(bases), "source": source_label, "row": row_index})
    return out


# --------------------------------------------------------------------------
# input discovery
# --------------------------------------------------------------------------

def discover_dir(path, what, max_files, max_depth):
    """Tabular files plus BED and markdown files under a directory, breadth-limited."""
    if not path:
        return [], [], []
    if not os.path.isdir(path):
        warn(what, path)
        return [], [], []
    tables, beds, others = [], [], []
    root_depth = os.path.normpath(path).count(os.sep)
    for current, dirnames, filenames in os.walk(path):
        if os.path.normpath(current).count(os.sep) - root_depth >= max_depth:
            dirnames[:] = []
        dirnames.sort()
        for name in sorted(filenames):
            full = os.path.join(current, name)
            lower = name.lower()
            if lower.endswith(".bed"):
                beds.append(full)
            elif lower.endswith(TABLE_EXT) or lower.endswith(".md"):
                tables.append(full)
            else:
                others.append(full)
    if len(tables) > max_files:
        print("WARN: {0} at {1} has {2} tabular files, parsing only the first {3}".format(
            ascii_only(what), ascii_only(path), len(tables), max_files))
        tables = tables[:max_files]
    return tables, beds, others


# --------------------------------------------------------------------------
# anonymisation
# --------------------------------------------------------------------------

def build_sample_map(long_records, schemas, extra_names):
    """Collect real sample names, then map them to S01..Snn sorted by real name."""
    real = set()
    for schema in schemas:
        for column in schema.get("sample_cols", []):
            for record in schema["records"]:
                value = str(record.get(column, "")).strip()
                if value and value != "*":
                    real.add(value)
    real.update(n for n in extra_names if n)

    ordered = sorted(real)
    width = max(2, len(str(len(ordered))))
    return {name: "S{0}".format(str(i + 1).zfill(width)) for i, name in enumerate(ordered)}


def scrub(text, sample_map):
    """Replace every real sample name with its anon ID. Longest names first."""
    out = str(text)
    for real in sorted(sample_map, key=len, reverse=True):
        if real and real in out:
            out = out.replace(real, sample_map[real])
    return out


# --------------------------------------------------------------------------
# trade-off analysis
# --------------------------------------------------------------------------

def classify_series(scope_pairs, metric):
    """Bucket a (scope, metric) series as mask cost, HERV cross-mapping, or other."""
    text = (" ".join(v for _, v in scope_pairs) + " " + metric).upper()
    metric_low = metric.lower()
    masky = any(w in metric_low for w in ("mask", "similar", "shared", "retain"))

    subject = ""
    for key, value in scope_pairs:
        low = key.lower()
        if "query" in low and ("categ" in low or "ref" in low):
            subject = value.upper()
            break
    if not subject:
        for key, value in scope_pairs:
            if norm(key) in ("category", "cat", "class", "referenceid", "reference", "ref"):
                subject = value.upper()
                break

    if subject.startswith("HIV") or subject.startswith("HTLV"):
        return MASK_COST
    if "HERV" in subject:
        return HERV_CROSS

    # A column named after a viral category is that category's read/count readout,
    # i.e. the sensitivity side of the trade-off.
    metric_norm = norm(metric)
    if metric_norm in ("hiv1", "hiv2", "htlv1", "htlv2", "hiv", "htlv"):
        return MASK_COST
    if metric_norm == "herv":
        return HERV_CROSS

    if masky and ("HIV" in text or "HTLV" in text):
        return MASK_COST
    if "HERV" in text:
        return HERV_CROSS
    if masky:
        return MASK_COST
    return OTHER


def build_series(long_records):
    """Group long records into series keyed by (bucket, source, scope, metric)."""
    series = {}
    for record in long_records:
        num = to_num(record["value"])
        if num is None:
            continue
        scope = scope_string(record["scope_pairs"])
        key = (classify_series(record["scope_pairs"], record["metric"]),
               record["source"], scope, record["metric"])
        series.setdefault(key, {}).setdefault(record["k"], []).append(num)
    return series


def series_direction(k_values, per_k):
    """monotone-up / monotone-down / mixed / flat across ascending k."""
    values = [sum(per_k[k]) / len(per_k[k]) for k in k_values]
    ups = downs = 0
    for before, after in zip(values, values[1:]):
        if after > before:
            ups += 1
        elif after < before:
            downs += 1
    if ups and downs:
        return "mixed"
    if ups:
        return "monotone-up"
    if downs:
        return "monotone-down"
    return "flat"


def render_table(headers, rows):
    """Fixed-width ASCII table."""
    if not rows:
        return ["(no rows)"]
    widths = [len(str(h)) for h in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(str(cell)))
    out = ["  ".join(str(h).ljust(widths[i]) for i, h in enumerate(headers)).rstrip()]
    out.append("  ".join("-" * widths[i] for i in range(len(headers))))
    for row in rows:
        out.append("  ".join(str(c).ljust(widths[i]) for i, c in enumerate(row)).rstrip())
    return out


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

def write_report(report_path, args, schemas, long_records, series, k_values, sample_map,
                 missing, bed_files, md_files, other_files):
    lines = []
    add = lines.append

    add("=" * 96)
    add("K-MER MASKING LADDER -- PARSED SUMMARY")
    add("=" * 96)
    add("generated            : {0}".format(time.strftime("%Y-%m-%d")))
    add("module               : a3_kmer_ladder.py")
    add("ladder summary       : {0}".format(ascii_only(args.ladder_summary)))
    add("ladder directory     : {0}".format(ascii_only(args.ladder_dir)))
    add("mask metrics dir     : {0}".format(ascii_only(args.mask_metrics_dir)))
    add("masked build dir     : {0}".format(ascii_only(args.masked_build_dir or "(not given)")))
    add("markdown tables      : {0}".format("parsed (--include-md)" if args.include_md
                                            else "not parsed (pass --include-md)"))
    add("tables parsed        : {0}".format(len(schemas)))
    add("long-form rows       : {0}".format(len(long_records)))
    add("k values found       : {0}".format(", ".join(k_values) if k_values else "(none)"))
    add("samples anonymised   : {0}".format(len(sample_map)))
    add("")

    if missing:
        add("-" * 96)
        add("MISSING INPUTS (skipped, not fatal)")
        add("-" * 96)
        for what, where in missing:
            add("  {0}: {1}".format(ascii_only(what), ascii_only(where)))
        add("")

    add("-" * 96)
    add("DETECTED SCHEMA PER FILE (the input schema was sniffed, not assumed)")
    add("-" * 96)
    if not schemas:
        add("  No tabular input could be parsed.")
    for schema in schemas:
        delim_name = {"\t": "TAB", ",": "COMMA", ";": "SEMICOLON", "|": "PIPE"}.get(
            schema["delimiter"], repr(schema["delimiter"]))
        add("  file      : {0}".format(scrub(ascii_only(schema["path"]), sample_map)))
        add("  delimiter : {0}{1}".format(delim_name,
                                          "  (markdown table)" if schema["markdown"] else ""))
        add("  header    : {0}".format("row 1" if schema["header_detected"]
                                       else "NOT DETECTED - columns named col1..colN"))
        add("  rows      : {0}".format(len(schema["records"])))
        add("  columns   : {0}".format(", ".join(scrub(c, sample_map) for c in schema["columns"])))
        add("  k column  : {0}".format(scrub(schema["k_col"], sample_map)
                                       if schema["k_col"]
                                       else "NONE - k taken from the file path"))
        # a wide table can carry sample names as column headers, so scrub these
        add("  id cols   : {0}".format(
            ", ".join(scrub(c, sample_map) for c in schema["id_cols"]) or "(none)"))
        add("  metrics   : {0}".format(
            ", ".join(scrub(c, sample_map) for c in schema["metric_cols"]) or "(none)"))
        if schema["metric_col"]:
            add("  long form : already long ({0} / {1})".format(schema["metric_col"],
                                                                schema["value_col"]))
        if schema["comments"]:
            add("  comments  : {0}".format(scrub(ascii_only(" | ".join(schema["comments"][:3])),
                                                 sample_map)))
        if args.preview_lines > 0:
            for raw in schema["preview"][:args.preview_lines]:
                add("  raw       : {0}".format(scrub(clean_cell(raw), sample_map)[:180]))
        add("")

    add("-" * 96)
    add("COLUMN AUDIT (what was looked for vs what is actually present)")
    add("-" * 96)
    present = set()
    for schema in schemas:
        present.update(norm(c) for c in schema["columns"])
    found = [c for c in EXPECTED_COLUMNS if norm(c) in present]
    absent = [c for c in EXPECTED_COLUMNS if norm(c) not in present]
    add("  Recognised columns (lineage of mask_shared_retro_regions.py) that ARE present:")
    add("    {0}".format(", ".join(found) if found else "(none)"))
    add("  Recognised columns that are ABSENT - no value was guessed for these:")
    add("    {0}".format(", ".join(absent) if absent else "(none)"))
    extra = sorted(set()) if not schemas else sorted(
        {scrub(c, sample_map) for s in schemas for c in s["columns"]
         if norm(c) not in {norm(e) for e in EXPECTED_COLUMNS}})
    add("  Columns present that were NOT in the recognised list (reported, not interpreted):")
    add("    {0}".format(", ".join(extra) if extra else "(none)"))
    add("")

    add("-" * 96)
    add("TRADE-OFF ACROSS k")
    add("-" * 96)
    if not series or len(k_values) < 2:
        add("  Fewer than two distinct k values carry numeric data, so no ladder trade-off")
        add("  can be computed. Everything parsed is still in the long-form TSV.")
        add("")
    else:
        shown_k = k_values[:args.max_k_columns]
        # (metric, scope) pairs that occur under more than one source would render as
        # duplicate rows, so those rows carry the source too.
        source_count = {}
        for key in series:
            source_count.setdefault((key[3], key[2]), set()).add(key[1])
        ambiguous = {pair for pair, sources in source_count.items() if len(sources) > 1}
        if len(k_values) > len(shown_k):
            add("  NOTE: {0} k values found, showing the first {1}. Full data is in the TSV.".format(
                len(k_values), len(shown_k)))
            add("")
        for bucket in (MASK_COST, HERV_CROSS, OTHER):
            rows = []
            for key in sorted(series):
                if key[0] != bucket:
                    continue
                per_k = series[key]
                ks = [k for k in k_values if k in per_k]
                if len(ks) < 2:
                    continue
                # Metric first so it is never the part lost to truncation.
                # series keys hold RAW values, so scrub before the label is shown.
                label = "{0} | {1}".format(key[3], key[2])
                if (key[3], key[2]) in ambiguous:
                    label = "{0} @{1}".format(label, key[1])
                label = scrub(label, sample_map)
                if len(label) > 76:
                    label = label[:73] + "..."
                row = [label]
                for k in shown_k:
                    if k in per_k:
                        values = per_k[k]
                        mean = sum(values) / len(values)
                        cell = fmt_report_num(mean)
                        if len(values) > 1:
                            cell = "{0}(n{1})".format(cell, len(values))
                    else:
                        cell = "."
                    row.append(cell)
                first = sum(per_k[ks[0]]) / len(per_k[ks[0]])
                last = sum(per_k[ks[-1]]) / len(per_k[ks[-1]])
                row.append(fmt_report_num(last - first))
                row.append(series_direction(ks, per_k))
                rows.append(row)
            if not rows:
                continue
            add("  {0}".format(BUCKET_TITLE[bucket]))
            add("")
            headers = ["metric | scope"] + ["k={0}".format(k) for k in shown_k] + \
                      ["delta(first->last k)", "direction"]
            for line in render_table(headers, rows):
                add("    " + line)
            add("")

        cost = [k for k in series if k[0] == MASK_COST]
        herv = [k for k in series if k[0] == HERV_CROSS]
        add("  READING THIS TABLE")
        add("    Mask-cost series: {0}. HERV-cross series: {1}.".format(len(cost), len(herv)))
        add("    A mask-cost value that falls as k rises means less HIV1/HTLV1 sequence is")
        add("    removed at large k, so more of the real retroviral genome is kept.")
        add("    A HERV-cross value that rises as k rises means more residual HERV-shared")
        add("    sequence survives the mask, so more scope for HERV reads to land on HIV1/HTLV1.")
        add("")

    add("-" * 96)
    add("NOTE ON MECHANISM (not measured by this module)")
    add("-" * 96)
    add("  k sets how long an exact match must be before a position counts as HERV-shared.")
    add("  Larger k is stricter: fewer shared k-mers, less HIV1/HTLV1 sequence masked, more")
    add("  detection sensitivity retained, and more residual cross-mapping risk. Smaller k is")
    add("  the opposite trade. The ladder exists to choose that operating point; this module")
    add("  reports where each k sits, it does not recommend one.")
    add("")

    add("-" * 96)
    add("CAVEATS")
    add("-" * 96)
    add("  1. These are exact shared-kmer masks. That is not the VirCAPP production masking")
    add("     rule, so numbers here are not interchangeable with a vendor-masked panel.")
    add("  2. Reference-space metrics (masked bp, similarity percent) are not read counts. A")
    add("     mask percent does not translate directly into lost reads.")
    add("  3. Numbers from different source files are only comparable when the same reference")
    add("     build and the same category definitions were used. The source column is kept in")
    add("     the long-form TSV so this stays checkable.")
    add("  4. Where a cell held several values for one (k, scope, metric), the report shows the")
    add("     mean and marks it (nN). The TSV keeps every original row.")
    if bed_files:
        add("  5. Mask BED files were summarised as interval count and masked bp per reference:")
        for bed_path in bed_files:
            add("     {0}".format(scrub(ascii_only(bed_path), sample_map)))
    if md_files and not args.include_md:
        add("  6. Markdown tables were found but not parsed. Rerun with --include-md to fold")
        add("     them in:")
        for md_path in md_files[:6]:
            add("     {0}".format(scrub(ascii_only(md_path), sample_map)))
    if other_files:
        add("  7. Non-tabular files present and ignored: {0}".format(len(other_files)))
    add("")
    if sample_map:
        add("  Sample-level fields were present. They are anonymised to S01..Snn here; the real")
        add("  mapping is in {0}_sample_key.tsv, which must not be committed or emailed.".format(
            args.prefix))
    else:
        add("  No sample-level field was found in any parsed table, so no sample key was")
        add("  written. Nothing in this report or in the TSV is a sample identifier.")
    add("")

    text = "\n".join(ascii_only(line) for line in lines) + "\n"
    with open(report_path, "w", encoding="ascii", errors="replace", newline="") as handle:
        handle.write(text)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Parse and summarise the already-computed k-mer masking ladder.")
    parser.add_argument("--ladder-summary", default=DEF_LADDER_SUMMARY,
                        help="ladder sweep table (default: %(default)s)")
    parser.add_argument("--ladder-dir", default=DEF_LADDER_DIR,
                        help="per-k ladder directory (default: %(default)s)")
    parser.add_argument("--mask-metrics-dir", default=DEF_MASK_METRICS_DIR,
                        help="mask metrics directory (default: %(default)s)")
    parser.add_argument("--masked-build-dir", default=DEF_MASKED_BUILD_DIR,
                        help="masked reference build dir; its tables are folded in "
                             "(default: %(default)s). Pass '' to skip.")
    parser.add_argument("--extra-table", action="append", default=[],
                        help="additional table to fold in; repeatable")
    parser.add_argument("--outdir", default="./a3_kmer_ladder_out",
                        help="output directory (default: %(default)s)")
    parser.add_argument("--prefix", default="kmer_ladder",
                        help="output filename prefix (default: %(default)s)")
    parser.add_argument("--include-md", action="store_true",
                        help="also parse markdown tables such as mask_report.md")
    parser.add_argument("--preview-lines", type=int, default=2,
                        help="raw input lines echoed per file in the report (default: %(default)s)")
    parser.add_argument("--max-k-columns", type=int, default=12,
                        help="k columns shown in the report tables (default: %(default)s)")
    parser.add_argument("--max-files", type=int, default=400,
                        help="cap on tabular files read per directory (default: %(default)s)")
    parser.add_argument("--max-depth", type=int, default=3,
                        help="directory recursion depth (default: %(default)s)")
    return parser


def main():
    args = build_parser().parse_args()

    missing = []
    tables, beds, md_files, other_files = [], [], [], []

    if args.ladder_summary:
        if os.path.isfile(args.ladder_summary):
            tables.append(args.ladder_summary)
        else:
            warn("ladder summary table", args.ladder_summary)
            missing.append(("ladder summary table", args.ladder_summary))

    for path, what in ((args.ladder_dir, "ladder directory"),
                       (args.mask_metrics_dir, "mask metrics directory"),
                       (args.masked_build_dir, "masked reference build directory")):
        if not path:
            continue
        if not os.path.isdir(path):
            warn(what, path)
            missing.append((what, path))
            continue
        found_tables, found_beds, found_other = discover_dir(
            path, what, args.max_files, args.max_depth)
        tables.extend(found_tables)
        beds.extend(found_beds)
        other_files.extend(found_other)

    for path in args.extra_table:
        if os.path.isfile(path):
            tables.append(path)
        else:
            warn("extra table", path)
            missing.append(("extra table", path))

    seen = set()
    ordered_tables = []
    for path in tables:
        key = os.path.normpath(path)
        if key not in seen:
            seen.add(key)
            ordered_tables.append(path)

    # ---- parse ----
    schemas, long_records, stem_samples = [], [], set()
    for path in ordered_tables:
        if path.lower().endswith(".md"):
            md_files.append(path)
            if not args.include_md:
                continue
        schema = sniff_table(path, include_md=args.include_md)
        if not schema:
            continue
        classify_columns(schema)
        stem = os.path.splitext(os.path.basename(path))[0]
        if looks_like_sample_stem(stem):
            stem_samples.add(stem)
            source_label = "STEM:" + stem
        else:
            source_label = os.path.basename(path)
        schemas.append(schema)
        long_records.extend(melt(schema, source_label))

    for path in beds:
        long_records.extend(bed_records(path, os.path.basename(path)))

    print("=" * 78)
    print("a3_kmer_ladder.py  --  {0}".format(time.strftime("%Y-%m-%d")))
    print("=" * 78)
    if not schemas and not long_records:
        print("No parsable ladder or mask-metric table was found. Nothing to summarise.")

    sample_map = build_sample_map(long_records, schemas, stem_samples)

    for schema in schemas:
        delim_name = {"\t": "TAB", ",": "COMMA", ";": "SEMICOLON", "|": "PIPE"}.get(
            schema["delimiter"], repr(schema["delimiter"]))
        print("")
        print("FILE   {0}".format(scrub(ascii_only(schema["path"]), sample_map)))
        print("  delimiter={0}  header={1}  rows={2}  cols={3}".format(
            delim_name, "yes" if schema["header_detected"] else "NO (synthesised names)",
            len(schema["records"]), len(schema["columns"])))
        print("  columns : {0}".format(", ".join(scrub(c, sample_map)
                                                 for c in schema["columns"])))
        print("  k column: {0}".format(scrub(schema["k_col"], sample_map)
                                       if schema["k_col"]
                                       else "none (k derived from path)"))

    k_values = sorted({r["k"] for r in long_records}, key=k_sort_key)
    series = build_series(long_records)

    # ---- write ----
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except OSError as exc:
        print("WARN: cannot create output directory at {0} ({1}), skipping all output".format(
            ascii_only(args.outdir), ascii_only(exc)))
        return 0

    long_path = os.path.join(args.outdir, "{0}_long.tsv".format(args.prefix))
    report_path = os.path.join(args.outdir, "{0}_summary_report.txt".format(args.prefix))
    key_path = os.path.join(args.outdir, "{0}_sample_key.tsv".format(args.prefix))

    # Inputs are immutable: never let an output path collide with a file we read.
    inputs_seen = {os.path.abspath(p) for p in ordered_tables} | {os.path.abspath(p) for p in beds}
    for out_path in (long_path, report_path, key_path):
        if os.path.abspath(out_path) in inputs_seen:
            print("WARN: output {0} would overwrite an input file, refusing to write "
                  "anything; choose a different --outdir or --prefix".format(
                      ascii_only(out_path)))
            return 0

    rows = []
    for record in long_records:
        scope = scrub(scope_string(record["scope_pairs"]), sample_map)
        rows.append([
            record["k"],
            ascii_only(scope),
            ascii_only(scrub(record["metric"], sample_map)),
            ascii_only(scrub(record["value"], sample_map)),
            fmt_num(to_num(record["value"])),
            ascii_only(scrub(record["source"], sample_map)),
            str(record["row"]),
        ])
    rows.sort(key=lambda r: (k_sort_key(r[0]), r[5], r[1], r[2], int(r[6])))

    with open(long_path, "w", encoding="ascii", errors="replace", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["k", "scope", "metric", "value", "value_num", "source", "row"])
        writer.writerows(rows)

    write_report(report_path, args, schemas, long_records, series, k_values, sample_map,
                 missing, beds, md_files, other_files)

    if sample_map:
        with open(key_path, "w", encoding="ascii", errors="replace", newline="") as handle:
            handle.write("# CONTAINS IDENTIFIERS - DO NOT COMMIT OR EMAIL\n")
            handle.write("# generated {0} by a3_kmer_ladder.py\n".format(
                time.strftime("%Y-%m-%d")))
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["real_sample", "anon_sample", "group"])
            for real in sorted(sample_map):
                writer.writerow([ascii_only(real), sample_map[real], group_label(real)])

    # ---- stdout summary ----
    print("")
    print("-" * 78)
    print("SUMMARY")
    print("-" * 78)
    print("tables parsed   : {0}".format(len(schemas)))
    print("mask BED files  : {0}".format(len(beds)))
    print("long-form rows  : {0}".format(len(rows)))
    print("k values        : {0}".format(", ".join(k_values) if k_values else "(none)"))
    print("numeric series  : {0}  (mask-cost {1}, HERV-cross {2}, other {3})".format(
        len(series),
        sum(1 for key in series if key[0] == MASK_COST),
        sum(1 for key in series if key[0] == HERV_CROSS),
        sum(1 for key in series if key[0] == OTHER)))
    print("missing inputs  : {0}".format(len(missing)))
    print("samples anon'd  : {0}".format(len(sample_map)))
    print("")
    print("wrote {0}".format(long_path))
    print("wrote {0}".format(report_path))
    if sample_map:
        print("wrote {0}   <- CONTAINS IDENTIFIERS, do not commit or email".format(key_path))
    else:
        print("no sample-level field found; no sample key written")
    print("")
    print("Read {0} for the detected schema, the per-k pivot, and the column audit.".format(
        os.path.basename(report_path)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
