from __future__ import annotations

import argparse
import csv
from pathlib import Path


VALID_STATUSES = {"positive", "negative", "unknown"}


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def ratio(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.6f}" if denominator else "NA"


def validate_labels(rows: list[dict[str, str]]) -> dict[str, str]:
    labels: dict[str, str] = {}
    for line_number, row in enumerate(rows, start=2):
        sample = row.get("sample_id", "").strip()
        status = row.get("status", "").strip().lower()
        if not sample:
            raise ValueError(f"Missing sample_id at label row {line_number}")
        if status not in VALID_STATUSES:
            raise ValueError(f"Invalid status for {sample}: {status!r}; expected positive, negative, or unknown")
        if sample in labels:
            raise ValueError(f"Duplicate label for sample_id: {sample}")
        labels[sample] = status
    return labels


def evaluate(
    count_rows: list[dict[str, str]], labels: dict[str, str], category: str, threshold: int
) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not count_rows or category not in count_rows[0]:
        raise ValueError(f"Count table does not contain category column: {category}")
    counts = {row["sample"]: int(float(row.get(category, "") or 0)) for row in count_rows}
    missing = sorted(set(labels) - set(counts))
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"{len(missing)} labelled samples are absent from counts: {preview}")

    rows: list[dict[str, object]] = []
    tp = fp = tn = fn = 0
    for sample, status in labels.items():
        count = counts[sample]
        detected = count >= threshold
        outcome = "not_evaluated"
        if status == "positive":
            outcome = "TP" if detected else "FN"
            tp += int(detected)
            fn += int(not detected)
        elif status == "negative":
            outcome = "FP" if detected else "TN"
            fp += int(detected)
            tn += int(not detected)
        rows.append(
            {
                "sample_id": sample,
                "known_status": status,
                "category": category,
                "count": count,
                "threshold": threshold,
                "detected": int(detected),
                "outcome": outcome,
            }
        )

    metrics: dict[str, object] = {
        "category": category,
        "threshold": threshold,
        "labelled_samples": len(labels),
        "positive_labels": sum(status == "positive" for status in labels.values()),
        "negative_labels": sum(status == "negative" for status in labels.values()),
        "unknown_labels": sum(status == "unknown" for status in labels.values()),
        "TP": tp,
        "FN": fn,
        "TN": tn,
        "FP": fp,
        "sensitivity": ratio(tp, tp + fn),
        "specificity": ratio(tn, tn + fp),
        "positive_predictive_value": ratio(tp, tp + fp),
        "negative_predictive_value": ratio(tn, tn + fn),
    }
    return rows, metrics


def write_report(path: Path, metrics: dict[str, object]) -> None:
    lines = [
        "# Known-Status Detection Performance",
        "",
        f"- Category: {metrics['category']}",
        f"- Operational threshold: >= {metrics['threshold']} coordinate-deduplicated alignments",
        f"- Labelled samples: {metrics['labelled_samples']}",
        f"- Positive / negative / unknown labels: {metrics['positive_labels']} / {metrics['negative_labels']} / {metrics['unknown_labels']}",
        f"- TP / FN / TN / FP: {metrics['TP']} / {metrics['FN']} / {metrics['TN']} / {metrics['FP']}",
        f"- Sensitivity: {metrics['sensitivity']}",
        f"- Specificity: {metrics['specificity']}",
        "",
        "The threshold is analytical and must not be described as clinically validated without an independent validation design.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate viral detection against authoritative sample-level labels.")
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True, help="CSV/TSV with sample_id,status")
    parser.add_argument("--category", required=True)
    parser.add_argument("--threshold", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.threshold < 1:
        parser.error("--threshold must be >= 1")

    labels = validate_labels(read_table(args.labels))
    rows, metrics = evaluate(read_table(args.counts), labels, args.category, args.threshold)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "sample_classification.tsv", rows, list(rows[0].keys()) if rows else [])
    write_tsv(args.output_dir / "performance_metrics.tsv", [metrics], list(metrics.keys()))
    write_report(args.output_dir / "performance_report.md", metrics)
    print(f"Wrote known-status performance outputs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
