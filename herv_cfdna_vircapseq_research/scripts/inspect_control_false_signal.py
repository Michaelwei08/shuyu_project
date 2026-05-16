from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


RETRO = {"hiv", "htlv", "herv", "line1"}


def load_control_samples(manifest_csv: Path) -> set[str]:
    controls: set[str] = set()
    with manifest_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("expected_group") == "healthy_negative":
                controls.add(row["sample_id"])
    return controls


def inspect(manifest_csv: Path, hits_csv: Path, output_csv: Path, score_delta: float) -> None:
    controls = load_control_samples(manifest_csv)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with hits_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["sample_id"] in controls:
                grouped[(row["sample_id"], row["read_id"])].append(row)

    fields = ["sample_id", "read_id", "top_category", "competing_categories", "top_score", "score_gap", "review_reason"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for (sample_id, read_id), rows in sorted(grouped.items()):
            ranked = sorted(rows, key=lambda row: float(row["alignment_score"]), reverse=True)
            top = ranked[0]
            if top["reference_category"].lower() not in {"hiv", "htlv"}:
                continue
            top_score = float(top["alignment_score"])
            competitors = [row for row in ranked[1:] if row["reference_category"].lower() in RETRO]
            if competitors:
                best_competitor = max(float(row["alignment_score"]) for row in competitors)
                gap = top_score - best_competitor
                reason = "near_tie_with_endogenous_retroelement" if gap <= score_delta else "viral_top_hit_in_control"
            else:
                gap = 999.0
                reason = "no_retroelement_competitor_observed"
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "read_id": read_id,
                    "top_category": top["reference_category"].lower(),
                    "competing_categories": "|".join(sorted({row["reference_category"].lower() for row in competitors})),
                    "top_score": f"{top_score:.3f}",
                    "score_gap": f"{gap:.3f}",
                    "review_reason": reason,
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect likely false HIV/HTLV signal in healthy controls.")
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("hits_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--score-delta", type=float, default=5)
    args = parser.parse_args()
    inspect(args.manifest_csv, args.hits_csv, args.output_csv, args.score_delta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
