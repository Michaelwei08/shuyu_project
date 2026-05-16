from __future__ import annotations

import argparse
import csv
import html
import re
from dataclasses import dataclass
from pathlib import Path


LARGE_GB = 40.0


@dataclass(frozen=True)
class DatasetRules:
    target_disease_match: str
    sample_type: str
    expected_virus_signal: str
    clinical_association_ready: str
    primary_use: str


RULES = {
    "PRJNA809536": DatasetRules(
        "partial",
        "patient_blood_or_pbmc",
        "EBV_positive_HL_context",
        "no_small_case_report",
        "HL_EBV_positive_control",
    ),
    "PRJNA172563": DatasetRules(
        "partial",
        "dlbcl_cell_line",
        "unknown_or_cell_line_viral_signal",
        "no_cell_lines",
        "DLBCL_WGS_method_smoke_test",
    ),
    "PRJNA772081": DatasetRules(
        "partial",
        "immunodeficiency_lymphoma_patient",
        "possible_EBV_or_anellovirus_context",
        "limited_small_mixed_modalities",
        "immunodeficiency_lymphoma_bridge",
    ),
    "PRJDB11215": DatasetRules(
        "outside_p2_yes_p1",
        "htlv1_infected_cell_line",
        "HTLV1_expected",
        "no_cell_lines",
        "HTLV_discrimination_control",
    ),
    "PRJNA1212064": DatasetRules(
        "partial",
        "ebv_positive_burkitt_cell_line",
        "EBV_expected",
        "no_cell_line_single_run",
        "EBV_mapping_positive_control",
    ),
    "PRJNA998843": DatasetRules(
        "partial",
        "hhv8_positive_pel_cell_line",
        "HHV8_expected",
        "no_cell_line_single_project",
        "herpesvirus_mapping_control",
    ),
    "PRJNA528941": DatasetRules(
        "partial",
        "raji_ebv_positive_cell_line",
        "EBV_expected",
        "no_cell_line_single_run",
        "EBV_mapping_positive_control",
    ),
    "PRJDB18873": DatasetRules(
        "no_host_wgs",
        "ebv_viral_genome",
        "EBV_expected",
        "no_viral_only",
        "EBV_reference_validation",
    ),
    "PRJNA551423": DatasetRules(
        "ptld_but_no_host_wgs",
        "ebv_amplicon",
        "EBV_expected_PTLD_context",
        "no_amplicon_only",
        "PTLD_EBV_variant_reference",
    ),
}


def normalize_project(value: str) -> str:
    match = re.search(r"(PRJ[A-Z]{0,2}\d+)", value)
    return match.group(1) if match else value.strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def runinfo_path(runinfo_dir: Path, project: str) -> Path | None:
    candidates = sorted(runinfo_dir.glob(f"ncbi_runinfo_{project}_direct_download_*.csv"))
    if not candidates:
        candidates = sorted(runinfo_dir.glob(f"ncbi_runinfo_{project}_2026-*.csv"))
    return candidates[-1] if candidates else None


def is_direct(row: dict[str, str]) -> bool:
    return row.get("download_path", "").startswith("http")


def is_wgs(row: dict[str, str]) -> bool:
    return row.get("LibraryStrategy", "").strip().upper() == "WGS"


def size_gb(row: dict[str, str]) -> float:
    try:
        return float(row.get("size_MB") or 0.0) / 1024.0
    except ValueError:
        return 0.0


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def score_dataset(
    project: str,
    rules: DatasetRules,
    direct_wgs_runs: int,
    large_direct_wgs_runs: int,
    direct_total_gb: float,
) -> tuple[int, str]:
    score = 0
    if direct_wgs_runs > 0:
        score += 1
    if large_direct_wgs_runs > 0 or direct_total_gb >= 100:
        score += 1
    if "HL" in rules.primary_use or "DLBCL" in rules.primary_use or "immunodeficiency_lymphoma" in rules.primary_use:
        score += 1
    if rules.expected_virus_signal.startswith(("EBV", "HHV8", "HTLV1")):
        score += 1
    if rules.clinical_association_ready.startswith("yes"):
        score += 1

    if project == "PRJDB18873":
        return 2, "usable_for_EBV_reference_only"
    if project == "PRJNA551423":
        return 1, "usable_for_PTLD_EBV_variant_only"
    if score >= 4:
        return score, "supports_smoke_test_not_full_hypothesis"
    if score >= 3:
        return score, "partial_support"
    return score, "weak_support"


def build_summary(candidate_csv: Path, runinfo_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for candidate in read_csv(candidate_csv):
        project = normalize_project(candidate["accession_or_project"])
        rules = RULES.get(
            project,
            DatasetRules("unknown", "unknown", "unknown", "unknown", "unclassified"),
        )
        path = runinfo_path(runinfo_dir, project)
        runinfo = read_csv(path) if path else []
        direct_rows = [row for row in runinfo if is_direct(row)]
        direct_wgs = [row for row in direct_rows if is_wgs(row)]
        direct_non_wgs = [row for row in direct_rows if not is_wgs(row)]
        large_direct_wgs = [row for row in direct_wgs if size_gb(row) >= LARGE_GB]
        direct_total_gb = sum(size_gb(row) for row in direct_wgs)
        score, verdict = score_dataset(
            project,
            rules,
            len(direct_wgs),
            len(large_direct_wgs),
            direct_total_gb,
        )
        rows.append(
            {
                "project": project,
                "dataset_name": candidate["dataset_name"],
                "runinfo_found": yes_no(path is not None),
                "manifest_direct_wgs_runs": candidate.get("direct_raw_wgs_runs", ""),
                "observed_direct_wgs_runs": len(direct_wgs),
                "observed_direct_non_wgs_runs": len(direct_non_wgs),
                "observed_direct_wgs_ge40gb": len(large_direct_wgs),
                "observed_direct_wgs_total_gb": round(direct_total_gb, 1),
                "target_disease_match": rules.target_disease_match,
                "sample_type": rules.sample_type,
                "expected_virus_signal": rules.expected_virus_signal,
                "clinical_association_ready": rules.clinical_association_ready,
                "primary_use": rules.primary_use,
                "hypothesis_score_0_5": score,
                "verdict": verdict,
            }
        )
    return rows


def build_hypotheses(summary: list[dict[str, object]]) -> list[dict[str, object]]:
    def any_row(predicate) -> bool:
        return any(predicate(row) for row in summary)

    return [
        {
            "hypothesis": "Direct raw WGS exists for a target HL/DLBCL/PTLD disease context",
            "status": "partial",
            "evidence": "PRJNA809536 supports HL/EBV but is a 3-run case-report-scale dataset; PRJNA172563 supports DLBCL/NHL but mainly cell lines.",
        },
        {
            "hypothesis": "Direct raw WGS is large enough for detector smoke testing",
            "status": "pass" if any_row(lambda row: int(row["observed_direct_wgs_ge40gb"]) > 0) else "fail",
            "evidence": "At least one retained dataset has direct WGS runs >=40 GB; PRJNA172563 has 12 and PRJNA809536 has 2.",
        },
        {
            "hypothesis": "Direct data can validate EBV/herpesvirus positive controls",
            "status": "pass",
            "evidence": "HL EBV context, Raji/EB-3 EBV-positive cell lines, BC-3 HHV8-positive PEL, and EBV viral genomes are directly downloadable.",
        },
        {
            "hypothesis": "Direct data can support anellovirus detection in lymphoma WGS",
            "status": "weak",
            "evidence": "No direct dataset has known anellovirus-positive lymphoma WGS annotation; anellovirus can only be explored opportunistically in host WGS.",
        },
        {
            "hypothesis": "Direct data can support clinical association replication",
            "status": "fail",
            "evidence": "Direct datasets are too small, cell-line-heavy, or viral-only; controlled datasets are still needed for clinical association replication.",
        },
        {
            "hypothesis": "Direct data can support HTLV/HIV/retroviral discrimination controls",
            "status": "partial",
            "evidence": "PRJDB11215 provides direct HTLV-1-infected cell-line WGS; no direct HIV cohort WGS with healthy controls was retained.",
        },
    ]


def render_svg(summary: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    width = 1100
    height = 620
    left = 260
    top = 80
    row_h = 46
    bar_w = 500
    max_gb = max(float(row["observed_direct_wgs_total_gb"]) for row in summary) or 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="1100" height="620" fill="#f8fafc"/>',
        '<text x="550" y="38" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="24" font-weight="700" fill="#0f172a">Direct-Download Dataset Smoke Test</text>',
        '<text x="550" y="62" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="13" fill="#475569">Bar length = direct WGS GB; right badge = hypothesis score out of 5</text>',
    ]
    for idx, row in enumerate(summary):
        y = top + idx * row_h
        gb = float(row["observed_direct_wgs_total_gb"])
        score = int(row["hypothesis_score_0_5"])
        fill = "#0f766e" if score >= 4 else "#d97706" if score >= 3 else "#64748b"
        label = f"{row['project']} ({row['primary_use']})"
        label = label[:45] + "..." if len(label) > 48 else label
        current_bar = max(2.0, (gb / max_gb) * bar_w)
        parts.extend(
            [
                f'<text x="{left - 14}" y="{y + 25}" text-anchor="end" font-family="Segoe UI, Arial, sans-serif" font-size="12" fill="#0f172a">{html.escape(label)}</text>',
                f'<rect x="{left}" y="{y + 8}" width="{bar_w}" height="22" rx="4" fill="#e2e8f0"/>',
                f'<rect x="{left}" y="{y + 8}" width="{current_bar:.1f}" height="22" rx="4" fill="{fill}"/>',
                f'<text x="{left + bar_w + 12}" y="{y + 24}" font-family="Segoe UI, Arial, sans-serif" font-size="12" fill="#334155">{gb:g} GB; {row["observed_direct_wgs_runs"]} WGS runs; {row["observed_direct_wgs_ge40gb"]} >=40GB</text>',
                f'<rect x="980" y="{y + 7}" width="58" height="24" rx="12" fill="{fill}"/>',
                f'<text x="1009" y="{y + 24}" text-anchor="middle" font-family="Segoe UI, Arial, sans-serif" font-size="12" font-weight="700" fill="#ffffff">{score}/5</text>',
            ]
        )
    parts.append("</svg>")
    output.write_text("".join(parts), encoding="utf-8")


def write_report(summary: list[dict[str, object]], hypotheses: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    best = sorted(summary, key=lambda row: (-int(row["hypothesis_score_0_5"]), -float(row["observed_direct_wgs_total_gb"])))
    lines = [
        "# Direct-Download Smoke Test",
        "",
        "## Bottom Line",
        "",
        "The directly downloadable datasets support method smoke testing and viral-positive controls, but they do not fully support the original clinical-cohort hypothesis. The missing piece is still a directly downloadable large patient-level HL/DLBCL/PTLD WGS cohort with usable clinical metadata.",
        "",
        "## Hypothesis Matrix",
        "",
        "| Hypothesis | Status | Evidence |",
        "|---|---:|---|",
    ]
    for row in hypotheses:
        lines.append(f"| {row['hypothesis']} | {row['status']} | {row['evidence']} |")
    lines.extend(["", "## Top Dataset Verdicts", ""])
    for row in best:
        lines.append(
            f"- `{row['project']}`: {row['verdict']}; {row['observed_direct_wgs_runs']} direct WGS runs, "
            f"{row['observed_direct_wgs_ge40gb']} >=40 GB, {row['observed_direct_wgs_total_gb']} GB total. "
            f"Primary use: `{row['primary_use']}`."
        )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Smoke-test direct-download candidate datasets using RunInfo metadata.")
    parser.add_argument("candidate_csv", type=Path)
    parser.add_argument("runinfo_dir", type=Path)
    parser.add_argument("results_csv", type=Path)
    parser.add_argument("hypothesis_csv", type=Path)
    parser.add_argument("report_md", type=Path)
    parser.add_argument("plot_svg", type=Path)
    args = parser.parse_args()

    summary = build_summary(args.candidate_csv, args.runinfo_dir)
    hypotheses = build_hypotheses(summary)
    write_csv(args.results_csv, summary, list(summary[0].keys()))
    write_csv(args.hypothesis_csv, hypotheses, ["hypothesis", "status", "evidence"])
    write_report(summary, hypotheses, args.report_md)
    render_svg(summary, args.plot_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
