from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


DISEASES = ["HL", "PTLD", "DLBCL"]
VIRAL_PHENOTYPES = ["ebv_load", "anellovirus_load"]
SUBTYPE_PHENOTYPES = ["v1_pan_viral", "v3_high_ebv_low_anello"]


def as_number(value: str) -> float | None:
    clean = value.strip()
    if clean in {"", "NA", "N/A", "unknown"}:
        return None
    if clean.lower() == "true":
        return 1.0
    if clean.lower() == "false":
        return 0.0
    return float(clean)


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def variance(values: list[float]) -> float:
    if len(values) < 2:
        return float("nan")
    center = mean(values)
    return sum((value - center) ** 2 for value in values) / (len(values) - 1)


def ci_from_groups(positive: list[float], reference: list[float]) -> tuple[float, float, float, float]:
    diff = mean(positive) - mean(reference)
    var_pos = variance(positive)
    var_ref = variance(reference)
    se = math.sqrt(var_pos / len(positive) + var_ref / len(reference)) if len(positive) > 1 and len(reference) > 1 else float("nan")
    if math.isnan(se):
        return diff, se, float("nan"), float("nan")
    return diff, se, diff - 1.96 * se, diff + 1.96 * se


def solve_linear_system(matrix: list[list[float]], vector: list[float]) -> list[float] | None:
    n = len(vector)
    aug = [row[:] + [value] for row, value in zip(matrix, vector)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if abs(aug[pivot][col]) < 1e-12:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [value / scale for value in aug[col]]
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            aug[row] = [value - factor * aug[col][idx] for idx, value in enumerate(aug[row])]
    return [row[-1] for row in aug]


def invert_matrix(matrix: list[list[float]]) -> list[list[float]] | None:
    n = len(matrix)
    inverse: list[list[float]] = []
    for idx in range(n):
        unit = [0.0] * n
        unit[idx] = 1.0
        solution = solve_linear_system(matrix, unit)
        if solution is None:
            return None
        inverse.append(solution)
    return [[inverse[col][row] for col in range(n)] for row in range(n)]


def adjusted_effect(rows: list[dict[str, str]], phenotype: str, positive_column: str, adjust: list[str]) -> tuple[float, float, float, float] | None:
    design: list[list[float]] = []
    y: list[float] = []
    batches = sorted({row.get("batch", "") for row in rows if row.get("batch", "")})
    batch_levels = batches[1:]
    for row in rows:
        value = as_number(row.get(phenotype, ""))
        positive = as_number(row.get(positive_column, ""))
        if value is None or positive is None:
            continue
        covariates = [1.0, positive]
        if "sequencing_depth" in adjust:
            depth = as_number(row.get("sequencing_depth", ""))
            if depth is None:
                continue
            covariates.append(depth)
        if "batch" in adjust:
            covariates.extend(1.0 if row.get("batch", "") == level else 0.0 for level in batch_levels)
        design.append(covariates)
        y.append(value)
    if len(design) <= len(design[0]):
        return None
    xtx = [[sum(row[i] * row[j] for row in design) for j in range(len(design[0]))] for i in range(len(design[0]))]
    xty = [sum(row[i] * value for row, value in zip(design, y)) for i in range(len(design[0]))]
    beta = solve_linear_system(xtx, xty)
    inv = invert_matrix(xtx)
    if beta is None or inv is None:
        return None
    residuals = [value - sum(coef * x for coef, x in zip(beta, row)) for row, value in zip(design, y)]
    dof = max(len(y) - len(beta), 1)
    sigma2 = sum(residual * residual for residual in residuals) / dof
    se = math.sqrt(max(sigma2 * inv[1][1], 0.0))
    effect = beta[1]
    return effect, se, effect - 1.96 * se, effect + 1.96 * se


def disease_rows(rows: list[dict[str, str]], disease: str) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for row in rows:
        if row["disease_label"] == disease:
            output = dict(row)
            output["is_positive"] = "1"
            selected.append(output)
        elif row["disease_label"] == "control":
            output = dict(row)
            output["is_positive"] = "0"
            selected.append(output)
    return selected


def write_result(writer: csv.DictWriter, test_name: str, model: str, phenotype: str, group: str, result: tuple[float, float, float, float] | None, n_positive: int, n_reference: int) -> None:
    if result is None:
        effect = se = ci_low = ci_high = float("nan")
    else:
        effect, se, ci_low, ci_high = result
    writer.writerow(
        {
            "test_name": test_name,
            "model": model,
            "phenotype": phenotype,
            "group": group,
            "n_positive": n_positive,
            "n_reference": n_reference,
            "effect_estimate": f"{effect:.6f}" if not math.isnan(effect) else "",
            "standard_error": f"{se:.6f}" if not math.isnan(se) else "",
            "ci95_low": f"{ci_low:.6f}" if not math.isnan(ci_low) else "",
            "ci95_high": f"{ci_high:.6f}" if not math.isnan(ci_high) else "",
            "note": "scaffold linear effect; replace example input with real cohort data",
        }
    )


def run(input_csv: Path, output_csv: Path) -> None:
    with input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    fields = ["test_name", "model", "phenotype", "group", "n_positive", "n_reference", "effect_estimate", "standard_error", "ci95_low", "ci95_high", "note"]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for disease in DISEASES:
            selected = disease_rows(rows, disease)
            n_positive = sum(1 for row in selected if row["is_positive"] == "1")
            n_reference = sum(1 for row in selected if row["is_positive"] == "0")
            for phenotype in [*VIRAL_PHENOTYPES, *SUBTYPE_PHENOTYPES]:
                positive = [as_number(row[phenotype]) for row in selected if row["is_positive"] == "1"]
                reference = [as_number(row[phenotype]) for row in selected if row["is_positive"] == "0"]
                positive_values = [value for value in positive if value is not None]
                reference_values = [value for value in reference if value is not None]
                write_result(writer, "disease_vs_control", "unadjusted", phenotype, disease, ci_from_groups(positive_values, reference_values), n_positive, n_reference)
                write_result(writer, "disease_vs_control", "depth_adjusted", phenotype, disease, adjusted_effect(selected, phenotype, "is_positive", ["sequencing_depth"]), n_positive, n_reference)
                write_result(writer, "disease_vs_control", "depth_batch_adjusted", phenotype, disease, adjusted_effect(selected, phenotype, "is_positive", ["sequencing_depth", "batch"]), n_positive, n_reference)

        disease_only = [row for row in rows if row["disease_label"] != "control"]
        for row in disease_only:
            row["is_adverse"] = "1" if row.get("clinical_outcome") == "adverse" else "0"
        n_adverse = sum(1 for row in disease_only if row["is_adverse"] == "1")
        n_non_adverse = len(disease_only) - n_adverse
        for phenotype in SUBTYPE_PHENOTYPES:
            positive = [as_number(row[phenotype]) for row in disease_only if row["is_adverse"] == "1"]
            reference = [as_number(row[phenotype]) for row in disease_only if row["is_adverse"] == "0"]
            write_result(writer, "clinical_outcome", "unadjusted", phenotype, "all_diseases", ci_from_groups([v for v in positive if v is not None], [v for v in reference if v is not None]), n_adverse, n_non_adverse)

        for disease in DISEASES:
            within = [row for row in disease_only if row["disease_label"] == disease]
            n_pos = sum(1 for row in within if row["is_adverse"] == "1")
            n_ref = len(within) - n_pos
            for phenotype in VIRAL_PHENOTYPES:
                if n_pos == 0 or n_ref == 0:
                    write_result(writer, "within_disease_outcome", "unadjusted", phenotype, disease, None, n_pos, n_ref)
                    continue
                positive = [as_number(row[phenotype]) for row in within if row["is_adverse"] == "1"]
                reference = [as_number(row[phenotype]) for row in within if row["is_adverse"] == "0"]
                write_result(writer, "within_disease_outcome", "unadjusted", phenotype, disease, ci_from_groups([v for v in positive if v is not None], [v for v in reference if v is not None]), n_pos, n_ref)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run scaffold clinical association grids for viral phenotypes.")
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    run(args.input_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
