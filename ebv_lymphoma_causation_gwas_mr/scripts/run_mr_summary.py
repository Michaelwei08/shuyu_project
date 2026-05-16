from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def weighted_median(ratios: list[tuple[float, float]]) -> float:
    ordered = sorted(ratios, key=lambda item: item[0])
    total_weight = sum(weight for _, weight in ordered)
    cumulative = 0.0
    for ratio, weight in ordered:
        cumulative += weight
        if cumulative >= total_weight / 2:
            return ratio
    return ordered[-1][0]


def normal_p(z_score: float) -> float:
    return math.erfc(abs(z_score) / math.sqrt(2))


def weighted_regression(xs: list[float], ys: list[float], ws: list[float], intercept: bool) -> tuple[float, float, float, float, float, float]:
    if not intercept:
        denom = sum(w * x * x for x, w in zip(xs, ws))
        slope = sum(w * x * y for x, y, w in zip(xs, ys, ws)) / denom
        se = math.sqrt(1 / denom) if denom > 0 else float("nan")
        p_value = normal_p(slope / se) if se > 0 else float("nan")
        return 0.0, slope, se, p_value, float("nan"), float("nan")
    sw = sum(ws)
    xbar = sum(w * x for x, w in zip(xs, ws)) / sw
    ybar = sum(w * y for y, w in zip(ys, ws)) / sw
    denom = sum(w * (x - xbar) ** 2 for x, w in zip(xs, ws))
    slope = sum(w * (x - xbar) * (y - ybar) for x, y, w in zip(xs, ys, ws)) / denom
    intercept_value = ybar - slope * xbar
    residuals = [y - (intercept_value + slope * x) for x, y in zip(xs, ys)]
    dof = max(len(xs) - 2, 1)
    residual_scale = sum(w * residual * residual for residual, w in zip(residuals, ws)) / dof
    slope_se = math.sqrt(residual_scale / denom) if denom > 0 else float("nan")
    p_value = normal_p(slope / slope_se) if slope_se > 0 else float("nan")
    intercept_se = math.sqrt(residual_scale * (1 / sw + (xbar**2 / denom))) if denom > 0 and sw > 0 else float("nan")
    intercept_p = normal_p(intercept_value / intercept_se) if intercept_se > 0 else float("nan")
    return intercept_value, slope, slope_se, p_value, intercept_se, intercept_p


def run(input_csv: Path, output_csv: Path) -> None:
    bx: list[float] = []
    by: list[float] = []
    weights: list[float] = []
    ratios: list[tuple[float, float]] = []
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"beta_exposure", "se_exposure", "beta_outcome", "se_outcome"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            exposure = float(row["beta_exposure"])
            outcome = float(row["beta_outcome"])
            se_outcome = float(row["se_outcome"])
            if exposure == 0 or se_outcome <= 0:
                continue
            weight = 1 / (se_outcome**2)
            bx.append(exposure)
            by.append(outcome)
            weights.append(weight)
            ratios.append((outcome / exposure, weight))

    ivw_intercept, ivw_slope, ivw_se, ivw_p, ivw_intercept_se, ivw_intercept_p = weighted_regression(bx, by, weights, intercept=False)
    egger_intercept, egger_slope, egger_se, egger_p, egger_intercept_se, egger_intercept_p = weighted_regression(bx, by, weights, intercept=True)
    median = weighted_median(ratios)
    # This is a lightweight placeholder summary for the ConMix slot; final analysis should use a validated MR package.
    conmix_proxy = median

    with output_csv.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "method",
                "estimate",
                "standard_error",
                "p_value",
                "intercept",
                "intercept_standard_error",
                "intercept_p_value",
                "n_instruments",
                "note",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "method": "IVW",
                "estimate": f"{ivw_slope:.6f}",
                "standard_error": f"{ivw_se:.6f}" if not math.isnan(ivw_se) else "",
                "p_value": f"{ivw_p:.6g}" if not math.isnan(ivw_p) else "",
                "intercept": f"{ivw_intercept:.6f}",
                "intercept_standard_error": f"{ivw_intercept_se:.6f}" if not math.isnan(ivw_intercept_se) else "",
                "intercept_p_value": f"{ivw_intercept_p:.6g}" if not math.isnan(ivw_intercept_p) else "",
                "n_instruments": len(bx),
                "note": "weighted regression through origin; lightweight normal approximation",
            }
        )
        writer.writerow(
            {
                "method": "MR-Egger",
                "estimate": f"{egger_slope:.6f}",
                "standard_error": f"{egger_se:.6f}" if not math.isnan(egger_se) else "",
                "p_value": f"{egger_p:.6g}" if not math.isnan(egger_p) else "",
                "intercept": f"{egger_intercept:.6f}",
                "intercept_standard_error": f"{egger_intercept_se:.6f}" if not math.isnan(egger_intercept_se) else "",
                "intercept_p_value": f"{egger_intercept_p:.6g}" if not math.isnan(egger_intercept_p) else "",
                "n_instruments": len(bx),
                "note": "weighted regression with intercept; lightweight normal approximation",
            }
        )
        writer.writerow({"method": "Weighted median", "estimate": f"{median:.6f}", "standard_error": "", "p_value": "", "intercept": "", "intercept_standard_error": "", "intercept_p_value": "", "n_instruments": len(bx), "note": "weighted median of ratio estimates; no scaffold SE"})
        writer.writerow({"method": "ConMix", "estimate": f"{conmix_proxy:.6f}", "standard_error": "", "p_value": "", "intercept": "", "intercept_standard_error": "", "intercept_p_value": "", "n_instruments": len(bx), "note": "placeholder proxy; use validated ConMix implementation for final analysis"})


def main() -> int:
    parser = argparse.ArgumentParser(description="Run lightweight MR method summaries on harmonized instruments.")
    parser.add_argument("harmonized_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    run(args.harmonized_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
