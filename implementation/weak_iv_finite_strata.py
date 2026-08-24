#!/usr/bin/env python3
"""Conservative weak-IV projection for a fixed finite-strata contrast.

This reference implementation matches Section E.3 of the guide. It assumes:

* one pre-event period and one post-event period;
* a flat no-event continuation;
* a fixed, finite set of strata;
* externally fixed reference weights; and
* valid asymptotically normal standard errors for each reduced form and first
  stage under the user's sampling and clustering design.

The procedure forms Bonferroni-simultaneous intervals for the primitive
moments, projects each reduced-form/first-stage rectangle through the ratio
map, and combines the resulting intervals with the fixed reference weights.
If a required first-stage interval contains zero, the final outer confidence
set is reported as the whole real line.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


Interval = Tuple[float, float]


REQUIRED_COLUMNS: Sequence[str] = (
    "stratum",
    "reference_weight",
    "rho_pre",
    "se_rho_pre",
    "pi_pre",
    "se_pi_pre",
    "rho_post",
    "se_rho_post",
    "pi_post",
    "se_pi_post",
)


def _finite_float(row: Mapping[str, str], name: str) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Column {name!r} must contain numeric values") from exc
    if not math.isfinite(value):
        raise ValueError(f"Column {name!r} must contain finite values")
    return value


def read_rows(path: Path) -> List[Dict[str, float | str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in REQUIRED_COLUMNS if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        rows: List[Dict[str, float | str]] = []
        for raw in reader:
            parsed: Dict[str, float | str] = {"stratum": raw["stratum"]}
            for name in REQUIRED_COLUMNS[1:]:
                parsed[name] = _finite_float(raw, name)
            for name in (
                "se_rho_pre",
                "se_pi_pre",
                "se_rho_post",
                "se_pi_post",
            ):
                if float(parsed[name]) < 0:
                    raise ValueError(f"{name} must be nonnegative")
            if float(parsed["reference_weight"]) <= 0:
                raise ValueError(
                    "reference_weight must be strictly positive; omit zero-weight rows"
                )
            rows.append(parsed)

    if not rows:
        raise ValueError("Input file contains no strata")

    total_weight = sum(float(row["reference_weight"]) for row in rows)
    if not math.isclose(total_weight, 1.0, rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError(f"reference_weight must sum to one; found {total_weight:.12g}")
    return rows


def normal_interval(estimate: float, standard_error: float, critical_value: float) -> Interval:
    radius = critical_value * standard_error
    return estimate - radius, estimate + radius


def ratio_interval(numerator: Interval, denominator: Interval) -> Interval:
    """Project a rectangle through n/d.

    Returning the whole real line when the denominator interval contains zero
    is deliberately conservative. It avoids silently turning weak first-stage
    information into a bounded ratio interval.
    """

    d_low, d_high = denominator
    if d_low <= 0.0 <= d_high:
        return -math.inf, math.inf

    n_low, n_high = numerator
    values = (
        n_low / d_low,
        n_low / d_high,
        n_high / d_low,
        n_high / d_high,
    )
    return min(values), max(values)


def subtract_intervals(left: Interval, right: Interval) -> Interval:
    if any(math.isinf(value) for value in (*left, *right)):
        return -math.inf, math.inf
    return left[0] - right[1], left[1] - right[0]


def weighted_sum(intervals: Iterable[Tuple[float, Interval]]) -> Interval:
    lower = 0.0
    upper = 0.0
    for weight, interval in intervals:
        if weight == 0:
            continue
        if any(math.isinf(value) for value in interval):
            return -math.inf, math.inf
        lower += weight * interval[0]
        upper += weight * interval[1]
    return lower, upper


def _json_bound(value: float) -> float | None:
    return None if math.isinf(value) else value


def project(rows: Sequence[Mapping[str, float | str]], alpha: float) -> Dict[str, object]:
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie strictly between zero and one")
    if not rows:
        raise ValueError("Input contains no strata")

    reference_weights = [float(row["reference_weight"]) for row in rows]
    if any(not math.isfinite(weight) or weight <= 0.0 for weight in reference_weights):
        raise ValueError(
            "reference_weight must be finite and strictly positive; omit zero-weight rows"
        )
    total_weight = sum(reference_weights)
    if not math.isclose(total_weight, 1.0, rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError(f"reference_weight must sum to one; found {total_weight:.12g}")

    primitive_count = 4 * len(rows)
    critical_value = NormalDist().inv_cdf(1.0 - alpha / (2.0 * primitive_count))
    stratum_results: List[Dict[str, object]] = []
    weighted_changes: List[Tuple[float, Interval]] = []
    first_stage_intervals: List[Interval] = []
    point_estimate = 0.0

    required_first_stages = [
        float(row[f"pi_{period}"])
        for row in rows
        for period in ("pre", "post")
    ]
    if all(value > 0.0 for value in required_first_stages):
        point_orientation_status = "positive"
    elif all(value < 0.0 for value in required_first_stages):
        point_orientation_status = "negative"
    else:
        point_orientation_status = "mixed_or_zero"

    for row in rows:
        intervals: Dict[str, Interval] = {}
        for period in ("pre", "post"):
            intervals[f"rho_{period}"] = normal_interval(
                float(row[f"rho_{period}"]),
                float(row[f"se_rho_{period}"]),
                critical_value,
            )
            intervals[f"pi_{period}"] = normal_interval(
                float(row[f"pi_{period}"]),
                float(row[f"se_pi_{period}"]),
                critical_value,
            )
            first_stage_intervals.append(intervals[f"pi_{period}"])

        ratio_pre = ratio_interval(intervals["rho_pre"], intervals["pi_pre"])
        ratio_post = ratio_interval(intervals["rho_post"], intervals["pi_post"])
        change = subtract_intervals(ratio_post, ratio_pre)
        weight = float(row["reference_weight"])
        weighted_changes.append((weight, change))

        pi_pre = float(row["pi_pre"])
        pi_post = float(row["pi_post"])
        if pi_pre == 0.0 or pi_post == 0.0:
            point_estimate = math.nan
        elif not math.isnan(point_estimate):
            point_estimate += weight * (
                float(row["rho_post"]) / pi_post - float(row["rho_pre"]) / pi_pre
            )

        stratum_results.append(
            {
                "stratum": str(row["stratum"]),
                "reference_weight": weight,
                "primitive_intervals": {
                    name: [_json_bound(bounds[0]), _json_bound(bounds[1])]
                    for name, bounds in intervals.items()
                },
                "ratio_pre_interval": [_json_bound(ratio_pre[0]), _json_bound(ratio_pre[1])],
                "ratio_post_interval": [_json_bound(ratio_post[0]), _json_bound(ratio_post[1])],
                "change_interval": [_json_bound(change[0]), _json_bound(change[1])],
            }
        )

    confidence_set = weighted_sum(weighted_changes)
    all_real = math.isinf(confidence_set[0]) or math.isinf(confidence_set[1])
    if any(bounds[1] <= 0.0 for bounds in first_stage_intervals):
        interval_assessment = "contradicts_declared_positive_orientation"
        causal_interpretation_stop = True
        orientation_message = (
            "At least one simultaneous first-stage interval lies wholly at or below zero, "
            "contradicting the declared positive orientation at this confidence level. Stop "
            "the common-orientation causal interpretation unless the design is redefined and "
            "defended; the script does not reorient the data."
        )
    elif all(bounds[0] > 0.0 for bounds in first_stage_intervals):
        interval_assessment = "consistent_with_declared_positive_orientation"
        causal_interpretation_stop = False
        orientation_message = (
            "Every simultaneous first-stage interval lies above zero under the declared "
            "positive orientation. This diagnostic does not establish monotonicity or "
            "exclusion."
        )
    else:
        interval_assessment = "inconclusive_because_an_interval_spans_zero"
        causal_interpretation_stop = False
        orientation_message = (
            "No required simultaneous first-stage interval lies wholly below zero, but at "
            "least one spans zero with a positive upper endpoint. Treat orientation evidence "
            "as inconclusive and weak "
            "identification as potentially severe; sampling noise alone is not a causal-"
            "assumption failure."
        )
    return {
        "alpha": alpha,
        "coverage_level": 1.0 - alpha,
        "primitive_moment_count": primitive_count,
        "bonferroni_critical_value": critical_value,
        "point_estimate": None if math.isnan(point_estimate) else point_estimate,
        "confidence_set": {
            "lower": _json_bound(confidence_set[0]),
            "upper": _json_bound(confidence_set[1]),
            "all_real": all_real,
            "display": "(-inf, inf)" if all_real else f"[{confidence_set[0]:.6g}, {confidence_set[1]:.6g}]",
        },
        "first_stage_orientation": {
            "point_estimate_status": point_orientation_status,
            "simultaneous_interval_assessment": interval_assessment,
            "causal_interpretation_stop": causal_interpretation_stop,
            "message": orientation_message,
        },
        "strata": stratum_results,
        "scope": (
            "Fixed finite strata, fixed reference weights, one pre period, one post period, "
            "flat continuation. Primitive standard errors and asymptotic normal approximations "
            "must already match the sampling and clustering design; the normal critical value "
            "is not a few-cluster correction."
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path, help="CSV containing the required primitive moments")
    parser.add_argument("--alpha", type=float, default=0.05, help="Tail probability (default: 0.05)")
    parser.add_argument("--output", type=Path, help="Optional JSON output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = project(read_rows(args.input_csv), args.alpha)
    rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
