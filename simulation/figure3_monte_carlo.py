#!/usr/bin/env python3
"""Monte Carlo experiment for Section 3.2 and Figure 3 of the guide.

The experiment studies common standardization, complier turnover, and the
no-event continuation. It uses no external data.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIGURE_DIR = ROOT / "figures"
RESULT_DIR = ROOT / "results"

THETA = 1.0
PERIODS = np.array([-2, -1, 0, 1])
UNIVERSE_PROB_X = 0.50

SCENARIOS = {
    "stable": "Stable compliers,\nflat path",
    "observed_shift": "Observed X\nshift",
    "turnover": "Unobserved\nturnover",
    "secular": "Linear no-event\ntrend",
}

ESTIMATORS = {
    "raw_flat": "Unstandardized, flat",
    "standardized_flat": "Common-standardized, flat",
    "standardized_linear": "Common-standardized, linear",
}


def arm_difference(values: np.ndarray, z: np.ndarray) -> float:
    """Return the difference in sample means between instrument arms."""
    return values[z == 1].mean() - values[z == 0].mean()


def period_sample(
    scenario: str, period: int, sample_size: int, rng: np.random.Generator
) -> dict[str, np.ndarray | float]:
    """Draw one repeated cross section."""
    if scenario == "observed_shift":
        observation_prob_0, observation_prob_1 = (
            (0.90, 0.30) if period <= 0 else (0.30, 0.90)
        )
    else:
        observation_prob_0, observation_prob_1 = 0.60, 0.60

    prob_x = (
        UNIVERSE_PROB_X
        * observation_prob_1
        / (
            UNIVERSE_PROB_X * observation_prob_1
            + (1 - UNIVERSE_PROB_X) * observation_prob_0
        )
    )

    x = rng.binomial(1, prob_x, sample_size)
    z = rng.binomial(1, 0.50, sample_size)

    if scenario == "turnover":
        # 0: stable, 1: entrant, 2: exiter, 3: never-complier.
        response_type = rng.choice(4, sample_size, p=[0.30, 0.20, 0.20, 0.30])
        if period <= 0:
            complier = np.isin(response_type, [0, 2])
        else:
            complier = np.isin(response_type, [0, 1])
        base_effect = np.choose(response_type, [1.0, 2.0, 0.0, 1.0])
    else:
        complier = rng.binomial(1, 0.35, sample_size).astype(bool)
        base_effect = 0.5 + 1.5 * x

    if scenario == "secular":
        base_effect = base_effect + 0.4 * period

    treatment_effect = base_effect + THETA * (period == 1)
    d = z * complier.astype(int)
    untreated_outcome = 0.3 * x + 0.15 * period + rng.normal(0, 1, sample_size)
    y = untreated_outcome + treatment_effect * d

    conditional_wald = np.empty(2)
    conditional_first_stage = np.empty(2)
    empirical_x_mass = np.empty(2)
    for value_x in (0, 1):
        within_x = x == value_x
        conditional_first_stage[value_x] = arm_difference(d[within_x], z[within_x])
        reduced_form = arm_difference(y[within_x], z[within_x])
        conditional_wald[value_x] = reduced_form / conditional_first_stage[value_x]
        empirical_x_mass[value_x] = within_x.mean()

    raw_first_stage = arm_difference(d, z)
    raw_wald = arm_difference(y, z) / raw_first_stage

    return {
        "conditional_wald": conditional_wald,
        "conditional_first_stage": conditional_first_stage,
        "empirical_x_mass": empirical_x_mass,
        "raw_wald": raw_wald,
    }


def one_replication(
    scenario: str, sample_size: int, rng: np.random.Generator
) -> dict[str, float]:
    """Compute all three estimators in one replication."""
    period_data = {
        int(period): period_sample(scenario, int(period), sample_size, rng)
        for period in PERIODS
    }

    reference = period_data[0]
    reference_numerator = (
        reference["conditional_first_stage"] * reference["empirical_x_mass"]
    )
    reference_weights = reference_numerator / reference_numerator.sum()

    raw_effects = np.array([period_data[int(t)]["raw_wald"] for t in PERIODS])
    standardized_effects = np.array(
        [
            np.dot(reference_weights, period_data[int(t)]["conditional_wald"])
            for t in PERIODS
        ]
    )

    pre_design = np.column_stack([np.ones(3), PERIODS[:3]])
    pre_coefficients = np.linalg.lstsq(
        pre_design, standardized_effects[:3], rcond=None
    )[0]
    predicted_no_event_post = pre_coefficients @ np.array([1.0, 1.0])

    return {
        "raw_flat": raw_effects[3] - raw_effects[2],
        "standardized_flat": standardized_effects[3] - standardized_effects[2],
        "standardized_linear": standardized_effects[3] - predicted_no_event_post,
    }


def run_experiment(
    replications: int = 2_000, sample_size: int = 2_500, seed: int = 20260817
) -> pd.DataFrame:
    """Run and summarize the four design scenarios."""
    rng = np.random.default_rng(seed)
    records: list[dict[str, float | str | int]] = []
    for scenario in SCENARIOS:
        for replication in range(replications):
            estimates = one_replication(scenario, sample_size, rng)
            for estimator, estimate in estimates.items():
                records.append(
                    {
                        "scenario": scenario,
                        "scenario_label": SCENARIOS[scenario].replace("\n", " "),
                        "replication": replication,
                        "estimator": estimator,
                        "estimator_label": ESTIMATORS[estimator],
                        "estimate": estimate,
                    }
                )

    draws = pd.DataFrame.from_records(records)
    summary = (
        draws.groupby(
            ["scenario", "scenario_label", "estimator", "estimator_label"],
            sort=False,
        )["estimate"]
        .agg(
            mean="mean",
            standard_deviation="std",
            q025=lambda value: value.quantile(0.025),
            median="median",
            q975=lambda value: value.quantile(0.975),
        )
        .reset_index()
    )
    summary["bias"] = summary["mean"] - THETA
    rmse = (
        draws.assign(squared_error=lambda frame: (frame["estimate"] - THETA) ** 2)
        .groupby(["scenario", "estimator"], sort=False)["squared_error"]
        .mean()
        .pow(0.5)
        .rename("rmse")
        .reset_index()
    )
    summary = summary.merge(rmse, on=["scenario", "estimator"], how="left")
    summary["replications"] = replications
    summary["sample_size_per_period"] = sample_size
    return summary


def validate_summary(summary: pd.DataFrame) -> None:
    """Stop if simulated means depart materially from their population values."""
    expected_means = {
        ("stable", "raw_flat"): 1.0,
        ("stable", "standardized_flat"): 1.0,
        ("stable", "standardized_linear"): 1.0,
        ("observed_shift", "raw_flat"): 1.75,
        ("observed_shift", "standardized_flat"): 1.0,
        ("observed_shift", "standardized_linear"): 1.0,
        ("turnover", "raw_flat"): 1.8,
        ("turnover", "standardized_flat"): 1.8,
        ("turnover", "standardized_linear"): 1.8,
        ("secular", "raw_flat"): 1.4,
        ("secular", "standardized_flat"): 1.4,
        ("secular", "standardized_linear"): 1.0,
    }
    indexed = summary.set_index(["scenario", "estimator"])
    for key, expected in expected_means.items():
        simulated = indexed.loc[key, "mean"]
        if abs(simulated - expected) > 0.035:
            raise RuntimeError(
                f"Monte Carlo mean for {key} is {simulated:.4f}; "
                f"expected {expected:.4f}."
            )


def draw_figure(summary: pd.DataFrame) -> None:
    """Plot mean estimates across the four design scenarios."""
    colors = ["#595959", "#2878B5", "#E07A1F"]
    markers = ["o", "s", "^"]
    offsets = [-0.22, 0.0, 0.22]
    figure, axis = plt.subplots(figsize=(9.2, 5.3))
    scenario_order = list(SCENARIOS)

    for color, marker, offset, estimator in zip(
        colors, markers, offsets, ESTIMATORS, strict=True
    ):
        values = (
            summary[summary["estimator"] == estimator]
            .set_index("scenario")
            .loc[scenario_order]
        )
        x_locations = np.arange(len(scenario_order)) + offset
        axis.scatter(
            x_locations,
            values["mean"],
            s=72,
            marker=marker,
            color=color,
            label=ESTIMATORS[estimator],
            zorder=3,
        )
        axis.vlines(
            x_locations,
            values["q025"],
            values["q975"],
            color=color,
            linewidth=1.4,
            alpha=0.65,
            zorder=2,
        )

    axis.axhline(
        THETA,
        color="black",
        linewidth=1.2,
        linestyle="--",
        label="True change in LATE",
    )
    axis.set_xticks(np.arange(len(scenario_order)))
    axis.set_xticklabels([SCENARIOS[name] for name in scenario_order])
    axis.set_ylabel("Estimated contrast")
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, ncol=2, loc="upper left")
    axis.set_ylim(-0.2, 3.0)
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(
            FIGURE_DIR / f"figure3_standardization_and_continuation.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(figure)


def main() -> None:
    """Run the experiment and write its fixed-seed outputs."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    summary = run_experiment()
    validate_summary(summary)
    summary.to_csv(RESULT_DIR / "standardization_summary.csv", index=False)

    metadata = {
        "seed": 20260817,
        "replications": 2000,
        "sample_size_per_period": 2500,
        "periods": PERIODS.tolist(),
        "true_effect_change": THETA,
    }
    (RESULT_DIR / "figure3_monte_carlo_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    draw_figure(summary)
    print(summary.to_string(index=False, float_format="%.4f"))


if __name__ == "__main__":
    main()
