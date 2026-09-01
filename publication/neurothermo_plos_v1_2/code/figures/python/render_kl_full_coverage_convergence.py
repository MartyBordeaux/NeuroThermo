#!/usr/bin/env python3
"""Render FigS3 from the frozen full-coverage KL convergence summaries."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo-mpl")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "kl_convergence_v1_0_1"
OUT = ROOT / "results" / "figures"


def main():
    ensemble = pd.read_csv(DATA / "ensemble_convergence_summary.csv")
    scenario = pd.read_csv(DATA / "scenario_convergence_gates.csv")
    cells = pd.read_csv(DATA / "cell_balanced_deltas.csv")

    variant = "seed_median_curve_isotonic"
    main = ensemble[(ensemble.marker_variant == variant) & ensemble.view.isin(["xyz", "xy"])]
    colors = {"xyz": "#54278f", "xy": "#e08214"}
    styles = {"marker_first": "-", "curve_first": "--"}

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), constrained_layout=True)

    for (view, order), group in main.groupby(["view", "aggregation_order"]):
        group = group.sort_values("dt_ms", ascending=False)
        x = np.arange(len(group))
        axes[0].plot(
            x,
            group.median_pair_delta,
            marker="o",
            color=colors[view],
            linestyle=styles[order],
            label=f"{view}, {order.replace('_', ' ')}",
        )
    axes[0].axhline(0, color="black", linewidth=0.7)
    labels = ["0.05", "0.025", "0.0125"]
    axes[0].set_xticks(range(len(labels)), labels)
    axes[0].set(xlabel="Integration step dt (ms)", ylabel=r"Median $p_{KL}-p_{firing}$")
    axes[0].legend(fontsize=7, frameon=False)

    pass_rates = scenario.groupby("view")["pass"].mean().reindex(["xyz", "xy", "z"])
    axes[1].bar(pass_rates.index, pass_rates.values, color=[colors["xyz"], colors["xy"], "#777777"])
    axes[1].axhline(0.80, color="black", linestyle="--", linewidth=0.8)
    axes[1].set(ylim=(0, 1), ylabel="Scenario convergence pass fraction", xlabel="Distributional view")

    fine = min(cells.dt_ms)
    selected = cells[
        (cells.dt_ms == fine)
        & (cells.marker_variant == variant)
        & (cells.aggregation_order == "marker_first")
    ]
    arrays = [
        selected.loc[selected.view == view, "kl_minus_firing_p"].to_numpy()
        for view in ("xyz", "xy")
    ]
    axes[2].boxplot(arrays, tick_labels=["xyz", "xy"], showfliers=False)
    axes[2].axhline(0, color="black", linewidth=0.7)
    axes[2].set(xlabel="Distributional view", ylabel=r"Cell-balanced $p_{KL}-p_{firing}$")

    for label, axis in zip("ABC", axes):
        axis.text(-0.13, 1.04, label, transform=axis.transAxes, fontweight="bold", fontsize=12)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "FigS3_KL_full_coverage_convergence.pdf", bbox_inches="tight")
    fig.savefig(OUT / "FigS3_KL_full_coverage_convergence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
