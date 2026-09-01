#!/usr/bin/env python3
"""Render marker-definition robustness for the full-coverage KL analysis."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo-mpl")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data"
OUTPUTS = (ROOT / "results" / "figures",)

COLORS = {"xyz": "#2B6CB0", "xy": "#C44E52", "z": "#3C8D73"}
VARIANTS = [
    ("seed_median_curve_isotonic", "Median curve\nisotonic"),
    ("seed_median_curve_first", "Median curve\nfirst"),
    ("seed_median_curve_persistent", "Median curve\npersistent"),
    ("median_seed_isotonic", "Median seed\nisotonic"),
    ("q25_seed_isotonic", "Q25 seed\nisotonic"),
    ("q75_seed_isotonic", "Q75 seed\nisotonic"),
]


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "legend.fontsize": 8,
        "savefig.facecolor": "white",
    })


def save(fig, stem):
    for directory in OUTPUTS:
        directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(directory / f"{stem}.pdf", bbox_inches="tight")
        fig.savefig(directory / f"{stem}.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


def marker_robustness():
    data = pd.read_csv(
        DATA / "kl_convergence_v1_0_1" / "ensemble_convergence_summary.csv"
    )
    data = data[(data.dt_ms == 0.025) & data.view.isin(["xyz", "xy"])]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7), sharey=True,
                             constrained_layout=True)
    x = np.arange(len(VARIANTS), dtype=float)
    for ax, view in zip(axes, ("xyz", "xy")):
        for aggregation, offset, linestyle, label in [
            ("marker_first", -0.08, "-", "marker first"),
            ("curve_first", 0.08, "--", "curve first"),
        ]:
            rows = data[(data.view == view)
                        & (data.aggregation_order == aggregation)].set_index("marker_variant")
            med = np.array([rows.loc[key, "median_pair_delta"] for key, _ in VARIANTS])
            lo = np.array([rows.loc[key, "q25_pair_delta"] for key, _ in VARIANTS])
            hi = np.array([rows.loc[key, "q75_pair_delta"] for key, _ in VARIANTS])
            ax.errorbar(x + offset, med, yerr=np.vstack([med - lo, hi - med]),
                        color=COLORS[view], linestyle=linestyle, marker="o",
                        markersize=3.8, linewidth=1.0, capsize=2, label=label)
        ax.axhline(0, color="#555555", linewidth=0.9)
        ax.axhspan(-0.55, 0, color="#EAF3EA", alpha=0.65, zorder=0)
        ax.set_xticks(x, [label for _, label in VARIANTS], rotation=25, ha="right")
        ax.set_title("Full xyz" if view == "xyz" else "Fast xy (supporting)")
        ax.set_xlabel("Marker definition")
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)
        ax.legend(frameon=False, loc="lower right")
    axes[0].set_ylabel(r"Median pair $p_{KL}-p_{firing}$ (IQR)")
    axes[0].set_ylim(-0.55, 0.22)
    save(fig, "FigS2_multiseed_marker_robustness")


def main():
    style()
    marker_robustness()


if __name__ == "__main__":
    main()
