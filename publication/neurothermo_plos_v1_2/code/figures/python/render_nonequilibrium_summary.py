#!/usr/bin/env python3
"""Render the focused slow-coordinate and time-reversal summary (main Fig 5)."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo-mpl")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "nonequilibrium_geometry_v1_0_1"
OUT = ROOT / "results" / "figures"

BLUE = "#2B6CB0"
RED = "#C44E52"
PURPLE = "#7A5195"
GREY = "#555555"


def summarize(frame, value):
    grouped = frame.groupby("p", sort=True)[value]
    return pd.DataFrame({
        "p": grouped.median().index,
        "median": grouped.median().values,
        "q25": grouped.quantile(0.25).values,
        "q75": grouped.quantile(0.75).values,
    })


def ribbon(ax, summary, color, label=None, alpha=0.18):
    x = summary["p"].to_numpy(float)
    y = summary["median"].to_numpy(float)
    lo = summary["q25"].to_numpy(float)
    hi = summary["q75"].to_numpy(float)
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, linewidth=0)
    ax.plot(x, y, color=color, linewidth=2.2, label=label)


def panel_label(ax, label):
    ax.text(0.015, 0.985, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold")


def main():
    balanced = pd.read_csv(DATA / "animal_balanced_geometry.csv")
    strata = pd.read_csv(DATA / "animal_pair_balanced_geometry.csv.gz")
    strata["slow_fraction"] = (
        strata["path_fi_conditional_z_given_xy"] / strata["path_fi_xyz"]
    )

    fi_xyz = balanced[["p", "path_fi_xyz_median", "path_fi_xyz_q25", "path_fi_xyz_q75"]].rename(
        columns={"path_fi_xyz_median": "median", "path_fi_xyz_q25": "q25", "path_fi_xyz_q75": "q75"}
    )
    fi_xy = balanced[["p", "path_fi_xy_median", "path_fi_xy_q25", "path_fi_xy_q75"]].rename(
        columns={"path_fi_xy_median": "median", "path_fi_xy_q25": "q25", "path_fi_xy_q75": "q75"}
    )
    slow = summarize(strata, "slow_fraction")
    db = balanced[["p", "markov_db_violation_median", "markov_db_violation_q25", "markov_db_violation_q75"]].rename(
        columns={"markov_db_violation_median": "median", "markov_db_violation_q25": "q25", "markov_db_violation_q75": "q75"}
    )

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    })

    fig = plt.figure(figsize=(10.5, 7.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9], hspace=0.34, wspace=0.28)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    ribbon(ax_a, fi_xyz, BLUE, r"Full $xyz$")
    ribbon(ax_a, fi_xy, RED, r"Fast $xy$")
    ax_a.set_yscale("log")
    ax_a.set_ylabel("Path Fisher information")
    ax_a.legend(loc="upper center")
    panel_label(ax_a, "A")

    ribbon(ax_b, slow, PURPLE)
    ax_b.axhline(0.5, color=GREY, linestyle="--", linewidth=1.0)
    ax_b.set_ylim(0.0, 1.0)
    ax_b.set_ylabel("Conditional slow-coordinate contribution")
    panel_label(ax_b, "B")

    ribbon(ax_c, db, BLUE)
    ax_c.axhline(0.05, color=GREY, linestyle="--", linewidth=1.0)
    ax_c.text(0.02, 0.072, "Detailed-balance threshold", color=GREY,
              ha="left", va="bottom", fontsize=8)
    ax_c.set_ylim(0.0, 0.86)
    ax_c.set_ylabel("Detailed-balance violation")
    panel_label(ax_c, "C")

    for ax in (ax_a, ax_b, ax_c):
        ax.set_xlim(-0.01, 1.01)
        ax.set_xlabel("Constructed path coordinate, p")
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "Fig5_nonequilibrium_geometry.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Fig5_nonequilibrium_geometry.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
