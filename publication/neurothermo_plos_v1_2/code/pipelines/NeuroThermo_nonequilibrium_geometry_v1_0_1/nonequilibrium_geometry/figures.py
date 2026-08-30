from __future__ import annotations

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo-matplotlib")

import matplotlib.pyplot as plt
import numpy as np


BLUE, RED, GREY = "#2f73b7", "#c94f52", "#555555"


def _band(axis, frame, prefix, color, label):
    p = frame["p"].to_numpy()
    median = frame[prefix + "_median"].to_numpy()
    q25 = frame[prefix + "_q25"].to_numpy()
    q75 = frame[prefix + "_q75"].to_numpy()
    axis.fill_between(p, q25, q75, color=color, alpha=0.18, linewidth=0)
    axis.plot(p, median, color=color, lw=2.2, label=label)


def make_figures(output, ensemble, protocol, fluctuation):
    output = Path(output)
    figure_dir = output / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), constrained_layout=True)
    _band(axes[0, 0], ensemble, "path_fi_xyz", BLUE, "FI, xyz")
    _band(axes[0, 0], ensemble, "path_fi_xy", RED, "FI, xy")
    axes[0, 0].set(ylabel="Path Fisher information")
    axes[0, 0].set_yscale("log")
    axes[0, 0].legend(frameon=False)
    _band(axes[0, 1], ensemble, "friction_metric", GREY, "Friction metric")
    axes[0, 1].set(ylabel="Integrated force covariance")
    axes[0, 1].set_yscale("log")
    _band(axes[1, 0], ensemble, "markov_db_violation", BLUE, "Coarse Markov")
    axes[1, 0].axhline(0.05, color="0.4", ls="--", lw=1)
    axes[1, 0].set(ylabel="Detailed-balance violation")
    _band(axes[1, 1], ensemble, "stationary_current_divergence_relative", GREY, "Divergence residual")
    axes[1, 1].axhline(0.25, color="0.4", ls="--", lw=1)
    axes[1, 1].set(ylabel="Relative residual")
    for label, axis in zip(("A", "B", "C", "D"), axes.flat):
        axis.text(0.01, 0.98, label, transform=axis.transAxes, va="top", ha="left", fontweight="bold")
        axis.set_xlabel("Constructed path coordinate, p")
        axis.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"Fig_nonequilibrium_geometry.{suffix}", dpi=300)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(5.8, 4.5), constrained_layout=True)
    styles = {"path_fi_xyz": (BLUE, "FI xyz"), "path_fi_xy": (RED, "FI xy"), "friction": (GREY, "Friction")}
    for metric, group in protocol.groupby("metric"):
        color, label = styles.get(metric, (GREY, metric))
        summary = group.groupby("normalized_time")["p"].agg(
            median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75)
        ).reset_index()
        axis.fill_between(summary["normalized_time"], summary["q25"], summary["q75"], color=color, alpha=0.15, linewidth=0)
        axis.plot(summary["normalized_time"], summary["median"], lw=2.2, color=color, label=label)
    axis.plot([0, 1], [0, 1], ls="--", color="0.7", label="Linear p")
    axis.set(xlabel="Normalized protocol time", ylabel="p")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"Fig_adaptive_protocols.{suffix}", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.1), constrained_layout=True)
    observed = set(fluctuation["schedule"])
    labels = [name for name in ("linear", "path_fi_xyz", "path_fi_xy", "friction") if name in observed]
    locations = np.arange(len(labels))
    colors = {"WT_to_SCA3": BLUE, "SCA3_to_WT": RED}
    width = 0.36
    for direction_index, direction in enumerate(("WT_to_SCA3", "SCA3_to_WT")):
        group = fluctuation.loc[fluctuation["direction"].eq(direction)]
        hs = group.groupby("schedule")["mean_exp_minus_Y"].median().reindex(labels)
        ift = group.groupby("schedule")["median_sigma"].median().reindex(labels)
        offset = (direction_index - 0.5) * width
        axes[0].bar(locations + offset, hs, width=width, color=colors[direction], label=direction)
        axes[1].bar(locations + offset, ift, width=width, color=colors[direction], label=direction)
    axes[0].axhline(1.0, color="0.25", ls="--", lw=1)
    for panel_label, axis in zip(("A", "B"), axes):
        axis.set_xticks(locations, labels, rotation=15)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(0.01, 0.98, panel_label, transform=axis.transAxes, va="top", ha="left", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("Sampled exponential average")
    axes[1].set_ylabel("Median log path ratio")
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"Fig_fluctuation_diagnostics.{suffix}", dpi=300)
    plt.close(fig)
