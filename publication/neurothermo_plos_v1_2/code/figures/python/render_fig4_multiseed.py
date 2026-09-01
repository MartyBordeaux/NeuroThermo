#!/usr/bin/env python3
"""Render Fig 4 from the full-coverage KL convergence analysis v1.0.1."""

from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo-mpl")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "kl_convergence_v1_0_1"
FIGURES = ROOT / "results" / "figures"
PRIMARY_VARIANT = "seed_median_curve_isotonic"
PRIMARY_DT = 0.025


def panel_label(axis, label):
    axis.text(-0.13, 1.075, label, transform=axis.transAxes,
              fontsize=12, fontweight="bold", va="top")


def main():
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    pairs = pd.read_csv(DATA / "pair_markers_both_orders.csv")
    cells = pd.read_csv(DATA / "cell_balanced_deltas.csv")
    ensemble = pd.read_csv(DATA / "ensemble_convergence_summary.csv")

    primary_pairs = pairs[
        (pairs.dt_ms == PRIMARY_DT)
        & (pairs.marker_variant == PRIMARY_VARIANT)
    ].copy()
    primary_cells = cells[
        (cells.dt_ms == PRIMARY_DT)
        & (cells.marker_variant == PRIMARY_VARIANT)
        & (cells.aggregation_order == "marker_first")
    ].copy()

    xyz = "#54278F"
    xy = "#E08214"
    grey = "#555555"
    pale_green = "#EAF3EA"
    rng = np.random.default_rng(20260829)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.35), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.flat

    # A: paired full-state markers versus the firing marker.
    a = primary_pairs[
        (primary_pairs.view == "xyz")
        & (primary_pairs.aggregation_order == "marker_first")
    ]
    lim = (0.05, 0.95)
    ax_a.plot(lim, lim, color=grey, linestyle="--", linewidth=0.9)
    negative = a.kl_minus_firing_p < 0
    ax_a.scatter(a.firing_balance_p[negative], a.kl_balance_p[negative],
                 color=xyz, s=24, alpha=0.78, label="KL at lower p")
    ax_a.scatter(a.firing_balance_p[~negative], a.kl_balance_p[~negative],
                 facecolor="white", edgecolor=xyz, s=25, linewidth=0.9,
                 label="KL not at lower p")
    ax_a.text(0.04, 0.96, f"{int(negative.sum())}/{len(a)} pairs below identity",
              transform=ax_a.transAxes, color=xyz, fontsize=7.2,
              ha="left", va="top")
    ax_a.set(xlabel="Firing-phenotype balance, p",
             ylabel="Full-state KL balance, p", xlim=lim, ylim=lim)
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.legend(loc="lower right", fontsize=6.8, frameon=True)
    panel_label(ax_a, "A")

    # B: primary pair-level deltas under both noncommuting aggregation orders.
    groups = [
        ("xyz", "marker_first", "Full xyz\nmarker first", xyz),
        ("xyz", "curve_first", "Full xyz\ncurve first", xyz),
        ("xy", "marker_first", "Fast xy\nmarker first", xy),
        ("xy", "curve_first", "Fast xy\ncurve first", xy),
    ]
    arrays = []
    for view, order, _, _ in groups:
        arrays.append(primary_pairs[
            (primary_pairs.view == view)
            & (primary_pairs.aggregation_order == order)
        ].kl_minus_firing_p.to_numpy())
    box = ax_b.boxplot(arrays, positions=np.arange(4), widths=0.55,
                       patch_artist=True, showfliers=False,
                       medianprops={"color": "black", "linewidth": 1.0})
    for patch, (_, _, _, color) in zip(box["boxes"], groups):
        patch.set_facecolor(color)
        patch.set_alpha(0.13)
        patch.set_edgecolor(color)
    for i, (values, (_, _, _, color)) in enumerate(zip(arrays, groups)):
        jitter = rng.uniform(-0.11, 0.11, size=len(values))
        ax_b.scatter(i + jitter, values, s=11, color=color, alpha=0.52,
                     linewidth=0)
    ax_b.axhline(0, color=grey, linewidth=0.9)
    ax_b.axvline(1.5, color="#BBBBBB", linewidth=0.7)
    ax_b.axhspan(-0.48, 0, color=pale_green, alpha=0.58, zorder=0)
    ax_b.set_xticks(np.arange(4), [g[2] for g in groups])
    ax_b.set_ylabel(r"$p_{KL}-p_{firing}$")
    ax_b.set_ylim(-0.48, 0.34)
    ax_b.text(0.02, 0.03, "KL balance at lower constructed p",
              transform=ax_b.transAxes, color="#3E6B3E", fontsize=7.0)
    panel_label(ax_b, "B")

    # C: endpoint-cell-balanced primary marker-first differences.
    order = list(primary_cells[primary_cells.view == "xyz"]
                 .sort_values(["endpoint_genotype", "endpoint_cell"])
                 .endpoint_cell)
    y_positions = np.arange(len(order))[::-1]
    ax_c.axvspan(-0.46, 0, color=pale_green, alpha=0.70, zorder=0)
    for y, cell in zip(y_positions, order):
        subset = primary_cells[primary_cells.endpoint_cell == cell].set_index("view")
        x_xyz = float(subset.loc["xyz", "kl_minus_firing_p"])
        x_xy = float(subset.loc["xy", "kl_minus_firing_p"])
        ax_c.plot([x_xyz, x_xy], [y, y], color="#BDBDBD", linewidth=0.8)
        marker = "o" if subset.iloc[0].endpoint_genotype == "WT" else "^"
        ax_c.scatter(x_xyz, y + 0.08, color=xyz, marker=marker, s=24,
                     label="Full xyz" if cell == order[0] else None, zorder=2)
        ax_c.scatter(x_xy, y - 0.08, color=xy, marker=marker, s=24,
                     label="Fast xy" if cell == order[0] else None, zorder=2)
    ax_c.axvline(0, color=grey, linestyle="--", linewidth=0.9)
    ax_c.set_yticks(y_positions, [cell.replace("_", " ") for cell in order])
    ax_c.set(xlabel=r"Cell-balanced $p_{KL}-p_{firing}$",
             ylabel="Endpoint cell", xlim=(-0.46, 0.09))
    ax_c.legend(loc="lower right", frameon=True)
    panel_label(ax_c, "C")

    # D: integration-step stability for both aggregation orders.
    central = ensemble[
        (ensemble.marker_variant == PRIMARY_VARIANT)
        & ensemble.view.isin(["xyz", "xy"])
    ].copy()
    dt_order = [0.05, 0.025, 0.0125]
    x = np.arange(3)
    for view, color in [("xyz", xyz), ("xy", xy)]:
        for aggregation, linestyle, label_suffix in [
            ("marker_first", "-", "marker first"),
            ("curve_first", "--", "curve first"),
        ]:
            z = central[(central.view == view)
                        & (central.aggregation_order == aggregation)].set_index("dt_ms")
            y = np.array([z.loc[v, "median_pair_delta"] for v in dt_order])
            ax_d.plot(x, y, color=color, linestyle=linestyle, marker="o",
                      markersize=4, linewidth=1.25,
                      label=f"{view}, {label_suffix}")
    ax_d.axhline(0, color=grey, linewidth=0.9)
    ax_d.axhspan(-0.20, 0, color=pale_green, alpha=0.58, zorder=0)
    ax_d.set_xticks(x, ["0.05", "0.025", "0.0125"])
    ax_d.set(xlabel="Integration step, dt (ms)",
             ylabel=r"Median pair $p_{KL}-p_{firing}$", ylim=(-0.20, 0.02))
    ax_d.legend(fontsize=6.7, loc="upper right", frameon=True)
    panel_label(ax_d, "D")

    for axis in axes.flat:
        axis.tick_params(direction="out", length=3)

    pdf_path = FIGURES / "Fig4_thermodynamic_information_geometry.pdf"
    png_path = FIGURES / "Fig4_thermodynamic_information_geometry.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
