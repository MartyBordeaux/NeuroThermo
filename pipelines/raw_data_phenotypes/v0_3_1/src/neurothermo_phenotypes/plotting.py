from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Union

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo_matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LABELS = {
    "sustained_rate_hz": "Sustained rate, Hz",
    "first_spike_latency_ms": "First-spike latency, ms",
    "work_per_spike_fJ": "External work/spike, fJ",
    "mean_power_signed_fW": "External mean power, fW",
    "permutation_entropy_norm": "Permutation entropy",
    "spectral_entropy_norm": "Spectral entropy",
    "predictive_information_nats": "Predictive information, nats",
    "path_kl_rate_excess_nats_s": "Surrogate-excess path KL, nats/s",
}


def plot_response_curves(
    summary: pd.DataFrame, output: Union[str, Path], x_label: str,
    shared_support_only: bool = False,
    inference: pd.DataFrame = None,
) -> None:
    if summary.empty:
        return
    features = [f for f in LABELS if f in set(summary.feature)]
    ncols, nrows = 2, math.ceil(len(features) / 2)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(11, 3.5 * nrows), squeeze=False, constrained_layout=True
    )
    colors = {"WT": "#2563a8", "SCA3": "#b43c39"}
    q_by_feature = {}
    if inference is not None and not inference.empty and "fdr_q_global" in inference:
        q_by_feature = dict(zip(inference["feature"], inference["fdr_q_global"]))
    for ax, feature in zip(axes.flat, features):
        panel = summary[summary.feature == feature]
        if shared_support_only and "shared_support" in panel:
            panel = panel[panel["shared_support"]]
        if panel.empty:
            ax.text(.5, .5, "No shared support\n(n >= 3 cells/group)", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(LABELS[feature])
            ax.set_xlabel(x_label)
            continue
        for group in ["WT", "SCA3"]:
            data = panel[panel.group == group].sort_values("axis_value")
            if data.empty:
                continue
            x = data.axis_value.to_numpy(float)
            ax.plot(x, data["median"], color=colors[group], lw=2, label=group)
            ax.fill_between(x, data.q25.to_numpy(float), data.q75.to_numpy(float), color=colors[group], alpha=.18)
        title = LABELS[feature]
        q_value = q_by_feature.get(feature)
        if q_value is not None and np.isfinite(q_value):
            title += f" (global q={q_value:.3g})"
        ax.set_title(title)
        unique_x = np.sort(panel["axis_value"].unique())
        count_parts = []
        for group in ["WT", "SCA3"]:
            values = panel.loc[panel.group == group, "n_cells"].to_numpy(int)
            if len(values):
                count_parts.append(f"{group} n={values.min()}–{values.max()}")
        if len(unique_x) <= 13:
            tick_labels = []
            for value in unique_x:
                level = panel[panel.axis_value == value]
                counts = {
                    row.group: int(row.n_cells) for row in level.itertuples(index=False)
                }
                tick_labels.append(f"{value:g}\n{counts.get('WT', 0)}/{counts.get('SCA3', 0)}")
            ax.set_xticks(unique_x, tick_labels)
            ax.set_xlabel(x_label + "\nsecond line: n WT/SCA3")
        else:
            ax.set_xlabel(x_label)
        if count_parts and len(unique_x) > 13:
            ax.text(.99, .02, "; ".join(count_parts), transform=ax.transAxes, fontsize=7, va="bottom", ha="right")
        ax.grid(alpha=.2)
    for ax in axes.flat[len(features):]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_effect_heatmap(comparisons: pd.DataFrame, output: Union[str, Path]) -> None:
    if comparisons.empty:
        return
    use = comparisons[comparisons.feature.isin(LABELS)].copy()
    if use.empty:
        return
    table = use.pivot(index="feature", columns="current_pA", values="cliffs_delta_SCA3_vs_WT")
    q_table = use.pivot(index="feature", columns="current_pA", values="fdr_q_global")
    table = table.reindex([f for f in LABELS if f in table.index])
    q_table = q_table.reindex(index=table.index, columns=table.columns)
    fig, ax = plt.subplots(figsize=(12, max(3, .55 * len(table))))
    image = ax.imshow(table.to_numpy(float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_yticks(range(len(table)), [LABELS.get(x, x) for x in table.index])
    ax.set_xticks(range(len(table.columns)), [f"{x:g}" for x in table.columns], rotation=45, ha="right")
    ax.set_xlabel("Injected current, pA")
    ax.set_title("Cliff's delta: SCA3 relative to WT")
    for row in range(len(table.index)):
        for col in range(len(table.columns)):
            if np.isfinite(q_table.iloc[row, col]) and q_table.iloc[row, col] < 0.05:
                ax.text(col, row, "*", ha="center", va="center", color="black", fontsize=11)
    fig.colorbar(image, ax=ax, label="Effect size")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_qc(features: pd.DataFrame, output: Union[str, Path]) -> None:
    counts = features.groupby(["group", "qc_status"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=["PASS", "WARN", "FAIL"], fill_value=0)
    status_rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
    cell_status = (
        features.assign(_rank=features["qc_status"].map(status_rank))
        .groupby(["group", "cell_id"])["_rank"].max().map({0: "PASS", 1: "WARN", 2: "FAIL"})
        .rename("qc_status").reset_index()
    )
    cell_counts = cell_status.groupby(["group", "qc_status"]).size().unstack(fill_value=0)
    cell_counts = cell_counts.reindex(index=counts.index, columns=["PASS", "WARN", "FAIL"], fill_value=0)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    colors = {"PASS": "#3b7f55", "WARN": "#d69e2e", "FAIL": "#a6a6a6"}
    labels = {"PASS": "Pass", "WARN": "Warning", "FAIL": "Fatal QC"}
    for ax, table, title, ylabel in [
        (axes[0], counts, "Sweep-level QC", "Sweeps"),
        (axes[1], cell_counts, "Cell-level worst QC", "Cells"),
    ]:
        x = np.arange(len(table))
        bottom = np.zeros(len(table))
        for status in ["PASS", "WARN", "FAIL"]:
            ax.bar(x, table[status], bottom=bottom, label=labels[status], color=colors[status])
            bottom += table[status].to_numpy(float)
        ax.set_xticks(x, table.index)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(.5, 1.06))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
