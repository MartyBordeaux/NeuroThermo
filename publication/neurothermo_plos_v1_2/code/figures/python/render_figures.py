#!/usr/bin/env python3
"""Render the frozen v1.1.8 NeuroThermo publication figures with Matplotlib.

The panel definitions, data columns, colors, scales, thresholds, and layouts
mirror the frozen R scripts.  This renderer exists so the integrated review
PDF can be built in environments without R; the source CSV tables remain the
single source of numerical values.
"""

from pathlib import Path
import os
import subprocess
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/neurothermo-mpl")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


RELEASE_ROOT = Path(__file__).resolve().parents[3]
DATA = RELEASE_ROOT / "data" / "figure_source"
OUT = RELEASE_ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

WT = "#2B6CB0"
SCA3 = "#C44E52"
EARLY = "#3C8D73"
COUPLED = "#4C5F9E"
LATE = "#B36A3C"
KAPPA = "#7A5195"
JCOL = "#EF5675"
COMBINED = "#003F5C"
INTERACTION = "#A0A0A0"
GROUP_COLORS = {"WT": WT, "SCA3": SCA3}
PATH_COLORS = {"Drive early": EARLY, "Coupled": COUPLED, "Drive late": LATE}
STAGE_COLORS = {"WT exit": "#4D4D4D", "Balance": "#8C6D31", "SCA3 entry": "#8B1A1A"}
STAGE_STYLES = {"WT exit": ":", "Balance": "--", "SCA3 entry": "-."}


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.framealpha": 0.9,
    "legend.edgecolor": "#CCCCCC",
    "savefig.facecolor": "white",
})


def read(name):
    return pd.read_csv(DATA / name)


def panel(ax, label):
    ax.text(0.005, 0.995, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold")


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


def stage_rows():
    s = read("fig2_primary_isi_staging.csv")
    s = s[(s.path_family == "coupled") & (s.subset == "core_secure_pairs")].copy()
    mapping = {"wt_exit_p_isi": "WT exit", "balance_p_isi": "Balance", "sca3_entry_p_isi": "SCA3 entry"}
    s["stage"] = s.metric.map(mapping)
    return s.dropna(subset=["stage"])


def add_stage_lines(ax, show_labels=False):
    for _, row in stage_rows().iterrows():
        ax.axvline(row["median"], color="#666666", lw=0.8,
                   ls=STAGE_STYLES[row.stage], alpha=0.6,
                   label=row.stage if show_labels else None)


def fig1():
    d = read("fig1_endpoint_cells.csv")
    fig = plt.figure(figsize=(7.2, 6.3), constrained_layout=True)
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.02])
    top = outer[0].subgridspec(1, 3, wspace=0.38)
    rng = np.random.default_rng(7123)
    specs = [
        ("rheobase_J_best", r"$J_{rheo}$ (pA/pF)", True, "A"),
        ("exp_q75_firing_rate_hz", "Experimental firing rate at q = 0.75 (Hz)", False, "B"),
        ("exp_q75_mean_isi_ms", "Experimental mean ISI at q = 0.75 (ms)", False, "C"),
    ]
    for i, (var, ylabel, logy, tag) in enumerate(specs):
        ax = fig.add_subplot(top[0, i])
        for xpos, group in enumerate(["WT", "SCA3"]):
            vals = d.loc[d.group == group, var].dropna().to_numpy(float)
            x = xpos + rng.uniform(-0.075, 0.075, len(vals))
            ax.scatter(x, vals, s=24, color=GROUP_COLORS[group], alpha=0.85,
                       edgecolor="none", zorder=3)
            q25, med, q75 = np.quantile(vals, [0.25, 0.5, 0.75])
            ax.errorbar(xpos, med, yerr=[[med-q25], [q75-med]], color="black",
                        capsize=5, lw=1.2, zorder=4)
            ax.plot([xpos-0.08, xpos+0.08], [med, med], color="black", lw=2.0, zorder=5)
        ax.set_xticks([0, 1], ["WT", "SCA3"])
        ax.set_ylabel(ylabel)
        if logy:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        panel(ax, tag)

    bottom = outer[1].subgridspec(1, 2, wspace=0.28)
    for i, (metric, ex, mo) in enumerate([
        ("Firing rate", "exp_q75_firing_rate_hz", "model_q75_firing_rate_hz"),
        ("Mean ISI", "exp_q75_mean_isi_ms", "model_q75_mean_isi_ms"),
    ]):
        ax = fig.add_subplot(bottom[0, i])
        for group in ["WT", "SCA3"]:
            g = d[d.group == group]
            ax.scatter(g[ex], g[mo], s=25, color=GROUP_COLORS[group], alpha=0.85,
                       edgecolor="none", label=group)
        vals = np.r_[d[ex].dropna().to_numpy(), d[mo].dropna().to_numpy()]
        lo, hi = max(0, vals.min()*0.9), vals.max()*1.05
        ax.plot([lo, hi], [lo, hi], ls="--", lw=0.9, color="#555555")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("Experiment")
        ax.set_ylabel("Best-fit HR model" if i == 0 else "")
        ax.set_title(metric, fontweight="bold", pad=5)
        if i == 0:
            panel(ax, "D")
            ax.legend(loc="lower left", bbox_to_anchor=(0.32, 0.06))
    save(fig, "Fig1_endpoint_phenotype")


def fig2():
    curves = read("fig2_core_secure_curves.csv")
    stages = read("fig2_primary_isi_staging.csv")
    stages = stages[stages.subset == "core_secure_pairs"].copy()
    path_map = {"drive_early": "Drive early", "coupled": "Coupled", "drive_late": "Drive late"}
    metric_map = {"wt_exit_p_isi": "WT exit", "balance_p_isi": "Balance", "sca3_entry_p_isi": "SCA3 entry"}
    curves["path"] = curves.path_family.map(path_map)
    stages["path"] = stages.path_family.map(path_map)
    stages["stage"] = stages.metric.map(metric_map)
    stages = stages.dropna(subset=["stage"])
    ref = read("fig2_projection_reference.csv")
    ref = ref[ref.projection == "isi_primary_v1_0_frozen"].iloc[0]
    thresholds = {"WT exit": ref.wt_exit_A_threshold, "Balance": 0.5, "SCA3 entry": ref.sca3_entry_A_threshold}

    fig = plt.figure(figsize=(7.2, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], hspace=0.18)
    ax = fig.add_subplot(gs[0, 0])
    p = np.linspace(0, 1, 201)
    prot = {"Drive early": 2*p-p*p, "Coupled": p, "Drive late": p*p}
    ax.plot(p, p, color="#BBBBBB", lw=0.7)
    for name, y in prot.items(): ax.plot(p, y, color=PATH_COLORS[name], lw=1.6, label=name)
    ax.scatter([0, 1], [0, 1], s=38, color=[WT, SCA3], edgecolor="white", linewidth=0.8, zorder=6)
    ax.text(0.18, 0.07, "WT", color=WT, fontsize=8, fontweight="bold",
            ha="center", va="bottom")
    ax.text(0.82, 0.93, "SCA3", color=SCA3, fontsize=8, fontweight="bold",
            ha="center", va="top")
    ax.set_aspect("equal", adjustable="box"); ax.set(xlabel=r"$p_{intrinsic}$", ylabel=r"$p_{drive}$", xlim=(0,1), ylim=(0,1))
    panel(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    for path in ["Drive early", "Coupled", "Drive late"]:
        g = curves[curves.path == path].sort_values("path_progress")
        x = g.path_progress.to_numpy(float); med = g.A_isi_median.to_numpy(float)
        lo = g.A_isi_q25.to_numpy(float); hi = g.A_isi_q75.to_numpy(float)
        ax.fill_between(x, lo, hi, color=PATH_COLORS[path], alpha=0.12, linewidth=0)
        ax.plot(x, med, color=PATH_COLORS[path], lw=1.5, label=path)
    for stage, y in thresholds.items(): ax.axhline(y, color="#4D4D4D", lw=0.8, ls=STAGE_STYLES[stage])
    wt_y = curves.loc[np.isclose(curves.path_progress, 0), "A_isi_median"].median()
    sca3_y = curves.loc[np.isclose(curves.path_progress, 1), "A_isi_median"].median()
    ax.scatter([0, 1], [wt_y, sca3_y], s=38, color=[WT, SCA3], edgecolor="white", linewidth=0.8, zorder=6)
    ax.annotate("WT", (0, wt_y), xytext=(7, 8), textcoords="offset points", color=WT,
                fontsize=8, fontweight="bold", ha="left", va="bottom")
    ax.text(0.80, 1.00, "SCA3", color=SCA3, fontsize=8, fontweight="bold",
            ha="center", va="top")
    ax.set(xlabel="Path progress, p", ylabel=r"$A_{ISI}$", xlim=(0,1))
    panel(ax, "B")

    ax = fig.add_subplot(gs[1, :])
    paths = ["Drive early", "Coupled", "Drive late"]
    stages_order = ["WT exit", "Balance", "SCA3 entry"]
    offsets = [-0.18, 0, 0.18]
    for j, stage in enumerate(stages_order):
        for i, path in enumerate(paths):
            row = stages[(stages.path == path) & (stages.stage == stage)].iloc[0]
            ax.errorbar(i+offsets[j], row["median"], yerr=[[row["median"]-row.q25],[row.q75-row["median"]]],
                        fmt="o", ms=5, capsize=3, lw=1.1, color=STAGE_COLORS[stage],
                        label=stage if i == 0 else None)
    ax.set_xticks(range(3), paths); ax.set_ylim(0,1); ax.set_ylabel("Stage location, p")
    ax.legend(ncol=3, loc="upper center")
    panel(ax, "C")
    save(fig, "Fig2_transition_staging")


def grid_matrix(df, column):
    tab = df.pivot(index="p_drive", columns="p_intrinsic", values=column).sort_index().sort_index(axis=1)
    return tab.columns.to_numpy(float), tab.index.to_numpy(float), tab.to_numpy(float)


def fig3():
    surf = read("fig3_drive_surface_core_secure.csv")
    comp = read("fig3_coupled_component_sensitivity_core_secure.csv")
    inter = read("fig3_interaction_at_boundaries_core_secure.csv")
    ref = read("fig2_projection_reference.csv")
    ref = ref[ref.projection == "isi_primary_v1_0_frozen"].iloc[0]
    thresholds = [ref.wt_exit_A_threshold, 0.5, ref.sca3_entry_A_threshold]
    srows = stage_rows()
    stage_pos = dict(zip(srows.stage, srows["median"]))
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.6), constrained_layout=True)

    ax = axs[0,0]
    x,y,z = grid_matrix(surf, "A_isi_median")
    lim = np.nanmax(np.abs(z-0.5))
    im = ax.pcolormesh(x,y,z, shading="auto", cmap="coolwarm", vmin=0.5-lim, vmax=0.5+lim)
    ax.contour(x,y,z, levels=thresholds, colors="white", linewidths=0.8)
    ax.scatter([0, 1], [0, 1], s=38, color=[WT, SCA3], edgecolor="white", linewidth=0.8, zorder=6)
    wt_label = ax.annotate("WT", (0, 0), xytext=(7, 8), textcoords="offset points", color=WT,
                           fontsize=8, fontweight="bold", ha="left", va="bottom")
    sca3_label = ax.annotate("SCA3", (1, 1), xytext=(-7, -8), textcoords="offset points", color=SCA3,
                             fontsize=8, fontweight="bold", ha="right", va="top")
    for endpoint_label in (wt_label, sca3_label):
        endpoint_label.set_path_effects([
            path_effects.Stroke(linewidth=2.5, foreground="white"),
            path_effects.Normal(),
        ])
    fig.colorbar(im, ax=ax, label=r"$A_{ISI}$", fraction=0.046, pad=0.04)
    ax.set(xlabel=r"$p_{intrinsic}$", ylabel=r"$p_{drive}$"); ax.set_aspect("equal")
    panel(ax,"A")

    ax = axs[0,1]
    x,y,z = grid_matrix(surf, "drive_dominance_fraction_median")
    im = ax.pcolormesh(x,y,z, shading="auto", cmap="coolwarm", vmin=0, vmax=1)
    ax.scatter([0, 1], [0, 1], s=38, color=[WT, SCA3], edgecolor="white", linewidth=0.8, zorder=6)
    ax.annotate("WT", (0, 0), xytext=(7, 8), textcoords="offset points", color=WT,
                fontsize=8, fontweight="bold", ha="left", va="bottom")
    ax.annotate("SCA3", (1, 1), xytext=(-7, -8), textcoords="offset points", color=SCA3,
                fontsize=8, fontweight="bold", ha="right", va="top")
    fig.colorbar(im, ax=ax, label="Drive dominance", fraction=0.046, pad=0.04)
    ax.set(xlabel=r"$p_{intrinsic}$", ylabel=r"$p_{drive}$"); ax.set_aspect("equal")
    panel(ax,"B")

    ax = axs[1,0]
    combined = surf[np.isclose(surf.p_intrinsic, surf.p_drive)].sort_values("p_intrinsic")
    ax.axvspan(0.0, 0.30, color=JCOL, alpha=0.035, linewidth=0)
    ax.axvspan(0.70, 1.0, color=KAPPA, alpha=0.035, linewidth=0)
    components = [
        ("Combined drive", combined, "dA_ddrive", COMBINED),
        (r"$\kappa_I$", comp, "dA_dkappa", KAPPA),
        ("Applied J", comp, "dA_dJ", JCOL),
    ]
    for label, source, prefix, color in components:
        x = source.p_intrinsic.to_numpy(float); med=source[f"{prefix}_median"].to_numpy(float)
        lo=source[f"{prefix}_q25"].to_numpy(float); hi=source[f"{prefix}_q75"].to_numpy(float)
        ok=np.isfinite(lo)&np.isfinite(hi)
        ax.fill_between(x[ok],lo[ok],hi[ok],color=color,alpha=.10,linewidth=0)
        ax.plot(x,med,color=color,lw=1.35,label=label)
    ax.axhline(0,color="#666666",lw=.7)
    for v in stage_pos.values(): ax.axvline(v,color="#777777",lw=.65,ls="--")
    ax.text(0.15, 0.985, "J-shaped", transform=ax.get_xaxis_transform(),
            color=JCOL, fontsize=6.8, ha="center", va="top")
    ax.text(0.85, 0.985, r"$\kappa_I$-shaped", transform=ax.get_xaxis_transform(),
            color=KAPPA, fontsize=6.8, ha="center", va="top")
    ax.set(xlabel="Coupled path progress, p", ylabel=r"Local sensitivity $dA_{ISI}/dp$", xlim=(0,1))
    ax.legend(loc="best"); panel(ax,"C")

    ax=axs[1,1]
    stage_key={"WT_exit":"WT exit","balance":"Balance","SCA3_entry":"SCA3 entry"}
    components=[("Combined","median_abs_delta_combined",COMBINED),(r"$\kappa_I$","median_abs_delta_kappa",KAPPA),
                ("Applied J","median_abs_delta_J",JCOL),("Interaction","median_abs_interaction",INTERACTION)]
    positions={"WT exit":0,"Balance":1,"SCA3 entry":2}; offsets=np.linspace(-.24,.24,4)
    for k,(label,prefix,color) in enumerate(components):
        for _,row in inter.iterrows():
            stage=stage_key[row.stage]; xpos=positions[stage]+offsets[k]
            med=row[f"{prefix}_median"]; lo=row[f"{prefix}_q25"]; hi=row[f"{prefix}_q75"]
            ax.errorbar(xpos,med,yerr=[[med-lo],[hi-med]],fmt="o",ms=4,capsize=2,lw=1,color=color,
                        label=label if stage=="WT exit" else None)
    ax.set_xticks([0,1,2],["WT exit","Balance","SCA3 entry"])
    ax.set_ylabel(r"Median $|\Delta A_{ISI}|$ near boundary")
    ax.legend(ncol=2,loc="upper right"); panel(ax,"D")
    save(fig,"Fig3_intrinsic_drive_decomposition")


def fig4_legacy():
    d=read("fig4_thermodynamic_curves.csv")
    fig=plt.figure(figsize=(7.2,6.5),constrained_layout=True)
    outer=fig.add_gridspec(2,1,height_ratios=[1,1],hspace=.18)
    top=outer[0].subgridspec(1,2,wspace=.28)
    bottom=outer[1].subgridspec(1,3,wspace=.35)
    x=d.p.to_numpy(float)

    ax=fig.add_subplot(top[0,0])
    for endpoint,prefix,color in [("WT","kl_to_wt",WT),("SCA3","kl_to_sca3",SCA3)]:
        med=d[f"{prefix}_median"].to_numpy(float); lo=d[f"{prefix}_q25"].to_numpy(float); hi=d[f"{prefix}_q75"].to_numpy(float)
        ax.fill_between(x,lo,hi,color=color,alpha=.12,linewidth=0)
        ax.plot(x,med,color=color,lw=1.4,label=endpoint)
    add_stage_lines(ax,True); ax.set(xlabel="p",ylabel=r"$D_{KL}$",xlim=(0,1)); panel(ax,"A")
    handles=[Line2D([0],[0],color=WT,lw=2,label="WT"),Line2D([0],[0],color=SCA3,lw=2,label="SCA3")]
    handles += [Line2D([0],[0],color="#666666",ls=STAGE_STYLES[s],label=s) for s in ["WT exit","Balance","SCA3 entry"]]
    ax.legend(handles=handles,ncol=2,loc="best")

    ax=fig.add_subplot(top[0,1])
    med=d.kl_balance_median.to_numpy(float);lo=d.kl_balance_q25.to_numpy(float);hi=d.kl_balance_q75.to_numpy(float)
    ax.axhline(0,color="#555555",ls=":",lw=.8);ax.fill_between(x,lo,hi,color="#7A7A7A",alpha=.15,linewidth=0);ax.plot(x,med,color="#333333",lw=1.4)
    add_stage_lines(ax);ax.set(xlabel="p",ylabel=r"$\Delta D_{KL}$",xlim=(0,1));panel(ax,"B")

    specs=[("fisher","Fisher information","#54278F","#6A51A3",True,"C"),
           ("entropy_delta_wt",r"$\Delta H$ relative to WT","#006D2C","#238B45",False,"D"),
           ("epr","Model EPR","#99000D","#CB181D",True,"E")]
    for i,(prefix,ylabel,linec,fillc,logy,tag) in enumerate(specs):
        ax=fig.add_subplot(bottom[0,i]);med=d[f"{prefix}_median"].to_numpy(float);lo=d[f"{prefix}_q25"].to_numpy(float);hi=d[f"{prefix}_q75"].to_numpy(float)
        if prefix=="entropy_delta_wt":ax.axhline(0,color="#666666",lw=.7)
        ax.fill_between(x,lo,hi,color=fillc,alpha=.16,linewidth=0);ax.plot(x,med,color=linec,lw=1.35);add_stage_lines(ax)
        if logy:ax.set_yscale("log")
        ax.set(xlabel="p",ylabel=ylabel,xlim=(0,1));panel(ax,tag)
    save(fig,"Fig4_thermodynamic_information_geometry")


def fig4():
    """Render the superseding full-cohort multiseed Figure 4."""
    subprocess.run(
        [sys.executable, str(Path(__file__).with_name("render_fig4_multiseed.py"))],
        cwd=RELEASE_ROOT,
        check=True,
    )


def figS1():
    g=read("supp_s1_dynamic_group.csv"); cells=read("supp_s1_dynamic_cells.csv"); sup=read("supp_s1_support_summary.csv")
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.3),constrained_layout=True)
    source_style={"experiment":("-","o","Experiment"),"model_best":((0,(3,2)),"s","HR model")}
    for ax,var,ylabel,tag in [(axs[0,0],"firing_rate_hz_median","Firing rate (Hz)","A"),(axs[0,1],"mean_isi_ms_median","Mean ISI (ms)","B")]:
        for group in ["WT","SCA3"]:
            for source,(ls,marker,label) in source_style.items():
                z=g[(g.group==group)&(g.source==source)].sort_values("q")
                ax.plot(z.q,z[var],color=GROUP_COLORS[group],ls=ls,marker=marker,ms=4,lw=1.2,
                        label=f"{group}, {label}")
        ax.set(xlabel="Support-restricted current coordinate, q",ylabel=ylabel);panel(ax,tag)
    axs[0,0].legend(ncol=2,fontsize=7)
    ax=axs[1,0]
    for group in ["WT","SCA3"]:
        z=cells[cells.group==group]
        for _,c in z.groupby("cell_id"):
            ax.plot(c.q,c.firing_rate_hz,color=GROUP_COLORS[group],alpha=.32,lw=.7)
            ax.scatter(c.q,c.firing_rate_hz,color=GROUP_COLORS[group],alpha=.55,s=8)
    ax.set(xlabel="q",ylabel="Experimental firing rate (Hz)");panel(ax,"C")
    ax.legend(handles=[Line2D([0],[0],color=WT,label="WT"),Line2D([0],[0],color=SCA3,label="SCA3")],loc="best")
    ax=axs[1,1];qs=[.25,.5,.75];width=.34
    for k,group in enumerate(["WT","SCA3"]):
        z=sup[sup.group==group].set_index("q").loc[qs]
        pos=np.arange(3)+(k-.5)*width
        bars=ax.bar(pos,z.fraction_supported,width=width*.92,color=GROUP_COLORS[group],label=group)
        for bar,(_,row) in zip(bars,z.iterrows()):
            ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.03,f"{int(row.n_supported)}/{int(row.n_cells)}",ha="center",va="bottom",fontsize=7)
    ax.set_xticks(range(3),["0.25","0.50","0.75"]);ax.set_ylim(0,1.13);ax.yaxis.set_major_formatter(ticker.PercentFormatter(1));ax.set(xlabel="q",ylabel="Cells within observed-current support")
    ax.legend();panel(ax,"D")
    save(fig,"FigS1_support_restricted_dynamics")


def figS2():
    loc=read("supp_thermo_noise_marker_locations.csv");shifts=read("supp_thermo_noise_marker_shifts.csv");dt=read("supp_thermo_dt_convergence.csv");align=read("fig4_marker_alignment_by_pair.csv")
    lab={"fisher_interior_peak_p_median":"Fisher peak","entropy_peak_p_median":"Entropy peak","kl_balance_p_median":"KL balance","epr_onset_p_median":"EPR onset","epr_peak_p_median":"EPR peak"}
    colors={v:c for v,c in zip(["KL balance","Fisher peak","Entropy peak","EPR onset","EPR peak"],["#333333","#54278F","#006D2C","#E6550D","#99000D"])}
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.4),constrained_layout=True)
    ax=axs[0,0]; conds=["D/2","D","2D"]; markers=["KL balance","Fisher peak","Entropy peak","EPR onset","EPR peak"]
    offsets=np.linspace(-.25,.25,len(markers))
    for k,m in enumerate(markers):
        z=loc[loc.marker.map(lab)==m].set_index("noise_condition")
        xs=[];ys=[]
        for i,c in enumerate(conds):
            if c not in z.index:continue
            row=z.loc[c];x=i+offsets[k];xs.append(x);ys.append(row["median"])
            ax.errorbar(x,row["median"],yerr=[[row["median"]-row.q25],[row.q75-row["median"]]],fmt="o",ms=3,capsize=2,lw=.8,color=colors[m],label=m if i==0 else None)
        ax.plot(xs,ys,color=colors[m],lw=.6)
    ax.set_xticks(range(3),conds);ax.set_ylim(0,1);ax.set(xlabel="Diffusion level",ylabel="Marker location, p");ax.legend(fontsize=6.7,ncol=2);panel(ax,"A")

    ax=axs[0,1];shifts["label"]=shifts.marker.map(lab);order=markers;conds2=["half","double"];offs=[-.16,.16]
    for k,c in enumerate(conds2):
        z=shifts[shifts.noise_condition==c].set_index("label")
        for i,m in enumerate(order):
            if m not in z.index:continue
            row=z.loc[m];ax.errorbar(row.median_abs_shift_p,i+offs[k],xerr=[[0],[row.q90_abs_shift_p-row.median_abs_shift_p]],fmt="o",ms=4,capsize=2,lw=.9,label={"half":"D/2","double":"2D"}[c] if i==0 else None)
    ax.set_yticks(range(len(order)),order);ax.set_xlabel(r"$|\Delta p|$: median to 90th percentile");ax.legend();panel(ax,"B")

    ax=axs[1,0];dtlab={k.replace("_median",""):v for k,v in lab.items()};dt["label"]=dt.marker.map(dtlab)
    for i,m in enumerate(order):
        z=dt[dt.label==m]
        if z.empty:continue
        row=z.iloc[0];ax.errorbar(row.median_abs_shift_p,i,xerr=[[0],[row.max_abs_shift_p-row.median_abs_shift_p]],fmt="o",color="#444444",capsize=2,lw=.9)
    ax.set_yticks(range(len(order)),order);ax.set_xlabel(r"dt 0.05 to 0.025 ms: $|\Delta p|$ median to maximum");panel(ax,"C")

    ax=axs[1,1]
    pairs=[("kl_balance_p_median","wt_exit_p_isi_weighted_median","KL balance vs WT exit"),("fisher_interior_peak_p_median","balance_p_isi_weighted_median","Fisher peak vs balance"),("entropy_peak_p_median","balance_p_isi_weighted_median","Entropy peak vs balance"),("epr_peak_p_median","sca3_entry_p_isi_weighted_median","EPR peak vs SCA3 entry")]
    mks=["o","s","^","D"]
    ax.plot([0,1],[0,1],ls="--",color="#777777",lw=.8)
    for (thermo,dyn,label),mk in zip(pairs,mks):ax.scatter(align[dyn],align[thermo],s=16,alpha=.55,label=label,marker=mk)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel="Frozen dynamical stage",ylabel="Thermodynamic marker");ax.set_aspect("equal");ax.legend(fontsize=6,loc="upper left");panel(ax,"D")
    save(fig,"FigS2_thermodynamic_robustness")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); figS1()
    print(f"Rendered four main and one supporting figure to {OUT}")
