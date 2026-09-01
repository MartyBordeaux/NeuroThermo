#!/usr/bin/env python3
"""Post-visual-QC morphology-calibrated WT--SCA3 phenotype analysis."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import platform
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences, peak_widths
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict


FEATURES = ["prominence_mV", "local_amplitude_mV", "max_dvdt_mV_per_ms",
            "abs_min_dvdt_mV_per_ms", "half_width_ms", "ahp_depth_mV"]
GROUPS = ("WT", "SCA3")


def normalize_qc_note(value):
    text = str(value).strip().lower().replace("ё", "е")
    text = text.replace("–", "-").replace("—", "-")
    text = text.replace("последниe", "последние")
    return re.sub(r"[ \t]+", " ", text)


def qc_cell_id(group, cell):
    return f"{str(group).upper()}_{int(cell):02d}"


def validate_fixed_qc(qc, sweeps):
    required = {"type", "cell", "sweep", "conclusion"}
    missing = required.difference(qc.columns)
    if missing:
        raise ValueError(f"fixed QC is missing columns: {sorted(missing)}")
    qc = qc.copy()
    qc["group"] = qc["type"].str.strip().str.upper()
    qc["cell_id"] = [qc_cell_id(g, c) for g, c in zip(qc["group"], qc["cell"])]
    qc["current_key_pA"] = qc["sweep"].astype(float).round().astype(int)
    qc["qc_note"] = qc["conclusion"].map(normalize_qc_note)
    if not set(qc.group).issubset(GROUPS):
        raise ValueError(f"unknown QC groups: {sorted(set(qc.group).difference(GROUPS))}")
    if qc.duplicated(["group", "cell_id", "current_key_pA"]).any():
        bad = qc.loc[qc.duplicated(["group", "cell_id", "current_key_pA"], keep=False),
                     ["group", "cell_id", "current_key_pA"]]
        raise ValueError(f"duplicate fixed-QC rows:\n{bad.to_string(index=False)}")
    sweep_keys = {(s["group"], s["cell_id"], int(round(s["current_pA"]))) for s in sweeps}
    qc_keys = set(map(tuple, qc[["group", "cell_id", "current_key_pA"]].itertuples(index=False, name=None)))
    missing_qc, extra_qc = sorted(sweep_keys-qc_keys), sorted(qc_keys-sweep_keys)
    if missing_qc or extra_qc:
        raise ValueError(f"fixed-QC coverage mismatch; missing={missing_qc[:10]}, extra={extra_qc[:10]}")
    return qc


def _nearest_index(block, mask, column, target, tolerance, instruction):
    candidates = block.index[np.asarray(mask, bool)]
    if not len(candidates):
        raise ValueError(f"QC instruction has no matching candidate: {instruction}")
    distance = (block.loc[candidates, column].astype(float)-float(target)).abs()
    idx = distance.idxmin()
    if float(distance.loc[idx]) > float(tolerance):
        raise ValueError(f"QC instruction is outside tolerance: {instruction}; nearest distance={distance.loc[idx]:.3g}")
    return idx


def apply_fixed_qc(events, sweeps, qc, cfg):
    """Apply frozen per-sweep visual decisions over classifier output."""
    qc = validate_fixed_qc(qc, sweeps)
    out = events.copy()
    out["algorithm_detected"] = out["detected"].astype(bool)
    out["fixed_qc_detected"] = out["algorithm_detected"]
    out["qc_note"] = ""
    out["qc_action"] = "unchanged"
    time_tol = float(cfg.get("manual_time_match_tolerance_ms", 75.0))

    for row in qc.itertuples(index=False):
        mask = ((out.group == row.group) & (out.cell_id == row.cell_id) &
                (out.current_pA.round().astype(int) == row.current_key_pA))
        idx = out.index[mask]
        note = row.qc_note
        out.loc[idx, "qc_note"] = note
        if note == "ok":
            continue
        block = out.loc[idx].sort_values("time_ms")
        instructions = [x.strip() for x in note.splitlines() if x.strip()]
        for instruction in instructions:
            current = out.loc[block.index, "fixed_qc_detected"].astype(bool)
            changed = []
            if instruction == "все rejected":
                changed = list(block.index[current])
                out.loc[block.index, "fixed_qc_detected"] = False
            elif re.fullmatch(r"(?:последние? 2|2 последних) detected - rejected", instruction):
                changed = list(block.index[current][-2:])
                if len(changed) != 2:
                    raise ValueError(f"fewer than two detected events for: {row.cell_id}, {row.sweep}, {instruction}")
                out.loc[changed, "fixed_qc_detected"] = False
            elif instruction == "последний detected - rejected":
                changed = list(block.index[current][-1:])
                if len(changed) != 1:
                    raise ValueError(f"no detected event for: {row.cell_id}, {row.sweep}, {instruction}")
                out.loc[changed, "fixed_qc_detected"] = False
            elif instruction == "последний rejected - detected":
                changed = list(block.index[~current][-1:])
                if len(changed) != 1:
                    raise ValueError(f"no rejected event for: {row.cell_id}, {row.sweep}, {instruction}")
                out.loc[changed, "fixed_qc_detected"] = True
            elif instruction == "detected около -20 мв и ниже - rejected":
                changed = list(block.index[current & (block.peak_voltage_mV <= -20.0)])
                out.loc[changed, "fixed_qc_detected"] = False
            elif (m := re.fullmatch(r"detected около (\d+(?:\.\d+)?) мс - rejected", instruction)):
                target = float(m.group(1))
                chosen = _nearest_index(block, current, "time_ms", target, time_tol, instruction)
                changed = [chosen]
                out.loc[chosen, "fixed_qc_detected"] = False
            elif (m := re.fullmatch(r"rejected около (\d+(?:\.\d+)?) мс - detected", instruction)):
                target = float(m.group(1))
                chosen = _nearest_index(block, ~current, "time_ms", target, time_tol, instruction)
                changed = [chosen]
                out.loc[chosen, "fixed_qc_detected"] = True
            elif (m := re.fullmatch(r"все после (\d+(?:\.\d+)?) мс - rejected", instruction)):
                target = float(m.group(1))
                changed = list(block.index[current & (block.time_ms > target)])
                out.loc[changed, "fixed_qc_detected"] = False
            elif instruction == "rejected после первых трех пиков - detected":
                accepted = block.index[current]
                rejected = block.index[~current]
                if len(accepted) < 3:
                    raise ValueError(f"fewer than three detected peaks for: {row.cell_id}, {row.sweep}")
                third_time = float(out.loc[accepted[2], "time_ms"])
                later = [i for i in rejected if float(out.loc[i, "time_ms"]) > third_time]
                if not later:
                    raise ValueError(f"no rejected event after first three peaks: {row.cell_id}, {row.sweep}")
                changed = [later[0]]
                out.loc[changed, "fixed_qc_detected"] = True
            else:
                raise ValueError(f"unsupported fixed-QC instruction: {instruction!r}")
            if changed:
                previous = out.loc[changed, "qc_action"].astype(str)
                out.loc[changed, "qc_action"] = [instruction if x == "unchanged" else f"{x}; {instruction}" for x in previous]
    out["qc_changed"] = out["algorithm_detected"] != out["fixed_qc_detected"]
    out["detected"] = out["fixed_qc_detected"]
    summary = qc[["group", "cell_id", "current_key_pA", "conclusion", "qc_note"]].copy()
    return out, summary


def stimulus_from_command(time_ms, command_pA):
    time_ms, command_pA = np.asarray(time_ms, float), np.asarray(command_pA, float)
    edge = max(1, min(len(command_pA) // 20, int(round(20 / np.median(np.diff(time_ms))))))
    baseline = float(np.median(command_pA[:edge]))
    deviation = np.abs(command_pA - baseline)
    mask = deviation > max(1e-9, .01 * float(np.nanmax(deviation)))
    idx = np.flatnonzero(mask)
    if not len(idx):
        raise ValueError("no command step")
    cuts = np.flatnonzero(np.diff(idx) > 1)
    runs = np.split(idx, cuts + 1)
    run = max(runs, key=len)
    start, stop = int(run[0]), int(run[-1] + 1)
    return float(time_ms[start]), float(time_ms[min(stop, len(time_ms)-1)]), float(np.median(command_pA[start:stop]) - baseline)


def cell_id(group, path):
    return f"{group}_{path.stem.split('_')[-1]}"


def find_cc(root):
    files = sorted({*root.rglob("cc_*.abf"), *root.rglob("CC_*.abf")})
    return [p for p in files if "capa" not in p.name.lower()]


def load_sweeps(group, path, cfg):
    import pyabf
    abf = pyabf.ABF(str(path))
    result = []
    for sweep in map(int, abf.sweepList):
        abf.setSweep(sweep, channel=int(cfg["command_channel"]))
        t = np.asarray(abf.sweepX, float) * 1000
        command = np.asarray(abf.sweepC, float)
        try:
            onset, offset, current = stimulus_from_command(t, command)
        except ValueError:
            abf.setSweep(int(abf.sweepList[-1]), channel=int(cfg["command_channel"]))
            rt, rc = np.asarray(abf.sweepX, float)*1000, np.asarray(abf.sweepC, float)
            onset, offset, _ = stimulus_from_command(rt, rc)
            abf.setSweep(sweep, channel=int(cfg["command_channel"]))
            t, command = np.asarray(abf.sweepX, float)*1000, np.asarray(abf.sweepC, float)
            pre, stim = command[t < onset], command[(t >= onset) & (t < offset)]
            current = float(np.median(stim) - np.median(pre))
        abf.setSweep(sweep, channel=int(cfg["voltage_channel"]))
        v = np.asarray(abf.sweepY, float)
        result.append(dict(group=group, cell_id=cell_id(group, path), abf_path=str(path.resolve()),
                           sweep_index=sweep, current_pA=round(current, 1), onset_ms=onset,
                           offset_ms=offset, time_ms=t, voltage_mV=v))
    return result


def window_stat(t, values, lo, hi, fn, default):
    x = values[(t >= lo) & (t <= hi)]
    return float(fn(x)) if len(x) else float(default)


def candidate_features(sweep, cfg):
    t, v = sweep["time_ms"], sweep["voltage_mV"]
    dt = float(np.median(np.diff(t)))
    distance = max(1, int(round(float(cfg["minimum_peak_distance_ms"]) / dt)))
    peaks, _ = find_peaks(v, prominence=float(cfg["candidate_prominence_mV"]),
                          height=float(cfg["candidate_height_mV"]), distance=distance)
    if not len(peaks):
        return pd.DataFrame()
    prom = peak_prominences(v, peaks)[0]
    widths = peak_widths(v, peaks, rel_height=.5)[0] * dt
    dvdt = np.gradient(v, t)
    rows = []
    for p, pr, width in zip(peaks, prom, widths):
        tp = float(t[p])
        baseline = window_stat(t, v, tp-5, tp-2, np.median, v[p])
        rise = window_stat(t, dvdt, tp-3, tp, np.max, 0)
        fall = window_stat(t, dvdt, tp, tp+4, np.min, 0)
        ahp = baseline - window_stat(t, v, tp+1, tp+12, np.min, baseline)
        row = {k: sweep[k] for k in ("group", "cell_id", "abf_path", "sweep_index", "current_pA", "onset_ms", "offset_ms")}
        row.update(peak_index=int(p), time_ms=tp, peak_voltage_mV=float(v[p]),
                   prominence_mV=float(pr), local_amplitude_mV=float(v[p]-baseline),
                   max_dvdt_mV_per_ms=float(rise), abs_min_dvdt_mV_per_ms=float(abs(fall)),
                   half_width_ms=float(width), ahp_depth_mV=float(ahp))
        rows.append(row)
    return pd.DataFrame(rows)


def marker_class(frame):
    strict = (frame.prominence_mV >= 10) & (frame.peak_voltage_mV >= -10)
    primary = (frame.prominence_mV >= 10) & (frame.peak_voltage_mV >= -20)
    relaxed = (frame.prominence_mV >= 8) & (frame.peak_voltage_mV >= -30)
    return np.select([strict, primary, relaxed], ["strict", "primary", "relaxed"], default="broad")


def build_training(events, rules):
    merged = events.merge(rules[["group", "cell_id", "training_rule"]], on=["group", "cell_id"], how="left")
    merged["marker_class"] = marker_class(merged)
    rank = {"strict": 3, "primary": 2, "relaxed": 1, "broad": 0}
    reviewed = merged.marker_class.ne("broad") & merged.training_rule.isin(["strict", "primary", "relaxed"])
    train = merged[reviewed].copy()
    train["label"] = [int(rank[m] >= rank[r]) for m, r in zip(train.marker_class, train.training_rule)]
    return train


def model_factory(cfg):
    return RandomForestClassifier(n_estimators=int(cfg["rf_trees"]), min_samples_leaf=3,
                                  class_weight="balanced", max_features="sqrt",
                                  random_state=int(cfg["random_seed"]), n_jobs=-1)


def calibrate(train, cfg):
    clean = train.dropna(subset=FEATURES + ["label"]).copy()
    if clean.label.nunique() != 2 or clean.cell_id.nunique() < 3:
        raise RuntimeError("manual QC anchors do not contain both classes across enough cells")
    logo = LeaveOneGroupOut()
    prob = cross_val_predict(model_factory(cfg), clean[FEATURES], clean.label.astype(int),
                             groups=clean.cell_id, cv=logo, method="predict_proba", n_jobs=1)[:, 1]
    rows = []
    for threshold in map(float, cfg["probability_threshold_grid"]):
        pred = prob >= threshold
        precision, recall, f1, _ = precision_recall_fscore_support(clean.label, pred, average="binary", zero_division=0)
        rows.append(dict(threshold=threshold, precision=precision, recall=recall, f1=f1))
    scores = pd.DataFrame(rows)
    best = scores.sort_values(["f1", "precision", "threshold"], ascending=[False, False, False]).iloc[0]
    model = model_factory(cfg).fit(clean[FEATURES], clean.label.astype(int))
    audit = {"selected_threshold": float(best.threshold), "n_events": int(len(clean)),
             "n_cells": int(clean.cell_id.nunique()), "positive_fraction": float(clean.label.mean()),
             "leave_one_cell_out_auc": float(roc_auc_score(clean.label, prob))}
    clean["cv_probability"] = prob
    return model, float(best.threshold), scores, clean, audit


def sustained_mask(times, max_isi_ms, minimum_spikes):
    times = np.asarray(times, float)
    out = np.zeros(len(times), bool)
    if len(times) < minimum_spikes:
        return out
    gaps = np.diff(times) > max_isi_ms
    starts = np.r_[0, np.flatnonzero(gaps)+1]
    stops = np.r_[np.flatnonzero(gaps)+1, len(times)]
    for start, stop in zip(starts, stops):
        if stop-start >= minimum_spikes:
            out[start:stop] = True
    return out


def apply_rules(sweeps, detected, rules, cfg):
    rule_map = {(r.group, r.cell_id): r for r in rules.itertuples(index=False)}
    rows = []
    for s in sweeps:
        rule = rule_map.get((s["group"], s["cell_id"]))
        lo, hi = s["onset_ms"] + cfg["stimulus_edge_trim_ms"], s["offset_ms"] - cfg["stimulus_edge_trim_ms"]
        ev = detected[(detected.abf_path == s["abf_path"]) & (detected.sweep_index == s["sweep_index"]) &
                      detected.detected & detected.time_ms.between(lo, hi, inclusive="left")].sort_values("time_ms")
        forced_zero = bool(rule and pd.notna(rule.valid_zero_below_pA) and s["current_pA"] < float(rule.valid_zero_below_pA))
        above_valid = bool(rule and pd.notna(rule.maximum_valid_current_pA) and s["current_pA"] > float(rule.maximum_valid_current_pA))
        # Fixed event-level QC is authoritative. Legacy forced-zero and valid-range
        # rules are retained as audit/sensitivity flags, not applied to primary counts.
        times = ev.time_ms.to_numpy(float)
        keep_train = sustained_mask(times, float(cfg["sustained_train_max_isi_ms"]), int(cfg["minimum_train_spikes"]))
        duration = max((hi-lo)/1000, 1e-9)
        early_stop, late_start = lo + (hi-lo)/3, hi - (hi-lo)/3
        early = int(np.sum((times >= lo) & (times < early_stop)))
        late = int(np.sum((times >= late_start) & (times < hi)))
        rows.append(dict(group=s["group"], cell_id=s["cell_id"], abf_path=s["abf_path"], sweep_index=s["sweep_index"],
                         current_pA=s["current_pA"], onset_ms=s["onset_ms"], offset_ms=s["offset_ms"],
                         main_cell_include=int(rule.main_include) if rule else 1,
                         phenotype=rule.phenotype if rule else "unspecified", forced_valid_zero=int(forced_zero),
                         excluded_above_valid_range=int(above_valid), primary_current_range=int(cfg["primary_current_min_pA"] <= s["current_pA"] <= cfg["primary_current_max_pA"]),
                         total_spike_count=int(len(times)), sustained_spike_count=int(keep_train.sum()),
                         total_rate_hz=float(len(times)/duration), sustained_rate_hz=float(keep_train.sum()/duration),
                         first_spike_latency_ms=float(times[0]-s["onset_ms"]) if len(times) else np.nan,
                         early_spike_count=early, late_spike_count=late,
                         depolarization_block_candidate=int(early >= 3 and late/max(early,1) <= .25)))
    return pd.DataFrame(rows)


def cohort_frame(features, cohort):
    base = features[features.primary_current_range == 1].copy()
    if cohort == "fixed_qc_all":
        return base
    if cohort == "conservative":
        return base[(base.main_cell_include == 1) & (base.excluded_above_valid_range == 0)]
    if cohort == "exclude_low_amplitude_cell":
        return base[(base.main_cell_include == 1) & (base.excluded_above_valid_range == 0) &
                    (base.cell_id != "SCA3_02")]
    raise ValueError(cohort)


def permutation_p(x, y, rng, reps):
    x, y = np.asarray(x,float), np.asarray(y,float)
    observed, pool, nx = np.median(y)-np.median(x), np.r_[x,y], len(x)
    total = math.comb(len(pool), nx)
    iterator = itertools.combinations(range(len(pool)), nx) if total <= reps else (rng.permutation(len(pool))[:nx] for _ in range(reps))
    extreme = n = 0
    for chosen in iterator:
        mask = np.zeros(len(pool), bool); mask[list(chosen)] = True
        stat = np.median(pool[~mask])-np.median(pool[mask])
        extreme += abs(stat) >= abs(observed)-1e-12; n += 1
    return extreme/n


def summarize_cohort(features, cohort, cfg):
    source = cohort_frame(features, cohort)
    cells = source.groupby(["group","cell_id"], as_index=False).agg(
        mean_sustained_rate_hz=("sustained_rate_hz","mean"),
        maximum_sustained_rate_hz=("sustained_rate_hz","max"),
        mean_total_rate_hz=("total_rate_hz","mean"),
        block_fraction=("depolarization_block_candidate","mean"))
    recruit = []
    for (group, cell), block in source.groupby(["group","cell_id"]):
        positive = block[block.sustained_spike_count >= 2]
        recruit.append(dict(group=group, cell_id=cell, recruitment_current_pA=float(positive.current_pA.min()) if len(positive) else np.nan))
    cells = cells.merge(pd.DataFrame(recruit), on=["group","cell_id"], how="left")
    rng = np.random.default_rng(int(cfg["random_seed"]))
    tests=[]
    for endpoint in ["mean_sustained_rate_hz","maximum_sustained_rate_hz","mean_total_rate_hz","recruitment_current_pA","block_fraction"]:
        wt=cells.loc[cells.group=="WT",endpoint].dropna().to_numpy(); sc=cells.loc[cells.group=="SCA3",endpoint].dropna().to_numpy()
        tests.append(dict(cohort=cohort, endpoint=endpoint, n_WT=len(wt), n_SCA3=len(sc),
                          WT_median=np.median(wt) if len(wt) else np.nan, SCA3_median=np.median(sc) if len(sc) else np.nan,
                          difference_SCA3_minus_WT=np.median(sc)-np.median(wt) if len(wt) and len(sc) else np.nan,
                          permutation_p=permutation_p(wt,sc,rng,int(cfg["permutation_replicates"])) if min(len(wt),len(sc))>=3 else np.nan))
    return cells, pd.DataFrame(tests), source


def count_model(source, cohort):
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    block=source.copy(); block["current_z"]=(block.current_pA-block.current_pA.mean())/block.current_pA.std(ddof=0)
    formula="sustained_spike_count ~ C(group, Treatment(reference='WT')) * (current_z + I(current_z ** 2))"
    pois=smf.glm(formula,data=block,family=sm.families.Poisson()).fit()
    dispersion=float(np.sum(pois.resid_pearson**2)/max(pois.df_resid,1)); alpha=max((dispersion-1)/max(pois.fittedvalues.mean(),1e-6),1e-6)
    nb=smf.glm(formula,data=block,family=sm.families.NegativeBinomial(alpha=alpha)).fit(cov_type="cluster",cov_kwds={"groups":block.cell_id})
    names=list(nb.params.index); idx=[i for i,n in enumerate(names) if "SCA3" in n]; R=np.zeros((len(idx),len(names)))
    for r,c in enumerate(idx): R[r,c]=1
    p=float(nb.wald_test(R,scalar=True).pvalue) if idx else np.nan
    return dict(cohort=cohort,n_sweeps=len(block),n_cells=block.cell_id.nunique(),poisson_dispersion=dispersion,negative_binomial_alpha=alpha,group_joint_p=p)


def figures(sweeps, events, features, output):
    trace_dir=output/"03_visual_audit"; trace_dir.mkdir(parents=True,exist_ok=True)
    lookup={(s["abf_path"],s["sweep_index"]):s for s in sweeps}
    for (group,cell), block in features.groupby(["group","cell_id"]):
        target=trace_dir/f"{cell}_detections.pdf"
        with PdfPages(target) as pdf:
            for _,f in block.sort_values(["current_pA","sweep_index"]).iterrows():
                s=lookup[(f.abf_path,int(f.sweep_index))]; ev=events[(events.abf_path==f.abf_path)&(events.sweep_index==f.sweep_index)]
                fig,ax=plt.subplots(figsize=(11,3)); ax.plot(s["time_ms"],s["voltage_mV"],color=".2",lw=.6)
                yes=ev[ev.detected]; no=ev[~ev.detected]
                ax.scatter(no.time_ms,no.peak_voltage_mV,marker="x",color="#999999",s=12,label="rejected")
                ax.scatter(yes.time_ms,yes.peak_voltage_mV,facecolors="none",edgecolors="#d73027",s=28,label="detected")
                add=ev[(~ev.algorithm_detected)&ev.fixed_qc_detected]
                remove=ev[ev.algorithm_detected&(~ev.fixed_qc_detected)]
                ax.scatter(add.time_ms,add.peak_voltage_mV,marker="+",color="#1b9e77",s=50,label="manual add")
                ax.scatter(remove.time_ms,remove.peak_voltage_mV,marker="s",facecolors="none",edgecolors="#7570b3",s=38,label="manual remove")
                ax.axvspan(s["onset_ms"],s["offset_ms"],color="#eeeeee",zorder=-3)
                ax.set(title=f"{cell} | {f.current_pA:g} pA | total={f.total_spike_count}, sustained={f.sustained_spike_count}",xlabel="Time (ms)",ylabel="mV"); ax.legend(fontsize=7)
                fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    main=cohort_frame(features,"fixed_qc_all")
    fig,ax=plt.subplots(figsize=(7,5))
    per=main.groupby(["group","cell_id","current_pA"],as_index=False).sustained_rate_hz.mean()
    colors={"WT":"#2166ac","SCA3":"#b2182b"}
    for group in GROUPS:
        for _,b in per[per.group==group].groupby("cell_id"):
            ax.plot(b.current_pA,b.sustained_rate_hz,color=colors[group],alpha=.25,lw=.8)
        g=per[per.group==group].groupby("current_pA",as_index=False).sustained_rate_hz.mean()
        ax.plot(g.current_pA,g.sustained_rate_hz,"o-",color=colors[group],lw=2,label=group)
    ax.set(xlabel="Injected current (pA)",ylabel="Sustained firing rate (Hz)"); ax.legend(); fig.tight_layout()
    fig.savefig(output/"03_visual_audit"/"FI_main_cohort.pdf"); fig.savefig(output/"03_visual_audit"/"FI_main_cohort.png",dpi=250); plt.close(fig)


def main():
    p=argparse.ArgumentParser(); p.add_argument("--wt-root",required=True); p.add_argument("--sca3-root",required=True); p.add_argument("--config",default="config.json"); p.add_argument("--qc-rules",default="qc_rules.csv"); p.add_argument("--fixed-qc",default="qc2.csv"); p.add_argument("--output",default="results_stage1_qc_fixed"); a=p.parse_args()
    cfg=json.loads(Path(a.config).read_text()); rules=pd.read_csv(a.qc_rules); output=Path(a.output); output.mkdir(parents=True,exist_ok=True)
    sweeps=[]
    for group,root in [("WT",Path(a.wt_root)),("SCA3",Path(a.sca3_root))]:
        files=find_cc(root)
        if not files: raise SystemExit(f"no current-clamp ABFs in {root}")
        for path in files: print(f"extract {group} {path.name}",flush=True); sweeps.extend(load_sweeps(group,path,cfg))
    tables=[candidate_features(s,cfg) for s in sweeps]; events=pd.concat([x for x in tables if len(x)],ignore_index=True)
    train=build_training(events,rules); model,threshold,grid,cv,audit=calibrate(train,cfg)
    events["spike_probability"]=model.predict_proba(events[FEATURES])[:,1]; events["detected"]=events.spike_probability>=threshold
    events,qc_summary=apply_fixed_qc(events,sweeps,pd.read_csv(a.fixed_qc),cfg)
    features=apply_rules(sweeps,events,rules,cfg)
    (output/"00_calibration").mkdir(exist_ok=True); (output/"01_features").mkdir(exist_ok=True); (output/"02_statistics").mkdir(exist_ok=True)
    events.to_csv(output/"00_calibration"/"candidate_events_with_predictions.csv",index=False); cv.to_csv(output/"00_calibration"/"training_events_leave_one_cell_out.csv",index=False); grid.to_csv(output/"00_calibration"/"threshold_selection.csv",index=False)
    features.to_csv(output/"01_features"/"sweep_features.csv",index=False); rules.to_csv(output/"01_features"/"formal_qc_rules.csv",index=False)
    qc_summary.to_csv(output/"01_features"/"fixed_manual_qc.csv",index=False)
    override_summary=(events.groupby(["group","cell_id"],as_index=False)
                      .agg(candidate_events=("detected","size"),
                           algorithm_spikes=("algorithm_detected","sum"),
                           fixed_qc_spikes=("fixed_qc_detected","sum"),
                           changed_events=("qc_changed","sum")))
    override_summary["manual_additions"]=(events[(~events.algorithm_detected)&events.fixed_qc_detected]
                                          .groupby(["group","cell_id"]).size()
                                          .reindex(pd.MultiIndex.from_frame(override_summary[["group","cell_id"]]),fill_value=0).to_numpy())
    override_summary["manual_removals"]=(events[events.algorithm_detected&(~events.fixed_qc_detected)]
                                         .groupby(["group","cell_id"]).size()
                                         .reindex(pd.MultiIndex.from_frame(override_summary[["group","cell_id"]]),fill_value=0).to_numpy())
    override_summary.to_csv(output/"01_features"/"qc_override_summary_by_cell.csv",index=False)
    all_cells=[]; all_tests=[]; models=[]
    for cohort in ["fixed_qc_all","conservative","exclude_low_amplitude_cell"]:
        cells,tests,source=summarize_cohort(features,cohort,cfg); cells["cohort"]=cohort; all_cells.append(cells); all_tests.append(tests)
        try: models.append(count_model(source,cohort))
        except Exception as exc: models.append(dict(cohort=cohort,error=f"{type(exc).__name__}: {exc}"))
    pd.concat(all_cells,ignore_index=True).to_csv(output/"01_features"/"cell_phenotypes_all_cohorts.csv",index=False); pd.concat(all_tests,ignore_index=True).to_csv(output/"02_statistics"/"scalar_tests_all_cohorts.csv",index=False); pd.DataFrame(models).to_csv(output/"02_statistics"/"negative_binomial_models.csv",index=False)
    figures(sweeps,events,features,output)
    audit.update(feature_columns=FEATURES,absolute_peak_voltage_used_for_classifier=False,
                 primary_cohort="fixed_qc_all",primary_cohort_exclusions=[],
                 sensitivity_cohorts=["conservative","exclude_low_amplitude_cell"],
                 fixed_qc_rows=int(len(qc_summary)),fixed_qc_changed_events=int(events.qc_changed.sum()),
                 fixed_qc_sha256=hashlib.sha256(Path(a.fixed_qc).read_bytes()).hexdigest(),
                 python=sys.version,platform=platform.platform())
    (output/"manifest.json").write_text(json.dumps(audit,indent=2)+"\n")
    report=["# Fixed-manual-QC Stage 1", "", f"Selected morphology probability threshold: **{threshold:.2f}**.", f"Leave-one-cell-out AUC: **{audit['leave_one_cell_out_auc']:.3f}**.", f"Manual event decisions changed: **{audit['fixed_qc_changed_events']}**.", "", "The primary cohort contains all 25 cells reviewed sweep by sweep in qc2.csv.", "The prior conservative cohort is retained as a prespecified sensitivity analysis.", "Primary inference remains cell-level until animal identifiers are supplied."]
    (output/"REPORT.md").write_text("\n".join(report)+"\n")
    print(f"completed: {output.resolve()}",flush=True)


if __name__ == "__main__": main()
