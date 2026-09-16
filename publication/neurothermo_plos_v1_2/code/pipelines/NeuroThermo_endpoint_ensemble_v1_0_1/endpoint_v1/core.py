from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple
import json
import math

import numpy as np
import pandas as pd
import yaml


def srd(a: float, b: float) -> float:
    den = abs(a) + abs(b)
    return 0.0 if den == 0 else 2.0 * abs(a - b) / den


def robust_center_scale(x: pd.Series) -> Tuple[float, float]:
    x = pd.to_numeric(x, errors="coerce").dropna().astype(float)
    center = float(x.median())
    mad = float((x - center).abs().median())
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 1e-12:
        q1, q3 = x.quantile([0.25, 0.75])
        scale = float((q3 - q1) / 1.349)
    if not np.isfinite(scale) or scale <= 1e-12:
        scale = float(x.std(ddof=0)) if len(x) > 1 else 1.0
    if not np.isfinite(scale) or scale <= 1e-12:
        scale = 1.0
    return center, scale


def weighted_quantile(values, weights, q):
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    m = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[m], weights[m]
    if len(values) == 0:
        return np.nan
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cum = np.cumsum(weights)
    target = float(q) * cum[-1]
    return float(values[np.searchsorted(cum, target, side="left")])


def cliffs_delta(x: Iterable[float], y: Iterable[float]) -> float:
    x = np.asarray(list(x), float); y = np.asarray(list(y), float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = sum(float(a > b) for a in x for b in y)
    lt = sum(float(a < b) for a in x for b in y)
    return (gt - lt) / (len(x) * len(y))


def load_inputs(base: Path) -> Dict[str, pd.DataFrame]:
    files = {
        "fit": "v3_9_cell_fit_summary.csv",
        "alts": "final_identifiability_alternatives.csv",
        "ident": "v3_9_identifiability.csv",
        "qexp": "q_interpolated_experiment.csv",
        "qmodel": "q_interpolated_model_all_solutions.csv",
        "qrob": "q_scalar_robustness_near_optimal.csv",
        "rheo": "rheobase_refinement_all_solutions.csv",
        "animals": "animal_id_accepted_cohort.csv",
    }
    out = {}
    for k, fn in files.items():
        p = base / fn
        if not p.exists():
            raise FileNotFoundError(p)
        out[k] = pd.read_csv(p)
    return out


def primary_cells(inp: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    fit = inp["fit"].copy()
    fit = fit[(fit["final_v3_9_decision"] == "ACCEPT") & (fit["primary_support"] == "MULTI_SWEEP")].copy()
    if len(fit) != 18:
        raise ValueError(f"Expected 18 primary multi-sweep cells, got {len(fit)}")
    animals = inp["animals"][["group", "cell_id", "animal_id", "animal_id_status", "experiment_day_code"]].copy()
    fit = fit.merge(animals, on=["group", "cell_id"], how="left", validate="one_to_one")
    return fit


def q_best_table(inp, q: float, source: str) -> pd.DataFrame:
    if source == "experiment":
        df = inp["qexp"].copy()
        df = df[(df["q"] == q) & df["supported_by_observed_current_range"]].copy()
        return df
    df = inp["qmodel"].copy()
    df = df[(df["q"] == q) & df["supported_by_observed_current_range"] & (df["solution"] == "best")].copy()
    return df


def build_secure_flags(inp, cells, q: float, threshold: float) -> pd.DataFrame:
    ident = dict(zip(inp["ident"].cell_id.astype(str), inp["ident"].identifiability.astype(str)))
    rob = inp["qrob"]
    rows = []
    for cell in cells.cell_id:
        is_id = ident.get(cell) == "IDENTIFIABLE"
        rec = {"cell_id": cell, "parameter_identifiable": is_id}
        for metric in ["firing_rate_hz", "mean_isi_ms"]:
            z = rob[(rob.cell_id == cell) & (rob.q == q) & (rob.metric == metric)]
            alt_eval = bool(len(z) and z.iloc[0].status == "EVALUATED")
            alt_stable = bool(alt_eval and bool(z.iloc[0].stable_20pct))
            rec[f"{metric}_secure_q{q:g}"] = bool(is_id or alt_stable)
            rec[f"{metric}_alt_evaluated_q{q:g}"] = alt_eval
        rows.append(rec)
    return pd.DataFrame(rows)


def rheobase_secure(inp, cells, threshold: float) -> pd.DataFrame:
    ident = dict(zip(inp["ident"].cell_id.astype(str), inp["ident"].identifiability.astype(str)))
    rows = []
    for cell in cells.cell_id:
        g = inp["rheo"][inp["rheo"].cell_id == cell]
        best = g[g.solution == "best"]
        if best.empty:
            raise ValueError(f"Missing best rheobase for {cell}")
        b = float(best.iloc[0].rheobase_J)
        vals = g[g.solution == "alternative"].rheobase_J.dropna().astype(float).tolist()
        srds = [srd(b, v) for v in vals]
        is_id = ident.get(cell) == "IDENTIFIABLE"
        stable = bool(len(srds) and max(srds) <= threshold)
        rows.append({
            "cell_id": cell,
            "rheobase_J_best": b,
            "rheobase_pA_best": float(best.iloc[0].rheobase_pA),
            "n_rheobase_alternatives": len(srds),
            "rheobase_max_srd": max(srds) if srds else np.nan,
            "rheobase_secure": bool(is_id or stable),
        })
    return pd.DataFrame(rows)


def build_cell_endpoints(inp, q_primary: float, q_ext: float, threshold: float):
    cells = primary_cells(inp)
    keep = ["group", "cell_id", "animal_id", "animal_id_status", "experiment_day_code", "capacitance_pF", "b", "r", "s", "kappa_I", "identifiability", "cell_loss"]
    core = cells[keep].copy()
    rb = rheobase_secure(inp, cells, threshold)
    core = core.merge(rb, on="cell_id", how="left", validate="one_to_one")

    for q, label in [(q_primary, "q75"), (q_ext, "q50")]:
        exp = q_best_table(inp, q, "experiment")
        mod = q_best_table(inp, q, "model")
        exp_cols = ["cell_id", "firing_rate_hz", "mean_isi_ms", "median_isi_ms", "spike_count", "train_duration_ms"]
        mod_cols = ["cell_id", "firing_rate_hz", "mean_isi_ms", "median_isi_ms", "spike_count", "train_duration_ms"]
        exp = exp[exp_cols].rename(columns={c: f"exp_{label}_{c}" for c in exp_cols if c != "cell_id"})
        mod = mod[mod_cols].rename(columns={c: f"model_{label}_{c}" for c in mod_cols if c != "cell_id"})
        core = core.merge(exp, on="cell_id", how="left", validate="one_to_one")
        core = core.merge(mod, on="cell_id", how="left", validate="one_to_one")
        sec = build_secure_flags(inp, cells, q, threshold)
        # rename q strings to stable labels
        ren = {c: c.replace(f"q{q:g}", label) for c in sec.columns if c != "cell_id"}
        sec = sec.rename(columns=ren)
        core = core.merge(sec, on="cell_id", how="left", validate="one_to_one")

    # Model-vs-experiment diagnostics.
    for label in ["q75", "q50"]:
        for metric in ["firing_rate_hz", "mean_isi_ms"]:
            a = core[f"exp_{label}_{metric}"]
            b = core[f"model_{label}_{metric}"]
            core[f"model_exp_srd_{label}_{metric}"] = [srd(x, y) if pd.notna(x) and pd.notna(y) else np.nan for x, y in zip(a, b)]

    core["core_q75_secure"] = core[["rheobase_secure", "firing_rate_hz_secure_q75", "mean_isi_ms_secure_q75"]].all(axis=1)
    core["extended_q50_supported"] = core["exp_q50_firing_rate_hz"].notna()
    core["extended_q50_secure"] = core["extended_q50_supported"] & core[["firing_rate_hz_secure_q50", "mean_isi_ms_secure_q50"]].all(axis=1)
    return core


def build_solution_support(inp, cell_endpoints, q_primary: float, q_ext: float):
    fit = inp["fit"].set_index("cell_id")
    alts = inp["alts"].copy()
    alts = alts[alts.near_optimal == True].copy()  # noqa: E712
    rheo = inp["rheo"].copy()
    qm = inp["qmodel"].copy()

    rows = []
    for _, c in cell_endpoints.iterrows():
        cell = c.cell_id
        # best member
        members = [{
            "solution_key": "best", "source": "best", "b": c.b, "r": c.r, "s": c.s, "kappa_I": c.kappa_I,
            "alt_loss": np.nan,
        }]
        for _, a in alts[alts.cell_id == cell].iterrows():
            source = f"{a.parameter}_{a.side}"
            members.append({
                "solution_key": source, "source": source, "b": a.b, "r": a.r, "s": a.s, "kappa_I": a.kappa_I,
                "alt_loss": a.loss,
            })
        nmem = len(members)
        for m in members:
            rr = rheo[(rheo.cell_id == cell) & (rheo.source == m["source"])]
            if m["source"] == "best":
                rr = rheo[(rheo.cell_id == cell) & (rheo.solution == "best")]
            if rr.empty:
                raise ValueError(f"No rheobase row for {cell} {m['source']}")
            rec = {
                "group": c.group, "cell_id": cell, "animal_id": c.animal_id,
                **m,
                "rheobase_J": float(rr.iloc[0].rheobase_J),
                "rheobase_pA": float(rr.iloc[0].rheobase_pA),
                "n_solution_members_cell": nmem,
                "within_cell_support_weight": 1.0 / nmem,
            }
            for q, label in [(q_primary, "q75"), (q_ext, "q50")]:
                if m["source"] == "best":
                    z = qm[(qm.cell_id == cell) & (qm.q == q) & (qm.solution == "best")]
                else:
                    z = qm[(qm.cell_id == cell) & (qm.q == q) & (qm.solution == "alternative") & (qm.source == m["source"])]
                if z.empty:
                    for metric in ["firing_rate_hz", "mean_isi_ms", "median_isi_ms", "spike_count", "train_duration_ms"]:
                        rec[f"{label}_{metric}"] = np.nan
                else:
                    z = z.iloc[0]
                    for metric in ["firing_rate_hz", "mean_isi_ms", "median_isi_ms", "spike_count", "train_duration_ms"]:
                        rec[f"{label}_{metric}"] = z.get(metric, np.nan)
            rows.append(rec)
    out = pd.DataFrame(rows)
    n_by_group = cell_endpoints.groupby("group").cell_id.nunique().to_dict()
    out["cell_weight_within_group"] = out.group.map({g: 1.0/n for g,n in n_by_group.items()})
    out["group_support_weight"] = out["cell_weight_within_group"] * out["within_cell_support_weight"]
    return out


def add_geometry(cell_endpoints, support):
    # Core geometry intentionally excludes mean ISI because it is almost reciprocal to rate.
    best = cell_endpoints.copy()
    best["log10_rheobase_J"] = np.log10(best.rheobase_J_best.astype(float))
    best["log10_firing_rate_q75"] = np.log10(best.exp_q75_firing_rate_hz.astype(float))
    transforms = []
    for col in ["log10_rheobase_J", "log10_firing_rate_q75"]:
        center, scale = robust_center_scale(best[col])
        transforms.append({"coordinate": col, "center": center, "scale": scale, "method": "pooled best-cell median / MAD"})
        best["z_" + col] = (best[col] - center) / scale
        if col == "log10_rheobase_J":
            support[col] = np.log10(support.rheobase_J.astype(float))
        else:
            support[col] = np.log10(support.q75_firing_rate_hz.astype(float))
        support["z_" + col] = (support[col] - center) / scale
    return best, support, pd.DataFrame(transforms)


def summaries(cell_endpoints, support):
    rows = []
    for metric in ["rheobase_J_best", "exp_q75_firing_rate_hz", "exp_q75_mean_isi_ms", "exp_q50_firing_rate_hz", "exp_q50_mean_isi_ms"]:
        for group, g in cell_endpoints.groupby("group"):
            x = pd.to_numeric(g[metric], errors="coerce").dropna()
            rows.append({"level": "experimental_best_cells", "group": group, "metric": metric, "n": len(x), "median": x.median(), "q25": x.quantile(.25), "q75": x.quantile(.75)})
    # descriptive Cliff delta, SCA3 relative to WT
    for metric in ["rheobase_J_best", "exp_q75_firing_rate_hz", "exp_q75_mean_isi_ms", "exp_q50_firing_rate_hz", "exp_q50_mean_isi_ms"]:
        s = pd.to_numeric(cell_endpoints[cell_endpoints.group=="SCA3"][metric], errors="coerce").dropna()
        w = pd.to_numeric(cell_endpoints[cell_endpoints.group=="WT"][metric], errors="coerce").dropna()
        rows.append({"level": "descriptive_effect", "group": "SCA3_vs_WT", "metric": metric, "n": len(s)+len(w), "median": np.nan, "q25": np.nan, "q75": np.nan, "cliffs_delta": cliffs_delta(s,w)})
    return pd.DataFrame(rows)


def uncertainty_summary(cell_endpoints, support):
    specs = [
        ("rheobase_J", "rheobase_J_best"),
        ("q75_firing_rate_hz", "exp_q75_firing_rate_hz"),
        ("q75_mean_isi_ms", "exp_q75_mean_isi_ms"),
        ("q50_firing_rate_hz", "exp_q50_firing_rate_hz"),
        ("q50_mean_isi_ms", "exp_q50_mean_isi_ms"),
    ]
    rows = []
    for sol_metric, exp_metric in specs:
        for group in ["WT", "SCA3"]:
            ce = cell_endpoints[cell_endpoints.group==group]
            bestvals = pd.to_numeric(ce[exp_metric], errors="coerce").dropna()
            if len(bestvals) >= 2 and (bestvals > 0).all():
                lv = np.log10(bestvals)
                between_iqr = float(lv.quantile(.75)-lv.quantile(.25))
            else:
                between_iqr = np.nan
            within = []
            for cell in ce.cell_id:
                vals = pd.to_numeric(support[support.cell_id==cell][sol_metric], errors="coerce").dropna()
                vals = vals[vals>0]
                if len(vals) >= 2:
                    within.append(float(np.log10(vals).max()-np.log10(vals).min()))
            med_within = float(np.median(within)) if within else 0.0
            rows.append({
                "group": group, "metric": sol_metric,
                "between_cell_best_log10_iqr": between_iqr,
                "median_within_cell_solution_log10_span": med_within,
                "within_to_between_ratio": med_within/between_iqr if np.isfinite(between_iqr) and between_iqr>0 else np.nan,
                "n_cells_with_multiple_solutions": len(within),
            })
    return pd.DataFrame(rows)


def group_geometry(cell_geom):
    rows=[]
    zcols=["z_log10_rheobase_J","z_log10_firing_rate_q75"]
    centroids={}
    for group,g in cell_geom.groupby('group'):
        vec=g[zcols].mean().values.astype(float)
        centroids[group]=vec
        rows.append({"group":group,"n_cells":len(g),"centroid_z_log10_rheobase_J":vec[0],"centroid_z_log10_firing_rate_q75":vec[1]})
    if 'WT' in centroids and 'SCA3' in centroids:
        d=float(np.linalg.norm(centroids['SCA3']-centroids['WT']))
        rows.append({"group":"SCA3_minus_WT","n_cells":18,"centroid_z_log10_rheobase_J":centroids['SCA3'][0]-centroids['WT'][0],"centroid_z_log10_firing_rate_q75":centroids['SCA3'][1]-centroids['WT'][1],"euclidean_centroid_distance":d})
    return pd.DataFrame(rows)


def sample_for_visualization(support, n_per_group: int, seed: int):
    rng=np.random.default_rng(seed)
    rows=[]
    for group,g in support.groupby('group'):
        cells=sorted(g.cell_id.unique())
        for i in range(n_per_group):
            cell=rng.choice(cells)
            cg=g[g.cell_id==cell]
            j=int(rng.integers(0,len(cg)))
            r=cg.iloc[j].to_dict()
            r['draw']=i
            rows.append(r)
    return pd.DataFrame(rows)


def run(config_path: Path):
    cfg=yaml.safe_load(config_path.read_text())
    package_root=config_path.resolve().parents[1]
    in_dir=package_root/cfg['input_dir']
    out_dir=package_root/cfg['output_dir']
    out_dir.mkdir(parents=True,exist_ok=True)
    inp=load_inputs(in_dir)
    cells=build_cell_endpoints(inp,float(cfg['primary_q']),float(cfg['extended_q']),float(cfg['robustness_srd_threshold']))
    support=build_solution_support(inp,cells,float(cfg['primary_q']),float(cfg['extended_q']))
    cell_geom,support_geom,transform=add_geometry(cells,support)
    summary=summaries(cells,support_geom)
    uncert=uncertainty_summary(cells,support_geom)
    geom=group_geometry(cell_geom)
    vis=sample_for_visualization(support_geom,int(cfg['n_visual_samples_per_group']),int(cfg['random_seed']))

    cells.to_csv(out_dir/'endpoint_cells_full_observable.csv',index=False)
    cell_geom.to_csv(out_dir/'endpoint_cells_transition_core.csv',index=False)
    support_geom.to_csv(out_dir/'endpoint_solution_support.csv',index=False)
    transform.to_csv(out_dir/'transition_core_transform.csv',index=False)
    summary.to_csv(out_dir/'group_endpoint_summary.csv',index=False)
    uncert.to_csv(out_dir/'uncertainty_decomposition.csv',index=False)
    geom.to_csv(out_dir/'endpoint_geometry.csv',index=False)
    vis.to_csv(out_dir/'balanced_support_draws_for_visualization.csv',index=False)

    # Transition-ready subset: all actual parameter support members, all cells retained.
    cols=["group","cell_id","animal_id","solution_key","source","b","r","s","kappa_I","rheobase_J","q75_firing_rate_hz","q75_mean_isi_ms","q50_firing_rate_hz","q50_mean_isi_ms","cell_weight_within_group","within_cell_support_weight","group_support_weight","z_log10_rheobase_J","z_log10_firing_rate_q75"]
    support_geom[cols].to_csv(out_dir/'transition_ready_endpoint_support.csv',index=False)

    run_summary={
        "version":"1.0.1",
        "primary_cells":int(cells.cell_id.nunique()),
        "WT_cells":int((cells.group=='WT').sum()),
        "SCA3_cells":int((cells.group=='SCA3').sum()),
        "q75_experimental_support":int(cells.exp_q75_firing_rate_hz.notna().sum()),
        "q50_experimental_support":int(cells.exp_q50_firing_rate_hz.notna().sum()),
        "solution_support_members":int(len(support_geom)),
        "cells_with_near_optimal_members":int((support_geom.groupby('cell_id').size()>1).sum()),
        "core_q75_secure_cells":int(cells.core_q75_secure.sum()),
        "rheobase_secure_cells":int(cells.rheobase_secure.sum()),
        "q75_rate_secure_cells":int(cells.firing_rate_hz_secure_q75.sum()),
        "q75_isi_secure_cells":int(cells.mean_isi_ms_secure_q75.sum()),
        "formal_animal_level_p_values":False,
        "solution_support_weights_are_probabilities":False,
        "q50_imputation":False,
        "primary_transition_observable_coordinates":["log10_rheobase_J","log10_firing_rate_q75"],
        "full_observable_record":["rheobase_J","firing_rate_q75","mean_isi_q75","firing_rate_q50","mean_isi_q50"],
    }
    (out_dir/'RUN_SUMMARY.json').write_text(json.dumps(run_summary,indent=2),encoding='utf-8')
    return run_summary


def validate(config_path: Path):
    cfg=yaml.safe_load(config_path.read_text())
    package_root=config_path.resolve().parents[1]
    inp=load_inputs(package_root/cfg['input_dir'])
    cells=primary_cells(inp)
    q75=q_best_table(inp,float(cfg['primary_q']),'experiment')
    q50=q_best_table(inp,float(cfg['extended_q']),'experiment')
    checks={
        "version":cfg['version'],
        "primary_cells":int(len(cells)),
        "WT_cells":int((cells.group=='WT').sum()),
        "SCA3_cells":int((cells.group=='SCA3').sum()),
        "q75_supported_cells":int(q75.cell_id.nunique()),
        "q50_supported_cells":int(q50.cell_id.nunique()),
        "q50_missing_cells":sorted(set(cells.cell_id)-set(q50.cell_id)),
    }
    if checks['primary_cells']!=18 or checks['q75_supported_cells']!=18 or checks['q50_supported_cells']!=17:
        raise ValueError(checks)
    return checks
