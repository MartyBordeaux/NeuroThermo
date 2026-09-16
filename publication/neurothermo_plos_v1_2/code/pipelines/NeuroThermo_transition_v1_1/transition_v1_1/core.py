from __future__ import annotations
from pathlib import Path
import json, hashlib, zipfile, tempfile, shutil
import numpy as np
import pandas as pd

from .geometry import fit_reference, reference_from_v1_tables, project_native, persistent_crossing, weighted_quantile, interp_at

REQUIRED = [
    "transition_paths.csv",
    "transition_pair_scenarios.csv",
    "transition_protocol_endpoint_states.csv",
    "transition_projection_reference.csv",
    "transition_projection_transform.csv",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def resolve_input(cfg, config_path):
    base = Path(config_path).resolve().parent.parent
    candidates = cfg["input"].get("candidates", [])
    checked = []
    for raw in candidates:
        p = (base / raw).resolve() if not Path(raw).is_absolute() else Path(raw)
        checked.append(str(p))
        if p.is_dir() and all((p/f).exists() for f in REQUIRED):
            return {"kind":"dir", "path":p, "root":p, "temporary":None}
        if p.is_file() and p.suffix.lower() == ".zip":
            td = Path(tempfile.mkdtemp(prefix="transition_v1_1_"))
            with zipfile.ZipFile(p) as z:
                z.extractall(td)
            roots = [td] + [x for x in td.iterdir() if x.is_dir()]
            for r in roots:
                if all((r/f).exists() for f in REQUIRED):
                    return {"kind":"zip", "path":p, "root":r, "temporary":td}
            shutil.rmtree(td, ignore_errors=True)
    raise FileNotFoundError("No valid transition v1.0 result directory/zip found. Checked: " + "; ".join(checked))


def cleanup_input(info):
    if info.get("temporary"):
        shutil.rmtree(info["temporary"], ignore_errors=True)


def load_csvs(root):
    return {name: pd.read_csv(root/name) for name in REQUIRED}


def build_refs(old, frozen_ref, cfg):
    # ISI is intentionally frozen from v1.0 because its observable definition was already latency-invariant and consistent.
    oldref = old["transition_projection_reference.csv"]
    oldtr = old["transition_projection_transform.csv"]
    rr = oldref[oldref.projection.eq("isi")].iloc[0]
    tt = oldtr[oldtr.projection.eq("isi")].copy()
    isi_ref = reference_from_v1_tables(rr, tt)

    ep = old["transition_protocol_endpoint_states.csv"]
    best_rheo = ep[ep.source.eq("best")][["cell_id","protocol_rheobase_J"]].drop_duplicates("cell_id")
    active_cells = frozen_ref.merge(best_rheo, on="cell_id", how="left")
    active_ref = fit_reference(
        active_cells, "protocol_rheobase_J", "exp_active_rate_q75_hz",
        secure_col="core_q75_secure",
        wt_q=float(cfg["staging"]["wt_exit_quantile"]),
        sca_q=float(cfg["staging"]["sca3_entry_quantile"]),
    )
    return isi_ref, active_ref, active_cells


def reproject_paths(paths, isi_ref, active_ref):
    d = paths.copy()
    # Semantic correction: support_rate_hz in v1.0 was computed in a fixed window starting at the first model spike.
    d["active_support_rate_hz"] = pd.to_numeric(d["support_rate_hz"], errors="coerce")
    Aact, Oact = project_native(d, active_ref, "rheobase_J", "active_support_rate_hz")
    Aisi, Oisi = project_native(d, isi_ref, "rheobase_J", "mean_isi_ms")
    d["A_active"] = Aact
    d["orth_active"] = Oact
    d["A_isi_v1_1"] = Aisi
    d["orth_isi_v1_1"] = Oisi
    return d


def marker_table(paths, scenarios, isi_ref, active_ref, cfg):
    persist = int(cfg["staging"]["persistence_points"])
    tol = float(cfg["staging"]["projection_agreement_tolerance"])
    smap = scenarios.set_index("scenario_id")
    rows=[]
    for (sid, fam), g in paths.groupby(["scenario_id","path_family"], sort=False):
        g=g.sort_values("path_progress"); p=g.path_progress.to_numpy(float)
        sc=smap.loc[int(sid)]
        row={"scenario_id":int(sid),"path_family":fam}
        for k in ["biological_pair_key","wt_cell_id","sca_cell_id","wt_solution_key","sca_solution_key","within_pair_support_weight","scenario_weight","biological_pair_weight"]:
            row[k]=sc[k]
        for nm,col,ref in [("isi","A_isi_v1_1",isi_ref),("active","A_active",active_ref)]:
            A=pd.to_numeric(g[col],errors="coerce").to_numpy(float)
            row[f"wt_exit_p_{nm}"]=persistent_crossing(p,A,float(ref["wt_exit_A_threshold"]),persist)
            row[f"balance_p_{nm}"]=persistent_crossing(p,A,0.5,persist)
            row[f"sca3_entry_p_{nm}"]=persistent_crossing(p,A,float(ref["sca3_entry_A_threshold"]),persist)
            row[f"start_A_{nm}"]=A[0] if len(A) else np.nan
            row[f"end_A_{nm}"]=A[-1] if len(A) else np.nan
            row[f"start_inside_wt_{nm}"]=bool(np.isfinite(A[0]) and A[0] <= ref["wt_exit_A_threshold"]) if len(A) else False
            row[f"end_inside_sca3_{nm}"]=bool(np.isfinite(A[-1]) and A[-1] >= ref["sca3_entry_A_threshold"]) if len(A) else False
        for stage in ["wt_exit","balance","sca3_entry"]:
            a=row[f"{stage}_p_active"]; i=row[f"{stage}_p_isi"]
            gap=abs(a-i) if np.isfinite(a) and np.isfinite(i) else np.nan
            row[f"{stage}_active_isi_gap"]=gap
            row[f"{stage}_dual_agrees"]=bool(np.isfinite(gap) and gap <= tol)
            row[f"{stage}_dual_consensus_p"]=0.5*(a+i) if row[f"{stage}_dual_agrees"] else np.nan
        # Occupancy is descriptive only; evaluate it at the primary ISI-derived stage positions.
        occ=pd.to_numeric(g["occupancy_fraction"],errors="coerce").to_numpy(float)
        for stage in ["wt_exit","balance","sca3_entry"]:
            row[f"occupancy_at_{stage}_isi"]=interp_at(p,occ,row[f"{stage}_p_isi"])
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(markers, frozen_ref):
    secure=set(frozen_ref.loc[frozen_ref.core_q75_secure.fillna(False).astype(bool),"cell_id"].astype(str))
    metrics=[]
    for s in ["wt_exit","balance","sca3_entry"]:
        metrics += [f"{s}_p_isi",f"{s}_p_active",f"{s}_dual_consensus_p"]
    rows=[]
    for (pair,fam),g in markers.groupby(["biological_pair_key","path_family"]):
        w=pd.to_numeric(g.within_pair_support_weight,errors="coerce").to_numpy(float)
        r={"biological_pair_key":pair,"path_family":fam,"wt_cell_id":g.wt_cell_id.iloc[0],"sca_cell_id":g.sca_cell_id.iloc[0],
           "both_core_secure":bool(str(g.wt_cell_id.iloc[0]) in secure and str(g.sca_cell_id.iloc[0]) in secure),"n_support_scenarios":len(g)}
        for m in metrics:
            vals=pd.to_numeric(g[m],errors="coerce").to_numpy(float)
            r[m+"_weighted_median"]=weighted_quantile(vals,w,.5)
            r[m+"_weighted_q25"]=weighted_quantile(vals,w,.25)
            r[m+"_weighted_q75"]=weighted_quantile(vals,w,.75)
            r[m+"_support_fraction"]=float(np.sum(w[np.isfinite(vals)])/np.sum(w)) if np.sum(w)>0 else np.nan
        for s in ["wt_exit","balance","sca3_entry"]:
            vals=pd.to_numeric(g[f"{s}_active_isi_gap"],errors="coerce").to_numpy(float)
            r[f"{s}_active_isi_gap_weighted_median"]=weighted_quantile(vals,w,.5)
            agree=g[f"{s}_dual_agrees"].astype(bool).to_numpy()
            r[f"{s}_dual_agreement_weight"]=float(np.sum(w[agree])/np.sum(w)) if np.sum(w)>0 else np.nan
            vals_occ=pd.to_numeric(g[f"occupancy_at_{s}_isi"],errors="coerce").to_numpy(float)
            r[f"occupancy_at_{s}_isi_weighted_median"]=weighted_quantile(vals_occ,w,.5)
        rows.append(r)
    pairdf=pd.DataFrame(rows)
    overall=[]
    for fam,g in pairdf.groupby("path_family"):
        for subset,x in [("all_pairs",g),("core_secure_pairs",g[g.both_core_secure])]:
            for m in metrics:
                v=pd.to_numeric(x[m+"_weighted_median"],errors="coerce").to_numpy(float); v=v[np.isfinite(v)]
                overall.append({"path_family":fam,"subset":subset,"metric":m,"n_biological_pairs_total":len(x),"n_pairs_with_marker":len(v),
                    "median":float(np.median(v)) if len(v) else np.nan,"q25":float(np.quantile(v,.25)) if len(v) else np.nan,"q75":float(np.quantile(v,.75)) if len(v) else np.nan})
    return pairdf,pd.DataFrame(overall)


def summarize_curves(paths):
    metrics=["A_isi_v1_1","A_active","rheobase_J","active_support_rate_hz","mean_isi_ms","occupancy_fraction"]
    pair=[]
    for (pk,fam,p),g in paths.groupby(["biological_pair_key","path_family","path_progress"]):
        w=pd.to_numeric(g.within_pair_support_weight,errors="coerce").to_numpy(float)
        r={"biological_pair_key":pk,"path_family":fam,"path_progress":p,"wt_cell_id":g.wt_cell_id.iloc[0],"sca_cell_id":g.sca_cell_id.iloc[0]}
        for m in metrics:
            r[m+"_weighted_median"]=weighted_quantile(pd.to_numeric(g[m],errors="coerce"),w,.5)
        pair.append(r)
    pair=pd.DataFrame(pair)
    ens=[]
    for (fam,p),g in pair.groupby(["path_family","path_progress"]):
        r={"path_family":fam,"path_progress":p,"n_biological_pairs":len(g)}
        for m in metrics:
            v=pd.to_numeric(g[m+"_weighted_median"],errors="coerce").to_numpy(float);v=v[np.isfinite(v)]
            r[m+"_median"]=float(np.median(v)) if len(v) else np.nan
            r[m+"_q25"]=float(np.quantile(v,.25)) if len(v) else np.nan
            r[m+"_q75"]=float(np.quantile(v,.75)) if len(v) else np.nan
        ens.append(r)
    return pair,pd.DataFrame(ens)


def reference_tables(isi_ref,active_ref):
    rows=[]; trs=[]
    for name,ref in [("isi_primary_v1_0_frozen",isi_ref),("active_rate_experimental_v2_1",active_ref)]:
        rows.append({"projection":name,"centroid_distance":ref["centroid_distance"],"wt_exit_A_threshold":ref["wt_exit_A_threshold"],"sca3_entry_A_threshold":ref["sca3_entry_A_threshold"],
            "cloud_overlap":ref["cloud_overlap"],"corridor_radius_q90":ref["corridor_radius_q90"],"wt_centroid_0":ref["cwt"][0],"wt_centroid_1":ref["cwt"][1],"sca3_centroid_0":ref["csc"][0],"sca3_centroid_1":ref["csc"][1]})
        coords=["log10_rheobase","log10_isi"] if name.startswith("isi") else ["log10_rheobase","log10_active_rate"]
        for i,c in enumerate(coords):
            trs.append({"projection":name,"coordinate":c,"center":ref["center"][i],"scale":ref["scale"][i]})
    return pd.DataFrame(rows),pd.DataFrame(trs)


def correction_audit(old, paths_new, active_cells, active_ref):
    rows=[]
    # ISI reprojection must be numerically identical to v1.0 A_isi.
    oldA=pd.to_numeric(old["transition_paths.csv"]["A_isi"],errors="coerce").to_numpy(float)
    newA=pd.to_numeric(paths_new["A_isi_v1_1"],errors="coerce").to_numpy(float)
    diff=np.abs(oldA-newA); finite=np.isfinite(diff)
    rows.append({"check":"ISI reprojection vs v1.0 stored A_isi","value":float(np.nanmax(diff)) if finite.any() else np.nan,"criterion":"max_abs_diff <= 1e-10","pass":bool(finite.any() and np.nanmax(diff)<=1e-10)})
    # Semantic equality: old support_rate is copied, not recomputed.
    a=pd.to_numeric(old["transition_paths.csv"]["support_rate_hz"],errors="coerce").to_numpy(float)
    b=pd.to_numeric(paths_new["active_support_rate_hz"],errors="coerce").to_numpy(float)
    eq=np.allclose(a,b,equal_nan=True,rtol=0,atol=0)
    rows.append({"check":"support_rate_hz renamed to active_support_rate_hz","value":1.0 if eq else 0.0,"criterion":"exact equality","pass":bool(eq)})
    # Endpoint validation against v2.1 model active rate for best solutions.
    ep=old["transition_protocol_endpoint_states.csv"]
    b_ep=ep[ep.source.eq("best")][["cell_id","support_rate_hz"]].merge(active_cells[["cell_id","model_active_rate_q75_hz"]],on="cell_id",how="left")
    x=b_ep.support_rate_hz.to_numpy(float); y=b_ep.model_active_rate_q75_hz.to_numpy(float)
    srd=2*np.abs(x-y)/(np.abs(x)+np.abs(y)+1e-12)
    rows.append({"check":"endpoint active-rate agreement with dynamic v2.1 best model","value":float(np.nanmedian(srd)),"criterion":"report median SRD; <=0.05 expected","pass":bool(np.nanmedian(srd)<=0.05)})
    rows.append({"check":"active reference cloud overlap","value":1.0 if active_ref["cloud_overlap"] else 0.0,"criterion":"must be false","pass":bool(not active_ref["cloud_overlap"])})
    return pd.DataFrame(rows)


def make_figures(out, path_summary, ensemble, active_ref, isi_ref):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []
    figdir=out/"figures"; figdir.mkdir(exist_ok=True)
    files=[]
    ss=path_summary[path_summary.subset.eq("core_secure_pairs")].copy()
    families=[x for x in ["drive_early","coupled","drive_late"] if x in set(ss.path_family)]
    # Primary ISI stage markers.
    plt.figure(figsize=(8.5,5.2))
    offs={"wt_exit_p_isi":-0.17,"balance_p_isi":0.0,"sca3_entry_p_isi":0.17}
    labs={"wt_exit_p_isi":"WT-exit","balance_p_isi":"balance","sca3_entry_p_isi":"SCA3-entry"}
    x=np.arange(len(families))
    for met,off in offs.items():
        z=ss[ss.metric.eq(met)].set_index("path_family").loc[families]
        y=z["median"].to_numpy(); err=np.vstack([y-z.q25.to_numpy(),z.q75.to_numpy()-y])
        plt.errorbar(x+off,y,yerr=err,fmt="o-",capsize=4,label=labs[met])
    plt.xticks(x,[f.replace("_"," ") for f in families]);plt.ylim(0,1);plt.ylabel("path progress p");plt.xlabel("path family")
    plt.title("Primary ISI staging after v1.1 correction");plt.legend();plt.grid(True,alpha=.3);plt.tight_layout()
    p=figdir/"primary_ISI_stage_markers.png";plt.savefig(p,dpi=180);plt.close();files.append(p)
    # Active vs ISI balance at path-family level.
    plt.figure(figsize=(6.5,5.2))
    for fam in families:
        z=ss[ss.path_family.eq(fam)].set_index("metric")
        xa=float(z.loc["balance_p_isi","median"]); ya=float(z.loc["balance_p_active","median"])
        plt.scatter([xa],[ya],s=60);plt.annotate(fam.replace("_"," "),(xa,ya),xytext=(5,5),textcoords="offset points")
    plt.plot([0,1],[0,1],linestyle="--");plt.xlim(0,1);plt.ylim(0,1);plt.xlabel("ISI-based balance p");plt.ylabel("active-rate-based balance p")
    plt.title("Independent projection check");plt.grid(True,alpha=.3);plt.tight_layout()
    p=figdir/"active_vs_ISI_balance.png";plt.savefig(p,dpi=180);plt.close();files.append(p)
    return files


def validate(cfg, config_path):
    info=resolve_input(cfg,config_path)
    try:
        old=load_csvs(info["root"])
        frozen=pd.read_csv(Path(config_path).resolve().parent.parent/cfg["frozen"]["q75_reference_cells"])
        if len(old["transition_paths.csv"]) != 121524:
            raise ValueError("unexpected v1.0 transition path row count")
        if old["transition_pair_scenarios.csv"].scenario_id.nunique()!=988:
            raise ValueError("unexpected v1.0 scenario count")
        if frozen.cell_id.nunique()!=18 or int(frozen.core_q75_secure.sum())!=12:
            raise ValueError("unexpected frozen reference cohort")
        isi_ref,active_ref,active_cells=build_refs(old,frozen,cfg)
        return {"version":"1.1.0","input_kind":info["kind"],"input_path":str(info["path"]),"path_rows":len(old["transition_paths.csv"]),"support_scenarios":old["transition_pair_scenarios.csv"].scenario_id.nunique(),
            "biological_pairs":old["transition_pair_scenarios.csv"].biological_pair_key.nunique(),"secure_cells":int(frozen.core_q75_secure.sum()),
            "isi_wt_exit_A":isi_ref["wt_exit_A_threshold"],"isi_sca3_entry_A":isi_ref["sca3_entry_A_threshold"],
            "active_wt_exit_A":active_ref["wt_exit_A_threshold"],"active_sca3_entry_A":active_ref["sca3_entry_A_threshold"]}
    finally:
        cleanup_input(info)


def run(cfg,config_path):
    info=resolve_input(cfg,config_path)
    try:
        old=load_csvs(info["root"])
        package_root=Path(config_path).resolve().parent.parent
        frozen=pd.read_csv(package_root/cfg["frozen"]["q75_reference_cells"])
        out=(package_root/cfg["output"]["dir"]).resolve();out.mkdir(parents=True,exist_ok=True)
        isi_ref,active_ref,active_cells=build_refs(old,frozen,cfg)
        paths=reproject_paths(old["transition_paths.csv"],isi_ref,active_ref)
        markers=marker_table(paths,old["transition_pair_scenarios.csv"],isi_ref,active_ref,cfg)
        pair,overall=summarize(markers,frozen)
        paircurve,ens=summarize_curves(paths)
        refs,trans=reference_tables(isi_ref,active_ref)
        audit=correction_audit(old,paths,active_cells,active_ref)
        active_ref_cells=active_ref["cells"].copy()
        active_ref_cells.to_csv(out/"active_rate_reference_cells.csv",index=False)
        refs.to_csv(out/"transition_projection_reference_v1_1.csv",index=False)
        trans.to_csv(out/"transition_projection_transform_v1_1.csv",index=False)
        paths.to_csv(out/"transition_paths_reprojected.csv",index=False)
        markers.to_csv(out/"scenario_stage_markers_v1_1.csv",index=False)
        pair.to_csv(out/"biological_pair_stage_summary_v1_1.csv",index=False)
        overall.to_csv(out/"path_family_stage_summary_v1_1.csv",index=False)
        paircurve.to_csv(out/"biological_pair_curve_summary_v1_1.csv",index=False)
        ens.to_csv(out/"ensemble_curve_summary_v1_1.csv",index=False)
        audit.to_csv(out/"CORRECTION_AUDIT.csv",index=False)
        # Primary stage-only compact table.
        primary=overall[(overall.subset.eq("core_secure_pairs")) & overall.metric.isin(["wt_exit_p_isi","balance_p_isi","sca3_entry_p_isi"])].copy()
        primary.to_csv(out/"PRIMARY_ISI_STAGING.csv",index=False)
        # Occupancy summaries at primary stages from pair-level values.
        orows=[]
        for fam,g in pair.groupby("path_family"):
            for subset,x in [("all_pairs",g),("core_secure_pairs",g[g.both_core_secure])]:
                for stage in ["wt_exit","balance","sca3_entry"]:
                    c=f"occupancy_at_{stage}_isi_weighted_median";v=pd.to_numeric(x[c],errors="coerce").dropna().to_numpy(float)
                    orows.append({"path_family":fam,"subset":subset,"stage":stage,"n_pairs":len(v),"median_occupancy":float(np.median(v)) if len(v) else np.nan,"q25":float(np.quantile(v,.25)) if len(v) else np.nan,"q75":float(np.quantile(v,.75)) if len(v) else np.nan})
        pd.DataFrame(orows).to_csv(out/"occupancy_at_primary_ISI_stages.csv",index=False)
        figures=make_figures(out,overall,ens,active_ref,isi_ref)
        summary={"version":"1.1.0","analysis":"WT_to_SCA3_transition_reprojection","new_HR_simulations":0,"input_path":str(info["path"]),"input_kind":info["kind"],
            "path_rows":int(len(paths)),"support_scenarios":int(old["transition_pair_scenarios.csv"].scenario_id.nunique()),"biological_pairs":int(old["transition_pair_scenarios.csv"].biological_pair_key.nunique()),
            "primary_projection":"ISI (frozen v1.0 consistent reference)","secondary_projection":"experimental q=.75 active rate (v2.1)",
            "active_wt_exit_A":active_ref["wt_exit_A_threshold"],"active_sca3_entry_A":active_ref["sca3_entry_A_threshold"],"isi_wt_exit_A":isi_ref["wt_exit_A_threshold"],"isi_sca3_entry_A":isi_ref["sca3_entry_A_threshold"],
            "audit_pass":bool(audit["pass"].all()),"figure_files":[str(p.relative_to(out)) for p in figures]}
        (out/"RUN_SUMMARY.json").write_text(json.dumps(summary,indent=2,sort_keys=True),encoding="utf-8")
        return summary
    finally:
        cleanup_input(info)
