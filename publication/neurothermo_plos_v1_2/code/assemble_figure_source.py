#!/usr/bin/env python3
"""Rebuild publication-facing Fig.1--3 source tables from upstream results.

Committed data/figure_source files are frozen references. This script writes a
clean recomputation under results/recomputed/figure_source and validates
scientific equivalence without overwriting the frozen release data.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pyabf

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "data" / "figure_source"
OUT = ROOT / "results" / "recomputed" / "figure_source"
OUT.mkdir(parents=True, exist_ok=True)


def truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1"])


def read(name: str, directory: Path) -> pd.DataFrame:
    return pd.read_csv(directory / name)


def date_label_from_abf(group: str, cell_id: str) -> str:
    nn = cell_id.split("_")[-1]
    p = ROOT / "data" / "raw" / group / f"cc_{nn}.abf"
    abf = pyabf.ABF(str(p), loadData=False)
    if abf.abfDateTime is None:
        raise RuntimeError(f"ABF datetime missing: {p}")
    return f"{group}_{abf.abfDateTime:%y%m%d}"


def build_fig1() -> pd.DataFrame:
    src = pd.read_csv(ROOT / "data" / "endpoint_ensemble_v1_0_results" / "endpoint_cells_full_observable.csv")
    cols = list(pd.read_csv(FROZEN / "fig1_endpoint_cells.csv", nrows=0).columns)
    z = src.copy()
    # Publication capacitance-aware normalization: preserve pA while changing J=I/Cm.
    z["capacitance_pF"] = pd.to_numeric(z["capacitance_pF"]) * 2.0
    z["kappa_I"] = pd.to_numeric(z["kappa_I"]) * 2.0
    z["rheobase_J_best"] = pd.to_numeric(z["rheobase_J_best"]) / 2.0

    labels = []
    for _, row in z.iterrows():
        status = str(row.get("animal_id_status", ""))
        day = row.get("experiment_day_code")
        if status == "RECOVERED_FROM_EXACT_ABF_MATCH_AND_DAY_CODE" and pd.notna(day):
            labels.append(f"{row['group']}_DD{int(float(day)):02d}")
        else:
            labels.append(date_label_from_abf(str(row["group"]), str(row["cell_id"])))
    z["animal_id"] = labels
    return z[cols].reset_index(drop=True)


def build_fig2_curves() -> pd.DataFrame:
    endpoint = pd.read_csv(ROOT / "data" / "endpoint_ensemble_v1_0_results" / "endpoint_cells_full_observable.csv")
    core = endpoint[truthy(endpoint["core_q75_secure"])].copy()
    wt = set(core.loc[core.group.eq("WT"), "cell_id"])
    sca = set(core.loc[core.group.eq("SCA3"), "cell_id"])
    if (len(wt), len(sca)) != (8, 4):
        raise RuntimeError(f"Expected 8 WT x 4 SCA3 core-secure cells, got {len(wt)} x {len(sca)}")

    pair = pd.read_csv(ROOT / "data" / "transition_v1_1_results" / "biological_pair_curve_summary_v1_1.csv")
    q = pair[pair.wt_cell_id.isin(wt) & pair.sca_cell_id.isin(sca)].copy()
    if q.biological_pair_key.nunique() != 32:
        raise RuntimeError(f"Expected 32 core-secure biological pairs, got {q.biological_pair_key.nunique()}")

    rows=[]
    metrics = {
        "A_isi": "A_isi_v1_1_weighted_median",
        "rheobase_J": "rheobase_J_weighted_median",
        "active_rate": "active_support_rate_hz_weighted_median",
        "mean_isi": "mean_isi_ms_weighted_median",
        "occupancy": "occupancy_fraction_weighted_median",
    }
    for (fam, prog), g in q.groupby(["path_family", "path_progress"], sort=True):
        row={"path_family":fam, "path_progress":prog, "n_pairs":int(g.biological_pair_key.nunique())}
        a=pd.to_numeric(g[metrics["A_isi"]], errors="coerce").dropna()
        row.update({"A_isi_median":a.median(), "A_isi_q25":a.quantile(.25), "A_isi_q75":a.quantile(.75)})
        # Publication capacitance-aware normalization halves current density J.
        r=pd.to_numeric(g[metrics["rheobase_J"]], errors="coerce").dropna()
        row["rheobase_J_median"] = r.median()/2.0
        for outcol,key in [("active_rate_median","active_rate"),("mean_isi_median","mean_isi"),("occupancy_median","occupancy")]:
            v=pd.to_numeric(g[metrics[key]], errors="coerce").dropna()
            row[outcol]=v.median()
        rows.append(row)
    return pd.DataFrame(rows)


def direct_projection(src_path: Path, frozen_name: str, mask=None) -> pd.DataFrame:
    src = pd.read_csv(src_path)
    if mask is not None:
        src = mask(src)
    cols = list(pd.read_csv(FROZEN/frozen_name, nrows=0).columns)
    return src[cols].reset_index(drop=True)


def validate(name: str, got: pd.DataFrame) -> dict:
    ref = pd.read_csv(FROZEN/name)
    if list(got.columns) != list(ref.columns):
        raise RuntimeError(f"{name}: column mismatch")
    if len(got) != len(ref):
        raise RuntimeError(f"{name}: row mismatch {len(got)} != {len(ref)}")
    maxdiff=0.0
    for c in ref.columns:
        if pd.api.types.is_numeric_dtype(ref[c]):
            a=pd.to_numeric(got[c],errors="coerce").to_numpy(float)
            b=pd.to_numeric(ref[c],errors="coerce").to_numpy(float)
            if not np.allclose(a,b,rtol=1e-12,atol=1e-12,equal_nan=True):
                finite=np.isfinite(a)&np.isfinite(b)
                d=float(np.max(np.abs(a[finite]-b[finite]))) if finite.any() else float("nan")
                raise RuntimeError(f"{name}: numeric mismatch in {c}; max_abs_diff={d}")
            finite=np.isfinite(a)&np.isfinite(b)
            if finite.any(): maxdiff=max(maxdiff,float(np.max(np.abs(a[finite]-b[finite]))))
        else:
            a=got[c].fillna("").astype(str).to_numpy(); b=ref[c].fillna("").astype(str).to_numpy()
            if not np.array_equal(a,b):
                raise RuntimeError(f"{name}: text mismatch in {c}")
    got.to_csv(OUT/name,index=False)
    return {"file":name,"rows":len(got),"max_abs_numeric_diff":maxdiff,"status":"PASS"}


def main() -> int:
    tables={
        "fig1_endpoint_cells.csv": build_fig1(),
        "fig2_core_secure_curves.csv": build_fig2_curves(),
        "fig2_primary_isi_staging.csv": direct_projection(
            ROOT/"data/transition_v1_1_results/PRIMARY_ISI_STAGING.csv", "fig2_primary_isi_staging.csv"),
        "fig2_projection_reference.csv": direct_projection(
            ROOT/"data/transition_v1_1_results/transition_projection_reference_v1_1.csv", "fig2_projection_reference.csv"),
        "fig3_drive_surface_core_secure.csv": direct_projection(
            ROOT/"data/transition_v1_2_1_results/drive_sensitivity_surface_v1_2_1.csv", "fig3_drive_surface_core_secure.csv",
            lambda d: d[d["subset"].eq("core_secure_pairs")]),
        "fig3_coupled_component_sensitivity_core_secure.csv": direct_projection(
            ROOT/"data/transition_v1_3_results/coupled_line_component_sensitivity_v1_3.csv", "fig3_coupled_component_sensitivity_core_secure.csv",
            lambda d: d[d["subset"].eq("core_secure_pairs") & d["projection"].eq("isi")]),
        "fig3_drive_sensitivity_at_boundaries_core_secure.csv": direct_projection(
            ROOT/"data/transition_v1_2_1_results/drive_sensitivity_at_stage_boundaries_v1_2_1.csv", "fig3_drive_sensitivity_at_boundaries_core_secure.csv",
            lambda d: d[d["subset"].eq("core_secure_pairs")]),
        "fig3_interaction_at_boundaries_core_secure.csv": direct_projection(
            ROOT/"data/transition_v1_3_results/interaction_at_stage_boundaries_v1_3.csv", "fig3_interaction_at_boundaries_core_secure.csv",
            lambda d: d[d["subset"].eq("core_secure_pairs") & d["projection"].eq("isi")]),
    }
    report=[]
    for name,df in tables.items():
        rec=validate(name,df); report.append(rec); print("PASS",name,"rows",rec["rows"],"maxdiff",rec["max_abs_numeric_diff"])
    meta={
        "status":"PASS",
        "tables":report,
        "note":"fig3_combined_drive_handoff_summary.csv is an auxiliary diagnostic and is not consumed by the canonical Fig3 renderer",
        "output_dir":str(OUT.relative_to(ROOT)),
    }
    (OUT/"FIGURE_SOURCE_ASSEMBLY.json").write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
    print("FIGURE_SOURCE_ASSEMBLY_PASS")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
