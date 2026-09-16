from __future__ import annotations
import numpy as np
import pandas as pd


def _mad_scale(x, eps=1e-12):
    a = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    a = a[np.isfinite(a)]
    if len(a) == 0:
        raise ValueError("cannot fit transform to empty coordinate")
    med = float(np.median(a))
    mad = float(np.median(np.abs(a - med)))
    if mad < eps:
        q25, q75 = np.quantile(a, [0.25, 0.75])
        mad = float((q75 - q25) / 1.349) if q75 > q25 else 1.0
    return med, mad


def fit_reference(cells, xcol, ycol, secure_col="core_q75_secure", wt_q=1.0, sca_q=0.0):
    d = cells[cells[secure_col].fillna(False).astype(bool)].copy()
    if (d.group == "WT").sum() < 2 or (d.group == "SCA3").sum() < 2:
        raise ValueError("insufficient secure endpoint cells")
    d["log10_x"] = np.log10(pd.to_numeric(d[xcol], errors="coerce"))
    d["log10_y"] = np.log10(pd.to_numeric(d[ycol], errors="coerce"))
    cx, sx = _mad_scale(d["log10_x"])
    cy, sy = _mad_scale(d["log10_y"])
    d["z0"] = (d["log10_x"] - cx) / sx
    d["z1"] = (d["log10_y"] - cy) / sy
    z = d[["z0", "z1"]].to_numpy(float)
    cwt = d.loc[d.group.eq("WT"), ["z0", "z1"]].mean().to_numpy(float)
    csc = d.loc[d.group.eq("SCA3"), ["z0", "z1"]].mean().to_numpy(float)
    delta = csc - cwt
    den = float(delta @ delta)
    if den <= 0:
        raise ValueError("degenerate endpoint centroids")
    A = ((z - cwt) @ delta) / den
    d["A"] = A
    wtA = d.loc[d.group.eq("WT"), "A"].to_numpy(float)
    scA = d.loc[d.group.eq("SCA3"), "A"].to_numpy(float)
    foot = cwt[None, :] + A[:, None] * delta[None, :]
    orth = np.linalg.norm(z - foot, axis=1)
    return {
        "cells": d,
        "center": np.array([cx, cy], float),
        "scale": np.array([sx, sy], float),
        "cwt": cwt,
        "csc": csc,
        "delta": delta,
        "den": den,
        "centroid_distance": float(np.sqrt(den)),
        "wt_exit_A_threshold": float(np.quantile(wtA, wt_q)),
        "sca3_entry_A_threshold": float(np.quantile(scA, sca_q)),
        "cloud_overlap": bool(np.quantile(wtA, wt_q) >= np.quantile(scA, sca_q)),
        "corridor_radius_q90": float(np.quantile(orth, 0.90)),
        "xcol": xcol,
        "ycol": ycol,
    }


def reference_from_v1_tables(ref_row, transform_rows):
    t = transform_rows.set_index("coordinate")
    coords = ["log10_rheobase", "log10_isi"]
    return {
        "center": np.array([float(t.loc[c, "center"]) for c in coords]),
        "scale": np.array([float(t.loc[c, "scale"]) for c in coords]),
        "cwt": np.array([float(ref_row.wt_centroid_0), float(ref_row.wt_centroid_1)]),
        "csc": np.array([float(ref_row.sca3_centroid_0), float(ref_row.sca3_centroid_1)]),
        "delta": np.array([float(ref_row.sca3_centroid_0-ref_row.wt_centroid_0), float(ref_row.sca3_centroid_1-ref_row.wt_centroid_1)]),
        "den": float(ref_row.centroid_distance) ** 2,
        "centroid_distance": float(ref_row.centroid_distance),
        "wt_exit_A_threshold": float(ref_row.wt_exit_A_threshold),
        "sca3_entry_A_threshold": float(ref_row.sca3_entry_A_threshold),
        "cloud_overlap": bool(ref_row.cloud_overlap),
        "corridor_radius_q90": float(ref_row.corridor_radius_q90),
        "xcol": "rheobase_J",
        "ycol": "mean_isi_ms",
    }


def project_native(df, ref, xcol, ycol):
    x = pd.to_numeric(df[xcol], errors="coerce").to_numpy(float)
    y = pd.to_numeric(df[ycol], errors="coerce").to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z0 = (np.log10(x) - ref["center"][0]) / ref["scale"][0]
        z1 = (np.log10(y) - ref["center"][1]) / ref["scale"][1]
    arr = np.column_stack([z0, z1])
    good = np.all(np.isfinite(arr), axis=1)
    A = np.full(len(df), np.nan)
    O = np.full(len(df), np.nan)
    if good.any():
        aa = ((arr[good] - ref["cwt"]) @ ref["delta"]) / ref["den"]
        foot = ref["cwt"][None, :] + aa[:, None] * ref["delta"][None, :]
        A[good] = aa
        O[good] = np.linalg.norm(arr[good] - foot, axis=1)
    return A, O


def persistent_crossing(x, y, thr, persistence=2):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    for i in range(len(x)):
        if not np.isfinite(y[i]) or y[i] <= thr:
            continue
        j = min(len(x), i + persistence)
        if np.all(np.isfinite(y[i:j])) and np.all(y[i:j] > thr):
            if i == 0:
                return float(x[0])
            if np.isfinite(y[i-1]) and y[i-1] <= thr:
                if y[i] == y[i-1]:
                    return float(x[i])
                return float(x[i-1] + (thr-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]))
            return float(x[i])
    return np.nan


def weighted_quantile(values, weights, q):
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[m], w[m]
    if len(v) == 0:
        return np.nan
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w) / np.sum(w)
    return float(np.interp(q, c, v))


def interp_at(x, y, x0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2 or not np.isfinite(x0):
        return np.nan
    x, y = x[m], y[m]
    o = np.argsort(x); x, y = x[o], y[o]
    if x0 < x[0] or x0 > x[-1]:
        return np.nan
    return float(np.interp(x0, x, y))
