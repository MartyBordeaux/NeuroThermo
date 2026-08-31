#!/usr/bin/env python3
"""Locate QC candidate peaks that lie numerically close to the frozen find_peaks thresholds."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "code" / "pipelines" / "NeuroThermo_stage1_qc_fixed"
CFG_PATH = PIPE / "config.json"
REF = ROOT / "data" / "calibration" / "candidate_events_with_predictions.csv"

spec = importlib.util.spec_from_file_location("qc_fixed", PIPE / "spike_qc_calibrated.py")
qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(qc)


def main() -> int:
    cfg = json.loads(CFG_PATH.read_text())
    prom_thr = float(cfg["candidate_prominence_mV"])
    height_thr = float(cfg["candidate_height_mV"])
    eps = 1e-9

    sweeps = []
    for group in ("WT", "SCA3"):
        root = ROOT / "data" / "raw" / group
        for path in qc.find_cc(root):
            sweeps.extend(qc.load_sweeps(group, path, cfg))

    rows = []
    relaxed_only = []
    for s in sweeps:
        t = s["time_ms"]
        v = s["voltage_mV"]
        dt = float(np.median(np.diff(t)))
        distance = max(1, int(round(float(cfg["minimum_peak_distance_ms"]) / dt)))
        exact, _ = find_peaks(v, prominence=prom_thr, height=height_thr, distance=distance)
        relaxed, _ = find_peaks(v, prominence=prom_thr-eps, height=height_thr-eps, distance=distance)
        exact_set = set(map(int, exact))
        rel_prom = peak_prominences(v, relaxed)[0] if len(relaxed) else np.array([])
        for p, pr in zip(relaxed, rel_prom):
            p = int(p)
            row = {
                "group": s["group"], "cell_id": s["cell_id"], "sweep_index": int(s["sweep_index"]),
                "current_pA": float(s["current_pA"]), "peak_index": p, "time_ms": float(t[p]),
                "peak_voltage_mV": float(v[p]), "prominence_mV": float(pr),
                "height_margin_mV": float(v[p]-height_thr), "prominence_margin_mV": float(pr-prom_thr),
                "exact_selected": p in exact_set,
            }
            rows.append(row)
            if p not in exact_set:
                relaxed_only.append(row)

    audit = pd.DataFrame(rows)
    exact_audit = audit[audit.exact_selected].copy()
    exact_count = int(exact_audit.shape[0])
    print(f"exact_candidate_count={exact_count}")
    print(f"relaxed_only_count={len(relaxed_only)}")

    exact_audit["boundary_margin_mV"] = exact_audit[["height_margin_mV", "prominence_margin_mV"]].min(axis=1)
    print("\nClosest exact candidates to either threshold:")
    print(exact_audit.sort_values("boundary_margin_mV").head(20).to_string(index=False))
    if relaxed_only:
        print("\nRelaxed-only candidates within 1e-9 mV:")
        print(pd.DataFrame(relaxed_only).to_string(index=False))

    if REF.is_file():
        ref = pd.read_csv(REF)
        key = ["group", "cell_id", "sweep_index", "peak_index"]
        frozen = ref[key + ["time_ms", "peak_voltage_mV", "prominence_mV"]].copy()
        got = exact_audit[key + ["time_ms", "peak_voltage_mV", "prominence_mV", "height_margin_mV", "prominence_margin_mV"]].copy()
        merged = frozen.merge(got, on=key, how="outer", indicator=True, suffixes=("_frozen", "_now"))
        missing = merged[merged._merge == "left_only"]
        extra = merged[merged._merge == "right_only"]
        print(f"\nfrozen_rows={len(frozen)} missing_from_exact={len(missing)} extra_vs_frozen={len(extra)}")
        if len(missing):
            print("Missing frozen identities:")
            print(missing.to_string(index=False))
        if len(extra):
            print("Extra current identities:")
            print(extra.to_string(index=False))

    out = ROOT / "results" / "recomputed" / "qc_peak_boundary_audit.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out, index=False)
    print(f"wrote={out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
