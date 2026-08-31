#!/usr/bin/env python3
"""Diagnose platform sensitivity of Stage-1 candidate peak extraction.

This intentionally stops before classifier fitting. It compares candidate identities and
raw peak descriptors produced from the publication ABFs against the frozen historical
candidate table.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "code" / "pipelines" / "NeuroThermo_stage1_qc_fixed"
RAW = ROOT / "data" / "raw"
REF = ROOT / "data" / "calibration" / "candidate_events_with_predictions.csv"
MOD_PATH = PIPE / "spike_qc_calibrated.py"


def load_module():
    spec = importlib.util.spec_from_file_location("stage1_qc_fixed", MOD_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {MOD_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def key_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["current_key_pA"] = pd.to_numeric(out["current_pA"]).round(6)
    out["time_key_ms"] = pd.to_numeric(out["time_ms"]).round(9)
    out["peak_index"] = pd.to_numeric(out["peak_index"]).astype(int)
    return out


def key_tuple(row):
    return (
        str(row.group), str(row.cell_id), int(row.sweep_index),
        float(row.current_key_pA), int(row.peak_index), float(row.time_key_ms),
    )


def main() -> int:
    mod = load_module()
    cfg = json.loads((PIPE / "config.json").read_text(encoding="utf-8"))
    sweeps = []
    for group in ("WT", "SCA3"):
        files = mod.find_cc(RAW / group)
        for path in files:
            sweeps.extend(mod.load_sweeps(group, path, cfg))
    tables = [mod.candidate_features(s, cfg) for s in sweeps]
    current = pd.concat([x for x in tables if len(x)], ignore_index=True)
    ref = pd.read_csv(REF)

    cur = key_frame(current)
    old = key_frame(ref)
    cur_map = {key_tuple(r): r for r in cur.itertuples(index=False)}
    old_map = {key_tuple(r): r for r in old.itertuples(index=False)}
    missing = sorted(set(old_map) - set(cur_map))
    extra = sorted(set(cur_map) - set(old_map))

    print(f"CURRENT_CANDIDATES={len(cur)}")
    print(f"FROZEN_CANDIDATES={len(old)}")
    print(f"MISSING_FROM_RECOMPUTE={len(missing)}")
    print(f"EXTRA_IN_RECOMPUTE={len(extra)}")
    print(f"CONFIG candidate_prominence_mV={cfg['candidate_prominence_mV']} candidate_height_mV={cfg['candidate_height_mV']} minimum_peak_distance_ms={cfg['minimum_peak_distance_ms']}")

    cols = [
        "group", "cell_id", "sweep_index", "current_pA", "peak_index", "time_ms",
        "peak_voltage_mV", "prominence_mV", "local_amplitude_mV", "max_dvdt_mV_per_ms",
        "abs_min_dvdt_mV_per_ms", "half_width_ms", "ahp_depth_mV",
    ]
    for label, keys, mp in (("MISSING", missing, old_map), ("EXTRA", extra, cur_map)):
        for i, k in enumerate(keys[:20], 1):
            r = mp[k]
            vals = {c: getattr(r, c) for c in cols if hasattr(r, c)}
            print(label, i, json.dumps(vals, default=str, sort_keys=True))

    # For identities present in both, quantify descriptor drift before classifier fitting.
    common = sorted(set(cur_map) & set(old_map))
    numeric = [
        "time_ms", "peak_voltage_mV", "prominence_mV", "local_amplitude_mV",
        "max_dvdt_mV_per_ms", "abs_min_dvdt_mV_per_ms", "half_width_ms", "ahp_depth_mV",
    ]
    for c in numeric:
        diffs = np.array([abs(float(getattr(cur_map[k], c)) - float(getattr(old_map[k], c))) for k in common], dtype=float)
        print(f"MAX_ABS_DIFF {c} {np.nanmax(diffs):.17g}")
        print(f"NONZERO_DIFF {c} {int(np.count_nonzero(diffs > 0))}")

    # Diagnostic only: do not fail merely because a platform-sensitive candidate was found.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
