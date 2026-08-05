from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Optional

import pandas as pd


REQUIRED = {
    "record_id", "cell_id", "group", "sweep_index", "current_pA",
    "capacitance_window_ms", "capacitance_pF", "current_density_pA_per_pF",
    "sustained_rate_hz", "first_spike_latency_ms",
}


def load_observations(path: Optional[str] = None, window_ms: int = 20) -> pd.DataFrame:
    if path:
        p = Path(path).expanduser().resolve()
        if p.is_dir():
            candidates = list(p.rglob("cell_current_observations.csv"))
            if not candidates:
                candidates = list(p.rglob("*observations*.csv"))
            if not candidates:
                raise FileNotFoundError("no observation CSV found below %s" % p)
            p = candidates[0]
        df = pd.read_csv(p)
    else:
        p = resources.files("neurothermo_per_cell").joinpath("data/frozen_v2_w20_observations.csv")
        df = pd.read_csv(p)
    missing = REQUIRED.difference(df.columns)
    if missing:
        raise ValueError("missing observation columns: %s" % sorted(missing))
    df = df[df["capacitance_window_ms"].astype(int) == int(window_ms)].copy()
    df = df.drop_duplicates(["group", "cell_id", "sweep_index"])
    df = df.sort_values(["group", "cell_id", "sweep_index"]).reset_index(drop=True)
    if df.empty:
        raise ValueError("no observations for capacitance window %d ms" % window_ms)
    return df


def validate_grid(df: pd.DataFrame) -> dict:
    counts = df.groupby(["group", "cell_id"]).size()
    if (counts < 8).any():
        raise ValueError("each cell must contain at least eight current levels")
    duplicate = df.duplicated(["group", "cell_id", "sweep_index"]).any()
    if duplicate:
        raise ValueError("duplicate cell/sweep rows")
    return {
        "rows": int(len(df)),
        "cells": int(df["cell_id"].nunique()),
        "groups": {str(k): int(v) for k, v in df.groupby("group")["cell_id"].nunique().items()},
        "levels_per_cell_min": int(counts.min()),
        "levels_per_cell_max": int(counts.max()),
    }
