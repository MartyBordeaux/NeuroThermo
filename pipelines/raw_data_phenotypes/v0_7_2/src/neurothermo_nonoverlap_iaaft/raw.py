from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

from .metrics import uniform_resample


def resolve_source_path(source_path: str, group: str, raw_root: Path) -> Path:
    original = Path(str(source_path)).expanduser()
    candidates = [original]
    basename = original.name
    candidates.extend([
        raw_root / str(group) / basename,
        raw_root / str(group).upper() / basename,
        raw_root / str(group).lower() / basename,
        raw_root / basename,
    ])
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "Raw trace not found for {}. Checked: {}".format(
            source_path, "; ".join(str(x) for x in candidates)
        )
    )


def _unit_scale(unit: str) -> float:
    return {"v": 1000.0, "mv": 1.0, "uv": 0.001, "µv": 0.001}.get(str(unit).strip().lower(), 1.0)


def _load_abf(path: Path, sweep_index: int) -> Tuple[np.ndarray, np.ndarray]:
    try:
        import pyabf
    except ImportError as exc:
        raise RuntimeError("ABF input requires pyabf. Install the package requirements first.") from exc
    abf = pyabf.ABF(str(path))
    if int(sweep_index) not in set(int(x) for x in abf.sweepList):
        raise IndexError("Sweep {} is absent from {}".format(sweep_index, path))
    abf.setSweep(int(sweep_index), channel=0)
    time_s = np.asarray(abf.sweepX, dtype=float).copy()
    voltage_mV = np.asarray(abf.sweepY, dtype=float).copy()
    voltage_mV *= _unit_scale(getattr(abf, "sweepUnitsY", "mV"))
    return time_s, voltage_mV


def _load_npz(path: Path, sweep_index: int) -> Tuple[np.ndarray, np.ndarray]:
    data = np.load(path, allow_pickle=False)
    time_s = np.asarray(data["time_s"], dtype=float)
    voltage_mV = np.asarray(data["voltage_mV"], dtype=float)
    if voltage_mV.ndim == 2:
        voltage_mV = voltage_mV[int(sweep_index)]
    if time_s.ndim == 2:
        time_s = time_s[int(sweep_index)]
    return time_s, voltage_mV


def _load_csv(path: Path, sweep_index: int) -> Tuple[np.ndarray, np.ndarray]:
    frame = pd.read_csv(path)
    if "sweep_index" in frame.columns:
        frame = frame[frame.sweep_index.astype(int) == int(sweep_index)]
    if not {"time_s", "voltage_mV"}.issubset(frame.columns) or frame.empty:
        raise ValueError("CSV trace requires time_s and voltage_mV columns: {}".format(path))
    return frame.time_s.to_numpy(float), frame.voltage_mV.to_numpy(float)


def stationary_trace(row: pd.Series, raw_root: Path, config) -> Tuple[np.ndarray, dict]:
    path = resolve_source_path(str(row.source_path), str(row.group), raw_root)
    suffix = path.suffix.lower()
    if suffix == ".abf":
        time_s, voltage_mV = _load_abf(path, int(row.sweep_index))
    elif suffix == ".npz":
        time_s, voltage_mV = _load_npz(path, int(row.sweep_index))
    elif suffix == ".csv":
        time_s, voltage_mV = _load_csv(path, int(row.sweep_index))
    else:
        raise ValueError("Unsupported raw trace format: {}".format(path))

    finite = np.isfinite(time_s) & np.isfinite(voltage_mV)
    time_s = time_s[finite]
    voltage_mV = voltage_mV[finite]
    stationary_start = float(row.stim_start_s) + float(config["surrogates"]["stationary_discard_ms"]) / 1000.0
    stationary_end = float(row.stim_end_s)
    mask = (time_s >= stationary_start) & (time_s <= stationary_end)
    minimum = int(config["surrogates"]["minimum_stationary_samples"])
    if int(mask.sum()) < minimum:
        raise ValueError(
            "Short stationary interval for {} at {} pA: {} samples".format(
                row.cell_id, row.current_pA, int(mask.sum())
            )
        )
    target_dt = float(config["surrogates"]["resample_dt_ms"]) / 1000.0
    _, values = uniform_resample(time_s[mask], voltage_mV[mask], target_dt)
    values = values - float(np.mean(values))
    metadata = {
        "resolved_source_path": str(path),
        "raw_stationary_samples": int(mask.sum()),
        "resampled_stationary_samples": int(len(values)),
        "stationary_start_s": stationary_start,
        "stationary_end_s": stationary_end,
        "target_dt_ms": target_dt * 1000.0,
    }
    return values, metadata
