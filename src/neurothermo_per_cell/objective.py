from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np
import pandas as pd

from .models import simulate


def predict(model: str, theta: Sequence[float], df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    currents = df["current_density_pA_per_pF"].to_numpy(float)
    rows = simulate(model, theta, currents, float(cfg["dt_ms"]),
                    float(cfg["duration_ms"]), float(cfg["sustained_start_ms"]))
    out = df.copy().reset_index(drop=True)
    for key in rows[0]:
        out[key] = [r[key] for r in rows]
    out["rate_residual_hz"] = out["pred_sustained_rate_hz"] - out["sustained_rate_hz"]
    return out


def metrics(pred: pd.DataFrame, cfg: dict) -> dict:
    obs = pred["sustained_rate_hz"].to_numpy(float)
    hat = pred["pred_sustained_rate_hz"].to_numpy(float)
    rate_scale = max(float(np.nanmax(obs)), float(cfg.get("rate_scale_floor_hz", 20.0)))
    rate_loss = float(np.mean(((hat - obs) / rate_scale) ** 2))
    obs_active, pred_active = obs > 0.0, hat > 0.0
    recruitment_loss = float(np.mean(obs_active != pred_active))
    obs_lat = pred["first_spike_latency_ms"].to_numpy(float)
    hat_lat = pred["pred_first_spike_latency_ms"].to_numpy(float)
    both = np.isfinite(obs_lat) & np.isfinite(hat_lat)
    if both.any():
        latency_rmse = float(np.sqrt(np.mean((hat_lat[both] - obs_lat[both]) ** 2)))
        latency_loss = float(np.mean(((hat_lat[both] - obs_lat[both]) / 200.0) ** 2))
    else:
        latency_rmse, latency_loss = np.nan, 1.0 if obs_active.any() else 0.0
    missing_latency = float(np.mean(np.isfinite(obs_lat) != np.isfinite(hat_lat)))
    weights = cfg["loss_weights"]
    total = (float(weights["rate"]) * rate_loss +
             float(weights["recruitment"]) * recruitment_loss +
             float(weights["latency"]) * (latency_loss + missing_latency))
    return {
        "loss_total": float(total),
        "rate_loss_component": rate_loss,
        "recruitment_loss_component": recruitment_loss,
        "latency_loss_component": latency_loss,
        "missing_latency_fraction": missing_latency,
        "fi_rmse_hz": float(np.sqrt(np.mean((hat - obs) ** 2))),
        "mean_rate_error_hz": float(np.mean(hat - obs)),
        "first_spike_latency_rmse_ms": latency_rmse,
        "balanced_recruitment_accuracy": _balanced_accuracy(obs_active, pred_active),
    }


def loss(model: str, theta: Sequence[float], df: pd.DataFrame, cfg: dict) -> float:
    try:
        value = metrics(predict(model, theta, df, cfg), cfg)["loss_total"]
        return value if np.isfinite(value) else 1e6
    except (FloatingPointError, OverflowError, ValueError):
        return 1e6


def _balanced_accuracy(y: np.ndarray, p: np.ndarray) -> float:
    scores = []
    for value in (False, True):
        mask = y == value
        if mask.any():
            scores.append(float(np.mean(p[mask] == value)))
    return float(np.mean(scores)) if scores else np.nan
