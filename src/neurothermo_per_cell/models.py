from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class ModelSpec:
    name: str
    names: Tuple[str, ...]
    bounds: Tuple[Tuple[float, float], ...]


SPECS = {
    "AdEx": ModelSpec(
        "AdEx",
        ("tau_m_ms", "V_T_mV", "tau_w_ms", "b_bar_mV_per_ms", "kappa_J"),
        ((5.0, 40.0), (-55.0, -35.0), (20.0, 1000.0), (0.0, 6.0), (0.2, 15.0)),
    ),
    "Izhikevich": ModelSpec(
        "Izhikevich",
        ("a", "b", "c_mV", "d", "kappa_J"),
        ((0.002, 0.2), (0.05, 0.35), (-75.0, -40.0), (0.0, 12.0), (0.2, 20.0)),
    ),
}


def as_dict(model: str, values: Sequence[float]) -> Dict[str, float]:
    return dict(zip(SPECS[model].names, map(float, values)))


def simulate(model: str, values: Sequence[float], currents: np.ndarray,
             dt_ms: float = 0.2, duration_ms: float = 1000.0,
             sustained_start_ms: float = 10.0) -> List[dict]:
    if model == "AdEx":
        return _simulate_adex(values, currents, dt_ms, duration_ms, sustained_start_ms)
    if model == "Izhikevich":
        return _simulate_izh(values, currents, dt_ms, duration_ms, sustained_start_ms)
    raise ValueError("unknown model: %s" % model)


def _summarize(spikes: List[float], duration_ms: float, sustained_start_ms: float) -> dict:
    sustained = [t for t in spikes if t >= sustained_start_ms]
    window_s = max((duration_ms - sustained_start_ms) / 1000.0, 1e-9)
    return {
        "pred_total_spike_count": len(spikes),
        "pred_sustained_spike_count": len(sustained),
        "pred_sustained_rate_hz": len(sustained) / window_s,
        "pred_first_spike_latency_ms": spikes[0] if spikes else np.nan,
    }


def _simulate_adex(values, currents, dt, duration, sustained_start):
    tau_m, v_t, tau_w, b_bar, kappa = map(float, values)
    e_l, delta_t, v_reset, v_peak = -65.0, 2.0, -58.0, 20.0
    n = int(round(duration / dt))
    output = []
    for j in currents:
        v, w, spikes = e_l, 0.0, []
        invalid = False
        for step in range(n):
            expo = np.exp(np.clip((v - v_t) / delta_t, -50.0, 30.0))
            dv = (-(v - e_l) + delta_t * expo - w + kappa * float(j)) / tau_m
            dw = -w / tau_w
            v += dt * dv
            w += dt * dw
            if not np.isfinite(v + w):
                invalid = True
                break
            if v >= v_peak:
                spikes.append((step + 1) * dt)
                v = v_reset
                w += b_bar
            if len(spikes) > 500:
                invalid = True
                break
        output.append(_summarize([] if invalid else spikes, duration, sustained_start))
    return output


def _simulate_izh(values, currents, dt, duration, sustained_start):
    a, b, c, d, kappa = map(float, values)
    n = int(round(duration / dt))
    output = []
    for j in currents:
        v, u, spikes = -65.0, b * -65.0, []
        invalid = False
        for step in range(n):
            # Two half-steps stabilize the standard Izhikevich voltage update.
            drive = kappa * float(j)
            for _ in range(2):
                v += 0.5 * dt * (0.04 * v * v + 5.0 * v + 140.0 - u + drive)
            u += dt * a * (b * v - u)
            if not np.isfinite(v + u):
                invalid = True
                break
            if v >= 30.0:
                spikes.append((step + 1) * dt)
                v, u = c, u + d
            if len(spikes) > 500:
                invalid = True
                break
        output.append(_summarize([] if invalid else spikes, duration, sustained_start))
    return output
