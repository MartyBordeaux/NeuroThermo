from __future__ import annotations
import numpy as np
from typing import Optional
from numba import njit

@njit(cache=True)
def _rhs(x, y, z, J, b, r, s, kappa_I, a, c, d, x_R):
    dx = y - a*x*x*x + b*x*x - z + kappa_I*J
    dy = c - d*x*x - y
    dz = r * (s*(x - x_R) - z)
    return dx, dy, dz

@njit(cache=True)
def _rk4_step(x, y, z, h, J, b, r, s, kappa_I, a, c, d, x_R):
    k1x, k1y, k1z = _rhs(x, y, z, J, b, r, s, kappa_I, a, c, d, x_R)
    k2x, k2y, k2z = _rhs(x+0.5*h*k1x, y+0.5*h*k1y, z+0.5*h*k1z, J, b, r, s, kappa_I, a, c, d, x_R)
    k3x, k3y, k3z = _rhs(x+0.5*h*k2x, y+0.5*h*k2y, z+0.5*h*k2z, J, b, r, s, kappa_I, a, c, d, x_R)
    k4x, k4y, k4z = _rhs(x+h*k3x, y+h*k3y, z+h*k3z, J, b, r, s, kappa_I, a, c, d, x_R)
    return (
        x + h*(k1x + 2*k2x + 2*k3x + k4x)/6.0,
        y + h*(k1y + 2*k2y + 2*k3y + k4y)/6.0,
        z + h*(k1z + 2*k2z + 2*k3z + k4z)/6.0,
    )

@njit(cache=True)
def _equilibrate(x, y, z, n, h, b, r, s, kappa_I, a, c, d, x_R):
    for _ in range(n):
        x, y, z = _rk4_step(x, y, z, h, 0.0, b, r, s, kappa_I, a, c, d, x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)):
            return x, y, z, False
        if abs(x) > 1e6 or abs(y) > 1e6 or abs(z) > 1e6:
            return x, y, z, False
    return x, y, z, True

@njit(cache=True)
def _simulate_spikes_core(
    J, stimulus_duration_ms, observation_end_ms, pre_ms, dt_ms, model_time_scale_ms,
    b, r, s, kappa_I, a, c, d, x_R, x0, y0, z0,
    threshold, refractory_ms, max_spikes,
):
    h = dt_ms / model_time_scale_ms
    n_pre = int(np.ceil(pre_ms / dt_ms))
    n_total = int(np.ceil(observation_end_ms / dt_ms))
    x, y, z, ok = _equilibrate(x0, y0, z0, n_pre, h, b, r, s, kappa_I, a, c, d, x_R)
    spikes = np.empty(max_spikes, dtype=np.float64)
    count = 0
    if not ok:
        return spikes[:0], False
    last_spike = -1e30
    x_prev = x
    for i in range(1, n_total + 1):
        t_prev = (i - 1) * dt_ms
        drive = J if t_prev < stimulus_duration_ms else 0.0
        x, y, z = _rk4_step(x, y, z, h, drive, b, r, s, kappa_I, a, c, d, x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)):
            return spikes[:count], False
        if abs(x) > 1e6 or abs(y) > 1e6 or abs(z) > 1e6:
            return spikes[:count], False
        t_ms = i * dt_ms
        if x_prev < threshold and x >= threshold and (t_ms - last_spike) >= refractory_ms:
            denom = x - x_prev
            frac = (threshold - x_prev) / denom if denom != 0.0 else 0.0
            t_cross = (i - 1 + frac) * dt_ms
            if t_cross <= observation_end_ms + 1e-9:
                if count < max_spikes:
                    spikes[count] = t_cross
                    count += 1
                else:
                    return spikes[:count], False
                last_spike = t_cross
        x_prev = x
    return spikes[:count], True

def simulate_spikes(theta, sweep, cfg, dt_ms: float, observation_end_ms: Optional[float] = None):
    m = cfg['model']
    stimulus_duration = float(sweep['stimulus_duration_ms'])
    if observation_end_ms is None:
        observation_end_ms = float(sweep['fit_end_ms'])
    observation_end_ms = float(observation_end_ms)
    max_spikes = max(1000, int(observation_end_ms / max(float(m['model_refractory_ms']), dt_ms)) + 10)
    spikes, ok = _simulate_spikes_core(
        float(sweep['J']), stimulus_duration, observation_end_ms,
        float(m['pre_ms']), float(dt_ms), float(m['model_time_scale_ms']),
        float(theta['b']), float(theta['r']), float(theta['s']), float(theta['kappa_I']),
        float(m['a']), float(m['c']), float(m['d']), float(m['x_R']),
        float(m['x0']), float(m['y0']), float(m['z0']),
        float(m['model_spike_threshold']), float(m['model_refractory_ms']), int(max_spikes),
    )
    return np.asarray(spikes, dtype=float), bool(ok)
