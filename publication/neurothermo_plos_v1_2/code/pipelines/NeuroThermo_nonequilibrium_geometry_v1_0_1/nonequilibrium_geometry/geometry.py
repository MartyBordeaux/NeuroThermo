from __future__ import annotations

import numpy as np

from .density import values_at_samples
from .model import drift_on_grid


def safe_masses(masses, cfg):
    masses = np.asarray(masses, float)
    axes = tuple(range(1, masses.ndim))
    floor = float(cfg["density"].get("mass_floor_relative", 1e-12)) * np.max(masses, axis=axes, keepdims=True)
    masses = np.maximum(masses, floor)
    return masses / masses.sum(axis=axes, keepdims=True)


def path_derivatives(masses, p_grid):
    masses = np.asarray(masses, float)
    p_grid = np.asarray(p_grid, float)
    log_mass = np.log(masses)
    dlog = np.gradient(log_mass, p_grid, axis=0, edge_order=2 if len(p_grid) >= 3 else 1)
    dphi = -dlog
    axes = tuple(range(1, masses.ndim))
    score_mean = np.sum(masses * dlog, axis=axes)
    centered = dlog - score_mean.reshape((-1,) + (1,) * (masses.ndim - 1))
    path_fi = np.sum(masses * centered * centered, axis=axes)
    return dphi, path_fi, score_mean


def local_kl_check(masses, p_grid, path_fi):
    rows = []
    eps = 1e-300
    for index in range(len(p_grid) - 1):
        left, right = masses[index], masses[index + 1]
        kl = float(np.sum(left * (np.log(left + eps) - np.log(right + eps))))
        dp = float(p_grid[index + 1] - p_grid[index])
        prediction = 0.25 * float(path_fi[index] + path_fi[index + 1]) * dp * dp
        rel = abs(kl - prediction) / max(abs(kl), abs(prediction), 1e-15)
        rows.append((index, kl, prediction, rel))
    return rows


def _positive_autocovariance_integral(series, dt, max_lag):
    series = np.asarray(series, float)
    series = series[np.isfinite(series)]
    if len(series) < 8:
        return np.nan, np.nan, np.nan
    centered = series - series.mean()
    variance = float(np.mean(centered * centered))
    if variance <= 0:
        return variance, 0.0, 0.0
    size = 1 << (2 * len(centered) - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size)
    acov = np.fft.irfft(spectrum * spectrum.conjugate(), n=size)[: len(centered)]
    acov /= np.arange(len(centered), 0, -1)
    upper = min(len(acov), int(max_lag) + 1)
    stop = upper
    for index in range(1, upper):
        if acov[index] <= 0:
            stop = index
            break
    integral = float(dt * (0.5 * acov[0] + acov[1:stop].sum()))
    tau = integral / variance
    return variance, tau, integral


def geometry_for_path(samples_by_p, masses, edges, centers, theta, current, p_grid, cfg):
    masses = safe_masses(masses, cfg)
    dphi, path_fi, score_mean = path_derivatives(masses, p_grid)
    # These marginals inherit positivity and normalization from the same full
    # density. Applying an independent floor here would break the Fisher chain
    # rule and can create a spurious negative conditional contribution.
    masses_xy = masses.sum(axis=3)
    masses_xy /= masses_xy.sum(axis=(1, 2), keepdims=True)
    masses_z = masses.sum(axis=(1, 2))
    masses_z /= masses_z.sum(axis=1, keepdims=True)
    _, path_fi_xy, score_mean_xy = path_derivatives(masses_xy, p_grid)
    _, path_fi_z, score_mean_z = path_derivatives(masses_z, p_grid)
    X, Y, Z = np.meshgrid(*centers, indexing="ij")
    spacing = [float(axis[1] - axis[0]) for axis in centers]
    D = np.asarray(cfg["noise"]["D"], float) * float(cfg["noise"].get("multiplier", 1.0))
    support_relative = float(cfg["geometry"]["support_mass_relative"])
    sample_dt = float(cfg["stationary"]["sample_stride_ms"])
    max_lag = int(round(float(cfg["geometry"]["friction_max_lag_ms"]) / sample_dt))
    rows = []
    for index, p in enumerate(p_grid):
        mass = masses[index]
        log_density = np.log(mass)
        gradients = np.gradient(log_density, *spacing, edge_order=2)
        drift = drift_on_grid(X, Y, Z, current[index], theta[index], cfg)
        circulation = tuple(drift[axis] - D[axis] * gradients[axis] for axis in range(3))
        density = mass / np.prod(spacing)
        probability_current = tuple(density * value for value in circulation)
        divergence = sum(
            np.gradient(probability_current[axis], spacing[axis], axis=axis, edge_order=2)
            for axis in range(3)
        )
        grad_power_grid = sum(D[axis] * gradients[axis] ** 2 for axis in range(3))
        epr_grid = sum(circulation[axis] ** 2 / D[axis] for axis in range(3))
        drift_power_grid = sum(drift[axis] ** 2 / D[axis] for axis in range(3))
        support = mass >= support_relative * mass.max()
        supported_mass = mass * support
        supported_mass /= supported_mass.sum()
        state_fi = float(np.sum(supported_mass * grad_power_grid))
        epr = float(np.sum(supported_mass * epr_grid))
        drift_power = float(np.sum(supported_mass * drift_power_grid))
        circulation_fraction = epr / max(epr + state_fi, 1e-300)
        current_l1 = float(np.sum(np.abs(divergence)) * np.prod(spacing))
        current_scale = float(
            sum(np.sum(np.abs(probability_current[axis])) * np.prod(spacing) / spacing[axis] for axis in range(3))
        )
        force = values_at_samples(dphi[index], samples_by_p[index], edges)
        sample_variance, tau, _ = _positive_autocovariance_integral(force, sample_dt, max_lag)
        # The trajectory supplies only the normalized correlation time. The
        # covariance amplitude must come from the same normalized density and
        # centered score used to define path FI; otherwise Var(X_p) and g(p)
        # are numerically inconsistent.
        metric_variance = float(path_fi[index])
        friction = float(metric_variance * tau) if np.isfinite(tau) else np.nan
        chain_difference = float(path_fi[index] - path_fi_xy[index])
        rows.append({
            "p": float(p),
            "path_fi": float(path_fi[index]),
            "path_fi_xyz": float(path_fi[index]),
            "path_fi_xy": float(path_fi_xy[index]),
            "path_fi_z": float(path_fi_z[index]),
            "path_fi_conditional_z_given_xy": max(chain_difference, 0.0),
            "path_fi_chain_rule_violation": max(-chain_difference, 0.0),
            "path_score_mean": float(score_mean[index]),
            "path_score_mean_xy": float(score_mean_xy[index]),
            "path_score_mean_z": float(score_mean_z[index]),
            "state_fi_D": state_fi,
            "epr_proxy": epr,
            "drift_power_Dinv": drift_power,
            "circulation_fraction": circulation_fraction,
            "stationary_current_divergence_relative": current_l1 / max(current_scale, 1e-300),
            "hs_force_mean": float(np.nanmean(force)),
            "hs_force_variance": metric_variance,
            "hs_force_variance_sample": sample_variance,
            "force_variance_consistency_ratio": (
                float(sample_variance / metric_variance) if metric_variance > 0 else np.nan
            ),
            "force_correlation_time_ms": tau,
            "friction_metric": friction,
        })
    return rows, local_kl_check(masses, p_grid, path_fi), dphi, masses


def cumulative_length(p, metric):
    p, metric = np.asarray(p, float), np.asarray(metric, float)
    speed = np.sqrt(np.maximum(metric, 0.0))
    increments = 0.5 * (speed[:-1] + speed[1:]) * np.diff(p)
    cumulative = np.concatenate(([0.0], np.cumsum(increments)))
    total = float(cumulative[-1])
    normalized = cumulative / total if total > 0 else np.linspace(0.0, 1.0, len(p))
    return cumulative, normalized, total


def adaptive_grid(p, metric, n_points):
    _, normalized, _ = cumulative_length(p, metric)
    target = np.linspace(0.0, 1.0, int(n_points))
    unique, indices = np.unique(normalized, return_index=True)
    return np.interp(target, unique, np.asarray(p)[indices])


def adaptive_indices(p, metric, n_points):
    """Select a fixed number of strictly increasing grid indices.

    Dynamic programming minimizes squared thermodynamic-coordinate mismatch
    while forcing both endpoints and equal unique-state counts across all
    schedules. This avoids the repeated nearest-neighbour positions in v1.0.0.
    """
    p, metric = np.asarray(p, float), np.asarray(metric, float)
    n_grid, n_points = len(p), int(n_points)
    if not 2 <= n_points <= n_grid:
        raise ValueError("protocol.n_points must be between 2 and path.n_p")
    _, normalized, _ = cumulative_length(p, metric)
    targets = np.linspace(0.0, 1.0, n_points)
    cost = (normalized[:, None] - targets[None, :]) ** 2
    infinity = np.inf
    dp = np.full((n_grid, n_points), infinity)
    previous = np.full((n_grid, n_points), -1, dtype=int)
    dp[0, 0] = cost[0, 0]
    for target_index in range(1, n_points):
        minimum_grid = target_index
        maximum_grid = n_grid - (n_points - target_index)
        for grid_index in range(minimum_grid, maximum_grid + 1):
            candidates = dp[:grid_index, target_index - 1]
            best = int(np.argmin(candidates))
            if np.isfinite(candidates[best]):
                dp[grid_index, target_index] = candidates[best] + cost[grid_index, target_index]
                previous[grid_index, target_index] = best
    indices = np.empty(n_points, dtype=int)
    indices[-1] = n_grid - 1
    if not np.isfinite(dp[indices[-1], -1]):
        raise RuntimeError("No strictly increasing adaptive protocol exists")
    for target_index in range(n_points - 1, 0, -1):
        indices[target_index - 1] = previous[indices[target_index], target_index]
    if indices[0] != 0 or indices[-1] != n_grid - 1 or len(np.unique(indices)) != n_points:
        raise RuntimeError("Adaptive protocol uniqueness invariant failed")
    return indices
