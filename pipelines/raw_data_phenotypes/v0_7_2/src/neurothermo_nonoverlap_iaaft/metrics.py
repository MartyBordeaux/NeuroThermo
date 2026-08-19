from __future__ import annotations

import itertools
import math
from typing import Dict, Sequence, Tuple

import numpy as np
from scipy.signal import resample_poly


def uniform_resample(time_s: np.ndarray, values: np.ndarray, target_dt_s: float) -> Tuple[np.ndarray, np.ndarray]:
    time_s = np.asarray(time_s, dtype=float)
    values = np.asarray(values, dtype=float)
    dt = float(np.median(np.diff(time_s)))
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Invalid time base")
    if abs(dt - target_dt_s) / target_dt_s < 0.02:
        return time_s, values
    ratio = target_dt_s / dt
    if ratio >= 1.0:
        down = max(1, int(round(ratio)))
        resampled = resample_poly(values, 1, down)
    else:
        up = max(1, int(round(1.0 / ratio)))
        resampled = resample_poly(values, up, 1)
    new_time = time_s[0] + np.arange(len(resampled)) * target_dt_s
    keep = new_time <= time_s[-1]
    return new_time[keep], resampled[keep]


def ordinal_pattern_codes(values: np.ndarray, order: int, delay: int = 1) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    width = (order - 1) * delay + 1
    if order < 2 or len(x) < width:
        return np.asarray([], dtype=np.int16)
    windows = np.lib.stride_tricks.sliding_window_view(x, width)[:, ::delay]
    permutations = list(itertools.permutations(range(order)))
    lookup = {pattern: index for index, pattern in enumerate(permutations)}
    ranks = np.argsort(windows, axis=1, kind="stable")
    return np.fromiter((lookup[tuple(row)] for row in ranks), dtype=np.int16, count=len(ranks))


def mutual_information_codes(first: np.ndarray, second: np.ndarray, n_states: int) -> float:
    first = np.asarray(first, dtype=np.int64)
    second = np.asarray(second, dtype=np.int64)
    if len(first) == 0 or len(first) != len(second):
        return np.nan
    joint = np.bincount(first * n_states + second, minlength=n_states * n_states).reshape(n_states, n_states)
    total = float(joint.sum())
    if total <= 0:
        return np.nan
    rows = joint.sum(axis=1)
    cols = joint.sum(axis=0)
    nz_i, nz_j = np.nonzero(joint)
    counts = joint[nz_i, nz_j].astype(float)
    return float(np.sum((counts / total) * np.log((counts * total) / (rows[nz_i] * cols[nz_j]))))


def ordinal_predictive_information(values: np.ndarray, order: int = 4, delay: int = 1, code_lag: int = 1) -> float:
    codes = ordinal_pattern_codes(values, order, delay)
    lag = int(code_lag)
    if lag < 1 or len(codes) <= lag:
        return np.nan
    return mutual_information_codes(codes[:-lag], codes[lag:], math.factorial(order))


def ordinal_pi_lags(values: np.ndarray, order: int, delay: int, code_lags: Sequence[int]) -> np.ndarray:
    codes = ordinal_pattern_codes(values, order, delay)
    states = math.factorial(order)
    output = []
    for lag in code_lags:
        lag = int(lag)
        output.append(mutual_information_codes(codes[:-lag], codes[lag:], states) if len(codes) > lag else np.nan)
    return np.asarray(output, dtype=float)


def shuffled_surrogate(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    return rng.permutation(np.asarray(values, dtype=float))


def _spectral_amplitude_nrmse(values: np.ndarray, target_amplitude: np.ndarray) -> float:
    amplitude = np.abs(np.fft.rfft(np.asarray(values, dtype=float)))
    denominator = float(np.linalg.norm(target_amplitude))
    if denominator <= np.finfo(float).eps:
        return 0.0
    return float(np.linalg.norm(amplitude - target_amplitude) / denominator)


def iaaft_surrogate(
    values: np.ndarray,
    rng: np.random.Generator,
    max_iterations: int = 200,
    improvement_tolerance: float = 1e-8,
    patience: int = 20,
) -> Tuple[np.ndarray, Dict[str, object]]:
    x = np.asarray(values, dtype=float)
    target_sorted = np.sort(x)
    target_amplitude = np.abs(np.fft.rfft(x))
    surrogate = rng.permutation(x)
    best_error = np.inf
    stagnant = 0
    stop_reason = "max_iterations"
    converged = False
    iterations = 0
    for iteration in range(1, int(max_iterations) + 1):
        spectrum = np.fft.rfft(surrogate)
        magnitude = np.abs(spectrum)
        phase = np.ones_like(spectrum, dtype=complex)
        nonzero = magnitude > np.finfo(float).eps
        phase[nonzero] = spectrum[nonzero] / magnitude[nonzero]
        projected = np.fft.irfft(target_amplitude * phase, n=len(x))
        order = np.argsort(projected, kind="stable")
        updated = np.empty_like(projected)
        updated[order] = target_sorted
        error = _spectral_amplitude_nrmse(updated, target_amplitude)
        iterations = iteration
        if np.array_equal(updated, surrogate):
            surrogate = updated
            stop_reason = "fixed_rank_order"
            converged = True
            break
        improvement = best_error - error
        scale = max(1.0, abs(best_error)) if np.isfinite(best_error) else 1.0
        if not np.isfinite(best_error) or improvement > float(improvement_tolerance) * scale:
            best_error = error
            stagnant = 0
        else:
            stagnant += 1
        surrogate = updated
        if stagnant >= int(patience):
            stop_reason = "spectral_error_plateau"
            converged = True
            break
    diagnostics = {
        "iterations": int(iterations),
        "converged": bool(converged),
        "stop_reason": stop_reason,
        "spectral_amplitude_nrmse": _spectral_amplitude_nrmse(surrogate, target_amplitude),
        "amplitude_sorted_max_abs_error": float(np.max(np.abs(np.sort(surrogate) - target_sorted))),
        "amplitude_distribution_exact": bool(np.array_equal(np.sort(surrogate), target_sorted)),
    }
    return surrogate, diagnostics


def normalized_autocorrelation(values: np.ndarray, maximum_lag: int) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x - float(np.mean(x))
    denominator = float(np.dot(x, x))
    if denominator <= np.finfo(float).eps:
        return np.full(int(maximum_lag), np.nan)
    return np.asarray([
        float(np.dot(x[:-lag], x[lag:]) / denominator)
        for lag in range(1, int(maximum_lag) + 1)
    ])


def surrogate_fidelity(values: np.ndarray, surrogate: np.ndarray, maximum_acf_lag: int) -> Dict[str, float]:
    x = np.asarray(values, dtype=float)
    y = np.asarray(surrogate, dtype=float)
    target_amplitude = np.abs(np.fft.rfft(x))
    target_psd = target_amplitude ** 2
    surrogate_psd = np.abs(np.fft.rfft(y)) ** 2
    if len(target_psd) > 1:
        target_psd = target_psd[1:]
        surrogate_psd = surrogate_psd[1:]
    positive = target_psd[target_psd > 0]
    epsilon = float(np.median(positive) * 1e-12) if len(positive) else np.finfo(float).tiny
    log_psd_rmse = float(np.sqrt(np.mean((np.log(target_psd + epsilon) - np.log(surrogate_psd + epsilon)) ** 2)))
    original_acf = normalized_autocorrelation(x, maximum_acf_lag)
    surrogate_acf = normalized_autocorrelation(y, maximum_acf_lag)
    acf_rmse = float(np.sqrt(np.nanmean((original_acf - surrogate_acf) ** 2)))
    return {
        "spectral_amplitude_nrmse": _spectral_amplitude_nrmse(y, target_amplitude),
        "log_psd_rmse": log_psd_rmse,
        "acf_rmse": acf_rmse,
        "amplitude_sorted_max_abs_error": float(np.max(np.abs(np.sort(y) - np.sort(x)))),
        "amplitude_distribution_exact": bool(np.array_equal(np.sort(y), np.sort(x))),
    }
