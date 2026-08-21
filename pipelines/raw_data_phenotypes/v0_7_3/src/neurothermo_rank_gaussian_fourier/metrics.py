from __future__ import annotations

import itertools
import math
from typing import Dict, Sequence, Tuple

import numpy as np
from scipy.signal import resample_poly
from scipy.stats import kurtosis, norm, skew


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


def ordinal_pi_lags(values: np.ndarray, order: int, delay: int, code_lags: Sequence[int]) -> np.ndarray:
    codes = ordinal_pattern_codes(values, order, delay)
    states = math.factorial(order)
    result = []
    for lag in code_lags:
        lag = int(lag)
        result.append(mutual_information_codes(codes[:-lag], codes[lag:], states) if len(codes) > lag else np.nan)
    return np.asarray(result, dtype=float)


def rank_gaussianize(values: np.ndarray) -> np.ndarray:
    """Stable empirical normal-score transform preserving the ordinal code sequence."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("rank_gaussianize requires a finite one-dimensional trace")
    order = np.argsort(x, kind="stable")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(len(x), dtype=float)
    probabilities = (ranks + 0.5) / float(len(x))
    return norm.ppf(probabilities)


def fourier_phase_surrogate(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Randomize Fourier phases while preserving every rFFT magnitude."""
    x = np.asarray(values, dtype=float)
    n = len(x)
    spectrum = np.fft.rfft(x)
    target_amplitude = np.abs(spectrum)
    phases = rng.uniform(0.0, 2.0 * np.pi, size=len(spectrum))
    randomized = target_amplitude * np.exp(1j * phases)
    randomized[0] = spectrum[0]
    if n % 2 == 0:
        randomized[-1] = target_amplitude[-1] * (1.0 if rng.random() < 0.5 else -1.0)
    return np.fft.irfft(randomized, n=n)


def normalized_linear_autocorrelation(values: np.ndarray, maximum_lag: int) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x - float(np.mean(x))
    denominator = float(np.dot(x, x))
    if denominator <= np.finfo(float).eps:
        return np.full(int(maximum_lag), np.nan)
    return np.asarray([
        float(np.dot(x[:-lag], x[lag:]) / denominator)
        for lag in range(1, int(maximum_lag) + 1)
    ])


def normalized_circular_autocorrelation(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x - float(np.mean(x))
    spectrum = np.fft.rfft(x)
    acf = np.fft.irfft(np.abs(spectrum) ** 2, n=len(x))
    if abs(float(acf[0])) <= np.finfo(float).eps:
        return np.full(len(x), np.nan)
    return np.asarray(acf / float(acf[0]), dtype=float)


def fourier_fidelity(values: np.ndarray, surrogate: np.ndarray, maximum_linear_acf_lag: int) -> Dict[str, float]:
    x = np.asarray(values, dtype=float)
    y = np.asarray(surrogate, dtype=float)
    target_amplitude = np.abs(np.fft.rfft(x))
    surrogate_amplitude = np.abs(np.fft.rfft(y))
    denominator = float(np.linalg.norm(target_amplitude))
    spectral_nrmse = 0.0 if denominator <= np.finfo(float).eps else float(
        np.linalg.norm(surrogate_amplitude - target_amplitude) / denominator
    )
    circular_difference = normalized_circular_autocorrelation(y) - normalized_circular_autocorrelation(x)
    linear_difference = (
        normalized_linear_autocorrelation(y, maximum_linear_acf_lag) -
        normalized_linear_autocorrelation(x, maximum_linear_acf_lag)
    )
    return {
        "spectral_amplitude_nrmse": spectral_nrmse,
        "spectral_amplitude_max_abs_error": float(np.max(np.abs(surrogate_amplitude - target_amplitude))),
        "circular_acf_rmse": float(np.sqrt(np.nanmean(circular_difference ** 2))),
        "circular_acf_max_abs_error": float(np.nanmax(np.abs(circular_difference))),
        "linear_acf_rmse": float(np.sqrt(np.nanmean(linear_difference ** 2))),
        "surrogate_mean": float(np.mean(y)),
        "surrogate_sd": float(np.std(y, ddof=1)),
        "surrogate_skewness": float(skew(y, bias=False)),
        "surrogate_excess_kurtosis": float(kurtosis(y, fisher=True, bias=False)),
    }
