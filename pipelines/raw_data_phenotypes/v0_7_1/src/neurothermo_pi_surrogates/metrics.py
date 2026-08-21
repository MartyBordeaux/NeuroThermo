from __future__ import annotations

import itertools
import math
from typing import Tuple

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


def ordinal_predictive_information(values: np.ndarray, order: int = 4, delay: int = 1) -> float:
    codes = ordinal_pattern_codes(values, order, delay)
    if len(codes) < 3:
        return np.nan
    return mutual_information_codes(codes[:-1], codes[1:], math.factorial(order))


def shuffled_surrogate(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    return rng.permutation(np.asarray(values, dtype=float))


def aaft_surrogate(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    n = len(x)
    ranks = np.empty(n, dtype=int)
    ranks[np.argsort(x, kind="stable")] = np.arange(n)
    gaussian = np.sort(rng.normal(size=n))[ranks]
    spectrum = np.fft.rfft(gaussian)
    phases = rng.uniform(0.0, 2.0 * np.pi, size=len(spectrum))
    phases[0] = 0.0
    if n % 2 == 0:
        phases[-1] = 0.0
    randomized = np.fft.irfft(np.abs(spectrum) * np.exp(1j * phases), n=n)
    randomized_ranks = np.empty(n, dtype=int)
    randomized_ranks[np.argsort(randomized, kind="stable")] = np.arange(n)
    return np.sort(x)[randomized_ranks]
