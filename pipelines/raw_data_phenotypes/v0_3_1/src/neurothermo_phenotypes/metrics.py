from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy.signal import resample_poly, welch
from scipy.stats import rankdata
from sklearn.metrics import mutual_info_score


def uniform_resample(time_s: np.ndarray, values: np.ndarray, target_dt_s: float) -> tuple[np.ndarray, np.ndarray]:
    time_s, values = np.asarray(time_s, float), np.asarray(values, float)
    dt = float(np.median(np.diff(time_s)))
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Invalid time base")
    if abs(dt - target_dt_s) / target_dt_s < 0.02:
        return time_s, values
    ratio = target_dt_s / dt
    if ratio >= 1.0:
        down = max(1, int(round(ratio)))
        y = resample_poly(values, 1, down)
    else:
        up = max(1, int(round(1.0 / ratio)))
        y = resample_poly(values, up, 1)
    t = time_s[0] + np.arange(len(y)) * target_dt_s
    keep = t <= time_s[-1]
    return t[keep], y[keep]


def ordinal_pattern_codes(values: np.ndarray, order: int, delay: int = 1) -> np.ndarray:
    x = np.asarray(values, float)
    width = (order - 1) * delay + 1
    if order < 2 or len(x) < width:
        return np.array([], dtype=int)
    permutations = list(__import__("itertools").permutations(range(order)))
    lookup = {p: i for i, p in enumerate(permutations)}
    codes = np.empty(len(x) - width + 1, dtype=int)
    for i in range(len(codes)):
        pattern = tuple(np.argsort(x[i : i + width : delay], kind="stable"))
        codes[i] = lookup[pattern]
    return codes


def permutation_entropy(values: np.ndarray, order: int = 4, delay: int = 1) -> float:
    codes = ordinal_pattern_codes(values, order, delay)
    if len(codes) == 0:
        return np.nan
    counts = np.bincount(codes, minlength=math.factorial(order)).astype(float)
    p = counts[counts > 0] / counts.sum()
    entropy = -float(np.sum(p * np.log(p)))
    return entropy / math.log(math.factorial(order))


def ordinal_predictive_information(values: np.ndarray, order: int = 4, delay: int = 1) -> float:
    codes = ordinal_pattern_codes(values, order, delay)
    if len(codes) < 3:
        return np.nan
    return float(mutual_info_score(codes[:-1], codes[1:]))


def normalized_spectral_entropy(values: np.ndarray, sample_rate_hz: float) -> float:
    x = np.asarray(values, float)
    if len(x) < 16 or np.nanstd(x) == 0:
        return 0.0
    x = x - np.nanmean(x)
    _, power = welch(x, fs=sample_rate_hz, nperseg=min(256, len(x)))
    power = power[np.isfinite(power) & (power > 0)]
    if len(power) < 2:
        return 0.0
    p = power / power.sum()
    return float(-np.sum(p * np.log(p)) / np.log(len(p)))


def _symbols(values: np.ndarray, n_bins: int) -> tuple[np.ndarray, int]:
    x = np.asarray(values, float)
    edges = np.unique(np.quantile(x, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 4:
        return np.array([], dtype=int), 0
    inner = edges[1:-1]
    return np.digitize(x, inner, right=False).astype(int), len(edges) - 1


def block_path_kl_rate(
    values: np.ndarray,
    dt_s: float,
    n_bins: int = 6,
    word_length: int = 3,
    delay: int = 1,
    pseudocount: float = 0.5,
) -> tuple[float, float, int]:
    symbols, actual_bins = _symbols(values, n_bins)
    width = (word_length - 1) * delay + 1
    if actual_bins < 2 or len(symbols) < width + 10:
        return np.nan, np.nan, 0
    n_states = actual_bins ** word_length
    counts = np.zeros(n_states, dtype=float)
    reverse_counts = np.zeros(n_states, dtype=float)
    powers = actual_bins ** np.arange(word_length - 1, -1, -1)
    n_words = len(symbols) - width + 1
    for start in range(n_words):
        word = symbols[start : start + width : delay]
        counts[int(np.dot(word, powers))] += 1
        reverse_counts[int(np.dot(word[::-1], powers))] += 1
    p = (counts + pseudocount) / (counts.sum() + pseudocount * n_states)
    q = (reverse_counts + pseudocount) / (reverse_counts.sum() + pseudocount * n_states)
    divergence = float(np.sum(p * np.log(p / q)))
    rate = divergence / max((word_length - 1) * delay * dt_s, np.finfo(float).eps)
    coverage = float(np.count_nonzero(counts) / n_states)
    return rate, coverage, n_words


def aaft_surrogate(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = np.asarray(values, float)
    n = len(x)
    ranks = rankdata(x, method="ordinal") - 1
    gaussian = np.sort(rng.normal(size=n))[ranks]
    spectrum = np.fft.rfft(gaussian)
    phases = rng.uniform(0.0, 2.0 * np.pi, size=len(spectrum))
    phases[0] = 0.0
    if n % 2 == 0:
        phases[-1] = 0.0
    randomized = np.fft.irfft(np.abs(spectrum) * np.exp(1j * phases), n=n)
    randomized_ranks = rankdata(randomized, method="ordinal") - 1
    return np.sort(x)[randomized_ranks]


def irreversibility_with_surrogates(values: np.ndarray, dt_s: float, cfg: dict, rng: np.random.Generator) -> dict:
    args = dict(
        n_bins=int(cfg["symbol_bins"]), word_length=int(cfg["word_length"]),
        delay=max(1, int(round((cfg["word_delay_ms"] / 1000.0) / dt_s))),
        pseudocount=float(cfg["pseudocount"]),
    )
    observed, coverage, n_words = block_path_kl_rate(values, dt_s, **args)
    if not np.isfinite(observed):
        return {
            "path_kl_rate_raw_nats_s": np.nan, "path_kl_surrogate_median_nats_s": np.nan,
            "path_kl_rate_bias_corrected_nats_s": np.nan, "path_kl_rate_excess_nats_s": np.nan,
            "path_kl_surrogate_p": np.nan, "path_word_coverage": coverage, "path_n_words": n_words,
        }
    surrogate_rates = []
    for _ in range(int(cfg["n_reversible_surrogates"])):
        surrogate = aaft_surrogate(values, rng)
        rate, _, _ = block_path_kl_rate(surrogate, dt_s, **args)
        if np.isfinite(rate):
            surrogate_rates.append(rate)
    if not surrogate_rates:
        median, corrected, p_value = np.nan, np.nan, np.nan
    else:
        s = np.asarray(surrogate_rates)
        median = float(np.median(s))
        corrected = float(observed - median)
        p_value = float((1 + np.sum(s >= observed)) / (1 + len(s)))
    return {
        "path_kl_rate_raw_nats_s": observed,
        "path_kl_surrogate_median_nats_s": median,
        "path_kl_rate_bias_corrected_nats_s": corrected,
        "path_kl_rate_excess_nats_s": max(0.0, corrected) if np.isfinite(corrected) else np.nan,
        "path_kl_surrogate_p": p_value,
        "path_word_coverage": coverage,
        "path_n_words": n_words,
    }


def electrical_work(
    time_s: np.ndarray,
    voltage_mV: np.ndarray,
    current_pA: float,
    baseline_voltage_mV: float,
    current_trace_pA: Optional[np.ndarray] = None,
) -> dict:
    t, v = np.asarray(time_s, float), np.asarray(voltage_mV, float)
    source = "protocol_scalar"
    i = np.full_like(v, float(current_pA))
    if current_trace_pA is not None and len(current_trace_pA) == len(v):
        cmd = np.asarray(current_trace_pA, float)
        baseline_n = max(3, len(cmd) // 20)
        cmd = cmd - np.median(cmd[:baseline_n])
        plateau = float(np.median(cmd))
        tolerance = max(5.0, 0.15 * abs(float(current_pA)))
        if np.isfinite(plateau) and abs(plateau - current_pA) <= tolerance:
            i, source = cmd, "recorded_or_command_trace"
    power_fW = i * (v - baseline_voltage_mV)
    signed_fJ = float(np.trapz(power_fW, t))
    positive_fJ = float(np.trapz(np.maximum(power_fW, 0.0), t))
    duration = float(t[-1] - t[0]) if len(t) > 1 else np.nan
    return {
        "external_work_signed_fJ": signed_fJ,
        "external_work_positive_fJ": positive_fJ,
        "mean_power_signed_fW": signed_fJ / duration if duration > 0 else np.nan,
        "work_current_source": source,
    }
