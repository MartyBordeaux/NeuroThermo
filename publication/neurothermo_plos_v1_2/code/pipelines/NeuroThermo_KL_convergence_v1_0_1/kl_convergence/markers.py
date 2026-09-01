from __future__ import annotations

import numpy as np


def weighted_quantile(values, weights, probability):
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return float("nan")
    values, weights = values[mask], weights[mask]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= weights.sum()
    return float(np.interp(float(probability), cumulative, values, left=values[0], right=values[-1]))


def isotonic_increasing(values):
    values = np.asarray(values, float)
    levels = []
    counts = []
    for value in values:
        levels.append(float(value))
        counts.append(1)
        while len(levels) >= 2 and levels[-2] > levels[-1]:
            count = counts[-2] + counts[-1]
            level = (levels[-2] * counts[-2] + levels[-1] * counts[-1]) / count
            levels[-2:] = [level]
            counts[-2:] = [count]
    output = np.empty(len(values), float)
    index = 0
    for level, count in zip(levels, counts):
        output[index:index + count] = level
        index += count
    return output


def crossing_count(values):
    values = np.asarray(values, float)
    signs = np.sign(values)
    for index in range(1, len(signs)):
        if signs[index] == 0:
            signs[index] = signs[index - 1]
    return int(np.sum(signs[:-1] * signs[1:] < 0))


def first_crossing(p, values):
    p, values = np.asarray(p, float), np.asarray(values, float)
    for index in range(len(values) - 1):
        left, right = values[index], values[index + 1]
        if left == 0:
            return float(p[index])
        if left < 0 <= right:
            if right == left:
                return float(p[index])
            return float(p[index] - left * (p[index + 1] - p[index]) / (right - left))
    if values[-1] == 0:
        return float(p[-1])
    return float("nan")


def persistent_crossing(p, values, persistence_points):
    p, values = np.asarray(p, float), np.asarray(values, float)
    persistence_points = int(persistence_points)
    for index in range(len(values) - 1):
        if values[index] < 0 <= values[index + 1]:
            stop = min(len(values), index + 1 + persistence_points)
            if np.all(values[index + 1:stop] >= 0):
                return first_crossing(p[index:index + 2], values[index:index + 2])
    return float("nan")


def curve_markers(p, values, persistence_points):
    iso = isotonic_increasing(values)
    return {
        "first": first_crossing(p, values),
        "persistent": persistent_crossing(p, values, persistence_points),
        "isotonic": first_crossing(p, iso),
        "crossing_count": crossing_count(values),
        "endpoint_direction": bool(values[0] < 0 and values[-1] > 0),
    }


def seed_ensemble_markers(p, curves, persistence_points):
    curves = np.asarray(curves, float)
    median_curve = np.median(curves, axis=0)
    curve = curve_markers(p, median_curve, persistence_points)
    seed_iso = np.asarray([curve_markers(p, row, persistence_points)["isotonic"] for row in curves], float)
    finite = seed_iso[np.isfinite(seed_iso)]
    quantiles = np.quantile(finite, [0.25, 0.5, 0.75]) if len(finite) else [np.nan] * 3
    return {
        "seed_median_curve_isotonic": curve["isotonic"],
        "seed_median_curve_first": curve["first"],
        "seed_median_curve_persistent": curve["persistent"],
        "q25_seed_isotonic": float(quantiles[0]),
        "median_seed_isotonic": float(quantiles[1]),
        "q75_seed_isotonic": float(quantiles[2]),
        "seed_isotonic_iqr": float(quantiles[2] - quantiles[0]),
        "median_curve_crossing_count": int(curve["crossing_count"]),
        "endpoint_direction_fraction": float(np.mean([row[0] < 0 and row[-1] > 0 for row in curves])),
    }


MARKER_VARIANTS = (
    "seed_median_curve_isotonic",
    "seed_median_curve_first",
    "seed_median_curve_persistent",
    "q25_seed_isotonic",
    "median_seed_isotonic",
    "q75_seed_isotonic",
)
