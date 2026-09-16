from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


def common_grid(samples_by_p, cfg):
    pooled = np.concatenate(samples_by_p, axis=0)
    dcfg = cfg["density"]
    low = np.quantile(pooled, float(dcfg["quantile_low"]), axis=0)
    high = np.quantile(pooled, float(dcfg["quantile_high"]), axis=0)
    span = np.maximum(high - low, 1e-8)
    margin = float(dcfg["margin_fraction"])
    low, high = low - margin * span, high + margin * span
    bins = int(dcfg["bins"])
    edges = [np.linspace(low[i], high[i], bins + 1) for i in range(3)]
    centers = [0.5 * (edge[:-1] + edge[1:]) for edge in edges]
    return edges, centers


def full_coverage_grid(extents, cfg):
    """Build one scenario grid from extrema of every dt--seed pilot path.

    ``extents`` is an iterable of ``(axis_min, axis_max)`` pairs.  A positive
    margin protects the second deterministic pass from platform-level floating
    point differences while retaining the same bin count as the frozen
    analysis.
    """
    extents = list(extents)
    if not extents:
        raise ValueError("At least one pilot extent is required")
    lows = np.vstack([np.asarray(item[0], float) for item in extents])
    highs = np.vstack([np.asarray(item[1], float) for item in extents])
    if lows.shape[1:] != (3,) or highs.shape != lows.shape:
        raise ValueError("Pilot extents must contain three-dimensional minima and maxima")
    if not np.isfinite(lows).all() or not np.isfinite(highs).all():
        raise ValueError("Pilot extents must be finite")
    low = lows.min(axis=0)
    high = highs.max(axis=0)
    span = np.maximum(high - low, 1e-8)
    margin = float(cfg["density"].get("coverage_margin_fraction", 0.02))
    if margin <= 0:
        raise ValueError("coverage_margin_fraction must be positive")
    low = low - margin * span
    high = high + margin * span
    bins = int(cfg["density"]["bins"])
    edges = [np.linspace(low[i], high[i], bins + 1) for i in range(3)]
    centers = [0.5 * (edge[:-1] + edge[1:]) for edge in edges]
    return edges, centers


def histogram_mass(samples, edges, cfg):
    histogram, _ = np.histogramdd(samples, bins=edges)
    sigma = float(cfg["density"].get("gaussian_sigma_bins", 1.0))
    if sigma > 0:
        histogram = gaussian_filter(histogram.astype(float), sigma=sigma, mode="nearest")
    histogram += float(cfg["density"].get("pseudocount", 1e-10))
    return histogram / histogram.sum()


def build_masses(samples_by_p, edges, cfg):
    masses, retained = [], []
    for samples in samples_by_p:
        masses.append(histogram_mass(samples, edges, cfg))
        inside = np.ones(len(samples), dtype=bool)
        for axis in range(3):
            inside &= samples[:, axis] >= edges[axis][0]
            inside &= samples[:, axis] <= edges[axis][-1]
        retained.append(float(np.mean(inside)))
    return np.asarray(masses), np.asarray(retained)


def bin_indices(samples, edges):
    indices = np.column_stack([np.searchsorted(edge, samples[:, axis], side="right") - 1 for axis, edge in enumerate(edges)])
    shape = tuple(len(edge) - 1 for edge in edges)
    valid = np.ones(len(samples), dtype=bool)
    for axis, length in enumerate(shape):
        valid &= (indices[:, axis] >= 0) & (indices[:, axis] < length)
    clipped = np.column_stack([np.clip(indices[:, axis], 0, shape[axis] - 1) for axis in range(3)])
    return clipped.astype(int), valid


def values_at_samples(field, samples, edges):
    indices, valid = bin_indices(samples, edges)
    values = field[indices[:, 0], indices[:, 1], indices[:, 2]].astype(float)
    values[~valid] = np.nan
    return values


def quantile_state_edges(samples_by_p, state_shape):
    pooled = np.concatenate(samples_by_p, axis=0)
    edges = []
    for axis, count in enumerate(state_shape):
        internal = np.quantile(pooled[:, axis], np.linspace(0.0, 1.0, int(count) + 1)[1:-1])
        edges.append(np.concatenate(([-np.inf], np.unique(internal), [np.inf])))
    return edges


def assign_states(samples, state_edges):
    dimensions = [len(edge) - 1 for edge in state_edges]
    indices = [np.searchsorted(edge, samples[:, axis], side="right") - 1 for axis, edge in enumerate(state_edges)]
    for axis in range(3):
        indices[axis] = np.clip(indices[axis], 0, dimensions[axis] - 1)
    return np.ravel_multi_index(tuple(indices), tuple(dimensions)), int(np.prod(dimensions))
