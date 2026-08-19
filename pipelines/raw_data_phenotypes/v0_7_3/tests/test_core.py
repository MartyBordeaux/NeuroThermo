from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from neurothermo_rank_gaussian_fourier.metrics import (
    fourier_fidelity,
    fourier_phase_surrogate,
    ordinal_pattern_codes,
    ordinal_pi_lags,
    rank_gaussianize,
)
from neurothermo_rank_gaussian_fourier.pipeline import _fourier_worker, _seed, _surrogate_statistics
from neurothermo_rank_gaussian_fourier.raw import resolve_source_path
from neurothermo_rank_gaussian_fourier.statistics import exact_lag_family_tests, label_masks


def test_rank_gaussianization_preserves_ordinal_codes_including_ties():
    values = np.asarray([2.0, 2.0, -1.0, 4.0, 0.5, 0.5, 3.0, -2.0, 1.0, 1.0])
    gaussianized = rank_gaussianize(values)
    assert np.array_equal(
        ordinal_pattern_codes(values, order=4, delay=1),
        ordinal_pattern_codes(gaussianized, order=4, delay=1),
    )


def test_rank_gaussianization_preserves_lagged_pi():
    rng = np.random.default_rng(9)
    values = rng.normal(size=901) ** 3
    transformed = rank_gaussianize(values)
    assert np.array_equal(
        ordinal_pi_lags(values, 4, 1, [4, 8, 16, 32]),
        ordinal_pi_lags(transformed, 4, 1, [4, 8, 16, 32]),
    )


def test_fourier_surrogate_preserves_spectrum_and_circular_acf():
    rng = np.random.default_rng(5)
    for length in [900, 901]:
        time = np.arange(length)
        values = np.sin(2 * np.pi * time / 31.0) + 0.3 * rng.normal(size=len(time))
        surrogate = fourier_phase_surrogate(values, rng)
        diagnostic = fourier_fidelity(values, surrogate, 32)
        assert diagnostic["spectral_amplitude_nrmse"] < 1e-12
        assert diagnostic["circular_acf_max_abs_error"] < 1e-12


def test_fourier_worker_is_deterministic():
    values = rank_gaussianize(np.random.default_rng(4).normal(size=160))
    payload = ("WT_01", 400, values, 4, 1, [4, 8], 4, 20260818, 16)
    first = _fourier_worker(payload)
    second = _fourier_worker(payload)
    assert np.array_equal(first[0], second[0])
    assert first[1] == second[1]
    assert _seed(1, "WT_01", 100) != _seed(1, "WT_01", 150)


def test_tail_probabilities_include_lower_upper_and_two_sided():
    result = _surrogate_statistics(0.0, np.asarray([1.0, 2.0, 3.0]), "x")
    assert np.isclose(result["x_p_lower"], 0.25)
    assert np.isclose(result["x_p_upper"], 1.0)
    assert 0.0 < result["x_p_two_sided"] <= 1.0


def test_exact_lag_tests_apply_maxT_and_keep_primary_flag():
    rows = []
    cells = ["w1", "w2", "s1", "s2"]
    groups = ["WT", "WT", "SCA3", "SCA3"]
    for cell, group, offset in zip(cells, groups, [2.0, 1.0, -1.0, -2.0]):
        for current in [100, 200]:
            rows.append({
                "cell_id": cell, "group": group, "current_pA": current,
                "shuffle_lag_4ms_centered_PI_nats": offset,
                "shuffle_lag_8ms_centered_PI_nats": offset / 2.0,
                "fourier_lag_4ms_centered_PI_nats": offset * 0.8,
                "fourier_lag_8ms_centered_PI_nats": offset * 0.4,
            })
    table = pd.DataFrame(rows)
    labels = np.asarray([False, False, True, True])
    masks = label_masks(4, 2, 10)
    families = {
        "shuffle": {4: "shuffle_lag_4ms_centered_PI_nats", 8: "shuffle_lag_8ms_centered_PI_nats"},
        "fourier": {4: "fourier_lag_4ms_centered_PI_nats", 8: "fourier_lag_8ms_centered_PI_nats"},
    }
    auc, current, cell_values = exact_lag_family_tests(
        table, families, cells, labels, [100, 200], masks, -1.0
    )
    assert len(auc) == 4
    assert len(current) == 8
    assert len(cell_values) == 16
    assert bool(auc[(auc.surrogate_family == "shuffle") & (auc.lag_ms == 4)].primary_prespecified.iloc[0])
    assert not bool(auc[(auc.surrogate_family == "fourier") & (auc.lag_ms == 4)].primary_prespecified.iloc[0])


def test_raw_path_fallback_uses_group_directory(tmp_path: Path):
    target = tmp_path / "SCA3" / "cc_01.npz"
    target.parent.mkdir()
    target.write_bytes(b"x")
    assert resolve_source_path("/old/cc_01.npz", "SCA3", tmp_path) == target.resolve()
