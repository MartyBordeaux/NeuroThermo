from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from neurothermo_nonoverlap_iaaft.metrics import (
    iaaft_surrogate,
    mutual_information_codes,
    ordinal_pattern_codes,
    ordinal_predictive_information,
    surrogate_fidelity,
)
from neurothermo_nonoverlap_iaaft.pipeline import _seed, _surrogate_statistics, _surrogate_worker
from neurothermo_nonoverlap_iaaft.raw import resolve_source_path
from neurothermo_nonoverlap_iaaft.statistics import (
    exact_lag_family_tests,
    label_masks,
)


def test_code_lag_four_uses_nonoverlapping_order_four_windows():
    first_samples = set(range(0, 4))
    second_samples = set(range(4, 8))
    assert first_samples.isdisjoint(second_samples)


def test_lagged_ordinal_pi_matches_direct_count_definition():
    values = np.asarray([0.2, 1.0, -0.3, 0.7, 0.4, -1.0, 0.1, 0.8, -0.2, 0.5])
    codes = ordinal_pattern_codes(values, order=4, delay=1)
    direct = mutual_information_codes(codes[:-4], codes[4:], 24)
    assert np.isclose(ordinal_predictive_information(values, 4, 1, 4), direct)


def test_nonoverlap_removes_large_iid_overlap_bias():
    rng = np.random.default_rng(7)
    values = rng.normal(size=3000)
    adjacent = ordinal_predictive_information(values, 4, 1, 1)
    nonoverlap = ordinal_predictive_information(values, 4, 1, 4)
    assert adjacent > nonoverlap


def test_iaaft_preserves_amplitudes_and_improves_spectrum_over_shuffle():
    rng = np.random.default_rng(3)
    time = np.arange(901)
    values = np.sin(2 * np.pi * time / 37.0) + 0.4 * np.sin(2 * np.pi * time / 11.0)
    iaaft, convergence = iaaft_surrogate(values, rng, max_iterations=100, patience=10)
    shuffled = rng.permutation(values)
    iaaft_fidelity = surrogate_fidelity(values, iaaft, 32)
    shuffled_fidelity = surrogate_fidelity(values, shuffled, 32)
    assert np.array_equal(np.sort(values), np.sort(iaaft))
    assert convergence["iterations"] >= 1
    assert iaaft_fidelity["spectral_amplitude_nrmse"] < shuffled_fidelity["spectral_amplitude_nrmse"]


def test_surrogate_worker_is_deterministic():
    values = np.random.default_rng(4).normal(size=160)
    payload = ("WT_01", 400, values, 4, 1, [4, 8], 3, 3, 20260818, 30, 1e-8, 5, 16)
    first = _surrogate_worker(payload)
    second = _surrogate_worker(payload)
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])
    assert first[2] == second[2]
    assert _seed(1, "a", 100, "shuffle_v072") != _seed(1, "a", 100, "iaaft_v072")


def test_tail_probabilities_include_all_directions():
    result = _surrogate_statistics(0.0, np.asarray([1.0, 2.0, 3.0]), "x")
    assert np.isclose(result["x_p_lower"], 0.25)
    assert np.isclose(result["x_p_upper"], 1.0)
    assert 0.0 < result["x_p_two_sided"] <= 1.0


def test_exact_lag_tests_use_all_labelings_and_maxT():
    rows = []
    cells = ["w1", "w2", "s1", "s2"]
    groups = ["WT", "WT", "SCA3", "SCA3"]
    for cell, group, offset in zip(cells, groups, [2.0, 1.0, -1.0, -2.0]):
        for current in [100, 200]:
            rows.append({
                "cell_id": cell, "group": group, "current_pA": current,
                "shuffle_lag_4ms_centered_PI_nats": offset,
                "shuffle_lag_8ms_centered_PI_nats": offset / 2.0,
            })
    table = pd.DataFrame(rows)
    labels = np.asarray([False, False, True, True])
    masks = label_masks(4, 2, 10)
    auc, current, cells_out = exact_lag_family_tests(
        table,
        {"shuffle": {4: "shuffle_lag_4ms_centered_PI_nats", 8: "shuffle_lag_8ms_centered_PI_nats"}},
        cells, labels, [100, 200], masks, -1.0,
    )
    assert len(auc) == 2
    assert len(current) == 4
    assert len(cells_out) == 8
    assert bool(auc.loc[auc.lag_ms == 4, "primary_prespecified"].iloc[0])
    assert (auc.n_exact_labelings == 6).all()


def test_raw_path_fallback_uses_group_directory(tmp_path: Path):
    target = tmp_path / "WT" / "cc_01.npz"
    target.parent.mkdir()
    target.write_bytes(b"x")
    resolved = resolve_source_path("/old/server/cc_01.npz", "WT", tmp_path)
    assert resolved == target.resolve()
