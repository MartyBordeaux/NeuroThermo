from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from neurothermo_pi_surrogates.metrics import (
    aaft_surrogate,
    mutual_information_codes,
    ordinal_pattern_codes,
    ordinal_predictive_information,
)
from neurothermo_pi_surrogates.pipeline import _seed, _surrogate_worker
from neurothermo_pi_surrogates.raw import resolve_source_path
from neurothermo_pi_surrogates.statistics import (
    cell_auc,
    exact_cell_difference,
    label_masks,
    ridge_crossfit_by_cell_current,
)


def test_ordinal_pi_matches_direct_count_definition():
    values = np.asarray([0.2, 1.0, -0.3, 0.7, 0.4, -1.0, 0.1, 0.8, -0.2])
    codes = ordinal_pattern_codes(values, order=4, delay=1)
    direct = mutual_information_codes(codes[:-1], codes[1:], 24)
    assert np.isclose(ordinal_predictive_information(values, 4, 1), direct)


def test_overlapping_ordinal_patterns_have_positive_shuffle_null():
    rng = np.random.default_rng(7)
    values = rng.normal(size=901)
    shuffled = rng.permutation(values)
    assert ordinal_predictive_information(shuffled, 4, 1) > 0.5


def test_aaft_preserves_amplitude_distribution_exactly():
    values = np.linspace(-2.0, 3.0, 101) ** 3
    surrogate = aaft_surrogate(values, np.random.default_rng(10))
    assert np.array_equal(np.sort(values), np.sort(surrogate))


def test_surrogate_worker_is_deterministic():
    values = np.random.default_rng(4).normal(size=120)
    payload = ("WT_01", 400, values, 4, 1, 5, 5, 20260818)
    first = _surrogate_worker(payload)
    second = _surrogate_worker(payload)
    assert np.array_equal(first[0], second[0])
    assert np.array_equal(first[1], second[1])
    assert _seed(1, "a", 100, "shuffle") != _seed(1, "a", 100, "aaft")


def test_exact_cell_difference_enumerates_all_labelings():
    values = np.asarray([0.0, 1.0, 5.0, 6.0])
    masks = label_masks(4, 2, 10)
    observed_labels = np.asarray([False, False, True, True])
    observed, p_expected, p_two, _ = exact_cell_difference(
        values, masks, observed_labels, expected_direction=1.0
    )
    assert observed == 5.0
    assert np.isclose(p_expected, 1.0 / 6.0)
    assert np.isclose(p_two, 2.0 / 6.0)


def test_cell_auc_is_normalized():
    matrix = np.asarray([[0.0, 1.0, 2.0], [4.0, 4.0, 4.0]])
    assert np.allclose(cell_auc(matrix, [100, 150, 200]), [1.0, 4.0])


def test_crossfit_excludes_target_cell():
    rows = []
    cells = ["c1", "c2", "c3", "c4"]
    for current in [100, 150]:
        for index, cell in enumerate(cells):
            rows.append({
                "cell_id": cell, "current_pA": current,
                "outcome": 2.0 + 0.5 * index, "x": float(index),
            })
    table = pd.DataFrame(rows)
    first, _, _ = ridge_crossfit_by_cell_current(
        table, cells, [100, 150], "outcome", ["x"], 1.0
    )
    changed = table.copy()
    changed.loc[changed.cell_id == "c1", "outcome"] = 2000.0
    second, _, _ = ridge_crossfit_by_cell_current(
        changed, cells, [100, 150], "outcome", ["x"], 1.0
    )
    target = table.cell_id == "c1"
    assert np.allclose(first[target], second[target])


def test_raw_path_fallback_uses_group_directory(tmp_path: Path):
    target = tmp_path / "WT" / "cc_01.npz"
    target.parent.mkdir()
    target.write_bytes(b"x")
    resolved = resolve_source_path("/old/server/cc_01.npz", "WT", tmp_path)
    assert resolved == target.resolve()
