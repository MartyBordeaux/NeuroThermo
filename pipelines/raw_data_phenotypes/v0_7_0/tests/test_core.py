from __future__ import annotations

import numpy as np
import pandas as pd

from neurothermo_predictive_validation.pipeline import (
    _cell_auc,
    _exact_cell_difference,
    _ridge_crossfit_by_cell_current,
)


def test_crossfit_excludes_entire_target_cell():
    rows = []
    cells = ["c1", "c2", "c3", "c4"]
    for current in [100, 150]:
        for index, cell in enumerate(cells):
            x = float(index)
            rows.append({
                "cell_id": cell,
                "current_pA": current,
                "log_predictive_information": 2.0 + 0.5 * x,
                "x": x,
            })
    table = pd.DataFrame(rows)
    first, _, _ = _ridge_crossfit_by_cell_current(table, cells, [100, 150], ["x"], 1.0)
    changed = table.copy()
    changed.loc[changed.cell_id == "c1", "log_predictive_information"] = 2000.0
    second, _, _ = _ridge_crossfit_by_cell_current(changed, cells, [100, 150], ["x"], 1.0)
    target = table.cell_id == "c1"
    assert np.allclose(first[target], second[target])


def test_cell_auc_is_normalized_by_current_span():
    values = np.asarray([[0.0, 1.0, 2.0], [4.0, 4.0, 4.0]])
    result = _cell_auc(values, [100, 150, 200])
    assert np.allclose(result, [1.0, 4.0])


def test_exact_cell_difference_uses_all_labelings():
    values = np.asarray([0.0, 1.0, 5.0, 6.0])
    masks = np.asarray([
        [True, True, False, False],
        [True, False, True, False],
        [True, False, False, True],
        [False, True, True, False],
        [False, True, False, True],
        [False, False, True, True],
    ])
    observed_labels = masks[-1]
    observed, p_expected, p_two = _exact_cell_difference(values, masks, observed_labels, 1.0)
    assert observed == 5.0
    assert np.isclose(p_expected, 1.0 / 6.0)
    assert np.isclose(p_two, 2.0 / 6.0)
