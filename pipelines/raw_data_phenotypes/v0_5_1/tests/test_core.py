import numpy as np
import pandas as pd

from neurothermo_cell_fingerprint.pipeline import _bh, _cliffs_delta, _transform, exact_label_masks


def test_transforms_and_direction_helpers():
    assert np.allclose(_transform([0, 1], "log1p"), [0, np.log(2)])
    assert np.isnan(_transform([0], "log")[0])
    assert _cliffs_delta([2, 3], [0, 1]) == 1.0


def test_bh_is_monotone_in_rank():
    q = _bh(np.array([0.03, 0.01, 0.02]))
    table = pd.DataFrame({"p": [0.03, 0.01, 0.02], "q": q}).sort_values("p")
    assert np.all(np.diff(table["q"]) >= -1e-12)


def test_exact_masks_enumerate_all_assignments():
    labels = np.array(["WT", "WT", "WT", "SCA3", "SCA3"])
    masks = exact_label_masks(labels, exact_max=100)
    assert masks.shape == (10, 5)
    assert np.all(masks.sum(axis=1) == 3)
