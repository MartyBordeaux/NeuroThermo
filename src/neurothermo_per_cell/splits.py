from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pandas as pd


def current_level_splits(cell: pd.DataFrame, n_folds: int, include_full: bool = True):
    order = cell.sort_values("current_density_pA_per_pF").index.to_numpy()
    positive = order[1:] if len(order) > 1 else order
    splits = []
    if include_full:
        splits.append(("full", order, np.array([], dtype=int)))
    for fold in range(n_folds):
        test = positive[fold::n_folds]
        train = np.array([i for i in order if i not in set(test)], dtype=int)
        splits.append(("cv%d" % fold, train, test))
    return splits
