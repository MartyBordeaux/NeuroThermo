import unittest

import numpy as np

from neurothermo_per_cell.data import load_observations, validate_grid
from neurothermo_per_cell.models import SPECS, simulate
from neurothermo_per_cell.splits import current_level_splits


class PipelineTests(unittest.TestCase):
    def test_frozen_data(self):
        audit = validate_grid(load_observations())
        self.assertEqual(audit["rows"], 312)
        self.assertEqual(audit["cells"], 24)
        self.assertEqual(audit["groups"], {"SCA3": 9, "WT": 15})

    def test_splits_do_not_overlap(self):
        cell = load_observations().query("cell_id == 'WT_02'")
        for name, train, test in current_level_splits(cell, 3):
            self.assertFalse(set(train).intersection(test))
            if name != "full":
                self.assertIn(cell.index[0], train)

    def test_models_are_finite(self):
        for name, spec in SPECS.items():
            midpoint = [(lo + hi) / 2 for lo, hi in spec.bounds]
            result = simulate(name, midpoint, np.array([0.0, 1.0]), dt_ms=0.5, duration_ms=20)
            self.assertEqual(len(result), 2)
            self.assertTrue(np.isfinite(result[0]["pred_sustained_rate_hz"]))


if __name__ == "__main__":
    unittest.main()
