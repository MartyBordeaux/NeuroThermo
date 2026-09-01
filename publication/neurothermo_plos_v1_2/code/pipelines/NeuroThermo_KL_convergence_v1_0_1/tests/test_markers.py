import unittest

import numpy as np

from kl_convergence.markers import curve_markers, isotonic_increasing, seed_ensemble_markers, weighted_quantile


class MarkerTests(unittest.TestCase):
    def test_pava_is_monotone(self):
        fitted = isotonic_increasing([0.0, 2.0, 1.0, 3.0])
        self.assertTrue(np.all(np.diff(fitted) >= 0))
        self.assertTrue(np.allclose(fitted, [0.0, 1.5, 1.5, 3.0]))

    def test_linear_crossing(self):
        result = curve_markers([0.0, 0.5, 1.0], [-1.0, 0.0, 1.0], 2)
        self.assertAlmostEqual(result["first"], 0.5)
        self.assertAlmostEqual(result["isotonic"], 0.5)

    def test_seed_summary(self):
        curves = np.asarray([[-1.0, 0.0, 1.0], [-2.0, -1.0, 1.0], [-1.0, 1.0, 2.0]])
        result = seed_ensemble_markers([0.0, 0.5, 1.0], curves, 2)
        self.assertTrue(np.isfinite(result["seed_median_curve_isotonic"]))
        self.assertGreaterEqual(result["seed_isotonic_iqr"], 0.0)

    def test_weighted_quantile(self):
        self.assertAlmostEqual(weighted_quantile([0, 10], [9, 1], .5), 1.0)


if __name__ == "__main__":
    unittest.main()
