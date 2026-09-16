import unittest

import numpy as np

from kl_convergence.density import build_masses, full_coverage_grid


class FullCoverageGridTests(unittest.TestCase):
    def config(self):
        return {
            "density": {
                "bins": 8,
                "coverage_margin_fraction": 0.02,
                "gaussian_sigma_bins": 0.0,
                "pseudocount": 1e-10,
            }
        }

    def test_grid_contains_all_pilot_samples(self):
        left = np.asarray([[-3.0, -1.0, 0.0], [-2.0, 2.0, 1.0]])
        right = np.asarray([[4.0, -4.0, -2.0], [5.0, 3.0, 6.0]])
        extents = [(left.min(axis=0), left.max(axis=0)),
                   (right.min(axis=0), right.max(axis=0))]
        edges, _ = full_coverage_grid(extents, self.config())
        _, retained = build_masses([left, right], edges, self.config())
        self.assertTrue(np.array_equal(retained, np.ones(2)))

    def test_positive_margin_is_required(self):
        cfg = self.config()
        cfg["density"]["coverage_margin_fraction"] = 0.0
        with self.assertRaises(ValueError):
            full_coverage_grid([([0, 0, 0], [1, 1, 1])], cfg)


if __name__ == "__main__":
    unittest.main()
