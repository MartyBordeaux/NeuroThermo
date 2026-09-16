import numpy as np
import unittest

from nonequilibrium_geometry.geometry import adaptive_grid, adaptive_indices, cumulative_length, path_derivatives


class GeometryTests(unittest.TestCase):
    def test_path_fisher_zero_for_constant_family(self):
        masses = np.ones((5, 3, 4, 2), dtype=float)
        masses /= masses.sum(axis=(1, 2, 3), keepdims=True)
        _, fisher, score = path_derivatives(masses, np.linspace(0, 1, 5))
        self.assertTrue(np.allclose(fisher, 0.0))
        self.assertTrue(np.allclose(score, 0.0))

    def test_bernoulli_path_fisher_midpoint(self):
        p = np.linspace(0.2, 0.8, 301)
        masses = np.column_stack((p, 1 - p))
        _, fisher, score = path_derivatives(masses, p)
        middle = len(p) // 2
        self.assertTrue(np.isclose(fisher[middle], 4.0, rtol=2e-4))
        self.assertLess(abs(score[middle]), 1e-8)

    def test_constant_metric_has_linear_adaptive_grid(self):
        p = np.linspace(0, 1, 11)
        _, normalized, total = cumulative_length(p, np.full_like(p, 4.0))
        self.assertTrue(np.isclose(total, 2.0))
        self.assertTrue(np.allclose(normalized, p))
        self.assertTrue(np.allclose(adaptive_grid(p, np.ones_like(p), 7), np.linspace(0, 1, 7)))

    def test_adaptive_indices_are_unique_and_include_endpoints(self):
        p = np.linspace(0, 1, 31)
        metric = 1.0 + 100.0 * np.exp(-((p - 0.2) / 0.05) ** 2)
        indices = adaptive_indices(p, metric, 15)
        self.assertEqual(len(indices), 15)
        self.assertEqual(indices[0], 0)
        self.assertEqual(indices[-1], 30)
        self.assertTrue(np.all(np.diff(indices) > 0))

    def test_centered_fisher_is_invariant_to_constant_score_offset(self):
        p = np.linspace(0.2, 0.8, 301)
        masses = np.column_stack((p, 1 - p))
        _, fisher, score = path_derivatives(masses, p)
        self.assertTrue(np.all(fisher >= 0))
        self.assertLess(np.max(np.abs(score[2:-2])), 5e-5)
