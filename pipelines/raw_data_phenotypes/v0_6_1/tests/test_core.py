import unittest

import numpy as np

from neurothermo_vulnerability.pipeline import _combine_excitability, _first_persistent, _robust_reference, _transform


class CoreTests(unittest.TestCase):
    def test_first_persistent_requires_adjacent_levels(self):
        currents = [100, 150, 200, 250]
        self.assertEqual(_first_persistent([True, False, True, True], currents, 2, 300), (200, False))
        self.assertEqual(_first_persistent([True, False, True, False], currents, 2, 300), (300, True))


    def test_transform_does_not_impute_invalid_values(self):
        values = np.asarray([0.0, 1.0, np.nan, -1.0])
        logged = _transform(values, "log")
        self.assertTrue(np.isnan(logged[0]))
        self.assertEqual(logged[1], 0.0)
        self.assertTrue(np.isnan(logged[2]))
        self.assertTrue(np.isnan(logged[3]))


    def test_zero_mad_uses_frozen_fallback(self):
        center, scale, n, fallback = _robust_reference(np.ones(6), 0.5)
        self.assertEqual(center, 1.0)
        self.assertEqual(scale, 0.5)
        self.assertEqual(n, 6)
        self.assertTrue(fallback)

    def test_excitability_combines_redundant_subfeatures_into_one_domain(self):
        firing = np.asarray([2.0, 4.0, np.nan])
        isi = np.asarray([4.0, np.nan, 3.0])
        domain, count = _combine_excitability(firing, isi)
        np.testing.assert_allclose(domain, [3.0, 4.0, 3.0])
        np.testing.assert_array_equal(count, [2, 1, 1])


if __name__ == "__main__":
    unittest.main()
