import unittest

import numpy as np

from kl_convergence.model import stationary_samples_nested


class NestedIntegratorTests(unittest.TestCase):
    def config(self):
        return {
            "model": {"a": 1.0, "c": 1.0, "d": 5.0, "x_R": -1.6,
                      "x0": -1.6, "y0": -10.0, "z0": 1.0, "model_time_scale_ms": 1.0},
            "noise": {"D": [0.0025, 0.01, 0.00025], "multiplier": 1.0},
            "stationary": {"burn_ms": 1.0, "sample_ms": 2.0, "sample_stride_ms": .5},
            "convergence": {"dt_ms": [.1, .05, .025]},
        }

    def test_reproducible(self):
        cfg = self.config()
        left = stationary_samples_nested(123, [3.0, .01, 4.0, .1], .5, cfg, .05)[0]
        right = stationary_samples_nested(123, [3.0, .01, 4.0, .1], .5, cfg, .05)[0]
        self.assertTrue(np.array_equal(left, right))
        self.assertEqual(left.shape, (4, 3))

    def test_all_steps_use_same_sample_clock(self):
        cfg = self.config()
        for dt in (.1, .05, .025):
            samples, _, ok = stationary_samples_nested(7, [3.0, .01, 4.0, .1], .5, cfg, dt)
            self.assertTrue(ok)
            self.assertEqual(samples.shape, (4, 3))


if __name__ == "__main__":
    unittest.main()
