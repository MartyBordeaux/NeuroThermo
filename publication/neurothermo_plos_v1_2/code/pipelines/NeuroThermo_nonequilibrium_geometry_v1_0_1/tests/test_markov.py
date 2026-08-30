import numpy as np
import unittest

from nonequilibrium_geometry.markov import detailed_balance_metrics, exact_hatano_sasa


CFG = {"markov": {"cycle_min_flux": 1e-12, "cycle_min_occupancy": 0.0, "cycle_max_states": 10}}


class MarkovTests(unittest.TestCase):
    def test_reversible_chain_passes_detailed_balance(self):
        matrix = np.array([[0.8, 0.2], [0.2, 0.8]])
        pi = np.array([0.5, 0.5])
        metrics, _ = detailed_balance_metrics(matrix, pi, CFG)
        self.assertLess(metrics["markov_db_violation"], 1e-14)
        self.assertLess(abs(metrics["markov_entropy_per_lag"]), 1e-14)

    def test_exact_hatano_sasa_identity(self):
        matrices, stationary = [], []
        for q in (0.25, 0.5, 0.75):
            pi = np.array([q, 1 - q])
            matrices.append(np.tile(pi, (2, 1)))
            stationary.append(pi)
        self.assertTrue(np.isclose(exact_hatano_sasa(matrices, stationary), 1.0, atol=1e-12))
