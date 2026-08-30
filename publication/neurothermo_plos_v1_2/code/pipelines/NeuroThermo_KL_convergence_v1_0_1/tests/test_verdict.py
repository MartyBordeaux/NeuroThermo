import unittest

import pandas as pd

from kl_convergence.markers import MARKER_VARIANTS
from kl_convergence.verdict import build_verdict


class VerdictTierTests(unittest.TestCase):
    def config(self):
        return {
            "convergence": {"dt_ms": [0.05, 0.025, 0.0125],
                            "primary_dt_ms": 0.025, "seeds": [1, 2, 3]},
            "validation": {"expected_scenarios": 1},
            "verdict": {"enforce": True},
            "gates": {
                "minimum_grid_retention": 1.0,
                "minimum_scenario_convergence_fraction": 0.8,
                "maximum_fine_median_pair_delta": -0.05,
                "minimum_fine_pair_fraction_negative": 0.6,
                "maximum_fine_median_cell_delta": -0.05,
                "minimum_fine_cell_fraction_negative": 0.75,
                "minimum_fine_animal_pair_fraction_negative": 0.625,
                "ensemble_max_shift_coarse_primary": 0.05,
                "ensemble_max_shift_primary_fine": 0.05,
                "maximum_aggregation_order_disagreement": 0.075,
            },
        }

    def frames(self):
        rows = []
        for dt in (0.05, 0.025, 0.0125):
            for view in ("xyz", "xy"):
                for variant in MARKER_VARIANTS:
                    for order in ("marker_first", "curve_first"):
                        rows.append({
                            "dt_ms": dt, "view": view, "marker_variant": variant,
                            "aggregation_order": order, "median_pair_delta": -0.2,
                            "fraction_pairs_negative": 0.9, "median_cell_delta": -0.2,
                            "fraction_cells_negative": 0.9, "q75_cell_delta": -0.1,
                            "fraction_animal_pairs_negative": 0.9,
                            "leave_one_animal_out_all_negative": True,
                        })
        scenario = pd.DataFrame({"view": ["xyz", "xy"], "pass": [True, True]})
        retention = pd.DataFrame({"retained_min": [1.0] * 9})
        return scenario, pd.DataFrame(rows), retention

    def test_all_gates_pass(self):
        scenario, ensemble, retention = self.frames()
        verdict = build_verdict(self.config(), scenario, ensemble, None, retention)
        self.assertEqual(verdict["decision"], "KEEP_AS_MAIN_RESULT")

    def test_supporting_failure_downgrades(self):
        scenario, ensemble, retention = self.frames()
        scenario.loc[scenario.view == "xy", "pass"] = False
        verdict = build_verdict(self.config(), scenario, ensemble, None, retention)
        self.assertEqual(verdict["decision"], "KEEP_AS_ENSEMBLE_RESULT_WITH_LIMITATIONS")
        self.assertTrue(verdict["fatal_gates_pass"])

    def test_primary_failure_removes(self):
        scenario, ensemble, retention = self.frames()
        mask = ((ensemble.dt_ms == 0.0125) & (ensemble.view == "xyz") &
                (ensemble.marker_variant == "seed_median_curve_isotonic") &
                (ensemble.aggregation_order == "marker_first"))
        ensemble.loc[mask, "median_pair_delta"] = 0.1
        verdict = build_verdict(self.config(), scenario, ensemble, None, retention)
        self.assertEqual(verdict["decision"], "REMOVE_KL_RESULT")
        self.assertFalse(verdict["fatal_gates_pass"])


if __name__ == "__main__":
    unittest.main()
