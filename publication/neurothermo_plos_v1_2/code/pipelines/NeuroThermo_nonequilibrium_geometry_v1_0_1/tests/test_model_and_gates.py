from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import unittest

from nonequilibrium_geometry.model import drift_vec, fixed_points
from nonequilibrium_geometry.verdict import formalism_verdict, validate_physical_mapping


ROOT = Path(__file__).resolve().parents[1]


def config():
    return yaml.safe_load((ROOT / "configs" / "smoke_nonequilibrium_geometry_v1_0_1.yaml").read_text())


class ModelAndGateTests(unittest.TestCase):
    def test_fixed_points_satisfy_drift(self):
        cfg = config()
        theta = [3.0, 0.01, 4.0, 0.1]
        for x, y, z, _ in fixed_points(theta, 0.5, cfg):
            residual = drift_vec(x, y, z, 0.5, *theta, 1.0, 1.0, 5.0, -1.6)
            self.assertLess(np.linalg.norm(residual), 1e-8)

    def test_incomplete_physical_mapping_is_rejected(self):
        cfg = config()
        cfg["physical_mapping"] = {"enabled": True, "beta": 1.0, "energy_definition": "candidate", "work_controls": ["J"]}
        with self.assertRaisesRegex(ValueError, "cannot omit morphed controls"):
            validate_physical_mapping(cfg)

    def test_markov_gate_can_classify_ness_when_continuous_qc_fails(self):
        cfg = config()
        geometry = pd.DataFrame({
            "circulation_fraction": [1.0, 1.0],
            "stationary_current_divergence_relative": [0.8, 0.9],
            "markov_db_violation": [0.4, 0.5],
        })
        verdict = formalism_verdict(geometry, cfg)
        self.assertEqual(verdict["stationary_formalism"], "NESS")
        self.assertEqual(verdict["continuous_current_status"], "DIAGNOSTIC_INVALID")
        self.assertFalse(verdict["continuous_current_used_for_formalism"])
