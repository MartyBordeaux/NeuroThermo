from pathlib import Path
import unittest

import pandas as pd

from nonequilibrium_geometry.aggregation import annotate_animal_pairs, balanced_mean


ROOT = Path(__file__).resolve().parents[1]


class AggregationTests(unittest.TestCase):
    def test_known_core_cells_map_to_animals(self):
        mapping = pd.read_csv(ROOT / "data" / "animal_mapping.csv")
        frame = pd.DataFrame({"biological_pair_key": ["WT_03__TO__SCA3_01", "WT_09__TO__SCA3_06"]})
        annotated = annotate_animal_pairs(frame, mapping)
        self.assertEqual(annotated.loc[0, "animal_pair_key"], "WT_150130__TO__SCA3_DD20")
        self.assertEqual(annotated.loc[1, "animal_pair_key"], "WT_DD09__TO__SCA3_DD24")

    def test_balanced_mean_prevents_cell_count_weighting(self):
        frame = pd.DataFrame({
            "animal_pair_key": ["a", "a", "b"],
            "seed": [1, 1, 1],
            "p": [0.5, 0.5, 0.5],
            "value": [0.0, 2.0, 10.0],
        })
        animal = balanced_mean(frame, ["animal_pair_key", "seed", "p"])
        self.assertEqual(animal.loc[animal["animal_pair_key"].eq("a"), "value"].item(), 1.0)
        self.assertEqual(animal.loc[animal["animal_pair_key"].eq("b"), "value"].item(), 10.0)

    def test_core_secure_mapping_has_expected_animal_counts(self):
        mapping = pd.read_csv(ROOT / "data" / "animal_mapping.csv")
        wt = ["WT_03", "WT_04", "WT_08", "WT_09", "WT_10", "WT_13", "WT_14", "WT_16"]
        sca3 = ["SCA3_01", "SCA3_04", "SCA3_06", "SCA3_09"]
        frame = pd.DataFrame({
            "biological_pair_key": [f"{left}__TO__{right}" for left in wt for right in sca3]
        })
        annotated = annotate_animal_pairs(frame, mapping)
        self.assertEqual(annotated["wt_animal_id"].nunique(), 4)
        self.assertEqual(annotated["sca3_animal_id"].nunique(), 2)
        self.assertEqual(annotated["animal_pair_key"].nunique(), 8)
