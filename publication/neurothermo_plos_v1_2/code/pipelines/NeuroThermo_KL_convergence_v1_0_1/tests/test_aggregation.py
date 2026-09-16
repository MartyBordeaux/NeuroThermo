import unittest

import pandas as pd

from kl_convergence.aggregation import animal_pair_balanced, cell_balanced, leave_one_animal_out


class AggregationTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame([
            {"biological_pair_key": "WT_1__TO__SCA3_1", "dt_ms": .025, "view": "xyz",
             "marker_variant": "m", "aggregation_order": "marker_first", "kl_minus_firing_p": -.2},
            {"biological_pair_key": "WT_2__TO__SCA3_1", "dt_ms": .025, "view": "xyz",
             "marker_variant": "m", "aggregation_order": "marker_first", "kl_minus_firing_p": -.1},
        ])

    def test_cell_expansion(self):
        cells = cell_balanced(self.frame())
        self.assertEqual(set(cells.endpoint_cell), {"WT_1", "WT_2", "SCA3_1"})

    def test_animal_and_leave_one_out(self):
        mapping = pd.DataFrame({"cell_id": ["WT_1", "WT_2", "SCA3_1"],
                                "animal_id": ["WA", "WB", "SA"]})
        animals = animal_pair_balanced(self.frame(), mapping)
        self.assertEqual(len(animals), 2)
        loo = leave_one_animal_out(animals)
        self.assertEqual(set(loo.omitted_animal_id), {"WA", "WB", "SA"})


if __name__ == "__main__":
    unittest.main()
