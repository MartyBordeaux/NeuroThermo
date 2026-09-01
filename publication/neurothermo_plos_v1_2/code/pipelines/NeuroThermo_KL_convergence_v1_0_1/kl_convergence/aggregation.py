from __future__ import annotations

import numpy as np
import pandas as pd

from .markers import weighted_quantile


def pair_cells(pair_key):
    left, right = str(pair_key).split("__TO__", 1)
    return left, right


def weighted_summary(values, weights):
    return {
        "median": weighted_quantile(values, weights, 0.5),
        "q25": weighted_quantile(values, weights, 0.25),
        "q75": weighted_quantile(values, weights, 0.75),
    }


def cell_balanced(pair_frame):
    rows = []
    for _, row in pair_frame.iterrows():
        wt, sca = pair_cells(row.biological_pair_key)
        for genotype, cell in (("WT", wt), ("SCA3", sca)):
            rows.append({**row.to_dict(), "endpoint_genotype": genotype, "endpoint_cell": cell})
    expanded = pd.DataFrame(rows)
    keys = ["dt_ms", "view", "marker_variant", "aggregation_order", "endpoint_genotype", "endpoint_cell"]
    output = []
    for key, group in expanded.groupby(keys, sort=True):
        output.append(dict(zip(keys, key), n_partner_pairs=len(group),
                           kl_minus_firing_p=float(group.kl_minus_firing_p.median())))
    return pd.DataFrame(output)


def attach_animals(pair_frame, mapping):
    animal = dict(zip(mapping.cell_id.astype(str), mapping.animal_id.astype(str)))
    output = pair_frame.copy()
    cells = output.biological_pair_key.map(pair_cells)
    output["wt_cell_id"] = [item[0] for item in cells]
    output["sca3_cell_id"] = [item[1] for item in cells]
    output["wt_animal_id"] = output.wt_cell_id.map(animal)
    output["sca3_animal_id"] = output.sca3_cell_id.map(animal)
    if output[["wt_animal_id", "sca3_animal_id"]].isna().any().any():
        missing = output.loc[output.wt_animal_id.isna() | output.sca3_animal_id.isna(), "biological_pair_key"].unique()
        raise ValueError("Missing animal mapping for: " + ", ".join(map(str, missing)))
    output["animal_pair_key"] = output.wt_animal_id + "__TO__" + output.sca3_animal_id
    return output


def animal_pair_balanced(pair_frame, mapping):
    frame = attach_animals(pair_frame, mapping)
    keys = ["dt_ms", "view", "marker_variant", "aggregation_order", "animal_pair_key", "wt_animal_id", "sca3_animal_id"]
    rows = []
    for key, group in frame.groupby(keys, sort=True):
        rows.append(dict(zip(keys, key), n_cell_pairs=len(group),
                         kl_minus_firing_p=float(group.kl_minus_firing_p.median())))
    return pd.DataFrame(rows)


def leave_one_animal_out(animal_pairs):
    rows = []
    base_keys = ["dt_ms", "view", "marker_variant", "aggregation_order"]
    for key, group in animal_pairs.groupby(base_keys, sort=True):
        animals = sorted(set(group.wt_animal_id) | set(group.sca3_animal_id))
        for omitted in animals:
            kept = group[(group.wt_animal_id != omitted) & (group.sca3_animal_id != omitted)]
            rows.append(dict(zip(base_keys, key), omitted_animal_id=omitted,
                             n_animal_pairs=len(kept),
                             median_kl_minus_firing_p=float(kept.kl_minus_firing_p.median())))
    return pd.DataFrame(rows)
