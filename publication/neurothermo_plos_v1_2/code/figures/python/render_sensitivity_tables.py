"""Generate reviewer-requested descriptive sensitivity tables.

The script uses only the frozen publication CSV files. It does not refit the
Hindmarsh--Rose model or rerun stochastic simulations.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "figure_source"


def partial_spearman(x, y, controls):
    """Rank-residualized partial Spearman coefficient."""
    design = np.column_stack(
        [np.ones(len(x))] + [rankdata(np.asarray(c)) for c in controls]
    )
    rank_x = rankdata(np.asarray(x))
    rank_y = rankdata(np.asarray(y))
    residual_x = rank_x - design @ np.linalg.lstsq(
        design, rank_x, rcond=None
    )[0]
    residual_y = rank_y - design @ np.linalg.lstsq(
        design, rank_y, rcond=None
    )[0]
    return float(np.corrcoef(residual_x, residual_y)[0, 1])


def main():
    endpoint = pd.read_csv(DATA / "fig1_endpoint_cells.csv")
    characterization = pd.read_csv(
        DATA / "supp_characterization_primary_cells.csv"
    )
    endpoint = endpoint.merge(
        characterization[
            ["cell_id", "rheobase_pA", "rheobase_J", "median_circle_f1"]
        ],
        on="cell_id",
        how="left",
    )

    resolved = endpoint[endpoint["animal_id"] != "NA_NOT_RECOVERABLE"]
    variables = [
        "capacitance_pF",
        "rheobase_pA",
        "rheobase_J",
        "exp_q75_firing_rate_hz",
        "exp_q75_mean_isi_ms",
        "kappa_I",
        "cell_loss",
        "median_circle_f1",
    ]
    animal = (
        resolved.groupby(["group", "animal_id"])
        .agg(
            n_cells=("cell_id", "size"),
            **{variable: (variable, "median") for variable in variables},
        )
        .reset_index()
    )
    animal.to_csv(DATA / "supp_animal_medians.csv", index=False)

    genotype = (endpoint["group"] == "SCA3").astype(float).to_numpy()
    capacitance = endpoint["capacitance_pF"].to_numpy(float)
    kappa = endpoint["kappa_I"].to_numpy(float)
    associations = pd.DataFrame(
        [
            {
                "association": "kappa_I_vs_capacitance",
                "control": "none",
                "coefficient": spearmanr(kappa, capacitance).statistic,
                "method": "Spearman rank correlation",
                "n_cells": len(endpoint),
            },
            {
                "association": "kappa_I_vs_capacitance",
                "control": "genotype",
                "coefficient": partial_spearman(
                    kappa, capacitance, [genotype]
                ),
                "method": "partial Spearman by rank residualization",
                "n_cells": len(endpoint),
            },
            {
                "association": "kappa_I_vs_genotype_SCA3",
                "control": "capacitance",
                "coefficient": partial_spearman(kappa, genotype, [capacitance]),
                "method": "partial Spearman by rank residualization",
                "n_cells": len(endpoint),
            },
        ]
    )
    associations.to_csv(
        DATA / "supp_kappa_capacitance_association.csv", index=False
    )

    pairs = pd.read_csv(DATA / "fig4_marker_alignment_by_pair.csv")
    pairs[["wt_cell", "sca3_cell"]] = pairs[
        "biological_pair_key"
    ].str.split("__TO__", expand=True)
    cell_tables = []
    for role, column in (("WT", "wt_cell"), ("SCA3", "sca3_cell")):
        cell_table = (
            pairs.groupby(column)
            .agg(
                n_pair_combinations=("biological_pair_key", "size"),
                median_kl_balance_p=("kl_balance_p_median", "median"),
                median_dynamical_balance_p=(
                    "balance_p_isi_weighted_median",
                    "median",
                ),
            )
            .reset_index()
            .rename(columns={column: "cell_id"})
        )
        cell_table.insert(0, "endpoint_role", role)
        cell_table["kl_at_lower_p_than_dynamical_balance"] = (
            cell_table["median_kl_balance_p"]
            < cell_table["median_dynamical_balance_p"]
        )
        cell_tables.append(cell_table)
    pd.concat(cell_tables, ignore_index=True).to_csv(
        DATA / "supp_cell_balanced_marker_alignment.csv", index=False
    )


if __name__ == "__main__":
    main()
