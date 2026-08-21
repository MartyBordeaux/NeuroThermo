from pathlib import Path

import numpy as np
import pandas as pd

from neurothermo_phenotypes.config import load_config
from neurothermo_phenotypes.disease_coordinate import (
    build_disease_coordinate,
    plot_disease_coordinate,
    write_disease_coordinate_report,
)


def test_disease_coordinate_retains_all_cells_without_imputation(tmp_path: Path):
    cfg = load_config(None)
    cfg["statistics"]["bootstrap_iterations"] = 100
    cfg["disease_coordinate"].update({
        "currents_pA": [400.0, 500.0],
        "stability_bootstrap_iterations": 50,
        "exact_max_labelings": 200000,
    })
    scalar_rows = []
    feature_rows = []
    for group in ["WT", "SCA3"]:
        for cell_index in range(4):
            cell_id = f"{group}_{cell_index}"
            capacitance = (
                300.0 + 10.0 * cell_index if group == "WT"
                else 80.0 + 5.0 * cell_index
            )
            scalar_rows.append({
                "group": group, "cell_id": cell_id,
                "capacitance_20ms_pF": capacitance,
            })
            for current in [400.0, 500.0]:
                missing = group == "WT" and cell_index == 0
                feature_rows.append({
                    "group": group, "cell_id": cell_id,
                    "animal_id": "NA_NOT_RECOVERABLE", "current_pA": current,
                    "qc_pass": not missing, "thermo_eligible": not missing,
                    "mean_isi_ms": np.nan if missing else (
                        12.0 + cell_index if group == "WT" else 30.0 + cell_index
                    ),
                    "predictive_information_nats": np.nan if missing else (
                        1.9 - 0.01 * cell_index if group == "WT" else 1.7 - 0.01 * cell_index
                    ),
                })
    scores, reference, validation, long_scores = build_disease_coordinate(
        pd.DataFrame(feature_rows), pd.DataFrame(scalar_rows), cfg
    )
    assert len(scores) == 8
    assert scores["disease_burden_z"].notna().all()
    wt_missing = scores.loc[scores["cell_id"] == "WT_0"].iloc[0]
    assert wt_missing["domains_available"] == 1
    assert wt_missing["dynamic_values_available"] == 0
    assert wt_missing["evidence_grade"] == "structural_only"
    assert wt_missing["coordinate_reliability"] == "structural_only"
    complete = scores.loc[scores["cell_id"] == "SCA3_0"].iloc[0]
    assert complete["dynamic_values_available"] == 4
    assert complete["dynamic_values_expected"] == 4
    assert complete["evidence_grade"] == "full_dynamic"
    assert scores["crossfit_disease_burden_z"].notna().all()
    assert validation.loc[0, "internal_cross_fitted_auc_SCA3_vs_WT"] == 1.0
    assert validation.loc[0, "crossfit_min_SCA3_minus_max_WT_margin"] > 0
    assert scores["stability_bootstrap_iterations_valid"].eq(50).all()
    for column in [
        "bootstrap_p_outside_WT_robust_boundary",
        "bootstrap_p_outside_observed_WT_envelope",
        "bootstrap_p_WT_exit_consensus",
    ]:
        assert scores[column].between(0.0, 1.0).all()
    assert validation.loc[0, "permutation_mode"] == "exact"
    assert validation.loc[0, "valid_labelings"] == 70
    assert validation.loc[0, "descriptive_auc_SCA3_vs_WT"] == 1.0
    assert set(reference["level"]) == {"feature", "domain", "composite"}
    assert len(long_scores) == 8 * 5
    assert not long_scores.loc[
        (long_scores["cell_id"] == "WT_0")
        & (long_scores["domain"] != "structure"), "available"
    ].any()
    plot_path = tmp_path / "coordinate.png"
    plot_disease_coordinate(scores, plot_path)
    assert plot_path.exists() and plot_path.stat().st_size > 0
    report_path = tmp_path / "README.md"
    write_disease_coordinate_report(scores, validation, report_path)
    assert "q is not a probability" in report_path.read_text(encoding="utf-8")
    assert "leave-one-WT-out AUC" in report_path.read_text(encoding="utf-8")
