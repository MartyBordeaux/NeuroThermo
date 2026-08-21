from __future__ import annotations

import copy
import gc
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from neurothermo_vulnerability import pipeline as engine

from . import __version__
from .config import REQUIRED_UPSTREAM, REQUIRED_V061


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_dump(path: Path, payload: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)


def _read_sweeps(upstream: Path) -> pd.DataFrame:
    frame = pd.read_csv(upstream / "sweep_features.csv")
    if "animal_id" in frame.columns:
        frame = frame.drop(columns=["animal_id"])
    return frame


def validate_inputs(config: Mapping[str, Any], upstream: Path, v061: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    missing_upstream = [str(upstream / x) for x in REQUIRED_UPSTREAM if not (upstream / x).is_file()]
    missing_v061 = [str(v061 / x) for x in REQUIRED_V061 if not (v061 / x).is_file()]
    rows.append({"check": "required_v031_files", "passed": not missing_upstream, "detail": "none" if not missing_upstream else "; ".join(missing_upstream)})
    rows.append({"check": "required_v061_files", "passed": not missing_v061, "detail": "none" if not missing_v061 else "; ".join(missing_v061)})
    if missing_upstream or missing_v061:
        return pd.DataFrame(rows)

    sweeps = _read_sweeps(upstream)
    currents = [int(x) for x in config["analysis"]["currents_pA"]]
    groups = [str(x) for x in config["analysis"]["groups"]]
    covariate_sources = {
        str(spec.get("source", name))
        for name, spec in config["covariates"].items()
        if str(spec.get("source", name)) != "mean_isi_missing"
    }
    required = {
        "group", "cell_id", "current_pA", "predictive_information_nats",
        "firing_rate_hz", "mean_isi_ms",
    } | covariate_sources
    absent = sorted(required - set(sweeps.columns))
    rows.append({"check": "required_columns", "passed": not absent, "detail": "none" if not absent else ", ".join(absent)})
    if absent:
        return pd.DataFrame(rows)

    domain = sweeps[sweeps.current_pA.isin(currents) & sweeps.group.isin(groups)].copy()
    duplicates = int(domain.duplicated(["cell_id", "current_pA"]).sum())
    counts = domain[["group", "cell_id"]].drop_duplicates().groupby("group").size().to_dict()
    expected = config["analysis"]["expected_cells"]
    rows.append({"check": "one_row_per_cell_current", "passed": duplicates == 0, "detail": "duplicates={}".format(duplicates)})
    rows.append({"check": "frozen_cell_counts", "passed": all(int(counts.get(g, 0)) == int(expected[g]) for g in groups), "detail": json.dumps(counts, sort_keys=True)})
    per_cell = domain.groupby("cell_id").current_pA.nunique()
    bad_grid = per_cell[per_cell != len(currents)].to_dict()
    rows.append({"check": "complete_current_grid", "passed": not bad_grid, "detail": "none" if not bad_grid else json.dumps(bad_grid, sort_keys=True)})
    one_group = domain.groupby("cell_id").group.nunique()
    rows.append({"check": "one_group_per_cell", "passed": bool((one_group == 1).all()), "detail": "max={}".format(int(one_group.max()))})
    positive_pi = pd.to_numeric(domain.predictive_information_nats, errors="coerce")
    rows.append({"check": "positive_finite_predictive_information", "passed": bool(np.isfinite(positive_pi).all() and (positive_pi > 0).all()), "detail": "finite={}/{}; min={:.6g}".format(int(np.isfinite(positive_pi).sum()), len(positive_pi), float(positive_pi.min()))})

    vsummary = pd.read_csv(v061 / "cell_vulnerability_summary.csv")
    upstream_cells = set(domain.cell_id.astype(str).unique())
    v061_cells = set(vsummary.cell_id.astype(str).unique())
    rows.append({"check": "v061_cell_identity", "passed": upstream_cells == v061_cells, "detail": "upstream={}; v061={}; symmetric_difference={}".format(len(upstream_cells), len(v061_cells), sorted(upstream_cells ^ v061_cells))})
    vmanifest = json.loads((v061 / "analysis_manifest.json").read_text(encoding="utf-8"))
    rows.append({"check": "v061_version", "passed": str(vmanifest.get("pipeline_version")) == "0.6.1", "detail": str(vmanifest.get("pipeline_version"))})
    rows.append({"check": "cell_level_only", "passed": "animal_id" not in domain.columns and not bool(vmanifest.get("animal_level_inference", True)), "detail": "animal_id_removed=True; v061_animal_level_inference={}".format(vmanifest.get("animal_level_inference"))})
    return pd.DataFrame(rows)


def _transform(values: pd.Series, name: str) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(x.shape, np.nan, dtype=float)
    finite = np.isfinite(x)
    if name == "identity":
        out[finite] = x[finite]
    elif name == "log":
        valid = finite & (x > 0)
        out[valid] = np.log(x[valid])
    elif name == "log1p":
        valid = finite & (x >= 0)
        out[valid] = np.log1p(x[valid])
    else:
        raise ValueError("Unknown covariate transform: {}".format(name))
    return out


def _prepare_adjustment_table(config: Mapping[str, Any], upstream: Path, cell_ids: Sequence[str], currents: Sequence[int]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sweeps = _read_sweeps(upstream)
    table = sweeps[sweeps.cell_id.astype(str).isin(cell_ids) & sweeps.current_pA.isin(currents)].copy()
    table["cell_id"] = table.cell_id.astype(str)
    order = pd.MultiIndex.from_product([cell_ids, currents], names=["cell_id", "current_pA"])
    table = table.set_index(["cell_id", "current_pA"]).reindex(order).reset_index()
    table["log_predictive_information"] = np.log(pd.to_numeric(table.predictive_information_nats, errors="coerce"))
    table["mean_isi_missing"] = table.mean_isi_ms.isna().astype(float)
    variation_rows = []
    for name, spec in config["covariates"].items():
        source = str(spec.get("source", name))
        if source == "mean_isi_missing":
            values = table["mean_isi_missing"].to_numpy(float)
        else:
            values = _transform(table[source], str(spec["transform"]))
        table[name] = values
        finite = values[np.isfinite(values)]
        variation_rows.append({
            "covariate": name,
            "source_column": source,
            "transform": str(spec["transform"]),
            "n_finite": int(len(finite)),
            "n_unique": int(len(np.unique(finite))) if len(finite) else 0,
            "minimum": float(np.min(finite)) if len(finite) else np.nan,
            "maximum": float(np.max(finite)) if len(finite) else np.nan,
            "standard_deviation": float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan,
            "used_in_any_adjustment": any(name in values for values in config["adjustment_modes"].values()),
        })
    keep = ["group", "cell_id", "current_pA", "predictive_information_nats", "log_predictive_information"] + list(config["covariates"])
    return table[keep], pd.DataFrame(variation_rows)


def _ridge_crossfit_by_cell_current(table: pd.DataFrame, cell_ids: Sequence[str], currents: Sequence[int], covariates: Sequence[str], ridge_lambda: float):
    prediction = np.full(len(table), np.nan, dtype=float)
    coefficients: List[Dict[str, Any]] = []
    for current in currents:
        current_mask = table.current_pA.to_numpy(int) == int(current)
        for cell in cell_ids:
            target_mask = current_mask & (table.cell_id.to_numpy(str) == str(cell))
            train_mask = current_mask & ~target_mask
            target_indices = np.flatnonzero(target_mask)
            if len(target_indices) != 1:
                raise ValueError("Expected one row for {} at {} pA".format(cell, current))
            train = table.loc[train_mask]
            test = table.loc[target_mask]
            y = train.log_predictive_information.to_numpy(float)
            if not np.isfinite(y).all():
                raise ValueError("Non-finite predictive information in residualization training data.")
            train_columns = []
            test_columns = []
            active_names = []
            metadata = []
            for name in covariates:
                x_train = train[name].to_numpy(float)
                x_test = test[name].to_numpy(float)
                finite = np.isfinite(x_train)
                median = float(np.nanmedian(x_train)) if finite.any() else 0.0
                x_train = np.where(finite, x_train, median)
                x_test = np.where(np.isfinite(x_test), x_test, median)
                center = float(np.mean(x_train))
                scale = float(np.std(x_train, ddof=1)) if len(x_train) > 1 else 0.0
                active = bool(np.isfinite(scale) and scale > 1e-12)
                metadata.append((name, median, center, scale, active))
                if active:
                    train_columns.append((x_train - center) / scale)
                    test_columns.append((x_test - center) / scale)
                    active_names.append(name)
            x_train_matrix = np.column_stack([np.ones(len(train))] + train_columns)
            x_test_matrix = np.column_stack([np.ones(len(test))] + test_columns)
            penalty = np.eye(x_train_matrix.shape[1], dtype=float) * float(ridge_lambda)
            penalty[0, 0] = 0.0
            beta = np.linalg.solve(x_train_matrix.T @ x_train_matrix + penalty, x_train_matrix.T @ y)
            prediction[target_indices] = x_test_matrix @ beta
            coefficients.append({
                "target_cell_id": str(cell), "current_pA": int(current), "term": "intercept",
                "standardized_coefficient": float(beta[0]), "active": True,
                "training_n_cells": int(len(train)), "ridge_lambda": float(ridge_lambda),
            })
            coefficient_by_name = dict(zip(active_names, beta[1:]))
            for name, median, center, scale, active in metadata:
                coefficients.append({
                    "target_cell_id": str(cell), "current_pA": int(current), "term": name,
                    "standardized_coefficient": float(coefficient_by_name[name]) if active else 0.0,
                    "active": active, "training_n_cells": int(len(train)),
                    "ridge_lambda": float(ridge_lambda), "training_imputation_median": median,
                    "training_center": center, "training_scale": scale,
                })
    residual = table.log_predictive_information.to_numpy(float) - prediction
    return prediction, residual, pd.DataFrame(coefficients)


def _residualize(config: Mapping[str, Any], table: pd.DataFrame, cell_ids: Sequence[str], currents: Sequence[int]):
    primary_lambda = float(config["residualization"]["ridge_lambda"])
    outputs = table.copy()
    coefficient_tables = []
    mode_residuals: Dict[str, np.ndarray] = {}
    for mode, covariates in config["adjustment_modes"].items():
        pred, residual, coefficients = _ridge_crossfit_by_cell_current(table, cell_ids, currents, covariates, primary_lambda)
        outputs[mode + "__prediction"] = pred
        outputs[mode + "__residual"] = residual
        coefficients.insert(0, "adjustment_mode", mode)
        coefficient_tables.append(coefficients)
        mode_residuals[mode] = residual
    lambda_residuals: Dict[float, np.ndarray] = {primary_lambda: mode_residuals["activity_technical"]}
    for value in config["residualization"]["ridge_lambda_sensitivity"]:
        value = float(value)
        if np.isclose(value, primary_lambda):
            continue
        _, residual, _ = _ridge_crossfit_by_cell_current(
            table, cell_ids, currents, config["adjustment_modes"]["activity_technical"], value
        )
        lambda_residuals[value] = residual
    return outputs, pd.concat(coefficient_tables, ignore_index=True), mode_residuals, lambda_residuals


def _safe_spearman(first: pd.Series, second: pd.Series) -> float:
    pair = pd.concat([first, second], axis=1).dropna()
    if len(pair) < 4 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
        return np.nan
    return float(pair.iloc[:, 0].corr(pair.iloc[:, 1], method="spearman"))


def _residual_diagnostics(config: Mapping[str, Any], table: pd.DataFrame):
    performance = []
    pooled = []
    currentwise = []
    y = table.log_predictive_information
    for mode, covariates in config["adjustment_modes"].items():
        pred = table[mode + "__prediction"]
        residual = table[mode + "__residual"]
        ss_res = float(np.sum(np.square(y - pred)))
        ss_tot = float(np.sum(np.square(y - y.mean())))
        performance.append({
            "adjustment_mode": mode, "n": int(len(table)),
            "crossfit_MAE_log_PI": float(np.mean(np.abs(residual))),
            "crossfit_RMSE_log_PI": float(np.sqrt(np.mean(np.square(residual)))),
            "crossfit_R2_log_PI": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
            "spearman_observed_vs_prediction": _safe_spearman(y, pred),
            "residual_median": float(np.median(residual)),
        })
        for covariate in covariates:
            before = _safe_spearman(table[covariate], y)
            after = _safe_spearman(table[covariate], residual)
            pooled.append({
                "adjustment_mode": mode, "covariate": covariate, "n": int(table[[covariate, "log_predictive_information"]].dropna().shape[0]),
                "spearman_raw_log_PI": before, "spearman_residual_log_PI": after,
                "absolute_dependence_reduction": abs(before) - abs(after) if np.isfinite(before) and np.isfinite(after) else np.nan,
            })
            for current, sub in table.groupby("current_pA", sort=True):
                currentwise.append({
                    "adjustment_mode": mode, "covariate": covariate, "current_pA": int(current),
                    "n": int(sub[[covariate, "log_predictive_information"]].dropna().shape[0]),
                    "spearman_raw_log_PI": _safe_spearman(sub[covariate], sub.log_predictive_information),
                    "spearman_residual_log_PI": _safe_spearman(sub[covariate], sub[mode + "__residual"]),
                })
    return pd.DataFrame(performance), pd.DataFrame(pooled), pd.DataFrame(currentwise)


def _matrix_from_rows(table: pd.DataFrame, values: np.ndarray, cell_ids: Sequence[str], currents: Sequence[int]) -> np.ndarray:
    work = table[["cell_id", "current_pA"]].copy()
    work["value"] = values
    return work.pivot(index="cell_id", columns="current_pA", values="value").reindex(index=cell_ids, columns=currents).to_numpy(float)


def _mode_config(config: Mapping[str, Any]) -> Dict[str, Any]:
    mode_config = copy.deepcopy(dict(config))
    mode_config["primary_features"]["predictive_information_nats"]["transform"] = "identity"
    mode_config["primary_features"]["predictive_information_nats"]["direction"] = -1
    mode_config["primary_features"]["predictive_information_nats"]["fallback_scale"] = float(config["residualization"]["residual_fallback_scale"])
    return mode_config


def _run_exact_mode(config, mode, residual_matrix, base_matrices, cell_ids, labels, currents, masks, run_conditional_bootstrap=False):
    mode_config = _mode_config(config)
    matrices = dict(base_matrices)
    matrices["predictive_information_nats"] = residual_matrix
    cache = engine.SubsetReferenceCache(mode_config, matrices, int((~labels).sum()))
    burden, excitability, predictive, strict, observed_index, observed = engine._exact_target_specific_scores(mode_config, labels, currents, cache, masks)
    scores, references = engine._observed_score_tables(mode_config, cell_ids, labels, currents, cache, observed)
    scores = scores.rename(columns={"predictive_information_nats__oriented_z": "residual_predictive_information__oriented_z"})
    scores["residual_predictive_dynamics_z"] = scores["predictive_dynamics_z"]
    scores.insert(0, "adjustment_mode", mode)
    references.insert(0, "adjustment_mode", mode)
    references["feature"] = references.feature.replace({"predictive_information_nats": "residual_log_predictive_information"})
    cell_summary = engine._cell_summary(mode_config, scores, cell_ids, labels, currents)
    cell_summary.insert(0, "adjustment_mode", mode)

    current_tables = []
    curve_tables = []
    for analysis, values in (
        ("independent_domain_burden", burden),
        ("residual_predictive_dynamics", predictive),
    ):
        differences = engine._group_differences(values, masks)
        current_table, curve_table = engine._exact_difference_tables(currents, differences, observed_index, analysis, 1.0)
        current_table.insert(0, "adjustment_mode", mode)
        curve_table["adjustment_mode"] = mode
        current_tables.append(current_table)
        curve_tables.append(pd.DataFrame([curve_table]))
    i_exit, _ = engine._exact_i_exit(mode_config, currents, strict, masks, observed_index)
    i_exit.insert(0, "adjustment_mode", mode)
    i_exit["analysis"] = "calibrated_residual_predictive_I_exit"

    conditional_bootstrap = None
    sensitivity = None
    if run_conditional_bootstrap:
        conditional_bootstrap = engine._bootstrap(mode_config, labels, cell_ids, currents, cache)
        conditional_bootstrap.insert(0, "adjustment_mode", mode)
        conditional_bootstrap["residualizer_refit_in_bootstrap"] = False
        conditional_bootstrap["interpretation"] = "conditional on label-blind cross-fitted residuals"
        sensitivity = engine._sensitivity(mode_config, scores, currents)
        sensitivity.insert(0, "adjustment_mode", mode)

    result = {
        "scores": scores,
        "references": references,
        "cell_summary": cell_summary,
        "current_tests": pd.concat(current_tables, ignore_index=True),
        "curve_tests": pd.concat(curve_tables, ignore_index=True),
        "i_exit": i_exit,
        "conditional_bootstrap": conditional_bootstrap,
        "sensitivity": sensitivity,
    }
    del burden, excitability, predictive, strict, observed, cache
    gc.collect()
    return result


def _observed_only_exits(config, residual_matrix, base_matrices, cell_ids, labels, currents):
    mode_config = _mode_config(config)
    matrices = dict(base_matrices)
    matrices["predictive_information_nats"] = residual_matrix
    cache = engine.SubsetReferenceCache(mode_config, matrices, int((~labels).sum()))
    observed_mask = labels.reshape(1, -1)
    _, _, _, _, _, observed = engine._exact_target_specific_scores(mode_config, labels, currents, cache, observed_mask)
    scores, _ = engine._observed_score_tables(mode_config, cell_ids, labels, currents, cache, observed)
    summary = engine._cell_summary(mode_config, scores, cell_ids, labels, currents)
    return summary


def _cell_auc(values: np.ndarray, currents: Sequence[int]) -> np.ndarray:
    x = np.asarray(currents, dtype=float)
    out = np.full(values.shape[0], np.nan)
    for i, row in enumerate(values):
        valid = np.isfinite(row)
        if valid.sum() >= 2:
            out[i] = np.trapz(row[valid], x[valid]) / (x[valid][-1] - x[valid][0])
    return out


def _exact_cell_difference(values: np.ndarray, masks: np.ndarray, observed_labels: np.ndarray, expected_direction: float):
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    case_sum = masks.astype(float) @ np.where(valid, values, 0.0)
    case_n = masks.astype(float) @ valid.astype(float)
    control_masks = ~masks
    control_sum = control_masks.astype(float) @ np.where(valid, values, 0.0)
    control_n = control_masks.astype(float) @ valid.astype(float)
    differences = case_sum / case_n - control_sum / control_n
    observed_index = int(np.flatnonzero(np.all(masks == observed_labels[None, :], axis=1))[0])
    observed = float(differences[observed_index])
    oriented = expected_direction * differences
    observed_oriented = expected_direction * observed
    return observed, float(np.mean(oriented >= observed_oriented - 1e-12)), float(np.mean(np.abs(differences) >= abs(observed) - 1e-12))


def _ridge_sensitivity(config, table, lambda_residuals, base_matrices, cell_ids, labels, currents, masks):
    rows = []
    for value in sorted(lambda_residuals):
        matrix = _matrix_from_rows(table, lambda_residuals[value], cell_ids, currents)
        auc = _cell_auc(matrix, currents)
        observed, p_expected, p_two = _exact_cell_difference(auc, masks, labels, -1.0)
        summary = _observed_only_exits(config, matrix, base_matrices, cell_ids, labels, currents)
        exiting = summary[~summary.I_exit_calibrated_censored]
        rows.append({
            "ridge_lambda": float(value),
            "observed_residual_AUC_difference_SCA3_minus_WT": observed,
            "observed_expected_direction_difference_WT_minus_SCA3": -observed,
            "exact_p_residual_AUC_expected_direction": p_expected,
            "exact_p_residual_AUC_two_sided": p_two,
            "n_WT_exit_by_600": int(((exiting.group == config["analysis"]["groups"][0])).sum()),
            "n_SCA3_exit_by_600": int(((exiting.group == config["analysis"]["groups"][1])).sum()),
            "SCA3_exiting_cells": ";".join(exiting.loc[exiting.group == config["analysis"]["groups"][1], "cell_id"].astype(str)),
            "note": "exit counts are observed-label sensitivity; exact exit p-values are primary-lambda only",
        })
    return pd.DataFrame(rows)


def _group_curves(scores: pd.DataFrame, residual_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (mode, group, current), sub in scores.groupby(["adjustment_mode", "group", "current_pA"], sort=False):
        for metric in ["independent_domain_burden", "excitability_timing_z", "residual_predictive_dynamics_z"]:
            x = sub[metric].dropna()
            rows.append({"adjustment_mode": mode, "metric": metric, "score_mode": "target_specific_nested_LOO", "group": group, "current_pA": int(current), "n": int(len(x)), "median": float(x.median()), "q25": float(x.quantile(.25)), "q75": float(x.quantile(.75))})
    for mode in [c[:-10] for c in residual_table.columns if c.endswith("__residual")]:
        metric = mode + "__residual"
        for (group, current), sub in residual_table.groupby(["group", "current_pA"], sort=False):
            x = sub[metric].dropna()
            rows.append({"adjustment_mode": mode, "metric": "residual_log_predictive_information", "score_mode": "label_blind_leave_one_cell_out", "group": group, "current_pA": int(current), "n": int(len(x)), "median": float(x.median()), "q25": float(x.quantile(.25)), "q75": float(x.quantile(.75))})
    return pd.DataFrame(rows)


def _comparison(v061: Path, curve_tests: pd.DataFrame, i_exit_tests: pd.DataFrame) -> pd.DataFrame:
    rows = []
    old_curve = pd.read_csv(v061 / "group_curve_exact_tests.csv")
    old_burden = old_curve[old_curve.analysis == "independent_domain_burden"].iloc[0]
    old_exit = pd.read_csv(v061 / "I_exit_exact_test.csv").iloc[0]
    rows.append({
        "analysis_version": "v0.6.1", "adjustment_mode": "unadjusted_predictive_information",
        "burden_curve_AUC_expected_difference": float(old_burden.observed_expected_direction_curve_auc_difference),
        "burden_curve_exact_p_expected": float(old_burden.exact_p_curve_auc_expected_direction),
        "burden_maxT_exact_p_expected": float(old_burden.exact_p_omnibus_maxT_expected_direction),
        "n_WT_exit_by_600": int(old_exit.observed_n_WT_exit_by_600), "n_SCA3_exit_by_600": int(old_exit.observed_n_SCA3_exit_by_600),
        "I_exit_restricted_mean_difference_WT_minus_SCA3_pA": float(old_exit.observed_restricted_mean_difference_WT_minus_SCA3_pA),
        "I_exit_exact_p_expected": float(old_exit.exact_p_restricted_mean_expected_direction),
    })
    for mode in curve_tests.adjustment_mode.unique():
        burden = curve_tests[(curve_tests.adjustment_mode == mode) & (curve_tests.analysis == "independent_domain_burden")].iloc[0]
        exit_row = i_exit_tests[i_exit_tests.adjustment_mode == mode].iloc[0]
        rows.append({
            "analysis_version": "v0.7.0", "adjustment_mode": mode,
            "burden_curve_AUC_expected_difference": float(burden.observed_expected_direction_curve_auc_difference),
            "burden_curve_exact_p_expected": float(burden.exact_p_curve_auc_expected_direction),
            "burden_maxT_exact_p_expected": float(burden.exact_p_omnibus_maxT_expected_direction),
            "n_WT_exit_by_600": int(exit_row.observed_n_WT_exit_by_600), "n_SCA3_exit_by_600": int(exit_row.observed_n_SCA3_exit_by_600),
            "I_exit_restricted_mean_difference_WT_minus_SCA3_pA": float(exit_row.observed_restricted_mean_difference_WT_minus_SCA3_pA),
            "I_exit_exact_p_expected": float(exit_row.exact_p_restricted_mean_expected_direction),
        })
    return pd.DataFrame(rows)


def _exit_concordance(v061: Path, primary_summary: pd.DataFrame) -> pd.DataFrame:
    old = pd.read_csv(v061 / "cell_vulnerability_summary.csv")[["group", "cell_id", "I_exit_calibrated_pA", "I_exit_calibrated_censored"]]
    old = old.rename(columns={"I_exit_calibrated_pA": "I_exit_v061_pA", "I_exit_calibrated_censored": "I_exit_v061_censored"})
    new = primary_summary[["group", "cell_id", "I_exit_calibrated_pA", "I_exit_calibrated_censored"]].rename(columns={"I_exit_calibrated_pA": "I_exit_v070_primary_pA", "I_exit_calibrated_censored": "I_exit_v070_primary_censored"})
    merged = old.merge(new, on=["group", "cell_id"], how="outer", validate="one_to_one")
    merged["exit_status_concordant"] = merged.I_exit_v061_censored == merged.I_exit_v070_primary_censored
    merged["I_exit_shift_pA"] = merged.I_exit_v070_primary_pA - merged.I_exit_v061_pA
    return merged


def _plots(output: Path, config: Mapping[str, Any], residual_table: pd.DataFrame, curves: pd.DataFrame, scores: pd.DataFrame, comparison: pd.DataFrame, concordance: pd.DataFrame, pooled_dependence: pd.DataFrame, bootstrap: pd.DataFrame):
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    groups = config["analysis"]["groups"]
    colors = {groups[0]: "#2166ac", groups[1]: "#b2182b"}
    primary = "activity_technical"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for group in groups:
        sub = residual_table[residual_table.group == group]
        raw = sub.groupby("current_pA").log_predictive_information.agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
        raw.columns = ["median", "q25", "q75"]
        axes[0].plot(raw.index, raw["median"], marker="o", color=colors[group], label=group)
        axes[0].fill_between(raw.index.to_numpy(float), raw.q25.to_numpy(float), raw.q75.to_numpy(float), color=colors[group], alpha=.18)
        res = sub.groupby("current_pA")[primary + "__residual"].agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
        res.columns = ["median", "q25", "q75"]
        axes[1].plot(res.index, res["median"], marker="o", color=colors[group], label=group)
        axes[1].fill_between(res.index.to_numpy(float), res.q25.to_numpy(float), res.q75.to_numpy(float), color=colors[group], alpha=.18)
    axes[0].set(title="Unadjusted predictive information", ylabel="log predictive information", xlabel="Injected current (pA)")
    axes[1].axhline(0, color="black", ls="--", lw=1)
    axes[1].set(title="Activity + technical residual", ylabel="cross-fitted residual log predictive information", xlabel="Injected current (pA)")
    axes[0].legend(frameon=False)
    fig.tight_layout(); fig.savefig(figures / "predictive_information_before_after_adjustment.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for ax, metric, title in zip(axes, ["residual_predictive_dynamics_z", "independent_domain_burden"], ["Residual predictive-dynamics domain", "Adjusted two-domain burden"]):
        for group in groups:
            sub = curves[(curves.adjustment_mode == primary) & (curves.metric == metric) & (curves.group == group)].sort_values("current_pA")
            ax.plot(sub.current_pA, sub["median"], marker="o", color=colors[group], label=group)
            ax.fill_between(sub.current_pA.to_numpy(float), sub.q25.to_numpy(float), sub.q75.to_numpy(float), color=colors[group], alpha=.18)
        ax.axhline(float(config["analysis"]["robust_z_threshold"]), color="black", ls="--", lw=1)
        ax.set(title=title, xlabel="Injected current (pA)")
    axes[0].set_ylabel("WT-oriented robust z"); axes[1].legend(frameon=False)
    fig.tight_layout(); fig.savefig(figures / "adjusted_predictive_domain_and_burden.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    primary_scores = scores[scores.adjustment_mode == primary]
    for ax, group in zip(axes, groups):
        for cell, sub in primary_scores[primary_scores.group == group].groupby("cell_id"):
            sub = sub.sort_values("current_pA")
            ax.plot(sub.current_pA, sub.independent_domain_burden, marker="o", lw=1.2, alpha=.75, label=cell)
        ax.axhline(float(config["analysis"]["robust_z_threshold"]), color="black", ls="--", lw=1)
        ax.set(title=group, xlabel="Injected current (pA)")
    axes[0].set_ylabel("Adjusted independent-domain burden")
    fig.suptitle("Primary adjusted cell trajectories")
    fig.tight_layout(); fig.savefig(figures / "primary_adjusted_cell_trajectories.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for group in groups:
        sub = concordance[concordance.group == group]
        ax.scatter(sub.I_exit_v061_pA, sub.I_exit_v070_primary_pA, color=colors[group], s=50, label=group)
        for _, row in sub.iterrows():
            if row.I_exit_v061_pA != row.I_exit_v070_primary_pA:
                ax.annotate(row.cell_id, (row.I_exit_v061_pA, row.I_exit_v070_primary_pA), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.plot([400, 650], [400, 650], color="gray", ls="--", lw=1)
    ax.set(xlabel="v0.6.1 I_exit (pA)", ylabel="v0.7.0 primary adjusted I_exit (pA)", title="Exit concordance before and after adjustment", xlim=(425, 660), ylim=(425, 660))
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(figures / "I_exit_v061_v070_concordance.png", dpi=180); plt.close(fig)

    show = pooled_dependence[pooled_dependence.adjustment_mode == primary].copy()
    x = np.arange(len(show)); width = .35
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(x - width/2, show.spearman_raw_log_PI.abs(), width, label="raw |rho|", color="#777777")
    ax.bar(x + width/2, show.spearman_residual_log_PI.abs(), width, label="residual |rho|", color="#4daf4a")
    ax.set_xticks(x); ax.set_xticklabels(show.covariate, rotation=25, ha="right")
    ax.set(ylabel="absolute pooled Spearman correlation", title="Covariate dependence before and after adjustment")
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(figures / "covariate_dependence_reduction.png", dpi=180); plt.close(fig)

    plot_comparison = comparison.copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
    axes[0].bar(plot_comparison.adjustment_mode, plot_comparison.burden_curve_AUC_expected_difference, color=["#777777", "#4daf4a", "#984ea3", "#ff7f00"][:len(plot_comparison)])
    axes[0].set(ylabel="burden curve AUC difference", title="Effect size")
    axes[1].bar(plot_comparison.adjustment_mode, -np.log10(plot_comparison.burden_curve_exact_p_expected), color=["#777777", "#4daf4a", "#984ea3", "#ff7f00"][:len(plot_comparison)])
    axes[1].axhline(-math.log10(.05), color="black", ls="--", lw=1)
    axes[1].set(ylabel="-log10 exact p", title="Exact cell-label test")
    for ax in axes:
        ax.tick_params(axis="x", rotation=30)
    fig.tight_layout(); fig.savefig(figures / "adjustment_mode_comparison.png", dpi=180); plt.close(fig)

    if bootstrap is not None:
        ordered = bootstrap.sort_values(["group", "p_exit_by_600"])
        fig, ax = plt.subplots(figsize=(8, 5))
        for group in groups:
            sub = ordered[ordered.group == group]
            ax.scatter(sub.p_exit_by_600, sub.cell_id, color=colors[group], label=group, s=42)
        ax.set(xlabel="Conditional bootstrap probability of adjusted exit by 600 pA", title="WT-calibration uncertainty conditional on residuals")
        ax.legend(frameon=False); fig.tight_layout(); fig.savefig(figures / "conditional_bootstrap_exit_stability.png", dpi=180); plt.close(fig)


def _write_summary(output: Path, comparison: pd.DataFrame, curve_tests: pd.DataFrame, i_exit_tests: pd.DataFrame, concordance: pd.DataFrame, ridge_sensitivity: pd.DataFrame):
    primary_cmp = comparison[comparison.adjustment_mode == "activity_technical"].iloc[0]
    pred = curve_tests[(curve_tests.adjustment_mode == "activity_technical") & (curve_tests.analysis == "residual_predictive_dynamics")].iloc[0]
    exits = concordance[(~concordance.I_exit_v070_primary_censored)][["cell_id", "I_exit_v070_primary_pA"]]
    exit_text = ", ".join("{}={} pA".format(row.cell_id, int(row.I_exit_v070_primary_pA)) for _, row in exits.iterrows()) or "none"
    lambda_text = "; ".join("lambda {}: WT {}, SCA3 {} ({})".format(row.ridge_lambda, int(row.n_WT_exit_by_600), int(row.n_SCA3_exit_by_600), row.SCA3_exiting_cells or "none") for _, row in ridge_sensitivity.iterrows())
    text = """# NeuroThermo v0.7.0 — predictive-dynamics mechanistic validation

## Analysis unit and cohort

Each cell is an independent analysis unit. Animal identifiers are removed before analysis. The frozen cohort contains 13 WT and 7 SCA3 cells on the 100–600 pA current grid.

## Primary adjustment

Log predictive information is predicted separately at each current by label-blind leave-one-cell-out ridge regression. The primary covariates are log1p firing rate, log mean ISI with an explicit missingness indicator, log baseline noise, and stationary-sample count. WT/SCA3 labels are never used by the residualizer. The primary ridge penalty is lambda=1.0 and was frozen before exact group tests.

## Primary results

- Adjusted two-domain burden curve AUC difference: {burden_auc:.6g}; exact p={burden_p:.6g}.
- Adjusted current-wise burden maxT exact p: {burden_maxt:.6g}.
- Residual predictive-dynamics curve AUC difference: {pred_auc:.6g}; exact p={pred_p:.6g}.
- Adjusted exits: WT {n_wt}/13; SCA3 {n_sca}/7.
- Exiting cells: {exit_text}.
- Adjusted I_exit restricted-mean difference WT minus SCA3: {exit_diff:.6g} pA; exact p={exit_p:.6g}.

## Ridge sensitivity

{lambda_text}.

## Interpretation boundary

Persistence after adjustment supports a predictive-dynamics contribution not explained by the included activity and technical covariates. It does not establish causal mechanism, statistical independence of domains, disease time, irreversibility, or a thermodynamic phase transition. I_exit remains a threshold under imposed current stress.
""".format(
        burden_auc=primary_cmp.burden_curve_AUC_expected_difference,
        burden_p=primary_cmp.burden_curve_exact_p_expected,
        burden_maxt=primary_cmp.burden_maxT_exact_p_expected,
        pred_auc=pred.observed_expected_direction_curve_auc_difference,
        pred_p=pred.exact_p_curve_auc_expected_direction,
        n_wt=int(primary_cmp.n_WT_exit_by_600), n_sca=int(primary_cmp.n_SCA3_exit_by_600),
        exit_text=exit_text, exit_diff=primary_cmp.I_exit_restricted_mean_difference_WT_minus_SCA3_pA,
        exit_p=primary_cmp.I_exit_exact_p_expected, lambda_text=lambda_text,
    )
    (output / "RUN_SUMMARY.md").write_text(text, encoding="utf-8")


def run_pipeline(config: Mapping[str, Any], upstream: Path, v061: Path, output: Path) -> None:
    checks = validate_inputs(config, upstream, v061)
    if not bool(checks.passed.all()):
        raise ValueError("Input validation failed:\n" + checks.to_string(index=False))
    output.mkdir(parents=True, exist_ok=True)
    checks.to_csv(output / "input_validation.csv", index=False)

    base_config = copy.deepcopy(dict(config))
    raw_v061, cell_ids, labels, currents, base_matrices = engine._prepare(base_config, upstream)
    adjustment_table, variation = _prepare_adjustment_table(config, upstream, cell_ids, currents)
    residual_table, coefficients, mode_residuals, lambda_residuals = _residualize(config, adjustment_table, cell_ids, currents)
    performance, pooled_dependence, current_dependence = _residual_diagnostics(config, residual_table)
    residual_table.to_csv(output / "crossfit_residualized_predictive_information.csv", index=False)
    coefficients.to_csv(output / "crossfit_ridge_coefficients.csv", index=False)
    variation.to_csv(output / "technical_and_activity_covariate_variation.csv", index=False)
    performance.to_csv(output / "residualization_performance.csv", index=False)
    pooled_dependence.to_csv(output / "covariate_dependence_pooled.csv", index=False)
    current_dependence.to_csv(output / "covariate_dependence_by_current.csv", index=False)

    masks = engine._label_masks(len(cell_ids), int(labels.sum()), int(config["analysis"]["exact_max_labelings"]))
    results = {}
    primary_mode = str(config["residualization"]["primary_mode"])
    mode_order = [primary_mode] + [x for x in config["adjustment_modes"] if x != primary_mode]
    for mode in mode_order:
        matrix = _matrix_from_rows(residual_table, mode_residuals[mode], cell_ids, currents)
        results[mode] = _run_exact_mode(
            config, mode, matrix, base_matrices, cell_ids, labels, currents, masks,
            run_conditional_bootstrap=(mode == primary_mode),
        )

    scores = pd.concat([results[x]["scores"] for x in mode_order], ignore_index=True)
    references = pd.concat([results[x]["references"] for x in mode_order], ignore_index=True)
    cell_summaries = pd.concat([results[x]["cell_summary"] for x in mode_order], ignore_index=True)
    current_tests = pd.concat([results[x]["current_tests"] for x in mode_order], ignore_index=True)
    curve_tests = pd.concat([results[x]["curve_tests"] for x in mode_order], ignore_index=True)
    i_exit_tests = pd.concat([results[x]["i_exit"] for x in mode_order], ignore_index=True)
    bootstrap = results[primary_mode]["conditional_bootstrap"]
    sensitivity = results[primary_mode]["sensitivity"]
    scores.to_csv(output / "mode_cell_current_scores.csv", index=False)
    references.to_csv(output / "mode_target_specific_references.csv", index=False)
    cell_summaries.to_csv(output / "mode_cell_vulnerability_summary.csv", index=False)
    current_tests.to_csv(output / "mode_currentwise_exact_tests.csv", index=False)
    curve_tests.to_csv(output / "mode_group_curve_exact_tests.csv", index=False)
    i_exit_tests.to_csv(output / "mode_I_exit_exact_tests.csv", index=False)
    bootstrap.to_csv(output / "primary_conditional_bootstrap_cell_vulnerability.csv", index=False)
    sensitivity.to_csv(output / "primary_vulnerability_sensitivity_map.csv", index=False)

    curves = _group_curves(scores, residual_table)
    curves.to_csv(output / "mode_group_dynamic_curves.csv", index=False)
    comparison = _comparison(v061, curve_tests, i_exit_tests)
    comparison.to_csv(output / "v061_v070_mode_comparison.csv", index=False)
    concordance = _exit_concordance(v061, results[primary_mode]["cell_summary"])
    concordance.to_csv(output / "v061_v070_exit_concordance.csv", index=False)
    ridge_sensitivity = _ridge_sensitivity(config, residual_table, lambda_residuals, base_matrices, cell_ids, labels, currents, masks)
    ridge_sensitivity.to_csv(output / "ridge_lambda_sensitivity.csv", index=False)

    rule = {
        "pipeline_version": __version__, "analysis_unit": "cell", "animal_level_inference": False,
        "currents_pA": [int(x) for x in currents],
        "residualizer": {
            "label_blind": True, "crossfit_unit": "entire cell", "fit_separately_by_current": True,
            "outcome": "log(predictive_information_nats)", "model": "ridge linear regression",
            "primary_mode": primary_mode, "primary_covariates": list(config["adjustment_modes"][primary_mode]),
            "ridge_lambda": float(config["residualization"]["ridge_lambda"]),
            "ridge_lambda_sensitivity": [float(x) for x in config["residualization"]["ridge_lambda_sensitivity"]],
            "missing_mean_ISI": "training-fold median imputation plus explicit missingness indicator",
            "constant_covariates": "reported but excluded automatically within each current/fold",
        },
        "post_adjustment_scoring": {
            "excitability_timing": "same v0.6.1 firing-rate/mean-ISI domain",
            "predictive_dynamics": "WT-oriented robust z of label-blind cross-fitted log-PI residual",
            "burden": "mean of excitability/timing and residual predictive-dynamics z",
            "target_specific_nested_LOO_WT_calibration": True,
            "strict_exit": {
                "both_domains_above_z": float(config["analysis"]["robust_z_threshold"]),
                "burden_above_z": float(config["analysis"]["robust_z_threshold"]),
                "rank_p_at_or_below": float(config["analysis"]["rank_alpha"]),
                "persistence_adjacent_current_levels": int(config["analysis"]["persistence_levels"]),
                "censor_current_pA": int(config["analysis"]["censor_current_pA"]),
            },
        },
        "exact_inference": {"cell_labelings": int(len(masks)), "recomputes_target_specific_WT_calibration": True},
        "bootstrap_boundary": "conditional on fixed label-blind cross-fitted residuals; WT calibration is recomputed, residualizer is not refit",
        "interpretation": "current-resolved adjusted vulnerability, not disease time or thermodynamic phase transition",
    }
    _json_dump(output / "frozen_adjustment_and_scoring_rule.json", rule)
    _plots(output, config, residual_table, curves, scores, comparison, concordance, pooled_dependence, bootstrap)
    _write_summary(output, comparison, curve_tests, i_exit_tests, concordance, ridge_sensitivity)

    input_hashes = {"v0.3.1/" + item: _sha256(upstream / item) for item in REQUIRED_UPSTREAM}
    input_hashes.update({"v0.6.1/" + item: _sha256(v061 / item) for item in REQUIRED_V061})
    output_files = sorted(x for x in output.rglob("*") if x.is_file() and x.name != "analysis_manifest.json")
    manifest = {
        "pipeline_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(), "python_minimum": "3.9",
        "upstream_dir": str(upstream), "v061_results_dir": str(v061), "output_dir": str(output),
        "input_sha256": input_hashes, "output_sha256": {str(x.relative_to(output)): _sha256(x) for x in output_files},
        "n_cells": int(len(cell_ids)), "n_WT_cells": int((~labels).sum()), "n_SCA3_cells": int(labels.sum()),
        "n_currents": int(len(currents)), "n_exact_labelings": int(len(masks)), "animal_level_inference": False,
    }
    _json_dump(output / "analysis_manifest.json", manifest)
