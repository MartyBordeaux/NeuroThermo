from __future__ import annotations

import itertools
import math
from typing import Dict, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


def label_masks(n_cells: int, n_sca3: int, maximum: int) -> np.ndarray:
    total = math.comb(n_cells, n_sca3)
    if total > int(maximum):
        raise ValueError("Exact labeling count {} exceeds limit {}".format(total, maximum))
    masks = np.zeros((total, n_cells), dtype=bool)
    for row, combination in enumerate(itertools.combinations(range(n_cells), n_sca3)):
        masks[row, list(combination)] = True
    return masks


def cell_auc(matrix: np.ndarray, currents: Sequence[int]) -> np.ndarray:
    values = np.asarray(matrix, dtype=float)
    x = np.asarray(currents, dtype=float)
    result = np.full(values.shape[0], np.nan, dtype=float)
    for index, row in enumerate(values):
        valid = np.isfinite(row)
        if valid.sum() >= 2:
            result[index] = np.trapz(row[valid], x[valid]) / float(x[valid][-1] - x[valid][0])
    return result


def exact_cell_difference(
    values: np.ndarray, masks: np.ndarray, observed_labels: np.ndarray, expected_direction: float
) -> Tuple[float, float, float, np.ndarray]:
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    numeric_masks = masks.astype(float)
    numeric_controls = (~masks).astype(float)
    case_sum = numeric_masks @ np.where(valid, values, 0.0)
    case_n = numeric_masks @ valid.astype(float)
    control_sum = numeric_controls @ np.where(valid, values, 0.0)
    control_n = numeric_controls @ valid.astype(float)
    differences = case_sum / case_n - control_sum / control_n
    matches = np.flatnonzero(np.all(masks == observed_labels[None, :], axis=1))
    if len(matches) != 1:
        raise RuntimeError("Observed labeling was not uniquely found.")
    observed = float(differences[int(matches[0])])
    oriented = float(expected_direction) * differences
    observed_oriented = float(expected_direction) * observed
    p_expected = float(np.mean(oriented >= observed_oriented - 1e-12))
    p_two_sided = float(np.mean(np.abs(differences) >= abs(observed) - 1e-12))
    return observed, p_expected, p_two_sided, differences


def exact_fixed_metric_tests(
    table: pd.DataFrame,
    metric_names: Sequence[str],
    cell_ids: Sequence[str],
    labels: np.ndarray,
    currents: Sequence[int],
    masks: np.ndarray,
    expected_direction: float = -1.0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    auc_rows = []
    current_rows = []
    for metric in metric_names:
        matrix = table.pivot(index="cell_id", columns="current_pA", values=metric).reindex(
            index=cell_ids, columns=currents
        ).to_numpy(float)
        auc_values = cell_auc(matrix, currents)
        observed, p_expected, p_two, _ = exact_cell_difference(
            auc_values, masks, labels, expected_direction
        )
        auc_rows.append({
            "metric": metric,
            "expected_direction_multiplier": float(expected_direction),
            "observed_raw_AUC_difference_SCA3_minus_WT": observed,
            "observed_expected_direction_AUC_difference": float(expected_direction) * observed,
            "exact_p_AUC_expected_direction": p_expected,
            "exact_p_AUC_two_sided": p_two,
            "n_exact_labelings": int(len(masks)),
            "current_min_pA": int(min(currents)),
            "current_max_pA": int(max(currents)),
        })

        differences = []
        observed_by_current = []
        p_expected_by_current = []
        p_two_by_current = []
        for k, current in enumerate(currents):
            obs, p_exp, p_two_current, null = exact_cell_difference(
                matrix[:, k], masks, labels, expected_direction
            )
            differences.append(null)
            observed_by_current.append(obs)
            p_expected_by_current.append(p_exp)
            p_two_by_current.append(p_two_current)
        difference_matrix = np.column_stack(differences)
        oriented_matrix = float(expected_direction) * difference_matrix
        max_expected = np.nanmax(oriented_matrix, axis=1)
        max_absolute = np.nanmax(np.abs(difference_matrix), axis=1)
        for k, current in enumerate(currents):
            observed = float(observed_by_current[k])
            observed_oriented = float(expected_direction) * observed
            current_rows.append({
                "metric": metric,
                "current_pA": int(current),
                "expected_direction_multiplier": float(expected_direction),
                "observed_raw_difference_SCA3_minus_WT": observed,
                "observed_expected_direction_difference": observed_oriented,
                "exact_p_expected_direction_unadjusted": float(p_expected_by_current[k]),
                "exact_p_two_sided_unadjusted": float(p_two_by_current[k]),
                "exact_p_expected_direction_maxT_adjusted": float(np.mean(max_expected >= observed_oriented - 1e-12)),
                "exact_p_two_sided_maxT_adjusted": float(np.mean(max_absolute >= abs(observed) - 1e-12)),
                "n_exact_labelings": int(len(masks)),
            })
    return pd.DataFrame(auc_rows), pd.DataFrame(current_rows)


def exact_lag_family_tests(
    table: pd.DataFrame,
    families: Mapping[str, Mapping[int, str]],
    cell_ids: Sequence[str],
    labels: np.ndarray,
    currents: Sequence[int],
    masks: np.ndarray,
    expected_direction: float = -1.0,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Exact AUC and current tests with maxT within each surrogate family."""
    auc_rows = []
    current_rows = []
    cell_auc_rows = []
    group_by_cell = table[["cell_id", "group"]].drop_duplicates().set_index("cell_id")["group"].to_dict()
    for family, lag_columns in families.items():
        auc_nulls = []
        auc_records = []
        current_nulls = []
        current_records = []
        for lag_ms, column in sorted(lag_columns.items()):
            matrix = table.pivot(index="cell_id", columns="current_pA", values=column).reindex(
                index=cell_ids, columns=currents
            ).to_numpy(float)
            auc_values = cell_auc(matrix, currents)
            observed, p_expected, p_two, null = exact_cell_difference(
                auc_values, masks, labels, expected_direction
            )
            auc_nulls.append(null)
            auc_records.append((lag_ms, column, observed, p_expected, p_two))
            for cell, value in zip(cell_ids, auc_values):
                cell_auc_rows.append({
                    "surrogate_family": family,
                    "lag_ms": int(lag_ms),
                    "metric": column,
                    "group": str(group_by_cell[cell]),
                    "cell_id": str(cell),
                    "cell_AUC_nats": float(value),
                })
            for k, current in enumerate(currents):
                current_observed, current_p, current_two, current_null = exact_cell_difference(
                    matrix[:, k], masks, labels, expected_direction
                )
                current_nulls.append(current_null)
                current_records.append((lag_ms, column, int(current), current_observed, current_p, current_two))

        auc_null_matrix = np.column_stack(auc_nulls)
        auc_max_expected = np.nanmax(float(expected_direction) * auc_null_matrix, axis=1)
        auc_max_absolute = np.nanmax(np.abs(auc_null_matrix), axis=1)
        for lag_ms, column, observed, p_expected, p_two in auc_records:
            oriented = float(expected_direction) * float(observed)
            auc_rows.append({
                "surrogate_family": family,
                "lag_ms": int(lag_ms),
                "metric": column,
                "primary_prespecified": bool(family == "shuffle" and int(lag_ms) == 4),
                "expected_direction_multiplier": float(expected_direction),
                "observed_raw_AUC_difference_SCA3_minus_WT": float(observed),
                "observed_expected_direction_AUC_difference_WT_minus_SCA3": oriented,
                "exact_p_AUC_expected_direction_unadjusted": float(p_expected),
                "exact_p_AUC_two_sided_unadjusted": float(p_two),
                "exact_p_AUC_expected_direction_maxT_across_lags": float(np.mean(auc_max_expected >= oriented - 1e-12)),
                "exact_p_AUC_two_sided_maxT_across_lags": float(np.mean(auc_max_absolute >= abs(observed) - 1e-12)),
                "n_exact_labelings": int(len(masks)),
            })

        current_null_matrix = np.column_stack(current_nulls)
        current_max_expected = np.nanmax(float(expected_direction) * current_null_matrix, axis=1)
        current_max_absolute = np.nanmax(np.abs(current_null_matrix), axis=1)
        for lag_ms, column, current, observed, p_expected, p_two in current_records:
            oriented = float(expected_direction) * float(observed)
            current_rows.append({
                "surrogate_family": family,
                "lag_ms": int(lag_ms),
                "current_pA": int(current),
                "metric": column,
                "expected_direction_multiplier": float(expected_direction),
                "observed_raw_difference_SCA3_minus_WT": float(observed),
                "observed_expected_direction_difference_WT_minus_SCA3": oriented,
                "exact_p_expected_direction_unadjusted": float(p_expected),
                "exact_p_two_sided_unadjusted": float(p_two),
                "exact_p_expected_direction_maxT_across_lags_and_currents": float(np.mean(current_max_expected >= oriented - 1e-12)),
                "exact_p_two_sided_maxT_across_lags_and_currents": float(np.mean(current_max_absolute >= abs(observed) - 1e-12)),
                "n_exact_labelings": int(len(masks)),
            })
    return pd.DataFrame(auc_rows), pd.DataFrame(current_rows), pd.DataFrame(cell_auc_rows)


def transform_covariate(values: pd.Series, transform: str) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    result = np.full(x.shape, np.nan, dtype=float)
    finite = np.isfinite(x)
    if transform == "identity":
        result[finite] = x[finite]
    elif transform == "log":
        valid = finite & (x > 0)
        result[valid] = np.log(x[valid])
    elif transform == "log1p":
        valid = finite & (x >= 0)
        result[valid] = np.log1p(x[valid])
    else:
        raise ValueError("Unknown transform: {}".format(transform))
    return result


def prepare_covariates(config: Mapping, table: pd.DataFrame) -> pd.DataFrame:
    output = table.copy()
    output["mean_isi_missing"] = output.mean_isi_ms.isna().astype(float)
    for name, specification in config["covariates"].items():
        source = str(specification.get("source", name))
        if source == "mean_isi_missing":
            output[name] = output["mean_isi_missing"].to_numpy(float)
        else:
            output[name] = transform_covariate(output[source], str(specification["transform"]))
    return output


def ridge_crossfit_by_cell_current(
    table: pd.DataFrame,
    cell_ids: Sequence[str],
    currents: Sequence[int],
    outcome: str,
    covariates: Sequence[str],
    ridge_lambda: float,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    prediction = np.full(len(table), np.nan, dtype=float)
    coefficients = []
    for current in currents:
        current_mask = table.current_pA.to_numpy(int) == int(current)
        for cell in cell_ids:
            target_mask = current_mask & (table.cell_id.to_numpy(str) == str(cell))
            target_indices = np.flatnonzero(target_mask)
            if len(target_indices) != 1:
                raise ValueError("Expected one row for {} at {} pA".format(cell, current))
            train = table.loc[current_mask & ~target_mask]
            test = table.loc[target_mask]
            y = train[outcome].to_numpy(float)
            if not np.isfinite(y).all():
                raise ValueError("Non-finite residualization outcome in training data.")
            train_columns = []
            test_columns = []
            metadata = []
            active_names = []
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
            beta = np.linalg.solve(
                x_train_matrix.T @ x_train_matrix + penalty, x_train_matrix.T @ y
            )
            prediction[target_indices] = x_test_matrix @ beta
            coefficients.append({
                "target_cell_id": str(cell), "current_pA": int(current),
                "term": "intercept", "standardized_coefficient": float(beta[0]),
                "active": True, "training_n_cells": int(len(train)),
                "ridge_lambda": float(ridge_lambda),
            })
            coefficient_by_name: Dict[str, float] = dict(zip(active_names, beta[1:]))
            for name, median, center, scale, active in metadata:
                coefficients.append({
                    "target_cell_id": str(cell), "current_pA": int(current),
                    "term": name,
                    "standardized_coefficient": float(coefficient_by_name[name]) if active else 0.0,
                    "active": active, "training_n_cells": int(len(train)),
                    "ridge_lambda": float(ridge_lambda),
                    "training_imputation_median": median,
                    "training_center": center, "training_scale": scale,
                })
    observed = table[outcome].to_numpy(float)
    residual = observed - prediction
    return prediction, residual, pd.DataFrame(coefficients)
