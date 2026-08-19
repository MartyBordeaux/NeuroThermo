from __future__ import annotations

from functools import lru_cache
from itertools import combinations
from math import comb
import warnings

import numpy as np
import pandas as pd


ELECTRICAL_WORK_FEATURES = [
    "external_work_signed_fJ",
    "external_work_positive_fJ",
    "mean_power_signed_fW",
]

SPIKE_CONDITIONAL_THERMO_FEATURES = [
    "work_per_spike_fJ",
    "permutation_entropy_norm",
    "spectral_entropy_norm",
    "predictive_information_nats",
    "path_kl_rate_raw_nats_s",
    "path_kl_rate_bias_corrected_nats_s",
    "path_kl_rate_excess_nats_s",
]

THERMO_FEATURES = [*ELECTRICAL_WORK_FEATURES, *SPIKE_CONDITIONAL_THERMO_FEATURES]

CONVENTIONAL_FEATURES = [
    "firing_rate_hz", "sustained_rate_hz", "first_spike_latency_ms",
    "mean_isi_ms", "cv_isi", "adaptation_ratio", "baseline_voltage_mV",
]

FIXED_DOMAIN_INTEGRATED_FEATURES = [
    *ELECTRICAL_WORK_FEATURES,
    "firing_rate_hz", "sustained_rate_hz", "baseline_voltage_mV",
]

FEATURE_FAMILIES = {
    **{name: "electrical_work" for name in ELECTRICAL_WORK_FEATURES},
    "work_per_spike_fJ": "electrical_efficiency",
    "permutation_entropy_norm": "information_dynamics",
    "spectral_entropy_norm": "information_dynamics",
    "predictive_information_nats": "information_dynamics",
    "path_kl_rate_raw_nats_s": "observable_irreversibility",
    "path_kl_rate_bias_corrected_nats_s": "observable_irreversibility",
    "path_kl_rate_excess_nats_s": "observable_irreversibility",
    **{name: "electrophysiology" for name in CONVENTIONAL_FEATURES},
}


def _requires_thermo_eligibility(feature: str) -> bool:
    return feature in SPIKE_CONDITIONAL_THERMO_FEATURES


def _add_fdr_columns(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    out["fdr_q_global"] = benjamini_hochberg(out["permutation_p"])
    out["fdr_q_family"] = np.nan
    for _, idx in out.groupby("family").groups.items():
        out.loc[idx, "fdr_q_family"] = benjamini_hochberg(out.loc[idx, "permutation_p"])
    out["fdr_q"] = out["fdr_q_global"]
    return out


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    return float(np.mean(a[:, None] > b[None, :]) - np.mean(a[:, None] < b[None, :]))


def _bootstrap_difference(a, b, n_iter, rng):
    diffs = np.empty(n_iter, float)
    for i in range(n_iter):
        aa = rng.choice(a, len(a), replace=True)
        bb = rng.choice(b, len(b), replace=True)
        diffs[i] = np.median(bb) - np.median(aa)
    return np.quantile(diffs, [0.025, 0.975])


@lru_cache(maxsize=None)
def _exact_group_masks(n_total: int, n_b: int) -> np.ndarray:
    total = comb(n_total, n_b)
    masks = np.zeros((total, n_total), dtype=bool)
    for row_index, b_indices in enumerate(combinations(range(n_total), n_b)):
        masks[row_index, list(b_indices)] = True
    masks.setflags(write=False)
    return masks


def _permutation_p(a, b, n_iter, rng, exact_max=200000):
    observed = abs(float(np.median(b) - np.median(a)))
    combined = np.r_[a, b]
    n_a = len(a)
    n_b = len(b)
    total = comb(len(combined), n_b)
    if total <= int(exact_max):
        masks = _exact_group_masks(len(combined), n_b)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            b_medians = np.nanmedian(
                np.where(masks, combined[None, :], np.nan), axis=1
            )
            a_medians = np.nanmedian(
                np.where(~masks, combined[None, :], np.nan), axis=1
            )
        differences = np.abs(b_medians - a_medians)
        return (
            float(np.mean(differences >= observed - 1e-12)),
            "exact", int(total),
        )
    count = 0
    for _ in range(n_iter):
        perm = rng.permutation(combined)
        diff = abs(float(np.median(perm[n_a:]) - np.median(perm[:n_a])))
        count += diff >= observed
    return float((count + 1) / (n_iter + 1)), "monte_carlo", int(n_iter)


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    p = p_values.to_numpy(float)
    q = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if len(valid) == 0:
        return q
    order = valid[np.argsort(p[valid])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def compare_groups_by_current(
    features: pd.DataFrame, cfg: dict, feature_names=None
) -> pd.DataFrame:
    stats_cfg = cfg["statistics"]
    rng = np.random.default_rng(int(cfg["seed"]))
    rows = []
    frame = features[features["qc_pass"]].copy()
    frame["current_key"] = frame["current_pA"].round(int(stats_cfg["current_round_decimals"]))
    selected = feature_names or [*THERMO_FEATURES, *CONVENTIONAL_FEATURES]
    for feature in selected:
        source = frame[frame["thermo_eligible"]].copy() if _requires_thermo_eligibility(feature) else frame
        if feature not in source:
            continue
        for current, level in source.groupby("current_key"):
            a = level.loc[level.group == "WT", feature].dropna().to_numpy(float)
            b = level.loc[level.group == "SCA3", feature].dropna().to_numpy(float)
            minimum = int(stats_cfg["minimum_cells_per_group"])
            if len(a) < minimum or len(b) < minimum:
                continue
            lo, hi = _bootstrap_difference(a, b, int(stats_cfg["bootstrap_iterations"]), rng)
            p_value, permutation_mode, valid_labelings = _permutation_p(
                a, b, int(stats_cfg["permutation_iterations"]), rng,
                int(stats_cfg.get("exact_max_labelings", 200000)),
            )
            rows.append({
                "feature": feature, "family": FEATURE_FAMILIES.get(feature, "other"),
                "current_pA": current,
                "n_WT": len(a), "n_SCA3": len(b),
                "median_WT": float(np.median(a)), "median_SCA3": float(np.median(b)),
                "median_difference_SCA3_minus_WT": float(np.median(b) - np.median(a)),
                "difference_ci95_low": float(lo), "difference_ci95_high": float(hi),
                "cliffs_delta_SCA3_vs_WT": cliffs_delta(b, a),
                "permutation_mode": permutation_mode,
                "valid_labelings": valid_labelings,
                "permutation_p": p_value,
            })
    result = pd.DataFrame(rows)
    return _add_fdr_columns(result)


def rheobase_brackets(features: pd.DataFrame) -> pd.DataFrame:
    rows = []
    valid = features[features["qc_pass"]].copy()
    for (group, cell_id), cell in valid.groupby(["group", "cell_id"]):
        if "threshold_probe_spiking" in cell:
            threshold_spiking = cell["threshold_probe_spiking"].fillna(False).astype(bool)
        else:
            threshold_spiking = cell["n_spikes"] > 0
        spiking = np.sort(cell.loc[threshold_spiking, "current_pA"].dropna().unique())
        silent = np.sort(cell.loc[~threshold_spiking, "current_pA"].dropna().unique())
        upper = float(spiking[0]) if len(spiking) else np.nan
        lower_candidates = silent[silent < upper] if np.isfinite(upper) else silent
        lower = float(lower_candidates[-1]) if len(lower_candidates) else np.nan
        rows.append({
            "group": group, "cell_id": cell_id,
            "rheobase_lower_nonspiking_pA": lower,
            "rheobase_upper_spiking_pA": upper,
            "rheobase_midpoint_pA": (lower + upper) / 2 if np.isfinite(lower) and np.isfinite(upper) else upper,
            "rheobase_bracket_width_pA": upper - lower if np.isfinite(lower) and np.isfinite(upper) else np.nan,
        })
    return pd.DataFrame(rows)


def integrated_cell_features(features: pd.DataFrame, cfg: dict, return_coverage: bool = False):
    rows = []
    coverage_rows = []
    valid = features[features["qc_pass"]].copy()
    grid = np.asarray(cfg["input"].get("analysis_currents_pA") or [], float)
    tolerance = float(cfg["input"].get("current_tolerance_pA", 1e-6))
    if len(grid) < 2:
        raise ValueError("Fixed-domain integration requires at least two analysis_currents_pA")
    domain_label = f"{grid.min():g}_{grid.max():g}pA"
    for (group, cell_id), cell in valid.groupby(["group", "cell_id"]):
        row = {"group": group, "cell_id": cell_id}
        for feature in FIXED_DOMAIN_INTEGRATED_FEATURES:
            curve = cell[["current_pA", feature]].dropna().groupby("current_pA", as_index=False).median()
            values = []
            missing = []
            for current in grid:
                matched = curve[np.abs(curve["current_pA"].to_numpy(float) - current) <= tolerance]
                if len(matched) == 1:
                    values.append(float(matched[feature].iloc[0]))
                else:
                    values.append(np.nan)
                    missing.append(float(current))
            complete = not missing
            coverage_rows.append({
                "group": group, "cell_id": cell_id, "feature": feature,
                "expected_levels": int(len(grid)), "observed_levels": int(len(grid) - len(missing)),
                "complete_common_domain": complete,
                "missing_currents_pA": ";".join(f"{x:g}" for x in missing),
                "domain_min_pA": float(grid.min()), "domain_max_pA": float(grid.max()),
            })
            if complete:
                row[f"{feature}__auc_mean_{domain_label}"] = float(
                    np.trapz(np.asarray(values, float), grid) / (grid.max() - grid.min())
                )
        rows.append(row)
    integrated = pd.DataFrame(rows)
    coverage = pd.DataFrame(coverage_rows)
    return (integrated, coverage) if return_coverage else integrated


def compare_integrated_cells(integrated: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rng = np.random.default_rng(int(cfg["seed"]) + 1)
    scfg = cfg["statistics"]
    rows = []
    for feature in [c for c in integrated.columns if "__auc_mean_" in c]:
        a = integrated.loc[integrated.group == "WT", feature].dropna().to_numpy(float)
        b = integrated.loc[integrated.group == "SCA3", feature].dropna().to_numpy(float)
        if min(len(a), len(b)) < int(scfg["minimum_cells_per_group"]):
            continue
        lo, hi = _bootstrap_difference(a, b, int(scfg["bootstrap_iterations"]), rng)
        p_value, permutation_mode, valid_labelings = _permutation_p(
            a, b, int(scfg["permutation_iterations"]), rng,
            int(scfg.get("exact_max_labelings", 200000)),
        )
        rows.append({
            "feature": feature, "family": FEATURE_FAMILIES.get(feature.split("__", 1)[0], "other"),
            "n_WT": len(a), "n_SCA3": len(b),
            "median_WT": float(np.median(a)), "median_SCA3": float(np.median(b)),
            "median_difference_SCA3_minus_WT": float(np.median(b) - np.median(a)),
            "difference_ci95_low": float(lo), "difference_ci95_high": float(hi),
            "cliffs_delta_SCA3_vs_WT": cliffs_delta(b, a),
            "permutation_mode": permutation_mode,
            "valid_labelings": valid_labelings,
            "permutation_p": p_value,
        })
    out = pd.DataFrame(rows)
    return _add_fdr_columns(out)


def response_summary(features: pd.DataFrame, axis: str, feature_names: list[str], cfg: dict) -> pd.DataFrame:
    rows = []
    frame = features[features["qc_pass"]].copy()
    if axis == "current_pA":
        frame["axis_value"] = frame[axis].round(int(cfg["statistics"]["current_round_decimals"]))
        grouped = frame.groupby(["group", "axis_value"])
        for (group, x), level in grouped:
            for feature in feature_names:
                eligible = level[level["thermo_eligible"]] if _requires_thermo_eligibility(feature) else level
                values = eligible[feature].dropna().to_numpy(float)
                if len(values):
                    rows.append({"group": group, "axis": axis, "axis_value": x, "feature": feature, "n_cells": len(values), "median": np.median(values), "q25": np.quantile(values, .25), "q75": np.quantile(values, .75)})
        return _annotate_shared_support(pd.DataFrame(rows), int(cfg["statistics"]["minimum_cells_per_group"]))
    step = float(cfg["statistics"]["J_grid_step_pA_per_pF"])
    finite = frame[axis].replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return pd.DataFrame()
    grid = np.arange(0.0, finite.max() + step / 2, step)
    for group, group_frame in frame.groupby("group"):
        for feature in feature_names:
            cell_curves = []
            for _, cell in group_frame.groupby("cell_id"):
                source = cell[cell["thermo_eligible"]] if _requires_thermo_eligibility(feature) else cell
                curve = source[[axis, feature]].replace([np.inf, -np.inf], np.nan).dropna().groupby(axis, as_index=False).median().sort_values(axis)
                if len(curve) < 3:
                    continue
                values = np.full(len(grid), np.nan)
                inside = (grid >= curve[axis].min()) & (grid <= curve[axis].max())
                values[inside] = np.interp(grid[inside], curve[axis], curve[feature])
                cell_curves.append(values)
            if not cell_curves:
                continue
            matrix = np.asarray(cell_curves)
            for j, x in enumerate(grid):
                values = matrix[:, j]
                values = values[np.isfinite(values)]
                if len(values):
                    rows.append({"group": group, "axis": axis, "axis_value": x, "feature": feature, "n_cells": len(values), "median": np.median(values), "q25": np.quantile(values, .25), "q75": np.quantile(values, .75)})
    return _annotate_shared_support(pd.DataFrame(rows), int(cfg["statistics"]["minimum_cells_per_group"]))


def _annotate_shared_support(summary: pd.DataFrame, minimum: int) -> pd.DataFrame:
    if summary.empty:
        return summary
    counts = summary.pivot_table(
        index=["axis", "axis_value", "feature"], columns="group", values="n_cells", aggfunc="max"
    ).fillna(0)
    counts["shared_support"] = (counts.get("WT", 0) >= minimum) & (counts.get("SCA3", 0) >= minimum)
    support = counts[["shared_support"]].reset_index()
    return summary.merge(support, on=["axis", "axis_value", "feature"], how="left")


def current_density_support(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if summary.empty:
        return pd.DataFrame(rows)
    for (axis, feature), part in summary.groupby(["axis", "feature"]):
        shared = np.sort(part.loc[part["shared_support"], "axis_value"].unique())
        rows.append({
            "axis": axis, "feature": feature,
            "shared_support_points": int(len(shared)),
            "shared_support_min": float(shared.min()) if len(shared) else np.nan,
            "shared_support_max": float(shared.max()) if len(shared) else np.nan,
            "shared_support_available": bool(len(shared)),
        })
    return pd.DataFrame(rows)


def compare_response_curves(
    features: pd.DataFrame, cfg: dict, feature_names=None
) -> pd.DataFrame:
    frame = features[features["qc_pass"]].copy()
    grid = np.asarray(cfg["input"].get("analysis_currents_pA") or [], float)
    minimum = int(cfg["statistics"]["minimum_cells_per_group"])
    iterations = int(cfg["statistics"]["permutation_iterations"])
    rng = np.random.default_rng(int(cfg["seed"]) + 2)
    rows = []
    selected = feature_names or [*THERMO_FEATURES, *CONVENTIONAL_FEATURES]
    for feature in selected:
        source = frame[frame["thermo_eligible"]].copy() if _requires_thermo_eligibility(feature) else frame.copy()
        if feature not in source:
            continue
        source["cell_key"] = source["group"].astype(str) + "|" + source["cell_id"].astype(str)
        cell_group = source.groupby("cell_key")["group"].first()
        table = source.pivot_table(index="cell_key", columns="current_pA", values=feature, aggfunc="median")
        levels = []
        for current in grid:
            if current not in table.columns:
                continue
            values = table[current]
            if values[cell_group.eq("WT")].notna().sum() >= minimum and values[cell_group.eq("SCA3")].notna().sum() >= minimum:
                levels.append(float(current))
        if not levels:
            continue
        matrix = table[levels].to_numpy(float)
        labels = cell_group.loc[table.index].to_numpy(str)

        def curve_statistics(group_labels):
            wt_matrix = matrix[group_labels == "WT"]
            sca_matrix = matrix[group_labels == "SCA3"]
            if np.any(np.sum(np.isfinite(wt_matrix), axis=0) == 0) or np.any(np.sum(np.isfinite(sca_matrix), axis=0) == 0):
                return None
            wt = np.nanmedian(wt_matrix, axis=0)
            sca = np.nanmedian(sca_matrix, axis=0)
            difference = sca - wt
            if len(levels) == 1:
                return float(difference[0]), float(abs(difference[0])), float(abs(difference[0]))
            span = levels[-1] - levels[0]
            signed = float(np.trapz(difference, levels) / span)
            absolute = float(np.trapz(np.abs(difference), levels) / span)
            return signed, absolute, float(np.max(np.abs(difference)))

        observed_statistics = curve_statistics(labels)
        if observed_statistics is None:
            continue
        signed, absolute, maximum = observed_statistics
        n_sca = int(np.sum(labels == "SCA3"))
        total_labelings = comb(len(labels), n_sca)
        exact_max = int(cfg["statistics"].get("exact_max_labelings", 200000))
        exceed, valid_labelings = 0, 0
        if total_labelings <= exact_max:
            permutation_mode = "exact"
            masks = _exact_group_masks(len(labels), n_sca)
            differences = np.full((len(masks), len(levels)), np.nan, dtype=float)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                for column in range(len(levels)):
                    values = matrix[:, column]
                    sca_median = np.nanmedian(
                        np.where(masks, values[None, :], np.nan), axis=1
                    )
                    wt_median = np.nanmedian(
                        np.where(~masks, values[None, :], np.nan), axis=1
                    )
                    differences[:, column] = sca_median - wt_median
            valid_rows = np.isfinite(differences).all(axis=1)
            valid_differences = differences[valid_rows]
            if len(levels) == 1:
                permuted_absolute = np.abs(valid_differences[:, 0])
            else:
                permuted_absolute = np.trapz(
                    np.abs(valid_differences), levels, axis=1
                ) / float(levels[-1] - levels[0])
            valid_labelings = int(len(permuted_absolute))
            exceed = int(np.sum(permuted_absolute >= absolute - 1e-12))
            p_value = float(exceed / valid_labelings) if valid_labelings else np.nan
        else:
            permutation_mode = "monte_carlo"
            attempts = 0
            while valid_labelings < iterations and attempts < iterations * 20:
                attempts += 1
                permuted = curve_statistics(rng.permutation(labels))
                if permuted is None:
                    continue
                _, permuted_absolute, _ = permuted
                exceed += permuted_absolute >= absolute
                valid_labelings += 1
            p_value = float((exceed + 1) / (valid_labelings + 1)) if valid_labelings else np.nan
        if valid_labelings == 0:
            continue
        rows.append({
            "feature": feature, "family": FEATURE_FAMILIES.get(feature, "other"),
            "n_WT_cells": int(np.sum(labels == "WT")), "n_SCA3_cells": int(np.sum(labels == "SCA3")),
            "n_shared_current_levels": int(len(levels)),
            "shared_current_min_pA": float(min(levels)), "shared_current_max_pA": float(max(levels)),
            "signed_curve_difference_SCA3_minus_WT": signed,
            "absolute_curve_difference": absolute, "maximum_absolute_difference": maximum,
            "permutation_mode": permutation_mode,
            "valid_labelings": int(valid_labelings),
            "permutation_p": p_value,
        })
    return _add_fdr_columns(pd.DataFrame(rows))


def _curve_triplet(difference: np.ndarray, levels: np.ndarray):
    difference = np.asarray(difference, float)
    levels = np.asarray(levels, float)
    if len(levels) == 1:
        value = float(difference[0])
        return value, abs(value), abs(value)
    span = float(levels[-1] - levels[0])
    signed = float(np.trapz(difference, levels) / span)
    absolute = float(np.trapz(np.abs(difference), levels) / span)
    maximum = float(np.max(np.abs(difference)))
    return signed, absolute, maximum


def _survival_probabilities(values: np.ndarray, exact: bool) -> np.ndarray:
    values = np.asarray(values, float)
    ordered = np.sort(values)
    exceed = len(values) - np.searchsorted(ordered, values, side="left")
    if exact:
        return exceed / len(values)
    return (exceed + 1) / (len(values) + 1)


def compare_two_part_all_cells(features: pd.DataFrame, cfg: dict):
    """All-cell availability + conditional-value curve tests.

    Every accepted cell enters the binary availability component. Conditional
    values are used only where physically defined; no non-spiking value is
    imputed and no cell is removed globally for an incomplete curve.
    """
    scfg = cfg["statistics"]
    configured_levels = scfg.get("two_part_currents_pA") or []
    primary = scfg.get("two_part_primary_features") or []
    secondary = scfg.get("two_part_secondary_features") or []
    selected = [(feature, "primary") for feature in primary]
    selected += [(feature, "secondary") for feature in secondary]
    comparison_columns = [
        "inference_role", "feature", "family", "n_WT_total", "n_SCA3_total",
        "n_WT_any_defined", "n_SCA3_any_defined", "current_levels_pA",
        "availability_signed_curve_difference_SCA3_minus_WT",
        "availability_absolute_curve_difference",
        "availability_maximum_absolute_difference", "availability_permutation_p",
        "conditional_signed_curve_difference_SCA3_minus_WT",
        "conditional_absolute_curve_difference",
        "conditional_maximum_absolute_difference", "conditional_permutation_p",
        "combined_fisher_statistic", "permutation_mode", "valid_labelings",
        "permutation_p", "status",
    ]
    coverage_columns = [
        "inference_role", "feature", "group", "total_cells",
        "cells_any_defined", "cells_all_levels_defined",
        "median_defined_levels", "minimum_defined_levels",
        "maximum_defined_levels", "required_current_levels_pA",
    ]
    cell_columns = [
        "inference_role", "feature", "group", "cell_id", "current_pA",
        "available", "value",
    ]
    if not configured_levels or not selected:
        return (
            pd.DataFrame(columns=comparison_columns),
            pd.DataFrame(columns=coverage_columns),
            pd.DataFrame(columns=cell_columns),
        )

    levels = np.asarray(configured_levels, float)
    level_label = ";".join(f"{value:g}" for value in levels)
    # The denominator is the accepted biological cohort, not the subset with
    # at least one usable sweep. A cell with no defined conditional value must
    # remain in the availability component as zero availability.
    cells = features[["group", "cell_id"]].drop_duplicates().sort_values(
        ["group", "cell_id"]
    )
    cell_index = pd.MultiIndex.from_frame(cells[["group", "cell_id"]])
    labels = cells["group"].astype(str).to_numpy()
    n_wt = int(np.sum(labels == "WT"))
    n_sca = int(np.sum(labels == "SCA3"))
    exact_max = int(scfg.get("exact_max_labelings", 200000))
    iterations = int(scfg["permutation_iterations"])
    rng = np.random.default_rng(int(cfg["seed"]) + 5)
    comparison_rows, coverage_rows, cell_rows = [], [], []
    total_labelings = comb(len(labels), n_sca)
    if total_labelings <= exact_max:
        permutation_mode = "exact"
        label_masks = _exact_group_masks(len(labels), n_sca)
    else:
        permutation_mode = "monte_carlo"
        label_masks = np.zeros((iterations, len(labels)), dtype=bool)
        for row_index in range(iterations):
            label_masks[
                row_index,
                rng.choice(len(labels), size=n_sca, replace=False),
            ] = True

    for feature, role in selected:
        if feature not in features:
            continue
        source = features[
            features["qc_pass"] & features["current_pA"].isin(levels)
        ].copy()
        if _requires_thermo_eligibility(feature):
            source.loc[~source["thermo_eligible"], feature] = np.nan
        table = source.pivot_table(
            index=["group", "cell_id"], columns="current_pA",
            values=feature, aggfunc="median", dropna=False,
        ).reindex(index=cell_index, columns=levels)
        matrix = table.to_numpy(float)
        availability = np.isfinite(matrix)

        for group in ["WT", "SCA3"]:
            group_mask = labels == group
            defined_levels = availability[group_mask].sum(axis=1)
            coverage_rows.append({
                "inference_role": role, "feature": feature, "group": group,
                "total_cells": int(group_mask.sum()),
                "cells_any_defined": int(np.sum(defined_levels > 0)),
                "cells_all_levels_defined": int(np.sum(defined_levels == len(levels))),
                "median_defined_levels": float(np.median(defined_levels)),
                "minimum_defined_levels": int(np.min(defined_levels)),
                "maximum_defined_levels": int(np.max(defined_levels)),
                "required_current_levels_pA": level_label,
            })
        for (group, cell_id), values in table.iterrows():
            for current, value in zip(levels, values.to_numpy(float)):
                cell_rows.append({
                    "inference_role": role, "feature": feature,
                    "group": group, "cell_id": cell_id,
                    "current_pA": float(current),
                    "available": bool(np.isfinite(value)),
                    "value": float(value) if np.isfinite(value) else np.nan,
                })

        def component_statistics(group_labels):
            wt_mask = group_labels == "WT"
            sca_mask = group_labels == "SCA3"
            availability_difference = (
                availability[sca_mask].mean(axis=0)
                - availability[wt_mask].mean(axis=0)
            )
            availability_stats = _curve_triplet(availability_difference, levels)
            conditional_difference = []
            for column in range(len(levels)):
                wt_values = matrix[wt_mask, column]
                sca_values = matrix[sca_mask, column]
                wt_values = wt_values[np.isfinite(wt_values)]
                sca_values = sca_values[np.isfinite(sca_values)]
                if len(wt_values) == 0 or len(sca_values) == 0:
                    return availability_stats, None
                conditional_difference.append(
                    float(np.median(sca_values) - np.median(wt_values))
                )
            conditional_stats = _curve_triplet(
                np.asarray(conditional_difference, float), levels
            )
            return availability_stats, conditional_stats

        observed_availability, observed_conditional = component_statistics(labels)
        if observed_conditional is None:
            comparison_rows.append({
                "inference_role": role, "feature": feature,
                "family": FEATURE_FAMILIES.get(feature, "other"),
                "n_WT_total": n_wt, "n_SCA3_total": n_sca,
                "n_WT_any_defined": int(np.sum(availability[labels == "WT"].any(axis=1))),
                "n_SCA3_any_defined": int(np.sum(availability[labels == "SCA3"].any(axis=1))),
                "current_levels_pA": level_label,
                "status": "conditional_component_undefined",
            })
            continue

        sca_available = label_masks.astype(np.int16) @ availability.astype(np.int16)
        total_available = availability.sum(axis=0)
        wt_available = total_available[None, :] - sca_available
        availability_difference = (
            sca_available / n_sca - wt_available / n_wt
        )
        if len(levels) == 1:
            availability_absolute_all = np.abs(availability_difference[:, 0])
        else:
            availability_absolute_all = np.trapz(
                np.abs(availability_difference), levels, axis=1
            ) / float(levels[-1] - levels[0])

        conditional_difference = np.full(
            (len(label_masks), len(levels)), np.nan, dtype=float
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for column in range(len(levels)):
                values = matrix[:, column]
                sca_median = np.nanmedian(
                    np.where(label_masks, values[None, :], np.nan), axis=1
                )
                wt_median = np.nanmedian(
                    np.where(~label_masks, values[None, :], np.nan), axis=1
                )
                conditional_difference[:, column] = sca_median - wt_median
        valid_rows = np.isfinite(conditional_difference).all(axis=1)
        if len(levels) == 1:
            conditional_absolute_all = np.abs(conditional_difference[:, 0])
        else:
            conditional_absolute_all = np.trapz(
                np.abs(conditional_difference), levels, axis=1
            ) / float(levels[-1] - levels[0])
        availability_null = availability_absolute_all[valid_rows]
        conditional_null = conditional_absolute_all[valid_rows]
        valid = len(availability_null)
        exact = permutation_mode == "exact"
        if valid == 0:
            comparison_rows.append({
                "inference_role": role, "feature": feature,
                "family": FEATURE_FAMILIES.get(feature, "other"),
                "n_WT_total": n_wt, "n_SCA3_total": n_sca,
                "n_WT_any_defined": int(np.sum(availability[labels == "WT"].any(axis=1))),
                "n_SCA3_any_defined": int(np.sum(availability[labels == "SCA3"].any(axis=1))),
                "current_levels_pA": level_label,
                "permutation_mode": permutation_mode, "valid_labelings": 0,
                "status": "no_valid_labelings",
            })
            continue

        availability_p_values = _survival_probabilities(availability_null, exact)
        conditional_p_values = _survival_probabilities(conditional_null, exact)
        if exact:
            availability_p = float(np.mean(
                availability_null >= observed_availability[1] - 1e-12
            ))
            conditional_p = float(np.mean(
                conditional_null >= observed_conditional[1] - 1e-12
            ))
        else:
            availability_p = float((
                np.sum(availability_null >= observed_availability[1]) + 1
            ) / (valid + 1))
            conditional_p = float((
                np.sum(conditional_null >= observed_conditional[1]) + 1
            ) / (valid + 1))
        combined_null = -2.0 * (
            np.log(availability_p_values) + np.log(conditional_p_values)
        )
        observed_combined = -2.0 * (
            np.log(availability_p) + np.log(conditional_p)
        )
        if exact:
            omnibus_p = float(np.mean(combined_null >= observed_combined - 1e-12))
        else:
            omnibus_p = float((
                np.sum(combined_null >= observed_combined) + 1
            ) / (valid + 1))
        comparison_rows.append({
            "inference_role": role, "feature": feature,
            "family": FEATURE_FAMILIES.get(feature, "other"),
            "n_WT_total": n_wt, "n_SCA3_total": n_sca,
            "n_WT_any_defined": int(np.sum(availability[labels == "WT"].any(axis=1))),
            "n_SCA3_any_defined": int(np.sum(availability[labels == "SCA3"].any(axis=1))),
            "current_levels_pA": level_label,
            "availability_signed_curve_difference_SCA3_minus_WT": observed_availability[0],
            "availability_absolute_curve_difference": observed_availability[1],
            "availability_maximum_absolute_difference": observed_availability[2],
            "availability_permutation_p": availability_p,
            "conditional_signed_curve_difference_SCA3_minus_WT": observed_conditional[0],
            "conditional_absolute_curve_difference": observed_conditional[1],
            "conditional_maximum_absolute_difference": observed_conditional[2],
            "conditional_permutation_p": conditional_p,
            "combined_fisher_statistic": observed_combined,
            "permutation_mode": permutation_mode,
            "valid_labelings": int(valid),
            "permutation_p": omnibus_p,
            "status": "tested",
        })

    comparisons = pd.DataFrame(comparison_rows, columns=comparison_columns)
    adjusted = []
    for role, part in comparisons.groupby("inference_role", sort=False):
        adjusted.append(_add_fdr_columns(part))
    comparisons = pd.concat(adjusted, ignore_index=True) if adjusted else comparisons
    return (
        comparisons,
        pd.DataFrame(coverage_rows, columns=coverage_columns),
        pd.DataFrame(cell_rows, columns=cell_columns),
    )


def cell_scalar_phenotypes(manifest: pd.DataFrame, rheobase: pd.DataFrame) -> pd.DataFrame:
    capacitance_columns = [
        name for name in ["capacitance_pF", "capacitance_10ms_pF", "capacitance_20ms_pF", "capacitance_50ms_pF"]
        if name in manifest
    ]
    cells = manifest[["group", "cell_id", *capacitance_columns]].drop_duplicates(["group", "cell_id"])
    out = cells.merge(rheobase, on=["group", "cell_id"], how="left", validate="one_to_one")
    for capacitance in capacitance_columns:
        suffix = capacitance.replace("capacitance_", "").replace("capacitance", "Cm")
        out[f"rheobase_midpoint_pA_per_{suffix}"] = out["rheobase_midpoint_pA"] / out[capacitance]
    return out


def compare_cell_scalars(cells: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    excluded = {"group", "cell_id"}
    features = [name for name in cells.columns if name not in excluded]
    rng = np.random.default_rng(int(cfg["seed"]) + 3)
    iterations = int(cfg["statistics"]["permutation_iterations"])
    minimum = int(cfg["statistics"]["minimum_cells_per_group"])
    rows = []
    for feature in features:
        a = cells.loc[cells.group == "WT", feature].dropna().to_numpy(float)
        b = cells.loc[cells.group == "SCA3", feature].dropna().to_numpy(float)
        if min(len(a), len(b)) < minimum:
            continue
        lo, hi = _bootstrap_difference(a, b, int(cfg["statistics"]["bootstrap_iterations"]), rng)
        family = "capacitance" if feature.startswith("capacitance") else "rheobase"
        p_value, permutation_mode, valid_labelings = _permutation_p(
            a, b, iterations, rng,
            int(cfg["statistics"].get("exact_max_labelings", 200000)),
        )
        rows.append({
            "feature": feature, "family": family, "n_WT": len(a), "n_SCA3": len(b),
            "median_WT": float(np.median(a)), "median_SCA3": float(np.median(b)),
            "median_difference_SCA3_minus_WT": float(np.median(b) - np.median(a)),
            "difference_ci95_low": float(lo), "difference_ci95_high": float(hi),
            "cliffs_delta_SCA3_vs_WT": cliffs_delta(b, a),
            "permutation_mode": permutation_mode,
            "valid_labelings": valid_labelings,
            "permutation_p": p_value,
        })
    return _add_fdr_columns(pd.DataFrame(rows))
