from __future__ import annotations

from math import comb
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .statistics import _exact_group_masks


DOMAIN_ORDER = ["structure", "spike_timing", "predictive_dynamics"]
DOMAIN_LABELS = {
    "structure": "Structure: low log(Cm)",
    "spike_timing": "Spike timing: long ISI",
    "predictive_dynamics": "Dynamics: low predictive information",
}


def _robust_scale(values, minimum_scale):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, "unavailable"
    center = float(np.median(values))
    mad_scale = float(1.4826 * np.median(np.abs(values - center)))
    if np.isfinite(mad_scale) and mad_scale >= minimum_scale:
        return center, mad_scale, "MAD"
    if len(values) >= 2:
        std_scale = float(np.std(values, ddof=1))
        if np.isfinite(std_scale) and std_scale >= minimum_scale:
            return center, std_scale, "SD_fallback"
    return center, float(minimum_scale), "minimum_scale_fallback"


def _prepare_matrix(features, cell_scalars, cfg):
    dcfg = cfg["disease_coordinate"]
    levels = np.asarray(dcfg["currents_pA"], float)
    capacitance = str(dcfg["capacitance_feature"])
    cells = cell_scalars[["group", "cell_id"]].drop_duplicates().copy()
    animal = features[["group", "cell_id", "animal_id"]].drop_duplicates(
        ["group", "cell_id"]
    )
    cells = cells.merge(animal, on=["group", "cell_id"], how="left")
    cells = cells.sort_values(["group", "cell_id"]).reset_index(drop=True)
    cell_index = pd.MultiIndex.from_frame(cells[["group", "cell_id"]])

    scalar = cell_scalars.set_index(["group", "cell_id"]).reindex(cell_index)
    cap_values = scalar[capacitance].to_numpy(float)
    if np.any(np.isfinite(cap_values) & (cap_values <= 0)):
        raise ValueError("Disease coordinate requires strictly positive capacitance")

    arrays = [np.log(cap_values)]
    metadata = [{
        "feature": "log_" + capacitance,
        "source_feature": capacitance,
        "domain": "structure",
        "current_pA": np.nan,
        "orientation": -1.0,
        "orientation_label": "lower_is_SCA3_like",
    }]
    accepted = features[features["qc_pass"] & features["current_pA"].isin(levels)].copy()
    dynamic_specs = [
        ("mean_isi_ms", "spike_timing", 1.0, "higher_is_SCA3_like", False),
        (
            "predictive_information_nats", "predictive_dynamics", -1.0,
            "lower_is_SCA3_like", True,
        ),
    ]
    for feature, domain, orientation, label, requires_thermo in dynamic_specs:
        source = accepted.copy()
        if requires_thermo:
            source.loc[~source["thermo_eligible"], feature] = np.nan
        table = source.pivot_table(
            index=["group", "cell_id"], columns="current_pA",
            values=feature, aggfunc="median", dropna=False,
        ).reindex(index=cell_index, columns=levels)
        for column, current in enumerate(levels):
            arrays.append(table.iloc[:, column].to_numpy(float))
            metadata.append({
                "feature": f"{feature}__{current:g}pA",
                "source_feature": feature,
                "domain": domain,
                "current_pA": float(current),
                "orientation": orientation,
                "orientation_label": label,
            })
    return cells, np.column_stack(arrays), pd.DataFrame(metadata)


def _score_from_reference(values, labels, metadata, cfg):
    dcfg = cfg["disease_coordinate"]
    minimum_scale = float(dcfg["minimum_scale"])
    z_clip = float(dcfg["z_clip"])
    wt_mask = labels == "WT"
    z_values = np.full(values.shape, np.nan, dtype=float)
    reference_rows = []
    for column, meta in metadata.iterrows():
        center, scale, method = _robust_scale(values[wt_mask, column], minimum_scale)
        if np.isfinite(center) and np.isfinite(scale):
            z_values[:, column] = np.clip(
                float(meta["orientation"]) * (values[:, column] - center) / scale,
                -z_clip, z_clip,
            )
        reference_rows.append({
            **meta.to_dict(),
            "WT_defined_cells": int(np.isfinite(values[wt_mask, column]).sum()),
            "WT_center": center,
            "WT_scale": scale,
            "scale_method": method,
            "z_clip": z_clip,
        })
    domain_scores = np.full((len(values), len(DOMAIN_ORDER)), np.nan, dtype=float)
    for domain_index, domain in enumerate(DOMAIN_ORDER):
        columns = np.flatnonzero(metadata["domain"].to_numpy() == domain)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            domain_scores[:, domain_index] = np.nanmedian(z_values[:, columns], axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        burden = np.nanmean(domain_scores, axis=1)
    return z_values, domain_scores, burden, pd.DataFrame(reference_rows)


def _permutation_statistics(values, metadata, observed_labels, cfg):
    dcfg = cfg["disease_coordinate"]
    exact_max = int(dcfg["exact_max_labelings"])
    iterations = int(dcfg["permutation_iterations"])
    minimum_scale = float(dcfg["minimum_scale"])
    z_clip = float(dcfg["z_clip"])
    orientation = metadata["orientation"].to_numpy(float)
    domain_columns = [
        np.flatnonzero(metadata["domain"].to_numpy() == domain)
        for domain in DOMAIN_ORDER
    ]
    n_sca = int(np.sum(observed_labels == "SCA3"))
    total = comb(len(observed_labels), n_sca)
    if total <= exact_max:
        masks = _exact_group_masks(len(observed_labels), n_sca)
        mode = "exact"
    else:
        rng = np.random.default_rng(int(cfg["seed"]) + 31)
        masks = np.zeros((iterations, len(observed_labels)), dtype=bool)
        for row in range(iterations):
            masks[row, rng.choice(len(observed_labels), n_sca, replace=False)] = True
        mode = "monte_carlo"

    null_statistics = np.full(len(masks), np.nan, dtype=float)
    chunk_size = int(dcfg.get("permutation_chunk_size", 256))
    for start in range(0, len(masks), chunk_size):
        stop = min(start + chunk_size, len(masks))
        sca_masks = masks[start:stop]
        wt_masks = ~sca_masks
        expanded = np.where(wt_masks[:, :, None], values[None, :, :], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            centers = np.nanmedian(expanded, axis=1)
            mad = 1.4826 * np.nanmedian(
                np.abs(expanded - centers[:, None, :]), axis=1
            )
            std = np.nanstd(expanded, axis=1, ddof=1)
        scales = np.where(np.isfinite(mad) & (mad >= minimum_scale), mad, std)
        scales = np.where(
            np.isfinite(scales) & (scales >= minimum_scale), scales, minimum_scale
        )
        z = orientation[None, None, :] * (
            values[None, :, :] - centers[:, None, :]
        ) / scales[:, None, :]
        z = np.clip(z, -z_clip, z_clip)
        domain_scores = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for columns in domain_columns:
                domain_scores.append(np.nanmedian(z[:, :, columns], axis=2))
            burden = np.nanmean(np.stack(domain_scores, axis=2), axis=2)
            sca_median = np.nanmedian(
                np.where(sca_masks, burden, np.nan), axis=1
            )
            wt_median = np.nanmedian(
                np.where(wt_masks, burden, np.nan), axis=1
            )
        null_statistics[start:stop] = sca_median - wt_median

    actual_mask = observed_labels == "SCA3"
    actual_rows = np.flatnonzero(np.all(masks == actual_mask[None, :], axis=1))
    if len(actual_rows) != 1:
        raise RuntimeError("Observed group assignment is absent from permutation set")
    observed = float(null_statistics[actual_rows[0]])
    valid = np.isfinite(null_statistics)
    null = null_statistics[valid]
    if mode == "exact":
        p_one = float(np.mean(null >= observed - 1e-12))
        p_two = float(np.mean(np.abs(null) >= abs(observed) - 1e-12))
    else:
        p_one = float((np.sum(null >= observed) + 1) / (len(null) + 1))
        p_two = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))
    return observed, p_one, p_two, mode, int(valid.sum())


def _auc(scores, labels):
    wt = scores[labels == "WT"]
    sca = scores[labels == "SCA3"]
    wt = wt[np.isfinite(wt)]
    sca = sca[np.isfinite(sca)]
    if len(wt) == 0 or len(sca) == 0:
        return np.nan
    return float(np.mean(sca[:, None] > wt[None, :]) + 0.5 * np.mean(
        sca[:, None] == wt[None, :]
    ))


def _nanmean_rows(values):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=1)


def _score_target_against_reference(target, reference_values, metadata, cfg):
    """Score one held-out cell and the reference rows using the same WT transform."""
    dcfg = cfg["disease_coordinate"]
    minimum_scale = float(dcfg["minimum_scale"])
    z_clip = float(dcfg["z_clip"])
    target_z = np.full(len(metadata), np.nan, dtype=float)
    reference_z = np.full(reference_values.shape, np.nan, dtype=float)
    for column, meta in metadata.iterrows():
        center, scale, _ = _robust_scale(reference_values[:, column], minimum_scale)
        if not (np.isfinite(center) and np.isfinite(scale)):
            continue
        orientation = float(meta["orientation"])
        reference_z[:, column] = np.clip(
            orientation * (reference_values[:, column] - center) / scale,
            -z_clip, z_clip,
        )
        if np.isfinite(target[column]):
            target_z[column] = np.clip(
                orientation * (target[column] - center) / scale,
                -z_clip, z_clip,
            )

    target_domains = np.full(len(DOMAIN_ORDER), np.nan, dtype=float)
    reference_domains = np.full(
        (len(reference_values), len(DOMAIN_ORDER)), np.nan, dtype=float
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for domain_index, domain in enumerate(DOMAIN_ORDER):
            columns = np.flatnonzero(metadata["domain"].to_numpy() == domain)
            target_domains[domain_index] = np.nanmedian(target_z[columns])
            reference_domains[:, domain_index] = np.nanmedian(
                reference_z[:, columns], axis=1
            )
        target_burden = float(np.nanmean(target_domains))
        reference_burden = np.nanmean(reference_domains, axis=1)
    return target_domains, target_burden, reference_domains, reference_burden


def _cross_fitted_scores(values, labels, metadata, cfg, full_burden):
    """Leave each WT cell out of its own reference; retain full-WT scoring for SCA3."""
    wt_indices = np.flatnonzero(labels == "WT")
    crossfit = np.asarray(full_burden, float).copy()
    for target_index in wt_indices:
        reference_indices = wt_indices[wt_indices != target_index]
        if len(reference_indices) == 0:
            crossfit[target_index] = np.nan
            continue
        _, target_burden, _, _ = _score_target_against_reference(
            values[target_index], values[reference_indices], metadata, cfg
        )
        crossfit[target_index] = target_burden
    return crossfit


def _evidence_grades(dynamic_counts, expected_dynamic_values):
    dynamic_counts = np.asarray(dynamic_counts, int)
    high_cutoff = int(np.ceil(0.60 * expected_dynamic_values))
    return np.select(
        [
            dynamic_counts == expected_dynamic_values,
            dynamic_counts >= high_cutoff,
            dynamic_counts >= 1,
        ],
        ["full_dynamic", "partial_high", "sparse_dynamic"],
        default="structural_only",
    )


def _bootstrap_cell_stability(values, labels, metadata, cfg):
    """Reference-resampling stability with leave-one-WT-out scoring for WT cells."""
    dcfg = cfg["disease_coordinate"]
    n_iter = int(dcfg.get(
        "stability_bootstrap_iterations",
        cfg["statistics"]["bootstrap_iterations"],
    ))
    threshold = float(dcfg["robust_z_threshold"])
    minimum_scale = float(dcfg["minimum_scale"])
    consensus_min = int(dcfg["consensus_min_domains"])
    rng = np.random.default_rng(int(cfg["seed"]) + 33)
    wt_indices = np.flatnonzero(labels == "WT")
    orientation = metadata["orientation"].to_numpy(float)
    domain_columns = [
        np.flatnonzero(metadata["domain"].to_numpy() == domain)
        for domain in DOMAIN_ORDER
    ]

    def scale_rows(matrix):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            center = np.nanmedian(matrix, axis=1)
            mad = 1.4826 * np.nanmedian(
                np.abs(matrix - center[:, None]), axis=1
            )
            std = np.nanstd(matrix, axis=1, ddof=1)
        scale = np.where(np.isfinite(mad) & (mad >= minimum_scale), mad, std)
        scale = np.where(
            np.isfinite(scale) & (scale >= minimum_scale), scale, minimum_scale
        )
        return center, scale

    rows = []
    for target_index, label in enumerate(labels):
        if label == "WT":
            pool = wt_indices[wt_indices != target_index]
            mode = "leave_one_WT_out_bootstrap"
        else:
            pool = wt_indices
            mode = "full_WT_bootstrap"
        sampled = rng.choice(pool, size=(n_iter, len(pool)), replace=True)
        reference_values = values[sampled]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            centers = np.nanmedian(reference_values, axis=1)
            mad = 1.4826 * np.nanmedian(
                np.abs(reference_values - centers[:, None, :]), axis=1
            )
            std = np.nanstd(reference_values, axis=1, ddof=1)
        scales = np.where(
            np.isfinite(mad) & (mad >= minimum_scale), mad, std
        )
        scales = np.where(
            np.isfinite(scales) & (scales >= minimum_scale),
            scales, minimum_scale,
        )
        target_z = orientation[None, :] * (
            values[target_index][None, :] - centers
        ) / scales
        reference_z = orientation[None, None, :] * (
            reference_values - centers[:, None, :]
        ) / scales[:, None, :]
        target_z = np.clip(target_z, -float(dcfg["z_clip"]), float(dcfg["z_clip"]))
        reference_z = np.clip(
            reference_z, -float(dcfg["z_clip"]), float(dcfg["z_clip"])
        )
        target_domains = np.full((n_iter, len(DOMAIN_ORDER)), np.nan)
        reference_domains = np.full(
            (n_iter, len(pool), len(DOMAIN_ORDER)), np.nan
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for domain_index, columns in enumerate(domain_columns):
                target_domains[:, domain_index] = np.nanmedian(
                    target_z[:, columns], axis=1
                )
                reference_domains[:, :, domain_index] = np.nanmedian(
                    reference_z[:, :, columns], axis=2
                )
            target_burden = np.nanmean(target_domains, axis=1)
            reference_burden = np.nanmean(reference_domains, axis=2)
            observed_maximum = np.nanmax(reference_burden, axis=1)

        composite_center, composite_scale = scale_rows(reference_burden)
        composite_boundary = composite_center + threshold * composite_scale
        robust_outside = target_burden > composite_boundary
        domain_hits = np.zeros(n_iter, dtype=int)
        for domain_index in range(len(DOMAIN_ORDER)):
            center, scale = scale_rows(reference_domains[:, :, domain_index])
            domain_hits += (
                np.isfinite(target_domains[:, domain_index])
                & np.isfinite(center)
                & (target_domains[:, domain_index] > center + threshold * scale)
            )
        valid_mask = (
            np.isfinite(target_burden)
            & np.isfinite(composite_center)
            & np.isfinite(observed_maximum)
        )
        valid = int(valid_mask.sum())
        robust_hits = int(np.sum(robust_outside & valid_mask))
        envelope_hits = int(np.sum(
            (target_burden > observed_maximum) & valid_mask
        ))
        consensus_hits = int(np.sum(
            robust_outside & (domain_hits >= consensus_min) & valid_mask
        ))
        denominator = float(valid) if valid else np.nan
        rows.append({
            "stability_reference_mode": mode,
            "stability_bootstrap_iterations_requested": n_iter,
            "stability_bootstrap_iterations_valid": valid,
            "bootstrap_p_outside_WT_robust_boundary": robust_hits / denominator,
            "bootstrap_p_outside_observed_WT_envelope": envelope_hits / denominator,
            "bootstrap_p_WT_exit_consensus": consensus_hits / denominator,
        })
    return pd.DataFrame(rows)


def build_disease_coordinate(features, cell_scalars, cfg, compute_stability=True):
    """Build a model-free WT-to-SCA3 endpoint coordinate without imputation."""
    cells, values, metadata = _prepare_matrix(features, cell_scalars, cfg)
    labels = cells["group"].astype(str).to_numpy()
    z_values, domain_values, burden, reference = _score_from_reference(
        values, labels, metadata, cfg
    )
    dcfg = cfg["disease_coordinate"]
    threshold = float(dcfg["robust_z_threshold"])
    wt_anchor = float(np.median(burden[labels == "WT"]))
    sca_anchor = float(np.median(burden[labels == "SCA3"]))
    anchor_span = sca_anchor - wt_anchor
    if not np.isfinite(anchor_span) or anchor_span <= 0:
        raise ValueError("Disease coordinate requires SCA3 burden above WT burden")
    q = (burden - wt_anchor) / anchor_span

    scores = cells.copy()
    for index, domain in enumerate(DOMAIN_ORDER):
        scores[f"domain_z__{domain}"] = domain_values[:, index]
    scores["disease_burden_z"] = burden
    scores["q_endpoint_unbounded"] = q
    scores["q_endpoint_clipped_0_1"] = np.clip(q, 0.0, 1.0)
    scores["domains_available"] = np.isfinite(domain_values).sum(axis=1)
    scores["dynamic_values_available"] = np.isfinite(values[:, 1:]).sum(axis=1)
    expected_dynamic_values = int(values.shape[1] - 1)
    scores["dynamic_values_expected"] = expected_dynamic_values
    scores["dynamic_coverage_fraction"] = (
        scores["dynamic_values_available"] / float(expected_dynamic_values)
    )
    scores["evidence_grade"] = _evidence_grades(
        scores["dynamic_values_available"].to_numpy(), expected_dynamic_values
    )
    # Kept as a compatibility alias, but the old domain-count categories are retired.
    scores["coordinate_reliability"] = scores["evidence_grade"]

    crossfit_burden = _cross_fitted_scores(
        values, labels, metadata, cfg, burden
    )
    crossfit_wt_anchor = float(np.nanmedian(crossfit_burden[labels == "WT"]))
    crossfit_sca_anchor = float(np.nanmedian(crossfit_burden[labels == "SCA3"]))
    crossfit_span = crossfit_sca_anchor - crossfit_wt_anchor
    scores["crossfit_disease_burden_z"] = crossfit_burden
    scores["crossfit_q_endpoint_unbounded"] = (
        (crossfit_burden - crossfit_wt_anchor) / crossfit_span
        if np.isfinite(crossfit_span) and crossfit_span > 0 else np.nan
    )
    scores["crossfit_score_mode"] = np.where(
        labels == "WT", "leave_one_WT_out", "full_WT_reference"
    )

    domain_boundaries = {}
    for index, domain in enumerate(DOMAIN_ORDER):
        center, scale, _ = _robust_scale(
            domain_values[labels == "WT", index], float(dcfg["minimum_scale"])
        )
        boundary = center + threshold * scale
        domain_boundaries[domain] = boundary
        scores[f"outside_WT__{domain}"] = domain_values[:, index] > boundary
    scores["domains_outside_WT"] = scores[
        [f"outside_WT__{domain}" for domain in DOMAIN_ORDER]
    ].sum(axis=1)
    composite_center, composite_scale, composite_method = _robust_scale(
        burden[labels == "WT"], float(dcfg["minimum_scale"])
    )
    composite_boundary = composite_center + threshold * composite_scale
    observed_wt_max = float(np.nanmax(burden[labels == "WT"]))
    scores["outside_WT_robust_boundary"] = burden > composite_boundary
    scores["outside_observed_WT_envelope"] = burden > observed_wt_max
    scores["WT_exit_consensus_marker"] = (
        scores["outside_WT_robust_boundary"]
        & (scores["domains_outside_WT"] >= int(dcfg["consensus_min_domains"]))
    )

    if compute_stability:
        stability = _bootstrap_cell_stability(values, labels, metadata, cfg)
        scores = pd.concat(
            [scores.reset_index(drop=True), stability.reset_index(drop=True)], axis=1
        )
    else:
        scores["stability_reference_mode"] = "not_computed"
        scores["stability_bootstrap_iterations_requested"] = 0
        scores["stability_bootstrap_iterations_valid"] = 0
        scores["bootstrap_p_outside_WT_robust_boundary"] = np.nan
        scores["bootstrap_p_outside_observed_WT_envelope"] = np.nan
        scores["bootstrap_p_WT_exit_consensus"] = np.nan

    feature_rows = []
    for cell_index, cell in cells.iterrows():
        for column, meta in metadata.iterrows():
            feature_rows.append({
                "group": cell["group"], "cell_id": cell["cell_id"],
                "animal_id": cell["animal_id"], **meta.to_dict(),
                "raw_value": values[cell_index, column],
                "oriented_WT_z": z_values[cell_index, column],
                "available": bool(np.isfinite(values[cell_index, column])),
            })
    feature_scores = pd.DataFrame(feature_rows)

    reference = reference.assign(level="feature")
    extra_reference = []
    for index, domain in enumerate(DOMAIN_ORDER):
        center, scale, method = _robust_scale(
            domain_values[labels == "WT", index], float(dcfg["minimum_scale"])
        )
        extra_reference.append({
            "level": "domain", "feature": domain, "source_feature": domain,
            "domain": domain, "current_pA": np.nan, "orientation": 1.0,
            "orientation_label": "higher_is_SCA3_like",
            "WT_defined_cells": int(np.isfinite(domain_values[labels == "WT", index]).sum()),
            "WT_center": center, "WT_scale": scale, "scale_method": method,
            "z_clip": float(dcfg["z_clip"]),
            "WT_robust_upper_boundary": domain_boundaries[domain],
        })
    extra_reference.append({
        "level": "composite", "feature": "disease_burden_z",
        "source_feature": "equal_domain_mean", "domain": "composite",
        "current_pA": np.nan, "orientation": 1.0,
        "orientation_label": "higher_is_SCA3_like",
        "WT_defined_cells": int(np.isfinite(burden[labels == "WT"]).sum()),
        "WT_center": composite_center, "WT_scale": composite_scale,
        "scale_method": composite_method, "z_clip": float(dcfg["z_clip"]),
        "WT_robust_upper_boundary": composite_boundary,
        "WT_observed_maximum": observed_wt_max,
        "WT_endpoint_anchor": wt_anchor, "SCA3_endpoint_anchor": sca_anchor,
    })
    reference = pd.concat([reference, pd.DataFrame(extra_reference)], ignore_index=True)

    observed, p_one, p_two, mode, valid = _permutation_statistics(
        values, metadata, labels, cfg
    )
    rng = np.random.default_rng(int(cfg["seed"]) + 32)
    bootstrap = []
    wt_scores = burden[labels == "WT"]
    sca_scores = burden[labels == "SCA3"]
    for _ in range(int(cfg["statistics"]["bootstrap_iterations"])):
        bootstrap.append(
            np.median(rng.choice(sca_scores, len(sca_scores), replace=True))
            - np.median(rng.choice(wt_scores, len(wt_scores), replace=True))
        )
    ci_low, ci_high = np.quantile(bootstrap, [0.025, 0.975])
    validation = pd.DataFrame([{
        "n_WT": int(np.sum(labels == "WT")),
        "n_SCA3": int(np.sum(labels == "SCA3")),
        "WT_median_burden_z": wt_anchor,
        "SCA3_median_burden_z": sca_anchor,
        "median_difference_SCA3_minus_WT": sca_anchor - wt_anchor,
        "difference_ci95_low": float(ci_low),
        "difference_ci95_high": float(ci_high),
        "descriptive_auc_SCA3_vs_WT": _auc(burden, labels),
        "descriptive_auc_scope": "internal_training_cohort",
        "internal_cross_fitted_auc_SCA3_vs_WT": _auc(crossfit_burden, labels),
        "cross_fitted_auc_scope": "leave_one_WT_out_internal_not_external",
        "crossfit_WT_median_burden_z": crossfit_wt_anchor,
        "crossfit_SCA3_median_burden_z": crossfit_sca_anchor,
        "crossfit_min_SCA3_minus_max_WT_margin": float(
            np.nanmin(crossfit_burden[labels == "SCA3"])
            - np.nanmax(crossfit_burden[labels == "WT"])
        ),
        "descriptive_auc_structure_only": _auc(domain_values[:, 0], labels),
        "descriptive_auc_spike_timing_only": _auc(domain_values[:, 1], labels),
        "descriptive_auc_predictive_dynamics_only": _auc(domain_values[:, 2], labels),
        "descriptive_auc_without_structure": _auc(
            _nanmean_rows(domain_values[:, 1:]), labels
        ),
        "descriptive_auc_without_spike_timing": _auc(
            _nanmean_rows(domain_values[:, [0, 2]]), labels
        ),
        "descriptive_auc_without_predictive_dynamics": _auc(
            _nanmean_rows(domain_values[:, :2]), labels
        ),
        "permutation_mode": mode,
        "valid_labelings": valid,
        "permutation_p_one_sided_predeclared_direction": p_one,
        "permutation_p_two_sided": p_two,
        "permutation_observed_recomputed_difference": observed,
        "WT_robust_upper_boundary": composite_boundary,
        "WT_observed_maximum": observed_wt_max,
        "SCA3_above_WT_robust_boundary": int(np.sum(
            scores.loc[scores["group"] == "SCA3", "outside_WT_robust_boundary"]
        )),
        "SCA3_above_observed_WT_envelope": int(np.sum(
            scores.loc[scores["group"] == "SCA3", "outside_observed_WT_envelope"]
        )),
        "SCA3_WT_exit_consensus": int(np.sum(
            scores.loc[scores["group"] == "SCA3", "WT_exit_consensus_marker"]
        )),
    }])
    return scores, reference, validation, feature_scores


def plot_disease_coordinate(scores, output_path):
    ordered = scores.sort_values("disease_burden_z").reset_index(drop=True)
    colors = ordered["group"].map({"WT": "#2b6cb0", "SCA3": "#c43c35"})
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    axes[0].scatter(
        ordered["q_endpoint_unbounded"], np.arange(len(ordered)),
        c=colors, s=55, edgecolor="white", linewidth=0.7,
    )
    axes[0].axvline(0, color="#2b6cb0", linestyle="--", linewidth=1.2, label="WT median")
    axes[0].axvline(1, color="#c43c35", linestyle="--", linewidth=1.2, label="SCA3 median")
    axes[0].set_yticks(np.arange(len(ordered)))
    axes[0].set_yticklabels(ordered["cell_id"])
    axes[0].set_xlabel("Endpoint disease coordinate q (unbounded)")
    axes[0].set_title("Cell-level WT to SCA3 coordinate")
    axes[0].grid(axis="x", alpha=0.25)
    axes[0].legend(frameon=False)

    matrix = ordered[[f"domain_z__{domain}" for domain in DOMAIN_ORDER]].to_numpy(float)
    image = axes[1].imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-5, vmax=5)
    axes[1].set_xticks(np.arange(len(DOMAIN_ORDER)))
    axes[1].set_xticklabels(["Structure", "Mean ISI", "Predictive info"], rotation=20, ha="right")
    axes[1].set_yticks(np.arange(len(ordered)))
    axes[1].set_yticklabels(ordered["cell_id"])
    axes[1].set_title("Oriented WT-reference domain scores")
    colorbar = fig.colorbar(image, ax=axes[1], fraction=0.046, pad=0.04)
    colorbar.set_label("Robust z; positive is SCA3-like")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def write_disease_coordinate_report(scores, validation, output_path):
    row = validation.iloc[0]
    animal_unavailable = (
        "animal_id" not in scores
        or scores["animal_id"].isna().any()
        or scores["animal_id"].astype(str).eq("NA_NOT_RECOVERABLE").any()
    )
    lines = [
        "# Model-free WT-to-SCA3 endpoint coordinate",
        "",
        "The coordinate combines equally weighted structure, spike-timing and predictive-dynamics domains. Positive robust z values point toward the observed SCA3 endpoint.",
        "",
        "q=0 is the WT endpoint median and q=1 is the SCA3 endpoint median. q is not a probability, degeneration time or thermodynamic reaction coordinate. Injected current is a perturbation axis, not disease time.",
        "",
        "Missing conditional values are not imputed. `evidence_grade` is based on current-level dynamic coverage rather than only on the number of represented domains.",
        "",
        "## Validation",
        "",
        f"- Cells: {int(row['n_WT'])} WT and {int(row['n_SCA3'])} SCA3.",
        f"- Median burden difference: {row['median_difference_SCA3_minus_WT']:.6g} (95% bootstrap interval {row['difference_ci95_low']:.6g} to {row['difference_ci95_high']:.6g}).",
        f"- Exact recomputed two-sided permutation p: {row['permutation_p_two_sided']:.6g} across {int(row['valid_labelings'])} valid labelings.",
        f"- Descriptive full-coordinate AUC: {row['descriptive_auc_SCA3_vs_WT']:.4f}.",
        f"- Internal leave-one-WT-out AUC: {row['internal_cross_fitted_auc_SCA3_vs_WT']:.4f}; minimum cross-fitted SCA3-minus-WT margin: {row['crossfit_min_SCA3_minus_max_WT_margin']:.6g}.",
        f"- Descriptive structure-only AUC: {row['descriptive_auc_structure_only']:.4f}; without-structure AUC: {row['descriptive_auc_without_structure']:.4f}.",
        f"- SCA3 above the observed WT envelope: {int(row['SCA3_above_observed_WT_envelope'])}/{int(row['n_SCA3'])}.",
        f"- SCA3 meeting the multi-domain WT-exit consensus: {int(row['SCA3_WT_exit_consensus'])}/{int(row['n_SCA3'])}.",
        "",
        "The descriptive AUC uses the cohort that constructs the WT reference. The leave-one-WT-out AUC removes WT self-normalization but remains internal because SCA3 cells and feature definitions are unchanged. Neither value is external or clinical performance.",
    ]
    if animal_unavailable:
        lines.extend([
            "",
            "Animal identifiers are unavailable for at least one cell. The analysis unit is therefore the cell; animal-level independence and generalization cannot be verified.",
        ])
    lines.extend([
        "",
        "## Cell scores",
        "",
        "| Group | Cell | q | Dynamic coverage | Evidence | Above observed WT | Bootstrap P(outside envelope) | WT-exit consensus | Bootstrap P(consensus) |",
        "| --- | --- | ---: | ---: | --- | --- | ---: | --- | ---: |",
    ])
    for cell in scores.sort_values("q_endpoint_unbounded").itertuples(index=False):
        lines.append(
            f"| {cell.group} | {cell.cell_id} | {cell.q_endpoint_unbounded:.4f} | "
            f"{int(cell.dynamic_values_available)}/{int(cell.dynamic_values_expected)} | "
            f"{cell.evidence_grade} | {bool(cell.outside_observed_WT_envelope)} | "
            f"{cell.bootstrap_p_outside_observed_WT_envelope:.3f} | "
            f"{bool(cell.WT_exit_consensus_marker)} | "
            f"{cell.bootstrap_p_WT_exit_consensus:.3f} |"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
