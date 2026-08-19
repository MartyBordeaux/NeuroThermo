import hashlib
import itertools
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "neurothermo-dependency-aware-mpl"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import __version__
from .config import REQUIRED_INPUTS


KEYS = ["group", "cell_id"]


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_value(value):
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _write_json(path, value):
    with Path(path).open("w", encoding="utf-8") as handle:
        json.dump(_json_value(value), handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def _truthy(series):
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin(("true", "1", "yes", "y"))


def _transform(values, mode):
    values = np.asarray(values, dtype=float)
    if mode == "identity":
        return values.copy()
    out = np.full(values.shape, np.nan, dtype=float)
    if mode == "log":
        mask = values > 0
        out[mask] = np.log(values[mask])
        return out
    if mode == "log1p":
        mask = values >= 0
        out[mask] = np.log1p(values[mask])
        return out
    raise ValueError("Unknown transform: {}".format(mode))


def _read_inputs(upstream):
    data = {
        "sweep": pd.read_csv(upstream / "sweep_features.csv"),
        "scalar": pd.read_csv(upstream / "cell_scalar_phenotypes.csv"),
        "integrated": pd.read_csv(upstream / "cell_integrated_phenotypes.csv"),
        "legacy": pd.read_csv(upstream / "disease_coordinate" / "cell_disease_coordinate.csv"),
    }
    for key, frame in data.items():
        animal_columns = [column for column in frame.columns if column.lower() == "animal_id"]
        if animal_columns:
            data[key] = frame.drop(columns=animal_columns)
    return data


def _required_columns(config):
    scalar = {"group", "cell_id"}
    sweep = {"group", "cell_id", "current_pA", "qc_pass", "thermo_eligible"}
    for spec in config["features"].values():
        target = scalar if spec["source"] == "scalar" else sweep
        target.add(spec["source_column"])
    return scalar, sweep


def validate_inputs(config, upstream):
    missing = [str(upstream / item) for item in REQUIRED_INPUTS if not (upstream / item).exists()]
    if missing:
        raise FileNotFoundError("Required inputs are missing:\n" + "\n".join(missing))
    data = _read_inputs(upstream)
    scalar_required, sweep_required = _required_columns(config)
    absent_scalar = sorted(scalar_required.difference(data["scalar"].columns))
    absent_sweep = sorted(sweep_required.difference(data["sweep"].columns))
    if absent_scalar or absent_sweep:
        raise ValueError("Missing columns. scalar={} sweep={}".format(absent_scalar, absent_sweep))
    expected_groups = set(config["analysis"]["groups"])
    found_groups = set(data["scalar"]["group"].dropna().astype(str))
    if found_groups != expected_groups:
        raise ValueError("Expected groups {}, found {}".format(sorted(expected_groups), sorted(found_groups)))
    if data["scalar"].duplicated(KEYS).any():
        raise ValueError("Duplicate group/cell_id rows in cell_scalar_phenotypes.csv")
    scalar_keys = set(map(tuple, data["scalar"][KEYS].astype(str).to_numpy()))
    sweep_keys = set(map(tuple, data["sweep"][KEYS].astype(str).drop_duplicates().to_numpy()))
    if scalar_keys != sweep_keys:
        raise ValueError("Scalar and sweep files do not contain the same cells")
    counts = data["scalar"].groupby("group")["cell_id"].nunique().to_dict()
    checks = [
        {"check": "required_files", "status": "PASS", "detail": str(len(REQUIRED_INPUTS))},
        {"check": "required_columns", "status": "PASS", "detail": "all present"},
        {"check": "cell_key_consistency", "status": "PASS", "detail": str(len(scalar_keys))},
        {"check": "inference_unit", "status": "PASS", "detail": "cell"},
        {"check": "animal_id_used", "status": "PASS", "detail": "false"},
    ]
    return {
        "upstream_dir": str(upstream),
        "n_cells": int(len(scalar_keys)),
        "n_sweeps": int(len(data["sweep"])),
        "group_cell_counts": {str(k): int(v) for k, v in counts.items()},
        "animal_id_used": False,
        "checks": checks,
    }


def build_cell_features(config, data, current_min=None, current_max=None):
    acfg = config["analysis"]
    current_min = float(acfg["high_current_min_pA"] if current_min is None else current_min)
    current_max = float(acfg["high_current_max_pA"] if current_max is None else current_max)
    scalar = data["scalar"].copy()
    scalar["group"] = scalar["group"].astype(str)
    scalar["cell_id"] = scalar["cell_id"].astype(str)
    out = scalar[KEYS].copy()
    for name, spec in config["features"].items():
        if spec["source"] == "scalar":
            values = pd.to_numeric(scalar[spec["source_column"]], errors="coerce")
            out[name] = values
            out[name + "__n_points"] = values.notna().astype(int)

    sweep = data["sweep"].copy()
    sweep["group"] = sweep["group"].astype(str)
    sweep["cell_id"] = sweep["cell_id"].astype(str)
    sweep = sweep[_truthy(sweep["qc_pass"])].copy()
    current = pd.to_numeric(sweep["current_pA"], errors="coerce")
    sweep = sweep[current.between(current_min, current_max)].copy()
    thermo_ok = _truthy(sweep["thermo_eligible"])
    min_points = int(acfg["minimum_conditional_points"])
    for name, spec in config["features"].items():
        if spec["source"] != "sweep":
            continue
        values = pd.to_numeric(sweep[spec["source_column"]], errors="coerce")
        eligible = values.notna()
        if bool(spec.get("conditional", False)):
            eligible &= thermo_ok
        temp = sweep.loc[eligible, KEYS].copy()
        temp["value"] = values.loc[eligible]
        summary = temp.groupby(KEYS)["value"].agg(["median", "count"]).reset_index()
        summary = summary.rename(columns={"median": name, "count": name + "__n_points"})
        if bool(spec.get("conditional", False)):
            summary.loc[summary[name + "__n_points"] < min_points, name] = np.nan
        out = out.merge(summary, on=KEYS, how="left", validate="one_to_one")
        out[name + "__n_points"] = out[name + "__n_points"].fillna(0).astype(int)
    counts = sweep.groupby(KEYS).size().rename("window_sweeps_available").reset_index()
    out = out.merge(counts, on=KEYS, how="left", validate="one_to_one")
    out["window_sweeps_available"] = out["window_sweeps_available"].fillna(0).astype(int)
    out["current_window_min_pA"] = current_min
    out["current_window_max_pA"] = current_max
    return out.sort_values(KEYS).reset_index(drop=True)


def _robust_center_scale(values, minimum_scale, fallback_scale=None):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan, "missing"
    center = float(np.median(values))
    scale = float(1.4826 * np.median(np.abs(values - center)))
    method = "1.4826_MAD"
    if not np.isfinite(scale) or scale < minimum_scale:
        scale = float(np.std(values, ddof=1)) if values.size > 1 else np.nan
        method = "sample_SD_fallback"
    if not np.isfinite(scale) or scale < minimum_scale:
        scale = float(fallback_scale) if fallback_scale is not None else minimum_scale
        method = "fixed_fallback"
    if not np.isfinite(scale) or scale < minimum_scale:
        scale = minimum_scale
    return center, scale, method


def fit_reference(features, labels, feature_names, config, fallback=None):
    labels = np.asarray(labels)
    minimum_scale = float(config["analysis"]["minimum_scale"])
    reference = {}
    for name in feature_names:
        spec = config["features"][name]
        values = _transform(features[name].to_numpy(), spec["transform"])
        fallback_scale = None if fallback is None or name not in fallback else fallback[name]["scale"]
        center, scale, method = _robust_center_scale(values[labels == "WT"], minimum_scale, fallback_scale)
        reference[name] = {
            "center": center,
            "scale": scale,
            "scale_method": method,
            "transform": spec["transform"],
            "direction": int(spec["direction"]),
        }
    return reference


def score_specification(features, reference, specification, config):
    z_clip = float(config["analysis"]["z_clip"])
    feature_scores = {}
    for name in specification:
        spec = config["features"][name]
        values = _transform(features[name].to_numpy(), spec["transform"])
        ref = reference[name]
        z = int(spec["direction"]) * (values - ref["center"]) / ref["scale"]
        feature_scores[name] = np.clip(z, -z_clip, z_clip)
    domains = []
    for domain in specification.values():
        if domain not in domains:
            domains.append(domain)
    scores = pd.DataFrame(index=features.index)
    for name, values in feature_scores.items():
        scores["feature_z__" + name] = values
    for domain in domains:
        names = [name for name, assigned in specification.items() if assigned == domain]
        matrix = np.column_stack([feature_scores[name] for name in names])
        counts = np.isfinite(matrix).sum(axis=1)
        scores["domain_z__" + domain] = np.divide(
            np.nansum(matrix, axis=1), counts,
            out=np.full(len(features), np.nan), where=counts > 0,
        )
    domain_columns = ["domain_z__" + domain for domain in domains]
    domain_matrix = scores[domain_columns].to_numpy(dtype=float)
    counts = np.isfinite(domain_matrix).sum(axis=1)
    scores["fingerprint_burden_z"] = np.divide(
        np.nansum(domain_matrix, axis=1), counts,
        out=np.full(len(features), np.nan), where=counts > 0,
    )
    scores["features_available"] = scores[["feature_z__" + x for x in specification]].notna().sum(axis=1)
    scores["domains_available"] = counts
    return scores, domains


def add_coordinate_and_boundary(features, scores, domains, config):
    out = pd.concat([features[KEYS].reset_index(drop=True), scores.reset_index(drop=True)], axis=1)
    wt = out["group"] == "WT"; sca = out["group"] == "SCA3"
    wt_anchor = float(out.loc[wt, "fingerprint_burden_z"].median())
    sca_anchor = float(out.loc[sca, "fingerprint_burden_z"].median())
    denominator = sca_anchor - wt_anchor
    if not np.isfinite(denominator) or denominator <= float(config["analysis"]["minimum_scale"]):
        raise RuntimeError("SCA3 endpoint is not above WT for this fingerprint specification")
    out["q_endpoint_unbounded"] = (out["fingerprint_burden_z"] - wt_anchor) / denominator
    out["q_endpoint_clipped_0_1"] = out["q_endpoint_unbounded"].clip(0, 1)
    threshold = float(config["analysis"]["robust_z_threshold"])
    domain_columns = ["domain_z__" + d for d in domains]
    out["domains_outside_WT"] = (out[domain_columns] > threshold).sum(axis=1)
    min_domains = int(config["analysis"]["minimum_domains_for_boundary"])
    consensus = int(config["analysis"]["consensus_min_domains"])
    eligible = out["domains_available"] >= min_domains
    out["outside_WT_robust_boundary"] = eligible & (out["fingerprint_burden_z"] > threshold)
    wt_max = float(out.loc[wt, "fingerprint_burden_z"].max())
    out["outside_observed_WT_envelope"] = eligible & (out["fingerprint_burden_z"] > wt_max)
    out["WT_exit_consensus_marker"] = eligible & out["outside_WT_robust_boundary"] & (out["domains_outside_WT"] >= consensus)
    n_domains = len(domains)
    out["evidence_grade"] = np.select(
        [out["domains_available"] == n_domains, out["domains_available"] == n_domains - 1,
         out["domains_available"] >= max(2, n_domains - 2)],
        ["full", "partial_high", "partial"], default="insufficient",
    )
    return out, {
        "wt_median_burden": wt_anchor, "sca3_median_burden": sca_anchor,
        "q_denominator": denominator, "wt_observed_max_burden": wt_max,
        "robust_z_threshold": threshold, "consensus_min_domains": consensus,
        "minimum_domains_for_boundary": min_domains,
    }


def _auc(labels, values):
    labels = np.asarray(labels); values = np.asarray(values, dtype=float)
    a = values[(labels == "SCA3") & np.isfinite(values)]
    b = values[(labels == "WT") & np.isfinite(values)]
    if not len(a) or not len(b):
        return np.nan
    return float(((a[:, None] > b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()) / (len(a) * len(b)))


def _cliffs_delta(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if not len(a) or not len(b):
        return np.nan
    return float(((a[:, None] > b[None, :]).sum() - (a[:, None] < b[None, :]).sum()) / (len(a) * len(b)))


def _bh(values):
    values = np.asarray(values, dtype=float); out = np.full(values.shape, np.nan)
    valid = np.where(np.isfinite(values))[0]
    if not len(valid):
        return out
    order = valid[np.argsort(values[valid])]
    adjusted = values[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[order] = np.clip(adjusted, 0, 1)
    return out


def exact_label_masks(labels, exact_max):
    labels = np.asarray(labels); n = len(labels); n_sca = int(np.sum(labels == "SCA3"))
    total = math.comb(n, n_sca)
    if total > exact_max:
        raise RuntimeError("Exact label set has {} assignments, above exact_max_labelings={}".format(total, exact_max))
    masks = np.zeros((total, n), dtype=bool)
    for row, combo in enumerate(itertools.combinations(range(n), n_sca)):
        masks[row, :] = True; masks[row, list(combo)] = False
    return masks


def permutation_feature_scores(features, feature_names, wt_masks, config):
    minimum_scale = float(config["analysis"]["minimum_scale"])
    z_clip = float(config["analysis"]["z_clip"])
    scores = {}
    for name in feature_names:
        spec = config["features"][name]
        values = _transform(features[name].to_numpy(), spec["transform"])
        reference_values = np.where(wt_masks, values[None, :], np.nan)
        center = np.nanmedian(reference_values, axis=1)
        scale = 1.4826 * np.nanmedian(np.abs(reference_values - center[:, None]), axis=1)
        fallback = np.nanstd(reference_values, axis=1, ddof=1)
        bad = ~np.isfinite(scale) | (scale < minimum_scale)
        scale[bad] = fallback[bad]
        scale[~np.isfinite(scale) | (scale < minimum_scale)] = minimum_scale
        z = int(spec["direction"]) * (values[None, :] - center[:, None]) / scale[:, None]
        scores[name] = np.clip(z, -z_clip, z_clip).astype(np.float32)
    return scores


def permutation_burden(feature_z, specification):
    domains = list(dict.fromkeys(specification.values()))
    domain_arrays = []
    for domain in domains:
        names = [name for name, assigned in specification.items() if assigned == domain]
        matrix = np.stack([feature_z[name] for name in names], axis=0)
        counts = np.isfinite(matrix).sum(axis=0)
        domain_arrays.append(np.divide(
            np.nansum(matrix, axis=0), counts,
            out=np.full(matrix.shape[1:], np.nan, dtype=np.float32), where=counts > 0,
        ))
    matrix = np.stack(domain_arrays, axis=0); counts = np.isfinite(matrix).sum(axis=0)
    return np.divide(
        np.nansum(matrix, axis=0), counts,
        out=np.full(matrix.shape[1:], np.nan, dtype=np.float32), where=counts > 0,
    )


def permutation_validation(feature_z, specification, wt_masks, observed_difference):
    burden = permutation_burden(feature_z, specification)
    statistics = (
        np.nanmedian(np.where(~wt_masks, burden, np.nan), axis=1)
        - np.nanmedian(np.where(wt_masks, burden, np.nan), axis=1)
    )
    valid = statistics[np.isfinite(statistics)]
    return {
        "permutation_mode": "exact", "valid_labelings": int(len(valid)),
        "permutation_p_one_sided": float(np.mean(valid >= observed_difference - 1e-12)),
        "permutation_p_two_sided": float(np.mean(np.abs(valid) >= abs(observed_difference) - 1e-12)),
    }


def bootstrap_difference(values, labels, iterations, rng):
    values = np.asarray(values, dtype=float); labels = np.asarray(labels)
    wt = values[(labels == "WT") & np.isfinite(values)]
    sca = values[(labels == "SCA3") & np.isfinite(values)]
    if not len(wt) or not len(sca):
        return np.nan, np.nan
    draws = np.empty(iterations)
    for i in range(iterations):
        draws[i] = (
            np.median(rng.choice(sca, len(sca), replace=True))
            - np.median(rng.choice(wt, len(wt), replace=True))
        )
    return tuple(np.quantile(draws, [0.025, 0.975]))


def primary_feature_comparisons(features, labels, wt_masks, config):
    labels = np.asarray(labels)
    rng = np.random.default_rng(int(config["analysis"]["random_seed"]) + 11)
    rows = []
    for name, spec in config["features"].items():
        raw = pd.to_numeric(features[name], errors="coerce").to_numpy(dtype=float)
        direction = int(spec["direction"])
        oriented = _transform(raw, spec["transform"])
        if direction != 0:
            oriented = direction * oriented
        wt = labels == "WT"; sca = labels == "SCA3"
        observed = np.nanmedian(oriented[sca]) - np.nanmedian(oriented[wt])
        perm = (
            np.nanmedian(np.where(~wt_masks, oriented[None, :], np.nan), axis=1)
            - np.nanmedian(np.where(wt_masks, oriented[None, :], np.nan), axis=1)
        )
        valid = perm[np.isfinite(perm)]
        ci = bootstrap_difference(raw, labels, int(config["analysis"]["bootstrap_iterations"]), rng)
        rows.append({
            "feature": name, "role": spec.get("role", "core"),
            "predeclared_SCA3_direction": (
                "higher" if direction > 0 else "lower" if direction < 0 else "none_diagnostic"
            ),
            "n_WT": int(np.isfinite(raw[wt]).sum()), "n_SCA3": int(np.isfinite(raw[sca]).sum()),
            "median_WT_raw": float(np.nanmedian(raw[wt])), "median_SCA3_raw": float(np.nanmedian(raw[sca])),
            "median_difference_SCA3_minus_WT_raw": float(np.nanmedian(raw[sca]) - np.nanmedian(raw[wt])),
            "raw_difference_ci95_low": float(ci[0]), "raw_difference_ci95_high": float(ci[1]),
            "oriented_transformed_median_difference": float(observed),
            "oriented_cliffs_delta": _cliffs_delta(oriented[sca], oriented[wt]),
            "valid_labelings": int(len(valid)),
            "permutation_p_one_sided": (
                float(np.mean(valid >= observed - 1e-12)) if direction != 0 else np.nan
            ),
            "permutation_p_two_sided": float(np.mean(np.abs(valid) >= abs(observed) - 1e-12)),
        })
    out = pd.DataFrame(rows)
    out["fdr_q_two_sided_all_features"] = _bh(out["permutation_p_two_sided"])
    out["fdr_q_two_sided_within_role"] = np.nan
    for _, idx in out.groupby("role").groups.items():
        out.loc[idx, "fdr_q_two_sided_within_role"] = _bh(out.loc[idx, "permutation_p_two_sided"].to_numpy())
    return out


def observed_validation(name, features, scores, labels, feature_z, specification, wt_masks, config):
    labels = np.asarray(labels); burden = scores["fingerprint_burden_z"].to_numpy(dtype=float)
    observed = float(np.nanmedian(burden[labels == "SCA3"]) - np.nanmedian(burden[labels == "WT"]))
    permutation = permutation_validation(feature_z, specification, wt_masks, observed)
    rng = np.random.default_rng(int(config["analysis"]["random_seed"]) + sum(map(ord, name)))
    ci = bootstrap_difference(burden, labels, int(config["analysis"]["bootstrap_iterations"]), rng)
    row = {
        "fingerprint": name, "n_features": len(specification),
        "n_domains": len(set(specification.values())),
        "n_complete_cells": int(features[list(specification)].notna().all(axis=1).sum()),
        "median_difference_SCA3_minus_WT": observed,
        "difference_ci95_low": float(ci[0]), "difference_ci95_high": float(ci[1]),
        "auc_in_sample_WT_referenced": _auc(labels, burden),
    }
    row.update(permutation)
    return row


def dependency_audit(features, config):
    names = list(config["features"])
    raw = features[names].apply(pd.to_numeric, errors="coerce")
    within = raw.copy()
    for name in names:
        within[name] = raw[name] - features.groupby("group")[name].transform("median")
    rows = []
    for mode, matrix in (("all_cells", raw), ("within_group_centered", within)):
        corr = matrix.corr(method="spearman")
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                rows.append({"mode": mode, "feature_a": a, "feature_b": b, "spearman_rho": corr.loc[a, b]})
    derived = [
        {
            "derived_feature": "rheobase_midpoint_pA_per_20ms_pF",
            "numerator_or_source": "rheobase_midpoint_pA",
            "denominator_or_linked_feature": "capacitance_20ms_pF",
            "mathematical_relation": "rheobase_midpoint_pA / capacitance_20ms_pF",
            "core_status": "excluded_from_primary_core",
        },
        {
            "derived_feature": "work_per_spike_high_current_fJ",
            "numerator_or_source": "external_work_positive_high_current_fJ",
            "denominator_or_linked_feature": "spike_count_high_current",
            "mathematical_relation": "sweep work divided by spike count before cell median",
            "core_status": "excluded_from_primary_core",
        },
        {
            "derived_feature": "spike_count_high_current",
            "numerator_or_source": "n_spikes",
            "denominator_or_linked_feature": "firing_rate_high_current_hz",
            "mathematical_relation": "firing rate equals spike count divided by the fixed stimulus duration",
            "core_status": "excluded_from_primary_core",
        },
    ]
    return pd.DataFrame(rows), pd.DataFrame(derived)


def leave_out_audit(features, labels, core_spec, feature_z, wt_masks, config):
    rows = []; variants = []
    for feature in core_spec:
        variants.append(("feature", feature, {k: v for k, v in core_spec.items() if k != feature}))
    for domain in dict.fromkeys(core_spec.values()):
        variants.append(("domain", domain, {k: v for k, v in core_spec.items() if v != domain}))
    for kind, removed, spec in variants:
        reference = fit_reference(features, labels, list(spec), config)
        scores, _ = score_specification(features, reference, spec, config)
        validation = observed_validation(
            "leave_one_{}_{}".format(kind, removed), features, scores, labels,
            feature_z, spec, wt_masks, config,
        )
        validation.update({"removal_type": kind, "removed": removed})
        rows.append(validation)
    return pd.DataFrame(rows)


def current_window_audit(data, labels, core_spec, wt_masks, config):
    rows = []
    for low, high in config["analysis"]["current_windows_pA"]:
        features = build_cell_features(config, data, low, high)
        reference = fit_reference(features, labels, list(core_spec), config)
        scores, domains = score_specification(features, reference, core_spec, config)
        scored, _ = add_coordinate_and_boundary(features, scores, domains, config)
        feature_z = permutation_feature_scores(features, list(core_spec), wt_masks, config)
        validation = observed_validation(
            "window_{}_{}".format(low, high), features, scores, labels,
            feature_z, core_spec, wt_masks, config,
        )
        validation.update({
            "current_min_pA": low, "current_max_pA": high,
            "WT_exit_WT_cells": int(scored.loc[scored.group == "WT", "WT_exit_consensus_marker"].sum()),
            "WT_exit_SCA3_cells": int(scored.loc[scored.group == "SCA3", "WT_exit_consensus_marker"].sum()),
        })
        rows.append(validation)
    return pd.DataFrame(rows)


def leave_one_cell_audit(features, labels, core_spec, full_reference, config):
    labels = np.asarray(labels); rows = []
    for idx in range(len(features)):
        keep = np.ones(len(features), dtype=bool); keep[idx] = False
        reduced = features.loc[keep].reset_index(drop=True); reduced_labels = labels[keep]
        reference = fit_reference(reduced, reduced_labels, list(core_spec), config, fallback=full_reference)
        scores, domains = score_specification(reduced, reference, core_spec, config)
        scored, _ = add_coordinate_and_boundary(reduced, scores, domains, config)
        burden = scored["fingerprint_burden_z"].to_numpy()
        rows.append({
            "excluded_group": labels[idx], "excluded_cell_id": features.iloc[idx]["cell_id"],
            "n_cells_remaining": int(keep.sum()),
            "median_difference_SCA3_minus_WT": float(
                np.median(burden[reduced_labels == "SCA3"]) - np.median(burden[reduced_labels == "WT"])
            ),
            "auc_in_sample_WT_referenced": _auc(reduced_labels, burden),
            "WT_exit_WT_cells": int(scored.loc[scored.group == "WT", "WT_exit_consensus_marker"].sum()),
            "WT_exit_SCA3_cells": int(scored.loc[scored.group == "SCA3", "WT_exit_consensus_marker"].sum()),
        })
    return pd.DataFrame(rows)


def threshold_audit(scored, domains, config):
    rows = []; domain_columns = ["domain_z__" + d for d in domains]
    for threshold in config["threshold_audit"]["robust_z_values"]:
        domain_hits = (scored[domain_columns] > float(threshold)).sum(axis=1)
        for min_domains in config["threshold_audit"]["minimum_domains_values"]:
            eligible = scored["domains_available"] >= int(min_domains)
            outside = eligible & (scored["fingerprint_burden_z"] > float(threshold))
            for consensus in config["threshold_audit"]["consensus_domains_values"]:
                calls = outside & (domain_hits >= int(consensus))
                wt_calls = int(calls[scored.group == "WT"].sum())
                sca_calls = int(calls[scored.group == "SCA3"].sum())
                rows.append({
                    "robust_z_threshold": threshold, "minimum_domains_available": min_domains,
                    "consensus_min_domains": consensus, "WT_calls": wt_calls, "SCA3_calls": sca_calls,
                    "WT_false_positive_fraction": wt_calls / int((scored.group == "WT").sum()),
                    "SCA3_call_fraction": sca_calls / int((scored.group == "SCA3").sum()),
                })
    return pd.DataFrame(rows)


def bootstrap_stability(features, labels, core_spec, full_reference, config):
    labels = np.asarray(labels); acfg = config["analysis"]
    iterations = int(acfg["bootstrap_iterations"])
    threshold = float(acfg["robust_z_threshold"])
    consensus = int(acfg["consensus_min_domains"])
    min_domains = int(acfg["minimum_domains_for_boundary"])
    rng = np.random.default_rng(int(acfg["random_seed"]) + 29)
    wt_indices = np.where(labels == "WT")[0]; rows = []
    for idx in range(len(features)):
        pool = wt_indices[wt_indices != idx] if labels[idx] == "WT" else wt_indices
        burdens = []; boundary_hits = 0; consensus_hits = 0
        for _ in range(iterations):
            sampled = rng.choice(pool, len(pool), replace=True)
            ref_features = features.iloc[sampled].reset_index(drop=True)
            ref_labels = np.full(len(sampled), "WT", dtype=object)
            reference = fit_reference(ref_features, ref_labels, list(core_spec), config, fallback=full_reference)
            scores, domains = score_specification(
                features.iloc[[idx]].reset_index(drop=True), reference, core_spec, config
            )
            row = scores.iloc[0]; burden = float(row["fingerprint_burden_z"])
            if not np.isfinite(burden):
                continue
            burdens.append(burden)
            outside = int(row["domains_available"]) >= min_domains and burden > threshold
            boundary_hits += int(outside)
            domain_hits = sum(
                float(row["domain_z__" + d]) > threshold
                for d in domains if np.isfinite(row["domain_z__" + d])
            )
            consensus_hits += int(outside and domain_hits >= consensus)
        values = np.asarray(burdens); valid = len(values)
        rows.append({
            "group": labels[idx], "cell_id": features.iloc[idx]["cell_id"],
            "bootstrap_reference_mode": "leave_one_WT_out" if labels[idx] == "WT" else "full_WT",
            "bootstrap_iterations_requested": iterations, "bootstrap_iterations_valid": valid,
            "bootstrap_burden_median": float(np.median(values)) if valid else np.nan,
            "bootstrap_burden_ci95_low": float(np.quantile(values, 0.025)) if valid else np.nan,
            "bootstrap_burden_ci95_high": float(np.quantile(values, 0.975)) if valid else np.nan,
            "bootstrap_p_outside_WT_robust_boundary": boundary_hits / valid if valid else np.nan,
            "bootstrap_p_WT_exit_consensus": consensus_hits / valid if valid else np.nan,
        })
    return pd.DataFrame(rows)


def cross_fitted_scores(features, labels, core_spec, full_reference, config, anchors):
    labels = np.asarray(labels); rows = []
    for idx in range(len(features)):
        ref_labels = labels.copy(); mode = "full_WT_reference"
        if labels[idx] == "WT":
            ref_labels[idx] = "EXCLUDED"; mode = "leave_one_WT_out"
        reference = fit_reference(features, ref_labels, list(core_spec), config, fallback=full_reference)
        scores, _ = score_specification(
            features.iloc[[idx]].reset_index(drop=True), reference, core_spec, config
        )
        row = scores.iloc[0].to_dict()
        row.update({
            "group": labels[idx], "cell_id": features.iloc[idx]["cell_id"], "crossfit_mode": mode,
            "crossfit_q_endpoint_unbounded": (
                float(row["fingerprint_burden_z"]) - anchors["wt_median_burden"]
            ) / anchors["q_denominator"],
        })
        rows.append(row)
    out = pd.DataFrame(rows)
    out["crossfit_q_endpoint_clipped_0_1"] = out["crossfit_q_endpoint_unbounded"].clip(0, 1)
    return out


def _plot_dependency_heatmap(dependencies, feature_specs, path):
    names = list(feature_specs); matrix = np.eye(len(names))
    subset = dependencies[dependencies["mode"] == "within_group_centered"]
    lookup = {(r.feature_a, r.feature_b): r.spearman_rho for r in subset.itertuples()}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i < j:
                matrix[i, j] = matrix[j, i] = lookup.get((a, b), lookup.get((b, a), np.nan))
    fig, ax = plt.subplots(figsize=(10, 8))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-1, vmax=1)
    short = [x.replace("_high_current", "").replace("_midpoint", "") for x in names]
    ax.set_xticks(range(len(names))); ax.set_yticks(range(len(names)))
    ax.set_xticklabels(short, rotation=60, ha="right", fontsize=8)
    ax.set_yticklabels(short, fontsize=8)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, "{:.2f}".format(matrix[i, j]), ha="center", va="center", fontsize=7)
    ax.set_title("Within-group Spearman dependence")
    fig.colorbar(image, ax=ax, label="Spearman rho")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def _plot_overview(scored, path, threshold):
    ordered = scored.sort_values("fingerprint_burden_z").reset_index(drop=True)
    colors = ordered.group.map({"WT": "#2a9d8f", "SCA3": "#e76f51"})
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(range(len(ordered)), ordered.fingerprint_burden_z, color=colors)
    ax.axhline(threshold, color="black", ls="--", lw=1, label="default robust boundary")
    ax.set_xticks(range(len(ordered))); ax.set_xticklabels(ordered.cell_id, rotation=90, fontsize=8)
    ax.set_ylabel("Dependency-reduced fingerprint burden")
    ax.legend(frameon=False); fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def _plot_sensitivity(variants, windows, path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    table = variants.sort_values("median_difference_SCA3_minus_WT")
    axes[0].barh(table.fingerprint, table.median_difference_SCA3_minus_WT, color="#457b9d")
    axes[0].set_xlabel("Median burden difference (SCA3 - WT)")
    axes[0].set_title("Fingerprint specifications")
    x = np.arange(len(windows))
    labels = ["{}–{}".format(int(a), int(b)) for a, b in zip(windows.current_min_pA, windows.current_max_pA)]
    axes[1].plot(x, windows.median_difference_SCA3_minus_WT, marker="o")
    axes[1].set_xticks(x); axes[1].set_xticklabels(labels, rotation=30)
    axes[1].set_xlabel("Current window (pA)"); axes[1].set_ylabel("Median burden difference")
    axes[1].set_title("Current-window sensitivity")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def _write_summary(path, validation, feature_comparisons, variant_validation, window_audit, leave_cell, scored):
    primary = variant_validation[variant_validation.fingerprint == "dependency_reduced_core"].iloc[0]
    absolute = feature_comparisons[feature_comparisons.feature == "rheobase_midpoint_pA"].iloc[0]
    density = feature_comparisons[
        feature_comparisons.feature == "rheobase_midpoint_pA_per_20ms_pF"
    ].iloc[0]
    work = feature_comparisons[feature_comparisons.feature == "work_per_spike_high_current_fJ"].iloc[0]
    lines = [
        "# NeuroThermo dependency-aware fingerprint v{}".format(__version__), "",
        "## Analysis unit", "",
        "- Recorded cells are independent observations: {} WT and {} SCA3.".format(
            validation["group_cell_counts"].get("WT", 0),
            validation["group_cell_counts"].get("SCA3", 0),
        ),
        "- Animal identity is removed at the input boundary and is not used.",
        "- Scope of inference: the recorded cells, not the animal population.", "",
        "## Dependency-reduced primary fingerprint", "",
        "- Median burden difference (SCA3 - WT): {:.6g}.".format(primary.median_difference_SCA3_minus_WT),
        "- Exact recomputed two-sided permutation p: {:.6g} across {} labelings.".format(
            primary.permutation_p_two_sided, int(primary.valid_labelings)
        ),
        "- In-sample WT-referenced AUC: {:.4f}.".format(primary.auc_in_sample_WT_referenced),
        "- Default WT-exit calls: {} WT and {} SCA3 cells.".format(
            int(scored.loc[scored.group == "WT", "WT_exit_consensus_marker"].sum()),
            int(scored.loc[scored.group == "SCA3", "WT_exit_consensus_marker"].sum()),
        ), "", "## Dependency audit", "",
        "- Absolute rheobase median difference: {:.6g} pA; exact p={:.6g}.".format(
            absolute.median_difference_SCA3_minus_WT_raw, absolute.permutation_p_two_sided
        ),
        "- Capacitance-normalized rheobase difference: {:.6g} pA/pF; exact p={:.6g}. It is a derived diagnostic, not independent core evidence.".format(
            density.median_difference_SCA3_minus_WT_raw, density.permutation_p_two_sided
        ),
        "- Work per spike difference: {:.6g} fJ; exact p={:.6g}. It remains a derived diagnostic linked to spike count.".format(
            work.median_difference_SCA3_minus_WT_raw, work.permutation_p_two_sided
        ), "", "## Robustness", "",
        "- Specification AUC range: {:.4f} to {:.4f}.".format(
            variant_validation.auc_in_sample_WT_referenced.min(),
            variant_validation.auc_in_sample_WT_referenced.max(),
        ),
        "- Current-window exact p range: {:.6g} to {:.6g}.".format(
            window_audit.permutation_p_two_sided.min(), window_audit.permutation_p_two_sided.max()
        ),
        "- Leave-one-cell AUC range: {:.4f} to {:.4f}.".format(
            leave_cell.auc_in_sample_WT_referenced.min(),
            leave_cell.auc_in_sample_WT_referenced.max(),
        ), "",
        "The endpoint q coordinate and WT-exit boundary are operational phenotype summaries. They are not disease probabilities, biological time, or evidence of irreversible degeneration.",
    ]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(config, upstream, output):
    validation = validate_inputs(config, upstream); data = _read_inputs(upstream)
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(validation["checks"]).to_csv(output / "input_validation.csv", index=False)
    features = build_cell_features(config, data); labels = features.group.to_numpy(dtype=object)
    wt_masks = exact_label_masks(labels, int(config["analysis"]["exact_max_labelings"]))
    all_spec_features = sorted({
        name for spec in config["fingerprints"].values() for name in spec["features"]
    })
    feature_z = permutation_feature_scores(features, all_spec_features, wt_masks, config)

    primary_names = [name for name, spec in config["fingerprints"].items() if spec.get("primary", False)]
    if len(primary_names) != 1:
        raise ValueError("Exactly one fingerprint specification must have primary: true")
    primary_name = primary_names[0]; primary_spec = config["fingerprints"][primary_name]["features"]
    primary_reference = fit_reference(features, labels, list(primary_spec), config)
    primary_scores, primary_domains = score_specification(
        features, primary_reference, primary_spec, config
    )
    scored, anchors = add_coordinate_and_boundary(features, primary_scores, primary_domains, config)

    variant_rows = []; specifications = []
    for name, item in config["fingerprints"].items():
        spec = item["features"]; reference = fit_reference(features, labels, list(spec), config)
        scores, _ = score_specification(features, reference, spec, config)
        variant_rows.append(observed_validation(
            name, features, scores, labels, feature_z, spec, wt_masks, config
        ))
        for feature, domain in spec.items():
            specifications.append({
                "fingerprint": name, "primary": bool(item.get("primary", False)),
                "feature": feature, "assigned_domain": domain,
                "feature_role": config["features"][feature].get("role", "core"),
            })
    variants = pd.DataFrame(variant_rows)
    leave_out = leave_out_audit(features, labels, primary_spec, feature_z, wt_masks, config)
    windows = current_window_audit(data, labels, primary_spec, wt_masks, config)
    leave_cell = leave_one_cell_audit(features, labels, primary_spec, primary_reference, config)
    thresholds = threshold_audit(scored, primary_domains, config)
    feature_comparisons = primary_feature_comparisons(features, labels, wt_masks, config)
    dependencies, derived = dependency_audit(features, config)
    bootstrap = bootstrap_stability(features, labels, primary_spec, primary_reference, config)
    crossfit = cross_fitted_scores(features, labels, primary_spec, primary_reference, config, anchors)

    long_rows = []
    for name in primary_spec:
        column = "feature_z__" + name
        for idx in range(len(features)):
            long_rows.append({
                "group": features.at[idx, "group"], "cell_id": features.at[idx, "cell_id"],
                "feature": name, "assigned_domain": primary_spec[name],
                "raw_value": features.at[idx, name],
                "oriented_robust_z": primary_scores.at[idx, column],
                "available": bool(np.isfinite(primary_scores.at[idx, column])),
            })
    long_scores = pd.DataFrame(long_rows)
    reference_rows = [{
        "feature": name, "assigned_domain": primary_spec[name],
        "role": config["features"][name].get("role", "core"), **ref,
    } for name, ref in primary_reference.items()]

    features.to_csv(output / "cell_features_dependency_audit.csv", index=False)
    scored.to_csv(output / "cell_dependency_reduced_scores.csv", index=False)
    long_scores.to_csv(output / "primary_feature_scores_long.csv", index=False)
    pd.DataFrame(reference_rows).to_csv(output / "primary_fingerprint_reference.csv", index=False)
    pd.DataFrame(specifications).to_csv(output / "fingerprint_specifications.csv", index=False)
    feature_comparisons.to_csv(output / "feature_comparisons.csv", index=False)
    variants.to_csv(output / "fingerprint_variant_validation.csv", index=False)
    leave_out.to_csv(output / "leave_one_feature_domain_audit.csv", index=False)
    windows.to_csv(output / "current_window_sensitivity.csv", index=False)
    leave_cell.to_csv(output / "leave_one_cell_audit.csv", index=False)
    thresholds.to_csv(output / "WT_exit_threshold_map.csv", index=False)
    dependencies.to_csv(output / "feature_dependency_correlations.csv", index=False)
    derived.to_csv(output / "derived_feature_audit.csv", index=False)
    bootstrap.to_csv(output / "bootstrap_cell_stability.csv", index=False)
    crossfit.to_csv(output / "cross_fitted_cell_scores.csv", index=False)

    _write_json(output / "frozen_dependency_reduced_scoring_rule.json", {
        "pipeline_version": __version__, "inference_unit": "cell", "animal_id_used": False,
        "primary_fingerprint": primary_name, "specification": primary_spec,
        "reference": primary_reference, "endpoint_anchors_and_boundaries": anchors,
        "analysis_configuration": config["analysis"],
    })
    _plot_dependency_heatmap(dependencies, config["features"], output / "feature_dependency_heatmap.png")
    _plot_overview(
        scored, output / "dependency_reduced_fingerprint_overview.png",
        float(config["analysis"]["robust_z_threshold"]),
    )
    _plot_sensitivity(variants, windows, output / "robustness_sensitivity.png")
    _write_summary(output / "RUN_SUMMARY.md", validation, feature_comparisons, variants, windows, leave_cell, scored)

    manifest = {
        "pipeline_name": "neurothermo-dependency-aware-fingerprint",
        "pipeline_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
        "inference_unit": "cell", "scope_of_inference": "recorded_cells", "animal_id_used": False,
        "upstream_dir": str(upstream), "output_dir": str(output), "n_cells": validation["n_cells"],
        "group_cell_counts": validation["group_cell_counts"], "primary_fingerprint": primary_name,
        "exact_cell_labelings": int(len(wt_masks)),
        "input_sha256": {item: _sha256(upstream / item) for item in REQUIRED_INPUTS},
        "generated_files": sorted(p.name for p in output.iterdir() if p.is_file()),
        "configuration": {
            k: config[k] for k in ("analysis", "threshold_audit", "features", "fingerprints")
        },
    }
    _write_json(output / "analysis_manifest.json", manifest)
    return output
