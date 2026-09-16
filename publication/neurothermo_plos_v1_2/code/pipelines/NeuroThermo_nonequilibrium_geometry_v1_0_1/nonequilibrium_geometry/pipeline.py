from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from pathlib import Path
import json
import os
import traceback
import hashlib

import numpy as np
import pandas as pd

from .data import REQUIRED, load_frozen, resolve_frozen_dir, scenario_arrays, select_scenarios, sha256, write_input_manifest
from .density import build_masses, common_grid, quantile_state_edges
from .geometry import adaptive_indices, cumulative_length, geometry_for_path
from .figures import make_figures
from .markov import detailed_balance_metrics, exact_hatano_sasa, simulate_hatano_sasa, simulate_path_probability_ift, transition_matrix
from .model import stationary_samples
from .preflight import run_preflight
from .verdict import formalism_verdict, validate_physical_mapping
from .aggregation import annotate_animal_pairs, balanced_mean, distribution_summary, load_animal_mapping


def _json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(type(value).__name__)


def _task_seed(base, scenario_id, p_index=0):
    return int((int(base) + 100003 * int(scenario_id) + 1009 * int(p_index)) % (2**31 - 1))


def _representative_indices(p_grid, requested):
    return sorted({int(np.argmin(np.abs(np.asarray(p_grid) - float(value)))) for value in requested})


def analyse_task(row_dict, base_seed, cfg, checkpoint_dir, fingerprint):
    row = pd.Series(row_dict)
    scenario_id = int(row.scenario_id)
    checkpoint = Path(checkpoint_dir) / f"scenario_{scenario_id:05d}_seed_{int(base_seed)}.json"
    cache_required = bool(cfg.get("markov", {}).get("save_cache", False))
    cache_dir = Path(checkpoint_dir).parent / "markov_cache"
    cache_path = cache_dir / f"scenario_{scenario_id:05d}_seed_{int(base_seed)}.npz"
    if cfg["output"].get("resume", True) and checkpoint.is_file():
        cached = json.loads(checkpoint.read_text(encoding="utf-8"))
        cache_ready = (not cache_required) or cache_path.is_file()
        if cached.get("_fingerprint") == fingerprint and cached.get("_status") == "OK" and cache_ready:
            return cached
    p_grid = np.linspace(0.0, 1.0, int(cfg["path"]["n_p"]))
    theta, current = scenario_arrays(row, p_grid)
    samples_by_p, start = [], None
    for index in range(len(p_grid)):
        samples, start, finite = stationary_samples(
            _task_seed(base_seed, scenario_id, index), theta[index], current[index], cfg, start=start
        )
        if not finite or len(samples) < int(cfg["stationary"]["minimum_saved_samples"]):
            raise RuntimeError(f"Non-finite/short stationary run: scenario={scenario_id}, seed={base_seed}, p={p_grid[index]:g}")
        samples_by_p.append(samples)
    edges, centers = common_grid(samples_by_p, cfg)
    masses, retained = build_masses(samples_by_p, edges, cfg)
    geometry_rows, kl_rows, dphi, safe_mass = geometry_for_path(samples_by_p, masses, edges, centers, theta, current, p_grid, cfg)

    state_edges = quantile_state_edges(samples_by_p, cfg["markov"]["state_shape"])
    lag_samples = max(1, int(round(float(cfg["markov"]["lag_ms"]) / float(cfg["stationary"]["sample_stride_ms"]))))
    transitions, stationary, cycle_rows = [], [], []
    for index, samples in enumerate(samples_by_p):
        matrix, invariant, occupancy, _ = transition_matrix(samples, state_edges, lag_samples, cfg["markov"]["pseudocount"])
        metrics, affinities = detailed_balance_metrics(matrix, invariant, cfg, empirical_occupancy=occupancy)
        transitions.append(matrix)
        stationary.append(invariant)
        geometry_rows[index].update(metrics)
        geometry_rows[index]["retained_sample_fraction"] = float(retained[index])
        for affinity in affinities:
            cycle_rows.append({
                "scenario_id": scenario_id,
                "biological_pair_key": str(row.biological_pair_key),
                "seed": int(base_seed),
                "p": float(p_grid[index]),
                "affinity": float(affinity),
                "abs_affinity": abs(float(affinity)),
            })

    if cache_required:
        cache_dir.mkdir(exist_ok=True)
        temporary_cache = cache_path.with_suffix(".tmp")
        with temporary_cache.open("wb") as handle:
            np.savez_compressed(
                handle,
                scientific_fingerprint=np.asarray(fingerprint),
                p_grid=p_grid,
                transitions=np.asarray(transitions, dtype=np.float32),
                stationary=np.asarray(stationary, dtype=np.float64),
            )
        os.replace(temporary_cache, cache_path)

    path_fi = np.asarray([item["path_fi"] for item in geometry_rows])
    path_fi_xy = np.asarray([item["path_fi_xy"] for item in geometry_rows])
    friction = np.asarray([item["friction_metric"] for item in geometry_rows])
    fi_cumulative, fi_normalized, fi_total = cumulative_length(p_grid, path_fi)
    fr_cumulative, fr_normalized, fr_total = cumulative_length(p_grid, np.nan_to_num(friction, nan=0.0))
    protocol_points = int(cfg["protocol"]["n_points"])
    adaptive_fi = adaptive_indices(p_grid, path_fi, protocol_points)
    adaptive_fi_xy = adaptive_indices(p_grid, path_fi_xy, protocol_points)
    adaptive_friction = adaptive_indices(p_grid, np.nan_to_num(friction, nan=0.0), protocol_points)
    linear_indices = adaptive_indices(p_grid, np.ones_like(p_grid), protocol_points)
    for index, item in enumerate(geometry_rows):
        item.update({
            "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
            "analysis_scenario_weight": float(row.analysis_scenario_weight),
            "analysis_within_pair_weight": float(row.analysis_within_pair_weight),
            "thermodynamic_length_path_fi_cumulative": float(fi_cumulative[index]),
            "thermodynamic_length_path_fi_fraction": float(fi_normalized[index]),
            "thermodynamic_length_friction_cumulative": float(fr_cumulative[index]),
            "thermodynamic_length_friction_fraction": float(fr_normalized[index]),
        })
    local_kl = [{
        "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
        "left_p": float(p_grid[index]), "right_p": float(p_grid[index + 1]),
        "kl_forward": kl, "fisher_quadratic": prediction, "relative_error": relative,
    } for index, kl, prediction, relative in kl_rows]
    protocol = []
    for kind, indices in (("path_fi_xyz", adaptive_fi), ("path_fi_xy", adaptive_fi_xy), ("friction", adaptive_friction)):
        for index, grid_index in enumerate(indices):
            protocol.append({
                "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
                "metric": kind, "protocol_index": index, "normalized_time": index / max(len(indices) - 1, 1),
                "path_index": int(grid_index), "p": float(p_grid[grid_index]),
            })

    fluctuation = []
    schedules = {
        "linear": linear_indices,
        "path_fi_xyz": adaptive_fi,
        "path_fi_xy": adaptive_fi_xy,
        "friction": adaptive_friction,
    }
    for schedule_name, indices in schedules.items():
        schedule_matrices = [transitions[index] for index in indices]
        schedule_stationary = [stationary[index] for index in indices]
        for direction, matrices, distributions in (
            ("WT_to_SCA3", schedule_matrices, schedule_stationary),
            ("SCA3_to_WT", schedule_matrices[::-1], schedule_stationary[::-1]),
        ):
            offset = sum(ord(character) for character in schedule_name + direction)
            mc = simulate_hatano_sasa(matrices, distributions, cfg["fluctuation"]["trajectories"], _task_seed(base_seed + 700001 + offset, scenario_id))
            ift = simulate_path_probability_ift(matrices, distributions, cfg["fluctuation"]["trajectories"], _task_seed(base_seed + 900001 + offset, scenario_id))
            fluctuation.append({
                "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
                "schedule": schedule_name, "direction": direction, "n_protocol_steps": len(indices),
                "n_unique_path_positions": int(len(np.unique(indices))),
                "exact_mean_exp_minus_Y": exact_hatano_sasa(matrices, distributions),
                **mc, **ift,
            })

    potential = []
    for index in _representative_indices(p_grid, cfg["output"]["potential_p"]):
        mass = safe_mass[index]
        marginals = (mass.sum(axis=(1, 2)), mass.sum(axis=(0, 2)), mass.sum(axis=(0, 1)))
        for axis, axis_mass, axis_centers in zip(("x", "y", "z"), marginals, centers):
            phi = -np.log(axis_mass)
            phi -= np.nanmin(phi)
            for coordinate, probability, potential_value in zip(axis_centers, axis_mass, phi):
                potential.append({
                    "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
                    "p": float(p_grid[index]), "axis": axis, "coordinate": float(coordinate),
                    "marginal_probability": float(probability), "relative_phi": float(potential_value),
                })
    sensitivity = []
    for variant in cfg.get("sensitivity", {}).get("density_variants", []):
        variant_cfg = deepcopy(cfg)
        variant_cfg["density"].update({key: value for key, value in variant.items() if key != "name"})
        variant_edges, variant_centers = common_grid(samples_by_p, variant_cfg)
        variant_masses, _ = build_masses(samples_by_p, variant_edges, variant_cfg)
        variant_rows, _, _, _ = geometry_for_path(
            samples_by_p, variant_masses, variant_edges, variant_centers, theta, current, p_grid, variant_cfg
        )
        for item in variant_rows:
            sensitivity.append({
                "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
                "variant_type": "density", "variant": str(variant["name"]), "p": item["p"],
                "path_fi_xyz": item["path_fi_xyz"], "path_fi_xy": item["path_fi_xy"],
                "circulation_fraction": item["circulation_fraction"],
                "stationary_current_divergence_relative": item["stationary_current_divergence_relative"],
            })
    for lag_ms in cfg.get("sensitivity", {}).get("markov_lag_ms", []):
        variant_lag = max(1, int(round(float(lag_ms) / float(cfg["stationary"]["sample_stride_ms"]))))
        for index, samples in enumerate(samples_by_p):
            matrix, invariant, occupancy, _ = transition_matrix(samples, state_edges, variant_lag, cfg["markov"]["pseudocount"])
            metrics, _ = detailed_balance_metrics(matrix, invariant, cfg, empirical_occupancy=occupancy)
            sensitivity.append({
                "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
                "variant_type": "markov_lag", "variant": f"lag_{float(lag_ms):g}_ms", "p": float(p_grid[index]),
                "markov_db_violation": metrics["markov_db_violation"],
                "markov_entropy_per_lag": metrics["markov_entropy_per_lag"],
                "markov_empirical_pi_l1": metrics["markov_empirical_pi_l1"],
            })
    summary = [{
        "scenario_id": scenario_id, "biological_pair_key": str(row.biological_pair_key), "seed": int(base_seed),
        "path_fi_length": fi_total, "friction_length": fr_total,
    }]
    payload = {"_status": "OK", "_fingerprint": fingerprint,
               "geometry": geometry_rows, "local_kl": local_kl, "cycles": cycle_rows, "protocol": protocol,
               "fluctuation": fluctuation, "potential": potential, "sensitivity": sensitivity, "summary": summary}
    temporary = checkpoint.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, default=_json_default), encoding="utf-8")
    os.replace(temporary, checkpoint)
    return payload


def _flatten(results, key):
    return pd.DataFrame([row for result in results for row in result[key]])


def _fingerprint(cfg, frozen):
    scientific_cfg = {key: value for key, value in cfg.items() if key not in ("output", "parallel")}
    payload = {
        "pipeline_version": "1.0.1",
        "config": scientific_cfg,
        "frozen_sha256": {name: sha256(frozen / name) for name in REQUIRED},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _pair_balance(geometry):
    excluded = {"scenario_id", "seed", "p", "analysis_scenario_weight", "analysis_within_pair_weight"}
    numeric = [column for column in geometry.select_dtypes(include=[np.number]).columns if column not in excluded]
    rows = []
    for (pair, seed, p), group in geometry.groupby(["biological_pair_key", "seed", "p"], sort=True):
        weights = group["analysis_within_pair_weight"].to_numpy(float)
        weights /= weights.sum()
        row = {"biological_pair_key": pair, "seed": int(seed), "p": float(p), "n_support_scenarios": len(group)}
        for column in numeric:
            values = group[column].to_numpy(float)
            valid = np.isfinite(values)
            row[column] = float(np.average(values[valid], weights=weights[valid])) if valid.any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _weighted_ensemble(pair_geometry):
    numeric = [column for column in pair_geometry.select_dtypes(include=[np.number]).columns if column not in ("seed", "p")]
    rows = []
    for p, group in pair_geometry.groupby("p", sort=True):
        row = {"p": p, "n_pair_seed": len(group), "n_pairs": group["biological_pair_key"].nunique()}
        for column in numeric:
            values = group[column].dropna()
            row[column + "_median"] = float(values.median()) if len(values) else np.nan
            row[column + "_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[column + "_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _balance_protocol_or_fluctuation(frame, scenarios, grouping):
    weights = scenarios[["scenario_id", "biological_pair_key", "analysis_within_pair_weight"]].drop_duplicates()
    merged = frame.merge(weights, on=["scenario_id", "biological_pair_key"], how="left", validate="many_to_one")
    if merged["analysis_within_pair_weight"].isna().any():
        raise RuntimeError("Missing scenario weights during balanced aggregation")
    merged = merged.drop(columns=["scenario_id"])
    return balanced_mean(
        merged,
        ["biological_pair_key", "seed"] + list(grouping),
        weight_column="analysis_within_pair_weight",
        count_name="n_support_scenarios",
    )


def _protocol_performance(fluctuation):
    rows = []
    for (schedule, direction), group in fluctuation.groupby(["schedule", "direction"], sort=True):
        rows.append({
            "schedule": schedule,
            "direction": direction,
            "n_animal_pair_seed": len(group),
            "mean_exp_minus_Y_median": float(group["mean_exp_minus_Y"].median()),
            "mean_exp_minus_Y_q25": float(group["mean_exp_minus_Y"].quantile(0.25)),
            "mean_exp_minus_Y_q75": float(group["mean_exp_minus_Y"].quantile(0.75)),
            "absolute_error_from_one_median": float((group["mean_exp_minus_Y"] - 1.0).abs().median()),
            "ess_fraction_median": float(group["ess_fraction"].median()),
            "median_Y_median": float(group["median_Y"].median()),
            "mean_exp_minus_sigma_median": float(group["mean_exp_minus_sigma"].median()),
            "median_sigma_median": float(group["median_sigma"].median()),
            "unique_positions_min": int(group["n_unique_path_positions"].min()),
            "unique_positions_max": int(group["n_unique_path_positions"].max()),
        })
    return pd.DataFrame(rows)


def _protocol_verdict(performance, n_points):
    comparisons = []
    supported = True
    for direction in sorted(performance["direction"].unique()):
        baseline = performance.loc[
            performance["schedule"].eq("linear") & performance["direction"].eq(direction)
        ]
        if len(baseline) != 1:
            raise RuntimeError(f"Missing unique linear baseline for {direction}")
        baseline = baseline.iloc[0]
        for schedule in ("path_fi_xyz", "path_fi_xy", "friction"):
            candidate = performance.loc[
                performance["schedule"].eq(schedule) & performance["direction"].eq(direction)
            ]
            if len(candidate) != 1:
                raise RuntimeError(f"Missing unique {schedule} result for {direction}")
            candidate = candidate.iloc[0]
            error_improved = bool(
                candidate["absolute_error_from_one_median"] < baseline["absolute_error_from_one_median"]
            )
            ess_improved = bool(candidate["ess_fraction_median"] > baseline["ess_fraction_median"])
            supported &= error_improved and ess_improved
            comparisons.append({
                "schedule": schedule,
                "direction": direction,
                "absolute_error_improved": error_improved,
                "ess_improved": ess_improved,
            })
    return {
        "adaptive_sampling_improvement": "SUPPORTED" if supported else "NOT_SUPPORTED",
        "comparison_basis": "equal_unique_path_positions",
        "required_unique_positions": int(n_points),
        "comparisons": comparisons,
    }


def _numerical_validation(geometry, local_kl, fluctuation, cfg):
    divergence = geometry["stationary_current_divergence_relative"]
    chain = geometry["path_fi_chain_rule_violation"]
    return {
        "centered_path_fisher": True,
        "friction_covariance_amplitude_source": "centered_density_path_fisher",
        "sample_to_metric_force_variance_ratio_median": float(
            geometry["force_variance_consistency_ratio"].replace([np.inf, -np.inf], np.nan).median()
        ),
        "fisher_chain_rule_violation_fraction": float((chain > 1e-10).mean()),
        "local_kl_relative_error_median": float(local_kl["relative_error"].median()),
        "local_kl_relative_error_le_0_25_fraction": float((local_kl["relative_error"] <= 0.25).mean()),
        "continuous_current_adequate_fraction": float(
            (divergence <= float(cfg["gates"]["max_current_divergence_relative"])).mean()
        ),
        "all_protocols_have_required_unique_positions": bool(
            (fluctuation["n_unique_path_positions"] == int(cfg["protocol"]["n_points"])).all()
        ),
    }


def run(cfg, frozen_dir):
    validate_physical_mapping(cfg)
    frozen = resolve_frozen_dir(frozen_dir)
    scenarios_all, pair_stage, _, _ = load_frozen(frozen)
    scenarios = select_scenarios(scenarios_all, pair_stage, cfg)
    animal_mapping = load_animal_mapping(cfg)
    mapping_required = bool(cfg.get("animal_mapping", {}).get("required", False))
    if animal_mapping is not None:
        annotate_animal_pairs(
            scenarios[["biological_pair_key"]].drop_duplicates(), animal_mapping, required=mapping_required
        )
    output = Path(cfg["output"]["dir"]).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    checkpoint = output / "checkpoints"
    checkpoint.mkdir(exist_ok=True)
    write_input_manifest(frozen, output)
    fingerprint = _fingerprint(cfg, frozen)
    fixed, continuation, endpoints, endpoint_membership = run_preflight(scenarios, cfg)
    fixed.to_csv(output / "preflight_fixed_points.csv", index=False)
    continuation.to_csv(output / "preflight_current_continuation.csv", index=False)
    endpoints.to_csv(output / "preflight_endpoint_summary.csv", index=False)
    endpoint_membership.to_csv(output / "preflight_endpoint_membership.csv", index=False)
    if not endpoints["finite"].all() or not continuation["finite"].all():
        raise RuntimeError("Deterministic preflight failed: at least one trajectory was non-finite.")

    tasks = [(row.to_dict(), int(seed), cfg, str(checkpoint), fingerprint) for _, row in scenarios.iterrows() for seed in cfg["multiseed"]["seeds"]]
    results, failures = [], []
    workers = int(cfg["parallel"]["workers"])
    if workers == 1:
        for task in tasks:
            try:
                results.append(analyse_task(*task))
            except Exception as error:
                failures.append({"scenario_id": task[0]["scenario_id"], "seed": task[1], "error": str(error), "traceback": traceback.format_exc()})
                if not cfg["output"].get("continue_on_error", False):
                    raise
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(analyse_task, *task): task for task in tasks}
            for future in as_completed(futures):
                task = futures[future]
                try:
                    results.append(future.result())
                except Exception as error:
                    failures.append({"scenario_id": task[0]["scenario_id"], "seed": task[1], "error": str(error), "traceback": traceback.format_exc()})
                    if not cfg["output"].get("continue_on_error", False):
                        raise
    if not results:
        raise RuntimeError("No scenario-seed task completed.")
    frames = {
        "stationary_geometry.csv": _flatten(results, "geometry"),
        "local_kl_fisher_check.csv": _flatten(results, "local_kl"),
        "markov_cycle_affinities.csv": _flatten(results, "cycles"),
        "adaptive_protocols.csv": _flatten(results, "protocol"),
        "fluctuation_relations.csv": _flatten(results, "fluctuation"),
        "potential_marginals.csv": _flatten(results, "potential"),
        "estimator_sensitivity.csv": _flatten(results, "sensitivity"),
        "thermodynamic_length_summary.csv": _flatten(results, "summary"),
    }
    for name, frame in frames.items():
        frame.to_csv(output / name, index=False)
    pair_geometry = _pair_balance(frames["stationary_geometry.csv"])
    pair_geometry.to_csv(output / "pair_balanced_geometry.csv", index=False)
    cell_ensemble = _weighted_ensemble(pair_geometry)
    cell_ensemble.to_csv(output / "ensemble_geometry.csv", index=False)

    pair_protocol = _balance_protocol_or_fluctuation(
        frames["adaptive_protocols.csv"], scenarios, ["metric", "protocol_index", "normalized_time"]
    )
    pair_protocol.to_csv(output / "pair_balanced_protocols.csv", index=False)
    pair_fluctuation = _balance_protocol_or_fluctuation(
        frames["fluctuation_relations.csv"], scenarios, ["schedule", "direction"]
    )
    pair_fluctuation.to_csv(output / "pair_balanced_fluctuation_relations.csv", index=False)

    figure_ensemble = cell_ensemble
    figure_protocol = pair_protocol
    figure_fluctuation = pair_fluctuation
    animal_pair_geometry = None
    animal_pair_count = 0
    wt_animal_count = 0
    sca3_animal_count = 0
    if animal_mapping is not None:
        annotated_geometry = annotate_animal_pairs(pair_geometry, animal_mapping, required=mapping_required)
        animal_pair_geometry = balanced_mean(
            annotated_geometry,
            ["animal_pair_key", "wt_animal_id", "sca3_animal_id", "seed", "p"],
            count_name="n_cell_pairs",
        )
        animal_pair_geometry.to_csv(output / "animal_pair_balanced_geometry.csv", index=False)
        animal_ensemble = distribution_summary(animal_pair_geometry, ["p"])
        animal_ensemble.insert(1, "n_animal_pair_seed", animal_pair_geometry.groupby("p").size().to_numpy())
        animal_ensemble.insert(2, "n_animal_pairs", animal_pair_geometry["animal_pair_key"].nunique())
        animal_ensemble.to_csv(output / "animal_balanced_geometry.csv", index=False)

        annotated_protocol = annotate_animal_pairs(pair_protocol, animal_mapping, required=mapping_required)
        animal_pair_protocol = balanced_mean(
            annotated_protocol,
            ["animal_pair_key", "wt_animal_id", "sca3_animal_id", "seed", "metric", "protocol_index", "normalized_time"],
            count_name="n_cell_pairs",
        )
        animal_pair_protocol.to_csv(output / "animal_pair_balanced_protocols.csv", index=False)

        annotated_fluctuation = annotate_animal_pairs(pair_fluctuation, animal_mapping, required=mapping_required)
        animal_pair_fluctuation = balanced_mean(
            annotated_fluctuation,
            ["animal_pair_key", "wt_animal_id", "sca3_animal_id", "seed", "schedule", "direction"],
            count_name="n_cell_pairs",
        )
        animal_pair_fluctuation.to_csv(output / "animal_pair_balanced_fluctuation_relations.csv", index=False)
        performance = _protocol_performance(animal_pair_fluctuation)
        performance.to_csv(output / "protocol_performance_summary.csv", index=False)

        selected_cells = sorted(
            set(annotated_geometry["wt_cell_id"]).union(set(annotated_geometry["sca3_cell_id"]))
        )
        animal_mapping.loc[animal_mapping["cell_id"].isin(selected_cells)].sort_values(
            ["genotype", "animal_id", "cell_id"]
        ).to_csv(output / "animal_mapping_used.csv", index=False)
        animal_pair_count = int(animal_pair_geometry["animal_pair_key"].nunique())
        wt_animal_count = int(animal_pair_geometry["wt_animal_id"].nunique())
        sca3_animal_count = int(animal_pair_geometry["sca3_animal_id"].nunique())
        figure_ensemble = animal_ensemble
        figure_protocol = animal_pair_protocol
        figure_fluctuation = animal_pair_fluctuation
    else:
        performance = _protocol_performance(pair_fluctuation)
        performance.to_csv(output / "protocol_performance_summary.csv", index=False)

    protocol_verdict = _protocol_verdict(performance, cfg["protocol"]["n_points"])
    (output / "PROTOCOL_VERDICT.json").write_text(json.dumps(protocol_verdict, indent=2) + "\n", encoding="utf-8")
    numerical_validation = _numerical_validation(
        frames["stationary_geometry.csv"],
        frames["local_kl_fisher_check.csv"],
        frames["fluctuation_relations.csv"],
        cfg,
    )
    (output / "NUMERICAL_VALIDATION.json").write_text(
        json.dumps(numerical_validation, indent=2) + "\n", encoding="utf-8"
    )

    endpoint_status = endpoint_membership.merge(
        endpoints[["endpoint_group_id", "oscillatory"]], on="endpoint_group_id", how="left", validate="many_to_one"
    )
    oscillatory_ids = set(
        endpoint_status.groupby("scenario_id")["oscillatory"].all().loc[lambda x: x].index.astype(int)
    )
    oscillatory_geometry = frames["stationary_geometry.csv"].loc[
        frames["stationary_geometry.csv"]["scenario_id"].astype(int).isin(oscillatory_ids)
    ].copy()
    oscillatory_pair_geometry = _pair_balance(oscillatory_geometry) if len(oscillatory_geometry) else pd.DataFrame()
    oscillatory_pair_geometry.to_csv(output / "oscillatory_endpoint_pair_balanced_geometry.csv", index=False)

    if cfg["output"].get("make_figures", True):
        make_figures(output, figure_ensemble, figure_protocol, figure_fluctuation)
    pd.DataFrame(failures).to_csv(output / "failures.csv", index=False)
    verdict = formalism_verdict(
        frames["stationary_geometry.csv"],
        cfg,
        animal_geometry=animal_pair_geometry,
        oscillatory_geometry=oscillatory_geometry,
    )
    (output / "FORMALISM_VERDICT.json").write_text(json.dumps(verdict, indent=2) + "\n", encoding="utf-8")
    validation = {
        "status": "PASS" if not failures else "PARTIAL",
        "scientific_fingerprint": fingerprint,
        "frozen_dir": str(frozen), "selected_scenarios": len(scenarios), "scenario_seed_tasks": len(tasks),
        "completed_tasks": len(results), "failed_tasks": len(failures),
        "preflight_all_finite": bool(endpoints["finite"].all() and continuation["finite"].all()),
        "preflight_oscillatory_endpoint_groups": int(endpoints["oscillatory"].sum()),
        "preflight_total_endpoint_groups": int(len(endpoints)),
        "animal_mapping_applied": animal_mapping is not None,
        "animal_pairs": animal_pair_count,
        "wt_animals": wt_animal_count,
        "sca3_animals": sca3_animal_count,
        "protocol": protocol_verdict,
        "numerical_validation": numerical_validation,
        "formalism": verdict,
    }
    (output / "RUN_SUMMARY.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    return output, validation
