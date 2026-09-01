from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import os
import time
import traceback

import numpy as np
import pandas as pd

from .aggregation import animal_pair_balanced, cell_balanced, leave_one_animal_out
from .data import REQUIRED, load_frozen, resolve_frozen_dir, scenario_arrays, select_scenarios, sha256, write_input_manifest
from .density import build_masses, full_coverage_grid
from .markers import MARKER_VARIANTS, seed_ensemble_markers, weighted_quantile
from .model import stationary_samples_nested
from .verdict import build_verdict


VIEWS = ("xyz", "xy", "z")


def _task_seed(base, scenario_id, p_index):
    return int((int(base) + 100003 * int(scenario_id) + 1009 * int(p_index)) % (2**31 - 1))


def _dt_label(value):
    return ("%.8g" % float(value)).replace(".", "p")


def _atomic_json(path, payload):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    os.replace(str(temporary), str(path))


def _fingerprint(cfg, frozen):
    scientific = deepcopy(cfg)
    scientific.pop("_config_path", None)
    scientific.pop("parallel", None)
    scientific.get("output", {}).pop("dir", None)
    scientific.get("output", {}).pop("resume", None)
    mapping_path = Path(scientific["animal_mapping"].pop("path"))
    scientific["animal_mapping"]["sha256"] = sha256(mapping_path)
    payload = {
        "config": scientific,
        "inputs": {name: sha256(Path(frozen) / name) for name in REQUIRED},
        "version": "1.0.1",
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _kl(left, right):
    left = np.asarray(left, float)
    right = np.asarray(right, float)
    return float(np.sum(left * (np.log(left) - np.log(right))))


def _view_masses(masses, view):
    if view == "xyz":
        return masses
    if view == "xy":
        return masses.sum(axis=3)
    if view == "z":
        return masses.sum(axis=(1, 2))
    raise ValueError("Unknown view: " + str(view))


def _endpoint_curves(masses):
    result = {}
    for view in VIEWS:
        current = _view_masses(masses, view)
        wt, sca = current[0], current[-1]
        kl_wt = np.asarray([_kl(item, wt) for item in current])
        kl_sca = np.asarray([_kl(item, sca) for item in current])
        result[view] = {
            "kl_wt": kl_wt.tolist(),
            "kl_sca3": kl_sca.tolist(),
            "delta_kl": (kl_wt - kl_sca).tolist(),
        }
    return result


def _simulate_path(row, seed, dt_ms, cfg):
    p_grid = np.linspace(0.0, 1.0, int(cfg["path"]["n_p"]))
    theta, current = scenario_arrays(row, p_grid)
    samples = []
    start = None
    for index in range(len(p_grid)):
        values, final, ok = stationary_samples_nested(
            _task_seed(seed, int(row.scenario_id), index), theta[index], current[index], cfg, dt_ms, start=start
        )
        if not ok or len(values) < 10:
            raise RuntimeError("Non-finite/short run: scenario=%s seed=%s dt=%s p=%s" %
                               (row.scenario_id, seed, dt_ms, p_grid[index]))
        samples.append(values)
        start = final
    return p_grid, samples


def _checkpoint_path(output, scenario_id, dt_ms, seed):
    return Path(output) / "checkpoints" / ("scenario_%05d_dt_%s_seed_%d.json" %
                                            (int(scenario_id), _dt_label(dt_ms), int(seed)))


def _extent_checkpoint_path(output, scenario_id, dt_ms, seed):
    return Path(output) / "extent_checkpoints" / ("scenario_%05d_dt_%s_seed_%d.json" %
                                                   (int(scenario_id), _dt_label(dt_ms), int(seed)))


def _grid_path(output, scenario_id):
    return Path(output) / "reference_grids" / ("scenario_%05d.npz" % int(scenario_id))


def _valid_checkpoint(path, fingerprint):
    if not Path(path).is_file():
        return None
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("_status") == "OK" and payload.get("_fingerprint") == fingerprint:
        return payload
    return None


def _valid_extent(path, fingerprint):
    payload = _valid_checkpoint(path, fingerprint)
    if payload is None:
        return None
    required = ("axis_min", "axis_max", "sample_count")
    if not all(key in payload for key in required):
        return None
    return payload


def _save_grid(path, edges, fingerprint):
    path = Path(path)
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, x=edges[0], y=edges[1], z=edges[2], fingerprint=np.asarray(fingerprint))
    os.replace(str(temporary), str(path))


def _load_grid(path, fingerprint):
    with np.load(path, allow_pickle=False) as data:
        if str(data["fingerprint"].item()) != fingerprint:
            raise RuntimeError("Reference-grid fingerprint mismatch: " + str(path))
        return [data["x"], data["y"], data["z"]]


def _payload(row, seed, dt_ms, p_grid, masses, retained, fingerprint):
    return {
        "_status": "OK",
        "_fingerprint": fingerprint,
        "scenario_id": int(row.scenario_id),
        "biological_pair_key": str(row.biological_pair_key),
        "seed": int(seed),
        "dt_ms": float(dt_ms),
        "p": p_grid.tolist(),
        "retained_min": float(np.min(retained)),
        "retained_median": float(np.median(retained)),
        "views": _endpoint_curves(masses),
    }


def _extent_payload(row, seed, dt_ms, samples, fingerprint):
    axis_min = np.min(np.vstack([np.min(item, axis=0) for item in samples]), axis=0)
    axis_max = np.max(np.vstack([np.max(item, axis=0) for item in samples]), axis=0)
    return {
        "_status": "OK",
        "_fingerprint": fingerprint,
        "scenario_id": int(row.scenario_id),
        "biological_pair_key": str(row.biological_pair_key),
        "seed": int(seed),
        "dt_ms": float(dt_ms),
        "axis_min": axis_min.tolist(),
        "axis_max": axis_max.tolist(),
        "sample_count": int(sum(len(item) for item in samples)),
    }


def pilot_extent_task(row_dict, seed, dt_ms, cfg, output, fingerprint):
    row = pd.Series(row_dict)
    checkpoint = _extent_checkpoint_path(output, row.scenario_id, dt_ms, seed)
    cached = _valid_extent(checkpoint, fingerprint)
    if cached is not None:
        return cached
    p_grid, samples = _simulate_path(row, seed, dt_ms, cfg)
    if len(p_grid) != len(samples):
        raise RuntimeError("Pilot path/sample length mismatch")
    payload = _extent_payload(row, seed, dt_ms, samples, fingerprint)
    _atomic_json(checkpoint, payload)
    return payload


def build_scenario_grid(row, cfg, output, fingerprint):
    extents = []
    for dt_ms in cfg["convergence"]["dt_ms"]:
        for seed in cfg["convergence"]["seeds"]:
            path = _extent_checkpoint_path(output, row.scenario_id, dt_ms, seed)
            payload = _valid_extent(path, fingerprint)
            if payload is None:
                raise RuntimeError("Missing valid pilot extent: " + str(path))
            extents.append((payload["axis_min"], payload["axis_max"]))
    edges, _ = full_coverage_grid(extents, cfg)
    _save_grid(_grid_path(output, row.scenario_id), edges, fingerprint)
    return {
        "scenario_id": int(row.scenario_id),
        "n_extent_tasks": len(extents),
        "grid_x_min": float(edges[0][0]), "grid_x_max": float(edges[0][-1]),
        "grid_y_min": float(edges[1][0]), "grid_y_max": float(edges[1][-1]),
        "grid_z_min": float(edges[2][0]), "grid_z_max": float(edges[2][-1]),
    }


def analyse_task(row_dict, seed, dt_ms, cfg, output, fingerprint):
    row = pd.Series(row_dict)
    checkpoint = _checkpoint_path(output, row.scenario_id, dt_ms, seed)
    cached = _valid_checkpoint(checkpoint, fingerprint)
    if cached is not None:
        return cached
    edges = _load_grid(_grid_path(output, row.scenario_id), fingerprint)
    p_grid, samples = _simulate_path(row, int(seed), float(dt_ms), cfg)
    masses, retained = build_masses(samples, edges, cfg)
    required = float(cfg["density"].get("required_retention", 1.0))
    tolerance = float(cfg["density"].get("retention_tolerance", 1e-12))
    if float(np.min(retained)) + tolerance < required:
        raise RuntimeError(
            "Full-coverage violation: scenario=%s seed=%s dt=%s retained_min=%.12g required=%.12g" %
            (row.scenario_id, seed, dt_ms, float(np.min(retained)), required)
        )
    payload = _payload(row, seed, dt_ms, p_grid, masses, retained, fingerprint)
    _atomic_json(checkpoint, payload)
    return payload


def _scenario_marker_tables(results, scenarios, cfg):
    p_grid = np.linspace(0.0, 1.0, int(cfg["path"]["n_p"]))
    persistence = int(cfg["markers"]["persistence_points"])
    by_key = {}
    for result in results:
        for view in VIEWS:
            key = (result["scenario_id"], result["biological_pair_key"], result["dt_ms"], view)
            by_key.setdefault(key, []).append((result["seed"], np.asarray(result["views"][view]["delta_kl"], float)))
    rows, curve_cache = [], {}
    for key, items in sorted(by_key.items()):
        items.sort(key=lambda item: item[0])
        curves = np.asarray([item[1] for item in items])
        summary = seed_ensemble_markers(p_grid, curves, persistence)
        curve_cache[key] = {item[0]: item[1] for item in items}
        for variant in MARKER_VARIANTS:
            rows.append({
                "scenario_id": key[0], "biological_pair_key": key[1], "dt_ms": key[2], "view": key[3],
                "marker_variant": variant, "kl_balance_p": summary[variant],
                "seed_isotonic_iqr": summary["seed_isotonic_iqr"],
                "median_curve_crossing_count": summary["median_curve_crossing_count"],
                "endpoint_direction_fraction": summary["endpoint_direction_fraction"],
            })
    frame = pd.DataFrame(rows)
    weights = scenarios[["scenario_id", "analysis_within_pair_weight"]].drop_duplicates()
    frame = frame.merge(weights, on="scenario_id", how="left", validate="many_to_one")
    return frame, curve_cache


def _pair_marker_first(scenario_markers, pair_stage):
    stage = pair_stage[pair_stage.path_family.eq("coupled")].set_index("biological_pair_key")
    keys = ["biological_pair_key", "dt_ms", "view", "marker_variant"]
    rows = []
    for key, group in scenario_markers.groupby(keys, sort=True):
        values, weights = group.kl_balance_p.to_numpy(float), group.analysis_within_pair_weight.to_numpy(float)
        balance = float(stage.loc[key[0], "balance_p_isi_weighted_median"])
        marker = weighted_quantile(values, weights, 0.5)
        rows.append(dict(zip(keys, key), aggregation_order="marker_first", n_support_scenarios=len(group),
                         kl_balance_p=marker,
                         kl_balance_p_q25=weighted_quantile(values, weights, 0.25),
                         kl_balance_p_q75=weighted_quantile(values, weights, 0.75),
                         firing_balance_p=balance, kl_minus_firing_p=marker - balance))
    return pd.DataFrame(rows)


def _pair_curve_first(curve_cache, scenarios, pair_stage, cfg):
    p_grid = np.linspace(0.0, 1.0, int(cfg["path"]["n_p"]))
    persistence = int(cfg["markers"]["persistence_points"])
    stage = pair_stage[pair_stage.path_family.eq("coupled")].set_index("biological_pair_key")
    scenario_rows = scenarios.set_index("scenario_id")
    accum = {}
    for (scenario_id, pair, dt_ms, view), seeds in curve_cache.items():
        weight = float(scenario_rows.loc[int(scenario_id), "analysis_within_pair_weight"])
        for seed, curve in seeds.items():
            key = (pair, dt_ms, view, seed)
            accum[key] = accum.get(key, np.zeros_like(curve)) + weight * curve
    grouped = {}
    for (pair, dt_ms, view, seed), curve in accum.items():
        grouped.setdefault((pair, dt_ms, view), []).append((seed, curve))
    rows = []
    for key, items in sorted(grouped.items()):
        items.sort(key=lambda item: item[0])
        summary = seed_ensemble_markers(p_grid, np.asarray([item[1] for item in items]), persistence)
        balance = float(stage.loc[key[0], "balance_p_isi_weighted_median"])
        for variant in MARKER_VARIANTS:
            marker = float(summary[variant])
            rows.append({"biological_pair_key": key[0], "dt_ms": key[1], "view": key[2],
                         "marker_variant": variant, "aggregation_order": "curve_first",
                         "n_support_scenarios": int((scenarios.biological_pair_key == key[0]).sum()),
                         "kl_balance_p": marker, "kl_balance_p_q25": np.nan, "kl_balance_p_q75": np.nan,
                         "firing_balance_p": balance, "kl_minus_firing_p": marker - balance})
    return pd.DataFrame(rows)


def _scenario_convergence(scenario_markers, cfg):
    gates = cfg["gates"]
    dt_values = sorted((float(value) for value in cfg["convergence"]["dt_ms"]), reverse=True)
    fine, primary, coarse = min(dt_values), float(cfg["convergence"]["primary_dt_ms"]), max(dt_values)
    primary_variant = "seed_median_curve_isotonic"
    index = scenario_markers.set_index(["scenario_id", "view", "dt_ms", "marker_variant"])
    rows = []
    for scenario_id in sorted(scenario_markers.scenario_id.unique()):
        for view in VIEWS:
            def value(dt, variant, column="kl_balance_p"):
                try:
                    return float(index.loc[(scenario_id, view, dt, variant), column])
                except KeyError:
                    return float("nan")
            coarse_value = value(coarse, primary_variant)
            primary_value = value(primary, primary_variant)
            fine_value = value(fine, primary_variant)
            method_values = [value(fine, item) for item in
                             ("seed_median_curve_isotonic", "seed_median_curve_first", "seed_median_curve_persistent")]
            finite_methods = [item for item in method_values if np.isfinite(item)]
            disagreement = max(finite_methods) - min(finite_methods) if finite_methods else np.inf
            seed_iqr = value(fine, primary_variant, "seed_isotonic_iqr")
            endpoint_fraction = value(fine, primary_variant, "endpoint_direction_fraction")
            pass_value = bool(
                np.isfinite([coarse_value, primary_value, fine_value]).all()
                and abs(coarse_value - primary_value) <= float(gates["scenario_max_shift_coarse_primary"])
                and abs(primary_value - fine_value) <= float(gates["scenario_max_shift_primary_fine"])
                and seed_iqr <= float(gates["scenario_max_fine_seed_iqr"])
                and disagreement <= float(gates["scenario_max_fine_method_disagreement"])
                and endpoint_fraction >= float(gates["scenario_min_endpoint_direction_fraction"])
            )
            rows.append({"scenario_id": int(scenario_id), "view": view,
                         "shift_coarse_primary": coarse_value - primary_value,
                         "shift_primary_fine": primary_value - fine_value,
                         "fine_seed_iqr": seed_iqr, "fine_method_disagreement": disagreement,
                         "fine_endpoint_direction_fraction": endpoint_fraction, "pass": pass_value})
    return pd.DataFrame(rows)


def _ensemble_summary(pair_markers, cells, animal_pairs, loo):
    keys = ["dt_ms", "view", "marker_variant", "aggregation_order"]
    rows = []
    for key, group in pair_markers.groupby(keys, sort=True):
        cell = cells
        animal = animal_pairs
        jack = loo
        for column, value in zip(keys, key):
            cell = cell[cell[column] == value]
            animal = animal[animal[column] == value]
            jack = jack[jack[column] == value]
        pair_values = group.kl_minus_firing_p.to_numpy(float)
        cell_values = cell.kl_minus_firing_p.to_numpy(float)
        animal_values = animal.kl_minus_firing_p.to_numpy(float)
        rows.append(dict(zip(keys, key), n_pairs=len(pair_values),
                         median_pair_delta=float(np.nanmedian(pair_values)),
                         q25_pair_delta=float(np.nanquantile(pair_values, .25)),
                         q75_pair_delta=float(np.nanquantile(pair_values, .75)),
                         fraction_pairs_negative=float(np.mean(pair_values < 0)),
                         n_cells=len(cell_values), median_cell_delta=float(np.nanmedian(cell_values)),
                         q25_cell_delta=float(np.nanquantile(cell_values, .25)),
                         q75_cell_delta=float(np.nanquantile(cell_values, .75)),
                         fraction_cells_negative=float(np.mean(cell_values < 0)),
                         n_animal_pairs=len(animal_values), median_animal_pair_delta=float(np.nanmedian(animal_values)),
                         fraction_animal_pairs_negative=float(np.mean(animal_values < 0)),
                         leave_one_animal_out_all_negative=bool(len(jack) and (jack.median_kl_minus_firing_p < 0).all())))
    return pd.DataFrame(rows)


def run(cfg, frozen_dir):
    frozen = resolve_frozen_dir(frozen_dir)
    scenarios_all, pair_stage, _, _ = load_frozen(frozen)
    scenarios = select_scenarios(scenarios_all, pair_stage, cfg)
    output = Path(cfg["output"]["dir"]).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "extent_checkpoints").mkdir(exist_ok=True)
    (output / "checkpoints").mkdir(exist_ok=True)
    (output / "reference_grids").mkdir(exist_ok=True)
    write_input_manifest(frozen, output)
    fingerprint = _fingerprint(cfg, frozen)
    (output / "SCIENTIFIC_FINGERPRINT.txt").write_text(fingerprint + "\n", encoding="utf-8")

    workers = int(cfg["parallel"]["workers"])
    failures = []
    tasks = [(row.to_dict(), int(seed), float(dt_ms))
             for _, row in scenarios.iterrows()
             for dt_ms in cfg["convergence"]["dt_ms"]
             for seed in cfg["convergence"]["seeds"]]
    print("Pilot pass: measuring extrema for %d scenario-dt-seed paths" % len(tasks), flush=True)
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(pilot_extent_task, row, seed, dt_ms, cfg, str(output), fingerprint):
                   (row["scenario_id"], seed, dt_ms) for row, seed, dt_ms in tasks}
        for count, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
            except Exception as error:
                scenario_id, seed, dt_ms = futures[future]
                failures.append({"stage": "pilot_extent", "scenario_id": scenario_id,
                                 "seed": seed, "dt_ms": dt_ms,
                                 "error": str(error), "traceback": traceback.format_exc()})
            if count % 25 == 0 or count == len(futures):
                elapsed = (time.time() - started) / 3600.0
                print("Pilot paths: %d/%d; elapsed %.2f h" % (count, len(futures), elapsed), flush=True)
    if failures:
        _atomic_json(output / "FAILURES.json", failures)
        raise RuntimeError("Pilot-extent stage failed")

    print("Building %d all-task-extrema reference grids" % len(scenarios), flush=True)
    grid_rows = []
    for _, row in scenarios.iterrows():
        try:
            grid_rows.append(build_scenario_grid(row, cfg, str(output), fingerprint))
        except Exception as error:
            failures.append({"stage": "reference_grid", "scenario_id": int(row.scenario_id),
                             "error": str(error), "traceback": traceback.format_exc()})
    if failures:
        _atomic_json(output / "FAILURES.json", failures)
        raise RuntimeError("Full-coverage reference-grid stage failed")

    print("Analysis pass: running %d scenario-dt-seed paths" % len(tasks), flush=True)
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(analyse_task, row, seed, dt_ms, cfg, str(output), fingerprint):
                   (row["scenario_id"], seed, dt_ms) for row, seed, dt_ms in tasks}
        for count, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
            except Exception as error:
                scenario_id, seed, dt_ms = futures[future]
                failures.append({"stage": "path", "scenario_id": scenario_id, "seed": seed, "dt_ms": dt_ms,
                                 "error": str(error), "traceback": traceback.format_exc()})
            if count % 25 == 0 or count == len(futures):
                elapsed = (time.time() - started) / 3600.0
                print("Analysis paths: %d/%d; elapsed %.2f h" % (count, len(futures), elapsed), flush=True)
    if failures:
        _atomic_json(output / "FAILURES.json", failures)
        raise RuntimeError("One or more path tasks failed")

    results = []
    for _, row in scenarios.iterrows():
        for dt_ms in cfg["convergence"]["dt_ms"]:
            for seed in cfg["convergence"]["seeds"]:
                path = _checkpoint_path(output, row.scenario_id, dt_ms, seed)
                payload = _valid_checkpoint(path, fingerprint)
                if payload is None:
                    raise RuntimeError("Missing valid checkpoint: " + str(path))
                results.append(payload)

    scenario_markers, curve_cache = _scenario_marker_tables(results, scenarios, cfg)
    pair_first = _pair_marker_first(scenario_markers, pair_stage)
    pair_curve = _pair_curve_first(curve_cache, scenarios, pair_stage, cfg)
    pair_markers = pd.concat([pair_first, pair_curve], ignore_index=True)
    cells = cell_balanced(pair_markers)
    mapping = pd.read_csv(cfg["animal_mapping"]["path"])
    animal_pairs = animal_pair_balanced(pair_markers, mapping)
    loo = leave_one_animal_out(animal_pairs)
    scenario_qc = _scenario_convergence(scenario_markers, cfg)
    ensemble = _ensemble_summary(pair_markers, cells, animal_pairs, loo)
    retention = pd.DataFrame([{"scenario_id": item["scenario_id"], "seed": item["seed"], "dt_ms": item["dt_ms"],
                               "retained_min": item["retained_min"], "retained_median": item["retained_median"]}
                              for item in results])
    pilot_rows = []
    for _, row in scenarios.iterrows():
        for dt_ms in cfg["convergence"]["dt_ms"]:
            for seed in cfg["convergence"]["seeds"]:
                payload = _valid_extent(
                    _extent_checkpoint_path(output, row.scenario_id, dt_ms, seed), fingerprint
                )
                if payload is None:
                    raise RuntimeError("Pilot extent disappeared before finalization")
                pilot_rows.append({
                    "scenario_id": payload["scenario_id"], "seed": payload["seed"], "dt_ms": payload["dt_ms"],
                    "x_min": payload["axis_min"][0], "x_max": payload["axis_max"][0],
                    "y_min": payload["axis_min"][1], "y_max": payload["axis_max"][1],
                    "z_min": payload["axis_min"][2], "z_max": payload["axis_max"][2],
                    "sample_count": payload["sample_count"],
                })
    pilot_extents = pd.DataFrame(pilot_rows)
    grid_audit = pd.DataFrame(grid_rows)

    tables = {
        "scenario_markers.csv": scenario_markers,
        "scenario_convergence_gates.csv": scenario_qc,
        "pair_markers_both_orders.csv": pair_markers,
        "cell_balanced_deltas.csv": cells,
        "animal_pair_balanced_deltas.csv": animal_pairs,
        "leave_one_animal_out.csv": loo,
        "ensemble_convergence_summary.csv": ensemble,
        "grid_retention.csv": retention,
        "pilot_extents.csv": pilot_extents,
        "reference_grid_bounds.csv": grid_audit,
    }
    for name, frame in tables.items():
        frame.to_csv(output / name, index=False)
    verdict = build_verdict(cfg, scenario_qc, ensemble, pair_markers, retention)
    _atomic_json(output / "KL_CONVERGENCE_VERDICT.json", verdict)
    required_retention = float(cfg["density"].get("required_retention", 1.0))
    retention_tolerance = float(cfg["density"].get("retention_tolerance", 1e-12))
    coverage_audit = {
        "version": "1.0.1",
        "grid_strategy": str(cfg["density"].get("grid_strategy", "all_task_extrema")),
        "pilot_extent_tasks": len(pilot_extents),
        "reference_grids": len(grid_audit),
        "analysis_tasks": len(retention),
        "minimum_retained_mass": float(retention.retained_min.min()),
        "required_retained_mass": required_retention,
        "retention_tolerance": retention_tolerance,
        "all_tasks_full_coverage": bool(
            (retention.retained_min + retention_tolerance >= required_retention).all()
        ),
    }
    _atomic_json(output / "GRID_COVERAGE_AUDIT.json", coverage_audit)
    status = {"KEEP_AS_MAIN_RESULT": "PASS",
              "KEEP_AS_ENSEMBLE_RESULT_WITH_LIMITATIONS": "LIMITED",
              "REMOVE_KL_RESULT": "FAIL",
              "SMOKE_ONLY": "SMOKE"}[verdict["decision"]]
    analysis_simulations = (len(scenarios) * len(cfg["convergence"]["dt_ms"]) *
                            len(cfg["convergence"]["seeds"]) * int(cfg["path"]["n_p"]))
    summary = {
        "version": "1.0.1", "status": status,
        "selected_scenarios": len(scenarios), "dependent_pairs": int(scenarios.biological_pair_key.nunique()),
        "dt_ms": [float(value) for value in cfg["convergence"]["dt_ms"]],
        "seeds": [int(value) for value in cfg["convergence"]["seeds"]],
        "path_positions": int(cfg["path"]["n_p"]),
        "pilot_stationary_state_simulations": analysis_simulations,
        "analysis_stationary_state_simulations": analysis_simulations,
        "stationary_state_simulations": 2 * analysis_simulations,
        "grid_coverage": coverage_audit,
        "scientific_fingerprint": fingerprint, "verdict": verdict,
    }
    _atomic_json(output / "RUN_SUMMARY.json", summary)
    return output, summary
