from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import time
import traceback

import numpy as np
import pandas as pd

from .aggregation import animal_pair_balanced, cell_balanced, leave_one_animal_out
from .cli import load_config, validate
from .data import load_frozen, resolve_frozen_dir, select_scenarios
from .density import full_coverage_grid
from .pipeline import (
    _atomic_json,
    _ensemble_summary,
    _extent_checkpoint_path,
    _fingerprint,
    _grid_path,
    _load_grid,
    _pair_curve_first,
    _pair_marker_first,
    _save_grid,
    _scenario_marker_tables,
    _valid_checkpoint,
    _valid_extent,
    analyse_task,
    pilot_extent_task,
)

VERSION = "1.0.0"
PRIMARY_MARKER = "seed_median_curve_isotonic"
PRIMARY_DT_EXPECTED = 0.025
SEEDS_EXPECTED = [20260818, 21260821, 22260823, 23260837, 24260855]


def _label(multiplier: float) -> str:
    if np.isclose(multiplier, 0.5):
        return "half"
    if np.isclose(multiplier, 2.0):
        return "double"
    return ("x%.6g" % float(multiplier)).replace(".", "p")


def _build_primary_grid(row, cfg, output, fingerprint):
    primary_dt = float(cfg["convergence"]["primary_dt_ms"])
    extents = []
    for seed in cfg["convergence"]["seeds"]:
        path = _extent_checkpoint_path(output, row.scenario_id, primary_dt, seed)
        payload = _valid_extent(path, fingerprint)
        if payload is None:
            raise RuntimeError("Missing valid primary-dt pilot extent: " + str(path))
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


def _resolve_mapping(cfg):
    path = Path(cfg["animal_mapping"]["path"])
    if not path.is_absolute():
        path = (Path(cfg["_config_path"]).parent.parent / path).resolve()
        cfg["animal_mapping"]["path"] = str(path)
    return path


def _primary_comparison(condition_summary, baseline_summary, primary_dt, multiplier, label):
    keys = ["view", "marker_variant", "aggregation_order"]
    base = baseline_summary.loc[np.isclose(baseline_summary.dt_ms.astype(float), primary_dt)].copy()
    cond = condition_summary.copy()
    base = base[keys + [
        "median_pair_delta", "q25_pair_delta", "q75_pair_delta", "fraction_pairs_negative",
        "median_cell_delta", "fraction_cells_negative", "median_animal_pair_delta",
        "fraction_animal_pairs_negative", "leave_one_animal_out_all_negative",
    ]].rename(columns={c: c + "__baseline" for c in [
        "median_pair_delta", "q25_pair_delta", "q75_pair_delta", "fraction_pairs_negative",
        "median_cell_delta", "fraction_cells_negative", "median_animal_pair_delta",
        "fraction_animal_pairs_negative", "leave_one_animal_out_all_negative",
    ]})
    cond = cond[keys + [
        "median_pair_delta", "q25_pair_delta", "q75_pair_delta", "fraction_pairs_negative",
        "median_cell_delta", "fraction_cells_negative", "median_animal_pair_delta",
        "fraction_animal_pairs_negative", "leave_one_animal_out_all_negative",
    ]].rename(columns={c: c + "__condition" for c in [
        "median_pair_delta", "q25_pair_delta", "q75_pair_delta", "fraction_pairs_negative",
        "median_cell_delta", "fraction_cells_negative", "median_animal_pair_delta",
        "fraction_animal_pairs_negative", "leave_one_animal_out_all_negative",
    ]})
    merged = base.merge(cond, on=keys, how="inner", validate="one_to_one")
    merged.insert(0, "noise_condition", label)
    merged.insert(1, "noise_multiplier", float(multiplier))
    merged.insert(2, "dt_ms", float(primary_dt))
    merged["median_pair_delta_shift"] = merged["median_pair_delta__condition"] - merged["median_pair_delta__baseline"]
    merged["fraction_pairs_negative_shift"] = merged["fraction_pairs_negative__condition"] - merged["fraction_pairs_negative__baseline"]
    return merged


def _verdict(condition_summary, comparison, cfg, multiplier, label):
    primary = condition_summary.loc[condition_summary.marker_variant.eq(PRIMARY_MARKER)].copy()
    rows = []
    for order in ("marker_first", "curve_first"):
        sub = primary.loc[primary.aggregation_order.eq(order)].set_index("view")
        rec = {"aggregation_order": order}
        for view in ("xyz", "xy", "z"):
            row = sub.loc[view]
            rec[view + "_median_pair_delta"] = float(row.median_pair_delta)
            rec[view + "_fraction_pairs_negative"] = float(row.fraction_pairs_negative)
            rec[view + "_median_cell_delta"] = float(row.median_cell_delta)
            rec[view + "_fraction_cells_negative"] = float(row.fraction_cells_negative)
            rec[view + "_median_animal_pair_delta"] = float(row.median_animal_pair_delta)
            rec[view + "_fraction_animal_pairs_negative"] = float(row.fraction_animal_pairs_negative)
            rec[view + "_loo_all_negative"] = bool(row.leave_one_animal_out_all_negative)
        rec["kl_before_firing_full_state"] = bool(rec["xyz_median_pair_delta"] < 0.0)
        rec["kl_before_firing_slow_state"] = bool(rec["z_median_pair_delta"] < 0.0)
        rec["slow_more_negative_than_fast"] = bool(rec["z_median_pair_delta"] < rec["xy_median_pair_delta"])
        rows.append(rec)
    return {
        "version": VERSION,
        "noise_condition": label,
        "noise_multiplier": float(multiplier),
        "base_D": [float(v) for v in cfg["noise"]["D"]],
        "effective_D": (np.asarray(cfg["noise"]["D"], float) * float(multiplier)).tolist(),
        "primary_dt_ms": float(cfg["convergence"]["primary_dt_ms"]),
        "seeds": [int(v) for v in cfg["convergence"]["seeds"]],
        "primary_marker": PRIMARY_MARKER,
        "ordering_definition": "median(kl_balance_p - firing_balance_p) < 0",
        "slow_state_dominance_definition": "z median pair delta < xy median pair delta",
        "grid_policy": "condition-specific five-seed primary-dt all-extrema grid with frozen 2% margin and 100% retention",
        "full_state_ordering_preserved_both_aggregation_orders": all(r["kl_before_firing_full_state"] for r in rows),
        "slow_state_ordering_preserved_both_aggregation_orders": all(r["kl_before_firing_slow_state"] for r in rows),
        "slow_state_dominance_preserved_both_aggregation_orders": all(r["slow_more_negative_than_fast"] for r in rows),
        "rows": rows,
        "comparison_rows": int(len(comparison)),
    }


def run(config_path, frozen_dir, baseline_dir, output_dir, multiplier):
    multiplier = float(multiplier)
    if multiplier <= 0:
        raise ValueError("Noise multiplier must be positive")

    cfg = load_config(config_path)
    checks = validate(cfg, frozen_dir)
    primary_dt = float(cfg["convergence"]["primary_dt_ms"])
    seeds = [int(v) for v in cfg["convergence"]["seeds"]]
    if not np.isclose(primary_dt, PRIMARY_DT_EXPECTED, rtol=0.0, atol=1e-15):
        raise ValueError("Frozen primary dt mismatch")
    if seeds != SEEDS_EXPECTED:
        raise ValueError("Frozen seed sequence mismatch")
    if checks["selected_scenarios"] != 264 or checks["selected_pairs"] != 32:
        raise ValueError("Frozen cohort mismatch")

    frozen = resolve_frozen_dir(frozen_dir)
    scenarios_all, pair_stage, _, _ = load_frozen(frozen)
    scenarios = select_scenarios(scenarios_all, pair_stage, cfg)
    mapping_path = _resolve_mapping(cfg)
    mapping = pd.read_csv(mapping_path)

    label = _label(multiplier)
    cfg["noise"]["multiplier"] = multiplier
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "extent_checkpoints").mkdir(exist_ok=True)
    (output / "checkpoints").mkdir(exist_ok=True)
    (output / "reference_grids").mkdir(exist_ok=True)
    fingerprint = _fingerprint(cfg, frozen)
    (output / "SCIENTIFIC_FINGERPRINT.txt").write_text(fingerprint + "\n", encoding="utf-8")

    workers = int(cfg["parallel"]["workers"])
    tasks = [(row.to_dict(), seed) for _, row in scenarios.iterrows() for seed in seeds]
    failures = []

    print("Noise %s pilot pass: %d scenario-seed paths" % (label, len(tasks)), flush=True)
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(pilot_extent_task, row, seed, primary_dt, cfg, str(output), fingerprint):
            (int(row["scenario_id"]), int(seed)) for row, seed in tasks
        }
        for count, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
            except Exception as error:
                sid, seed = futures[future]
                failures.append({"stage": "pilot", "scenario_id": sid, "seed": seed,
                                 "error": str(error), "traceback": traceback.format_exc()})
            if count % 50 == 0 or count == len(futures):
                print("Pilot %d/%d; elapsed %.2f h" % (count, len(futures), (time.time()-started)/3600.0), flush=True)
    if failures:
        _atomic_json(output / "FAILURES.json", failures)
        raise RuntimeError("Pilot stage failed")

    grid_rows = []
    for _, row in scenarios.iterrows():
        grid_rows.append(_build_primary_grid(row, cfg, str(output), fingerprint))
    pd.DataFrame(grid_rows).to_csv(output / "reference_grid_bounds.csv", index=False)

    print("Noise %s analysis pass: %d scenario-seed paths" % (label, len(tasks)), flush=True)
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(analyse_task, row, seed, primary_dt, cfg, str(output), fingerprint):
            (int(row["scenario_id"]), int(seed)) for row, seed in tasks
        }
        for count, future in enumerate(as_completed(futures), 1):
            try:
                future.result()
            except Exception as error:
                sid, seed = futures[future]
                failures.append({"stage": "analysis", "scenario_id": sid, "seed": seed,
                                 "error": str(error), "traceback": traceback.format_exc()})
            if count % 50 == 0 or count == len(futures):
                print("Analysis %d/%d; elapsed %.2f h" % (count, len(futures), (time.time()-started)/3600.0), flush=True)
    if failures:
        _atomic_json(output / "FAILURES.json", failures)
        raise RuntimeError("Analysis stage failed")

    results = []
    retention_rows = []
    for _, row in scenarios.iterrows():
        for seed in seeds:
            path = Path(output) / "checkpoints" / ("scenario_%05d_dt_0p025_seed_%d.json" % (int(row.scenario_id), int(seed)))
            payload = _valid_checkpoint(path, fingerprint)
            if payload is None:
                raise RuntimeError("Missing valid checkpoint: " + str(path))
            results.append(payload)
            retention_rows.append({"scenario_id": int(row.scenario_id), "seed": int(seed),
                                   "retained_min": float(payload["retained_min"]),
                                   "retained_median": float(payload["retained_median"])})
    pd.DataFrame(retention_rows).to_csv(output / "grid_retention.csv", index=False)

    scenario_markers, curve_cache = _scenario_marker_tables(results, scenarios, cfg)
    pair_first = _pair_marker_first(scenario_markers, pair_stage)
    pair_curve = _pair_curve_first(curve_cache, scenarios, pair_stage, cfg)
    pair_markers = pd.concat([pair_first, pair_curve], ignore_index=True)
    cells = cell_balanced(pair_markers)
    animal_pairs = animal_pair_balanced(pair_markers, mapping)
    loo = leave_one_animal_out(animal_pairs)
    ensemble = _ensemble_summary(pair_markers, cells, animal_pairs, loo)

    scenario_markers.to_csv(output / "scenario_markers.csv", index=False)
    pair_markers.to_csv(output / "pair_markers.csv", index=False)
    cells.to_csv(output / "cell_balanced_deltas.csv", index=False)
    animal_pairs.to_csv(output / "animal_pair_balanced_deltas.csv", index=False)
    loo.to_csv(output / "leave_one_animal_out.csv", index=False)
    ensemble.to_csv(output / "ensemble_summary.csv", index=False)

    baseline_path = Path(baseline_dir).expanduser().resolve() / "ensemble_convergence_summary.csv"
    baseline = pd.read_csv(baseline_path)
    comparison = _primary_comparison(ensemble, baseline, primary_dt, multiplier, label)
    comparison.to_csv(output / "condition_vs_baseline.csv", index=False)

    verdict = _verdict(ensemble, comparison, cfg, multiplier, label)
    verdict.update({
        "selected_scenarios": 264,
        "selected_pairs": 32,
        "pilot_stationary_state_simulations": int(len(tasks) * int(cfg["path"]["n_p"])),
        "analysis_stationary_state_simulations": int(len(tasks) * int(cfg["path"]["n_p"])),
        "total_stationary_state_simulations": int(2 * len(tasks) * int(cfg["path"]["n_p"])),
        "minimum_grid_retention": float(pd.DataFrame(retention_rows).retained_min.min()),
    })
    _atomic_json(output / "NOISE_SENSITIVITY_CONDITION.json", verdict)
    print(json.dumps(verdict, indent=2, allow_nan=True), flush=True)
    return output, verdict


def main(argv=None):
    parser = argparse.ArgumentParser(description="Primary-dt single-condition diffusion sensitivity")
    parser.add_argument("--config", required=True)
    parser.add_argument("--frozen-dir", required=True)
    parser.add_argument("--baseline-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--multiplier", required=True, type=float)
    args = parser.parse_args(argv)
    run(args.config, args.frozen_dir, args.baseline_dir, args.out, args.multiplier)


if __name__ == "__main__":
    main()
