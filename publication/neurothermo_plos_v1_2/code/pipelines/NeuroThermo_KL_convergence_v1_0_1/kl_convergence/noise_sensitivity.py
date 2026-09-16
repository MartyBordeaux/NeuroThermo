from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import traceback

import numpy as np
import pandas as pd

from .aggregation import animal_pair_balanced, cell_balanced, leave_one_animal_out
from .cli import load_config, validate
from .data import REQUIRED, load_frozen, resolve_frozen_dir, select_scenarios, sha256
from .density import build_masses
from .model import stationary_samples_nested
from .pipeline import (
    VIEWS,
    _endpoint_curves,
    _ensemble_summary,
    _pair_curve_first,
    _pair_marker_first,
    _scenario_marker_tables,
    _task_seed,
)

VERSION = "1.0.0"
CONDITIONS = (("half", 0.5), ("double", 2.0))
PRIMARY_MARKER = "seed_median_curve_isotonic"


def _atomic_json(path: Path, payload: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    os.replace(str(temporary), str(path))


def _condition_fingerprint(cfg: dict, frozen_dir: Path, grid_bounds_path: Path,
                           mapping_path: Path, label: str, multiplier: float) -> str:
    payload = {
        "version": VERSION,
        "label": str(label),
        "multiplier": float(multiplier),
        "primary_dt_ms": float(cfg["convergence"]["primary_dt_ms"]),
        "seeds": [int(v) for v in cfg["convergence"]["seeds"]],
        "path_n_p": int(cfg["path"]["n_p"]),
        "model": cfg["model"],
        "noise_D": [float(v) for v in cfg["noise"]["D"]],
        "stationary": cfg["stationary"],
        "density": cfg["density"],
        "markers": cfg["markers"],
        "inputs": {name: sha256(frozen_dir / name) for name in REQUIRED},
        "grid_bounds_sha256": sha256(grid_bounds_path),
        "animal_mapping_sha256": sha256(mapping_path),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _checkpoint_path(output: Path, label: str, scenario_id: int, seed: int) -> Path:
    return output / "checkpoints" / label / f"scenario_{int(scenario_id):05d}_seed_{int(seed)}.json"


def _load_checkpoint(path: Path, fingerprint: str):
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if payload.get("_status") != "OK" or payload.get("_fingerprint") != fingerprint:
        return None
    return payload


def _edges_from_bounds(bounds_row: pd.Series, bins: int):
    return [
        np.linspace(float(bounds_row.grid_x_min), float(bounds_row.grid_x_max), bins + 1),
        np.linspace(float(bounds_row.grid_y_min), float(bounds_row.grid_y_max), bins + 1),
        np.linspace(float(bounds_row.grid_z_min), float(bounds_row.grid_z_max), bins + 1),
    ]


def _simulate_task(row_dict: dict, seed: int, multiplier: float, label: str,
                   cfg: dict, bounds_dict: dict, output: str, fingerprint: str):
    row = pd.Series(row_dict)
    scenario_id = int(row.scenario_id)
    output_path = Path(output)
    checkpoint = _checkpoint_path(output_path, label, scenario_id, int(seed))
    cached = _load_checkpoint(checkpoint, fingerprint)
    if cached is not None:
        return cached

    local_cfg = deepcopy(cfg)
    local_cfg["noise"]["multiplier"] = float(multiplier)
    primary_dt = float(local_cfg["convergence"]["primary_dt_ms"])
    p_grid = np.linspace(0.0, 1.0, int(local_cfg["path"]["n_p"]))

    from .data import scenario_arrays
    theta, current = scenario_arrays(row, p_grid)
    samples = []
    start = None
    for p_index in range(len(p_grid)):
        values, final, ok = stationary_samples_nested(
            _task_seed(int(seed), scenario_id, p_index),
            theta[p_index], current[p_index], local_cfg, primary_dt, start=start,
        )
        if not ok or len(values) < 10:
            raise RuntimeError(
                f"Non-finite/short run: condition={label} scenario={scenario_id} "
                f"seed={seed} p={p_grid[p_index]:.8g}"
            )
        samples.append(values)
        start = final

    edges = _edges_from_bounds(pd.Series(bounds_dict), int(local_cfg["density"]["bins"]))
    masses, retained = build_masses(samples, edges, local_cfg)
    required = float(local_cfg["density"].get("required_retention", 1.0))
    tolerance = float(local_cfg["density"].get("retention_tolerance", 1e-12))
    retained_min = float(np.min(retained))
    if retained_min + tolerance < required:
        raise RuntimeError(
            f"Frozen-grid retention failure: condition={label} scenario={scenario_id} "
            f"seed={seed} retained_min={retained_min:.12g} required={required:.12g}. "
            "Do not clip or silently enlarge the grid; construct a prespecified union-grid sensitivity instead."
        )

    payload = {
        "_status": "OK",
        "_fingerprint": fingerprint,
        "condition": label,
        "noise_multiplier": float(multiplier),
        "scenario_id": scenario_id,
        "biological_pair_key": str(row.biological_pair_key),
        "seed": int(seed),
        "dt_ms": primary_dt,
        "retained_min": retained_min,
        "retained_median": float(np.median(retained)),
        "views": _endpoint_curves(masses),
    }
    _atomic_json(checkpoint, payload)
    return payload


def _resolve_mapping(cfg: dict) -> Path:
    path = Path(cfg["animal_mapping"]["path"])
    if not path.is_absolute():
        config_path = Path(cfg["_config_path"])
        path = (config_path.parent.parent / path).resolve()
        cfg["animal_mapping"]["path"] = str(path)
    return path


def _summarize_condition(results, scenarios, pair_stage, mapping, cfg, label, multiplier):
    scenario_markers, curve_cache = _scenario_marker_tables(results, scenarios, cfg)
    pair_first = _pair_marker_first(scenario_markers, pair_stage)
    pair_curve = _pair_curve_first(curve_cache, scenarios, pair_stage, cfg)
    pair_markers = pd.concat([pair_first, pair_curve], ignore_index=True)
    cells = cell_balanced(pair_markers)
    animal_pairs = animal_pair_balanced(pair_markers, mapping)
    loo = leave_one_animal_out(animal_pairs)
    ensemble = _ensemble_summary(pair_markers, cells, animal_pairs, loo)
    for frame in (scenario_markers, pair_markers, cells, animal_pairs, loo, ensemble):
        frame.insert(0, "noise_condition", label)
        frame.insert(1, "noise_multiplier", float(multiplier))
    return scenario_markers, pair_markers, cells, animal_pairs, loo, ensemble


def _primary_rows(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.loc[
        (frame["marker_variant"] == PRIMARY_MARKER)
        & frame["view"].isin(VIEWS)
    ].copy()


def _make_verdict(summary: pd.DataFrame) -> dict:
    primary = _primary_rows(summary)
    records = []
    for order in ("marker_first", "curve_first"):
        subset = primary[primary.aggregation_order.eq(order)].set_index(["noise_multiplier", "view"])
        for multiplier in (0.5, 1.0, 2.0):
            record = {"noise_multiplier": multiplier, "aggregation_order": order}
            for view in VIEWS:
                row = subset.loc[(multiplier, view)]
                record[f"{view}_median_pair_delta"] = float(row.median_pair_delta)
                record[f"{view}_fraction_pairs_negative"] = float(row.fraction_pairs_negative)
            record["kl_before_firing_full_state"] = bool(record["xyz_median_pair_delta"] < 0.0)
            record["kl_before_firing_slow_state"] = bool(record["z_median_pair_delta"] < 0.0)
            record["slow_more_negative_than_fast"] = bool(
                record["z_median_pair_delta"] < record["xy_median_pair_delta"]
            )
            records.append(record)
    return {
        "version": VERSION,
        "primary_dt_ms": float(primary.dt_ms.iloc[0]),
        "noise_multipliers": [0.5, 1.0, 2.0],
        "primary_marker": PRIMARY_MARKER,
        "ordering_definition": "median(kl_balance_p - firing_balance_p) < 0",
        "kl_before_firing_preserved_all_noise_and_orders": all(
            item["kl_before_firing_full_state"] for item in records
        ),
        "slow_state_ordering_preserved_all_noise_and_orders": all(
            item["kl_before_firing_slow_state"] for item in records
        ),
        "slow_state_dominance_operational_definition": "z median pair delta < xy median pair delta",
        "slow_state_dominance_preserved_all_noise_and_orders": all(
            item["slow_more_negative_than_fast"] for item in records
        ),
        "rows": records,
    }


def run(config_path, frozen_dir, baseline_dir, output_dir):
    cfg = load_config(config_path)
    checks = validate(cfg, frozen_dir)
    primary_dt = float(cfg["convergence"]["primary_dt_ms"])
    seeds = [int(value) for value in cfg["convergence"]["seeds"]]
    if primary_dt != 0.025:
        raise ValueError(f"Frozen primary dt mismatch: {primary_dt} != 0.025 ms")
    if seeds != [20260818, 21260821, 22260823, 23260837, 24260855]:
        raise ValueError("Frozen seed sequence mismatch")

    frozen = resolve_frozen_dir(frozen_dir)
    scenarios_all, pair_stage, _, _ = load_frozen(frozen)
    scenarios = select_scenarios(scenarios_all, pair_stage, cfg)
    if len(scenarios) != 264 or scenarios.biological_pair_key.nunique() != 32:
        raise ValueError("Noise sensitivity must use the frozen 264-scenario / 32-pair cohort")

    baseline_dir = Path(baseline_dir).expanduser().resolve()
    grid_bounds_path = baseline_dir / "reference_grid_bounds.csv"
    baseline_summary_path = baseline_dir / "ensemble_convergence_summary.csv"
    if not grid_bounds_path.is_file() or not baseline_summary_path.is_file():
        raise FileNotFoundError(
            "Baseline KL results must contain reference_grid_bounds.csv and ensemble_convergence_summary.csv"
        )

    bounds = pd.read_csv(grid_bounds_path).set_index("scenario_id")
    missing_bounds = sorted(set(scenarios.scenario_id.astype(int)) - set(bounds.index.astype(int)))
    if missing_bounds:
        raise ValueError(f"Missing frozen reference-grid bounds for scenario IDs: {missing_bounds}")

    mapping_path = _resolve_mapping(cfg)
    mapping = pd.read_csv(mapping_path)
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    all_ensemble = []
    failures = []
    workers = int(cfg["parallel"]["workers"])

    for label, multiplier in CONDITIONS:
        fingerprint = _condition_fingerprint(
            cfg, frozen, grid_bounds_path, mapping_path, label, multiplier
        )
        tasks = [(row.to_dict(), seed) for _, row in scenarios.iterrows() for seed in seeds]
        results = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            future_map = {
                pool.submit(
                    _simulate_task,
                    row_dict,
                    seed,
                    multiplier,
                    label,
                    cfg,
                    bounds.loc[int(row_dict["scenario_id"])].to_dict(),
                    str(output),
                    fingerprint,
                ): (int(row_dict["scenario_id"]), seed)
                for row_dict, seed in tasks
            }
            for future in as_completed(future_map):
                scenario_id, seed = future_map[future]
                try:
                    results.append(future.result())
                except Exception as error:
                    failures.append({
                        "noise_condition": label,
                        "noise_multiplier": multiplier,
                        "scenario_id": scenario_id,
                        "seed": seed,
                        "error": str(error),
                        "traceback": traceback.format_exc(),
                    })
        if failures:
            _atomic_json(output / "FAILURES.json", failures)
            raise RuntimeError(
                "Noise sensitivity failed. Inspect FAILURES.json; no verdict was generated."
            )

        results.sort(key=lambda item: (item["scenario_id"], item["seed"]))
        scenario_markers, pair_markers, cells, animal_pairs, loo, ensemble = _summarize_condition(
            results, scenarios, pair_stage, mapping, cfg, label, multiplier
        )
        scenario_markers.to_csv(output / f"scenario_markers_{label}.csv", index=False)
        pair_markers.to_csv(output / f"pair_markers_{label}.csv", index=False)
        cells.to_csv(output / f"cell_balanced_{label}.csv", index=False)
        animal_pairs.to_csv(output / f"animal_pair_balanced_{label}.csv", index=False)
        loo.to_csv(output / f"leave_one_animal_out_{label}.csv", index=False)
        ensemble.to_csv(output / f"ensemble_summary_{label}.csv", index=False)
        all_ensemble.append(ensemble)

    baseline = pd.read_csv(baseline_summary_path)
    baseline = baseline.loc[np.isclose(baseline.dt_ms.astype(float), primary_dt)].copy()
    baseline.insert(0, "noise_condition", "reference")
    baseline.insert(1, "noise_multiplier", 1.0)
    combined = pd.concat([baseline] + all_ensemble, ignore_index=True, sort=False)
    combined.to_csv(output / "noise_sensitivity_ensemble_summary.csv", index=False)

    primary = _primary_rows(combined).sort_values(
        ["aggregation_order", "noise_multiplier", "view"]
    )
    primary.to_csv(output / "noise_sensitivity_primary_ordering.csv", index=False)
    verdict = _make_verdict(combined)
    verdict.update({
        "selected_scenarios": int(checks["selected_scenarios"]),
        "selected_pairs": int(checks["selected_pairs"]),
        "seeds": seeds,
        "base_D": [float(v) for v in cfg["noise"]["D"]],
        "effective_D": {
            "half": (np.asarray(cfg["noise"]["D"], float) * 0.5).tolist(),
            "reference": np.asarray(cfg["noise"]["D"], float).tolist(),
            "double": (np.asarray(cfg["noise"]["D"], float) * 2.0).tolist(),
        },
        "grid_policy": "reuse frozen reference_grid_bounds.csv; require 100% sample retention",
        "common_random_numbers": True,
        "new_stationary_state_simulations": int(
            len(scenarios) * len(seeds) * len(CONDITIONS) * int(cfg["path"]["n_p"])
        ),
    })
    _atomic_json(output / "NOISE_SENSITIVITY_VERDICT.json", verdict)
    return output, verdict


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Minimal primary-dt diffusion sensitivity for the frozen KL-before-firing result"
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--frozen-dir", required=True)
    parser.add_argument("--baseline-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    output, verdict = run(args.config, args.frozen_dir, args.baseline_dir, args.out)
    print(json.dumps(verdict, indent=2, allow_nan=True))
    print(f"Results written to {output}")


if __name__ == "__main__":
    main()
