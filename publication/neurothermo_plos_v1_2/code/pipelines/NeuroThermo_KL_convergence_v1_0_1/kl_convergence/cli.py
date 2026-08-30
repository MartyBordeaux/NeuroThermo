from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .data import load_frozen, resolve_frozen_dir, select_scenarios
from .figures import make_figures
from .pipeline import run


def load_config(path):
    path = Path(path).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    cfg["_config_path"] = str(path)
    return cfg


def validate(cfg, frozen_dir):
    root = resolve_frozen_dir(frozen_dir)
    scenarios, pair_stage, _, _ = load_frozen(root)
    selected = select_scenarios(scenarios, pair_stage, cfg)
    expected_scenarios = int(cfg["validation"]["expected_scenarios"])
    expected_pairs = int(cfg["validation"]["expected_pairs"])
    if len(selected) != expected_scenarios:
        raise ValueError("Selected %d scenarios; expected %d" % (len(selected), expected_scenarios))
    pairs = int(selected.biological_pair_key.nunique())
    if pairs != expected_pairs:
        raise ValueError("Selected %d pairs; expected %d" % (pairs, expected_pairs))
    dt_values = [float(value) for value in cfg["convergence"]["dt_ms"]]
    if len(dt_values) != len(set(dt_values)) or min(dt_values) <= 0:
        raise ValueError("convergence.dt_ms must contain distinct positive values")
    primary = float(cfg["convergence"]["primary_dt_ms"])
    if primary not in dt_values:
        raise ValueError("primary_dt_ms must be included in convergence.dt_ms")
    finest = min(dt_values)
    for value in dt_values:
        ratio = value / finest
        if not np.isclose(ratio, round(ratio), rtol=0.0, atol=1e-12):
            raise ValueError("Every dt must be an integer multiple of the finest dt")
    seeds = [int(value) for value in cfg["convergence"]["seeds"]]
    if len(seeds) < 3 or len(seeds) != len(set(seeds)):
        raise ValueError("At least three distinct seeds are required")
    if int(cfg["path"]["n_p"]) < 9:
        raise ValueError("At least nine path positions are required")
    for value in cfg["noise"]["D"]:
        if float(value) <= 0:
            raise ValueError("All diffusion coefficients must be positive")
    density = cfg["density"]
    if density.get("grid_strategy") != "all_task_extrema":
        raise ValueError("density.grid_strategy must be all_task_extrema")
    if float(density.get("coverage_margin_fraction", 0.0)) <= 0:
        raise ValueError("density.coverage_margin_fraction must be positive")
    if not np.isclose(float(density.get("required_retention", 0.0)), 1.0,
                      rtol=0.0, atol=0.0):
        raise ValueError("density.required_retention must be exactly 1.0")
    mapping_path = Path(cfg["animal_mapping"]["path"])
    if not mapping_path.is_absolute():
        mapping_path = Path(cfg["_config_path"]).parent.parent / mapping_path
        cfg["animal_mapping"]["path"] = str(mapping_path.resolve())
    mapping = pd.read_csv(mapping_path)
    cells = set()
    for key in selected.biological_pair_key.unique():
        left, right = str(key).split("__TO__", 1)
        cells.update((left, right))
    missing = sorted(cells - set(mapping.cell_id.astype(str)))
    if missing:
        raise ValueError("Animal mapping lacks selected cells: " + ", ".join(missing))
    if cfg.get("analysis", {}).get("role") == "server_full_ensemble":
        frozen_values = {
            "path.n_p": (int(cfg["path"]["n_p"]), 31),
            "stationary.burn_ms": (float(cfg["stationary"]["burn_ms"]), 2400.0),
            "stationary.sample_ms": (float(cfg["stationary"]["sample_ms"]), 6000.0),
            "stationary.sample_stride_ms": (float(cfg["stationary"]["sample_stride_ms"]), 0.5),
            "density.bins": (int(cfg["density"]["bins"]), 22),
        }
        for label, pair in frozen_values.items():
            if not np.isclose(pair[0], pair[1], rtol=0.0, atol=1e-15):
                raise ValueError("Frozen setting mismatch for %s: %s != %s" % (label, pair[0], pair[1]))
        if sorted(dt_values) != [0.0125, 0.025, 0.05]:
            raise ValueError("Server dt grid must be [0.05, 0.025, 0.0125] ms")
        expected_seeds = [20260818, 21260821, 22260823, 23260837, 24260855]
        if seeds != expected_seeds:
            raise ValueError("Frozen multiseed sequence mismatch")
    return {"status": "PASS", "frozen_dir": str(root), "selected_scenarios": len(selected),
            "selected_pairs": pairs, "dt_ms": dt_values, "seeds": seeds,
            "n_p": int(cfg["path"]["n_p"]), "workers": int(cfg["parallel"]["workers"])}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Full-ensemble convergence test for WT--SCA3 KL-balance ordering")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("validate", "run"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--config", required=True)
        sub.add_argument("--frozen-dir", required=True)
    fig = subparsers.add_parser("figures")
    fig.add_argument("--results-dir", required=True)
    args = parser.parse_args(argv)
    if args.command == "figures":
        paths = make_figures(Path(args.results_dir))
        print(json.dumps([str(path) for path in paths], indent=2))
        return
    cfg = load_config(args.config)
    checks = validate(cfg, args.frozen_dir)
    if args.command == "validate":
        print(json.dumps(checks, indent=2))
        return
    output, summary = run(cfg, args.frozen_dir)
    make_figures(output)
    print(json.dumps({
        "version": summary["version"], "status": summary["status"],
        "selected_scenarios": summary["selected_scenarios"],
        "dependent_pairs": summary["dependent_pairs"],
        "stationary_state_simulations": summary["stationary_state_simulations"],
        "decision": summary["verdict"]["decision"],
        "failed_gate_count": summary["verdict"]["n_failed_gates"],
    }, indent=2))
    print("Results written to %s" % output)
