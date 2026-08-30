from __future__ import annotations

import argparse
from pathlib import Path
import json

import yaml
import numpy as np

from .data import load_frozen, resolve_frozen_dir, select_scenarios
from .pipeline import run
from .verdict import validate_physical_mapping
from .aggregation import annotate_animal_pairs, load_animal_mapping


def load_config(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    return cfg


def validate(cfg, frozen_dir):
    validate_physical_mapping(cfg)
    root = resolve_frozen_dir(frozen_dir)
    scenarios, pair_stage, _, _ = load_frozen(root)
    selected = select_scenarios(scenarios, pair_stage, cfg)
    validation_cfg = cfg.get("validation", {})
    expected_scenarios = validation_cfg.get("expected_scenarios")
    expected_pairs = validation_cfg.get("expected_pairs")
    if expected_scenarios is not None and len(selected) != int(expected_scenarios):
        raise ValueError(f"Selected {len(selected)} scenarios; expected {int(expected_scenarios)}")
    observed_pairs = selected["biological_pair_key"].nunique()
    if expected_pairs is not None and observed_pairs != int(expected_pairs):
        raise ValueError(f"Selected {observed_pairs} biological pairs; expected {int(expected_pairs)}")
    mapping = load_animal_mapping(cfg)
    animal_counts = {"wt_animals": 0, "sca3_animals": 0, "animal_pairs": 0}
    if mapping is not None:
        annotated = annotate_animal_pairs(
            selected[["biological_pair_key"]].drop_duplicates(),
            mapping,
            required=bool(cfg.get("animal_mapping", {}).get("required", False)),
        )
        animal_counts = {
            "wt_animals": int(annotated["wt_animal_id"].nunique()),
            "sca3_animals": int(annotated["sca3_animal_id"].nunique()),
            "animal_pairs": int(annotated["animal_pair_key"].nunique()),
        }
        for label, key in (
            ("WT animals", "expected_wt_animals"),
            ("SCA3 animals", "expected_sca3_animals"),
            ("animal pairs", "expected_animal_pairs"),
        ):
            expected = validation_cfg.get(key)
            observed = animal_counts[key.replace("expected_", "")]
            if expected is not None and observed != int(expected):
                raise ValueError(f"Selected {observed} {label}; expected {int(expected)}")
    checks = {
        "status": "PASS", "frozen_dir": str(root), "selected_scenarios": len(selected), "selected_pairs": int(observed_pairs),
        "seeds": [int(value) for value in cfg["multiseed"]["seeds"]],
        "n_p": int(cfg["path"]["n_p"]), "workers": int(cfg["parallel"]["workers"]),
        **animal_counts,
    }
    if int(cfg["path"]["n_p"]) < 3:
        raise ValueError("path.n_p must be at least 3 for a central Fisher derivative.")
    protocol_points = int(cfg["protocol"]["n_points"])
    if not 2 <= protocol_points <= int(cfg["path"]["n_p"]):
        raise ValueError("protocol.n_points must be between 2 and path.n_p")
    if any(float(value) <= 0 for value in cfg["noise"]["D"]):
        raise ValueError("All diffusion coefficients must be strictly positive.")
    if int(cfg["density"]["bins"]) < 6:
        raise ValueError("density.bins must be at least 6.")
    if cfg.get("analysis", {}).get("name") == "nonequilibrium_geometry":
        frozen_values = {
            "path.n_p": (int(cfg["path"]["n_p"]), 31),
            "stationary.dt_ms": (float(cfg["stationary"]["dt_ms"]), 0.025),
            "stationary.burn_ms": (float(cfg["stationary"]["burn_ms"]), 2400.0),
            "stationary.sample_ms": (float(cfg["stationary"]["sample_ms"]), 6000.0),
            "stationary.sample_stride_ms": (float(cfg["stationary"]["sample_stride_ms"]), 0.5),
            "density.bins": (int(cfg["density"]["bins"]), 22),
        }
        for label, (observed, expected) in frozen_values.items():
            if not np.isclose(observed, expected, rtol=0, atol=1e-15):
                raise ValueError(f"Frozen setting mismatch for {label}: {observed} != {expected}")
        expected_seeds = [20260818, 21260821, 22260823, 23260837, 24260855]
        if [int(value) for value in cfg["multiseed"]["seeds"]] != expected_seeds:
            raise ValueError("Frozen multiseed sequence mismatch.")
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser(description="Nonequilibrium geometry of the frozen WT-to-SCA3 HR path")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("validate", "run"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--config", required=True)
        sub.add_argument("--frozen-dir", required=True, help="Exact directory containing the four frozen v1.0.1 CSV files")
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    checks = validate(cfg, args.frozen_dir)
    if args.command == "validate":
        print(json.dumps(checks, indent=2))
        return
    output, summary = run(cfg, args.frozen_dir)
    print(json.dumps(summary, indent=2))
    print(f"Results written to {output}")
