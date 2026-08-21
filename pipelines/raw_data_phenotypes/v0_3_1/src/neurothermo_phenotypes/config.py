from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Optional, Union

import yaml


DEFAULT_CONFIG = {
    "seed": 20260814,
    "input": {
        "data_root": ".",
        "voltage_channel": 0,
        "current_source": "configured",
        "protocol_currents_pA": [float(x) for x in range(0, 601, 50)],
        "analysis_currents_pA": [float(x) for x in range(0, 601, 50)],
        "enforce_common_current_grid": False,
        "require_complete_current_grid": False,
        "current_tolerance_pA": 1e-6,
        "sweep_overrides_csv": None,
        "curated_events_csv": None,
        "require_curated_events": False,
        "events_include_column": "fixed_qc_detected",
        "curated_sweeps_manifest_csv": None,
        "restrict_curated_events_to_frozen_sweeps": False,
        "curated_peak_overrides_csv": None,
        "curated_threshold_brackets_csv": None,
        "curated_hash_manifest_json": None,
        "strict_curated_metadata": True,
        "curated_metadata_tolerance_ms": 0.05,
        "curated_event_boundary_tolerance_ms": 0.0,
        "expected_curated_sweeps": None,
        "expected_curated_spike_events": None,
        "expected_curated_events_used": None,
    },
    "cohort": {},
    "stimulus": {
        "detection": "command_or_fixed",
        "fixed_start_s": 0.1,
        "fixed_end_s": 1.1,
        "command_threshold_pA": 5.0,
        "minimum_duration_s": 0.5,
    },
    "spikes": {
        "peak_threshold_mV": -20.0,
        "peak_prominence_mV": 15.0,
        "refractory_ms": 1.5,
        "steady_discard_ms": 10.0,
        "thermo_min_spikes": 3,
    },
    "qc": {
        "baseline_guard_ms": 5.0,
        "baseline_duration_ms": 75.0,
        "baseline_min_mV": -100.0,
        "baseline_max_mV": -30.0,
        "baseline_noise_warn_mV": 5.0,
        "minimum_sampling_hz": 2000.0,
        "minimum_finite_fraction": 0.999,
        "clipping_fraction_warn": 0.005,
    },
    "thermodynamics": {
        "stationary_discard_ms": 100.0,
        "resample_dt_ms": 1.0,
        "permutation_order": 4,
        "permutation_delay_ms": 1.0,
        "symbol_bins": 6,
        "word_length": 3,
        "word_delay_ms": 2.0,
        "pseudocount": 0.5,
        "n_reversible_surrogates": 100,
        "minimum_stationary_samples": 300,
    },
    "statistics": {
        "minimum_cells_per_group": 3,
        "bootstrap_iterations": 2000,
        "permutation_iterations": 10000,
        "current_round_decimals": 6,
        "J_grid_step_pA_per_pF": 0.1,
        "primary_inference_features": None,
        "diagnostic_inference_features": [],
        "exact_max_labelings": 200000,
        "two_part_currents_pA": None,
        "two_part_primary_features": [
            "mean_isi_ms", "predictive_information_nats",
        ],
        "two_part_secondary_features": [
            "work_per_spike_fJ", "path_kl_rate_excess_nats_s",
        ],
    },
    "disease_coordinate": {
        "enabled": True,
        "capacitance_feature": "capacitance_20ms_pF",
        "capacitance_sensitivity_features": [
            "capacitance_10ms_pF", "capacitance_20ms_pF", "capacitance_50ms_pF",
        ],
        "currents_pA": [400.0, 450.0, 500.0, 550.0, 600.0],
        "z_clip": 5.0,
        "minimum_scale": 1e-9,
        "robust_z_threshold": 2.5,
        "consensus_min_domains": 2,
        "stability_bootstrap_iterations": 2000,
        "exact_max_labelings": 200000,
        "permutation_iterations": 10000,
        "permutation_chunk_size": 256,
    },
    "runtime": {"workers": 2},
}


def _deep_update(base: dict, update: dict) -> dict:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
    return base


def load_config(path: Optional[Union[str, Path]]) -> dict:
    cfg = deepcopy(DEFAULT_CONFIG)
    if path is None:
        return cfg
    config_path = Path(path).expanduser()
    supplied = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    resolved = _deep_update(cfg, supplied)
    for key in [
        "curated_sweeps_manifest_csv", "curated_peak_overrides_csv",
        "curated_threshold_brackets_csv", "curated_hash_manifest_json",
    ]:
        value = resolved["input"].get(key)
        if value:
            value_path = Path(value).expanduser()
            if not value_path.is_absolute():
                value_path = (config_path.resolve().parent / value_path).resolve()
            resolved["input"][key] = str(value_path)
    return resolved
