from __future__ import annotations
from copy import deepcopy
from pathlib import Path
from typing import Any, Union
import yaml

DEFAULTS: dict[str, Any] = {
    'data': {
        'events_file': None,
        'events_include_column': 'fixed_qc_detected',
        'frozen_sweeps_manifest': '../calibration/frozen_accepted_spiking_sweeps_v3_5.csv',
        'peak_overrides_file': '../calibration/frozen_peak_overrides_v3_5.csv',
        'baseline_cell_summary_file': '../calibration/frozen_v3_1_cell_fit_summary.csv',
        'baseline_sweep_summary_file': '../calibration/frozen_v3_1_sweep_fit_summary.csv',
        'baseline_identifiability_file': '../calibration/frozen_v3_1_identifiability.csv',
        'seed_cell_summary_file': '../calibration/seed_cell_summary_v3_9.csv',
        'threshold_brackets_file': '../calibration/frozen_threshold_brackets_v3_5.csv',
        'strict_metadata_match': True,
        'metadata_tolerance_ms': 0.05,
    },
    'model': {
        'a': 1.0, 'c': 1.0, 'd': 5.0, 'x_R': -1.6,
        'x0': -1.6, 'y0': -10.0, 'z0': 1.0,
        'pre_ms': 500.0,
        'model_time_scale_ms': 1.0,
        'model_spike_threshold': 0.0,
        'model_refractory_ms': 1.0,
    },
    'bounds': {
        'b': {'min': 0.5, 'max': 7.0, 'scale': 'linear'},
        'r': {'min': 0.0001, 'max': 0.10, 'scale': 'log'},
        's': {'min': 0.05, 'max': 15.0, 'scale': 'linear'},
        'kappa_I': {'min': 0.0002, 'max': 2.0, 'scale': 'log'},
    },
    'loss': {
        'vp_tau_ms': 10.0,
        'normalize': True,
        'simulation_failure_loss': 1.0e6,
        'count_penalty_weight': 0.25,
        'sweep_weighting': 'equal',
        # No explicit ISI or adaptation term. Timing is constrained through VP on selected spike times.
        'fit_window': {
            'local_isi_count': 5,
            'post_last_spike_guard_isi_multiplier': 1.0,
            'min_guard_ms': 5.0,
            'max_guard_ms': 30.0,
            'single_spike_guard_ms': 15.0,
            'max_post_stimulus_ms': 15.0,
        },
        'latency_alignment': {
            'enabled': True,
            'method': 'exact_first_spike',
            'large_shift_warning_ms': 200.0,
        },
    },
    'threshold_constraint': {
        'enabled': True,
        # Binary rheobase bracket only; voltage plateaus and threshold-sweep spike counts are not fitted.
        'nonspiking_violation_penalty': 1.0,
        'first_spiking_violation_penalty': 1.0,
        'require_pass_for_auto_accept': True,
    },
    'optimization': {
        'seed': 20260815,
        'dt_search_ms': 0.15,
        'dt_refine_ms': 0.10,
        'dt_fine_ms': 0.05,
        'search_vp_tau_ms': 30.0,
        'prior_radius_fraction': 0.40,
        'prior_de_popsize': 14,
        'prior_de_maxiter': 85,
        'n_prior_starts': 1,
        'global_rescue_enabled': True,
        'global_search_always': True,
        'global_rescue_loss_threshold': 1.0,
        'n_global_starts': 3,
        'global_de_popsize': 18,
        'global_de_maxiter': 150,
        'refine_radius_fraction': 0.12,
        'refine_de_popsize': 18,
        'refine_de_maxiter': 130,
        'de_tol': 0.001,
        'n_jobs': 2,
    },
    'sweep_qc': {
        'circle_match_tolerance_ms': 10.0,
        'leading_extra_tolerance_ms': 5.0,
        'accept_max_count_error_fraction': 0.15,
        'accept_min_circle_f1': 0.80,
        'accept_max_leading_extra_spikes': 3,
        'bad_count_error_fraction': 0.50,
        'bad_min_circle_f1': 0.50,
    },
    'cell_qc': {
        'accept_min_fraction_sweeps_accept': 0.60,
        'accept_min_fraction_sweeps_nonbad': 0.80,
        'accept_min_median_circle_f1': 0.75,
        'bad_min_fraction_sweeps_bad': 0.50,
        'bad_max_median_circle_f1': 0.50,
        'primary_min_sweeps': 2,
    },
    'identifiability': {
        'enabled': True,
        'min_spiking_sweeps': 2,
        'reference_separation_fraction': 0.15,
        'reference_bounds': {
            'b': {'min': 1.0, 'max': 5.0, 'scale': 'linear'},
            'r': {'min': 0.0005, 'max': 0.05, 'scale': 'log'},
            's': {'min': 0.5, 'max': 10.0, 'scale': 'linear'},
            'kappa_I': {'min': 0.001, 'max': 1.0, 'scale': 'log'},
        },
        'near_optimal_absolute_loss': 0.05,
        'near_optimal_relative_loss': 0.25,
        'dt_ms': 0.10,
        'de_popsize': 10,
        'de_maxiter': 65,
        'n_jobs': 2,
    },
    'output': {
        'dir': '../results_cellfit_v3_9',
        'resume': True,
        'make_plots': True,
        'make_audit_pdf': True,
    },
}


def _merge(base, override):
    out = deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def _resolve(base: Path, value):
    if not value:
        return value
    p = Path(value).expanduser()
    return str(p if p.is_absolute() else (base / p).resolve())


def load_config(path: Union[str, Path]):
    path = Path(path).expanduser().resolve()
    with path.open('r', encoding='utf-8') as f:
        user = yaml.safe_load(f) or {}
    cfg = _merge(DEFAULTS, user)
    base = path.parent
    for key in (
        'events_file', 'frozen_sweeps_manifest', 'peak_overrides_file',
        'baseline_cell_summary_file', 'baseline_sweep_summary_file',
        'baseline_identifiability_file', 'seed_cell_summary_file', 'threshold_brackets_file',
    ):
        cfg['data'][key] = _resolve(base, cfg['data'].get(key))
    cfg['output']['dir'] = _resolve(base, cfg['output'].get('dir'))
    return cfg, path
