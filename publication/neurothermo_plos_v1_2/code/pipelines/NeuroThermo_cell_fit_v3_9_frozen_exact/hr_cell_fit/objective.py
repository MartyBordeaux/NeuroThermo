from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from typing import Optional
from .model import simulate_spikes
from .params import unit_to_theta
from .latency import align_first_spike


@dataclass
class SweepEval:
    sweep_id: str
    vp_loss: float
    raw_vp_loss: float
    count_error_fraction: float
    raw_count_error_fraction: float
    composite_loss: float
    model_spikes: np.ndarray
    raw_model_spikes: np.ndarray
    latency_shift_ms: float
    latency_alignment_applied: bool
    count_preserved_by_alignment: bool
    ok: bool


@dataclass
class ThresholdEval:
    nonspiking_model_spikes: np.ndarray
    first_spiking_model_spikes: np.ndarray
    nonspiking_violation: bool
    first_spiking_violation: bool
    nonspiking_penalty: float
    first_spiking_penalty: float
    total_penalty: float
    pass_constraint: bool
    ok: bool


@dataclass
class CellEval:
    loss: float
    spike_train_loss: float
    threshold_loss: float
    sweep_evals: list[SweepEval]
    threshold_eval: Optional[ThresholdEval]
    ok: bool


def _weights(cell, cfg):
    mode = str(cfg['loss'].get('sweep_weighting', 'equal')).lower()
    n = np.asarray([max(1, len(s['exp_spike_times_ms'])) for s in cell['sweeps']], dtype=float)
    if mode == 'equal':
        w = np.ones_like(n)
    elif mode == 'sqrt_spikes':
        w = np.sqrt(n)
    elif mode == 'spikes':
        w = n
    else:
        raise ValueError('Unknown loss.sweep_weighting=%r' % mode)
    return w / np.mean(w)


def _evaluate_threshold(theta, cell, cfg, dt_ms):
    tc = cfg.get('threshold_constraint', {})
    if not bool(tc.get('enabled', True)):
        return None
    bracket = cell.get('threshold_bracket')
    if not bracket:
        raise ValueError('%s has no threshold_bracket' % cell.get('cell_id', '<cell>'))

    low = bracket['nonspiking_sweep']
    high = bracket['first_spiking_sweep']
    low_end = float(low['stimulus_duration_ms'])
    high_end = float(high['stimulus_duration_ms'])
    low_mod, ok_low = simulate_spikes(theta, low, cfg, float(dt_ms), observation_end_ms=low_end)
    high_mod, ok_high = simulate_spikes(theta, high, cfg, float(dt_ms), observation_end_ms=high_end)
    ok = bool(ok_low and ok_high)
    if not ok:
        failure = float(cfg['loss']['simulation_failure_loss'])
        return ThresholdEval(
            np.asarray(low_mod, dtype=float), np.asarray(high_mod, dtype=float), True, True,
            failure, failure, failure, False, False,
        )

    # Binary rheobase information only. Latency alignment is NEVER applied to threshold probes.
    low_violation = len(low_mod) > 0
    high_violation = len(high_mod) == 0
    low_pen = float(tc.get('nonspiking_violation_penalty', 1.0)) if low_violation else 0.0
    high_pen = float(tc.get('first_spiking_violation_penalty', 1.0)) if high_violation else 0.0
    total = low_pen + high_pen
    return ThresholdEval(
        np.asarray(low_mod, dtype=float), np.asarray(high_mod, dtype=float),
        bool(low_violation), bool(high_violation), float(low_pen), float(high_pen),
        float(total), bool(not low_violation and not high_violation), True,
    )



def evaluate_theta(theta, cell, cfg, dt_ms: float, vp_tau_ms: Optional[float] = None,
                   latency_alignment_stage: str = 'fine') -> CellEval:
    vp_tau = float(cfg['loss']['vp_tau_ms'] if vp_tau_ms is None else vp_tau_ms)
    failure = float(cfg['loss']['simulation_failure_loss'])
    count_weight = float(cfg['loss'].get('count_penalty_weight', 0.0))
    ws = _weights(cell, cfg)
    out = []
    losses = []
    all_ok = True

    for w, sweep in zip(ws, cell['sweeps']):
        fit_end = float(sweep['fit_end_ms'])
        # v3.6 freezes the raw model train in the original fit window BEFORE latency alignment.
        mod_raw, ok = simulate_spikes(theta, sweep, cfg, float(dt_ms), observation_end_ms=fit_end)
        exp = np.asarray(sweep['exp_spike_times_ms'], dtype=float)
        exp = exp[(exp >= 0) & (exp <= fit_end + 1e-9)]
        mod_raw = np.asarray(mod_raw, dtype=float)
        if not ok:
            vp = failure
            raw_vp = failure
            count_frac = 1.0
            raw_count_frac = 1.0
            comp = failure
            aligned = np.asarray([], dtype=float)
            raw = np.asarray([], dtype=float)
            shift = 0.0
            alignment_applied = False
            count_preserved = True
            all_ok = False
        else:
            aln = align_first_spike(
                exp, mod_raw, fit_end, vp_tau, bool(cfg['loss']['normalize']), cfg,
                stage=str(latency_alignment_stage),
            )
            vp = float(aln.vp_loss)
            raw_vp = float(aln.raw_vp_loss)
            # Count penalty is ALWAYS based on the pre-alignment raw train.
            count_frac = float(aln.raw_count_error_fraction)
            raw_count_frac = float(aln.raw_count_error_fraction)
            aligned = np.asarray(aln.aligned_spikes, dtype=float)
            raw = np.asarray(aln.raw_spikes, dtype=float)
            shift = float(aln.shift_ms)
            alignment_applied = bool(aln.alignment_applied)
            count_preserved = bool(aln.count_preserved)
            if len(aligned) != len(raw):
                raise RuntimeError('v3.6 invariant violated: alignment changed spike count')
            comp = vp + count_weight * count_frac
        out.append(SweepEval(
            sweep['sweep_id'], float(vp), float(raw_vp), float(count_frac), float(raw_count_frac),
            float(comp), aligned, raw, float(shift), bool(alignment_applied), bool(count_preserved), bool(ok),
        ))
        losses.append(float(w) * float(comp))

    spike_loss = float(np.mean(losses))
    threshold_eval = _evaluate_threshold(theta, cell, cfg, dt_ms)
    threshold_loss = 0.0 if threshold_eval is None else float(threshold_eval.total_penalty)
    if threshold_eval is not None and not threshold_eval.ok:
        all_ok = False
    total = spike_loss + threshold_loss
    return CellEval(float(total), spike_loss, threshold_loss, out, threshold_eval, bool(all_ok))


def objective_unit(u, cell, cfg, dt_ms: float, vp_tau_ms: float, latency_alignment_stage: str = 'fine') -> float:
    theta = unit_to_theta(u, cfg)
    return evaluate_theta(theta, cell, cfg, dt_ms, vp_tau_ms, latency_alignment_stage).loss
