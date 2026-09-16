from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .vp import victor_purpura


@dataclass
class LatencyAlignmentResult:
    shift_ms: float
    vp_loss: float
    count_error_fraction: float
    aligned_spikes: np.ndarray
    raw_vp_loss: float
    raw_count_error_fraction: float
    raw_spikes: np.ndarray
    alignment_applied: bool
    count_preserved: bool


def align_first_spike(exp_spikes, model_spikes, fit_end_ms: float, vp_tau_ms: float,
                      normalize: bool, cfg=None, stage: str = 'fine') -> LatencyAlignmentResult:
    """Exactly align the first model spike to the first experimental spike.

    v3.6 rules:
      * the HR simulation itself is never shifted;
      * raw model spikes are defined in the original, unshifted fit window;
      * tau = t_exp_first - t_model_first whenever both trains are non-empty;
      * the same tau is added to every raw model spike;
      * aligned spikes are NOT clipped back to [0, fit_end_ms];
      * spike count and every model ISI are therefore exactly preserved;
      * no extra latency loss and no optimization over tau are used.

    ``stage`` is accepted only to keep the objective/optimizer interface identical
    across search, refinement, fine evaluation, and identifiability.
    """
    del stage
    exp = np.sort(np.asarray(exp_spikes, dtype=float))
    model = np.sort(np.asarray(model_spikes, dtype=float))
    fit_end = float(fit_end_ms)

    # The raw train is frozen BEFORE alignment. This is the key v3.6 invariant.
    exp = np.asarray(exp[(exp >= -1e-9) & (exp <= fit_end + 1e-9)], dtype=float)
    raw = np.asarray(model[(model >= -1e-9) & (model <= fit_end + 1e-9)], dtype=float)

    raw_vp = float(victor_purpura(exp, raw, tau_ms=float(vp_tau_ms), normalize=bool(normalize)))
    raw_count = abs(len(raw) - len(exp)) / max(len(exp), 1)

    enabled = True
    if cfg is not None:
        enabled = bool(cfg.get('loss', {}).get('latency_alignment', {}).get('enabled', True))

    if (not enabled) or exp.size == 0 or raw.size == 0:
        return LatencyAlignmentResult(
            shift_ms=0.0,
            vp_loss=raw_vp,
            count_error_fraction=float(raw_count),
            aligned_spikes=raw.copy(),
            raw_vp_loss=raw_vp,
            raw_count_error_fraction=float(raw_count),
            raw_spikes=raw.copy(),
            alignment_applied=False,
            count_preserved=True,
        )

    shift = float(exp[0] - raw[0])
    # Deliberately no clipping after the shift. Spikes may lie beyond fit_end_ms.
    aligned = np.asarray(raw + shift, dtype=float)
    vp = float(victor_purpura(exp, aligned, tau_ms=float(vp_tau_ms), normalize=bool(normalize)))
    count_frac = abs(len(raw) - len(exp)) / max(len(exp), 1)
    preserved = len(aligned) == len(raw)
    if not preserved:
        raise RuntimeError('v3.6 latency alignment changed model spike count')

    return LatencyAlignmentResult(
        shift_ms=shift,
        vp_loss=vp,
        count_error_fraction=float(count_frac),
        aligned_spikes=aligned,
        raw_vp_loss=raw_vp,
        raw_count_error_fraction=float(raw_count),
        raw_spikes=raw,
        alignment_applied=True,
        count_preserved=True,
    )
