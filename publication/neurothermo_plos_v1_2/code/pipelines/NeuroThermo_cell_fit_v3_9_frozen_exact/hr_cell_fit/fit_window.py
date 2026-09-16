from __future__ import annotations
import numpy as np


def compute_fit_window(exp_spikes_ms, stimulus_duration_ms: float, cfg) -> dict:
    """Define the VP observation window for v2.1.

    The window always starts at current-step onset (t=0) so latency is identifiable.
    It ends shortly after the last reviewed experimental spike. Any later quiescent
    plateau is not part of the loss. The current protocol itself is unchanged:
    J is applied until the experimental step offset and is zero afterwards.
    """
    s = np.sort(np.asarray(exp_spikes_ms, dtype=float))
    if s.size == 0:
        raise ValueError("compute_fit_window requires at least one experimental spike")
    if np.any(~np.isfinite(s)):
        raise ValueError("Experimental spike times must be finite")
    if np.any(s < 0.0):
        raise ValueError("Experimental spike times before stimulus onset are not supported")

    wcfg = cfg["loss"].get("fit_window", {})
    mode = str(wcfg.get("mode", "spiking_region")).lower()
    if mode != "spiking_region":
        raise ValueError(f"Unsupported loss.fit_window.mode={mode!r}; v2.1 requires 'spiking_region'")

    local_n = max(1, int(wcfg.get("local_isi_count", 5)))
    multiplier = float(wcfg.get("post_last_spike_guard_isi_multiplier", 1.0))
    min_guard = float(wcfg.get("min_guard_ms", 5.0))
    max_guard = float(wcfg.get("max_guard_ms", 30.0))
    single_guard = float(wcfg.get("single_spike_guard_ms", 15.0))
    max_post_stim = float(wcfg.get("max_post_stimulus_ms", 15.0))

    if min_guard < 0 or max_guard < min_guard or max_post_stim < 0:
        raise ValueError("Invalid fit-window guard configuration")

    local_isi = np.nan
    if s.size >= 2:
        isi = np.diff(s)
        local = isi[-min(local_n, isi.size):]
        local = local[np.isfinite(local) & (local > 0)]
        if local.size:
            local_isi = float(np.median(local))
            raw_guard = multiplier * local_isi
        else:
            raw_guard = single_guard
    else:
        raw_guard = single_guard

    guard = float(np.clip(raw_guard, min_guard, max_guard))
    last_spike = float(s[-1])
    protocol_end = float(stimulus_duration_ms) + max_post_stim
    fit_end = min(last_spike + guard, protocol_end)
    # Always retain the last selected event numerically, even under extreme config.
    fit_end = max(fit_end, last_spike + 1e-9)

    excluded_plateau = max(0.0, float(stimulus_duration_ms) - fit_end)
    return {
        "fit_start_ms": 0.0,
        "fit_end_ms": float(fit_end),
        "last_exp_spike_ms": last_spike,
        "guard_ms": guard,
        "local_isi_ms": local_isi,
        "stimulus_duration_ms": float(stimulus_duration_ms),
        "excluded_plateau_ms": float(excluded_plateau),
        "current_off_inside_fit_window": bool(fit_end > float(stimulus_duration_ms)),
    }
