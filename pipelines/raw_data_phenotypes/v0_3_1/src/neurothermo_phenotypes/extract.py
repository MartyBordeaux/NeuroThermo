from __future__ import annotations

import hashlib
from typing import Optional, Tuple

import numpy as np
from scipy.signal import find_peaks

from .io import Trace
from .metrics import (
    electrical_work,
    irreversibility_with_surrogates,
    normalized_spectral_entropy,
    ordinal_predictive_information,
    permutation_entropy,
    uniform_resample,
)


def _longest_true_run(mask: np.ndarray) -> Optional[Tuple[int, int]]:
    padded = np.r_[False, mask.astype(bool), False]
    changes = np.diff(padded.astype(int))
    starts, ends = np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)
    if len(starts) == 0:
        return None
    lengths = ends - starts
    k = int(np.argmax(lengths))
    return int(starts[k]), int(ends[k] - 1)


def stimulus_window(trace: Trace, cfg: dict) -> tuple[float, float, str]:
    if np.isfinite(trace.stim_start_override_s) and np.isfinite(trace.stim_end_override_s):
        return trace.stim_start_override_s, trace.stim_end_override_s, "override"
    stim_cfg = cfg["stimulus"]
    command = trace.current_trace_pA
    if stim_cfg["detection"] in {"command", "command_or_fixed"} and command is not None:
        n0 = max(10, len(command) // 20)
        baseline = float(np.nanmedian(command[:n0]))
        threshold = max(float(stim_cfg["command_threshold_pA"]), 8.0 * float(np.nanmedian(np.abs(command[:n0] - baseline))))
        run = _longest_true_run(np.abs(command - baseline) > threshold)
        if run is not None:
            start, end = float(trace.time_s[run[0]]), float(trace.time_s[run[1]])
            if end - start >= float(stim_cfg["minimum_duration_s"]):
                return start, end, "command"
    return float(stim_cfg["fixed_start_s"]), float(stim_cfg["fixed_end_s"]), "fixed"


def detect_spikes(time_s: np.ndarray, voltage_mV: np.ndarray, cfg: dict) -> np.ndarray:
    if len(time_s) < 3:
        return np.array([], dtype=float)
    dt_ms = float(np.median(np.diff(time_s))) * 1000.0
    distance = max(1, int(round(float(cfg["refractory_ms"]) / dt_ms)))
    peaks, _ = find_peaks(
        voltage_mV,
        height=float(cfg["peak_threshold_mV"]),
        prominence=float(cfg["peak_prominence_mV"]),
        distance=distance,
    )
    return np.asarray(time_s[peaks], float)


def _adaptation(isi_ms: np.ndarray) -> tuple[float, float]:
    if len(isi_ms) < 4:
        return np.nan, np.nan
    half = len(isi_ms) // 2
    early = float(np.median(isi_ms[:half]))
    late = float(np.median(isi_ms[-half:]))
    ratio = late / early if early > 0 else np.nan
    x = np.arange(len(isi_ms), dtype=float)
    slope = float(np.polyfit(x, isi_ms, 1)[0])
    return ratio, slope


def _seed(base_seed: int, trace: Trace) -> int:
    text = f"{base_seed}|{trace.group}|{trace.cell_id}|{trace.sweep_index}".encode()
    return int.from_bytes(hashlib.sha256(text).digest()[:8], "little") % (2**32 - 1)


def analyse_trace(trace: Trace, cfg: dict) -> dict:
    t, v = np.asarray(trace.time_s, float), np.asarray(trace.voltage_mV, float)
    order = np.argsort(t)
    t, v = t[order], v[order]
    command = None if trace.current_trace_pA is None else np.asarray(trace.current_trace_pA, float)[order]
    finite = np.isfinite(t) & np.isfinite(v)
    finite_fraction = float(np.mean(finite)) if len(finite) else 0.0
    t, v = t[finite], v[finite]
    command = None if command is None else command[finite]
    start, end, stim_source = stimulus_window(trace, cfg)
    qc_cfg, spike_cfg, thermo_cfg = cfg["qc"], cfg["spikes"], cfg["thermodynamics"]
    guard = float(qc_cfg["baseline_guard_ms"]) / 1000.0
    baseline_duration = float(qc_cfg["baseline_duration_ms"]) / 1000.0
    baseline_mask = (t >= max(t[0], start - guard - baseline_duration)) & (t < start - guard)
    stim_mask = (t >= start) & (t <= end)
    flags, fatal = [], False
    if baseline_mask.sum() < 20:
        flags.append("insufficient_baseline")
        fatal = True
    if stim_mask.sum() < 20:
        flags.append("insufficient_stimulus")
        fatal = True
    baseline_v = float(np.nanmedian(v[baseline_mask])) if baseline_mask.any() else np.nan
    baseline_noise = float(1.4826 * np.nanmedian(np.abs(v[baseline_mask] - baseline_v))) if baseline_mask.any() else np.nan
    sample_rate = 1.0 / float(np.median(np.diff(t))) if len(t) > 2 else np.nan
    if finite_fraction < float(qc_cfg["minimum_finite_fraction"]):
        flags.append("nonfinite_samples")
        fatal = True
    if not np.isfinite(sample_rate) or sample_rate < float(qc_cfg["minimum_sampling_hz"]) * (1.0 - 1e-6):
        flags.append("low_sampling_rate")
        fatal = True
    if np.isfinite(baseline_v) and not (float(qc_cfg["baseline_min_mV"]) <= baseline_v <= float(qc_cfg["baseline_max_mV"])):
        flags.append("baseline_out_of_range")
    if np.isfinite(baseline_noise) and baseline_noise > float(qc_cfg["baseline_noise_warn_mV"]):
        flags.append("high_baseline_noise")
    if stim_mask.any():
        stim_v = v[stim_mask]
        clip_fraction = max(float(np.mean(stim_v == np.nanmin(stim_v))), float(np.mean(stim_v == np.nanmax(stim_v))))
    else:
        clip_fraction = np.nan
    if np.isfinite(clip_fraction) and clip_fraction > float(qc_cfg["clipping_fraction_warn"]):
        flags.append("possible_clipping")
    if fatal:
        row = _base_row(trace, start, end, stim_source, finite_fraction, sample_rate, baseline_v, baseline_noise, flags, False)
        row.update(_empty_metrics("fatal_qc"))
        return row

    ts, vs = t[stim_mask], v[stim_mask]
    cmd_stim = None if command is None else command[stim_mask]
    curated = trace.metadata.get("curated_spikes_s")
    if curated is None:
        spikes = detect_spikes(ts, vs, spike_cfg)
        spike_source = "automatic_voltage_peaks"
        event_tolerance_s = 0.0
    else:
        event_tolerance_s = float(
            cfg["input"].get("curated_event_boundary_tolerance_ms", 0.0)
        ) / 1000.0
        loaded_curated_spikes = np.asarray(curated, float)
        spikes = loaded_curated_spikes
        spikes = spikes[
            (spikes >= start - event_tolerance_s)
            & (spikes <= end + event_tolerance_s)
        ]
        spike_source = "curated_events"
    threshold_probe = trace.metadata.get("threshold_probe_spikes_s")
    if threshold_probe is None:
        threshold_spikes = spikes
        threshold_probe_loaded = len(spikes)
    else:
        loaded_threshold_spikes = np.asarray(threshold_probe, float)
        threshold_probe_loaded = len(loaded_threshold_spikes)
        threshold_spikes = loaded_threshold_spikes[
            (loaded_threshold_spikes >= start - event_tolerance_s)
            & (loaded_threshold_spikes <= end + event_tolerance_s)
        ]
    steady_start = start + float(spike_cfg["steady_discard_ms"]) / 1000.0
    steady_spikes = spikes[spikes >= steady_start]
    isi_ms = np.diff(spikes) * 1000.0
    steady_isi_ms = np.diff(steady_spikes) * 1000.0
    duration = end - start
    steady_duration = end - steady_start
    adaptation_ratio, adaptation_slope = _adaptation(isi_ms)
    work = electrical_work(ts, vs, trace.current_pA, baseline_v, cmd_stim)
    n_spikes = int(len(spikes))
    row = _base_row(trace, start, end, stim_source, finite_fraction, sample_rate, baseline_v, baseline_noise, flags, True)
    row.update({
        "n_spikes": n_spikes,
        "n_sustained_spikes": int(len(steady_spikes)),
        "is_spiking": n_spikes > 0,
        "firing_rate_hz": n_spikes / duration if duration > 0 else np.nan,
        "sustained_rate_hz": len(steady_spikes) / steady_duration if steady_duration > 0 else np.nan,
        "first_spike_latency_ms": (spikes[0] - start) * 1000.0 if n_spikes else np.nan,
        "mean_isi_ms": float(np.mean(isi_ms)) if len(isi_ms) else np.nan,
        "cv_isi": float(np.std(isi_ms, ddof=1) / np.mean(isi_ms)) if len(isi_ms) >= 2 and np.mean(isi_ms) > 0 else np.nan,
        "steady_mean_isi_ms": float(np.mean(steady_isi_ms)) if len(steady_isi_ms) else np.nan,
        "adaptation_ratio": adaptation_ratio,
        "adaptation_slope_ms_per_interval": adaptation_slope,
        "median_voltage_stim_mV": float(np.median(vs)),
        "depolarization_area_mV_s": float(np.trapz(vs - baseline_v, ts)),
        "spike_source": spike_source,
        "curated_peak_override": trace.metadata.get("curated_peak_override", ""),
        "curated_frozen_sweep": bool(trace.metadata.get("frozen_sweep_membership", False)),
        "curated_events_loaded": int(trace.metadata.get("curated_events_loaded", len(spikes))) if curated is not None else 0,
        "curated_events_used": int(len(spikes)) if curated is not None else 0,
        "curated_events_outside_stimulus_window": int(
            trace.metadata.get("curated_events_loaded", len(spikes)) - len(spikes)
        ) if curated is not None else 0,
        "threshold_probe_events_loaded": int(threshold_probe_loaded),
        "threshold_probe_event_count": int(len(threshold_spikes)),
        "threshold_probe_spiking": bool(len(threshold_spikes) > 0),
        **work,
    })
    row["work_per_spike_fJ"] = work["external_work_signed_fJ"] / n_spikes if n_spikes else np.nan

    stationary_start = start + float(thermo_cfg["stationary_discard_ms"]) / 1000.0
    stationary_mask = (t >= stationary_start) & (t <= end)
    thermo_eligible = n_spikes >= int(spike_cfg["thermo_min_spikes"]) and stationary_mask.sum() >= int(thermo_cfg["minimum_stationary_samples"])
    row["thermo_eligible"] = bool(thermo_eligible)
    row["thermo_exclusion_reason"] = "" if thermo_eligible else ("nonspiking_or_too_few_spikes" if n_spikes < int(spike_cfg["thermo_min_spikes"]) else "short_stationary_window")
    if stationary_mask.sum() >= int(thermo_cfg["minimum_stationary_samples"]):
        target_dt = float(thermo_cfg["resample_dt_ms"]) / 1000.0
        rt, rv = uniform_resample(t[stationary_mask], v[stationary_mask], target_dt)
        rv = rv - np.mean(rv)
        perm_delay = max(1, int(round((float(thermo_cfg["permutation_delay_ms"]) / 1000.0) / target_dt)))
        row.update({
            "permutation_entropy_norm": permutation_entropy(rv, int(thermo_cfg["permutation_order"]), perm_delay),
            "spectral_entropy_norm": normalized_spectral_entropy(rv, 1.0 / target_dt),
            "predictive_information_nats": ordinal_predictive_information(rv, int(thermo_cfg["permutation_order"]), perm_delay),
            "stationary_samples": int(len(rv)),
        })
        rng = np.random.default_rng(_seed(int(cfg["seed"]), trace))
        row.update(irreversibility_with_surrogates(rv, target_dt, thermo_cfg, rng))
    else:
        for name in [
            "permutation_entropy_norm", "spectral_entropy_norm", "predictive_information_nats",
            "path_kl_rate_raw_nats_s", "path_kl_surrogate_median_nats_s",
            "path_kl_rate_bias_corrected_nats_s", "path_kl_rate_excess_nats_s",
            "path_kl_surrogate_p", "path_word_coverage", "path_n_words",
        ]:
            row[name] = np.nan
        row["stationary_samples"] = int(stationary_mask.sum())
    return row


def _base_row(trace, start, end, source, finite_fraction, sample_rate, baseline_v, baseline_noise, flags, qc_pass):
    def density(cap):
        return trace.current_pA / cap if np.isfinite(cap) and cap > 0 else np.nan
    return {
        "group": trace.group, "cell_id": trace.cell_id, "animal_id": trace.animal_id,
        "record_id": trace.record_id, "sweep_index": trace.sweep_index,
        "source_path": trace.source_path, "current_pA": trace.current_pA,
        "capacitance_pF": trace.capacitance_pF,
        "capacitance_10ms_pF": trace.capacitance_10ms_pF,
        "capacitance_20ms_pF": trace.capacitance_20ms_pF,
        "capacitance_50ms_pF": trace.capacitance_50ms_pF,
        "J_pA_per_pF": density(trace.capacitance_pF),
        "J_10ms_pA_per_pF": density(trace.capacitance_10ms_pF),
        "J_20ms_pA_per_pF": density(trace.capacitance_20ms_pF),
        "J_50ms_pA_per_pF": density(trace.capacitance_50ms_pF),
        "stim_start_s": start, "stim_end_s": end, "stimulus_window_source": source,
        "finite_fraction": finite_fraction, "sampling_rate_hz": sample_rate,
        "baseline_voltage_mV": baseline_v, "baseline_noise_mV": baseline_noise,
        "qc_flags": ";".join(flags),
        "qc_pass": bool(qc_pass),
        "qc_fatal": not bool(qc_pass),
        "qc_warning": bool(qc_pass) and bool(flags),
        "qc_status": "FAIL" if not qc_pass else ("WARN" if flags else "PASS"),
    }


def _empty_metrics(reason: str) -> dict:
    row = {
        "n_spikes": 0, "n_sustained_spikes": 0, "is_spiking": False,
        "thermo_eligible": False, "thermo_exclusion_reason": reason,
        "stationary_samples": 0, "work_current_source": "not_computed",
        "spike_source": "not_computed", "curated_peak_override": "",
        "curated_frozen_sweep": False,
        "curated_events_loaded": 0, "curated_events_used": 0,
        "curated_events_outside_stimulus_window": 0,
        "threshold_probe_events_loaded": 0,
        "threshold_probe_event_count": 0, "threshold_probe_spiking": False,
    }
    for name in [
        "firing_rate_hz", "sustained_rate_hz", "first_spike_latency_ms",
        "mean_isi_ms", "cv_isi", "steady_mean_isi_ms", "adaptation_ratio",
        "adaptation_slope_ms_per_interval", "median_voltage_stim_mV",
        "depolarization_area_mV_s", "external_work_signed_fJ",
        "external_work_positive_fJ", "mean_power_signed_fW", "work_per_spike_fJ",
        "permutation_entropy_norm", "spectral_entropy_norm", "predictive_information_nats",
        "path_kl_rate_raw_nats_s", "path_kl_surrogate_median_nats_s",
        "path_kl_rate_bias_corrected_nats_s", "path_kl_rate_excess_nats_s",
        "path_kl_surrogate_p", "path_word_coverage", "path_n_words",
    ]:
        row[name] = np.nan
    return row
