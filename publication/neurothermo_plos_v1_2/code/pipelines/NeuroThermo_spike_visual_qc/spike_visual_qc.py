#!/usr/bin/env python3
"""Standalone visual audit of NeuroThermo current-clamp spike detection."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences


GROUP_COLORS = {"WT": "#2166ac", "SCA3": "#b2182b"}


def load_config(path: Path) -> dict:
    cfg = json.loads(path.read_text(encoding="utf-8"))
    required = ["voltage_channel", "command_channel", "prominences_mV", "heights_mV",
                "primary_prominence_mV", "primary_height_mV", "minimum_peak_distance_ms"]
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"missing config keys: {missing}")
    return cfg


def longest_true_run(mask: np.ndarray) -> tuple[int, int]:
    idx = np.flatnonzero(mask)
    if not len(idx):
        raise ValueError("no non-baseline command interval")
    cuts = np.flatnonzero(np.diff(idx) > 1)
    starts = np.r_[0, cuts + 1]
    stops = np.r_[cuts + 1, len(idx)]
    k = int(np.argmax(stops - starts))
    run = idx[starts[k]:stops[k]]
    return int(run[0]), int(run[-1] + 1)


def stimulus_from_command(time_ms: np.ndarray, command_pA: np.ndarray) -> tuple[float, float, float]:
    dt = float(np.median(np.diff(time_ms)))
    edge = max(1, min(len(command_pA) // 20, int(round(20.0 / dt))))
    baseline = float(np.median(command_pA[:edge]))
    deviation = np.abs(command_pA - baseline)
    tolerance = max(1e-9, 0.01 * float(np.nanmax(deviation)))
    start, stop = longest_true_run(deviation > tolerance)
    onset = float(time_ms[start])
    offset = float(time_ms[min(stop, len(time_ms) - 1)])
    current = float(np.median(command_pA[start:stop]) - baseline)
    return onset, offset, current


def cell_id_from_path(group: str, path: Path) -> str:
    suffix = path.stem.split("_")[-1]
    return f"{group}_{suffix}"


def find_cc_files(root: Path) -> list[Path]:
    files = sorted({*root.rglob("cc_*.abf"), *root.rglob("CC_*.abf")})
    if not files:
        files = sorted(root.rglob("*.abf"))
    return [p for p in files if "capa" not in p.name.lower()]


def load_abf_sweeps(group: str, path: Path, cfg: dict) -> list[dict]:
    import pyabf
    abf = pyabf.ABF(str(path))
    cell_id = cell_id_from_path(group, path)
    rows = []
    for sweep_index in map(int, abf.sweepList):
        abf.setSweep(sweep_index, channel=int(cfg["command_channel"]))
        time_ms = np.asarray(abf.sweepX, float) * 1000.0
        command = np.asarray(abf.sweepC, float)
        try:
            onset, offset, current = stimulus_from_command(time_ms, command)
        except ValueError:
            abf.setSweep(int(abf.sweepList[-1]), channel=int(cfg["command_channel"]))
            ref_t = np.asarray(abf.sweepX, float) * 1000.0
            ref_c = np.asarray(abf.sweepC, float)
            onset, offset, _ = stimulus_from_command(ref_t, ref_c)
            abf.setSweep(sweep_index, channel=int(cfg["command_channel"]))
            time_ms = np.asarray(abf.sweepX, float) * 1000.0
            command = np.asarray(abf.sweepC, float)
            baseline_mask = time_ms < onset
            baseline = float(np.median(command[baseline_mask])) if baseline_mask.any() else float(command[0])
            stim_mask = (time_ms >= onset) & (time_ms < offset)
            current = float(np.median(command[stim_mask]) - baseline)
        abf.setSweep(sweep_index, channel=int(cfg["voltage_channel"]))
        voltage = np.asarray(abf.sweepY, float)
        rows.append({"group": group, "cell_id": cell_id, "abf_path": str(path.resolve()),
                     "sweep_index": sweep_index, "time_ms": time_ms, "voltage_mV": voltage,
                     "current_pA": current, "onset_ms": onset, "offset_ms": offset})
    return rows


def candidate_peaks(time_ms: np.ndarray, voltage: np.ndarray, cfg: dict) -> pd.DataFrame:
    dt = float(np.median(np.diff(time_ms)))
    distance = max(1, int(round(float(cfg["minimum_peak_distance_ms"]) / dt)))
    indices, _ = find_peaks(voltage, distance=distance)
    prominences = peak_prominences(voltage, indices)[0] if len(indices) else np.array([])
    floor_p = min(map(float, cfg["prominences_mV"]))
    floor_h = min(map(float, cfg["heights_mV"]))
    keep = (prominences >= floor_p) & (voltage[indices] >= floor_h)
    indices, prominences = indices[keep], prominences[keep]
    return pd.DataFrame({"peak_index": indices.astype(int), "time_ms": time_ms[indices],
                         "peak_voltage_mV": voltage[indices], "prominence_mV": prominences})


def tag_decisions(events: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    out = events.copy()
    for p in map(float, cfg["prominences_mV"]):
        for h in map(float, cfg["heights_mV"]):
            key = f"pass_p{p:g}_h{h:g}"
            out[key] = ((out.prominence_mV >= p) & (out.peak_voltage_mV >= h)).astype(int)
    pp, ph = float(cfg["primary_prominence_mV"]), float(cfg["primary_height_mV"])
    out["pass_primary"] = ((out.prominence_mV >= pp) & (out.peak_voltage_mV >= ph)).astype(int)
    out["pass_strict_height"] = ((out.prominence_mV >= pp) & (out.peak_voltage_mV >= -10.0)).astype(int)
    out["event_class"] = np.select(
        [out.pass_strict_height.eq(1), out.pass_primary.eq(1)],
        ["primary_and_strict", "primary_lost_at_minus10"],
        default="relaxed_only")
    return out


def analyze_sweep(sweep: dict, cfg: dict) -> tuple[pd.DataFrame, list[dict]]:
    events = tag_decisions(candidate_peaks(sweep["time_ms"], sweep["voltage_mV"], cfg), cfg)
    for key in ("group", "cell_id", "abf_path", "sweep_index", "current_pA", "onset_ms", "offset_ms"):
        events[key] = sweep[key]
    trim = float(cfg.get("stimulus_edge_trim_ms", 5.0))
    stim_start, stim_stop = sweep["onset_ms"] + trim, sweep["offset_ms"] - trim
    steady_cfg = [float(x) for x in cfg.get("steady_window_absolute_ms", [stim_start, stim_stop])]
    steady_start, steady_stop = max(stim_start, steady_cfg[0]), min(stim_stop, steady_cfg[1])
    if steady_stop <= steady_start:
        steady_start, steady_stop = stim_start, stim_stop
    events["inside_stimulus_window"] = ((events.time_ms >= stim_start) & (events.time_ms < stim_stop)).astype(int)
    events["inside_steady_window"] = ((events.time_ms >= steady_start) & (events.time_ms < steady_stop)).astype(int)
    duration_s = max((steady_stop - steady_start) / 1000.0, 1e-12)
    counts = []
    for p in map(float, cfg["prominences_mV"]):
        for h in map(float, cfg["heights_mV"]):
            col = f"pass_p{p:g}_h{h:g}"
            n = int(((events[col] == 1) & (events.inside_steady_window == 1)).sum())
            counts.append({"group": sweep["group"], "cell_id": sweep["cell_id"],
                           "abf_path": sweep["abf_path"], "sweep_index": sweep["sweep_index"],
                           "current_pA": sweep["current_pA"], "prominence_mV": p,
                           "height_mV": h, "spike_count": n, "rate_hz": n / duration_s,
                           "steady_start_ms": steady_start, "steady_stop_ms": steady_stop,
                           "analysis_duration_s": duration_s})
    return events, counts


def plot_trace_pages(sweeps: list[dict], event_table: pd.DataFrame, cfg: dict, output: Path) -> None:
    npp = int(cfg.get("trace_pages_sweeps_per_page", 4))
    for (group, cell_id), cell_sweeps in pd.DataFrame([
        {"group": s["group"], "cell_id": s["cell_id"], "i": i} for i, s in enumerate(sweeps)
    ]).groupby(["group", "cell_id"], sort=True):
        selected = [sweeps[int(i)] for i in cell_sweeps.i]
        selected.sort(key=lambda s: (s["current_pA"], s["sweep_index"]))
        target = output / group / f"{cell_id}_all_sweeps.pdf"
        target.parent.mkdir(parents=True, exist_ok=True)
        with PdfPages(target) as pdf:
            for page_start in range(0, len(selected), npp):
                page = selected[page_start:page_start + npp]
                fig, axes = plt.subplots(len(page), 1, figsize=(12, 2.6 * len(page)), squeeze=False)
                for ax, sweep in zip(axes[:, 0], page):
                    ev = event_table[(event_table.abf_path == sweep["abf_path"]) &
                                     (event_table.sweep_index == sweep["sweep_index"])]
                    ax.plot(sweep["time_ms"], sweep["voltage_mV"], color="0.2", lw=0.65)
                    ax.axvspan(sweep["onset_ms"], sweep["offset_ms"], color="#eeeeee", zorder=-5)
                    for h, color in [(-30, "#9ecae1"), (-20, "#fcae91"), (-10, "#fdd0a2")]:
                        ax.axhline(h, color=color, ls=":", lw=0.8)
                    styles = {
                        "primary_and_strict": dict(marker=".", color="black", s=18, zorder=5),
                        "primary_lost_at_minus10": dict(marker="o", facecolors="none", edgecolors="#e66101", s=34, zorder=6),
                        "relaxed_only": dict(marker="x", color="#0571b0", s=25, zorder=5),
                    }
                    for cls, style in styles.items():
                        part = ev[ev.event_class == cls]
                        ax.scatter(part.time_ms, part.peak_voltage_mV, **style)
                    n_primary = int(((ev.pass_primary == 1) & (ev.inside_steady_window == 1)).sum())
                    n_strict = int(((ev.pass_strict_height == 1) & (ev.inside_steady_window == 1)).sum())
                    ax.set_title(f"{group} {cell_id} | sweep {sweep['sweep_index']} | {sweep['current_pA']:.1f} pA | primary={n_primary}, strict={n_strict}", loc="left", fontsize=9)
                    ax.set_ylabel("mV")
                axes[-1, 0].set_xlabel("Time (ms)")
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)


def balanced_zoom_selection(events: pd.DataFrame, maximum: int) -> pd.DataFrame:
    ambiguous = events[(events.inside_stimulus_window == 1) &
                       (events.event_class != "primary_and_strict")].copy()
    if len(ambiguous) <= maximum:
        return ambiguous.sort_values(["group", "cell_id", "current_pA", "time_ms"])
    ambiguous["stratum"] = (ambiguous.group.astype(str) + "|" + ambiguous.cell_id.astype(str) + "|" +
                             ambiguous.event_class.astype(str))
    pieces = []
    strata = ambiguous.stratum.unique()
    quota = max(1, maximum // len(strata))
    for _, block in ambiguous.groupby("stratum", sort=True):
        positions = np.linspace(0, len(block) - 1, min(quota, len(block))).round().astype(int)
        pieces.append(block.sort_values(["current_pA", "time_ms"]).iloc[np.unique(positions)])
    return pd.concat(pieces).head(maximum)


def plot_event_zooms(selected: pd.DataFrame, sweeps: list[dict], cfg: dict, target: Path) -> None:
    lookup = {(s["abf_path"], s["sweep_index"]): s for s in sweeps}
    target.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(target) as pdf:
        if selected.empty:
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.text(.5, .5, "No threshold-disputed events found", ha="center", va="center")
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
        for _, event in selected.iterrows():
            s = lookup[(event.abf_path, int(event.sweep_index))]
            before, after = float(cfg.get("zoom_before_ms", 4)), float(cfg.get("zoom_after_ms", 7))
            mask = (s["time_ms"] >= event.time_ms - before) & (s["time_ms"] <= event.time_ms + after)
            fig, ax = plt.subplots(figsize=(8, 4.5))
            ax.plot(s["time_ms"][mask] - event.time_ms, s["voltage_mV"][mask], color="black", lw=1.2)
            ax.scatter([0], [event.peak_voltage_mV], color="#e66101" if event.pass_primary else "#0571b0", s=45, zorder=5)
            for h in (-30, -20, -10):
                ax.axhline(h, color="0.65", ls=":", lw=0.8)
            ax.set(title=(f"{event.group} {event.cell_id} | {event.current_pA:.1f} pA | "
                          f"peak={event.peak_voltage_mV:.1f} mV | prominence={event.prominence_mV:.1f} mV | {event.event_class}"),
                   xlabel="Time from candidate peak (ms)", ylabel="Voltage (mV)")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def plot_fi(counts: pd.DataFrame, events: pd.DataFrame, cfg: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    pp, ph = float(cfg["primary_prominence_mV"]), float(cfg["primary_height_mV"])
    current_min, current_max = float(cfg.get("current_min_pA", 0)), float(cfg.get("current_max_pA", 600))
    primary = counts[(counts.prominence_mV == pp) & (counts.height_mV == ph) &
                     counts.current_pA.between(current_min, current_max)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, group in zip(axes, ("WT", "SCA3")):
        block = primary[primary.group == group]
        for cell, cell_df in block.groupby("cell_id"):
            curve = cell_df.groupby("current_pA", as_index=False).rate_hz.mean().sort_values("current_pA")
            ax.plot(curve.current_pA, curve.rate_hz, marker="o", ms=3, lw=1, alpha=.75, label=cell)
        ax.set(title=f"{group}: every cell, primary detector", xlabel="Injected current (pA)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=7, ncol=2)
    axes[0].set_ylabel("Firing rate (Hz)")
    fig.tight_layout()
    fig.savefig(output / "fi_every_cell_primary.pdf")
    plt.close(fig)

    height_block = counts[(counts.prominence_mV == pp) & counts.current_pA.between(current_min, current_max)].copy()
    per_cell = height_block.groupby(["group", "cell_id", "current_pA", "height_mV"], as_index=False).rate_hz.mean()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    palette = {-30.0: "#2166ac", -20.0: "#7f7f7f", -10.0: "#b2182b"}
    for ax, group in zip(axes, ("WT", "SCA3")):
        block = per_cell[per_cell.group == group]
        for h, hdf in block.groupby("height_mV"):
            summary = hdf.groupby("current_pA").rate_hz.agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)]).reset_index()
            summary.columns = ["current_pA", "median", "q25", "q75"]
            color = palette.get(float(h), "black")
            ax.plot(summary.current_pA, summary["median"], marker="o", color=color, label=f"height {h:g} mV")
            ax.fill_between(summary.current_pA, summary.q25, summary.q75, color=color, alpha=.15)
        ax.set(title=f"{group}: median and cell IQR", xlabel="Injected current (pA)")
        ax.grid(alpha=.2); ax.legend()
    axes[0].set_ylabel("Firing rate (Hz)")
    fig.tight_layout(); fig.savefig(output / "fi_height_sensitivity.pdf"); plt.close(fig)

    inside = events[(events.inside_stimulus_window == 1) & (events.prominence_mV >= min(cfg["prominences_mV"]))]
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.arange(-30, 42, 2)
    for group in ("WT", "SCA3"):
        ax.hist(inside.loc[inside.group == group, "peak_voltage_mV"], bins=bins, histtype="step", lw=2,
                density=True, color=GROUP_COLORS[group], label=group)
    for h in (-20, -10): ax.axvline(h, color="0.4", ls="--")
    ax.set(xlabel="Candidate peak voltage (mV)", ylabel="Density", title="Peak-height distribution of relaxed detector candidates")
    ax.legend(); fig.tight_layout(); fig.savefig(output / "candidate_peak_voltage_distribution.pdf"); plt.close(fig)


def write_report(sweeps: list[dict], events: pd.DataFrame, counts: pd.DataFrame, output: Path, cfg: dict) -> None:
    pp, ph = float(cfg["primary_prominence_mV"]), float(cfg["primary_height_mV"])
    primary = counts[(counts.prominence_mV == pp) & (counts.height_mV == ph)]
    strict = counts[(counts.prominence_mV == pp) & (counts.height_mV == -10.0)]
    lines = ["# Spike visual QC", "", f"Sweeps: {len(sweeps)}", f"Candidate events: {len(events)}", "",
             "## Counts inside analysis window", ""]
    joined = primary.merge(strict, on=["group", "cell_id", "abf_path", "sweep_index", "current_pA"], suffixes=("_primary", "_strict"))
    summary = joined.groupby("group").agg(primary_events=("spike_count_primary", "sum"), strict_events=("spike_count_strict", "sum"))
    summary["fraction_lost_at_minus10"] = 1 - summary.strict_events / summary.primary_events.replace(0, np.nan)
    lines.extend(["```", summary.to_string(), "```", "", "No inferential tests or p-values are computed here.",
                  "Use trace pages and event zooms to decide whether the orange WT events are genuine action potentials."])
    (output / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wt-root", type=Path, required=True)
    parser.add_argument("--sca3-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config.json"))
    parser.add_argument("--output", type=Path, default=Path("results_spike_visual_qc"))
    args = parser.parse_args()
    cfg = load_config(args.config)
    args.output.mkdir(parents=True, exist_ok=True)
    sweeps = []
    for group, root in (("WT", args.wt_root), ("SCA3", args.sca3_root)):
        files = find_cc_files(root)
        if not files:
            raise SystemExit(f"No current-clamp ABFs found under {root}")
        for i, path in enumerate(files, 1):
            print(f"[{group} {i}/{len(files)}] {path}", flush=True)
            sweeps.extend(load_abf_sweeps(group, path, cfg))
    event_parts, count_rows = [], []
    for s in sweeps:
        ev, ct = analyze_sweep(s, cfg)
        event_parts.append(ev); count_rows.extend(ct)
    events = pd.concat(event_parts, ignore_index=True) if event_parts else pd.DataFrame()
    counts = pd.DataFrame(count_rows)
    tables = args.output / "01_tables"; tables.mkdir(exist_ok=True)
    events.to_csv(tables / "spike_candidates.csv", index=False)
    counts.to_csv(tables / "detector_counts.csv", index=False)
    inventory = pd.DataFrame([{k: s[k] for k in ("group", "cell_id", "abf_path", "sweep_index", "current_pA", "onset_ms", "offset_ms")} for s in sweeps])
    inventory.to_csv(tables / "sweep_inventory.csv", index=False)
    plot_trace_pages(sweeps, events, cfg, args.output / "02_trace_pages")
    selected = balanced_zoom_selection(events, int(cfg.get("maximum_zoom_events", 300)))
    selected.to_csv(tables / "zoomed_ambiguous_events.csv", index=False)
    plot_event_zooms(selected, sweeps, cfg, args.output / "03_event_zooms" / "ambiguous_event_zooms.pdf")
    plot_fi(counts, events, cfg, args.output / "04_fi_visual")
    write_report(sweeps, events, counts, args.output, cfg)
    print(f"Completed: {args.output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
