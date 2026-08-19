from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from collections import Counter
from typing import Iterator, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


@dataclass
class Trace:
    group: str
    cell_id: str
    record_id: str
    sweep_index: int
    source_path: str
    time_s: np.ndarray
    voltage_mV: np.ndarray
    current_trace_pA: Optional[np.ndarray]
    current_pA: float
    capacitance_pF: float
    capacitance_10ms_pF: float = np.nan
    capacitance_20ms_pF: float = np.nan
    capacitance_50ms_pF: float = np.nan
    animal_id: str = ""
    stim_start_override_s: float = np.nan
    stim_end_override_s: float = np.nan
    metadata: dict = field(default_factory=dict)


REQUIRED_MANIFEST = {"group", "cell_id", "path", "capacitance_pF"}


def _truthy(value) -> bool:
    if pd.isna(value):
        return True
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() not in {"0", "false", "no", "exclude", "fail"}


def _included_events(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    return series.fillna("").astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes", "y", "pass", "include", "included"}
    )


def read_manifest(path: Union[str, Path]) -> pd.DataFrame:
    path = Path(path).expanduser()
    frame = pd.read_csv(path)
    missing = sorted(REQUIRED_MANIFEST.difference(frame.columns))
    if missing:
        raise ValueError(f"Manifest lacks required columns: {', '.join(missing)}")
    if "include" in frame:
        frame = frame[frame["include"].map(_truthy)].copy()
    frame["group"] = frame["group"].astype(str).str.upper()
    frame["cell_id"] = frame["cell_id"].astype(str)
    if "record_id" not in frame:
        frame["record_id"] = frame["cell_id"] + "_CC"
    return frame.reset_index(drop=True)


def make_manifest_from_qc(
    qc_xlsx: Union[str, Path],
    raw_root: Union[str, Path],
    cohort_csv: Optional[Union[str, Path]] = None,
    expected_included: Optional[int] = None,
) -> pd.DataFrame:
    qc_xlsx = Path(qc_xlsx).expanduser()
    raw_root = Path(raw_root).expanduser()
    provenance = pd.read_excel(qc_xlsx, sheet_name="Provenance")
    capacitance = pd.read_excel(qc_xlsx, sheet_name="Capacitance QC")
    provenance.columns = [str(c).strip() for c in provenance.columns]
    capacitance.columns = [str(c).strip() for c in capacitance.columns]
    keep_cap = [
        "group", "cell_id", "capacitance_pF", "capacitance_10ms_pF",
        "capacitance_20ms_pF", "capacitance_50ms_pF", "capacitance_qc",
    ]
    keep_cap = [c for c in keep_cap if c in capacitance.columns]
    cap = capacitance[keep_cap].drop_duplicates(["group", "cell_id"])
    prov = provenance.drop_duplicates(["group", "cell_id"]).copy()
    out = prov.merge(cap, on=["group", "cell_id"], how="left", suffixes=("", "_cap"))
    out["path"] = out["cc_file"].astype(str)
    out["path_exists"] = out["path"].map(lambda p: (raw_root / p).exists())
    overall = out.get("overall_qc", pd.Series("PASS", index=out.index)).astype(str).str.upper()
    in_benchmark = out.get("in_new_benchmark", pd.Series(True, index=out.index)).map(_truthy)
    out["include"] = overall.eq("PASS") & in_benchmark & out["path_exists"]
    out["cohort_membership"] = True
    if cohort_csv:
        cohort = pd.read_csv(Path(cohort_csv).expanduser())
        if not {"group", "cell_id"}.issubset(cohort.columns):
            raise ValueError("Cohort CSV requires group and cell_id")
        cohort_keys = set(zip(cohort["group"].astype(str).str.upper(), cohort["cell_id"].astype(str)))
        available_keys = set(zip(out["group"].astype(str).str.upper(), out["cell_id"].astype(str)))
        missing = sorted(cohort_keys - available_keys)
        if missing:
            raise ValueError(f"Frozen cohort cells absent from QC workbook: {missing}")
        out["cohort_membership"] = [
            (str(g).upper(), str(c)) in cohort_keys for g, c in zip(out["group"], out["cell_id"])
        ]
        out["include"] &= out["cohort_membership"]
    if expected_included is not None and int(out["include"].sum()) != int(expected_included):
        raise ValueError(
            f"Included cohort has {int(out['include'].sum())} cells; expected {int(expected_included)}"
        )
    if "record_id" not in out:
        out["record_id"] = out["cell_id"].astype(str) + "_CC"
    columns = [
        "group", "cell_id", "record_id", "animal_id", "path",
        "capacitance_pF", "capacitance_10ms_pF", "capacitance_20ms_pF",
        "capacitance_50ms_pF", "overall_qc", "capacitance_qc",
        "include", "cohort_membership", "path_exists", "exclusion_reason", "notes",
    ]
    for col in columns:
        if col not in out:
            out[col] = "" if col not in {"capacitance_pF", "capacitance_10ms_pF", "capacitance_20ms_pF", "capacitance_50ms_pF"} else np.nan
    return out[columns].sort_values(["group", "cell_id"]).reset_index(drop=True)


def read_sweep_overrides(path: Optional[Union[str, Path]]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    frame = pd.read_csv(Path(path).expanduser())
    required = {"cell_id", "sweep_index"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Sweep overrides lack: {', '.join(sorted(missing))}")
    return frame


def _unit_scale(unit: str, target: str) -> float:
    unit = str(unit).strip().lower()
    if target == "mV":
        return {"v": 1000.0, "mv": 1.0, "uv": 0.001, "µv": 0.001}.get(unit, 1.0)
    return {"a": 1e12, "ma": 1e9, "ua": 1e6, "µa": 1e6, "na": 1e3, "pa": 1.0}.get(unit, 1.0)


def _protocol_current(cfg: dict, sweep_index: int) -> float:
    values = cfg["input"].get("protocol_currents_pA")
    if values is None or sweep_index >= len(values):
        return np.nan
    return float(values[sweep_index])


def _load_abf(path: Path, row: pd.Series, cfg: dict) -> Iterator[Trace]:
    try:
        import pyabf
    except ImportError as exc:
        raise RuntimeError("ABF input requires installation with: pip install -e '.[abf]'") from exc
    abf = pyabf.ABF(str(path))
    channel = int(cfg["input"].get("voltage_channel", 0))
    for sweep_index in abf.sweepList:
        abf.setSweep(int(sweep_index), channel=channel)
        time_s = np.asarray(abf.sweepX, dtype=float).copy()
        voltage = np.asarray(abf.sweepY, dtype=float).copy()
        voltage *= _unit_scale(getattr(abf, "sweepUnitsY", "mV"), "mV")
        try:
            command = np.asarray(abf.sweepC, dtype=float).copy()
            command *= _unit_scale(getattr(abf, "sweepUnitsC", "pA"), "pA")
            if command.shape != voltage.shape:
                command = None
        except Exception:
            command = None
        configured = _protocol_current(cfg, int(sweep_index))
        current = configured if cfg["input"].get("current_source") == "configured" else np.nan
        if not np.isfinite(current) and command is not None:
            n0 = max(10, len(command) // 20)
            baseline = float(np.median(command[:n0]))
            current = float(command[np.argmax(np.abs(command - baseline))] - baseline)
        yield _trace_from_arrays(row, path, int(sweep_index), time_s, voltage, command, current)


def _load_csv(path: Path, row: pd.Series, cfg: dict) -> Iterator[Trace]:
    frame = pd.read_csv(path)
    if not {"time_s", "voltage_mV"}.issubset(frame.columns):
        raise ValueError(f"CSV must be long-form with time_s and voltage_mV: {path}")
    if "sweep_index" not in frame:
        frame["sweep_index"] = 0
    for sweep_index, sweep in frame.groupby("sweep_index", sort=True):
        command = sweep["current_trace_pA"].to_numpy(float) if "current_trace_pA" in sweep else None
        current = float(sweep["current_pA"].median()) if "current_pA" in sweep else _protocol_current(cfg, int(sweep_index))
        yield _trace_from_arrays(
            row, path, int(sweep_index), sweep["time_s"].to_numpy(float),
            sweep["voltage_mV"].to_numpy(float), command, current,
        )


def _load_npz(path: Path, row: pd.Series, cfg: dict) -> Iterator[Trace]:
    data = np.load(path, allow_pickle=False)
    time = np.asarray(data["time_s"], dtype=float)
    voltage = np.asarray(data["voltage_mV"], dtype=float)
    if voltage.ndim == 1:
        voltage = voltage[None, :]
    command = np.asarray(data["current_trace_pA"], dtype=float) if "current_trace_pA" in data else None
    currents = np.asarray(data["current_pA"], dtype=float) if "current_pA" in data else None
    for sweep_index in range(voltage.shape[0]):
        t = time if time.ndim == 1 else time[sweep_index]
        cmd = None if command is None else (command if command.ndim == 1 else command[sweep_index])
        current = float(currents[sweep_index]) if currents is not None and currents.ndim else (
            float(currents) if currents is not None else _protocol_current(cfg, sweep_index)
        )
        yield _trace_from_arrays(row, path, sweep_index, t, voltage[sweep_index], cmd, current)


def _number(row: pd.Series, name: str) -> float:
    value = row.get(name, np.nan)
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _trace_from_arrays(row, path, sweep_index, time, voltage, command, current) -> Trace:
    return Trace(
        group=str(row["group"]), cell_id=str(row["cell_id"]),
        record_id=str(row.get("record_id", row["cell_id"])), sweep_index=int(sweep_index),
        source_path=str(path), time_s=np.asarray(time, float), voltage_mV=np.asarray(voltage, float),
        current_trace_pA=None if command is None else np.asarray(command, float), current_pA=float(current),
        capacitance_pF=_number(row, "capacitance_pF"),
        capacitance_10ms_pF=_number(row, "capacitance_10ms_pF"),
        capacitance_20ms_pF=_number(row, "capacitance_20ms_pF"),
        capacitance_50ms_pF=_number(row, "capacitance_50ms_pF"),
        animal_id=str(row.get("animal_id", "")),
    )


def load_traces(
    manifest: pd.DataFrame,
    manifest_path: Union[str, Path],
    cfg: dict,
    return_event_audit: bool = False,
) -> Union[List[Trace], Tuple[List[Trace], pd.DataFrame]]:
    manifest_path = Path(manifest_path).expanduser()
    root = Path(cfg["input"].get("data_root", ".")).expanduser()
    if not root.is_absolute():
        root = manifest_path.resolve().parent / root
    overrides = read_sweep_overrides(cfg["input"].get("sweep_overrides_csv"))
    override_map = {}
    for _, row in overrides.iterrows():
        override_map[(str(row.cell_id), int(row.sweep_index))] = row
    (
        event_map, threshold_event_map, peak_override_map, frozen_sweep_keys,
        frozen_window_map, curated_mode, event_audit,
    ) = _load_curated_spikes(
        cfg["input"], manifest_path.resolve().parent
    )
    traces: list[Trace] = []
    loaded_keys = set()
    for _, row in manifest.iterrows():
        path = Path(str(row["path"])).expanduser()
        if not path.is_absolute():
            path = root / path
        if not path.exists():
            raise FileNotFoundError(f"Trace file not found: {path}")
        suffix = path.suffix.lower()
        loader = {".abf": _load_abf, ".csv": _load_csv, ".npz": _load_npz}.get(suffix)
        if loader is None:
            raise ValueError(f"Unsupported input format {suffix}: {path}")
        for trace in loader(path, row, cfg):
            trace_key = (trace.group.upper(), trace.cell_id, trace.sweep_index)
            loaded_keys.add(trace_key)
            override = override_map.get((trace.cell_id, trace.sweep_index))
            if override is not None:
                if "include" in override and not _truthy(override["include"]):
                    continue
                if "current_pA" in override and pd.notna(override["current_pA"]):
                    trace.current_pA = float(override["current_pA"])
                if "stim_start_s" in override and pd.notna(override["stim_start_s"]):
                    trace.stim_start_override_s = float(override["stim_start_s"])
                if "stim_end_s" in override and pd.notna(override["stim_end_s"]):
                    trace.stim_end_override_s = float(override["stim_end_s"])
            if curated_mode:
                spikes = np.sort(np.asarray(event_map.get(trace_key, np.array([], dtype=float)), float))
                action = peak_override_map.get(trace_key, "")
                trace.metadata["curated_spikes_s"] = spikes
                trace.metadata["curated_peak_override"] = action
                trace.metadata["frozen_sweep_membership"] = trace_key in frozen_sweep_keys
                trace.metadata["curated_events_loaded"] = int(len(spikes))
                threshold_spikes = np.sort(np.asarray(
                    threshold_event_map.get(trace_key, np.array([], dtype=float)), float
                ))
                trace.metadata["threshold_probe_spikes_s"] = threshold_spikes
                trace.metadata["threshold_probe_events_loaded"] = int(len(threshold_spikes))
                if trace_key in frozen_window_map:
                    trace.stim_start_override_s, trace.stim_end_override_s = frozen_window_map[trace_key]
            traces.append(trace)
    missing_frozen = sorted(frozen_sweep_keys - loaded_keys)
    if missing_frozen:
        raise ValueError(f"Frozen curated sweeps are absent from loaded recordings: {missing_frozen}")
    if return_event_audit:
        return traces, event_audit
    return traces


def _resolve_input_path(value, base_dir: Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else base_dir / path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_curated_hashes(input_cfg: dict, base_dir: Path) -> None:
    value = input_cfg.get("curated_hash_manifest_json")
    if not value:
        return
    manifest_path = _resolve_input_path(value, base_dir)
    expected = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key in [
        "curated_sweeps_manifest_csv", "curated_peak_overrides_csv",
        "curated_threshold_brackets_csv",
    ]:
        path_value = input_cfg.get(key)
        if not path_value:
            continue
        path = _resolve_input_path(path_value, base_dir)
        expected_hash = expected.get(path.name)
        if expected_hash is None:
            raise ValueError(f"No frozen SHA-256 entry for {path.name}")
        actual_hash = _sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(f"Frozen calibration hash mismatch: {path}")


def _apply_peak_action(spikes: np.ndarray, action: str, sweep_id: str) -> np.ndarray:
    selected = np.sort(np.asarray(spikes, float))
    for part in action.split("+") if action else []:
        part = part.strip().upper()
        if part == "DROP_FIRST":
            if len(selected) < 2:
                raise ValueError(f"DROP_FIRST would remove the only spike in {sweep_id}")
            selected = selected[1:]
        elif part == "DROP_LAST":
            if len(selected) < 2:
                raise ValueError(f"DROP_LAST would remove the only spike in {sweep_id}")
            selected = selected[:-1]
        elif part not in {"", "NONE"}:
            raise ValueError(f"Unknown curated peak override {part!r} for {sweep_id}")
    return selected


def _load_curated_spikes(
    input_cfg: dict, base_dir: Path
) -> Tuple[dict, dict, dict, set, dict, bool, pd.DataFrame]:
    _validate_curated_hashes(input_cfg, base_dir)
    events_path = input_cfg.get("curated_events_csv")
    if not events_path:
        if input_cfg.get("require_curated_events", False):
            raise FileNotFoundError("Production analysis requires input.curated_events_csv")
        return {}, {}, {}, set(), {}, False, pd.DataFrame()
    events_path = _resolve_input_path(events_path, base_dir)
    if not events_path.exists():
        raise FileNotFoundError(f"Curated Stage-1 event table not found: {events_path}")
    events = pd.read_csv(events_path)
    required = {"group", "cell_id", "sweep_index", "time_ms"}
    missing = required.difference(events.columns)
    if missing:
        raise ValueError(f"Curated events lack: {', '.join(sorted(missing))}")
    include_col = input_cfg.get("events_include_column", "fixed_qc_detected")
    if include_col not in events:
        raise ValueError(f"Curated events lack include column: {include_col}")
    events = events[_included_events(events[include_col])].copy()
    events["group"] = events["group"].astype(str).str.upper()
    events["cell_id"] = events["cell_id"].astype(str)
    events["sweep_index"] = events["sweep_index"].astype(int)
    events["event_row_id"] = events.index.astype(int)
    raw_event_map = {
        (str(g).upper(), str(c), int(s)): np.sort(part["time_ms"].to_numpy(float)) / 1000.0
        for (g, c, s), part in events.groupby(["group", "cell_id", "sweep_index"])
    }
    event_map = {key: values.copy() for key, values in raw_event_map.items()}

    frozen_keys = set()
    frozen_all_keys = set()
    frozen_ids = {}
    frozen_window_map = {}
    frozen_path = input_cfg.get("curated_sweeps_manifest_csv")
    if frozen_path:
        frozen = pd.read_csv(_resolve_input_path(frozen_path, base_dir))
        required_frozen = {
            "sweep_id", "group", "cell_id", "sweep_index", "current_pA",
            "stim_onset_ms", "stim_offset_ms", "final_audit_decision",
        }
        missing = required_frozen.difference(frozen.columns)
        if missing:
            raise ValueError(f"Frozen sweep manifest lacks: {', '.join(sorted(missing))}")
        frozen = frozen[frozen["final_audit_decision"].astype(str).str.upper().eq("ACCEPT")].copy()
        frozen_all_keys = set(zip(
            frozen["group"].astype(str).str.upper(),
            frozen["cell_id"].astype(str),
            frozen["sweep_index"].astype(int),
        ))
        allowed = input_cfg.get("analysis_currents_pA")
        if allowed is not None:
            tolerance = float(input_cfg.get("current_tolerance_pA", 1e-6))
            allowed_values = np.asarray(allowed, float)
            frozen = frozen[
                frozen["current_pA"].astype(float).map(
                    lambda value: bool(np.any(np.abs(allowed_values - value) <= tolerance))
                )
            ].copy()
        if frozen["sweep_id"].duplicated().any():
            raise ValueError("Duplicate sweep_id in frozen curated sweep manifest")
        expected_sweeps = input_cfg.get("expected_curated_sweeps")
        if expected_sweeps is not None and len(frozen) != int(expected_sweeps):
            raise ValueError(f"Frozen current-matched sweep count is {len(frozen)}; expected {expected_sweeps}")
        strict = bool(input_cfg.get("strict_curated_metadata", True))
        metadata_tolerance = float(input_cfg.get("curated_metadata_tolerance_ms", 0.05))
        for _, frozen_row in frozen.iterrows():
            key = (
                str(frozen_row["group"]).upper(), str(frozen_row["cell_id"]),
                int(frozen_row["sweep_index"]),
            )
            if key not in event_map:
                raise ValueError(f"Frozen curated sweep has no selected Stage-1 events: {key}")
            frozen_keys.add(key)
            frozen_ids[key] = str(frozen_row["sweep_id"])
            frozen_window_map[key] = (
                float(frozen_row["stim_onset_ms"]) / 1000.0,
                float(frozen_row["stim_offset_ms"]) / 1000.0,
            )
            if strict:
                selected_rows = events[
                    events["group"].astype(str).str.upper().eq(key[0])
                    & events["cell_id"].astype(str).eq(key[1])
                    & events["sweep_index"].astype(int).eq(key[2])
                ]
                checks = []
                if "current_pA" in selected_rows:
                    checks.append(("current_pA", float(selected_rows["current_pA"].iloc[0]), float(frozen_row["current_pA"]), 1e-6))
                if "onset_ms" in selected_rows:
                    checks.append(("onset_ms", float(selected_rows["onset_ms"].iloc[0]), float(frozen_row["stim_onset_ms"]), metadata_tolerance))
                if "offset_ms" in selected_rows:
                    checks.append(("offset_ms", float(selected_rows["offset_ms"].iloc[0]), float(frozen_row["stim_offset_ms"]), metadata_tolerance))
                for name, observed, expected, tolerance in checks:
                    if abs(observed - expected) > tolerance:
                        raise ValueError(
                            f"{frozen_row['sweep_id']}: Stage-1 {name}={observed} disagrees with frozen {expected}"
                        )

    peak_map = {}
    peak_path = input_cfg.get("curated_peak_overrides_csv")
    if peak_path:
        peak_path = _resolve_input_path(peak_path, base_dir)
        peaks = pd.read_csv(peak_path)
        if not {"sweep_id", "action"}.issubset(peaks.columns):
            raise ValueError("Curated peak override CSV requires sweep_id and action")
        actions_by_id = dict(zip(peaks["sweep_id"].astype(str), peaks["action"].fillna("").astype(str)))
        for key, sweep_id in frozen_ids.items():
            action = actions_by_id.get(sweep_id, "").strip().upper()
            if action:
                event_map[key] = _apply_peak_action(event_map[key], action, sweep_id)
                peak_map[key] = action

    expected_events = input_cfg.get("expected_curated_spike_events")
    if expected_events is not None:
        actual_events = sum(len(event_map[key]) for key in frozen_keys)
        if actual_events != int(expected_events):
            raise ValueError(
                f"Frozen current-matched spike-event count is {actual_events}; expected {expected_events}"
            )
    threshold_event_map = {key: values.copy() for key, values in event_map.items()}
    restrict = bool(input_cfg.get("restrict_curated_events_to_frozen_sweeps", False))
    audit_rows = []
    for key, part in events.groupby(["group", "cell_id", "sweep_index"], sort=True):
        normalized_key = (str(key[0]).upper(), str(key[1]), int(key[2]))
        retained = Counter(np.round(event_map.get(normalized_key, np.array([], float)), 12))
        action = peak_map.get(normalized_key, "")
        for _, event in part.sort_values(["time_ms", "event_row_id"]).iterrows():
            time_s = float(event["time_ms"]) / 1000.0
            rounded = round(time_s, 12)
            if normalized_key in frozen_keys:
                if retained[rounded] > 0:
                    status = "accepted_frozen_after_override"
                    retained[rounded] -= 1
                else:
                    status = "excluded_peak_override"
            elif normalized_key in frozen_all_keys:
                status = "excluded_outside_common_current_grid"
            else:
                status = "excluded_not_frozen_sweep"
            audit_rows.append({
                "event_row_id": int(event["event_row_id"]),
                "group": normalized_key[0], "cell_id": normalized_key[1],
                "sweep_index": normalized_key[2],
                "current_pA": float(event["current_pA"]) if "current_pA" in event and pd.notna(event["current_pA"]) else np.nan,
                "time_ms": float(event["time_ms"]),
                "frozen_common_domain_sweep": normalized_key in frozen_keys,
                "peak_override_action": action,
                "event_status": status,
                "used_in_analysis": False,
            })
    event_audit = pd.DataFrame(audit_rows)
    if restrict:
        event_map = {key: event_map[key] for key in frozen_keys}
    return (
        event_map, threshold_event_map, peak_map, frozen_keys,
        frozen_window_map, True, event_audit,
    )
