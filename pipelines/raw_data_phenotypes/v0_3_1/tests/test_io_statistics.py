from pathlib import Path

import numpy as np
import pandas as pd

from neurothermo_phenotypes.config import load_config
from neurothermo_phenotypes.extract import analyse_trace
from neurothermo_phenotypes.io import Trace, load_traces, read_manifest
from neurothermo_phenotypes.pipeline import _filter_common_protocol, _finalize_event_audit
from neurothermo_phenotypes.statistics import (
    compare_two_part_all_cells,
    integrated_cell_features,
    rheobase_brackets,
)


def test_csv_loader_and_rheobase(tmp_path: Path):
    t = np.linspace(0, 1.2, 1201)
    rows = []
    for sweep, current in enumerate([0, 100]):
        v = np.full_like(t, -70.0)
        if current:
            for center in [0.3, 0.6, 0.9]:
                v += 90 * np.exp(-0.5 * ((t - center) / .001) ** 2)
        rows.append(pd.DataFrame({"time_s": t, "voltage_mV": v, "current_pA": current, "sweep_index": sweep}))
    trace_path = tmp_path / "cell.csv"
    pd.concat(rows).to_csv(trace_path, index=False)
    manifest_path = tmp_path / "manifest.csv"
    pd.DataFrame([{"group": "WT", "cell_id": "WT_01", "path": "cell.csv", "capacitance_pF": 100, "include": True}]).to_csv(manifest_path, index=False)
    cfg = load_config(None)
    cfg["input"]["data_root"] = "."
    traces = load_traces(read_manifest(manifest_path), manifest_path, cfg)
    assert len(traces) == 2
    features = pd.DataFrame({"group": ["WT", "WT"], "cell_id": ["WT_01", "WT_01"], "current_pA": [0, 100], "n_spikes": [0, 3], "qc_pass": [True, True]})
    bracket = rheobase_brackets(features).iloc[0]
    assert bracket.rheobase_lower_nonspiking_pA == 0
    assert bracket.rheobase_upper_spiking_pA == 100


def test_user_data_root_expands_tilde(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    data_root = tmp_path / "neurothermo"
    wt_root = data_root / "WT"
    wt_root.mkdir(parents=True)

    trace_path = wt_root / "cc_01.csv"
    pd.DataFrame({
        "time_s": [0.0, 0.001, 0.002],
        "voltage_mV": [-70.0, -69.0, -70.0],
        "current_pA": [100.0, 100.0, 100.0],
        "sweep_index": [0, 0, 0],
    }).to_csv(trace_path, index=False)

    manifest_path = data_root / "manifest.csv"
    pd.DataFrame([{
        "group": "WT", "cell_id": "WT_01", "path": "WT/cc_01.csv",
        "capacitance_pF": 100.0, "include": True,
    }]).to_csv(manifest_path, index=False)

    cfg = load_config(None)
    cfg["input"]["data_root"] = "~/neurothermo"
    traces = load_traces(read_manifest("~/neurothermo/manifest.csv"), "~/neurothermo/manifest.csv", cfg)
    assert len(traces) == 1
    assert Path(traces[0].source_path) == trace_path


def test_curated_events_and_frozen_invariants(tmp_path: Path):
    trace_path = tmp_path / "cell.csv"
    time = np.linspace(0.0, 1.2, 1201)
    pd.concat([
        pd.DataFrame({
            "time_s": time, "voltage_mV": np.full(len(time), -70.0),
            "current_pA": np.full(len(time), current),
            "sweep_index": np.full(len(time), sweep, int),
        })
        for sweep, current in [(0, 100.0), (1, 200.0)]
    ]).to_csv(trace_path, index=False)
    manifest_path = tmp_path / "manifest.csv"
    pd.DataFrame([{
        "group": "WT", "cell_id": "WT_01", "path": trace_path.name,
        "capacitance_pF": 100.0, "include": True,
    }]).to_csv(manifest_path, index=False)
    events_path = tmp_path / "events.csv"
    pd.DataFrame([
        {
            "group": "WT", "cell_id": "WT_01", "sweep_index": 0,
            "current_pA": 100.0, "onset_ms": 100.0, "offset_ms": 1100.0,
            "time_ms": 300.0, "fixed_qc_detected": True,
        },
        {
            "group": "WT", "cell_id": "WT_01", "sweep_index": 1,
            "current_pA": 200.0, "onset_ms": 100.0, "offset_ms": 1100.0,
            "time_ms": 400.0, "fixed_qc_detected": True,
        },
    ]).to_csv(events_path, index=False)
    frozen_path = tmp_path / "frozen.csv"
    pd.DataFrame([{
        "sweep_id": "WT_01__sw000__I100pA", "group": "WT", "cell_id": "WT_01",
        "sweep_index": 0, "current_pA": 100.0, "stim_onset_ms": 100.0,
        "stim_offset_ms": 1100.0, "final_audit_decision": "ACCEPT",
    }]).to_csv(frozen_path, index=False)
    cfg = load_config(None)
    cfg["input"].update({
        "data_root": str(tmp_path), "curated_events_csv": str(events_path),
        "require_curated_events": True, "curated_sweeps_manifest_csv": str(frozen_path),
        "analysis_currents_pA": [100.0, 200.0], "expected_curated_sweeps": 1,
        "expected_curated_spike_events": 1,
        "restrict_curated_events_to_frozen_sweeps": True,
    })
    traces, audit = load_traces(
        read_manifest(manifest_path), manifest_path, cfg, return_event_audit=True
    )
    assert len(traces) == 2
    assert np.allclose(traces[0].metadata["curated_spikes_s"], [0.3])
    assert traces[0].metadata["frozen_sweep_membership"] is True
    assert traces[0].stim_start_override_s == 0.1
    assert traces[0].stim_end_override_s == 1.1
    assert len(traces[1].metadata["curated_spikes_s"]) == 0
    assert traces[1].metadata["frozen_sweep_membership"] is False
    assert np.allclose(traces[1].metadata["threshold_probe_spikes_s"], [0.4])
    statuses = audit.set_index("event_row_id")["event_status"].to_dict()
    assert statuses == {
        0: "accepted_frozen_after_override",
        1: "excluded_not_frozen_sweep",
    }
    features = pd.DataFrame([analyse_trace(trace, cfg) for trace in traces])
    finalized = _finalize_event_audit(audit, features, cfg).set_index("event_row_id")
    assert bool(finalized.loc[0, "used_in_analysis"])
    assert bool(finalized.loc[0, "used_for_threshold_probe"])
    assert not bool(finalized.loc[1, "used_in_analysis"])
    assert bool(finalized.loc[1, "used_for_threshold_probe"])


def test_curated_events_outside_stimulus_are_counted():
    time = np.linspace(0.0, 1.2, 2401)
    trace = Trace(
        group="WT", cell_id="WT_01", record_id="WT_01", sweep_index=0,
        source_path="cell.csv", time_s=time,
        voltage_mV=np.full(len(time), -70.0), current_trace_pA=None,
        current_pA=100.0, capacitance_pF=100.0,
        metadata={
            "curated_spikes_s": np.array([0.05, 0.5, 1.15]),
            "curated_events_loaded": 3,
            "frozen_sweep_membership": True,
            "threshold_probe_spikes_s": np.array([0.05, 0.5, 1.15]),
        },
    )
    cfg = load_config(None)
    cfg["stimulus"].update({"detection": "fixed", "fixed_start_s": 0.1, "fixed_end_s": 1.1})
    row = analyse_trace(trace, cfg)
    assert row["curated_events_loaded"] == 3
    assert row["curated_events_used"] == 1
    assert row["curated_events_outside_stimulus_window"] == 2
    assert row["threshold_probe_events_loaded"] == 3
    assert row["threshold_probe_event_count"] == 1
    assert row["n_spikes"] == 1


def test_curated_boundary_tolerance_keeps_frozen_peak():
    time = np.linspace(0.0, 1.2, 2401)
    trace = Trace(
        group="WT", cell_id="WT_01", record_id="WT_01", sweep_index=0,
        source_path="cell.csv", time_s=time,
        voltage_mV=np.full(len(time), -70.0), current_trace_pA=None,
        current_pA=100.0, capacitance_pF=100.0,
        metadata={
            "curated_spikes_s": np.array([1.1008]),
            "curated_events_loaded": 1,
            "frozen_sweep_membership": True,
            "threshold_probe_spikes_s": np.array([1.1008]),
        },
    )
    cfg = load_config(None)
    cfg["stimulus"].update({"detection": "fixed", "fixed_start_s": 0.1, "fixed_end_s": 1.1})
    cfg["input"]["curated_event_boundary_tolerance_ms"] = 1.0
    row = analyse_trace(trace, cfg)
    assert row["curated_events_used"] == 1
    assert row["curated_events_outside_stimulus_window"] == 0
    assert row["threshold_probe_event_count"] == 1


def test_common_protocol_filter_audits_extra_sweeps():
    def trace(sweep, current):
        return Trace(
            group="WT", cell_id="WT_01", record_id="WT_01", sweep_index=sweep,
            source_path="cell.csv", time_s=np.array([0.0, 1.0]),
            voltage_mV=np.array([-70.0, -70.0]), current_trace_pA=None,
            current_pA=current, capacitance_pF=100.0,
        )
    cfg = load_config(None)
    cfg["input"].update({
        "analysis_currents_pA": [0.0, 50.0], "enforce_common_current_grid": True,
        "require_complete_current_grid": True,
    })
    kept, excluded = _filter_common_protocol([trace(0, 0), trace(1, 50), trace(2, 100)], cfg)
    assert [x.current_pA for x in kept] == [0.0, 50.0]
    assert excluded["current_pA"].tolist() == [100]


def test_fixed_domain_integration_requires_complete_grid():
    cfg = load_config(None)
    cfg["input"]["analysis_currents_pA"] = [0.0, 50.0, 100.0]
    rows = []
    for current, rate in [(0.0, 0.0), (50.0, 10.0), (100.0, 20.0)]:
        rows.append({
            "group": "WT", "cell_id": "WT_01", "current_pA": current,
            "qc_pass": True, "external_work_signed_fJ": current,
            "external_work_positive_fJ": current, "mean_power_signed_fW": current,
            "firing_rate_hz": rate, "sustained_rate_hz": rate,
            "baseline_voltage_mV": -70.0,
        })
    integrated, coverage = integrated_cell_features(pd.DataFrame(rows), cfg, return_coverage=True)
    assert integrated["firing_rate_hz__auc_mean_0_100pA"].iloc[0] == 10.0
    assert coverage["complete_common_domain"].all()


def test_two_part_all_cells_test_is_exact():
    cfg = load_config(None)
    cfg["statistics"].update({
        "minimum_cells_per_group": 3,
        "two_part_currents_pA": [400.0, 500.0],
        "two_part_primary_features": ["predictive_information_nats"],
        "two_part_secondary_features": [],
        "exact_max_labelings": 200000,
    })
    rows = []
    for group, offset in [("WT", 0.0), ("SCA3", 10.0)]:
        for cell in range(3):
            value = offset + float(cell)
            for current in [400.0, 500.0]:
                rows.append({
                    "group": group, "cell_id": f"{group}_{cell}",
                    "current_pA": current, "qc_pass": True,
                    "thermo_eligible": True, "curated_frozen_sweep": True,
                    "predictive_information_nats": (
                        np.nan if group == "SCA3" and cell == 2 else value
                    ),
                })
        for current in [400.0, 500.0]:
            rows.append({
                "group": group, "cell_id": f"{group}_all_fatal",
                "current_pA": current, "qc_pass": False,
                "thermo_eligible": False, "curated_frozen_sweep": False,
                "predictive_information_nats": np.nan,
            })
    comparisons, coverage, cells = compare_two_part_all_cells(pd.DataFrame(rows), cfg)
    result = comparisons.iloc[0]
    assert result["permutation_mode"] == "exact"
    assert result["valid_labelings"] == 70
    assert result["n_WT_total"] == 4
    assert result["n_SCA3_total"] == 4
    assert result["n_WT_any_defined"] == 3
    assert result["n_SCA3_any_defined"] == 2
    assert coverage["total_cells"].tolist() == [4, 4]
    assert len(cells) == 16
