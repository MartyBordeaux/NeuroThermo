import numpy as np
import pandas as pd
from spike_qc_calibrated import apply_fixed_qc, normalize_qc_note, stimulus_from_command, sustained_mask


def test_stimulus():
    t = np.arange(0, 1000, 0.1)
    c = np.zeros_like(t)
    c[(t >= 100) & (t < 900)] = 350
    onset, offset, current = stimulus_from_command(t, c)
    assert abs(onset - 100) < 0.2
    assert abs(offset - 900) < 0.2
    assert abs(current - 350) < 1e-9


def test_sustained_train_removes_isolated_event():
    times = np.array([10., 400., 410., 420., 800.])
    got = sustained_mask(times, max_isi_ms=25, minimum_spikes=2)
    assert got.tolist() == [False, True, True, True, False]


def test_normalize_qc_note_mixed_alphabet():
    assert normalize_qc_note("Последниe 2 detected – rejected ") == "последние 2 detected - rejected"


def test_fixed_qc_manual_remove_and_add():
    events = pd.DataFrame({
        "group": ["WT"]*5, "cell_id": ["WT_01"]*5,
        "current_pA": [500.0]*5, "time_ms": [100., 200., 300., 400., 500.],
        "peak_voltage_mV": [-5., -4., -3., -2., -1.],
        "detected": [True, True, True, False, True],
    })
    sweeps = [{"group": "WT", "cell_id": "WT_01", "current_pA": 500.0}]
    qc = pd.DataFrame({"type": ["wt"], "cell": [1], "sweep": [500],
                       "conclusion": ["последний detected – rejected\nrejected после первых трех пиков – detected"]})
    fixed, _ = apply_fixed_qc(events, sweeps, qc, {"manual_time_match_tolerance_ms": 75})
    assert fixed.algorithm_detected.tolist() == [True, True, True, False, True]
    assert fixed.fixed_qc_detected.tolist() == [True, True, True, True, False]
