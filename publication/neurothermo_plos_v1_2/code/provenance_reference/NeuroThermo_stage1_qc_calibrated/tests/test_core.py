import numpy as np
from spike_qc_calibrated import stimulus_from_command, sustained_mask


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
