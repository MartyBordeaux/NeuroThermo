from __future__ import annotations
import numpy as np


def fit_window_end_ms(exp_spikes_ms, stim_duration_ms, cfg):
    spikes=np.asarray(exp_spikes_ms,dtype=float)
    if spikes.size==0:return 0.0
    last=float(spikes[-1]);local_n=int(cfg.get('local_isi_count',5))
    if spikes.size>=2:
        isi=np.diff(spikes);local=isi[-min(local_n,isi.size):];typ=float(np.median(local));guard=float(cfg.get('post_last_spike_guard_isi_multiplier',1.0))*typ
    else:guard=float(cfg.get('single_spike_guard_ms',15.0))
    guard=max(float(cfg.get('min_guard_ms',5.0)),min(float(cfg.get('max_guard_ms',30.0)),guard))
    return min(float(stim_duration_ms)+float(cfg.get('max_post_stimulus_ms',15.0)),last+guard)
