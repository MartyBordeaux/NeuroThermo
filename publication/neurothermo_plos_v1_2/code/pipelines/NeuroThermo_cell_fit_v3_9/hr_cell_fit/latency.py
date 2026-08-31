from __future__ import annotations
import numpy as np
from .vp import victor_purpura


def _raw_spikes_in_model_interval(model_spikes_ms, max_ms):
    model=np.asarray(model_spikes_ms,dtype=float)
    return np.sort(model[(model>=0.0)&(model<=float(max_ms)+1e-9)])


def exact_first_spike_align(exp_spikes_ms, model_spikes_ms, max_ms, tau_ms=10.0, normalize=True, count_penalty_weight=0.25):
    exp=np.sort(np.asarray(exp_spikes_ms,dtype=float));raw=_raw_spikes_in_model_interval(model_spikes_ms,max_ms)
    if exp.size==0:raise ValueError('exact_first_spike_align expects a spiking experimental sweep')
    if raw.size==0:
        count_pen=float(count_penalty_weight)*abs(len(exp)-0)/max(len(exp),1)
        return {'tau_ms':float('nan'),'aligned_model_spike_times_ms':np.asarray([],dtype=float),'raw_model_spike_times_ms':raw,
                'vp_loss':float('inf'),'count_penalty':count_pen,'loss':float('inf'),'n_raw_spikes':0,'n_aligned_spikes':0,
                'raw_first_ms':np.nan,'aligned_first_ms':np.nan,'aligned_last_ms':np.nan,'exp_last_ms':float(exp[-1]),'aligned_last_spike_error_ms':np.nan,
                'raw_vp_loss':victor_purpura(exp,raw,tau_ms=tau_ms,normalize=normalize)}
    tau=float(exp[0]-raw[0]);aligned=raw+tau
    vp= victor_purpura(exp,aligned,tau_ms=tau_ms,normalize=normalize)
    raw_vp= victor_purpura(exp,raw,tau_ms=tau_ms,normalize=normalize)
    count_pen=float(count_penalty_weight)*abs(len(exp)-len(raw))/max(len(exp),1)
    return {'tau_ms':tau,'aligned_model_spike_times_ms':aligned,'raw_model_spike_times_ms':raw,'vp_loss':vp,'raw_vp_loss':raw_vp,
            'count_penalty':count_pen,'loss':vp+count_pen,'n_raw_spikes':int(raw.size),'n_aligned_spikes':int(aligned.size),'raw_first_ms':float(raw[0]),
            'aligned_first_ms':float(aligned[0]),'aligned_last_ms':float(aligned[-1]),'exp_last_ms':float(exp[-1]),
            'aligned_last_spike_error_ms':float(aligned[-1]-exp[-1])}
