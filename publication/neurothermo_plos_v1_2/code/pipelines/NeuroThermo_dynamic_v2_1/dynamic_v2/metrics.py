from __future__ import annotations
import numpy as np


def spike_metrics(spikes, duration_ms):
    s=np.asarray(spikes,float)
    n=len(s)
    out={'spike_count':n,'firing_rate_hz':n/(duration_ms/1000.0) if duration_ms>0 else np.nan,
         'active_rate_hz':np.nan,'latency_ms':np.nan,'train_duration_ms':np.nan,
         'mean_isi_ms':np.nan,'median_isi_ms':np.nan,'cv_isi':np.nan,
         'early_isi_ms':np.nan,'late_isi_ms':np.nan,'adaptation_index':np.nan}
    if n:
        out['latency_ms']=float(s[0])
    if n>=2:
        isi=np.diff(s)
        train=float(s[-1]-s[0])
        out['train_duration_ms']=train
        out['mean_isi_ms']=float(np.mean(isi)); out['median_isi_ms']=float(np.median(isi))
        if train>0: out['active_rate_hz']=float((n-1)/(train/1000.0))
        if np.mean(isi)>0: out['cv_isi']=float(np.std(isi,ddof=0)/np.mean(isi))
        k=min(3,len(isi)); early=float(np.mean(isi[:k])); late=float(np.mean(isi[-k:]))
        out['early_isi_ms']=early; out['late_isi_ms']=late
        den=late+early
        if den != 0: out['adaptation_index']=float((late-early)/den)
    return out


def symmetric_relative_difference(a,b,eps=1e-12):
    if not (np.isfinite(a) and np.isfinite(b)): return np.nan
    return float(2.0*abs(a-b)/(abs(a)+abs(b)+eps))


def align_first_spike(exp_spikes, model_spikes):
    e=np.asarray(exp_spikes,float); m=np.asarray(model_spikes,float)
    if len(e)==0 or len(m)==0:
        return m.copy(), np.nan
    tau=float(e[0]-m[0])
    return m+tau, tau


def restrict_spikes(spikes,start_ms,end_ms):
    s=np.asarray(spikes,float)
    return s[(s>=float(start_ms)) & (s<=float(end_ms))]
