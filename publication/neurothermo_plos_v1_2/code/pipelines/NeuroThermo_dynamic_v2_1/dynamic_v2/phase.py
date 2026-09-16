from __future__ import annotations
import numpy as np
from .model import divergence, speed

PHASE_KEYS=['x','y','z','speed','divergence']


def phase_profile(t,x,y,z,spikes,J,theta,cfg):
    grid=np.linspace(0.0,1.0,int(cfg['phase']['n_bins']))
    sp=np.asarray(spikes,float)
    if len(sp)<3:
        return None
    intervals=list(range(len(sp)-1))
    skip=int(cfg['phase']['skip_initial_cycles'])
    if len(intervals)>2*skip:
        intervals=intervals[skip:]
    cycles=[]
    v=speed(x,y,z,J,theta,cfg); div=divergence(x,theta,cfg)
    for i in intervals:
        a,b=sp[i],sp[i+1]
        mask=(t>=a)&(t<=b)
        if mask.sum()<5 or b<=a: continue
        ph=(t[mask]-a)/(b-a)
        cycles.append({'x':np.interp(grid,ph,x[mask]),'y':np.interp(grid,ph,y[mask]),
                       'z':np.interp(grid,ph,z[mask]),'speed':np.interp(grid,ph,v[mask]),
                       'divergence':np.interp(grid,ph,div[mask]),'period_ms':b-a})
    if len(cycles)<int(cfg['phase']['min_cycles']): return None
    prof={'phase':grid,'n_cycles':len(cycles),'period_ms':float(np.median([c['period_ms'] for c in cycles]))}
    for key in PHASE_KEYS:
        prof[key]=np.median(np.stack([c[key] for c in cycles]),axis=0)
    return prof


def supported_phase_profile(t,x,y,z,raw_model_spikes,aligned_model_spikes,exp_spikes,J,theta,cfg):
    """Phase profile from model cycles whose aligned boundaries lie inside the experimental train.

    The state trajectory is never shifted or rescaled.  Alignment is used only to decide which raw
    model spikes correspond to the experimentally supported temporal region.
    """
    exp=np.asarray(exp_spikes,float); raw=np.asarray(raw_model_spikes,float); ali=np.asarray(aligned_model_spikes,float)
    if len(exp)<2 or len(raw)!=len(ali): return None
    lo=float(exp[0]); hi=float(exp[-1])
    keep=(ali>=lo)&(ali<=hi)
    supported_raw=raw[keep]
    return phase_profile(t,x,y,z,supported_raw,J,theta,cfg)


def phase_descriptors(prof):
    if prof is None:
        return {'phase_cycles':0,'cycle_period_ms':np.nan,'x_range':np.nan,'y_range':np.nan,'z_range':np.nan,
                'mean_speed':np.nan,'peak_speed':np.nan,'mean_divergence':np.nan,'min_divergence':np.nan,
                'max_divergence':np.nan,'fraction_positive_divergence':np.nan,'cycle_arc_length':np.nan,
                'z_peak_phase':np.nan,'speed_peak_phase':np.nan}
    p=prof['phase']; x=prof['x']; y=prof['y']; z=prof['z']; v=prof['speed']; d=prof['divergence']
    arc=float(np.sum(np.sqrt(np.diff(x)**2+np.diff(y)**2+np.diff(z)**2)))
    return {'phase_cycles':int(prof['n_cycles']),'cycle_period_ms':float(prof['period_ms']),
            'x_range':float(np.ptp(x)),'y_range':float(np.ptp(y)),'z_range':float(np.ptp(z)),
            'mean_speed':float(np.mean(v)),'peak_speed':float(np.max(v)),
            'mean_divergence':float(np.mean(d)),'min_divergence':float(np.min(d)),
            'max_divergence':float(np.max(d)),'fraction_positive_divergence':float(np.mean(d>0)),
            'cycle_arc_length':arc,'z_peak_phase':float(p[int(np.argmax(z))]),
            'speed_peak_phase':float(p[int(np.argmax(v))])}


def profile_nrmse(best,alt,key,eps=1e-12):
    if best is None or alt is None: return np.nan
    a=np.asarray(best[key],float); b=np.asarray(alt[key],float)
    scale=float(np.ptp(a))
    if scale<eps: scale=max(float(np.std(a)),eps)
    return float(np.sqrt(np.mean((a-b)**2))/scale)
