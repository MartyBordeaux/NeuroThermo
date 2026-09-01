from __future__ import annotations
import math
import numpy as np
from numba import njit

@njit(cache=True)
def rhs(x, y, z, J, b, r, s, kappa_I, a, c, d, x_R):
    dx = y - a*x*x*x + b*x*x - z + kappa_I*J
    dy = c - d*x*x - y
    dz = r * (s*(x - x_R) - z)
    return dx, dy, dz

@njit(cache=True)
def rk4_step(x, y, z, h, J, b, r, s, kappa_I, a, c, d, x_R):
    k1x, k1y, k1z = rhs(x,y,z,J,b,r,s,kappa_I,a,c,d,x_R)
    k2x, k2y, k2z = rhs(x+.5*h*k1x,y+.5*h*k1y,z+.5*h*k1z,J,b,r,s,kappa_I,a,c,d,x_R)
    k3x, k3y, k3z = rhs(x+.5*h*k2x,y+.5*h*k2y,z+.5*h*k2z,J,b,r,s,kappa_I,a,c,d,x_R)
    k4x, k4y, k4z = rhs(x+h*k3x,y+h*k3y,z+h*k3z,J,b,r,s,kappa_I,a,c,d,x_R)
    return (x+h*(k1x+2*k2x+2*k3x+k4x)/6.0,
            y+h*(k1y+2*k2y+2*k3y+k4y)/6.0,
            z+h*(k1z+2*k2z+2*k3z+k4z)/6.0)

@njit(cache=True)
def pre_relax_core(pre_ms, dt_ms, time_scale_ms,
                   b,r,s,kappa_I,a,c,d,x_R,x0,y0,z0):
    h=dt_ms/time_scale_ms
    n=int(math.ceil(pre_ms/dt_ms))
    x,y,z=x0,y0,z0
    for _ in range(n):
        x,y,z=rk4_step(x,y,z,h,0.0,b,r,s,kappa_I,a,c,d,x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)):
            return x,y,z,False
    return x,y,z,True

@njit(cache=True)
def has_spike_from_state_core(x0,y0,z0,J,stim_ms,dt_ms,time_scale_ms,
                              b,r,s,kappa_I,a,c,d,x_R,threshold,refractory_ms):
    h=dt_ms/time_scale_ms
    n=int(math.ceil(stim_ms/dt_ms))
    x,y,z=x0,y0,z0
    xprev=x
    last=-1e30
    for i in range(1,n+1):
        x,y,z=rk4_step(x,y,z,h,J,b,r,s,kappa_I,a,c,d,x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:
            return False, False
        ti=i*dt_ms
        if xprev < threshold and x >= threshold and ti-last >= refractory_ms:
            return True, True
        xprev=x
    return False, True

@njit(cache=True)
def supported_spike_metrics_core(x0,y0,z0,J,active_window_ms,max_onset_ms,dt_ms,time_scale_ms,
                                 b,r,s,kappa_I,a,c,d,x_R,threshold,refractory_ms,max_spikes):
    """Measure suprathreshold dynamics in a window that starts at the first model spike.

    Current is held constant while searching for the first spike and throughout the
    active window. This intentionally removes absolute onset latency from the dynamic
    coordinate; absolute excitability is represented separately by rheobase.
    """
    h=dt_ms/time_scale_ms
    total_ms=max_onset_ms+active_window_ms
    n=int(math.ceil(total_ms/dt_ms))
    x,y,z=x0,y0,z0
    xprev=x
    last_cross=-1e30
    first=np.nan
    spikes=np.empty(max_spikes,dtype=np.float64)
    ns=0
    for i in range(1,n+1):
        x,y,z=rk4_step(x,y,z,h,J,b,r,s,kappa_I,a,c,d,x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:
            return 0,np.nan,np.nan,np.nan,np.nan,False
        ti=i*dt_ms
        if xprev < threshold and x >= threshold and ti-last_cross >= refractory_ms:
            den=x-xprev
            frac=(threshold-xprev)/den if den != 0 else 0.0
            tc=(i-1+frac)*dt_ms
            last_cross=tc
            if not np.isfinite(first):
                if tc > max_onset_ms:
                    return 0,np.nan,np.nan,np.nan,np.nan,True
                first=tc
            rel=tc-first
            if rel <= active_window_ms + 1e-9:
                if ns < max_spikes:
                    spikes[ns]=rel
                    ns+=1
            else:
                break
        xprev=x
        if np.isfinite(first) and ti-first > active_window_ms+dt_ms:
            break
    if ns==0:
        return 0,0.0,np.nan,0.0,np.nan,True
    rate=ns/(active_window_ms/1000.0) if active_window_ms>0 else np.nan
    if ns>=2:
        ssum=0.0
        for j in range(1,ns):
            ssum += spikes[j]-spikes[j-1]
        mean_isi=ssum/(ns-1)
        train=spikes[ns-1]-spikes[0]
        occupancy=train/active_window_ms if active_window_ms>0 else np.nan
    else:
        mean_isi=np.nan
        train=0.0
        occupancy=0.0
    return ns,rate,mean_isi,occupancy,first,True


def _pars(theta,cfg):
    m=cfg['model']
    return (float(theta['b']),float(theta['r']),float(theta['s']),float(theta['kappa_I']),
            float(m['a']),float(m['c']),float(m['d']),float(m['x_R']))


def pre_relax(theta,cfg,dt_ms=None):
    m=cfg['model']
    if dt_ms is None: dt_ms=float(cfg['simulation']['dt_ms'])
    b,r,s,k,a,c,d,xr=_pars(theta,cfg)
    return pre_relax_core(float(m['pre_ms']),float(dt_ms),float(m['model_time_scale_ms']),
                          b,r,s,k,a,c,d,xr,float(m['x0']),float(m['y0']),float(m['z0']))


def has_spike_from_state(state,theta,J,cfg,dt_ms=None):
    m=cfg['model']; sim=cfg['rheobase']
    if dt_ms is None: dt_ms=float(sim['dt_ms'])
    b,r,s,k,a,c,d,xr=_pars(theta,cfg)
    sp,ok=has_spike_from_state_core(float(state[0]),float(state[1]),float(state[2]),float(J),
        float(sim['stimulus_ms']),float(dt_ms),float(m['model_time_scale_ms']),
        b,r,s,k,a,c,d,xr,float(m['model_spike_threshold']),float(m['model_refractory_ms']))
    return bool(sp),bool(ok)


def refine_rheobase(theta,cfg,guess=None,pre_state=None):
    """Adaptive global rheobase search in current density J (pA/pF)."""
    rcfg=cfg['rheobase']
    if pre_state is None:
        x,y,z,ok=pre_relax(theta,cfg,dt_ms=float(rcfg['dt_ms']))
        if not ok: return {'rheobase_J':np.nan,'status':'PRE_RELAX_FAIL','iterations':0}
        state=(x,y,z)
    else:
        state=pre_state
    jmin=float(rcfg['search_min_J']); jmax=float(rcfg['search_max_J'])
    mult=float(rcfg['bracket_multiplier']); add=float(rcfg.get('bracket_add_J',0.02))
    if guess is None or not np.isfinite(guess): guess=max(jmin+add,0.5*(jmin+jmax))
    guess=min(max(float(guess),jmin+1e-12),jmax)
    sp,ok=has_spike_from_state(state,theta,guess,cfg)
    if not ok: return {'rheobase_J':np.nan,'status':'SIM_FAIL','iterations':0}
    probes=1
    if sp:
        hi=guess; lo=max(jmin,guess/mult-add)
        while lo>jmin+1e-12:
            spl,ok=has_spike_from_state(state,theta,lo,cfg); probes+=1
            if not ok: return {'rheobase_J':np.nan,'status':'SIM_FAIL','iterations':probes}
            if not spl: break
            hi=lo; lo=max(jmin,lo/mult-add)
        if lo<=jmin+1e-12:
            spl,ok=has_spike_from_state(state,theta,jmin,cfg); probes+=1
            if spl: return {'rheobase_J':jmin,'status':'AT_LOWER_SEARCH_BOUND','iterations':probes}
            lo=jmin
    else:
        lo=guess; hi=min(jmax,guess*mult+add)
        while True:
            sph,ok=has_spike_from_state(state,theta,hi,cfg); probes+=1
            if not ok: return {'rheobase_J':np.nan,'status':'SIM_FAIL','iterations':probes}
            if sph: break
            if hi>=jmax-1e-12:
                return {'rheobase_J':np.nan,'status':'NO_SPIKE_WITHIN_SEARCH_RANGE','iterations':probes}
            lo=hi; hi=min(jmax,hi*mult+add)
    tol=float(rcfg['tolerance_J']); maxit=int(rcfg['max_iterations'])
    it=0
    while hi-lo>tol and it<maxit:
        mid=.5*(lo+hi)
        spm,ok=has_spike_from_state(state,theta,mid,cfg); probes+=1
        if not ok: return {'rheobase_J':np.nan,'status':'SIM_FAIL','iterations':probes}
        if spm: hi=mid
        else: lo=mid
        it+=1
    return {'rheobase_J':float(hi),'status':'OK','iterations':probes}


def supported_metrics(theta,J,active_window_ms,cfg,pre_state=None):
    m=cfg['model']; scfg=cfg['simulation']
    if pre_state is None:
        x,y,z,ok=pre_relax(theta,cfg,dt_ms=float(scfg['dt_ms']))
        if not ok:
            return {'spike_count':0,'support_rate_hz':np.nan,'mean_isi_ms':np.nan,'occupancy_fraction':np.nan,'first_spike_ms':np.nan,'simulation_ok':False}
        state=(x,y,z)
    else:
        state=pre_state
    b,r,s,k,a,c,d,xr=_pars(theta,cfg)
    max_spikes=max(64,int((float(active_window_ms)+float(scfg['max_onset_ms']))/max(float(m['model_refractory_ms']),float(scfg['dt_ms'])))+20)
    n,rate,isi,occ,first,ok=supported_spike_metrics_core(float(state[0]),float(state[1]),float(state[2]),float(J),
        float(active_window_ms),float(scfg['max_onset_ms']),float(scfg['dt_ms']),float(m['model_time_scale_ms']),
        b,r,s,k,a,c,d,xr,float(m['model_spike_threshold']),float(m['model_refractory_ms']),int(max_spikes))
    return {'spike_count':int(n),'support_rate_hz':float(rate) if np.isfinite(rate) else np.nan,
            'mean_isi_ms':float(isi) if np.isfinite(isi) else np.nan,
            'occupancy_fraction':float(occ) if np.isfinite(occ) else np.nan,
            'first_spike_ms':float(first) if np.isfinite(first) else np.nan,'simulation_ok':bool(ok)}
