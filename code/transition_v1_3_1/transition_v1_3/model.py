from __future__ import annotations
import math
import numpy as np
from numba import njit

@njit(cache=True)
def deriv(x,y,z,J,b,r,s,kappa_I,a,c,d,x_R):
    dx=y-a*x*x*x+b*x*x-z+kappa_I*J
    dy=c-d*x*x-y
    dz=r*(s*(x-x_R)-z)
    return dx,dy,dz

@njit(cache=True)
def rk4_step(x,y,z,h,J,b,r,s,kappa_I,a,c,d,x_R):
    k1=deriv(x,y,z,J,b,r,s,kappa_I,a,c,d,x_R)
    k2=deriv(x+.5*h*k1[0],y+.5*h*k1[1],z+.5*h*k1[2],J,b,r,s,kappa_I,a,c,d,x_R)
    k3=deriv(x+.5*h*k2[0],y+.5*h*k2[1],z+.5*h*k2[2],J,b,r,s,kappa_I,a,c,d,x_R)
    k4=deriv(x+h*k3[0],y+h*k3[1],z+h*k3[2],J,b,r,s,kappa_I,a,c,d,x_R)
    return (x+h*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,
            y+h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,
            z+h*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6)

@njit(cache=True)
def pre_relax_core(pre_ms,dt_ms,time_scale_ms,b,r,s,k,a,c,d,xR,x0,y0,z0):
    h=dt_ms/time_scale_ms; n=int(round(pre_ms/dt_ms)); x,y,z=x0,y0,z0
    for _ in range(n):
        x,y,z=rk4_step(x,y,z,h,0.0,b,r,s,k,a,c,d,xR)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:
            return x,y,z,False
    return x,y,z,True

@njit(cache=True)
def has_spike_from_state_core(x,y,z,J,stim_ms,dt_ms,time_scale_ms,b,r,s,k,a,c,d,xR,threshold,refractory_ms):
    h=dt_ms/time_scale_ms; n=int(math.ceil(stim_ms/dt_ms)); xp=x; last=-1e30
    for i in range(1,n+1):
        x,y,z=rk4_step(x,y,z,h,J,b,r,s,k,a,c,d,xR)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:
            return False,False
        if xp < threshold and x >= threshold:
            den=x-xp; frac=(threshold-xp)/den if den!=0 else 0.0
            tc=(i-1+frac)*dt_ms
            if tc-last>=refractory_ms:return True,True
            last=tc
        xp=x
    return False,True

@njit(cache=True)
def supported_spike_metrics_core(x0,y0,z0,J,active_window_ms,max_onset_ms,dt_ms,time_scale_ms,b,r,s,k,a,c,d,xR,threshold,refractory_ms,max_spikes):
    h=dt_ms/time_scale_ms; total=max_onset_ms+active_window_ms; n=int(math.ceil(total/dt_ms))
    x,y,z=x0,y0,z0; xp=x; last=-1e30; first=np.nan; spikes=np.empty(max_spikes,np.float64); ns=0
    for i in range(1,n+1):
        x,y,z=rk4_step(x,y,z,h,J,b,r,s,k,a,c,d,xR)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:return 0,np.nan,np.nan,np.nan,np.nan,False
        ti=i*dt_ms
        if xp<threshold and x>=threshold and ti-last>=refractory_ms:
            den=x-xp; frac=(threshold-xp)/den if den!=0 else 0.0; tc=(i-1+frac)*dt_ms; last=tc
            if not np.isfinite(first):
                if tc>max_onset_ms:return 0,np.nan,np.nan,np.nan,np.nan,True
                first=tc
            rel=tc-first
            if rel<=active_window_ms+1e-9:
                if ns<max_spikes:spikes[ns]=rel;ns+=1
            else:break
        xp=x
        if np.isfinite(first) and ti-first>active_window_ms+dt_ms:break
    if ns==0:return 0,0.0,np.nan,0.0,np.nan,True
    rate=ns/(active_window_ms/1000.0) if active_window_ms>0 else np.nan
    if ns>=2:
        ss=0.0
        for j in range(1,ns):ss+=spikes[j]-spikes[j-1]
        isi=ss/(ns-1); train=spikes[ns-1]-spikes[0]; occ=train/active_window_ms if active_window_ms>0 else np.nan
    else:isi=np.nan;occ=0.0
    return ns,rate,isi,occ,first,True


def _pars(theta,cfg):
    m=cfg['model']; return float(theta['b']),float(theta['r']),float(theta['s']),float(theta['kappa_I']),float(m['a']),float(m['c']),float(m['d']),float(m['x_R'])

def pre_relax(theta,cfg,dt_ms=None):
    m=cfg['model']; dt_ms=float(cfg['simulation']['dt_ms']) if dt_ms is None else float(dt_ms)
    b,r,s,k,a,c,d,xr=_pars(theta,cfg)
    return pre_relax_core(float(m['pre_ms']),dt_ms,float(m['model_time_scale_ms']),b,r,s,k,a,c,d,xr,float(m['x0']),float(m['y0']),float(m['z0']))

def has_spike_from_state(state,theta,J,cfg,dt_ms=None):
    m=cfg['model']; rj=cfg['rheobase']; dt_ms=float(rj['dt_ms']) if dt_ms is None else float(dt_ms)
    b,r,s,k,a,c,d,xr=_pars(theta,cfg)
    return has_spike_from_state_core(float(state[0]),float(state[1]),float(state[2]),float(J),float(rj['stimulus_ms']),dt_ms,float(m['model_time_scale_ms']),b,r,s,k,a,c,d,xr,float(m['model_spike_threshold']),float(m['model_refractory_ms']))

def refine_rheobase(theta,cfg,guess=None,pre_state=None):
    rcfg=cfg['rheobase']
    if pre_state is None:
        x,y,z,ok=pre_relax(theta,cfg,dt_ms=float(rcfg['dt_ms']))
        if not ok:return {'rheobase_J':np.nan,'status':'PRE_RELAX_FAIL'}
        state=(x,y,z)
    else:state=pre_state
    lo=float(rcfg['search_min_J']); hi=float(rcfg['search_max_J'])
    mult=float(rcfg['bracket_multiplier']); add=float(rcfg.get('bracket_add_J',0.02))
    g=float(guess) if guess is not None and np.isfinite(guess) else .5*(lo+hi)
    g=min(max(g,lo+1e-12),hi)
    sp,ok=has_spike_from_state(state,theta,g,cfg)
    if not ok:return {'rheobase_J':np.nan,'status':'SIM_FAIL'}
    if sp:
        upper=g; lower=max(lo,g/mult-add)
        while lower>lo+1e-12:
            spl,ok=has_spike_from_state(state,theta,lower,cfg)
            if not ok:return {'rheobase_J':np.nan,'status':'SIM_FAIL'}
            if not spl:break
            upper=lower;lower=max(lo,lower/mult-add)
        if lower<=lo+1e-12:
            spl,ok=has_spike_from_state(state,theta,lo,cfg)
            if spl:return {'rheobase_J':lo,'status':'AT_LOWER_SEARCH_BOUND'}
            lower=lo
    else:
        lower=g;upper=min(hi,g*mult+add)
        while True:
            sph,ok=has_spike_from_state(state,theta,upper,cfg)
            if not ok:return {'rheobase_J':np.nan,'status':'SIM_FAIL'}
            if sph:break
            if upper>=hi-1e-12:return {'rheobase_J':np.nan,'status':'NO_SPIKE_WITHIN_SEARCH_RANGE'}
            lower=upper;upper=min(hi,upper*mult+add)
    for _ in range(int(rcfg['max_iterations'])):
        if upper-lower<=float(rcfg['tolerance_J']):break
        mid=.5*(lower+upper); sm,ok=has_spike_from_state(state,theta,mid,cfg)
        if not ok:return {'rheobase_J':np.nan,'status':'SIM_FAIL'}
        if sm:upper=mid
        else:lower=mid
    return {'rheobase_J':float(upper),'status':'OK'}

def supported_metrics(theta,J,active_window_ms,cfg,pre_state=None):
    m=cfg['model']; scfg=cfg['simulation']
    if pre_state is None:
        x,y,z,ok=pre_relax(theta,cfg,dt_ms=float(scfg['dt_ms']))
        if not ok:return {'spike_count':0,'support_rate_hz':np.nan,'mean_isi_ms':np.nan,'occupancy_fraction':np.nan,'first_spike_ms':np.nan,'simulation_ok':False}
        state=(x,y,z)
    else:state=pre_state
    b,r,s,k,a,c,d,xr=_pars(theta,cfg); max_spikes=2048
    n,rate,isi,occ,first,ok=supported_spike_metrics_core(float(state[0]),float(state[1]),float(state[2]),float(J),float(active_window_ms),float(scfg['max_onset_ms']),float(scfg['dt_ms']),float(m['model_time_scale_ms']),b,r,s,k,a,c,d,xr,float(m['model_spike_threshold']),float(m['model_refractory_ms']),max_spikes)
    return {'spike_count':int(n),'support_rate_hz':float(rate) if np.isfinite(rate) else np.nan,'mean_isi_ms':float(isi) if np.isfinite(isi) else np.nan,'occupancy_fraction':float(occ) if np.isfinite(occ) else np.nan,'first_spike_ms':float(first) if np.isfinite(first) else np.nan,'simulation_ok':bool(ok)}
