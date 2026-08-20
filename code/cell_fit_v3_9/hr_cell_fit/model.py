from __future__ import annotations
import numpy as np
from numba import njit

@njit(cache=True)
def _rhs(x,y,z,J,b,r,s,kappa_I,a,c,d,x_R):
    return (y-a*x*x*x+b*x*x-z+kappa_I*J,c-d*x*x-y,r*(s*(x-x_R)-z))

@njit(cache=True)
def _rk4(x,y,z,h,J,b,r,s,kappa_I,a,c,d,x_R):
    k1=_rhs(x,y,z,J,b,r,s,kappa_I,a,c,d,x_R)
    k2=_rhs(x+0.5*h*k1[0],y+0.5*h*k1[1],z+0.5*h*k1[2],J,b,r,s,kappa_I,a,c,d,x_R)
    k3=_rhs(x+0.5*h*k2[0],y+0.5*h*k2[1],z+0.5*h*k2[2],J,b,r,s,kappa_I,a,c,d,x_R)
    k4=_rhs(x+h*k3[0],y+h*k3[1],z+h*k3[2],J,b,r,s,kappa_I,a,c,d,x_R)
    return (x+h*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,y+h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6,z+h*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6)

@njit(cache=True)
def simulate_spikes_core(J,stim_ms,post_ms,dt_ms,b,r,s,kappa_I,a,c,d,x_R,x0,y0,z0,pre_ms,time_scale_ms,threshold,refractory_ms,max_spikes):
    h=dt_ms/time_scale_ms;x=x0;y=y0;z=z0
    npre=int(round(pre_ms/dt_ms))
    for _ in range(npre):
        x,y,z=_rk4(x,y,z,h,0.0,b,r,s,kappa_I,a,c,d,x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:return np.empty(0),False
    total_ms=stim_ms+post_ms;n=int(np.ceil(total_ms/dt_ms));out=np.empty(max_spikes,dtype=np.float64);nsp=0;xp=x;last=-1e30
    for i in range(1,n+1):
        t=(i-1)*dt_ms;drive=J if t<=stim_ms+1e-12 else 0.0
        x,y,z=_rk4(x,y,z,h,drive,b,r,s,kappa_I,a,c,d,x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:return out[:nsp],False
        if xp<threshold and x>=threshold:
            denom=x-xp;frac=(threshold-xp)/denom if denom!=0 else 0.0;tc=(i-1+frac)*dt_ms
            if tc-last>=refractory_ms:
                if nsp<max_spikes:out[nsp]=tc;nsp+=1
                last=tc
        xp=x
    return out[:nsp],True


def simulate_spikes(J,stim_ms,post_ms,dt_ms,params,model_cfg):
    max_spikes=max(128,int((stim_ms+post_ms)/max(model_cfg['model_refractory_ms'],dt_ms))+10)
    return simulate_spikes_core(float(J),float(stim_ms),float(post_ms),float(dt_ms),float(params['b']),float(params['r']),float(params['s']),float(params['kappa_I']),
        float(model_cfg['a']),float(model_cfg['c']),float(model_cfg['d']),float(model_cfg['x_R']),float(model_cfg['x0']),float(model_cfg['y0']),float(model_cfg['z0']),
        float(model_cfg['pre_ms']),float(model_cfg['model_time_scale_ms']),float(model_cfg['model_spike_threshold']),float(model_cfg['model_refractory_ms']),int(max_spikes))
