from __future__ import annotations
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
def simulate_trace_core(J, stimulus_duration_ms, observation_end_ms, pre_ms, dt_ms, time_scale_ms,
                        b,r,s,kappa_I,a,c,d,x_R,x0,y0,z0,
                        threshold,refractory_ms,max_spikes):
    h=dt_ms/time_scale_ms
    npre=int(np.ceil(pre_ms/dt_ms))
    x,y,z=x0,y0,z0
    for _ in range(npre):
        x,y,z=rk4_step(x,y,z,h,0.0,b,r,s,kappa_I,a,c,d,x_R)
    n=int(np.ceil(observation_end_ms/dt_ms))+1
    t=np.empty(n,dtype=np.float64); xs=np.empty(n); ys=np.empty(n); zs=np.empty(n)
    spikes=np.empty(max_spikes,dtype=np.float64); ns=0
    t[0]=0.0; xs[0]=x; ys[0]=y; zs[0]=z
    last=-1e30; xprev=x
    ok=True
    for i in range(1,n):
        drive=J if (i-1)*dt_ms < stimulus_duration_ms else 0.0
        x,y,z=rk4_step(x,y,z,h,drive,b,r,s,kappa_I,a,c,d,x_R)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6:
            ok=False
            return t[:i],xs[:i],ys[:i],zs[:i],spikes[:ns],ok
        ti=min(i*dt_ms,observation_end_ms)
        t[i]=ti; xs[i]=x; ys[i]=y; zs[i]=z
        if xprev < threshold and x >= threshold and ti-last >= refractory_ms:
            den=x-xprev
            frac=(threshold-xprev)/den if den != 0 else 0.0
            tc=(i-1+frac)*dt_ms
            if ns < max_spikes:
                spikes[ns]=tc; ns+=1; last=tc
        xprev=x
    return t,xs,ys,zs,spikes[:ns],ok

def simulate_trace(theta, J, cfg, duration_ms=None, dt_ms=None, observation_end_ms=None):
    m=cfg['model']; a=float(m['a']); c=float(m['c']); d=float(m['d']); x_R=float(m['x_R'])
    if duration_ms is None: duration_ms=float(cfg['analysis']['stimulus_duration_ms'])
    if dt_ms is None: dt_ms=float(cfg['analysis']['dt_ms'])
    if observation_end_ms is None: observation_end_ms=float(duration_ms)
    max_spikes=max(1000,int(observation_end_ms/max(float(m['model_refractory_ms']),dt_ms))+20)
    return simulate_trace_core(float(J),float(duration_ms),float(observation_end_ms),float(m['pre_ms']),float(dt_ms),float(m['model_time_scale_ms']),
        float(theta['b']),float(theta['r']),float(theta['s']),float(theta['kappa_I']),a,c,d,x_R,
        float(m['x0']),float(m['y0']),float(m['z0']),float(m['model_spike_threshold']),float(m['model_refractory_ms']),max_spikes)

def has_spike(theta,J,cfg,dt_ms=None):
    *_,spikes,ok=simulate_trace(theta,J,cfg,dt_ms=dt_ms)
    return bool(ok and len(spikes)>0)

def divergence(x, theta, cfg):
    m=cfg['model']
    return -3.0*float(m['a'])*np.asarray(x)**2 + 2.0*float(theta['b'])*np.asarray(x) - 1.0 - float(theta['r'])

def speed(x,y,z,J,theta,cfg):
    m=cfg['model']; out=[]
    for xi,yi,zi in zip(np.asarray(x),np.asarray(y),np.asarray(z)):
        dx,dy,dz=rhs(float(xi),float(yi),float(zi),float(J),float(theta['b']),float(theta['r']),float(theta['s']),float(theta['kappa_I']),float(m['a']),float(m['c']),float(m['d']),float(m['x_R']))
        out.append((dx*dx+dy*dy+dz*dz)**0.5)
    return np.asarray(out,float)
