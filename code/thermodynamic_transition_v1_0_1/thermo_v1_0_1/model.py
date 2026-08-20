from __future__ import annotations
import math
import numpy as np
from numba import njit

@njit(cache=True)
def drift_vec(x,y,z,J,b,r,s,k,a,c,d,xR):
    return (y-a*x*x*x+b*x*x-z+k*J,c-d*x*x-y,r*(s*(x-xR)-z))

@njit(cache=True)
def em_step(x,y,z,h,J,b,r,s,k,a,c,d,xR,Dx,Dy,Dz):
    fx,fy,fz=drift_vec(x,y,z,J,b,r,s,k,a,c,d,xR)
    return (x+fx*h+math.sqrt(2.0*Dx*h)*np.random.randn(),y+fy*h+math.sqrt(2.0*Dy*h)*np.random.randn(),z+fz*h+math.sqrt(2.0*Dz*h)*np.random.randn())

@njit(cache=True)
def stationary_samples_core(seed,x0,y0,z0,burn_steps,sample_steps,stride,h,J,b,r,s,k,a,c,d,xR,Dx,Dy,Dz):
    np.random.seed(seed); x,y,z=x0,y0,z0
    for _ in range(burn_steps):
        x,y,z=em_step(x,y,z,h,J,b,r,s,k,a,c,d,xR,Dx,Dy,Dz)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6: return np.empty((0,3),np.float64),x,y,z,False
    nsave=sample_steps//stride; out=np.empty((nsave,3),np.float64); j=0
    for i in range(sample_steps):
        x,y,z=em_step(x,y,z,h,J,b,r,s,k,a,c,d,xR,Dx,Dy,Dz)
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or abs(x)>1e6 or abs(y)>1e6 or abs(z)>1e6: return out[:j],x,y,z,False
        if (i+1)%stride==0 and j<nsave: out[j,0]=x; out[j,1]=y; out[j,2]=z; j+=1
    return out[:j],x,y,z,True


def stationary_samples(seed,theta,J,cfg,start=None):
    m=cfg['model']; s=cfg['stationary']; noise=cfg['noise']; dt=float(s['dt_ms']); h=dt/float(m['model_time_scale_ms'])
    if start is None: start=(float(m['x0']),float(m['y0']),float(m['z0']))
    burn_steps=int(round(float(s['burn_ms'])/dt)); sample_steps=int(round(float(s['sample_ms'])/dt)); stride=max(1,int(round(float(s['sample_stride_ms'])/dt)))
    D=np.asarray(noise['D'],float)*float(noise.get('multiplier',1.0))
    arr,x,y,z,ok=stationary_samples_core(int(seed),*map(float,start),burn_steps,sample_steps,stride,h,float(J),float(theta[0]),float(theta[1]),float(theta[2]),float(theta[3]),float(m['a']),float(m['c']),float(m['d']),float(m['x_R']),float(D[0]),float(D[1]),float(D[2]))
    return arr,(x,y,z),bool(ok)


def drift_on_grid(X,Y,Z,J,theta,cfg):
    m=cfg['model']; b,r,s,k=map(float,theta); a=float(m['a']); c=float(m['c']); d=float(m['d']); xr=float(m['x_R'])
    return Y-a*X**3+b*X**2-Z+k*float(J), c-d*X**2-Y, r*(s*(X-xr)-Z)
