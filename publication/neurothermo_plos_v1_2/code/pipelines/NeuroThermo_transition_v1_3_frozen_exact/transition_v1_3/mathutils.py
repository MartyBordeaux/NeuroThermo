from __future__ import annotations
import numpy as np


def persistent_crossing(x, y, thr, persistence=2):
    x=np.asarray(x,float); y=np.asarray(y,float)
    n=len(x)
    for i in range(n):
        if not np.isfinite(y[i]) or y[i] <= thr:
            continue
        j=min(n, i+int(persistence))
        if np.all(np.isfinite(y[i:j])) and np.all(y[i:j] > thr):
            if i == 0:
                return float(x[0])
            if np.isfinite(y[i-1]) and y[i-1] <= thr:
                if y[i] == y[i-1]:
                    return float(x[i])
                return float(x[i-1] + (thr-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]))
            return float(x[i])
    return np.nan


def weighted_censored_quantile(values, weights, q):
    """Right-censored quantile for crossings on [0,1].

    Non-finite crossing values mean "no crossing by 1" and retain their weight as
    right-censored mass. A quantile is finite only if enough total support crossed.
    """
    v=np.asarray(values,float); w=np.asarray(weights,float)
    good_w=np.isfinite(w)&(w>0)
    v=v[good_w]; w=w[good_w]
    if len(w)==0 or not np.isfinite(w.sum()) or w.sum()<=0:
        return np.nan
    total=float(w.sum())
    finite=np.isfinite(v)
    vf=v[finite]; wf=w[finite]
    if len(vf)==0 or wf.sum() + 1e-15 < float(q)*total:
        return np.nan
    o=np.argsort(vf); vf=vf[o]; wf=wf[o]
    c=np.cumsum(wf)
    target=float(q)*total
    # Linear interpolation in cumulative total-weight space; prepend first value at zero.
    xp=np.concatenate([[0.0],c])
    fp=np.concatenate([[vf[0]],vf])
    return float(np.interp(target,xp,fp))


def crossing_support_weight(values, weights):
    v=np.asarray(values,float); w=np.asarray(weights,float)
    good=np.isfinite(w)&(w>0)
    if not np.any(good): return np.nan
    den=float(np.sum(w[good]))
    return float(np.sum(w[good & np.isfinite(v)]) / den) if den>0 else np.nan


def drive_curve(p, family):
    p=np.asarray(p,float)
    if family=='coupled': return p
    if family=='drive_early': return 1.0-(1.0-p)**2
    if family=='drive_late': return p*p
    raise ValueError(family)


def bilinear_track(pis,pds,M,x,y):
    pis=np.asarray(pis,float); pds=np.asarray(pds,float); M=np.asarray(M,float)
    x=np.clip(np.asarray(x,float),pis[0],pis[-1]); y=np.clip(np.asarray(y,float),pds[0],pds[-1])
    i=np.searchsorted(pis,x,side='right')-1; j=np.searchsorted(pds,y,side='right')-1
    i=np.clip(i,0,len(pis)-2); j=np.clip(j,0,len(pds)-2)
    x0=pis[i]; x1=pis[i+1]; y0=pds[j]; y1=pds[j+1]
    q00=M[i,j]; q01=M[i,j+1]; q10=M[i+1,j]; q11=M[i+1,j+1]
    tx=np.divide(x-x0,x1-x0,out=np.zeros_like(x),where=(x1!=x0))
    ty=np.divide(y-y0,y1-y0,out=np.zeros_like(y),where=(y1!=y0))
    out=(1-tx)*(1-ty)*q00+(1-tx)*ty*q01+tx*(1-ty)*q10+tx*ty*q11
    valid=np.isfinite(q00)&np.isfinite(q01)&np.isfinite(q10)&np.isfinite(q11)
    out=np.where(valid,out,np.nan)
    return out


def weighted_finite_quantile(values, weights, q):
    """v1.1-compatible weighted quantile over finite markers only.

    Non-crossing support is reported separately through crossing_support_weight.
    """
    v=np.asarray(values,float); w=np.asarray(weights,float)
    m=np.isfinite(v)&np.isfinite(w)&(w>0)
    v,w=v[m],w[m]
    if len(v)==0:
        return np.nan
    o=np.argsort(v);v,w=v[o],w[o]
    c=np.cumsum(w)/np.sum(w)
    return float(np.interp(float(q),c,v))
