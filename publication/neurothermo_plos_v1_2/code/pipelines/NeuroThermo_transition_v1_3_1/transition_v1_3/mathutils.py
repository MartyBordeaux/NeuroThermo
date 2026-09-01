from __future__ import annotations
import numpy as np

def weighted_quantile(values, weights, q):
    v=np.asarray(values,float); w=np.asarray(weights,float)
    m=np.isfinite(v)&np.isfinite(w)&(w>0)
    v=v[m]; w=w[m]
    if len(v)==0:return np.nan
    o=np.argsort(v); v=v[o]; w=w[o]; cw=np.cumsum(w); t=float(q)*cw[-1]
    return float(v[np.searchsorted(cw,t,side='left')])

def weighted_median(values,weights):
    return weighted_quantile(values,weights,0.5)

def crossing_x(x,y,thr,persistence=2):
    x=np.asarray(x,float); y=np.asarray(y,float)
    for i in range(len(x)):
        if not np.isfinite(y[i]) or y[i]<thr: continue
        ok=True
        for k in range(persistence):
            if i+k>=len(y) or not np.isfinite(y[i+k]) or y[i+k]<thr: ok=False; break
        if not ok:continue
        if i==0:return float(x[0])
        j=i-1
        while j>=0 and not np.isfinite(y[j]):j-=1
        if j<0:return float(x[i])
        if y[j]>=thr:return float(x[i])
        den=y[i]-y[j]
        if abs(den)<1e-15:return float(x[i])
        return float(x[j]+(thr-y[j])*(x[i]-x[j])/den)
    return np.nan
