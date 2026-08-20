from __future__ import annotations
import numpy as np
from scipy.signal import find_peaks

def weighted_quantile(values,weights,q):
    v=np.asarray(values,float); w=np.asarray(weights,float); m=np.isfinite(v)&np.isfinite(w)&(w>0); v=v[m]; w=w[m]
    if v.size==0:return np.nan
    o=np.argsort(v);v=v[o];w=w[o];cw=np.cumsum(w);cw/=cw[-1];return float(np.interp(float(q),cw,v))
def weighted_mean(values,weights):
    v=np.asarray(values,float);w=np.asarray(weights,float);m=np.isfinite(v)&np.isfinite(w)&(w>0);return float(np.sum(v[m]*w[m])/np.sum(w[m])) if np.any(m) else np.nan
def finite_weight_fraction(values,weights):
    v=np.asarray(values,float);w=np.asarray(weights,float);good=np.isfinite(v)&np.isfinite(w)&(w>0);denom=np.sum(w[np.isfinite(w)&(w>0)]);return float(np.sum(w[good])/denom) if denom>0 else np.nan
def crossing_x(x,y,level=0.0):
    x=np.asarray(x,float);y=np.asarray(y,float)-float(level);m=np.isfinite(x)&np.isfinite(y);x=x[m];y=y[m]
    if len(x)<2:return np.nan
    for i in range(len(x)-1):
        if y[i]==0:return float(x[i])
        if y[i]*y[i+1]<0 or y[i+1]==0:
            den=y[i+1]-y[i];return float(x[i]) if den==0 else float(x[i]+(-y[i])*(x[i+1]-x[i])/den)
    return np.nan
def first_persistent_threshold(x,y,baseline_mask,z=3.0,persistence=2):
    x=np.asarray(x,float);y=np.asarray(y,float);bm=np.asarray(baseline_mask,bool);b=y[bm&np.isfinite(y)]
    if len(b)<2:return np.nan
    med=float(np.median(b));mad=float(np.median(np.abs(b-med)));scale=max(1.4826*mad,1e-12);flag=(y>med+z*scale)&np.isfinite(y);p=max(int(persistence),1)
    for i in range(0,len(flag)-p+1):
        if np.all(flag[i:i+p]):return float(x[i])
    return np.nan
def _moving_average_reflect(y,window):
    y=np.asarray(y,float);w=max(int(window),1)
    if w<=1:return y.copy()
    if w%2==0:w+=1
    pad=w//2;yp=np.pad(y,(pad,pad),mode='edge');return np.convolve(yp,np.ones(w)/w,mode='valid')
def interior_persistent_peak(x,y,exclude_edge_points=2,smooth_window=3,prominence_fraction=0.10,min_width_points=1.0):
    x=np.asarray(x,float);y=np.asarray(y,float);m=np.isfinite(x)&np.isfinite(y)
    if np.sum(m)<5:return (np.nan,np.nan,np.nan)
    yy=y.copy()
    if not np.all(np.isfinite(yy)):
        good=np.isfinite(yy)
        if np.sum(good)<5:return (np.nan,np.nan,np.nan)
        yy=np.interp(np.arange(len(yy)),np.where(good)[0],yy[good])
    sm=_moving_average_reflect(yy,smooth_window);q05,q95=np.quantile(sm,[.05,.95]);scale=max(float(q95-q05),1e-12);prom=max(float(prominence_fraction)*scale,1e-12);peaks,props=find_peaks(sm,prominence=prom,width=max(float(min_width_points),1.0));edge=max(int(exclude_edge_points),1);keep=[j for j,p in enumerate(peaks) if edge<=p<len(x)-edge]
    if not keep:return (np.nan,np.nan,np.nan)
    best=max(keep,key=lambda j:(float(props['prominences'][j]),float(sm[peaks[j]])));idx=int(peaks[best]);return float(x[idx]),float(props['prominences'][best]),float(props['widths'][best])
def expmean_and_ess(values):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    if len(v)==0:return np.nan,np.nan,np.nan
    a=-v;m=np.max(a);w=np.exp(a-m);logm=float(m+np.log(np.mean(w)));mean=float(np.exp(logm)) if logm<700 else np.inf;ess=float((w.sum()**2)/np.sum(w*w)) if np.sum(w*w)>0 else np.nan;return mean,logm,ess
def bootstrap_mean_ci(values,weights=None,n_boot=500,seed=0,alpha=0.05):
    v=np.asarray(values,float);weights=np.ones(len(v),float) if weights is None else weights;w=np.asarray(weights,float);m=np.isfinite(v)&np.isfinite(w)&(w>0);v=v[m];w=w[m]
    if len(v)==0:return (np.nan,np.nan,np.nan)
    w=w/w.sum();rng=np.random.default_rng(int(seed));point=float(np.sum(w*v));reps=np.empty(int(n_boot),float);idx=np.arange(len(v))
    for b in range(int(n_boot)):
        draw=rng.choice(idx,size=len(idx),replace=True,p=w);reps[b]=float(np.mean(v[draw]))
    return point,float(np.quantile(reps,alpha/2)),float(np.quantile(reps,1-alpha/2))
