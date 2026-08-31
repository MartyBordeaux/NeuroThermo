from __future__ import annotations
import math
import numpy as np
import pandas as pd


def mad_scale(x,eps=1e-12):
    a=np.asarray(x,float); a=a[np.isfinite(a)]
    med=float(np.median(a)); mad=float(np.median(np.abs(a-med)))
    if mad<eps:
        q75,q25=np.quantile(a,[.75,.25]); mad=float((q75-q25)/1.349) if q75>q25 else 1.0
    return med,mad


def robust_transform_fit(df,coords):
    rows=[]
    for c in coords:
        center,scale=mad_scale(df[c])
        rows.append({'coordinate':c,'center':center,'scale':scale,'method':'pooled best-cell median / MAD'})
    return pd.DataFrame(rows)


def apply_transform(df,transform,coords,prefix='z_'):
    out=df.copy()
    t=transform.set_index('coordinate')
    for c in coords:
        out[prefix+c]=(pd.to_numeric(out[c],errors='coerce')-float(t.loc[c,'center']))/float(t.loc[c,'scale'])
    return out


def projection_reference(best_cells,coords,secure_mask=None,wt_exit_q=.90,sca_entry_q=.10):
    d=best_cells.copy()
    if secure_mask is not None:
        ref=d[secure_mask(d)].copy()
        if ref.group.nunique()<2 or (ref.group=='WT').sum()<2 or (ref.group=='SCA3').sum()<2:
            ref=d.copy(); subset='all_best_cells_fallback'
        else: subset='core_secure_best_cells'
    else:
        ref=d.copy(); subset='all_best_cells'
    transform=robust_transform_fit(ref,coords)
    rz=apply_transform(ref,transform,coords)
    zcols=['z_'+c for c in coords]
    cwt=rz[rz.group.eq('WT')][zcols].mean().to_numpy(float)
    csc=rz[rz.group.eq('SCA3')][zcols].mean().to_numpy(float)
    delta=csc-cwt; den=float(np.dot(delta,delta))
    if den<=0: raise ValueError('degenerate endpoint centroids')
    z=rz[zcols].to_numpy(float)
    A=((z-cwt)@delta)/den
    rz=rz.assign(A=A)
    wtA=rz.loc[rz.group.eq('WT'),'A'].to_numpy(float)
    scA=rz.loc[rz.group.eq('SCA3'),'A'].to_numpy(float)
    wt_thr=float(np.quantile(wtA,wt_exit_q)); sc_thr=float(np.quantile(scA,sca_entry_q))
    # Orthogonal corridor from reference endpoints.
    foot=cwt[None,:]+A[:,None]*delta[None,:]
    orth=np.linalg.norm(z-foot,axis=1)
    corridor=float(np.quantile(orth,.90)) if len(orth) else np.nan
    return {
        'subset':subset,'transform':transform,'ref_cells':rz,
        'cwt':cwt,'csc':csc,'delta':delta,'den':den,
        'wt_exit_A_threshold':wt_thr,'sca3_entry_A_threshold':sc_thr,
        'cloud_overlap':bool(wt_thr>=sc_thr),
        'corridor_radius_q90':corridor,
        'centroid_distance':float(np.sqrt(den)),
    }


def project_points(df,ref,coords):
    z=apply_transform(df,ref['transform'],coords)
    zcols=['z_'+c for c in coords]
    arr=z[zcols].to_numpy(float)
    good=np.all(np.isfinite(arr),axis=1)
    A=np.full(len(z),np.nan); O=np.full(len(z),np.nan)
    if good.any():
        aa=((arr[good]-ref['cwt'])@ref['delta'])/ref['den']
        foot=ref['cwt'][None,:]+aa[:,None]*ref['delta'][None,:]
        A[good]=aa; O[good]=np.linalg.norm(arr[good]-foot,axis=1)
    z['A']=A; z['orthogonal_distance']=O
    return z


def _interp_cross(x,y,thr):
    x=np.asarray(x,float); y=np.asarray(y,float)
    for i in range(1,len(x)):
        if not (np.isfinite(y[i-1]) and np.isfinite(y[i])): continue
        if y[i-1] < thr <= y[i]:
            if y[i]==y[i-1]: return float(x[i])
            return float(x[i-1]+(thr-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]))
    return np.nan


def persistent_crossing(x,y,thr,persistence=2):
    x=np.asarray(x,float); y=np.asarray(y,float)
    n=len(x)
    for i in range(n):
        if not np.isfinite(y[i]) or y[i]<=thr: continue
        j=min(n,i+persistence)
        # Near the endpoint there may be fewer than `persistence` future grid
        # points. If every remaining point is beyond the threshold, accept the
        # crossing rather than making SCA3-entry impossible at p=1.
        if np.all(np.isfinite(y[i:j])) and np.all(y[i:j]>thr):
            if i==0: return float(x[0])
            if np.isfinite(y[i-1]) and y[i-1]<=thr:
                if y[i]==y[i-1]: return float(x[i])
                return float(x[i-1]+(thr-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]))
            return float(x[i])
    return np.nan


def count_reversals(y,tol=1e-4):
    y=np.asarray(y,float); y=y[np.isfinite(y)]
    if len(y)<3: return 0
    dy=np.diff(y); dy=dy[np.abs(dy)>tol]
    if len(dy)<2: return 0
    s=np.sign(dy)
    return int(np.sum(s[1:]!=s[:-1]))


def inverse_morph(p,name):
    p=float(np.clip(p,0,1))
    if name=='linear': return p
    if name=='quadratic': return math.sqrt(p)
    if name=='sqrt': return p*p
    if name=='smoothstep':
        lo,hi=0.0,1.0
        for _ in range(60):
            m=.5*(lo+hi); v=3*m*m-2*m*m*m
            if v<p: lo=m
            else: hi=m
        return .5*(lo+hi)
    raise ValueError(name)
