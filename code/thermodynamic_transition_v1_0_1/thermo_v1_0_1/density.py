from __future__ import annotations
import numpy as np
from scipy.ndimage import gaussian_filter
from .model import drift_on_grid


def grid_from_samples(samples_by_p,cfg):
    dc=cfg['density']; n=int(dc['bins']); pooled=np.concatenate(samples_by_p,axis=0)
    qlo=float(dc.get('quantile_low',0.001)); qhi=float(dc.get('quantile_high',0.999)); margin=float(dc.get('margin_fraction',0.15))
    lo=np.quantile(pooled,qlo,axis=0); hi=np.quantile(pooled,qhi,axis=0); span=np.maximum(hi-lo,1e-6); lo=lo-margin*span; hi=hi+margin*span
    edges=[np.linspace(lo[i],hi[i],n+1) for i in range(3)]; centers=[0.5*(e[:-1]+e[1:]) for e in edges]
    return lo,hi,edges,centers


def histogram_mass(samples,edges,cfg):
    H,_=np.histogramdd(samples,bins=edges); sigma=float(cfg['density'].get('gaussian_sigma_bins',1.0))
    if sigma>0: H=gaussian_filter(H.astype(float),sigma=sigma,mode='nearest')
    H += float(cfg['density'].get('pseudocount',1e-12)); H /= H.sum(); return H


def js_divergence(p,q,eps=1e-300):
    p=np.asarray(p,float); q=np.asarray(q,float); m=.5*(p+q)
    return .5*np.sum(p*np.log((p+eps)/(m+eps)))+.5*np.sum(q*np.log((q+eps)/(m+eps)))


def build_density_stack(samples_by_p,p_grid,params,J,cfg):
    lo,hi,edges,centers=grid_from_samples(samples_by_p,cfg); masses=[]; split_js=[]
    for s in samples_by_p:
        masses.append(histogram_mass(s,edges,cfg)); a=s[::2]; b=s[1::2]; pa=histogram_mass(a,edges,cfg); pb=histogram_mass(b,edges,cfg); split_js.append(js_divergence(pa,pb))
    masses=np.asarray(masses,float); metrics=stationary_metrics(masses,centers,np.asarray(p_grid,float),params,J,cfg); metrics['split_js']=np.asarray(split_js,float)
    return {'masses':masses,'lo':np.asarray(lo,float),'hi':np.asarray(hi,float),'centers':centers,'metrics':metrics}


def stationary_metrics(masses,centers,p_grid,params,J,cfg):
    floor_rel=float(cfg['density'].get('mass_floor_relative',1e-12)); pm=np.maximum(masses,floor_rel*np.max(masses,axis=(1,2,3),keepdims=True)); pm/=pm.sum(axis=(1,2,3),keepdims=True); logp=np.log(pm)
    Hdisc=-np.sum(pm*logp,axis=(1,2,3)); dx=[float(c[1]-c[0]) for c in centers]; vol=dx[0]*dx[1]*dx[2]; Hdiff=Hdisc+np.log(vol)
    p0=pm[0]; pN=pm[-1]; kl_wt=np.sum(pm*(logp-np.log(np.maximum(p0,floor_rel*np.max(p0)))),axis=(1,2,3)); kl_sca=np.sum(pm*(logp-np.log(np.maximum(pN,floor_rel*np.max(pN)))),axis=(1,2,3)); dlog=np.gradient(logp,p_grid,axis=0,edge_order=1); fisher=np.sum(pm*dlog*dlog,axis=(1,2,3))
    epr_cutoffs=[float(x) for x in cfg['epr']['rho_floor_relative_sensitivity']]; epr=np.empty((len(p_grid),len(epr_cutoffs)),float); retained=np.empty_like(epr); X,Y,Z=np.meshgrid(centers[0],centers[1],centers[2],indexing='ij'); D=np.asarray(cfg['noise']['D'],float)*float(cfg['noise'].get('multiplier',1.0))
    for ip in range(len(p_grid)):
        rho=pm[ip]/vol; grads=np.gradient(rho,dx[0],dx[1],dx[2],edge_order=1); Fx,Fy,Fz=drift_on_grid(X,Y,Z,float(J[ip]),params[ip],cfg); Js=[Fx*rho-D[0]*grads[0],Fy*rho-D[1]*grads[1],Fz*rho-D[2]*grads[2]]
        for ic,cut in enumerate(epr_cutoffs):
            mask=rho>=cut*np.max(rho); retained[ip,ic]=float(np.sum(pm[ip][mask])); denom=np.maximum(rho,1e-300); integ=np.zeros_like(rho)
            for k in range(3): integ += (Js[k]*Js[k])/(D[k]*denom)
            epr[ip,ic]=float(np.sum(integ[mask])*vol)
    primary_cut=float(cfg['epr']['primary_rho_floor_relative']); idx=int(np.argmin(np.abs(np.asarray(epr_cutoffs)-primary_cut)))
    return {'entropy_discrete':Hdisc,'entropy_differential':Hdiff,'kl_to_wt':kl_wt,'kl_to_sca3':kl_sca,'kl_balance':kl_wt-kl_sca,'fisher':fisher,'epr_all':epr,'epr_retained_mass_all':retained,'epr':epr[:,idx],'epr_retained_mass':retained[:,idx],'epr_cutoffs':np.asarray(epr_cutoffs,float)}
