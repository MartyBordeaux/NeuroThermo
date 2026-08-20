from __future__ import annotations
import numpy as np
from numba import njit

def _assign(X,C):
    out=np.empty(len(X),dtype=np.int32);chunk=5000
    for a in range(0,len(X),chunk):
        q=X[a:a+chunk];d=((q[:,None,:]-C[None,:,:])**2).sum(axis=2);out[a:a+chunk]=np.argmin(d,axis=1)
    return out
def common_partition(samples_by_p,cfg,seed):
    mc=cfg['markov'];K=int(mc['n_states']);rng=np.random.default_rng(int(seed));pooled=np.concatenate(samples_by_p,axis=0).astype(float);mu=pooled.mean(axis=0);sd=pooled.std(axis=0);sd=np.where(sd>1e-9,sd,1.0);Z=(pooled-mu)/sd;maxn=int(mc.get('kmeans_max_samples',20000));Zfit=Z[rng.choice(len(Z),size=maxn,replace=False)] if len(Z)>maxn else Z
    if len(Zfit)<K:K=max(4,len(Zfit)//2)
    C=Zfit[rng.choice(len(Zfit),size=K,replace=False)].copy()
    for _ in range(int(mc.get('kmeans_iterations',20))):
        lab=_assign(Zfit,C);new=C.copy()
        for k in range(K):
            q=Zfit[lab==k];new[k]=q.mean(axis=0) if len(q) else Zfit[rng.integers(0,len(Zfit))]
        if np.max(np.abs(new-C))<1e-5:C=new;break
        C=new
    return [_assign((s-mu)/sd,C) for s in samples_by_p],C,mu,sd
def transition_models(labels_by_p,dwell_ms,cfg):
    K=max(int(np.max(x)) for x in labels_by_p)+1;stride=float(cfg['stationary']['sample_stride_ms']);lag=max(1,int(round(float(dwell_ms)/stride)));alpha=float(cfg['markov'].get('transition_pseudocount',0.05));pis=[];Ts=[];stationarity=[];eprs=[]
    for lab in labels_by_p:
        n=len(lab);C=np.full((K,K),alpha,dtype=float)
        for t in range(n):C[lab[t],lab[(t+lag)%n]]+=1.0
        row=C.sum(axis=1);pi=row/row.sum();T=C/row[:,None];pis.append(pi);Ts.append(T);stationarity.append(float(np.max(np.abs(pi@T-pi))));flux=pi[:,None]*T;rev=flux.T;eprs.append(float(np.sum(flux*np.log(flux/rev)))/max(float(dwell_ms),1e-12))
    return np.asarray(pis),np.asarray(Ts),np.asarray(stationarity),np.asarray(eprs),lag
def exact_hatano_sasa(pis,Ts,direction):
    n=len(pis)
    if direction=='forward':
        v=pis[0].copy()
        for i in range(n-1):v=v*(pis[i+1]/pis[i]);v=v@Ts[i+1]
    else:
        v=pis[-1].copy()
        for i in range(n-1,0,-1):v=v*(pis[i-1]/pis[i]);v=v@Ts[i-1]
    return float(v.sum())
@njit(cache=True)
def _sample_cat(cdf,u):
    lo=0;hi=len(cdf)-1
    while lo<hi:
        mid=(lo+hi)//2
        if u<=cdf[mid]:hi=mid
        else:lo=mid+1
    return lo
@njit(cache=True)
def math_log(x):return np.log(max(x,1e-300))
@njit(cache=True)
def _simulate_discrete(seed,direction,ntraj,pis,Ts,cpi,cT):
    np.random.seed(seed);nlam=pis.shape[0];Y=np.empty(ntraj);S=np.empty(ntraj)
    for q in range(ntraj):
        states=np.empty(nlam,np.int64)
        if direction==1:
            s=_sample_cat(cpi[0],np.random.rand());states[0]=s;y=0.0;lf=math_log(pis[0,s])
            for i in range(nlam-1):y+=math_log(pis[i,s])-math_log(pis[i+1,s]);sn=_sample_cat(cT[i+1,s],np.random.rand());lf+=math_log(Ts[i+1,s,sn]);s=sn;states[i+1]=s
            lr=math_log(pis[-1,states[-1]])
            for i in range(nlam-2,-1,-1):lr+=math_log(Ts[i,states[i+1],states[i]])
            Y[q]=y;S[q]=lf-lr
        else:
            s=_sample_cat(cpi[-1],np.random.rand());states[-1]=s;y=0.0;lr=math_log(pis[-1,s])
            for i in range(nlam-1,0,-1):y+=math_log(pis[i,s])-math_log(pis[i-1,s]);sn=_sample_cat(cT[i-1,s],np.random.rand());lr+=math_log(Ts[i-1,s,sn]);s=sn;states[i-1]=s
            lf=math_log(pis[0,states[0]])
            for i in range(nlam-1):lf+=math_log(Ts[i+1,states[i],states[i+1]])
            Y[q]=y;S[q]=lr-lf
    return Y,S
def simulate_markov_fluctuations(seed,pis,Ts,cfg):
    n=int(cfg['markov']['n_trajectories']);pis=pis.astype(float);Ts=Ts.astype(float);cpi=np.cumsum(pis,axis=1);cT=np.cumsum(Ts,axis=2);yf,sf=_simulate_discrete(int(seed),1,n,pis,Ts,cpi,cT);yr,sr=_simulate_discrete(int(seed)+900001,-1,n,pis,Ts,cpi,cT)
    return {'Y_forward':yf,'Sigma_forward':sf,'Y_reverse':yr,'Sigma_reverse':sr,'HS_exact_forward':exact_hatano_sasa(pis,Ts,'forward'),'HS_exact_reverse':exact_hatano_sasa(pis,Ts,'reverse')}
