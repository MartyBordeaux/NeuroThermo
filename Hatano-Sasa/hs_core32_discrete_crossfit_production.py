#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Production-scale discrete Hatano-Sasa analysis for all 32 core-secure WT x SCA3
best-fit cell pairs, using the exact frozen q=0.75 current support from NeuroThermo.

Key design:
- exact q=.75 endpoint current: Jq75 = Jrheo_best + .75*(Jmax_observed-Jrheo_best)
- 31-point coupled path: b,s,J linear; r,kappa log
- stochastic HR: dt=.025 ms, D=(.0025,.01,.00025), burn=2400 ms,
  sample=6000 ms, retained every .5 ms
- train and holdout stationary ensembles are fully seed-disjoint
- train only defines the 4x4x3 partition, stationary pi and 50-ms transition matrices
- independent holdout only audits stationarity and defines endpoint start distributions
- finite-state HS expectation is propagated exactly (no rare-event Monte Carlo)
- deterministic x(t) is checked at all 31 path points for every pair
- local Y is summarized by path edge and around the pair-specific coarse KL balance

This is the best-fit layer. It does NOT replace the 264 retained-solution scenario
uncertainty analysis; it is the first production scaling layer across all 32 pairs.
"""
from __future__ import annotations
import argparse, json, math, os, time, warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from numba import njit, set_num_threads
from scipy.stats import t as student_t

# Frozen numerical settings
A=C=1.0; DPAR=5.0; XR=-1.6
DT=0.025
DX,DY,DZ=0.0025,0.01,0.00025
BURN_MS=12000.0
SAMPLE_MS=6000.0
RETAIN_MS=0.5
N_PATH=31
PGRID=np.linspace(0.,1.,N_PATH)
N_TRAIN_SEEDS=2
N_CAL_SEEDS=1
N_TEST_SEEDS=1
LAG_MS=50.0
TRANS_ALPHA=0.05
PI_FLOOR=1e-12
# previously frozen stationarity tolerances from strict pilot
MEAN_SHIFT_TOL=0.10
LOGVAR_TOL=math.log(1.2)
ACF_TOL=0.10
ACF_LAGS_MS=(50.0,100.0,250.0)
# seed families separated by 100M to guarantee no overlap across families/pairs/positions
TRAIN_BASE=100_000_000
CAL_BASE=300_000_000
TEST_BASE=500_000_000
DT_DET=0.05; PRE_MS=500.0; STIM_MS=1000.0

WT_ANIMAL={
 'WT_03':'WT_150130','WT_04':'WT_150305','WT_08':'WT_DD09','WT_09':'WT_DD09','WT_10':'WT_DD09',
 'WT_13':'WT_DD10','WT_14':'WT_DD10','WT_16':'WT_DD10'}
SCA_ANIMAL={'SCA3_01':'SCA3_DD20','SCA3_04':'SCA3_DD20','SCA3_06':'SCA3_DD24','SCA3_09':'SCA3_DD24'}

@njit(cache=True)
def stationary_trace(b,r,s,kappa,J,seed):
    np.random.seed(seed)
    x,y,z=-1.6,-10.0,1.0
    burn=int(round(BURN_MS/DT)); samp=int(round(SAMPLE_MS/DT)); retain=int(round(RETAIN_MS/DT))
    n=samp//retain
    out=np.empty((n,3),np.float64)
    sx=math.sqrt(2*DX*DT); sy=math.sqrt(2*DY*DT); sz=math.sqrt(2*DZ*DT)
    q=0
    for i in range(burn+samp):
        fx=y-x*x*x+b*x*x-z+kappa*J
        fy=1.0-5.0*x*x-y
        fz=r*(s*(x+1.6)-z)
        x+=fx*DT+sx*np.random.randn(); y+=fy*DT+sy*np.random.randn(); z+=fz*DT+sz*np.random.randn()
        if i>=burn and ((i-burn+1)%retain==0):
            out[q,0]=x; out[q,1]=y; out[q,2]=z; q+=1
    return out

@njit(cache=True)
def rk4_path_rate(path):
    out=np.empty((path.shape[0],5),np.float64)
    for pidx in range(path.shape[0]):
        b,r,s,kappa,J=path[pidx]
        x,y,z=-1.6,-10.0,1.0
        npre=int(round(PRE_MS/DT_DET))
        for _ in range(npre):
            # RK4 at J=0
            k1x=y-x*x*x+b*x*x-z; k1y=1-5*x*x-y; k1z=r*(s*(x+1.6)-z)
            x2=x+.5*DT_DET*k1x; y2=y+.5*DT_DET*k1y; z2=z+.5*DT_DET*k1z
            k2x=y2-x2*x2*x2+b*x2*x2-z2; k2y=1-5*x2*x2-y2; k2z=r*(s*(x2+1.6)-z2)
            x3=x+.5*DT_DET*k2x; y3=y+.5*DT_DET*k2y; z3=z+.5*DT_DET*k2z
            k3x=y3-x3*x3*x3+b*x3*x3-z3; k3y=1-5*x3*x3-y3; k3z=r*(s*(x3+1.6)-z3)
            x4=x+DT_DET*k3x; y4=y+DT_DET*k3y; z4=z+DT_DET*k3z
            k4x=y4-x4*x4*x4+b*x4*x4-z4; k4y=1-5*x4*x4-y4; k4z=r*(s*(x4+1.6)-z4)
            x+=DT_DET*(k1x+2*k2x+2*k3x+k4x)/6; y+=DT_DET*(k1y+2*k2y+2*k3y+k4y)/6; z+=DT_DET*(k1z+2*k2z+2*k3z+k4z)/6
        n=int(round(STIM_MS/DT_DET)); lastx=x; count=0; xmin=x; xmax=x; prev_t=-1e9; sumisi=0.; nisi=0
        for ii in range(n):
            # RK4 at J
            k1x=y-x*x*x+b*x*x-z+kappa*J; k1y=1-5*x*x-y; k1z=r*(s*(x+1.6)-z)
            x2=x+.5*DT_DET*k1x; y2=y+.5*DT_DET*k1y; z2=z+.5*DT_DET*k1z
            k2x=y2-x2*x2*x2+b*x2*x2-z2+kappa*J; k2y=1-5*x2*x2-y2; k2z=r*(s*(x2+1.6)-z2)
            x3=x+.5*DT_DET*k2x; y3=y+.5*DT_DET*k2y; z3=z+.5*DT_DET*k2z
            k3x=y3-x3*x3*x3+b*x3*x3-z3+kappa*J; k3y=1-5*x3*x3-y3; k3z=r*(s*(x3+1.6)-z3)
            x4=x+DT_DET*k3x; y4=y+DT_DET*k3y; z4=z+DT_DET*k3z
            k4x=y4-x4*x4*x4+b*x4*x4-z4+kappa*J; k4y=1-5*x4*x4-y4; k4z=r*(s*(x4+1.6)-z4)
            xn=x+DT_DET*(k1x+2*k2x+2*k3x+k4x)/6; yn=y+DT_DET*(k1y+2*k2y+2*k3y+k4y)/6; zn=z+DT_DET*(k1z+2*k2z+2*k3z+k4z)/6
            if lastx<0 and xn>=0:
                frac=(-lastx)/(xn-lastx); tt=(ii+frac)*DT_DET
                if tt-prev_t>=1.0:
                    if count>0: sumisi+=tt-prev_t; nisi+=1
                    prev_t=tt; count+=1
            x,y,z=xn,yn,zn; lastx=xn
            if x<xmin:xmin=x
            if x>xmax:xmax=x
        out[pidx,0]=count; out[pidx,1]=count/(STIM_MS/1000.0); out[pidx,2]=sumisi/nisi if nisi>0 else np.nan; out[pidx,3]=xmin; out[pidx,4]=xmax
    return out

def seed(base,pair_idx,pidx,sidx):
    return int(base + pair_idx*1_000_000 + pidx*10_000 + sidx*101)

def build_path(w,s):
    p=PGRID; out=np.empty((N_PATH,5),float)
    out[:,0]=(1-p)*w['b']+p*s['b']; out[:,2]=(1-p)*w['s']+p*s['s']; out[:,4]=(1-p)*w['J_q75']+p*s['J_q75']
    out[:,1]=np.exp((1-p)*np.log(w['r'])+p*np.log(s['r'])); out[:,3]=np.exp((1-p)*np.log(w['kappa_I'])+p*np.log(s['kappa_I']))
    return out

def acf(x,lag):
    x=np.asarray(x,float); x=x-x.mean(); den=np.dot(x,x)
    if den<=0 or lag>=len(x):return np.nan
    return float(np.dot(x[:-lag],x[lag:])/den)

def partition_edges(pool):
    ex=np.r_[-np.inf,np.quantile(pool[:,0],[.25,.5,.75]),np.inf]
    ey=np.r_[-np.inf,np.quantile(pool[:,1],[.25,.5,.75]),np.inf]
    ez=np.r_[-np.inf,np.quantile(pool[:,2],[1/3,2/3]),np.inf]
    return ex,ey,ez

def label_states(X,edges):
    ex,ey,ez=edges
    ix=np.searchsorted(ex,X[:,0],side='right')-1; iy=np.searchsorted(ey,X[:,1],side='right')-1; iz=np.searchsorted(ez,X[:,2],side='right')-1
    ix=np.clip(ix,0,3); iy=np.clip(iy,0,3); iz=np.clip(iz,0,2)
    return (ix*12+iy*3+iz).astype(np.int32)

def stationary_pi(T,init=None):
    K=T.shape[0]
    v=np.full(K,1/K) if init is None else np.asarray(init,float)/np.sum(init)
    for _ in range(100000):
        vn=v@T
        if np.max(np.abs(vn-v))<1e-14: v=vn; break
        v=vn
    v=np.maximum(v,PI_FLOOR); return v/v.sum()

def transition_model(label_lists,lag,K=48):
    # Circular lag counts: row and column marginals are identical, so the
    # empirical pi derived from row sums is stationary for T to numerical precision.
    Cnt=np.full((K,K),TRANS_ALPHA,float)
    for lab in label_lists:
        n=len(lab)
        for t in range(n):
            Cnt[lab[t],lab[(t+lag)%n]] += 1.0
    row=Cnt.sum(axis=1); pi=row/row.sum(); T=Cnt/row[:,None]
    return pi,T

def kl(p,q):
    p=np.asarray(p);q=np.asarray(q);return float(np.sum(p*np.log(np.maximum(p,PI_FLOOR)/np.maximum(q,PI_FLOOR))))

def js(p,q):
    m=.5*(p+q);return .5*kl(p,m)+.5*kl(q,m)

def pava(y):
    y=np.asarray(y,float); vals=[]; weights=[]; starts=[]; ends=[]
    for i,v in enumerate(y):
        vals.append(float(v)); weights.append(1.0); starts.append(i); ends.append(i)
        while len(vals)>=2 and vals[-2]>vals[-1]:
            w=weights[-2]+weights[-1]; vv=(vals[-2]*weights[-2]+vals[-1]*weights[-1])/w
            st=starts[-2]; en=ends[-1]
            vals[-2:]=[vv]; weights[-2:]=[w]; starts[-2:]=[st]; ends[-2:]=[en]
    out=np.empty(len(y))
    for v,st,en in zip(vals,starts,ends): out[st:en+1]=v
    return out

def zero_crossing(x,y):
    for i in range(len(y)-1):
        if y[i]==0:return float(x[i])
        if y[i]<0<=y[i+1]:
            if y[i+1]==y[i]:return float(x[i])
            return float(x[i]+(0-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]))
    return np.nan

def exact_hs(pi,Ts,qstart,direction,cvals=None):
    if cvals is None: cvals=np.zeros(N_PATH-1,float)
    if direction=='forward':
        v=np.asarray(qstart,float).copy()
        for i in range(N_PATH-1):
            # y = log(pi_i/pi_i+1) + c_i; exp(-y)=pi_i+1/pi_i * exp(-c_i)
            v=v*(pi[i+1]/pi[i])*math.exp(-float(cvals[i])); v=v@Ts[i+1]
    else:
        v=np.asarray(qstart,float).copy()
        for i in range(N_PATH-1,0,-1):
            j=i-1
            # reverse y = log(pi_i/pi_i-1) - c_j; exp(-y)=pi_i-1/pi_i * exp(+c_j)
            v=v*(pi[i-1]/pi[i])*math.exp(float(cvals[j])); v=v@Ts[i-1]
    return float(v.sum())

def local_y_stats(pi,Ts,qstart,direction,cvals):
    rows=[]; q=np.asarray(qstart,float).copy(); mids=.5*(PGRID[:-1]+PGRID[1:])
    if direction=='forward':
        for i in range(N_PATH-1):
            y=np.log(pi[i]/pi[i+1])+cvals[i]; mu=float(np.sum(q*y)); var=float(np.sum(q*(y-mu)**2)); expm=float(np.sum(q*np.exp(-y)))
            rows.append((i,mids[i],mu,math.sqrt(max(var,0)),expm)); q=q@Ts[i+1]
    else:
        for i in range(N_PATH-1,0,-1):
            j=i-1; y=np.log(pi[i]/pi[i-1])-cvals[j]; mu=float(np.sum(q*y)); var=float(np.sum(q*(y-mu)**2)); expm=float(np.sum(q*np.exp(-y)))
            rows.append((j,mids[j],mu,math.sqrt(max(var,0)),expm)); q=q@Ts[i-1]
        rows.sort(key=lambda z:z[0])
    return rows

def mean_ci_abs(values):
    x=np.asarray(values,float); n=len(x); m=float(np.mean(x)); se=float(np.std(x,ddof=1)/math.sqrt(n)) if n>1 else np.nan
    crit=float(student_t.ppf(.975,n-1)) if n>1 else np.nan
    lo=m-crit*se; hi=m+crit*se
    return m,lo,hi,max(abs(lo),abs(hi))

def pair_worker(arg):
    pair_idx,w,s,outroot=arg; set_num_threads(1); t0=time.time()
    pair=f"{w['cell_id']}__{s['cell_id']}"; out=Path(outroot)/pair; out.mkdir(parents=True,exist_ok=True)
    path=build_path(w,s)
    det=rk4_path_rate(path)
    pd.DataFrame({'p':PGRID,'n_spikes':det[:,0].astype(int),'rate_Hz':det[:,1],'mean_ISI_ms':det[:,2],'x_min':det[:,3],'x_max':det[:,4]}).to_csv(out/'deterministic_path.csv',index=False)
    train_by_p=[]; cal_by_p=[]; test_by_p=[]; stat_rows=[]
    for j,pars in enumerate(path):
        tr=[]
        for si in range(N_TRAIN_SEEDS): tr.append(stationary_trace(*map(float,pars),seed(TRAIN_BASE,pair_idx,j,si)))
        train_by_p.append(tr)
        cc=[]
        for si in range(N_CAL_SEEDS): cc.append(stationary_trace(*map(float,pars),seed(CAL_BASE,pair_idx,j,si)))
        cal_by_p.append(cc)
        tt=[]
        for si in range(N_TEST_SEEDS):
            z=stationary_trace(*map(float,pars),seed(TEST_BASE,pair_idx,j,si)); tt.append(z); h=len(z)//2
            for vi,var in enumerate(('x','y','z')):
                a=z[:h,vi]; b=z[h:,vi]; v1=a.var(ddof=1); v2=b.var(ddof=1); pooled=math.sqrt(max(1e-300,.5*(v1+v2)))
                rr={'p_index':j,'p':PGRID[j],'seed_index':si,'variable':var,'mean_shift_sd':(b.mean()-a.mean())/pooled,'log_var_ratio':math.log(max(v2,1e-300)/max(v1,1e-300))}
                for lm in ACF_LAGS_MS: rr[f'acf_diff_{int(lm)}ms']=acf(b,int(round(lm/RETAIN_MS)))-acf(a,int(round(lm/RETAIN_MS)))
                stat_rows.append(rr)
        test_by_p.append(tt)
    stat=pd.DataFrame(stat_rows); stat.to_csv(out/'stationarity_raw.csv',index=False)
    srows=[]
    for (j,p,var),g in stat.groupby(['p_index','p','variable']):
        row={'p_index':j,'p':p,'variable':var}
        for col,tol in [('mean_shift_sd',MEAN_SHIFT_TOL),('log_var_ratio',LOGVAR_TOL),('acf_diff_50ms',ACF_TOL),('acf_diff_100ms',ACF_TOL),('acf_diff_250ms',ACF_TOL)]:
            vals=g[col].values.astype(float); m=float(np.mean(vals)); edge=float(np.max(np.abs(vals)))
            row[col+'_mean']=m; row[col+'_abs']=edge; row[col+'_pass']=edge<=tol
        row['pass_all']=all(row[c+'_pass'] for c in ['mean_shift_sd','log_var_ratio','acf_diff_50ms','acf_diff_100ms','acf_diff_250ms'])
        srows.append(row)
    ssum=pd.DataFrame(srows); ssum.to_csv(out/'stationarity_summary.csv',index=False)
    pool=np.vstack([x for plist in train_by_p for x in plist])
    edges=partition_edges(pool)
    lag=int(round(LAG_MS/RETAIN_MS)); pis=[]; Ts=[]; calpis=[]; testpis=[]
    for j in range(N_PATH):
        labs=[label_states(x,edges) for x in train_by_p[j]]
        pi,T=transition_model(labs,lag); pis.append(pi); Ts.append(T)
        cl=np.concatenate([label_states(x,edges) for x in cal_by_p[j]])
        cp=np.bincount(cl,minlength=48).astype(float)+TRANS_ALPHA; cp/=cp.sum(); calpis.append(cp)
        tl=np.concatenate([label_states(x,edges) for x in test_by_p[j]])
        tp=np.bincount(tl,minlength=48).astype(float)+TRANS_ALPHA; tp/=tp.sum(); testpis.append(tp)
    pis=np.asarray(pis); Ts=np.asarray(Ts); calpis=np.asarray(calpis); testpis=np.asarray(testpis)
    # Cross-fit scalar normalization learned only on calibration data.
    cvals=np.empty(N_PATH-1,float); qcrow=[]
    for j in range(N_PATH-1):
        lr=np.log(pis[j]/pis[j+1])
        Af=float(np.sum(calpis[j]*np.exp(-lr))); Ar=float(np.sum(calpis[j+1]*np.exp(lr)))
        c=.5*(math.log(Af)-math.log(Ar)); cvals[j]=c
        tf=float(np.sum(testpis[j]*np.exp(-(lr+c)))); tr=float(np.sum(testpis[j+1]*np.exp(lr+c)))
        qcrow.append({'edge':j,'p0':PGRID[j],'p1':PGRID[j+1],'cal_c':c,'cal_raw_logA_forward':math.log(Af),'cal_raw_logA_reverse':math.log(Ar),'test_cal_logA_forward':math.log(tf),'test_cal_logA_reverse':math.log(tr)})
    pd.DataFrame(qcrow).to_csv(out/'crossfit_ratio_qc.csv',index=False)
    qF=testpis[0]; qR=testpis[-1]
    Zf_raw=exact_hs(pis,Ts,qF,'forward',None); Zr_raw=exact_hs(pis,Ts,qR,'reverse',None)
    Zf=exact_hs(pis,Ts,qF,'forward',cvals); Zr=exact_hs(pis,Ts,qR,'reverse',cvals)
    Zf_exact=exact_hs(pis,Ts,pis[0],'forward',None); Zr_exact=exact_hs(pis,Ts,pis[-1],'reverse',None)
    loc=[]
    for direction,q in [('forward',qF),('reverse',qR)]:
        for edge,pm,mu,sdv,expm in local_y_stats(pis,Ts,q,direction,cvals):
            loc.append({'pair':pair,'wt_cell':w['cell_id'],'sca3_cell':s['cell_id'],'direction':direction,'edge':edge,'p_mid':pm,'mean_y':mu,'sd_y':sdv,'local_expmean':expm,'mean_abs_proxy':abs(mu)})
    locdf=pd.DataFrame(loc); locdf.to_csv(out/'local_Y_by_edge.csv',index=False)
    # coarse endpoint-relative KL marker from stationary pi
    dwt=np.array([kl(p,pis[0]) for p in pis]); dsc=np.array([kl(p,pis[-1]) for p in pis]); delta=dwt-dsc; iso=pava(delta); pkl=zero_crossing(PGRID,iso)
    pd.DataFrame({'p':PGRID,'KL_to_WT':dwt,'KL_to_SCA3':dsc,'delta':delta,'delta_isotonic':iso}).to_csv(out/'coarse_KL_curve.csv',index=False)
    # untouched test-vs-train stationary distribution audit
    tv=.5*np.sum(np.abs(testpis-pis),axis=1); jsv=np.array([js(testpis[i]+PI_FLOOR, pis[i]+PI_FLOOR) for i in range(N_PATH)])
    stat_res=np.max(np.abs(np.einsum('pk,pkj->pj',pis,Ts)-pis),axis=1)
    pd.DataFrame({'p':PGRID,'test_TV':tv,'test_JS':jsv,'pi_stationarity_residual':stat_res}).to_csv(out/'discrete_stationary_distribution_qc.csv',index=False)
    # Y region summaries
    regrows=[]
    for direction,g in locdf.groupby('direction'):
        for name,mask in [('early',g.p_mid<.30),('global_KL_window',(g.p_mid>=.35)&(g.p_mid<=.55)),('late',g.p_mid>.70),('pair_KL_window',np.abs(g.p_mid-pkl)<=.10 if np.isfinite(pkl) else np.zeros(len(g),bool))]:
            gg=g.loc[mask]
            regrows.append({'direction':direction,'region':name,'n_edges':len(gg),'mean_of_mean_y':gg.mean_y.mean() if len(gg) else np.nan,'mean_abs_mean_y':gg.mean_y.abs().mean() if len(gg) else np.nan,'mean_sd_y':gg.sd_y.mean() if len(gg) else np.nan})
    pd.DataFrame(regrows).to_csv(out/'Y_region_summary.csv',index=False)
    # summary
    zss=ssum[ssum.variable=='z']; fwd=locdf[locdf.direction=='forward']; rev=locdf[locdf.direction=='reverse']
    sf={
        'pair':pair,'wt_cell':w['cell_id'],'sca3_cell':s['cell_id'],'wt_animal':WT_ANIMAL[w['cell_id']],'sca3_animal_day':SCA_ANIMAL[s['cell_id']],
        'wt_J_q75':w['J_q75'],'sca3_J_q75':s['J_q75'],
        'det_all_spiking':bool((det[:,0]>0).all()),'det_n_nonspiking_points':int(np.sum(det[:,0]==0)),'det_min_rate_Hz':float(np.min(det[:,1])),'det_max_rate_Hz':float(np.max(det[:,1])),
        'stationarity_all_pass_fraction':float(ssum.pass_all.mean()),'stationarity_z_pass_fraction':float(zss.pass_all.mean()),
        'stationarity_z_max_mean_abs':float(zss.mean_shift_sd_abs.max()),'stationarity_z_max_logvar_abs':float(zss.log_var_ratio_abs.max()),
        'stationarity_z_max_acf_abs':float(zss[[f'acf_diff_{int(x)}ms_abs' for x in ACF_LAGS_MS]].max(axis=1).max()),
        'holdout_TV_median':float(np.median(tv)),'holdout_TV_max':float(np.max(tv)),'holdout_JS_median':float(np.median(jsv)),
        'pi_stationarity_residual_max':float(np.max(stat_res)),
        'HS_exact_forward':Zf_exact,'HS_exact_reverse':Zr_exact,'HS_independent_start_forward_raw':Zf_raw,'HS_independent_start_reverse_raw':Zr_raw,
        'HS_independent_start_forward':Zf,'HS_independent_start_reverse':Zr,'log_HS_independent_start_forward':math.log(Zf),'log_HS_independent_start_reverse':math.log(Zr),
        'test_sum_logA_forward':float(pd.DataFrame(qcrow).test_cal_logA_forward.sum()),'test_sum_logA_reverse':float(pd.DataFrame(qcrow).test_cal_logA_reverse.sum()),
        'test_max_abs_logA':float(max(pd.DataFrame(qcrow).test_cal_logA_forward.abs().max(),pd.DataFrame(qcrow).test_cal_logA_reverse.abs().max())),
        'coarse_KL_balance_p':pkl,
        'forward_meanY':float(fwd.mean_y.sum()),'reverse_meanY':float(rev.mean_y.sum()),
        'forward_peak_abs_meanY_p':float(fwd.loc[fwd.mean_y.abs().idxmax(),'p_mid']),'reverse_peak_abs_meanY_p':float(rev.loc[rev.mean_y.abs().idxmax(),'p_mid']),
        'runtime_s':time.time()-t0}
    pd.DataFrame([sf]).to_csv(out/'pair_summary.csv',index=False)
    return sf

def qsummary(x):
    a=np.asarray(x,float); return {'median':float(np.nanmedian(a)),'q25':float(np.nanquantile(a,.25)),'q75':float(np.nanquantile(a,.75)),'min':float(np.nanmin(a)),'max':float(np.nanmax(a))}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--endpoints',default=str(Path(__file__).with_name('core_endpoints_best_q75.csv'))); ap.add_argument('--outdir',default=str(Path(__file__).with_name('results_core32_discrete_q75'))); ap.add_argument('--workers',type=int,default=8); args=ap.parse_args()
    ep=pd.read_csv(args.endpoints); ep['J_q75']=ep.rheobase_J+.75*(ep.J_max_observed-ep.rheobase_J)
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True); ep.to_csv(out/'core_endpoints_exact_q75.csv',index=False)
    cfg={'analysis':'discrete finite-state Hatano-Sasa production best-fit layer','n_pairs':32,'n_p':N_PATH,'dt_ms':DT,'D':[DX,DY,DZ],'burn_ms':BURN_MS,'sample_ms':SAMPLE_MS,'retain_ms':RETAIN_MS,
         'train_seeds_per_p':N_TRAIN_SEEDS,'calibration_seeds_per_p':N_CAL_SEEDS,'test_seeds_per_p':N_TEST_SEEDS,'partition':'4x4x3 pooled-train quantiles','lag_dwell_ms':LAG_MS,'transition_pseudocount':TRANS_ALPHA,
         'stationarity_tolerances':{'mean_shift_sd':MEAN_SHIFT_TOL,'abs_log_var_ratio':LOGVAR_TOL,'acf_difference':ACF_TOL},'stationarity_note':'single untouched 6-s test trace per p after 12-s burn; direct first-vs-second 3-s halves; calibration is separate','q':.75,'best_fit_only':True}
    (out/'config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    # warmup numba
    _=stationary_trace(3,.01,1,.5,1,1); _=rk4_path_rate(np.tile(np.array([3,.01,1,.5,1.]),(2,1)))
    wt=ep[ep.group=='WT'].sort_values('cell_id'); sc=ep[ep.group=='SCA3'].sort_values('cell_id'); jobs=[]; idx=0
    for _,w in wt.iterrows():
        for _,s in sc.iterrows():
            pair=f"{w.cell_id}__{s.cell_id}"
            if (out/pair/'pair_summary.csv').exists():
                idx+=1; continue
            jobs.append((idx,w.to_dict(),s.to_dict(),str(out))); idx+=1
    rows=[]
    # preserve already completed pairs when resuming after an interrupted long run
    for f in out.glob('*__*/pair_summary.csv'):
        try: rows.append(pd.read_csv(f).iloc[0].to_dict())
        except Exception: pass
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        fut={ex.submit(pair_worker,j):j[0] for j in jobs}
        for f in as_completed(fut):
            r=f.result(); rows.append(r); print(f"[{len(rows):02d}/32] {r['pair']} Zf={r['HS_independent_start_forward']:.4f} Zr={r['HS_independent_start_reverse']:.4f} pKL={r['coarse_KL_balance_p']:.3f} statZ={r['stationarity_z_pass_fraction']:.2f}",flush=True)
    ps=pd.DataFrame(rows).sort_values(['wt_cell','sca3_cell']); ps.to_csv(out/'all_pair_summary.csv',index=False)
    # collect local/region tables
    loc=[]; regs=[]; stat=[]
    for pair in ps.pair:
        loc.append(pd.read_csv(out/pair/'local_Y_by_edge.csv')); x=pd.read_csv(out/pair/'Y_region_summary.csv'); x['pair']=pair; regs.append(x)
        x=pd.read_csv(out/pair/'stationarity_summary.csv'); x['pair']=pair; stat.append(x)
    loc=pd.concat(loc,ignore_index=True); regs=pd.concat(regs,ignore_index=True); stat=pd.concat(stat,ignore_index=True)
    loc.to_csv(out/'all_local_Y_by_edge.csv',index=False); regs.to_csv(out/'all_Y_region_summary.csv',index=False); stat.to_csv(out/'all_stationarity_summary.csv',index=False)
    ens=loc.groupby(['direction','edge','p_mid']).agg(median_mean_y=('mean_y','median'),q25_mean_y=('mean_y',lambda x:np.quantile(x,.25)),q75_mean_y=('mean_y',lambda x:np.quantile(x,.75)),median_sd_y=('sd_y','median'),median_local_expmean=('local_expmean','median')).reset_index()
    ens.to_csv(out/'ensemble_local_Y_by_edge.csv',index=False)
    # cell-balanced summary: each endpoint cell average across partners first
    cb=[]
    for direction in ['forward','reverse']:
        col='log_HS_independent_start_forward' if direction=='forward' else 'log_HS_independent_start_reverse'
        for c in sorted(set(ps.wt_cell)|set(ps.sca3_cell)):
            g=ps[(ps.wt_cell==c)|(ps.sca3_cell==c)]; cb.append({'direction':direction,'cell':c,'mean_log_HS':float(g[col].mean()),'median_log_HS':float(g[col].median()),'n_partners':len(g)})
    cb=pd.DataFrame(cb); cb.to_csv(out/'cell_balanced_HS.csv',index=False)
    # 4 WT animals x2 SCA3 animal-days =8 crossed strata
    strata=[]
    for (wa,sa),g in ps.groupby(['wt_animal','sca3_animal_day']):
        strata.append({'wt_animal':wa,'sca3_animal_day':sa,'n_pairs':len(g),'forward_mean_log_HS':float(g.log_HS_independent_start_forward.mean()),'reverse_mean_log_HS':float(g.log_HS_independent_start_reverse.mean()),'median_pKL':float(g.coarse_KL_balance_p.median())})
    pd.DataFrame(strata).to_csv(out/'animal_day_strata_summary.csv',index=False)
    # overall summary json
    summ={
      'n_pairs':32,
      'det_all_spiking_pairs':int(ps.det_all_spiking.sum()),
      'stationarity_all_points_all_vars_pass_pairs':int((ps.stationarity_all_pass_fraction==1).sum()),
      'stationarity_all_z_points_pass_pairs':int((ps.stationarity_z_pass_fraction==1).sum()),
      'HS_exact_forward_abs_error_max':float(np.max(np.abs(ps.HS_exact_forward-1))),
      'HS_exact_reverse_abs_error_max':float(np.max(np.abs(ps.HS_exact_reverse-1))),
      'log_HS_independent_start_forward':qsummary(ps.log_HS_independent_start_forward),
      'log_HS_independent_start_reverse':qsummary(ps.log_HS_independent_start_reverse),
      'pairs_abs_log_HS_le_0p1_forward':int((ps.log_HS_independent_start_forward.abs()<=.1).sum()),
      'pairs_abs_log_HS_le_0p1_reverse':int((ps.log_HS_independent_start_reverse.abs()<=.1).sum()),
      'coarse_KL_balance_p':qsummary(ps.coarse_KL_balance_p),
      'holdout_TV_median_across_pairs':float(ps.holdout_TV_median.median()),
      'best_fit_only':True,
      'interpretation':'finite-state coarse-grained HS; retained-fit alternatives not yet propagated'}
    # local Y region aggregate, pair-level values
    regagg={}
    for (d,rn),g in regs.groupby(['direction','region']): regagg[f'{d}:{rn}']={k:v for k,v in qsummary(g.mean_abs_mean_y).items()}
    summ['mean_abs_local_Y_by_region']=regagg
    (out/'RUN_SUMMARY.json').write_text(json.dumps(summ,indent=2),encoding='utf-8')
    print(json.dumps(summ,indent=2),flush=True)

if __name__=='__main__': main()
