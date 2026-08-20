from __future__ import annotations
import math
import numpy as np
from numba import njit

@njit(cache=True)
def _drift(x,y,z,J,b,r,s,k,a,c,d,xR):
    return (y-a*x*x*x+b*x*x-z+k*J,c-d*x*x-y,r*(s*(x-xR)-z))
@njit(cache=True)
def _path_values(p,params,J):
    b=(1-p)*params[0,0]+p*params[-1,0];s=(1-p)*params[0,2]+p*params[-1,2];j=(1-p)*J[0]+p*J[-1];r=math.exp((1-p)*math.log(max(params[0,1],1e-300))+p*math.log(max(params[-1,1],1e-300)));k=math.exp((1-p)*math.log(max(params[0,3],1e-300))+p*math.log(max(params[-1,3],1e-300)));return b,r,s,k,j
@njit(cache=True)
def _interp_logmass(x,y,z,grid,lo,hi):
    nx,ny,nz=grid.shape;oob=x<lo[0] or x>hi[0] or y<lo[1] or y>hi[1] or z<lo[2] or z>hi[2]
    def coord(v,l,h,n):
        t=(v-l)/(h-l)*n-0.5;t=0.0 if t<0 else (n-1.0 if t>n-1 else t);i0=int(math.floor(t));i1=min(i0+1,n-1);return i0,i1,t-i0
    i0,i1,fx=coord(x,lo[0],hi[0],nx);j0,j1,fy=coord(y,lo[1],hi[1],ny);k0,k1,fz=coord(z,lo[2],hi[2],nz)
    c000=grid[i0,j0,k0];c100=grid[i1,j0,k0];c010=grid[i0,j1,k0];c110=grid[i1,j1,k0];c001=grid[i0,j0,k1];c101=grid[i1,j0,k1];c011=grid[i0,j1,k1];c111=grid[i1,j1,k1];c00=c000*(1-fx)+c100*fx;c10=c010*(1-fx)+c110*fx;c01=c001*(1-fx)+c101*fx;c11=c011*(1-fx)+c111*fx;c0=c00*(1-fy)+c10*fy;c1=c01*(1-fy)+c11*fy;return c0*(1-fz)+c1*fz,oob
@njit(cache=True)
def _hs_one(seed,direction,params,J,logmass,lo,hi,D,h,dwell_steps,endpoint_samples,a,c,d,xR):
    np.random.seed(seed);nlam=params.shape[0];nseg=nlam-1;pick=np.random.randint(0,endpoint_samples.shape[0]);x=endpoint_samples[pick,0];y=endpoint_samples[pick,1];z=endpoint_samples[pick,2];Y=0.0;oob=0;current=0 if direction==1 else nlam-1
    for seg in range(nseg):
        target=seg+1 if direction==1 else nlam-2-seg;lm0,ob0=_interp_logmass(x,y,z,logmass[current],lo,hi);lm1,ob1=_interp_logmass(x,y,z,logmass[target],lo,hi);Y+=lm0-lm1;oob+=int(ob0)+int(ob1);b,r,s,k=params[target];jj=J[target]
        for _ in range(dwell_steps):
            fx,fy,fz=_drift(x,y,z,jj,b,r,s,k,a,c,d,xR);x=x+fx*h+math.sqrt(2*D[0]*h)*np.random.randn();y=y+fy*h+math.sqrt(2*D[1]*h)*np.random.randn();z=z+fz*h+math.sqrt(2*D[2]*h)*np.random.randn()
        current=target
    return Y,float(oob)/(2.0*nseg)
@njit(cache=True)
def _log_kernel_cont(x0,y0,z0,x1,y1,z1,p,params,J,h,D,a,c,d,xR):
    b,r,s,k,jj=_path_values(p,params,J);fx,fy,fz=_drift(x0,y0,z0,jj,b,r,s,k,a,c,d,xR);rx=x1-x0-fx*h;ry=y1-y0-fy*h;rz=z1-z0-fz*h;return -0.25/h*(rx*rx/D[0]+ry*ry/D[1]+rz*rz/D[2])
@njit(cache=True)
def _sigma_one(seed,direction,params,J,logmass,lo,hi,D,h,total_steps,endpoint_samples,a,c,d,xR):
    np.random.seed(seed);states=np.empty((total_steps+1,3),np.float64);pick=np.random.randint(0,endpoint_samples.shape[0]);states[0]=endpoint_samples[pick];x,y,z=states[0,0],states[0,1],states[0,2]
    for n in range(total_steps):
        q=(n+0.5)/total_steps;p=q if direction==1 else 1.0-q;b,r,s,k,jj=_path_values(p,params,J);fx,fy,fz=_drift(x,y,z,jj,b,r,s,k,a,c,d,xR);x=x+fx*h+math.sqrt(2*D[0]*h)*np.random.randn();y=y+fy*h+math.sqrt(2*D[1]*h)*np.random.randn();z=z+fz*h+math.sqrt(2*D[2]*h)*np.random.randn();states[n+1]=[x,y,z]
    if direction==1:
        l0,ob0=_interp_logmass(*states[0],logmass[0],lo,hi);lN,obN=_interp_logmass(*states[-1],logmass[-1],lo,hi);log_num=l0;log_den=lN
        for n in range(total_steps):
            p=(n+0.5)/total_steps;log_num+=_log_kernel_cont(*states[n],*states[n+1],p,params,J,h,D,a,c,d,xR);rn=total_steps-1-n;log_den+=_log_kernel_cont(*states[rn+1],*states[rn],p,params,J,h,D,a,c,d,xR)
    else:
        lN,ob0=_interp_logmass(*states[0],logmass[-1],lo,hi);l0,obN=_interp_logmass(*states[-1],logmass[0],lo,hi);log_num=lN;log_den=l0
        for n in range(total_steps):
            p=1.0-(n+0.5)/total_steps;log_num+=_log_kernel_cont(*states[n],*states[n+1],p,params,J,h,D,a,c,d,xR);rn=total_steps-1-n;log_den+=_log_kernel_cont(*states[rn+1],*states[rn],p,params,J,h,D,a,c,d,xR)
    return log_num-log_den,0.5*(int(ob0)+int(obN))
@njit(cache=True)
def _ensemble(seed0,direction,ntraj,params,J,logmass,lo,hi,D,h,dwell_steps,endpoint_samples,a,c,d,xR):
    Y=np.empty(ntraj,np.float64);S=np.empty(ntraj,np.float64);O=np.empty(ntraj,np.float64);total_steps=(params.shape[0]-1)*dwell_steps
    for i in range(ntraj):Y[i],oh=_hs_one(seed0+7919*i,direction,params,J,logmass,lo,hi,D,h,dwell_steps,endpoint_samples,a,c,d,xR);S[i],os=_sigma_one(seed0+2000003+7919*i,direction,params,J,logmass,lo,hi,D,h,total_steps,endpoint_samples,a,c,d,xR);O[i]=max(oh,os)
    return Y,S,O
def simulate_protocols(seed,params,J,masses,lo,hi,endpoint_start,endpoint_end,cfg,dwell_ms):
    m=cfg['model'];pc=cfg['protocol'];noise=cfg['noise'];dt=float(pc['dt_ms']);h=dt/float(m['model_time_scale_ms']);dwell=max(1,int(round(float(dwell_ms)/dt)));floor=float(cfg['density'].get('mass_floor_relative',1e-12));mm=np.maximum(masses,floor*np.max(masses,axis=(1,2,3),keepdims=True));mm/=mm.sum(axis=(1,2,3),keepdims=True);logmass=np.log(mm);D=np.asarray(noise['D'],float)*float(noise.get('multiplier',1.0));n=int(pc['n_trajectories']);const=(float(m['a']),float(m['c']),float(m['d']),float(m['x_R']));YF,SF,OF=_ensemble(int(seed),1,n,params.astype(float),J.astype(float),logmass.astype(float),np.asarray(lo,float),np.asarray(hi,float),D,float(h),int(dwell),np.asarray(endpoint_start,float),*const);YR,SR,OR=_ensemble(int(seed)+1000003,-1,n,params.astype(float),J.astype(float),logmass.astype(float),np.asarray(lo,float),np.asarray(hi,float),D,float(h),int(dwell),np.asarray(endpoint_end,float),*const);return {'Y_forward':YF,'Sigma_forward':SF,'oob_forward':OF,'Y_reverse':YR,'Sigma_reverse':SR,'oob_reverse':OR}
@njit(cache=True)
def _systematic_resample(weights):
    n=len(weights);out=np.empty(n,np.int64);c=np.cumsum(weights);u0=np.random.rand()/n;j=0
    for i in range(n):
        u=u0+i/n
        while j<n-1 and u>c[j]:j+=1
        out[i]=j
    return out
@njit(cache=True)
def _hs_smc_once(seed,direction,params,J,logmass,lo,hi,D,h,dwell_steps,endpoint_samples,nparticles,a,c,d,xR):
    np.random.seed(seed);nlam=params.shape[0];nseg=nlam-1;P=np.empty((nparticles,3),np.float64)
    for q in range(nparticles):P[q]=endpoint_samples[np.random.randint(0,endpoint_samples.shape[0])]
    logZ=0.0;miness=1.0;sumess=0.0;oob=0;current=0 if direction==1 else nlam-1
    for seg in range(nseg):
        target=seg+1 if direction==1 else nlam-2-seg;lw=np.empty(nparticles,np.float64)
        for q in range(nparticles):
            l0,ob0=_interp_logmass(*P[q],logmass[current],lo,hi);l1,ob1=_interp_logmass(*P[q],logmass[target],lo,hi);lw[q]=l1-l0;oob+=int(ob0)+int(ob1)
        ma=np.max(lw);w=np.exp(lw-ma);meanw=np.mean(w)*math.exp(ma)
        if not np.isfinite(meanw) or meanw<=0:return np.nan,0.0,0.0,1.0
        logZ+=math.log(meanw);w=w/np.sum(w);ess=1.0/(np.sum(w*w)*nparticles);miness=min(miness,ess);sumess+=ess;idx=_systematic_resample(w);Q=np.empty_like(P)
        for q in range(nparticles):Q[q]=P[idx[q]]
        P=Q;b,r,s,k=params[target];jj=J[target]
        for _ in range(dwell_steps):
            for q in range(nparticles):
                x,y,z=P[q];fx,fy,fz=_drift(x,y,z,jj,b,r,s,k,a,c,d,xR);P[q,0]=x+fx*h+math.sqrt(2*D[0]*h)*np.random.randn();P[q,1]=y+fy*h+math.sqrt(2*D[1]*h)*np.random.randn();P[q,2]=z+fz*h+math.sqrt(2*D[2]*h)*np.random.randn()
        current=target
    return logZ,miness,sumess/nseg,float(oob)/(2.0*nseg*nparticles)
def simulate_hs_smc(seed,params,J,masses,lo,hi,endpoint_start,endpoint_end,cfg,dwell_ms):
    m=cfg['model'];pc=cfg['protocol'];noise=cfg['noise'];sc=pc['hs_smc'];dt=float(pc['dt_ms']);h=dt/float(m['model_time_scale_ms']);dwell=max(1,int(round(float(dwell_ms)/dt)));floor=float(cfg['density'].get('mass_floor_relative',1e-12));mm=np.maximum(masses,floor*np.max(masses,axis=(1,2,3),keepdims=True));mm/=mm.sum(axis=(1,2,3),keepdims=True);logmass=np.log(mm);D=np.asarray(noise['D'],float)*float(noise.get('multiplier',1.0));npart=int(sc['particles']);reps=int(sc['replicates']);const=(float(m['a']),float(m['c']),float(m['d']),float(m['x_R']));out={}
    for direction,name,ep,off in [(1,'forward',endpoint_start,0),(-1,'reverse',endpoint_end,10000019)]:
        logz=np.empty(reps);mine=np.empty(reps);meane=np.empty(reps);oob=np.empty(reps)
        for rix in range(reps):logz[rix],mine[rix],meane[rix],oob[rix]=_hs_smc_once(int(seed)+off+rix*65537,direction,params.astype(float),J.astype(float),logmass.astype(float),np.asarray(lo,float),np.asarray(hi,float),D,float(h),int(dwell),np.asarray(ep,float),npart,*const)
        out[f'logZ_{name}']=logz;out[f'minESS_{name}']=mine;out[f'meanESS_{name}']=meane;out[f'oob_{name}']=oob
    return out
