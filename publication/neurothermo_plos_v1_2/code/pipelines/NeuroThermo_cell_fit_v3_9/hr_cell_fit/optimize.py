from __future__ import annotations
import time
import numpy as np
from scipy.optimize import differential_evolution
from .params import pack_params,params_to_z,z_to_params
from .objective import evaluate_cell


def _box(center,radius):
    center=np.asarray(center,float);lo=np.maximum(0.0,center-radius);hi=np.minimum(1.0,center+radius);return list(zip(lo,hi))

def _de(obj,bounds,popsize,maxiter,tol,seed,workers):
    updating='deferred' if workers!=1 else 'immediate'
    return differential_evolution(obj,bounds,popsize=popsize,maxiter=maxiter,tol=tol,seed=seed,polish=True,workers=workers,updating=updating)

def _candidate_record(z,ev,source):return {'z':np.asarray(z,float),'eval':ev,'source':source}


def fit_cell(cell,seed_row,cfg,threshold_bracket=None):
    ocfg=cfg['optimization'];bounds_cfg=cfg['bounds'];rng=np.random.default_rng(int(ocfg['seed'])+abs(hash(cell['cell_id']))%1000000);candidates=[];t0=time.time()
    seed_params={p:float(seed_row[p]) for p in ('b','r','s','kappa_I')};seed_z=params_to_z(seed_params,bounds_cfg)
    seed_ev=evaluate_cell(cell,seed_z,cfg,dt_ms=ocfg['dt_refine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket)
    candidates.append(_candidate_record(seed_z,seed_ev,'v3_8_seed'))
    def obj_search(z):return evaluate_cell(cell,z,cfg,dt_ms=ocfg['dt_search_ms'],search_tau_ms=ocfg['search_vp_tau_ms'],threshold_bracket=threshold_bracket)['loss']
    for i in range(int(ocfg['n_prior_starts'])):
        res=_de(obj_search,_box(seed_z,float(ocfg['prior_radius_fraction'])),int(ocfg['prior_de_popsize']),int(ocfg['prior_de_maxiter']),float(ocfg['de_tol']),int(rng.integers(0,2**31-1)),int(ocfg['n_jobs']))
        ev=evaluate_cell(cell,res.x,cfg,dt_ms=ocfg['dt_refine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket);candidates.append(_candidate_record(res.x,ev,f'prior_{i}'))
    best=min(candidates,key=lambda x:x['eval']['loss'])
    if ocfg.get('global_search_always',False) or (ocfg['global_rescue_enabled'] and best['eval']['loss']>ocfg['global_rescue_loss_threshold']):
        for i in range(int(ocfg['n_global_starts'])):
            res=_de(obj_search,[(0,1)]*4,int(ocfg['global_de_popsize']),int(ocfg['global_de_maxiter']),float(ocfg['de_tol']),int(rng.integers(0,2**31-1)),int(ocfg['n_jobs']))
            ev=evaluate_cell(cell,res.x,cfg,dt_ms=ocfg['dt_refine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket);candidates.append(_candidate_record(res.x,ev,f'global_{i}'))
    best=min(candidates,key=lambda x:x['eval']['loss']);ref_center=best['z']
    def obj_refine(z):return evaluate_cell(cell,z,cfg,dt_ms=ocfg['dt_refine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket)['loss']
    res=_de(obj_refine,_box(ref_center,float(ocfg['refine_radius_fraction'])),int(ocfg['refine_de_popsize']),int(ocfg['refine_de_maxiter']),float(ocfg['de_tol']),int(rng.integers(0,2**31-1)),int(ocfg['n_jobs']))
    ev=evaluate_cell(cell,res.x,cfg,dt_ms=ocfg['dt_fine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket);candidates.append(_candidate_record(res.x,ev,'final_refine'))
    fine=[]
    for c in sorted(candidates,key=lambda x:x['eval']['loss'])[:5]:
        fev=evaluate_cell(cell,c['z'],cfg,dt_ms=ocfg['dt_fine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket);fine.append(_candidate_record(c['z'],fev,c['source']+'_fine'))
    best=min(fine,key=lambda x:x['eval']['loss']);return {'best_z':best['z'],'best_eval':best['eval'],'source':best['source'],'elapsed_s':time.time()-t0,'all_candidates':candidates}
