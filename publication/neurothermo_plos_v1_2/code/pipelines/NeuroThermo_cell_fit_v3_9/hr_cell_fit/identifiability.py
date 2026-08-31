from __future__ import annotations
import numpy as np
from scipy.optimize import differential_evolution
from .params import z_to_params, params_to_z, pack_params
from .objective import evaluate_cell

PARAMS=('b','r','s','kappa_I')

def _ref_coord(value,bound):
    lo=float(bound['min']);hi=float(bound['max'])
    if str(bound.get('scale','linear')).lower()=='log':return (np.log(value)-np.log(lo))/(np.log(hi)-np.log(lo))
    return (value-lo)/(hi-lo)

def _parameter_reference_distance(a,b,p,ref_bounds):return abs(_ref_coord(a[p],ref_bounds[p])-_ref_coord(b[p],ref_bounds[p]))


def profile_identifiability(cell,best_z,best_eval,cfg,threshold_bracket=None):
    icfg=cfg['identifiability']
    if len(cell['sweeps'])<int(icfg['min_spiking_sweeps']):return {'cell_id':cell['cell_id'],'identifiability':'INSUFFICIENT_SPIKING_SWEEPS'}
    best_p=best_eval['params'];best_loss=best_eval['loss'];sep=float(icfg['reference_separation_fraction']);ref=icfg['reference_bounds']
    alt_limit=best_loss+max(float(icfg['near_optimal_absolute_loss']),float(icfg['near_optimal_relative_loss'])*max(best_loss,1e-12))
    rng=np.random.default_rng(abs(hash(cell['cell_id']))%(2**31-1)); rows=[]; per_param={}
    for pi,p in enumerate(PARAMS):
        def constrained_obj(z):
            pp=z_to_params(z,cfg['bounds'])
            d=_parameter_reference_distance(pp,best_p,p,ref)
            if d<sep:return 1e4+1e3*(sep-d)
            ev=evaluate_cell(cell,z,cfg,dt_ms=icfg['dt_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket,identifiability_mode=True)
            return ev['loss']
        candidates=[]
        for side in (-1,1):
            seed=np.array(best_z,float);seed[pi]=np.clip(seed[pi]+side*sep,0,1);candidates.append(seed)
        bounds=[(0,1)]*4
        res=differential_evolution(constrained_obj,bounds,popsize=int(icfg['de_popsize']),maxiter=int(icfg['de_maxiter']),tol=cfg['optimization']['de_tol'],seed=int(rng.integers(0,2**31-1)),workers=1,polish=True,updating='immediate')
        candidates.append(np.asarray(res.x,float));evaluated=[]
        for z in candidates:
            pp=z_to_params(z,cfg['bounds']);dist=_parameter_reference_distance(pp,best_p,p,ref)
            if dist+1e-9<sep:continue
            ev=evaluate_cell(cell,z,cfg,dt_ms=icfg['dt_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=threshold_bracket,identifiability_mode=True)
            evaluated.append((ev['loss'],dist,z,ev))
        if not evaluated:per_param[p]='SEARCH_FAILED';continue
        loss,dist,z,ev=min(evaluated,key=lambda x:x[0]);alt_p=ev['params'];near=loss<=alt_limit;per_param[p]='NON_IDENTIFIABLE' if near else 'IDENTIFIABLE'
        row={'cell_id':cell['cell_id'],'parameter_tested':p,'best_loss':best_loss,'alt_loss':loss,'near_optimal_limit':alt_limit,'reference_parameter_distance':dist,'reference_separation_required':sep,'parameter_status':per_param[p],'latency_alignment_method':'exact_first_spike','latency_realigned_for_every_alternative':True}
        for k in PARAMS:row['best_'+k]=best_p[k];row['alt_'+k]=alt_p[k]
        row['best_median_abs_latency_shift_ms']=float(np.nanmedian([abs(s.get('latency_shift_ms',np.nan)) for s in best_eval['sweeps']]))
        row['alt_median_abs_latency_shift_ms']=float(np.nanmedian([abs(s.get('latency_shift_ms',np.nan)) for s in ev['sweeps']]))
        rows.append(row)
    final='IDENTIFIABLE' if all(per_param.get(p)=='IDENTIFIABLE' for p in PARAMS) else 'NON_IDENTIFIABLE'
    return {'cell_id':cell['cell_id'],'identifiability':final,'per_parameter':per_param,'alternative_rows':rows}
