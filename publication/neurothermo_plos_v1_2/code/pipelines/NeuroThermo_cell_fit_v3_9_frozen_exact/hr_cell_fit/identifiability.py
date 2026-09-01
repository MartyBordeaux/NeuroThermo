from __future__ import annotations
import math
import numpy as np
from scipy.optimize import differential_evolution
from .params import active_params, unit_to_theta
from .objective import objective_unit


def _value_to_unit(value, spec):
    lo=float(spec['min']); hi=float(spec['max']); x=float(value)
    if spec.get('scale','linear')=='log':
        return (math.log(x)-math.log(lo))/(math.log(hi)-math.log(lo))
    return (x-lo)/(hi-lo)


def _separated_physical_limits(best_value, ref_spec, frac):
    lo=float(ref_spec['min']); hi=float(ref_spec['max']); x=float(best_value)
    if ref_spec.get('scale','linear')=='log':
        d=float(frac)*(math.log(hi)-math.log(lo))
        return x*math.exp(-d), x*math.exp(d)
    d=float(frac)*(hi-lo)
    return x-d, x+d


def practical_identifiability(best_theta,best_loss,cell,cfg,seed_offset=0):
    ic=cfg['identifiability']
    if not bool(ic.get('enabled',True)):
        return {'overall':'SKIPPED','reason':'disabled'}
    if len(cell['sweeps']) < int(ic.get('min_spiking_sweeps',2)):
        return {'overall':'INSUFFICIENT_SPIKING_SWEEPS','reason':'fewer_than_min_spiking_sweeps'}
    names=active_params(cfg)
    frac=float(ic.get('reference_separation_fraction',0.15))
    ref=ic.get('reference_bounds') or cfg['bounds']
    dt=float(ic['dt_ms']); vp_tau=float(cfg['loss']['vp_tau_ms'])
    abs_tol=float(ic['near_optimal_absolute_loss']); rel_tol=float(ic['near_optimal_relative_loss'])
    tolerance=max(abs_tol,rel_tol*max(float(best_loss),1e-12)); threshold=float(best_loss)+tolerance
    seed=int(cfg['optimization']['seed'])+900000+int(seed_offset)
    results={}; any_nonid=False; alternatives=[]

    for pi,name in enumerate(names):
        low_lim,high_lim=_separated_physical_limits(best_theta[name],ref[name],frac)
        spec=cfg['bounds'][name]; wlo=float(spec['min']); whi=float(spec['max']); domains=[]
        if low_lim>wlo:
            uhi=min(1.0,max(0.0,_value_to_unit(min(low_lim,whi),spec)))
            if uhi>1e-9:
                b=[(0.0,1.0)]*len(names); b[pi]=(0.0,uhi); domains.append(('LOW',b))
        if high_lim<whi:
            ulo=min(1.0,max(0.0,_value_to_unit(max(high_lim,wlo),spec)))
            if ulo<1.0-1e-9:
                b=[(0.0,1.0)]*len(names); b[pi]=(ulo,1.0); domains.append(('HIGH',b))
        side_results=[]
        for si,(side,bounds) in enumerate(domains):
            # Every alternative recomputes first-spike alignment from scratch.
            res=differential_evolution(
                objective_unit,bounds=bounds,args=(cell,cfg,dt,vp_tau,'identifiability'),
                seed=seed+pi*1009+si*97,popsize=int(ic['de_popsize']),maxiter=int(ic['de_maxiter']),
                tol=float(cfg['optimization'].get('de_tol',0.001)),atol=0.0,polish=False,updating='immediate',workers=1,
            )
            alt_loss=float(res.fun); alt_theta=unit_to_theta(np.clip(res.x,0,1),cfg); near=alt_loss<=threshold
            rec={'side':side,'loss':alt_loss,'near_optimal':bool(near),'theta':alt_theta,'u':np.asarray(res.x).tolist()}
            side_results.append(rec); alternatives.append({'parameter':name,'side':side,'loss':alt_loss,'near_optimal':bool(near),**alt_theta})
        nonid=any(x['near_optimal'] for x in side_results); any_nonid=any_nonid or nonid
        results[name]={
            'status':'NONIDENTIFIABLE' if nonid else 'IDENTIFIABLE',
            'reference_separation_fraction':frac,
            'reference_bounds':ref[name],
            'best_separated_loss':min([x['loss'] for x in side_results],default=np.nan),
            'near_optimal_threshold':threshold,'sides':side_results,
        }
    return {
        'overall':'NONIDENTIFIABLE' if any_nonid else 'IDENTIFIABLE',
        'method':'wide_bound_separated_alternative_reoptimization_with_original_range_separation_and_exact_first_spike_alignment',
        'latency_realigned_for_every_alternative':True,'separation_uses_original_v3_6_reference_range':True,
        'loss_tolerance':tolerance,'near_optimal_threshold':threshold,'parameter_status':results,'alternatives':alternatives,
    }
