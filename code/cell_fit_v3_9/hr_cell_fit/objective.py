from __future__ import annotations
import math
import numpy as np
from .params import z_to_params
from .model import simulate_spikes
from .latency import exact_first_spike_align


def _count_spikes_at_current(cell,params,cfg,current_pA,stim_ms,post_ms=0.0,dt_ms=None):
    if dt_ms is None:dt_ms=cfg['optimization']['dt_refine_ms']
    J=float(current_pA)/float(cell['capacitance_pF']);sp,ok=simulate_spikes(J,float(stim_ms),float(post_ms),float(dt_ms),params,cfg['model'])
    return int(len(sp)) if ok else None


def evaluate_threshold_constraint(cell,params,cfg,threshold_bracket,dt_ms=None):
    if not cfg['threshold_constraint']['enabled'] or threshold_bracket is None:
        return {'threshold_loss':0.0,'threshold_pass':True,'model_spikes_at_nonspiking_current':np.nan,'model_spikes_at_first_spiking_current':np.nan}
    stim_ms=float(threshold_bracket['stimulus_duration_ms']);i0=float(threshold_bracket['nonspiking_current_pA']);i1=float(threshold_bracket['first_spiking_current_pA'])
    n0=_count_spikes_at_current(cell,params,cfg,i0,stim_ms,0.0,dt_ms);n1=_count_spikes_at_current(cell,params,cfg,i1,stim_ms,0.0,dt_ms)
    if n0 is None or n1 is None:return {'threshold_loss':cfg['loss']['simulation_failure_loss'],'threshold_pass':False,'model_spikes_at_nonspiking_current':-1,'model_spikes_at_first_spiking_current':-1}
    q=cfg['threshold_constraint'];loss=(q['nonspiking_violation_penalty'] if n0>0 else 0.0)+(q['first_spiking_violation_penalty'] if n1==0 else 0.0)
    return {'threshold_loss':float(loss),'threshold_pass':bool(n0==0 and n1>0),'model_spikes_at_nonspiking_current':n0,'model_spikes_at_first_spiking_current':n1}


def evaluate_cell(cell,z,cfg,dt_ms=None,search_tau_ms=None,threshold_bracket=None,identifiability_mode=False):
    if dt_ms is None:dt_ms=cfg['optimization']['dt_refine_ms']
    params=z_to_params(z,cfg['bounds']);sw=[];losses=[];align_cfg=cfg['loss']['latency_alignment']
    if not align_cfg.get('enabled',True) or str(align_cfg.get('method',''))!='exact_first_spike':raise ValueError('v3.9 requires exact_first_spike latency alignment')
    vp_tau=float(search_tau_ms if search_tau_ms is not None else cfg['loss']['vp_tau_ms'])
    for s in cell['sweeps']:
        sim_end=float(s['fit_end_ms']);stim_ms=min(float(s['stimulus_duration_ms']),sim_end);post_ms=max(0.0,sim_end-stim_ms)
        model_spikes,ok=simulate_spikes(float(s['J']),stim_ms,post_ms,float(dt_ms),params,cfg['model'])
        if not ok:
            return {'loss':cfg['loss']['simulation_failure_loss'],'spike_train_loss':cfg['loss']['simulation_failure_loss'],'threshold_loss':0.0,'params':params,'sweeps':[],'failed':True}
        exp=s['exp_spike_times_ms'];ali=exact_first_spike_align(exp,model_spikes,sim_end,tau_ms=vp_tau,normalize=cfg['loss']['normalize'],count_penalty_weight=cfg['loss']['count_penalty_weight'])
        loss=float(ali['loss']) if np.isfinite(ali['loss']) else cfg['loss']['simulation_failure_loss']
        losses.append(loss);sw.append({'sweep':s,'model_spikes':ali['aligned_model_spike_times_ms'],'raw_model_spikes':ali['raw_model_spike_times_ms'],'latency_shift_ms':ali['tau_ms'],'vp_loss':ali['vp_loss'],'raw_vp_loss':ali['raw_vp_loss'],'count_penalty':ali['count_penalty'],'loss':loss,'aligned_last_spike_error_ms':ali['aligned_last_spike_error_ms']})
    spike_loss=float(np.mean(losses)) if cfg['loss']['sweep_weighting']=='equal' else float(np.sum(losses))
    thr=evaluate_threshold_constraint(cell,params,cfg,threshold_bracket,dt_ms=dt_ms)
    return {'loss':spike_loss+thr['threshold_loss'],'spike_train_loss':spike_loss,'threshold_loss':thr['threshold_loss'],'threshold_pass':thr['threshold_pass'],
            'model_spikes_at_nonspiking_current':thr['model_spikes_at_nonspiking_current'],'model_spikes_at_first_spiking_current':thr['model_spikes_at_first_spiking_current'],
            'params':params,'sweeps':sw,'failed':False}
