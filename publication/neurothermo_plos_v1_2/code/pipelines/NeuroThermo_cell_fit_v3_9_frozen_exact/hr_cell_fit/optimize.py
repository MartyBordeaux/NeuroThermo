from __future__ import annotations
import time
import numpy as np
from scipy.optimize import differential_evolution
from .params import active_params, theta_to_unit, unit_to_theta
from .objective import objective_unit, evaluate_theta


def baseline_center(cell, cfg):
    return theta_to_unit(cell['baseline_theta'], cfg)


def _de(cell, cfg, dt, tau, bounds, seed, popsize, maxiter, alignment_stage):
    return differential_evolution(
        objective_unit,
        bounds=bounds,
        args=(cell, cfg, float(dt), float(tau), str(alignment_stage)),
        seed=int(seed),
        popsize=int(popsize),
        maxiter=int(maxiter),
        tol=float(cfg['optimization'].get('de_tol', 0.001)),
        atol=0.0,
        polish=False,
        updating='immediate',
        workers=1,
    )


def fit_cell(cell, cfg, seed_offset=0):
    t0 = time.time()
    opt = cfg['optimization']
    npar = len(active_params(cfg))
    seed = int(opt['seed']) + int(seed_offset)
    target_tau = float(cfg['loss']['vp_tau_ms'])
    search_tau = float(opt['search_vp_tau_ms'])
    dt_search = float(opt['dt_search_ms'])
    dt_refine = float(opt['dt_refine_ms'])
    dt_fine = float(opt['dt_fine_ms'])
    center = np.clip(baseline_center(cell, cfg), 0, 1)
    candidates = []

    ev0 = evaluate_theta(unit_to_theta(center, cfg), cell, cfg, dt_refine, target_tau, 'refine')
    candidates.append({
        'label':'frozen_best_seed', 'u':center.copy(), 'search_loss':np.nan,
        'target_coarse_loss':float(ev0.loss), 'target_coarse_spike_loss':float(ev0.spike_train_loss),
        'target_coarse_threshold_loss':float(ev0.threshold_loss), 'nfev':0,
    })

    radius = float(opt['prior_radius_fraction'])
    local_bounds = [(max(0.0, x-radius), min(1.0, x+radius)) for x in center]
    for k in range(int(opt.get('n_prior_starts', 2))):
        res = _de(cell, cfg, dt_search, search_tau, local_bounds, seed+1009*k,
                  opt['prior_de_popsize'], opt['prior_de_maxiter'], 'search')
        u = np.clip(res.x, 0, 1)
        ev = evaluate_theta(unit_to_theta(u, cfg), cell, cfg, dt_refine, target_tau, 'refine')
        candidates.append({
            'label':'prior_de_%d' % (k+1), 'u':u, 'search_loss':float(res.fun),
            'target_coarse_loss':float(ev.loss), 'target_coarse_spike_loss':float(ev.spike_train_loss),
            'target_coarse_threshold_loss':float(ev.threshold_loss), 'nfev':int(res.nfev),
        })

    best_pre = min(candidates, key=lambda c: c['target_coarse_loss'])
    rescue_used = False
    if bool(opt.get('global_rescue_enabled', True)) and (bool(opt.get('global_search_always', False)) or best_pre['target_coarse_loss'] > float(opt['global_rescue_loss_threshold'])):
        rescue_used = True
        for k in range(int(opt.get('n_global_starts', 1))):
            res = _de(cell, cfg, dt_search, search_tau, [(0.0,1.0)]*npar, seed+30011+4099*k,
                      opt['global_de_popsize'], opt['global_de_maxiter'], 'search')
            u = np.clip(res.x, 0, 1)
            ev = evaluate_theta(unit_to_theta(u, cfg), cell, cfg, dt_refine, target_tau, 'refine')
            candidates.append({
                'label':'global_wide_%d' % (k+1), 'u':u, 'search_loss':float(res.fun),
                'target_coarse_loss':float(ev.loss), 'target_coarse_spike_loss':float(ev.spike_train_loss),
                'target_coarse_threshold_loss':float(ev.threshold_loss), 'nfev':int(res.nfev),
            })
        best_pre = min(candidates, key=lambda c: c['target_coarse_loss'])

    rr = float(opt['refine_radius_fraction'])
    c = np.asarray(best_pre['u'])
    refine_bounds = [(max(0.0, x-rr), min(1.0, x+rr)) for x in c]
    res = _de(cell, cfg, dt_refine, target_tau, refine_bounds, seed+60013,
              opt['refine_de_popsize'], opt['refine_de_maxiter'], 'refine')
    u = np.clip(res.x, 0, 1)
    ev = evaluate_theta(unit_to_theta(u, cfg), cell, cfg, dt_refine, target_tau, 'refine')
    candidates.append({
        'label':'target_refine', 'u':u, 'search_loss':float(res.fun),
        'target_coarse_loss':float(ev.loss), 'target_coarse_spike_loss':float(ev.spike_train_loss),
        'target_coarse_threshold_loss':float(ev.threshold_loss), 'nfev':int(res.nfev),
    })

    final = []
    for cnd in candidates:
        th = unit_to_theta(cnd['u'], cfg)
        ev = evaluate_theta(th, cell, cfg, dt_fine, target_tau, 'fine')
        final.append({
            **cnd, 'theta':th, 'final_loss':float(ev.loss), 'spike_train_loss':float(ev.spike_train_loss),
            'threshold_loss':float(ev.threshold_loss), 'sweep_evals':ev.sweep_evals,
            'threshold_eval':ev.threshold_eval, 'ok':bool(ev.ok),
        })
    final.sort(key=lambda c: c['final_loss'])
    best = final[0]
    return {
        'theta':best['theta'], 'u':best['u'], 'loss':best['final_loss'],
        'spike_train_loss':best['spike_train_loss'], 'threshold_loss':best['threshold_loss'],
        'sweep_evals':best['sweep_evals'], 'threshold_eval':best['threshold_eval'],
        'best_label':best['label'], 'rescue_used':rescue_used, 'elapsed_s':float(time.time()-t0),
        'baseline_center_u':center,
        'solutions':[
            {
                'label':c['label'], 'final_loss':c['final_loss'], 'spike_train_loss':c['spike_train_loss'],
                'threshold_loss':c['threshold_loss'], 'target_coarse_loss':c['target_coarse_loss'],
                'target_coarse_spike_loss':c['target_coarse_spike_loss'],
                'target_coarse_threshold_loss':c['target_coarse_threshold_loss'],
                'search_loss':c['search_loss'], 'nfev':c['nfev'], 'u':np.asarray(c['u']).tolist(), **c['theta']
            }
            for c in final
        ],
    }
