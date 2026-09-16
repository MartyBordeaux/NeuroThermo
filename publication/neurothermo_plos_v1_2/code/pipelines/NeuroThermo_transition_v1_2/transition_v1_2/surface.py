from __future__ import annotations
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, math, os
import numpy as np
import pandas as pd

from .data import load_inputs, build_support_with_anchors, enumerate_scenarios, frozen_manifest
from .model import pre_relax, refine_rheobase, supported_metrics
from .geometry import load_geometry, project_scalar, stage_from_A, weighted_quantile, persistent_crossing


def _interp_state(sc, p_intrinsic, p_drive):
    pi = float(p_intrinsic)
    pdri = float(p_drive)
    b = (1-pi)*float(sc['wt_b']) + pi*float(sc['sca_b'])
    s = (1-pi)*float(sc['wt_s']) + pi*float(sc['sca_s'])
    r = math.exp((1-pi)*math.log(float(sc['wt_r'])) + pi*math.log(float(sc['sca_r'])))
    k = math.exp((1-pdri)*math.log(float(sc['wt_kappa_I'])) + pdri*math.log(float(sc['sca_kappa_I'])))
    J = (1-pdri)*float(sc['wt_J_q75']) + pdri*float(sc['sca_J_q75'])
    # Preserve the v1.0 coupled-path convention: the active support window follows intrinsic progress.
    window = (1-pi)*float(sc['wt_active_support_ms']) + pi*float(sc['sca_active_support_ms'])
    return {'b': b, 'r': r, 's': s, 'kappa_I': k}, float(J), float(window)


def _safe_log10(x):
    return math.log10(float(x)) if np.isfinite(x) and float(x) > 0 else np.nan


def _checkpoint_path(out, scenario_id):
    return out / 'checkpoints' / f'scenario_{int(scenario_id):04d}.csv.gz'


def _write_checkpoint(path, df):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + '.tmp')
    df.to_csv(tmp, index=False, compression='gzip')
    os.replace(tmp, path)


def _read_checkpoint(path):
    return pd.read_csv(path, compression='gzip')


def _scenario_surface_job(sc, cfg, geom):
    ni = int(cfg['surface']['n_intrinsic'])
    nd = int(cfg['surface']['n_drive'])
    pis = np.linspace(0.0, 1.0, ni)
    pds = np.linspace(0.0, 1.0, nd)
    rows = []
    prev_row_rb = np.full(nd, np.nan)

    for ii, pi in enumerate(pis):
        # At J=0, kappa_I does not affect the pre-relaxed state, so one pre-relaxation
        # is sufficient for the whole drive row at fixed intrinsic progress.
        theta0, _, window = _interp_state(sc, pi, 0.0)
        x, y, z, ok_pre = pre_relax(theta0, cfg, dt_ms=float(cfg['simulation']['dt_ms']))
        pre_state = (x, y, z) if ok_pre else None
        current_row_rb = np.full(nd, np.nan)
        previous_drive_rb = np.nan

        for jj, pdri in enumerate(pds):
            theta, J, window = _interp_state(sc, pi, pdri)
            if ok_pre:
                if np.isfinite(previous_drive_rb):
                    guess = previous_drive_rb
                elif np.isfinite(prev_row_rb[jj]):
                    guess = prev_row_rb[jj]
                else:
                    guess = (1-pi)*float(sc['wt_rheobase_J_endpoint']) + pi*float(sc['sca_rheobase_J_endpoint'])
                rb = refine_rheobase(theta, cfg, guess=guess, pre_state=pre_state)
                if np.isfinite(rb['rheobase_J']):
                    previous_drive_rb = float(rb['rheobase_J'])
                    current_row_rb[jj] = previous_drive_rb
                met = supported_metrics(theta, J, window, cfg, pre_state=pre_state)
            else:
                rb = {'rheobase_J': np.nan, 'status': 'PRE_RELAX_FAIL', 'iterations': 0}
                met = {'spike_count': 0, 'support_rate_hz': np.nan, 'mean_isi_ms': np.nan,
                       'occupancy_fraction': np.nan, 'first_spike_ms': np.nan, 'simulation_ok': False}

            log_rb = _safe_log10(rb['rheobase_J'])
            log_isi = _safe_log10(met['mean_isi_ms'])
            log_active = _safe_log10(met['support_rate_hz'])
            A_isi, orth_isi = project_scalar(log_rb, log_isi, geom['isi'])
            A_active, orth_active = project_scalar(log_rb, log_active, geom['active'])
            stage_isi = stage_from_A(A_isi, geom['isi'])
            stage_active = stage_from_A(A_active, geom['active'])

            rows.append({
                'scenario_id': int(sc['scenario_id']),
                'biological_pair_key': sc['biological_pair_key'],
                'wt_cell_id': sc['wt_cell_id'], 'sca_cell_id': sc['sca_cell_id'],
                'wt_solution_key': sc['wt_solution_key'], 'sca_solution_key': sc['sca_solution_key'],
                'intrinsic_index': ii, 'drive_index': jj,
                'p_intrinsic': float(pi), 'p_drive': float(pdri),
                'b': theta['b'], 'r': theta['r'], 's': theta['s'], 'kappa_I': theta['kappa_I'],
                'J_protocol': J, 'active_support_ms': window,
                'rheobase_J': rb['rheobase_J'], 'rheobase_status': rb['status'],
                'rheobase_iterations': rb['iterations'],
                'J_over_rheobase': float(J/rb['rheobase_J']) if np.isfinite(rb['rheobase_J']) and rb['rheobase_J'] > 0 else np.nan,
                'spike_count': met['spike_count'],
                'active_support_rate_hz': met['support_rate_hz'],
                'mean_isi_ms': met['mean_isi_ms'],
                'occupancy_fraction': met['occupancy_fraction'],
                'first_spike_ms': met['first_spike_ms'],
                'simulation_ok': met['simulation_ok'],
                'A_isi': A_isi, 'orth_isi': orth_isi,
                'stage_isi': stage_isi,
                'inside_isi_corridor': bool(np.isfinite(orth_isi) and orth_isi <= geom['isi']['corridor_radius_q90']),
                'A_active': A_active, 'orth_active': orth_active,
                'stage_active': stage_active,
                'inside_active_corridor': bool(np.isfinite(orth_active) and orth_active <= geom['active']['corridor_radius_q90']),
                'within_pair_support_weight': float(sc['within_pair_support_weight']),
                'scenario_weight': float(sc['scenario_weight']),
                'biological_pair_weight': float(sc['biological_pair_weight']),
            })
        prev_row_rb = current_row_rb
    return pd.DataFrame(rows)


def _pair_surface_for_pair(pair_scenarios, out, secure_cells):
    dfs = []
    for sid in pair_scenarios.scenario_id:
        dfs.append(_read_checkpoint(_checkpoint_path(out, sid)))
    d = pd.concat(dfs, ignore_index=True)
    rows = []
    numeric_metrics = ['A_isi','orth_isi','A_active','orth_active','rheobase_J',
                       'active_support_rate_hz','mean_isi_ms','occupancy_fraction','J_over_rheobase']
    for (ii, jj, pi, pdri), g in d.groupby(['intrinsic_index','drive_index','p_intrinsic','p_drive'], sort=True):
        w = pd.to_numeric(g.within_pair_support_weight, errors='coerce').to_numpy(float)
        sw = np.nansum(w)
        row = {
            'biological_pair_key': g.biological_pair_key.iloc[0],
            'wt_cell_id': g.wt_cell_id.iloc[0], 'sca_cell_id': g.sca_cell_id.iloc[0],
            'both_core_secure': bool(str(g.wt_cell_id.iloc[0]) in secure_cells and str(g.sca_cell_id.iloc[0]) in secure_cells),
            'intrinsic_index': int(ii), 'drive_index': int(jj),
            'p_intrinsic': float(pi), 'p_drive': float(pdri),
            'n_support_scenarios': int(len(g)),
        }
        for m in numeric_metrics:
            row[m + '_weighted_median'] = weighted_quantile(g[m], w, 0.5)
        valid_isi = np.isfinite(pd.to_numeric(g.A_isi, errors='coerce').to_numpy(float))
        valid_active = np.isfinite(pd.to_numeric(g.A_active, errors='coerce').to_numpy(float))
        row['valid_isi_weight'] = float(np.nansum(w[valid_isi]) / sw) if sw > 0 else np.nan
        row['valid_active_weight'] = float(np.nansum(w[valid_active]) / sw) if sw > 0 else np.nan
        for stage, col in [('WT_like','p_WT_like'),('TRANSITION','p_TRANSITION'),('SCA3_like','p_SCA3_like')]:
            mask = g.stage_isi.astype(str).eq(stage).to_numpy()
            row[col + '_isi'] = float(np.nansum(w[mask]) / sw) if sw > 0 else np.nan
            maska = g.stage_active.astype(str).eq(stage).to_numpy()
            row[col + '_active'] = float(np.nansum(w[maska]) / sw) if sw > 0 else np.nan
        row['isi_corridor_weight'] = float(np.nansum(w[g.inside_isi_corridor.astype(bool).to_numpy()]) / sw) if sw > 0 else np.nan
        row['active_corridor_weight'] = float(np.nansum(w[g.inside_active_corridor.astype(bool).to_numpy()]) / sw) if sw > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _ensemble_surface(pair_surface):
    rows = []
    metrics = [c for c in pair_surface.columns if c.endswith('_weighted_median')]
    probs = [c for c in pair_surface.columns if c.startswith('p_') and (c.endswith('_isi') or c.endswith('_active'))]
    probs += ['valid_isi_weight','valid_active_weight','isi_corridor_weight','active_corridor_weight']
    for subset, d in [('all_pairs', pair_surface), ('core_secure_pairs', pair_surface[pair_surface.both_core_secure])]:
        for (ii, jj, pi, pdri), g in d.groupby(['intrinsic_index','drive_index','p_intrinsic','p_drive'], sort=True):
            row = {'subset': subset, 'intrinsic_index': int(ii), 'drive_index': int(jj),
                   'p_intrinsic': float(pi), 'p_drive': float(pdri),
                   'n_biological_pairs': int(g.biological_pair_key.nunique())}
            for m in metrics:
                x = pd.to_numeric(g[m], errors='coerce').to_numpy(float)
                x = x[np.isfinite(x)]
                row[m.replace('_weighted_median','') + '_median'] = float(np.median(x)) if len(x) else np.nan
                row[m.replace('_weighted_median','') + '_q25'] = float(np.quantile(x, .25)) if len(x) else np.nan
                row[m.replace('_weighted_median','') + '_q75'] = float(np.quantile(x, .75)) if len(x) else np.nan
                row[m.replace('_weighted_median','') + '_n'] = int(len(x))
            for p in probs:
                x = pd.to_numeric(g[p], errors='coerce').to_numpy(float)
                x = x[np.isfinite(x)]
                row[p + '_mean'] = float(np.mean(x)) if len(x) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def _crossing_rows(pair_surface, geom, cfg, scan='drive'):
    persistence = int(cfg['staging']['persistence_points'])
    rows = []
    stages = [('WT_exit','wt_exit'), ('balance', None), ('SCA3_entry','sca3_entry')]
    for pair, gp in pair_surface.groupby('biological_pair_key'):
        secure = bool(gp.both_core_secure.iloc[0])
        if scan == 'drive':
            outer = 'p_intrinsic'; xcol = 'p_drive'
        else:
            outer = 'p_drive'; xcol = 'p_intrinsic'
        for fixed, gfix in gp.groupby(outer):
            gfix = gfix.sort_values(xcol)
            x = gfix[xcol].to_numpy(float)
            for proj, acol, refkey in [('isi','A_isi_weighted_median','isi'),('active','A_active_weighted_median','active')]:
                y = pd.to_numeric(gfix[acol], errors='coerce').to_numpy(float)
                ref = geom[refkey]
                for stage, attr in stages:
                    thr = 0.5 if attr is None else float(ref[attr])
                    cross = persistent_crossing(x, y, thr, persistence=persistence)
                    rows.append({
                        'biological_pair_key': pair,
                        'wt_cell_id': gp.wt_cell_id.iloc[0], 'sca_cell_id': gp.sca_cell_id.iloc[0],
                        'both_core_secure': secure,
                        'scan': scan, 'fixed_coordinate': outer, 'fixed_value': float(fixed),
                        'projection': proj, 'stage': stage, 'A_threshold': thr,
                        'crossing_coordinate': xcol, 'crossing_value': cross,
                    })
    return pd.DataFrame(rows)


def _summarize_crossings(crossings):
    rows = []
    for subset, d in [('all_pairs', crossings), ('core_secure_pairs', crossings[crossings.both_core_secure])]:
        for keys, g in d.groupby(['scan','fixed_coordinate','fixed_value','projection','stage','A_threshold','crossing_coordinate'], sort=True):
            x = pd.to_numeric(g.crossing_value, errors='coerce').to_numpy(float)
            finite = x[np.isfinite(x)]
            rows.append({
                'subset': subset,
                'scan': keys[0], 'fixed_coordinate': keys[1], 'fixed_value': float(keys[2]),
                'projection': keys[3], 'stage': keys[4], 'A_threshold': float(keys[5]),
                'crossing_coordinate': keys[6],
                'n_biological_pairs_total': int(len(g)), 'n_pairs_with_crossing': int(len(finite)),
                'support_fraction': float(len(finite)/len(g)) if len(g) else np.nan,
                'median': float(np.median(finite)) if len(finite) else np.nan,
                'q25': float(np.quantile(finite,.25)) if len(finite) else np.nan,
                'q75': float(np.quantile(finite,.75)) if len(finite) else np.nan,
            })
    return pd.DataFrame(rows)


def _drive_sensitivity(ensemble, geom, cfg):
    d = ensemble[ensemble.subset.eq('core_secure_pairs')].copy()
    pis = np.sort(d.p_intrinsic.unique())
    pds = np.sort(d.p_drive.unique())
    piv = d.pivot(index='p_intrinsic', columns='p_drive', values='A_isi_median').reindex(index=pis, columns=pds)
    A = piv.to_numpy(float)
    if A.shape[0] < 3 or A.shape[1] < 3:
        return pd.DataFrame(), pd.DataFrame()
    dpi = float(np.median(np.diff(pis)))
    dpd = float(np.median(np.diff(pds)))
    dAi, dAd = np.gradient(A, dpi, dpd)
    denom = np.abs(dAi) + np.abs(dAd)
    frac = np.divide(np.abs(dAd), denom, out=np.full_like(denom, np.nan), where=np.isfinite(denom)&(denom>0))
    grad = np.sqrt(dAi*dAi + dAd*dAd)
    rows = []
    for i, pi in enumerate(pis):
        for j, pdri in enumerate(pds):
            rows.append({'p_intrinsic':float(pi),'p_drive':float(pdri),'A_isi_median':A[i,j],
                         'dA_dintrinsic':dAi[i,j],'dA_ddrive':dAd[i,j],
                         'drive_dominance_fraction':frac[i,j],'gradient_magnitude':grad[i,j]})
    surf = pd.DataFrame(rows)
    band = float(cfg['sensitivity'].get('boundary_A_band', 0.05))
    sumrows = []
    for stage, thr in [('WT_exit',geom['isi']['wt_exit']),('balance',0.5),('SCA3_entry',geom['isi']['sca3_entry'])]:
        x = surf[np.isfinite(surf.A_isi_median) & (np.abs(surf.A_isi_median-float(thr)) <= band)]
        sumrows.append({
            'stage':stage,'A_threshold':float(thr),'A_band':band,'n_grid_points':int(len(x)),
            'median_drive_dominance_fraction':float(np.nanmedian(x.drive_dominance_fraction)) if len(x) else np.nan,
            'median_abs_dA_ddrive':float(np.nanmedian(np.abs(x.dA_ddrive))) if len(x) else np.nan,
            'median_abs_dA_dintrinsic':float(np.nanmedian(np.abs(x.dA_dintrinsic))) if len(x) else np.nan,
            'median_gradient_magnitude':float(np.nanmedian(x.gradient_magnitude)) if len(x) else np.nan,
        })
    return surf, pd.DataFrame(sumrows)


def _drive_curve(p, family):
    p = float(p)
    if family == 'coupled': return p
    if family == 'drive_early': return 1.0-(1.0-p)**2
    if family == 'drive_late': return p*p
    raise ValueError(family)


def _bilinear(pis, pds, M, x, y):
    x = float(np.clip(x, pis[0], pis[-1])); y = float(np.clip(y, pds[0], pds[-1]))
    i = min(np.searchsorted(pis, x, side='right')-1, len(pis)-2); i = max(i,0)
    j = min(np.searchsorted(pds, y, side='right')-1, len(pds)-2); j = max(j,0)
    x0,x1=pis[i],pis[i+1]; y0,y1=pds[j],pds[j+1]
    q00,q01,q10,q11=M[i,j],M[i,j+1],M[i+1,j],M[i+1,j+1]
    if not np.all(np.isfinite([q00,q01,q10,q11])): return np.nan
    tx=0 if x1==x0 else (x-x0)/(x1-x0); ty=0 if y1==y0 else (y-y0)/(y1-y0)
    return float((1-tx)*(1-ty)*q00 + (1-tx)*ty*q01 + tx*(1-ty)*q10 + tx*ty*q11)


def _legacy_surface_recovery(pair_surface, legacy, geom, cfg):
    pgrid = np.linspace(0,1,int(cfg['validation']['legacy_track_points']))
    pairrows=[]
    for pair,gp in pair_surface.groupby('biological_pair_key'):
        if not bool(gp.both_core_secure.iloc[0]): continue
        pis=np.sort(gp.p_intrinsic.unique()); pds=np.sort(gp.p_drive.unique())
        M=gp.pivot(index='p_intrinsic',columns='p_drive',values='A_isi_weighted_median').reindex(index=pis,columns=pds).to_numpy(float)
        for fam in ['drive_early','coupled','drive_late']:
            yy=np.array([_bilinear(pis,pds,M,p,_drive_curve(p,fam)) for p in pgrid])
            row={'biological_pair_key':pair,'path_family':fam}
            for stage,thr in [('wt_exit',geom['isi']['wt_exit']),('balance',0.5),('sca3_entry',geom['isi']['sca3_entry'])]:
                row[stage+'_p_surface']=persistent_crossing(pgrid,yy,float(thr),int(cfg['staging']['persistence_points']))
            pairrows.append(row)
    pairdf=pd.DataFrame(pairrows)
    rows=[]
    legacy2=legacy[(legacy.subset=='core_secure_pairs') & legacy.metric.isin(['wt_exit_p_isi','balance_p_isi','sca3_entry_p_isi'])].copy()
    for fam,g in pairdf.groupby('path_family'):
        for stage,metric in [('wt_exit','wt_exit_p_isi'),('balance','balance_p_isi'),('sca3_entry','sca3_entry_p_isi')]:
            x=pd.to_numeric(g[stage+'_p_surface'],errors='coerce').dropna().to_numpy(float)
            old=legacy2[(legacy2.path_family==fam)&(legacy2.metric==metric)]
            oldmed=float(old.iloc[0]['median']) if len(old) else np.nan
            med=float(np.median(x)) if len(x) else np.nan
            rows.append({'path_family':fam,'stage':stage,'n_pairs':len(x),'surface_recovered_median':med,
                         'frozen_v1_1_median':oldmed,'absolute_difference':abs(med-oldmed) if np.isfinite(med) and np.isfinite(oldmed) else np.nan})
    return pairdf,pd.DataFrame(rows)


def _make_figures(out, ensemble, crossing_summary, sensitivity, legacy, geom):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fdir=out/'figures'; fdir.mkdir(exist_ok=True)
    core=ensemble[ensemble.subset.eq('core_secure_pairs')].copy()
    pis=np.sort(core.p_intrinsic.unique()); pds=np.sort(core.p_drive.unique())

    def matrix(col):
        return core.pivot(index='p_drive',columns='p_intrinsic',values=col).reindex(index=pds,columns=pis).to_numpy(float)

    # Primary ISI stage map.
    M=matrix('A_isi_median')
    fig,ax=plt.subplots(figsize=(8.5,6.5))
    mesh=ax.pcolormesh(pis,pds,M,shading='auto')
    fig.colorbar(mesh,ax=ax,label='median A_ISI')
    cs=ax.contour(pis,pds,M,levels=[geom['isi']['wt_exit'],0.5,geom['isi']['sca3_entry']],linewidths=1.8)
    ax.clabel(cs,fmt={geom['isi']['wt_exit']:'WT-exit',0.5:'balance',geom['isi']['sca3_entry']:'SCA3-entry'})
    p=np.linspace(0,1,200)
    for fam in ['drive_early','coupled','drive_late']:
        ax.plot(p,[_drive_curve(x,fam) for x in p],linestyle='--',linewidth=1,label=fam.replace('_',' '))
    leg=legacy[(legacy.subset=='core_secure_pairs') & legacy.metric.isin(['wt_exit_p_isi','balance_p_isi','sca3_entry_p_isi'])]
    for _,r in leg.iterrows():
        pp=float(r['median']); ax.scatter([pp],[_drive_curve(pp,r['path_family'])],s=28)
    ax.set(xlabel='intrinsic progress p_intrinsic',ylabel='drive progress p_drive',
           title='WT→SCA3 intrinsic × drive map: primary ISI projection')
    ax.legend(loc='lower right',fontsize=8)
    fig.tight_layout(); fig.savefig(fdir/'01_primary_ISI_stage_map_core_secure.png',dpi=220); plt.close(fig)

    # Required drive for each boundary at fixed intrinsic progress.
    d=crossing_summary[(crossing_summary.subset=='core_secure_pairs')&(crossing_summary.scan=='drive')&(crossing_summary.projection=='isi')]
    fig,ax=plt.subplots(figsize=(8.5,5.8))
    for stage in ['WT_exit','balance','SCA3_entry']:
        x=d[d.stage==stage].sort_values('fixed_value')
        ax.plot(x.fixed_value,x['median'],label=stage.replace('_','-'))
        ax.fill_between(x.fixed_value.to_numpy(float),x.q25.to_numpy(float),x.q75.to_numpy(float),alpha=.15)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='fixed intrinsic progress',ylabel='drive progress required to cross boundary',
           title='Drive required for stage crossing (core-secure biological pairs)')
    ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(fdir/'02_required_drive_boundaries_core_secure.png',dpi=220); plt.close(fig)

    # Drive dominance.
    if len(sensitivity):
        S=sensitivity.pivot(index='p_drive',columns='p_intrinsic',values='drive_dominance_fraction').reindex(index=pds,columns=pis).to_numpy(float)
        fig,ax=plt.subplots(figsize=(8.5,6.2))
        mesh=ax.pcolormesh(pis,pds,S,shading='auto',vmin=0,vmax=1)
        fig.colorbar(mesh,ax=ax,label='|dA/d drive| / (|dA/d drive| + |dA/d intrinsic|)')
        ax.contour(pis,pds,M,levels=[geom['isi']['wt_exit'],0.5,geom['isi']['sca3_entry']],linewidths=1.2)
        ax.set(xlabel='intrinsic progress p_intrinsic',ylabel='drive progress p_drive',title='Local drive dominance of the ISI transition coordinate')
        fig.tight_layout(); fig.savefig(fdir/'03_drive_dominance_map_core_secure.png',dpi=220); plt.close(fig)

    # Secondary active-rate map.
    A=matrix('A_active_median')
    fig,ax=plt.subplots(figsize=(8.5,6.5))
    mesh=ax.pcolormesh(pis,pds,A,shading='auto')
    fig.colorbar(mesh,ax=ax,label='median A_active')
    cs=ax.contour(pis,pds,A,levels=[geom['active']['wt_exit'],0.5,geom['active']['sca3_entry']],linewidths=1.8)
    ax.clabel(cs,fmt={geom['active']['wt_exit']:'WT-exit',0.5:'balance',geom['active']['sca3_entry']:'SCA3-entry'})
    ax.set(xlabel='intrinsic progress p_intrinsic',ylabel='drive progress p_drive',title='Secondary active-rate transition map')
    fig.tight_layout(); fig.savefig(fdir/'04_secondary_active_rate_stage_map_core_secure.png',dpi=220); plt.close(fig)


def validate(cfg):
    support,cells,anchors,refs,transforms,boundaries,legacy=load_inputs(cfg)
    support=build_support_with_anchors(support,anchors)
    scenarios=enumerate_scenarios(support)
    geom=load_geometry(refs,transforms)
    out={
        'version':'1.2.0','analysis':'intrinsic_x_drive_transition_surface',
        'primary_cells':int(cells.cell_id.nunique()),'WT_cells':int((cells.group=='WT').sum()),'SCA3_cells':int((cells.group=='SCA3').sum()),
        'support_states':int(len(support)),'biological_pairs':int(scenarios.biological_pair_key.nunique()),
        'support_pair_scenarios':int(len(scenarios)),'scenario_weight_sum':float(scenarios.scenario_weight.sum()),
        'n_intrinsic':int(cfg['surface']['n_intrinsic']),'n_drive':int(cfg['surface']['n_drive']),
        'states_full_profile':int(len(scenarios)*int(cfg['surface']['n_intrinsic'])*int(cfg['surface']['n_drive'])),
        'primary_ISI_wt_exit_A':float(geom['isi']['wt_exit']),'primary_ISI_balance_A':0.5,
        'primary_ISI_sca3_entry_A':float(geom['isi']['sca3_entry']),
        'active_wt_exit_A':float(geom['active']['wt_exit']),'active_sca3_entry_A':float(geom['active']['sca3_entry']),
        'all_q75_supported':bool(anchors.q75_supported.all()),'frozen_sha256':frozen_manifest(cfg['data']['root'])}
    if (out['primary_cells'],out['WT_cells'],out['SCA3_cells']) != (18,12,6): raise ValueError('unexpected primary cohort')
    if out['support_states']!=64 or out['support_pair_scenarios']!=988: raise ValueError('unexpected support-state counts')
    if abs(out['scenario_weight_sum']-1.0)>1e-9: raise ValueError('scenario weights do not sum to one')
    if not out['all_q75_supported']: raise ValueError('q=.75 not supported for every primary cell')
    if geom['isi']['cloud_overlap'] or geom['active']['cloud_overlap']: raise ValueError('reference endpoint clouds overlap')
    if abs(out['primary_ISI_wt_exit_A']-0.135829)>1e-4 or abs(out['primary_ISI_sca3_entry_A']-0.797856)>1e-4:
        raise ValueError('unexpected frozen v1.1 ISI boundaries')
    return out


def run_all(cfg, resume=True):
    out=Path(cfg['output']['dir']); out.mkdir(parents=True,exist_ok=True)
    support,cells,anchors,refs,transforms,boundaries,legacy=load_inputs(cfg)
    support=build_support_with_anchors(support,anchors)
    scenarios=enumerate_scenarios(support)
    mode=str(cfg['surface'].get('scenario_mode','all_support'))
    if mode=='best_only': scenarios=scenarios[(scenarios.wt_source=='best')&(scenarios.sca_source=='best')].copy()
    maxs=cfg['surface'].get('max_scenarios')
    if maxs is not None: scenarios=scenarios.head(int(maxs)).copy()
    geom=load_geometry(refs,transforms)
    scenarios.to_csv(out/'transition_pair_scenarios_v1_2.csv',index=False)
    boundaries.to_csv(out/'staging_boundary_definitions_frozen_v1_1.csv',index=False)

    jobs=[]
    for _,r in scenarios.iterrows():
        sc=r.to_dict(); cp=_checkpoint_path(out,sc['scenario_id'])
        if not (resume and cp.exists()): jobs.append((sc,cp))
    workers=int(cfg['parallel']['workers'])
    if jobs:
        if workers<=1:
            for sc,cp in jobs:
                _write_checkpoint(cp,_scenario_surface_job(sc,cfg,geom))
        else:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                fut={ex.submit(_scenario_surface_job,sc,cfg,geom):(cp,sc['scenario_id']) for sc,cp in jobs}
                for f in as_completed(fut):
                    cp,sid=fut[f]; _write_checkpoint(cp,f.result())

    # Stream full state table without keeping the full ~1M-row server result in RAM.
    if bool(cfg['output'].get('save_full_states',True)):
        full=out/'transition_surface_states.csv.gz'
        first=True
        for sid in scenarios.scenario_id:
            d=_read_checkpoint(_checkpoint_path(out,sid))
            d.to_csv(full,index=False,compression='gzip',mode='w' if first else 'a',header=first)
            first=False

    secure_cells=set(cells.loc[cells.core_q75_secure.fillna(False).astype(bool),'cell_id'].astype(str))
    pair_parts=[]
    for pair,g in scenarios.groupby('biological_pair_key',sort=True):
        pair_parts.append(_pair_surface_for_pair(g,out,secure_cells))
    pair_surface=pd.concat(pair_parts,ignore_index=True)
    pair_surface.to_csv(out/'biological_pair_surface_summary.csv.gz',index=False,compression='gzip')
    ensemble=_ensemble_surface(pair_surface)
    ensemble.to_csv(out/'ensemble_surface_summary.csv',index=False)

    cross_drive=_crossing_rows(pair_surface,geom,cfg,scan='drive')
    cross_intr=_crossing_rows(pair_surface,geom,cfg,scan='intrinsic')
    crossings=pd.concat([cross_drive,cross_intr],ignore_index=True)
    crossings.to_csv(out/'biological_pair_boundary_crossings.csv.gz',index=False,compression='gzip')
    crossing_summary=_summarize_crossings(crossings)
    crossing_summary.to_csv(out/'boundary_curve_summary.csv',index=False)

    sens,sens_summary=_drive_sensitivity(ensemble,geom,cfg)
    sens.to_csv(out/'drive_sensitivity_surface.csv',index=False)
    sens_summary.to_csv(out/'drive_sensitivity_at_stage_boundaries.csv',index=False)

    legacy_pair,legacy_summary=_legacy_surface_recovery(pair_surface,legacy,geom,cfg)
    legacy_pair.to_csv(out/'legacy_path_surface_recovery_by_pair.csv',index=False)
    legacy_summary.to_csv(out/'legacy_path_surface_recovery_summary.csv',index=False)

    _make_figures(out,ensemble,crossing_summary,sens,legacy,geom)

    core=ensemble[ensemble.subset.eq('core_secure_pairs')]
    summary={
        'version':'1.2.0','analysis':'WT_to_SCA3_intrinsic_x_drive_surface',
        'scenario_mode':mode,'n_scenarios':int(scenarios.scenario_id.nunique()),
        'n_biological_pairs':int(scenarios.biological_pair_key.nunique()),
        'n_intrinsic':int(cfg['surface']['n_intrinsic']),'n_drive':int(cfg['surface']['n_drive']),
        'n_state_rows':int(len(scenarios)*int(cfg['surface']['n_intrinsic'])*int(cfg['surface']['n_drive'])),
        'core_secure_pairs':int(pair_surface[pair_surface.both_core_secure].biological_pair_key.nunique()),
        'primary_ISI_boundaries':{'WT_exit':geom['isi']['wt_exit'],'balance':0.5,'SCA3_entry':geom['isi']['sca3_entry']},
        'secondary_active_boundaries':{'WT_exit':geom['active']['wt_exit'],'balance':0.5,'SCA3_entry':geom['active']['sca3_entry']},
        'median_ISI_valid_weight_core_surface':float(np.nanmedian(core.valid_isi_weight_mean)) if len(core) else np.nan,
        'median_ISI_corridor_weight_core_surface':float(np.nanmedian(core.isi_corridor_weight_mean)) if len(core) else np.nan,
        'legacy_surface_recovery_max_abs_difference':float(np.nanmax(legacy_summary.absolute_difference)) if len(legacy_summary) else np.nan,
        'frozen_sha256':frozen_manifest(cfg['data']['root'])}
    (out/'RUN_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True),encoding='utf-8')
    return summary
