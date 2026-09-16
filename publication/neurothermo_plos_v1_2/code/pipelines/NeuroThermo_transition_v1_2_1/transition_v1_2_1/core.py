from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .io import V12Results
from .mathutils import persistent_crossing, weighted_finite_quantile, crossing_support_weight, drive_curve, bilinear_track

STAGES = [('WT_exit','wt_exit'),('balance','balance'),('SCA3_entry','sca3_entry')]
PROJS = [('isi','A_isi'),('active','A_active')]
SENS_METRICS = ['A_isi','dA_dintrinsic','dA_ddrive','drive_dominance_fraction','gradient_magnitude']


def _resolve(cfg):
    return V12Results(cfg['input'].get('directory'), cfg['input'].get('archive'))


def _thresholds(run):
    return {
      'isi': {'wt_exit':float(run['primary_ISI_boundaries']['WT_exit']), 'balance':0.5, 'sca3_entry':float(run['primary_ISI_boundaries']['SCA3_entry'])},
      'active': {'wt_exit':float(run['secondary_active_boundaries']['WT_exit']), 'balance':0.5, 'sca3_entry':float(run['secondary_active_boundaries']['SCA3_entry'])},
    }


def validate(cfg):
    src=_resolve(cfg)
    run=src.read_json('RUN_SUMMARY.json')
    scenarios=src.read_csv('transition_pair_scenarios_v1_2.csv')
    pair_old=src.read_csv('biological_pair_surface_summary.csv.gz',compression='gzip',usecols=['biological_pair_key','both_core_secure']).drop_duplicates()
    cps=src.checkpoint_ids()
    out={
      'version':'1.2.1','analysis':'scenario_first_uncertainty_postprocessing',
      'source':src.source_description,'source_version':run.get('version'),
      'n_scenarios':int(scenarios.scenario_id.nunique()),
      'n_biological_pairs':int(scenarios.biological_pair_key.nunique()),
      'core_secure_pairs':int(pair_old.loc[pair_old.both_core_secure.astype(bool),'biological_pair_key'].nunique()),
      'n_checkpoints':len(cps), 'scenario_weight_sum':float(scenarios.scenario_weight.sum()),
      'n_intrinsic':int(run['n_intrinsic']),'n_drive':int(run['n_drive']),
      'primary_ISI_boundaries':run['primary_ISI_boundaries'],
      'secondary_active_boundaries':run['secondary_active_boundaries'],
      'new_HR_simulations_required':False,
    }
    src.close()
    if out['source_version']!='1.2.0': raise ValueError('v1.2.0 results required')
    if out['n_scenarios']!=988 or out['n_biological_pairs']!=72 or out['core_secure_pairs']!=32: raise ValueError('unexpected v1.2 cohort/scenario counts')
    if out['n_checkpoints']!=out['n_scenarios']: raise ValueError('checkpoint count does not match scenarios')
    if abs(out['scenario_weight_sum']-1)>1e-9: raise ValueError('scenario weights do not sum to one')
    if abs(float(run['primary_ISI_boundaries']['WT_exit'])-0.1358293233470019)>1e-9: raise ValueError('unexpected ISI WT boundary')
    if abs(float(run['primary_ISI_boundaries']['SCA3_entry'])-0.7978563373093712)>1e-9: raise ValueError('unexpected ISI SCA3 boundary')
    return out


def _scenario_crossings(d, scrow, secure, thresholds, persistence):
    rows=[]
    common={'scenario_id':int(scrow.scenario_id),'biological_pair_key':scrow.biological_pair_key,
            'wt_cell_id':scrow.wt_cell_id,'sca_cell_id':scrow.sca_cell_id,
            'both_core_secure':bool(secure),'within_pair_support_weight':float(scrow.within_pair_support_weight)}
    for scan,outer,xcol in [('drive','p_intrinsic','p_drive'),('intrinsic','p_drive','p_intrinsic')]:
        for fixed,g in d.groupby(outer,sort=True):
            g=g.sort_values(xcol); x=g[xcol].to_numpy(float)
            for proj,acol in PROJS:
                y=pd.to_numeric(g[acol],errors='coerce').to_numpy(float)
                for stage,key in STAGES:
                    thr=thresholds[proj][key]
                    rows.append({**common,'scan':scan,'fixed_coordinate':outer,'fixed_value':float(fixed),
                                 'projection':proj,'stage':stage,'A_threshold':float(thr),
                                 'crossing_coordinate':xcol,'crossing_value':persistent_crossing(x,y,thr,persistence)})
    return rows


def _finite_median(x):
    a=np.asarray(x,float); a=a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else np.nan


def _scenario_sensitivity(d, scrow, secure, thresholds, band):
    pis=np.sort(d.p_intrinsic.unique()); pds=np.sort(d.p_drive.unique())
    M=d.pivot(index='p_intrinsic',columns='p_drive',values='A_isi').reindex(index=pis,columns=pds).to_numpy(float)
    dpi=float(np.median(np.diff(pis))); dpd=float(np.median(np.diff(pds)))
    dAi,dAd=np.gradient(M,dpi,dpd)
    denom=np.abs(dAi)+np.abs(dAd)
    frac=np.divide(np.abs(dAd),denom,out=np.full_like(denom,np.nan),where=np.isfinite(denom)&(denom>0))
    grad=np.sqrt(dAi*dAi+dAd*dAd)
    arrays={'A_isi':M,'dA_dintrinsic':dAi,'dA_ddrive':dAd,'drive_dominance_fraction':frac,'gradient_magnitude':grad}
    common={'scenario_id':int(scrow.scenario_id),'biological_pair_key':scrow.biological_pair_key,
            'both_core_secure':bool(secure),'within_pair_support_weight':float(scrow.within_pair_support_weight)}
    stage_rows=[]
    for stage,key in STAGES:
        thr=thresholds['isi'][key]; mask=np.isfinite(M)&(np.abs(M-thr)<=band)
        stage_rows.append({**common,'stage':stage,'A_threshold':float(thr),'A_band':float(band),
          'median_drive_dominance_fraction':_finite_median(frac[mask]),
          'median_abs_dA_ddrive':_finite_median(np.abs(dAd[mask])),
          'median_abs_dA_dintrinsic':_finite_median(np.abs(dAi[mask])),
          'median_gradient_magnitude':_finite_median(grad[mask]),
          'n_grid_points':int(np.sum(mask))})
    return pis,pds,arrays,stage_rows


def _legacy_scenario(d, scrow, secure, thresholds, persistence, ntrack):
    pis=np.sort(d.p_intrinsic.unique()); pds=np.sort(d.p_drive.unique())
    M=d.pivot(index='p_intrinsic',columns='p_drive',values='A_isi').reindex(index=pis,columns=pds).to_numpy(float)
    p=np.linspace(0,1,int(ntrack)); rows=[]
    for fam in ['drive_early','coupled','drive_late']:
        yy=bilinear_track(pis,pds,M,p,drive_curve(p,fam))
        row={'scenario_id':int(scrow.scenario_id),'biological_pair_key':scrow.biological_pair_key,
             'both_core_secure':bool(secure),'within_pair_support_weight':float(scrow.within_pair_support_weight),'path_family':fam}
        for stage,key in STAGES:
            row[stage+'_p']=persistent_crossing(p,yy,thresholds['isi'][key],persistence)
        rows.append(row)
    return rows


def _weighted_finite_matrix(V,w):
    """Weighted finite median per matrix element; V shape (n_solutions,n_i,n_d)."""
    V=np.asarray(V,float); w=np.asarray(w,float)
    n,ni,nd=V.shape; X=V.reshape(n,ni*nd)
    valid=np.isfinite(X)&np.isfinite(w[:,None])&(w[:,None]>0)
    vals=np.where(valid,X,np.inf)
    order=np.argsort(vals,axis=0)
    sv=np.take_along_axis(vals,order,axis=0)
    wb=np.broadcast_to(w[:,None],X.shape)
    sw=np.take_along_axis(np.where(valid,wb,0.0),order,axis=0)
    cs=np.cumsum(sw,axis=0); total=np.sum(sw,axis=0); target=.5*total
    idx=np.argmax(cs>=target[None,:],axis=0)
    med=sv[idx,np.arange(X.shape[1])]
    med[(total<=0)|~np.isfinite(med)]=np.nan
    return med.reshape(ni,nd)


def _pair_sensitivity_rows(pair,secure,pis,pds,scenario_arrays,weights):
    rows=[]
    stack={m:np.stack([a[m] for a in scenario_arrays],axis=0) for m in SENS_METRICS}
    med={m:_weighted_finite_matrix(stack[m],weights) for m in SENS_METRICS}
    for i,pi in enumerate(pis):
        for j,pdri in enumerate(pds):
            row={'biological_pair_key':pair,'both_core_secure':bool(secure),'p_intrinsic':float(pi),'p_drive':float(pdri)}
            for m in SENS_METRICS: row[m+'_weighted_median']=med[m][i,j]
            rows.append(row)
    return rows


def _aggregate_pair_crossing_frames(frames, weights):
    """Aggregate deterministic scenario crossing frames within one biological pair."""
    if not frames:
        return pd.DataFrame()
    base=frames[0].drop(columns=['scenario_id','crossing_value','within_pair_support_weight']).copy()
    V=np.stack([f.crossing_value.to_numpy(float) for f in frames],axis=0)
    w=np.asarray(weights,float)
    base['n_support_scenarios']=len(frames)
    base['crossing_support_weight']=[crossing_support_weight(V[:,j],w) for j in range(V.shape[1])]
    base['q25_weighted']=[weighted_finite_quantile(V[:,j],w,.25) for j in range(V.shape[1])]
    base['median_weighted']=[weighted_finite_quantile(V[:,j],w,.50) for j in range(V.shape[1])]
    base['q75_weighted']=[weighted_finite_quantile(V[:,j],w,.75) for j in range(V.shape[1])]
    return base


def _ensemble_crossings(pair):
    erows=[]
    gkeys=['scan','fixed_coordinate','fixed_value','projection','stage','A_threshold','crossing_coordinate']
    for subset,d in [('all_pairs',pair),('core_secure_pairs',pair[pair.both_core_secure.astype(bool)])]:
        for k,g in d.groupby(gkeys,sort=True,dropna=False):
            x=pd.to_numeric(g.median_weighted,errors='coerce').to_numpy(float); finite=x[np.isfinite(x)]
            sw=pd.to_numeric(g.crossing_support_weight,errors='coerce').to_numpy(float); sw=sw[np.isfinite(sw)]
            maj=g[pd.to_numeric(g.crossing_support_weight,errors='coerce')>=0.5]
            xm=pd.to_numeric(maj.median_weighted,errors='coerce').dropna().to_numpy(float)
            erows.append({'subset':subset,**dict(zip(gkeys,k)),'n_biological_pairs_total':int(len(g)),
              'n_pairs_with_marker':int(len(finite)),'pair_marker_fraction':float(len(finite)/len(g)) if len(g) else np.nan,
              'n_pairs_support_ge_0_5':int(len(xm)),'pair_majority_support_fraction':float(len(xm)/len(g)) if len(g) else np.nan,
              'median_pair_crossing_support_weight':float(np.median(sw)) if len(sw) else np.nan,
              'median':float(np.median(finite)) if len(finite) else np.nan,'q25':float(np.quantile(finite,.25)) if len(finite) else np.nan,'q75':float(np.quantile(finite,.75)) if len(finite) else np.nan,
              'median_majority_support':float(np.median(xm)) if len(xm) else np.nan,'q25_majority_support':float(np.quantile(xm,.25)) if len(xm) else np.nan,'q75_majority_support':float(np.quantile(xm,.75)) if len(xm) else np.nan})
    return pd.DataFrame(erows)


def _aggregate_sensitivity_grid(pair):
    parts=[]
    cols=[m+'_weighted_median' for m in SENS_METRICS]
    for subset,d in [('all_pairs',pair),('core_secure_pairs',pair[pair.both_core_secure.astype(bool)])]:
        g=d.groupby(['p_intrinsic','p_drive'],sort=True)
        med=g[cols].median().rename(columns={c:c.replace('_weighted_median','_median') for c in cols})
        q25=g[cols].quantile(.25).rename(columns={c:c.replace('_weighted_median','_q25') for c in cols})
        q75=g[cols].quantile(.75).rename(columns={c:c.replace('_weighted_median','_q75') for c in cols})
        n=g.size().rename('n_biological_pairs')
        z=pd.concat([n,med,q25,q75],axis=1).reset_index();z.insert(0,'subset',subset);parts.append(z)
    return pd.concat(parts,ignore_index=True)


def _weighted_finite(values,weights,q=.5):
    v=np.asarray(values,float); w=np.asarray(weights,float); m=np.isfinite(v)&np.isfinite(w)&(w>0)
    if not np.any(m): return np.nan
    v=v[m];w=w[m];o=np.argsort(v);v=v[o];w=w[o];c=np.cumsum(w)/np.sum(w)
    return float(np.interp(q,np.concatenate([[0],c]),np.concatenate([[v[0]],v])))


def _figures(out, ensemble_boundary, ens_sens_surface, legacy_summary):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fdir=out/'figures';fdir.mkdir(exist_ok=True)
    d=ensemble_boundary[(ensemble_boundary.subset=='core_secure_pairs')&(ensemble_boundary.scan=='drive')&(ensemble_boundary.projection=='isi')]
    fig,ax=plt.subplots(figsize=(8.8,6.0))
    for stage in ['WT_exit','balance','SCA3_entry']:
        x=d[d.stage==stage].sort_values('fixed_value')
        ax.plot(x.fixed_value,x['median'],linestyle=':',alpha=.45,label=stage.replace('_','-')+' (any marker)')
        ax.plot(x.fixed_value,x['median_majority_support'],label=stage.replace('_','-')+' (support>=0.5)')
        ax.fill_between(x.fixed_value.to_numpy(float),x.q25_majority_support.to_numpy(float),x.q75_majority_support.to_numpy(float),alpha=.15)
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='fixed intrinsic progress',ylabel='drive progress required for crossing',title='Scenario-first WT→SCA3 boundaries (core-secure pairs)');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(fdir/'01_corrected_required_drive_boundaries_core_secure.png',dpi=220);plt.close(fig)
    s=ens_sens_surface[ens_sens_surface.subset=='core_secure_pairs'];pis=np.sort(s.p_intrinsic.unique());pds=np.sort(s.p_drive.unique())
    M=s.pivot(index='p_drive',columns='p_intrinsic',values='drive_dominance_fraction_median').reindex(index=pds,columns=pis).to_numpy(float)
    fig,ax=plt.subplots(figsize=(8.8,6.2));mesh=ax.pcolormesh(pis,pds,M,shading='auto',vmin=0,vmax=1);fig.colorbar(mesh,ax=ax,label='median scenario-first drive dominance');ax.set(xlabel='intrinsic progress',ylabel='drive progress',title='Scenario-aware drive dominance (core-secure pairs)');fig.tight_layout();fig.savefig(fdir/'02_scenario_aware_drive_dominance_core_secure.png',dpi=220);plt.close(fig)
    order=['drive_early','coupled','drive_late'];stages=['WT_exit','balance','SCA3_entry'];fig,ax=plt.subplots(figsize=(9,5.8));x=np.arange(len(order));offs=[-.2,0,.2]
    for off,st in zip(offs,stages):
        vals=[]; olds=[]
        for fam in order:
            q=legacy_summary[(legacy_summary.path_family==fam)&(legacy_summary.stage==st)].iloc[0];vals.append(q.scenario_first_surface_median);olds.append(q.frozen_v1_1_median)
        ax.plot(x+off,vals,'o-',label=f'{st} surface');ax.scatter(x+off,olds,marker='x',s=50,label=f'{st} frozen v1.1')
    ax.set_xticks(x,[q.replace('_',' ') for q in order]);ax.set_ylim(0,1);ax.set_ylabel('path progress p');ax.set_title('Recovery of frozen v1.1 staging from v1.2 surfaces');ax.grid(alpha=.25);ax.legend(ncol=2,fontsize=8);fig.tight_layout();fig.savefig(fdir/'03_legacy_path_recovery_scenario_first.png',dpi=220);plt.close(fig)


def run_all(cfg):
    src=_resolve(cfg); run=src.read_json('RUN_SUMMARY.json'); thresholds=_thresholds(run); scenarios=src.read_csv('transition_pair_scenarios_v1_2.csv')
    pair_old=src.read_csv('biological_pair_surface_summary.csv.gz',compression='gzip',usecols=['biological_pair_key','both_core_secure']).drop_duplicates(); secure_map=dict(zip(pair_old.biological_pair_key,pair_old.both_core_secure.astype(bool)))
    out=Path(cfg['output']['dir']).expanduser().resolve();out.mkdir(parents=True,exist_ok=True)
    persistence=int(cfg['staging']['persistence_points']);ntrack=int(cfg['staging']['legacy_track_points']);band=float(cfg['sensitivity']['boundary_A_band'])
    pair_cross_frames=[]; pair_sens_rows=[]; pair_stage_rows=[]; pair_legacy_rows=[]
    for pair,gsc in scenarios.groupby('biological_pair_key',sort=True):
        secure=secure_map.get(pair,False); pair_arrays=[]; pair_weights=[]; pis0=pds0=None; cross_frames=[]; stage_frames=[]; legacy_frames=[]
        for _,sc in gsc.iterrows():
            d=src.read_checkpoint(int(sc.scenario_id))
            cross_frames.append(pd.DataFrame(_scenario_crossings(d,sc,secure,thresholds,persistence)))
            pis,pds,arrays,sr=_scenario_sensitivity(d,sc,secure,thresholds,band);stage_frames.append(pd.DataFrame(sr))
            legacy_frames.append(pd.DataFrame(_legacy_scenario(d,sc,secure,thresholds,persistence,ntrack)))
            if pis0 is None: pis0,pds0=pis,pds
            pair_arrays.append(arrays);pair_weights.append(float(sc.within_pair_support_weight))
        w=np.asarray(pair_weights,float)
        pair_cross_frames.append(_aggregate_pair_crossing_frames(cross_frames,w))
        pair_sens_rows.extend(_pair_sensitivity_rows(pair,secure,pis0,pds0,pair_arrays,w))
        S=pd.concat(stage_frames,ignore_index=True)
        for stage,g in S.groupby('stage',sort=False):
            row={'biological_pair_key':pair,'both_core_secure':bool(secure),'stage':stage,'A_threshold':float(g.A_threshold.iloc[0])}
            for m in ['median_drive_dominance_fraction','median_abs_dA_ddrive','median_abs_dA_dintrinsic','median_gradient_magnitude']:
                row[m+'_weighted_median']=_weighted_finite(g[m],w,.5)
            den=w.sum(); row['scenario_support_with_boundary_band']=float(np.sum(w[pd.to_numeric(g.n_grid_points,errors='coerce').to_numpy(float)>0])/den) if den>0 else np.nan
            pair_stage_rows.append(row)
        L=pd.concat(legacy_frames,ignore_index=True)
        for fam,g in L.groupby('path_family',sort=False):
            row={'biological_pair_key':pair,'both_core_secure':bool(secure),'path_family':fam}
            for stage,_ in STAGES:
                vals=pd.to_numeric(g[stage+'_p'],errors='coerce').to_numpy(float);row[stage+'_support_weight']=crossing_support_weight(vals,w);row[stage+'_p_weighted_median']=weighted_finite_quantile(vals,w,.5)
            pair_legacy_rows.append(row)
    src.close()
    pair_cross=pd.concat(pair_cross_frames,ignore_index=True);ens_cross=_ensemble_crossings(pair_cross)
    pair_cross.to_csv(out/'biological_pair_boundary_crossings_v1_2_1.csv',index=False);ens_cross.to_csv(out/'boundary_curve_summary_v1_2_1.csv',index=False)
    pair_sens_grid=pd.DataFrame(pair_sens_rows);ens_sens_grid=_aggregate_sensitivity_grid(pair_sens_grid)
    pair_sens_grid.to_csv(out/'biological_pair_drive_sensitivity_surface_v1_2_1.csv',index=False);ens_sens_grid.to_csv(out/'drive_sensitivity_surface_v1_2_1.csv',index=False)
    pair_sens_stage=pd.DataFrame(pair_stage_rows);ens_sens_stage=_ensemble_stage_sensitivity_from_pair(pair_sens_stage)
    pair_sens_stage.to_csv(out/'biological_pair_drive_sensitivity_at_boundaries_v1_2_1.csv',index=False);ens_sens_stage.to_csv(out/'drive_sensitivity_at_stage_boundaries_v1_2_1.csv',index=False)
    legacy_pair=pd.DataFrame(pair_legacy_rows);frozen_path=Path(__file__).resolve().parent.parent/'frozen'/'PRIMARY_ISI_STAGING.csv';legacy_summary=_legacy_summary_from_pair(legacy_pair,frozen_path)
    legacy_pair.to_csv(out/'biological_pair_legacy_path_recovery_v1_2_1.csv',index=False);legacy_summary.to_csv(out/'legacy_path_recovery_summary_v1_2_1.csv',index=False)
    _figures(out,ens_cross,ens_sens_grid,legacy_summary)
    core_sens=ens_sens_stage[ens_sens_stage.subset=='core_secure_pairs']; summary={'version':'1.2.1','analysis':'scenario_first_uncertainty_postprocessing','new_HR_simulations':0,'source_v1_2_state_rows':int(run['n_state_rows']),'n_scenarios':int(scenarios.scenario_id.nunique()),'n_biological_pairs':int(scenarios.biological_pair_key.nunique()),'core_secure_pairs':32,'primary_ISI_boundaries':thresholds['isi'],'secondary_active_boundaries':thresholds['active'],'legacy_recovery_max_abs_difference':float(np.nanmax(legacy_summary.absolute_difference)),'scenario_first_sensitivity_core_secure':core_sens.to_dict(orient='records')}
    (out/'RUN_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True),encoding='utf-8');_write_report(out,legacy_summary,core_sens);return summary


def _ensemble_stage_sensitivity_from_pair(pair):
    metrics=['median_drive_dominance_fraction','median_abs_dA_ddrive','median_abs_dA_dintrinsic','median_gradient_magnitude','scenario_support_with_boundary_band']
    rows=[]
    for subset,d in [('all_pairs',pair),('core_secure_pairs',pair[pair.both_core_secure.astype(bool)])]:
        for (st,thr),g in d.groupby(['stage','A_threshold'],sort=True):
            row={'subset':subset,'stage':st,'A_threshold':float(thr),'n_biological_pairs':int(len(g))}
            for m in metrics:
                col=m+'_weighted_median' if m!='scenario_support_with_boundary_band' else m
                x=pd.to_numeric(g[col],errors='coerce').dropna().to_numpy(float);row[m+'_median']=float(np.median(x)) if len(x) else np.nan;row[m+'_q25']=float(np.quantile(x,.25)) if len(x) else np.nan;row[m+'_q75']=float(np.quantile(x,.75)) if len(x) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def _legacy_summary_from_pair(pair,frozen_path):
    frozen=pd.read_csv(frozen_path);rows=[]
    for fam,g in pair[pair.both_core_secure.astype(bool)].groupby('path_family'):
        for stage,metric in [('WT_exit','wt_exit_p_isi'),('balance','balance_p_isi'),('SCA3_entry','sca3_entry_p_isi')]:
            x=pd.to_numeric(g[stage+'_p_weighted_median'],errors='coerce').dropna().to_numpy(float);med=float(np.median(x)) if len(x) else np.nan
            old=frozen[(frozen.path_family==fam)&(frozen.metric==metric)];oldmed=float(old.iloc[0]['median']) if len(old) else np.nan
            sup=pd.to_numeric(g[stage+'_support_weight'],errors='coerce').dropna().to_numpy(float)
            rows.append({'path_family':fam,'stage':stage,'n_core_pairs_total':int(len(g)),'n_pairs_with_marker':int(len(x)),
              'median_pair_crossing_support_weight':float(np.median(sup)) if len(sup) else np.nan,
              'scenario_first_surface_median':med,'frozen_v1_1_median':oldmed,'absolute_difference':abs(med-oldmed) if np.isfinite(med) and np.isfinite(oldmed) else np.nan})
    return pd.DataFrame(rows)


def _write_report(out,legacy,core_sens):
    lines=['# Transition ensemble v1.2.1','','No new Hindmarsh–Rose simulations are performed. The frozen v1.2 checkpoint surfaces are reprocessed with scenario-first uncertainty propagation.','','Correct order: scenario crossing -> v1.1-compatible within-pair weighted marker distribution with crossing-support weight retained separately -> equal-weight biological-pair ensemble.','','## Recovery of frozen v1.1 one-dimensional paths','']
    for _,r in legacy.iterrows(): lines.append(f"- {r.path_family} / {r.stage}: surface={r.scenario_first_surface_median:.6f}, frozen v1.1={r.frozen_v1_1_median:.6f}, |delta|={r.absolute_difference:.6f}.")
    lines += ['','## Scenario-aware drive sensitivity at ISI stage boundaries','']
    for _,r in core_sens.iterrows(): lines.append(f"- {r.stage}: drive dominance median={r.median_drive_dominance_fraction_median:.4f}; |dA/ddrive|={r.median_abs_dA_ddrive_median:.4f}; |dA/dintrinsic|={r.median_abs_dA_dintrinsic_median:.4f}.")
    lines += ['','Parameter-support states without a marker are represented by the separately reported crossing-support weight; marker quantiles use the same finite-marker weighting convention as frozen v1.1 for direct comparability.',''];(out/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
