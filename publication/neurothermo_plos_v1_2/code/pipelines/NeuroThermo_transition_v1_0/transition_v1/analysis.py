from __future__ import annotations
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip, json, math, os
import numpy as np
import pandas as pd

from .data import load_inputs, build_support_with_anchors, enumerate_scenarios, theta_from_row, frozen_manifest
from .model import pre_relax, refine_rheobase, supported_metrics
from .geometry import projection_reference, project_points, persistent_crossing, count_reversals, inverse_morph


def _drive_progress(p,family):
    p=float(p)
    if family=='coupled': return p
    if family=='drive_early': return 1.0-(1.0-p)**2
    if family=='drive_late': return p*p
    raise ValueError(f'unknown path family {family}')


def _interp_theta(sc,p,family):
    pdri=_drive_progress(p,family)
    # b and s are interpolated linearly; positive wide-range r and kappa_I in log space.
    b=(1-p)*float(sc['wt_b'])+p*float(sc['sca_b'])
    s=(1-p)*float(sc['wt_s'])+p*float(sc['sca_s'])
    r=math.exp((1-p)*math.log(float(sc['wt_r']))+p*math.log(float(sc['sca_r'])))
    k=math.exp((1-pdri)*math.log(float(sc['wt_kappa_I']))+pdri*math.log(float(sc['sca_kappa_I'])))
    J=(1-pdri)*float(sc['wt_J_q75'])+pdri*float(sc['sca_J_q75'])
    win=(1-p)*float(sc['wt_active_support_ms'])+p*float(sc['sca_active_support_ms'])
    return {'b':b,'r':r,'s':s,'kappa_I':k},J,win,pdri


def _scenario_path_job(sc,cfg,family):
    n=int(cfg['transition']['n_grid'])
    ps=np.linspace(0.0,1.0,n)
    rows=[]
    rb_guess=float(sc['wt_rheobase_J_endpoint'])
    for idx,p in enumerate(ps):
        theta,J,window,pdrive=_interp_theta(sc,float(p),family)
        x,y,z,ok=pre_relax(theta,cfg,dt_ms=float(cfg['simulation']['dt_ms']))
        if ok:
            rb=refine_rheobase(theta,cfg,guess=rb_guess,pre_state=(x,y,z))
            if np.isfinite(rb['rheobase_J']): rb_guess=float(rb['rheobase_J'])
            met=supported_metrics(theta,J,window,cfg,pre_state=(x,y,z))
        else:
            rb={'rheobase_J':np.nan,'status':'PRE_RELAX_FAIL','iterations':0}
            met={'spike_count':0,'support_rate_hz':np.nan,'mean_isi_ms':np.nan,'occupancy_fraction':np.nan,'first_spike_ms':np.nan,'simulation_ok':False}
        rows.append({
            'scenario_id':int(sc['scenario_id']),'biological_pair_key':sc['biological_pair_key'],
            'wt_cell_id':sc['wt_cell_id'],'sca_cell_id':sc['sca_cell_id'],
            'wt_solution_key':sc['wt_solution_key'],'sca_solution_key':sc['sca_solution_key'],
            'path_family':family,'grid_index':idx,'path_progress':float(p),'drive_progress':float(pdrive),
            'b':theta['b'],'r':theta['r'],'s':theta['s'],'kappa_I':theta['kappa_I'],
            'J_protocol':float(J),'active_support_ms':float(window),
            'rheobase_J':rb['rheobase_J'],'rheobase_status':rb['status'],'rheobase_iterations':rb['iterations'],
            'J_over_rheobase':float(J/rb['rheobase_J']) if np.isfinite(rb['rheobase_J']) and rb['rheobase_J']>0 else np.nan,
            **met,
            'within_pair_support_weight':float(sc['within_pair_support_weight']),
            'scenario_weight':float(sc['scenario_weight']),
            'biological_pair_weight':float(sc['biological_pair_weight']),
        })
    return rows


def _checkpoint_path(out,scenario_id,family):
    return out/'checkpoints'/f'scenario_{int(scenario_id):04d}__{family}.json.gz'


def _write_checkpoint(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    with gzip.open(tmp,'wt',encoding='utf-8') as f: json.dump(rows,f,separators=(',',':'))
    os.replace(tmp,path)


def _read_checkpoint(path):
    with gzip.open(path,'rt',encoding='utf-8') as f: return json.load(f)


def build_protocol_endpoints(support,cells,cfg):
    rows=[]
    for _,r in support.iterrows():
        theta=theta_from_row(r)
        x,y,z,ok=pre_relax(theta,cfg,dt_ms=float(cfg['simulation']['dt_ms']))
        if ok:
            rr=refine_rheobase(theta,cfg,guess=float(r.rheobase_J),pre_state=(x,y,z))
            met=supported_metrics(theta,float(r.J_q75),float(r.active_support_q75_ms),cfg,pre_state=(x,y,z))
        else:
            rr={'rheobase_J':np.nan,'status':'PRE_RELAX_FAIL','iterations':0}
            met={'spike_count':0,'support_rate_hz':np.nan,'mean_isi_ms':np.nan,'occupancy_fraction':np.nan,'first_spike_ms':np.nan,'simulation_ok':False}
        rows.append({
            'group':r.group,'cell_id':r.cell_id,'animal_id':r.animal_id,'solution_key':r.solution_key,'source':r.source,
            'b':r.b,'r':r.r,'s':r.s,'kappa_I':r.kappa_I,
            'stored_rheobase_J':r.rheobase_J,'protocol_rheobase_J':rr['rheobase_J'],'rheobase_status':rr['status'],
            'rheobase_relative_error':abs(rr['rheobase_J']-r.rheobase_J)/max(abs(r.rheobase_J),1e-12) if np.isfinite(rr['rheobase_J']) else np.nan,
            'J_q75':r.J_q75,'active_support_q75_ms':r.active_support_q75_ms,
            **met,'within_cell_support_weight':r.within_cell_support_weight,'group_support_weight':r.group_support_weight
        })
    ep=pd.DataFrame(rows)
    meta=cells[['cell_id','core_q75_secure','exp_q75_firing_rate_hz','exp_q75_mean_isi_ms','rheobase_J_best']].copy()
    ep=ep.merge(meta,on='cell_id',how='left',validate='many_to_one')
    return ep


def _safe_log10(x):
    x=pd.to_numeric(x,errors='coerce').to_numpy(float)
    out=np.full(len(x),np.nan); m=np.isfinite(x)&(x>0); out[m]=np.log10(x[m]); return out


def build_references(endpoint_protocol,cfg):
    best=endpoint_protocol[endpoint_protocol.source.eq('best')].copy()
    best['log10_rheobase']=_safe_log10(best.protocol_rheobase_J)
    best['log10_rate']=_safe_log10(best.support_rate_hz)
    best['log10_isi']=_safe_log10(best.mean_isi_ms)
    secure=lambda d: d.core_q75_secure.fillna(False).astype(bool)
    rate_ref=projection_reference(best,['log10_rheobase','log10_rate'],secure,
                                  float(cfg['staging']['wt_exit_quantile']),float(cfg['staging']['sca3_entry_quantile']))
    isi_ref=projection_reference(best,['log10_rheobase','log10_isi'],secure,
                                 float(cfg['staging']['wt_exit_quantile']),float(cfg['staging']['sca3_entry_quantile']))
    return best,rate_ref,isi_ref


def reference_tables(rate_ref,isi_ref):
    rows=[]; trans=[]
    for name,ref in [('rate',rate_ref),('isi',isi_ref)]:
        rows.append({'projection':name,'reference_subset':ref['subset'],'centroid_distance':ref['centroid_distance'],
                     'wt_exit_A_threshold':ref['wt_exit_A_threshold'],'sca3_entry_A_threshold':ref['sca3_entry_A_threshold'],
                     'cloud_overlap':ref['cloud_overlap'],'corridor_radius_q90':ref['corridor_radius_q90'],
                     'wt_centroid_0':ref['cwt'][0],'wt_centroid_1':ref['cwt'][1],
                     'sca3_centroid_0':ref['csc'][0],'sca3_centroid_1':ref['csc'][1]})
        t=ref['transform'].copy(); t.insert(0,'projection',name); trans.append(t)
    return pd.DataFrame(rows),pd.concat(trans,ignore_index=True)


def endpoint_validation(endpoint_protocol,rate_ref,isi_ref,cfg):
    b=endpoint_protocol[endpoint_protocol.source.eq('best')].copy()
    out=[]
    # Group directions in native coordinates.
    for metric,expcol,modcol in [('rate','exp_q75_firing_rate_hz','support_rate_hz'),('isi','exp_q75_mean_isi_ms','mean_isi_ms')]:
        for g in ['WT','SCA3']:
            x=b[b.group.eq(g)]
            out.append({'kind':'group_median','metric':metric,'group':g,'experiment':float(np.nanmedian(x[expcol])),
                        'protocol_model':float(np.nanmedian(x[modcol])),'n_cells':int(x[modcol].notna().sum())})
        wt=b[b.group.eq('WT')]; sc=b[b.group.eq('SCA3')]
        exp_delta=float(np.nanmedian(sc[expcol])-np.nanmedian(wt[expcol]))
        mod_delta=float(np.nanmedian(sc[modcol])-np.nanmedian(wt[modcol]))
        out.append({'kind':'direction','metric':metric,'group':'SCA3_minus_WT','experiment':exp_delta,'protocol_model':mod_delta,
                    'direction_agrees':bool(exp_delta*mod_delta>0),'n_cells':int(b[modcol].notna().sum())})
    return pd.DataFrame(out)


def project_all_paths(paths,rate_ref,isi_ref):
    d=paths.copy()
    d['log10_rheobase']=_safe_log10(d.rheobase_J)
    d['log10_rate']=_safe_log10(d.support_rate_hz)
    d['log10_isi']=_safe_log10(d.mean_isi_ms)
    pr=project_points(d[['log10_rheobase','log10_rate']].copy(),rate_ref,['log10_rheobase','log10_rate'])
    pi=project_points(d[['log10_rheobase','log10_isi']].copy(),isi_ref,['log10_rheobase','log10_isi'])
    d['A_rate']=pr.A.to_numpy(); d['orth_rate']=pr.orthogonal_distance.to_numpy()
    d['A_isi']=pi.A.to_numpy(); d['orth_isi']=pi.orthogonal_distance.to_numpy()
    return d


def _spearman_xy(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); m=np.isfinite(x)&np.isfinite(y); x=x[m]; y=y[m]
    if len(x)<3:return np.nan
    rx=pd.Series(x).rank(method='average').to_numpy(); ry=pd.Series(y).rank(method='average').to_numpy()
    if np.std(rx)==0 or np.std(ry)==0:return np.nan
    return float(np.corrcoef(rx,ry)[0,1])


def stage_markers(paths,scenarios,rate_ref,isi_ref,cfg):
    persist=int(cfg['staging']['persistence_points'])
    sched=cfg['transition']['schedule_morphs']
    smap=scenarios.set_index('scenario_id')
    rows=[]
    for (sid,fam),g in paths.groupby(['scenario_id','path_family'],sort=False):
        g=g.sort_values('path_progress'); p=g.path_progress.to_numpy(float)
        row={'scenario_id':int(sid),'path_family':fam}
        sc=smap.loc[int(sid)]
        for k in ['biological_pair_key','wt_cell_id','sca_cell_id','wt_solution_key','sca_solution_key','within_pair_support_weight','scenario_weight','biological_pair_weight']:
            row[k]=sc[k]
        for name,ref,col,orth in [('rate',rate_ref,'A_rate','orth_rate'),('isi',isi_ref,'A_isi','orth_isi')]:
            A=g[col].to_numpy(float)
            row[f'wt_exit_p_{name}']=persistent_crossing(p,A,float(ref['wt_exit_A_threshold']),persist)
            row[f'balance_p_{name}']=persistent_crossing(p,A,0.5,persist)
            row[f'sca3_entry_p_{name}']=persistent_crossing(p,A,float(ref['sca3_entry_A_threshold']),persist)
            row[f'start_A_{name}']=A[0] if len(A) else np.nan; row[f'end_A_{name}']=A[-1] if len(A) else np.nan
            row[f'start_inside_wt_{name}']=bool(np.isfinite(A[0]) and A[0]<=ref['wt_exit_A_threshold']) if len(A) else False
            row[f'end_inside_sca3_{name}']=bool(np.isfinite(A[-1]) and A[-1]>=ref['sca3_entry_A_threshold']) if len(A) else False
            row[f'spearman_p_A_{name}']=_spearman_xy(p,A)
            row[f'reversals_{name}']=count_reversals(A,float(cfg['staging']['reversal_tolerance']))
            oo=pd.to_numeric(g[orth],errors='coerce').to_numpy(float)
            row[f'max_orth_{name}']=float(np.nanmax(oo)) if np.isfinite(oo).any() else np.nan
            rad=float(ref['corridor_radius_q90'])
            row[f'fraction_in_corridor_{name}']=float(np.mean(oo[np.isfinite(oo)]<=rad)) if np.isfinite(oo).any() and np.isfinite(rad) else np.nan
        row['silent_fraction']=float(np.mean(pd.to_numeric(g.support_rate_hz,errors='coerce').fillna(0).to_numpy()<=0))
        row['isi_missing_fraction']=float(g.mean_isi_ms.isna().mean())
        for stage in ['wt_exit','balance','sca3_entry']:
            a=row[f'{stage}_p_rate']; b=row[f'{stage}_p_isi']
            row[f'{stage}_projection_gap']=abs(a-b) if np.isfinite(a) and np.isfinite(b) else np.nan
            row[f'{stage}_consensus_p']=0.5*(a+b) if np.isfinite(a) and np.isfinite(b) and abs(a-b)<=float(cfg['staging']['projection_consensus_tolerance']) else np.nan
            for morph in sched:
                for proj in ['rate','isi','consensus']:
                    key=f'{stage}_p_{proj}' if proj!='consensus' else f'{stage}_consensus_p'
                    val=row.get(key,np.nan)
                    row[f'{stage}_u_{proj}_{morph}']=inverse_morph(val,morph) if np.isfinite(val) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def weighted_quantile(values,weights,q):
    v=np.asarray(values,float); w=np.asarray(weights,float); m=np.isfinite(v)&np.isfinite(w)&(w>0); v=v[m]; w=w[m]
    if not len(v): return np.nan
    o=np.argsort(v); v=v[o]; w=w[o]; c=np.cumsum(w); c=c/c[-1]
    return float(np.interp(q,c,v))


def summarize_pairs(markers,cells,cfg):
    secure=set(cells.loc[cells.core_q75_secure.fillna(False).astype(bool),'cell_id'].astype(str))
    rows=[]
    metrics=[f'{s}_p_{p}' for s in ['wt_exit','balance','sca3_entry'] for p in ['rate','isi']]+[f'{s}_consensus_p' for s in ['wt_exit','balance','sca3_entry']]
    for (pair,fam),g in markers.groupby(['biological_pair_key','path_family']):
        w=g.within_pair_support_weight.to_numpy(float)
        r={'biological_pair_key':pair,'path_family':fam,'wt_cell_id':g.wt_cell_id.iloc[0],'sca_cell_id':g.sca_cell_id.iloc[0],
           'both_core_secure':bool(str(g.wt_cell_id.iloc[0]) in secure and str(g.sca_cell_id.iloc[0]) in secure),
           'n_support_scenarios':len(g)}
        for m in metrics:
            r[m+'_weighted_median']=weighted_quantile(g[m],w,.5)
            r[m+'_weighted_q25']=weighted_quantile(g[m],w,.25)
            r[m+'_weighted_q75']=weighted_quantile(g[m],w,.75)
            r[m+'_support_fraction']=float(np.sum(w[np.isfinite(g[m].to_numpy(float))])/np.sum(w)) if np.sum(w)>0 else np.nan
        r['silent_fraction_weighted']=float(np.average(g.silent_fraction,weights=w))
        r['projection_gap_balance_weighted_median']=weighted_quantile(g.balance_projection_gap,w,.5)
        rows.append(r)
    pairdf=pd.DataFrame(rows)
    overall=[]
    for fam,g in pairdf.groupby('path_family'):
        for subset_name,x in [('all_pairs',g),('core_secure_pairs',g[g.both_core_secure])]:
            for m in metrics:
                vals=pd.to_numeric(x[m+'_weighted_median'],errors='coerce').to_numpy(float); vals=vals[np.isfinite(vals)]
                overall.append({'path_family':fam,'subset':subset_name,'metric':m,'n_biological_pairs_total':len(x),'n_pairs_with_marker':len(vals),
                                'median':float(np.median(vals)) if len(vals) else np.nan,
                                'q25':float(np.quantile(vals,.25)) if len(vals) else np.nan,'q75':float(np.quantile(vals,.75)) if len(vals) else np.nan})
    return pairdf,pd.DataFrame(overall)



def summarize_curves(paths):
    metrics=['A_rate','A_isi','rheobase_J','support_rate_hz','mean_isi_ms','J_over_rheobase','orth_rate','orth_isi']
    pair_rows=[]
    for (pair,fam,p),g in paths.groupby(['biological_pair_key','path_family','path_progress']):
        w=pd.to_numeric(g.within_pair_support_weight,errors='coerce').to_numpy(float)
        row={'biological_pair_key':pair,'path_family':fam,'path_progress':p,'wt_cell_id':g.wt_cell_id.iloc[0],'sca_cell_id':g.sca_cell_id.iloc[0]}
        for m in metrics:
            row[m+'_weighted_median']=weighted_quantile(g[m],w,.5)
        pair_rows.append(row)
    pair=pd.DataFrame(pair_rows)
    ens=[]
    for (fam,p),g in pair.groupby(['path_family','path_progress']):
        row={'path_family':fam,'path_progress':p,'n_biological_pairs':len(g)}
        for m in metrics:
            x=pd.to_numeric(g[m+'_weighted_median'],errors='coerce').to_numpy(float); x=x[np.isfinite(x)]
            row[m+'_median']=float(np.median(x)) if len(x) else np.nan
            row[m+'_q25']=float(np.quantile(x,.25)) if len(x) else np.nan
            row[m+'_q75']=float(np.quantile(x,.75)) if len(x) else np.nan
            row[m+'_n']=len(x)
        ens.append(row)
    return pair,pd.DataFrame(ens)

def validate(cfg):
    support,cells,anchors,transform,geometry=load_inputs(cfg)
    support=build_support_with_anchors(support,anchors)
    scenarios=enumerate_scenarios(support)
    out={
        'version':'1.0.0','primary_cells':int(cells.cell_id.nunique()),'WT_cells':int((cells.group=='WT').sum()),'SCA3_cells':int((cells.group=='SCA3').sum()),
        'support_states':int(len(support)),'WT_support_states':int((support.group=='WT').sum()),'SCA3_support_states':int((support.group=='SCA3').sum()),
        'biological_pairs':int(scenarios.biological_pair_key.nunique()),'support_pair_scenarios':int(len(scenarios)),
        'scenario_weight_sum':float(scenarios.scenario_weight.sum()),'all_q75_supported':bool(anchors.q75_supported.all()),
        'path_families':list(cfg['transition']['path_families']),'n_grid':int(cfg['transition']['n_grid']),
        'frozen_sha256':frozen_manifest(cfg['data']['root'])
    }
    if out['primary_cells']!=18 or out['WT_cells']!=12 or out['SCA3_cells']!=6: raise ValueError('unexpected primary cohort')
    if out['support_states']!=64 or out['support_pair_scenarios']!=988: raise ValueError('unexpected support-state counts')
    if abs(out['scenario_weight_sum']-1.0)>1e-9: raise ValueError('scenario weights do not sum to one')
    if not out['all_q75_supported']: raise ValueError('q=.75 not supported for all primary cells')
    return out


def run_all(cfg,resume=True):
    out=Path(cfg['output']['dir']); out.mkdir(parents=True,exist_ok=True)
    support,cells,anchors,exp_transform,exp_geometry=load_inputs(cfg)
    support=build_support_with_anchors(support,anchors)
    scenarios=enumerate_scenarios(support)
    mode=str(cfg['transition'].get('scenario_mode','all_support'))
    if mode=='best_only': scenarios=scenarios[(scenarios.wt_source=='best')&(scenarios.sca_source=='best')].copy()
    maxs=cfg['transition'].get('max_scenarios')
    if maxs is not None: scenarios=scenarios.head(int(maxs)).copy()
    endpoint_protocol=build_protocol_endpoints(support,cells,cfg)
    best,rate_ref,isi_ref=build_references(endpoint_protocol,cfg)
    ref_summary,ref_transform=reference_tables(rate_ref,isi_ref)
    epval=endpoint_validation(endpoint_protocol,rate_ref,isi_ref,cfg)
    endpoint_protocol.to_csv(out/'transition_protocol_endpoint_states.csv',index=False)
    ref_summary.to_csv(out/'transition_projection_reference.csv',index=False)
    ref_transform.to_csv(out/'transition_projection_transform.csv',index=False)
    epval.to_csv(out/'endpoint_protocol_validation.csv',index=False)
    scenarios.to_csv(out/'transition_pair_scenarios.csv',index=False)

    families=list(cfg['transition']['path_families']); jobs=[]; allrows=[]
    for _,scrow in scenarios.iterrows():
        sc=scrow.to_dict()
        for fam in families:
            cp=_checkpoint_path(out,sc['scenario_id'],fam)
            if resume and cp.exists(): allrows.extend(_read_checkpoint(cp))
            else: jobs.append((sc,fam,cp))
    workers=int(cfg['parallel']['workers'])
    if jobs:
        if workers<=1:
            for sc,fam,cp in jobs:
                rr=_scenario_path_job(sc,cfg,fam); _write_checkpoint(cp,rr); allrows.extend(rr)
        else:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                fut={ex.submit(_scenario_path_job,sc,cfg,fam):(cp,sc['scenario_id'],fam) for sc,fam,cp in jobs}
                for f in as_completed(fut):
                    cp,sid,fam=fut[f]; rr=f.result(); _write_checkpoint(cp,rr); allrows.extend(rr)
    paths=pd.DataFrame(allrows)
    if paths.empty: raise RuntimeError('no transition paths produced')
    paths=paths.sort_values(['scenario_id','path_family','grid_index']).reset_index(drop=True)
    paths=project_all_paths(paths,rate_ref,isi_ref)
    paths.to_csv(out/'transition_paths.csv',index=False)
    markers=stage_markers(paths,scenarios,rate_ref,isi_ref,cfg)
    markers.to_csv(out/'scenario_stage_markers.csv',index=False)
    pair_summary,overall=summarize_pairs(markers,cells,cfg)
    pair_summary.to_csv(out/'biological_pair_stage_summary.csv',index=False)
    overall.to_csv(out/'path_family_stage_summary.csv',index=False)
    pair_curve,ensemble_curve=summarize_curves(paths)
    pair_curve.to_csv(out/'biological_pair_curve_summary.csv',index=False)
    ensemble_curve.to_csv(out/'ensemble_curve_summary.csv',index=False)

    summary={
        'version':'1.0.0','analysis':'WT_to_SCA3_transition_ensemble',
        'scenario_mode':mode,'n_scenarios':int(scenarios.scenario_id.nunique()),'n_path_families':len(families),'n_grid':int(cfg['transition']['n_grid']),
        'n_path_rows':int(len(paths)),'biological_pairs':int(scenarios.biological_pair_key.nunique()),
        'support_pair_scenarios_available_full':988,'path_families':families,
        'rate_reference_centroid_distance':rate_ref['centroid_distance'],'isi_reference_centroid_distance':isi_ref['centroid_distance'],
        'rate_cloud_overlap':rate_ref['cloud_overlap'],'isi_cloud_overlap':isi_ref['cloud_overlap'],
        'silent_path_row_fraction':float(np.mean(pd.to_numeric(paths.support_rate_hz,errors='coerce').fillna(0).to_numpy()<=0)),
        'missing_isi_row_fraction':float(paths.mean_isi_ms.isna().mean()),
        'endpoint_max_rheobase_relative_error':float(np.nanmax(endpoint_protocol.rheobase_relative_error)),
        'frozen_sha256':frozen_manifest(cfg['data']['root'])
    }
    (out/'RUN_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True),encoding='utf-8')
    return summary
