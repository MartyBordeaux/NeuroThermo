from __future__ import annotations
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, math, os
import numpy as np
import pandas as pd

from .v12io import V12Results
from .model import pre_relax, refine_rheobase, supported_metrics
from .geometry import load_geometry, project_scalar
from .mathutils import persistent_crossing, weighted_finite_quantile, crossing_support_weight, bilinear_track

STAGES=[('WT_exit','wt_exit'),('balance','balance'),('SCA3_entry','sca3_entry')]
PROJS=[('isi','A_isi'),('active','A_active')]
MODES=['combined','kappa_only','J_only']


def _resolve(cfg):
    return V12Results(cfg.get('input_combined',{}).get('directory'), cfg.get('input_combined',{}).get('archive'))


def _geom(root):
    f=Path(root)/'frozen'
    refs=pd.read_csv(f/'transition_projection_reference_v1_1.csv')
    trans=pd.read_csv(f/'transition_projection_transform_v1_1.csv')
    return load_geometry(refs,trans)


def _thresholds(run):
    return {
      'isi':{'wt_exit':float(run['primary_ISI_boundaries']['WT_exit']),'balance':0.5,'sca3_entry':float(run['primary_ISI_boundaries']['SCA3_entry'])},
      'active':{'wt_exit':float(run['secondary_active_boundaries']['WT_exit']),'balance':0.5,'sca3_entry':float(run['secondary_active_boundaries']['SCA3_entry'])},
    }


def _interp_theta(sc,pi,pk):
    pi=float(pi); pk=float(pk)
    return {
      'b':(1-pi)*float(sc['wt_b'])+pi*float(sc['sca_b']),
      's':(1-pi)*float(sc['wt_s'])+pi*float(sc['sca_s']),
      'r':math.exp((1-pi)*math.log(float(sc['wt_r']))+pi*math.log(float(sc['sca_r']))),
      'kappa_I':math.exp((1-pk)*math.log(float(sc['wt_kappa_I']))+pk*math.log(float(sc['sca_kappa_I']))),
    }


def _interp_J(sc,pj):
    pj=float(pj)
    return (1-pj)*float(sc['wt_J_q75'])+pj*float(sc['sca_J_q75'])


def _window(sc,pi):
    return (1-float(pi))*float(sc['wt_active_support_ms'])+float(pi)*float(sc['sca_active_support_ms'])


def _safe_log10(x):
    return math.log10(float(x)) if np.isfinite(x) and float(x)>0 else np.nan


def _checkpoint_path(out,mode,sid):
    return out/'checkpoints'/mode/f'scenario_{int(sid):04d}.csv.gz'


def _write_checkpoint(path,df):
    path.parent.mkdir(parents=True,exist_ok=True); tmp=Path(str(path)+'.tmp')
    df.to_csv(tmp,index=False,compression='gzip'); os.replace(tmp,path)


def _read_checkpoint(path): return pd.read_csv(path,compression='gzip')


def _simulate_mode(sc,cfg,geom,mode):
    ni=int(cfg['surface']['n_intrinsic']); nq=int(cfg['surface']['n_component'])
    pis=np.linspace(0,1,ni); qs=np.linspace(0,1,nq); rows=[]
    prev_rb=np.full(nq,np.nan)
    for ii,pi in enumerate(pis):
        theta_pre=_interp_theta(sc,pi,pi)
        x,y,z,okpre=pre_relax(theta_pre,cfg,dt_ms=float(cfg['simulation']['dt_ms']))
        pre=(x,y,z) if okpre else None
        current=np.full(nq,np.nan); prevq=np.nan
        # J-only: theta and rheobase are identical across q at fixed intrinsic progress.
        fixed_rb=None
        if okpre and mode=='J_only':
            theta_fixed=_interp_theta(sc,pi,pi)
            guess=(1-pi)*float(sc['wt_rheobase_J_endpoint'])+pi*float(sc['sca_rheobase_J_endpoint'])
            fixed_rb=refine_rheobase(theta_fixed,cfg,guess=guess,pre_state=pre)
        for jj,q in enumerate(qs):
            if mode=='kappa_only': pk,pj=q,pi
            elif mode=='J_only': pk,pj=pi,q
            else: raise ValueError(mode)
            theta=_interp_theta(sc,pi,pk); J=_interp_J(sc,pj); window=_window(sc,pi)
            if not okpre:
                rb={'rheobase_J':np.nan,'status':'PRE_RELAX_FAIL','iterations':0}
                met={'spike_count':0,'support_rate_hz':np.nan,'mean_isi_ms':np.nan,'occupancy_fraction':np.nan,'first_spike_ms':np.nan,'simulation_ok':False}
            else:
                if mode=='J_only':
                    rb=fixed_rb
                else:
                    if np.isfinite(prevq): guess=prevq
                    elif np.isfinite(prev_rb[jj]): guess=prev_rb[jj]
                    else: guess=(1-pi)*float(sc['wt_rheobase_J_endpoint'])+pi*float(sc['sca_rheobase_J_endpoint'])
                    rb=refine_rheobase(theta,cfg,guess=guess,pre_state=pre)
                    if np.isfinite(rb['rheobase_J']): prevq=float(rb['rheobase_J']); current[jj]=prevq
                met=supported_metrics(theta,J,window,cfg,pre_state=pre)
            lrb=_safe_log10(rb['rheobase_J']); lisi=_safe_log10(met['mean_isi_ms']); lact=_safe_log10(met['support_rate_hz'])
            Ai,oi=project_scalar(lrb,lisi,geom['isi']); Aa,oa=project_scalar(lrb,lact,geom['active'])
            rows.append({
              'scenario_id':int(sc['scenario_id']),'biological_pair_key':sc['biological_pair_key'],'wt_cell_id':sc['wt_cell_id'],'sca_cell_id':sc['sca_cell_id'],
              'mode':mode,'intrinsic_index':ii,'component_index':jj,'p_intrinsic':float(pi),'p_component':float(q),'p_kappa':float(pk),'p_J':float(pj),
              'b':theta['b'],'r':theta['r'],'s':theta['s'],'kappa_I':theta['kappa_I'],'J_protocol':J,'active_support_ms':window,
              'rheobase_J':rb['rheobase_J'],'rheobase_status':rb['status'],'rheobase_iterations':rb['iterations'],
              'spike_count':met['spike_count'],'active_support_rate_hz':met['support_rate_hz'],'mean_isi_ms':met['mean_isi_ms'],'occupancy_fraction':met['occupancy_fraction'],'first_spike_ms':met['first_spike_ms'],'simulation_ok':met['simulation_ok'],
              'A_isi':Ai,'orth_isi':oi,'A_active':Aa,'orth_active':oa,
              'within_pair_support_weight':float(sc['within_pair_support_weight']),'scenario_weight':float(sc['scenario_weight']),'biological_pair_weight':float(sc['biological_pair_weight'])})
        prev_rb=current
    return pd.DataFrame(rows)


def _scenario_job(sc,cfg,geom):
    return int(sc['scenario_id']), _simulate_mode(sc,cfg,geom,'kappa_only'), _simulate_mode(sc,cfg,geom,'J_only')


def _sample_combined(d,pis,qs,cols=('A_isi','A_active')):
    opi=np.sort(d.p_intrinsic.unique()); opd=np.sort(d.p_drive.unique())
    rows=[]
    X,Y=np.meshgrid(pis,qs,indexing='ij'); xf=X.ravel(); yf=Y.ravel()
    mats={}
    for col in cols:
        M=d.pivot(index='p_intrinsic',columns='p_drive',values=col).reindex(index=opi,columns=opd).to_numpy(float)
        mats[col]=bilinear_track(opi,opd,M,xf,yf).reshape(len(pis),len(qs))
    for i,pi in enumerate(pis):
        for j,q in enumerate(qs):
            rows.append({'p_intrinsic':float(pi),'p_component':float(q),'A_isi':mats['A_isi'][i,j],'A_active':mats['A_active'][i,j]})
    return pd.DataFrame(rows)


def _scenario_crossings(mode_df,sc,secure,thresholds,persistence,mode):
    rows=[]; common={'scenario_id':int(sc.scenario_id),'biological_pair_key':sc.biological_pair_key,'wt_cell_id':sc.wt_cell_id,'sca_cell_id':sc.sca_cell_id,'both_core_secure':bool(secure),'within_pair_support_weight':float(sc.within_pair_support_weight),'mode':mode}
    for pi,g in mode_df.groupby('p_intrinsic',sort=True):
        g=g.sort_values('p_component'); x=g.p_component.to_numpy(float)
        for proj,col in PROJS:
            y=pd.to_numeric(g[col],errors='coerce').to_numpy(float)
            for stage,key in STAGES:
                rows.append({**common,'p_intrinsic':float(pi),'projection':proj,'stage':stage,'A_threshold':thresholds[proj][key],'crossing_value':persistent_crossing(x,y,thresholds[proj][key],persistence)})
    return pd.DataFrame(rows)


def _aggregate_pair_crossings(frames,weights):
    if not frames:return pd.DataFrame()
    keys=['biological_pair_key','wt_cell_id','sca_cell_id','both_core_secure','mode','p_intrinsic','projection','stage','A_threshold']
    base=frames[0][keys].copy(); V=np.stack([f.crossing_value.to_numpy(float) for f in frames]); w=np.asarray(weights,float)
    base['n_support_scenarios']=len(frames)
    base['crossing_support_weight']=[crossing_support_weight(V[:,j],w) for j in range(V.shape[1])]
    base['q25_weighted']=[weighted_finite_quantile(V[:,j],w,.25) for j in range(V.shape[1])]
    base['median_weighted']=[weighted_finite_quantile(V[:,j],w,.5) for j in range(V.shape[1])]
    base['q75_weighted']=[weighted_finite_quantile(V[:,j],w,.75) for j in range(V.shape[1])]
    return base


def _ensemble_crossings(pair):
    rows=[]; keys=['mode','p_intrinsic','projection','stage','A_threshold']
    for subset,d in [('all_pairs',pair),('core_secure_pairs',pair[pair.both_core_secure.astype(bool)])]:
        for k,g in d.groupby(keys,sort=True):
            x=pd.to_numeric(g.median_weighted,errors='coerce').dropna().to_numpy(float); sw=pd.to_numeric(g.crossing_support_weight,errors='coerce').dropna().to_numpy(float)
            maj=g[pd.to_numeric(g.crossing_support_weight,errors='coerce')>=.5]; xm=pd.to_numeric(maj.median_weighted,errors='coerce').dropna().to_numpy(float)
            rows.append({'subset':subset,**dict(zip(keys,k)),'n_biological_pairs_total':len(g),'n_pairs_with_marker':len(x),'pair_marker_fraction':len(x)/len(g) if len(g) else np.nan,'n_pairs_support_ge_0_5':len(xm),'pair_majority_support_fraction':len(xm)/len(g) if len(g) else np.nan,'median_pair_crossing_support_weight':float(np.median(sw)) if len(sw) else np.nan,'median':float(np.median(x)) if len(x) else np.nan,'q25':float(np.quantile(x,.25)) if len(x) else np.nan,'q75':float(np.quantile(x,.75)) if len(x) else np.nan,'median_majority_support':float(np.median(xm)) if len(xm) else np.nan})
    return pd.DataFrame(rows)


def _scenario_effects(comb,kdf,jdf,sc,secure,thresholds,band):
    keys=['p_intrinsic','p_component']; z=comb.merge(kdf[keys+['A_isi','A_active']],on=keys,suffixes=('_combined','_kappa')).merge(jdf[keys+['A_isi','A_active']],on=keys)
    z=z.rename(columns={'A_isi':'A_isi_J','A_active':'A_active_J'})
    # Coupled baseline at p_component = p_intrinsic from combined surface, interpolated row-wise.
    base_i={}; base_a={}
    for pi,g in comb.groupby('p_intrinsic'):
        g=g.sort_values('p_component'); base_i[pi]=float(np.interp(pi,g.p_component,g.A_isi)); base_a[pi]=float(np.interp(pi,g.p_component,g.A_active))
    z['A_isi_base']=z.p_intrinsic.map(base_i); z['A_active_base']=z.p_intrinsic.map(base_a)
    for proj in ['isi','active']:
        C=z[f'A_{proj}_combined']-z[f'A_{proj}_base']; K=z[f'A_{proj}_kappa']-z[f'A_{proj}_base']; J=z[f'A_{proj}_J']-z[f'A_{proj}_base']
        z[f'delta_{proj}_combined']=C; z[f'delta_{proj}_kappa']=K; z[f'delta_{proj}_J']=J; z[f'interaction_{proj}']=C-K-J
    z['scenario_id']=int(sc.scenario_id); z['biological_pair_key']=sc.biological_pair_key; z['both_core_secure']=bool(secure); z['within_pair_support_weight']=float(sc.within_pair_support_weight)
    stage=[]
    for proj in ['isi','active']:
        for st,key in STAGES:
            thr=thresholds[proj][key]; mask=np.isfinite(z[f'A_{proj}_combined'])&(np.abs(z[f'A_{proj}_combined']-thr)<=band)
            row={'scenario_id':int(sc.scenario_id),'biological_pair_key':sc.biological_pair_key,'both_core_secure':bool(secure),'within_pair_support_weight':float(sc.within_pair_support_weight),'projection':proj,'stage':st,'A_threshold':thr,'n_grid_points':int(mask.sum())}
            for m in ['combined','kappa','J']:
                a=pd.to_numeric(z.loc[mask,f'delta_{proj}_{m}'],errors='coerce').dropna().to_numpy(float); row[f'median_abs_delta_{m}']=float(np.median(np.abs(a))) if len(a) else np.nan
            a=pd.to_numeric(z.loc[mask,f'interaction_{proj}'],errors='coerce').dropna().to_numpy(float); row['median_abs_interaction']=float(np.median(np.abs(a))) if len(a) else np.nan; row['median_signed_interaction']=float(np.median(a)) if len(a) else np.nan
            stage.append(row)
    return z,stage


def _weighted_matrix(V,w):
    V=np.asarray(V,float); w=np.asarray(w,float); n,ni,nj=V.shape; X=V.reshape(n,-1); out=np.full(X.shape[1],np.nan)
    for j in range(X.shape[1]): out[j]=weighted_finite_quantile(X[:,j],w,.5)
    return out.reshape(ni,nj)


def _pair_effect_surface(pair,scenario_effects,weights):
    pis=np.sort(scenario_effects[0].p_intrinsic.unique()); qs=np.sort(scenario_effects[0].p_component.unique())
    metrics=['delta_isi_combined','delta_isi_kappa','delta_isi_J','interaction_isi','delta_active_combined','delta_active_kappa','delta_active_J','interaction_active']
    mats={}
    for m in metrics:
        mats[m]=_weighted_matrix(np.stack([d.pivot(index='p_intrinsic',columns='p_component',values=m).reindex(index=pis,columns=qs).to_numpy(float) for d in scenario_effects]),weights)
    rows=[]
    secure=bool(scenario_effects[0].both_core_secure.iloc[0])
    for i,pi in enumerate(pis):
        for j,q in enumerate(qs):
            row={'biological_pair_key':pair,'both_core_secure':secure,'p_intrinsic':pi,'p_component':q}
            for m in metrics: row[m+'_weighted_median']=mats[m][i,j]
            rows.append(row)
    return rows


def _ensemble_effect(pair_effect):
    rows=[]; metrics=[c for c in pair_effect if c.endswith('_weighted_median')]
    for subset,d in [('all_pairs',pair_effect),('core_secure_pairs',pair_effect[pair_effect.both_core_secure.astype(bool)])]:
        for (pi,q),g in d.groupby(['p_intrinsic','p_component'],sort=True):
            row={'subset':subset,'p_intrinsic':pi,'p_component':q,'n_biological_pairs':g.biological_pair_key.nunique()}
            for m in metrics:
                x=pd.to_numeric(g[m],errors='coerce').dropna().to_numpy(float); stem=m.replace('_weighted_median',''); row[stem+'_median']=float(np.median(x)) if len(x) else np.nan; row[stem+'_q25']=float(np.quantile(x,.25)) if len(x) else np.nan; row[stem+'_q75']=float(np.quantile(x,.75)) if len(x) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def _aggregate_stage_effects(pair_stage):
    # First aggregate scenario values within biological pair, then pairs equally.
    metrics=['median_abs_delta_combined','median_abs_delta_kappa','median_abs_delta_J','median_abs_interaction','median_signed_interaction']
    prows=[]
    for (pair,proj,stage,thr),g in pair_stage.groupby(['biological_pair_key','projection','stage','A_threshold']):
        w=pd.to_numeric(g.within_pair_support_weight,errors='coerce').to_numpy(float); row={'biological_pair_key':pair,'both_core_secure':bool(g.both_core_secure.iloc[0]),'projection':proj,'stage':stage,'A_threshold':thr}
        for m in metrics: row[m+'_weighted_median']=weighted_finite_quantile(pd.to_numeric(g[m],errors='coerce'),w,.5)
        den=np.sum(w[np.isfinite(w)&(w>0)]); row['scenario_support_with_band']=float(np.sum(w[(pd.to_numeric(g.n_grid_points,errors='coerce').to_numpy(float)>0)&np.isfinite(w)&(w>0)])/den) if den>0 else np.nan
        prows.append(row)
    pair=pd.DataFrame(prows); erows=[]
    for subset,d in [('all_pairs',pair),('core_secure_pairs',pair[pair.both_core_secure.astype(bool)])]:
        for (proj,stage,thr),g in d.groupby(['projection','stage','A_threshold']):
            row={'subset':subset,'projection':proj,'stage':stage,'A_threshold':thr,'n_biological_pairs':len(g)}
            for m in metrics:
                x=pd.to_numeric(g[m+'_weighted_median'],errors='coerce').dropna().to_numpy(float); row[m+'_median']=float(np.median(x)) if len(x) else np.nan; row[m+'_q25']=float(np.quantile(x,.25)) if len(x) else np.nan; row[m+'_q75']=float(np.quantile(x,.75)) if len(x) else np.nan
            x=pd.to_numeric(g.scenario_support_with_band,errors='coerce').dropna().to_numpy(float); row['scenario_support_with_band_median']=float(np.median(x)) if len(x) else np.nan
            erows.append(row)
    return pair,pd.DataFrame(erows)


def _coupled_line_sensitivity(comb,kdf,jdf,sc,secure):
    rows=[]
    pis=np.sort(kdf.p_intrinsic.unique()); qs=np.sort(kdf.p_component.unique()); dq=float(np.median(np.diff(qs)))
    for proj in ['isi','active']:
        mats={}
        for name,d,col in [('combined',comb,f'A_{proj}'),('kappa',kdf,f'A_{proj}'),('J',jdf,f'A_{proj}')]:
            M=d.pivot(index='p_intrinsic',columns='p_component',values=col).reindex(index=pis,columns=qs).to_numpy(float); mats[name]=np.gradient(M,dq,axis=1)
        for i,pi in enumerate(pis):
            j=int(np.argmin(np.abs(qs-pi)))
            sk=mats['kappa'][i,j]; sj=mats['J'][i,j]; scb=mats['combined'][i,j]; inter=scb-sk-sj
            rows.append({'scenario_id':int(sc.scenario_id),'biological_pair_key':sc.biological_pair_key,'both_core_secure':bool(secure),'within_pair_support_weight':float(sc.within_pair_support_weight),'projection':proj,'p_intrinsic':pi,'dA_dcombined':scb,'dA_dkappa':sk,'dA_dJ':sj,'interaction_gradient_residual':inter})
    return rows


def _aggregate_line(line):
    metrics=['dA_dcombined','dA_dkappa','dA_dJ','interaction_gradient_residual']; prows=[]
    for (pair,proj,pi),g in line.groupby(['biological_pair_key','projection','p_intrinsic']):
        w=pd.to_numeric(g.within_pair_support_weight,errors='coerce').to_numpy(float); row={'biological_pair_key':pair,'both_core_secure':bool(g.both_core_secure.iloc[0]),'projection':proj,'p_intrinsic':pi}
        for m in metrics: row[m+'_weighted_median']=weighted_finite_quantile(pd.to_numeric(g[m],errors='coerce'),w,.5)
        prows.append(row)
    pair=pd.DataFrame(prows); erows=[]
    for subset,d in [('all_pairs',pair),('core_secure_pairs',pair[pair.both_core_secure.astype(bool)])]:
        for (proj,pi),g in d.groupby(['projection','p_intrinsic']):
            row={'subset':subset,'projection':proj,'p_intrinsic':pi,'n_biological_pairs':len(g)}
            for m in metrics:
                x=pd.to_numeric(g[m+'_weighted_median'],errors='coerce').dropna().to_numpy(float); row[m+'_median']=float(np.median(x)) if len(x) else np.nan; row[m+'_q25']=float(np.quantile(x,.25)) if len(x) else np.nan; row[m+'_q75']=float(np.quantile(x,.75)) if len(x) else np.nan
            erows.append(row)
    return pair,pd.DataFrame(erows)


def validate(cfg):
    src=_resolve(cfg); run=src.read_json('RUN_SUMMARY.json'); scenarios=src.read_csv('transition_pair_scenarios_v1_2.csv'); old=src.read_csv('biological_pair_surface_summary.csv.gz',compression='gzip',usecols=['biological_pair_key','both_core_secure']).drop_duplicates(); cps=src.checkpoint_ids()
    mode=cfg['surface'].get('scenario_mode','all_support')
    if mode=='best_only': scenarios=scenarios[(scenarios.wt_source=='best')&(scenarios.sca_source=='best')]
    maxn=cfg['surface'].get('max_scenarios');
    if maxn: scenarios=scenarios.head(int(maxn))
    out={'version':'1.3.0','analysis':'factorized_drive_decomposition','combined_source_version':run.get('version'),'source':src.source_description,'n_selected_scenarios':int(len(scenarios)),'n_all_support_scenarios':int(src.read_csv('transition_pair_scenarios_v1_2.csv').shape[0]),'n_biological_pairs_selected':int(scenarios.biological_pair_key.nunique()),'core_secure_pairs_available':int(old.loc[old.both_core_secure.astype(bool),'biological_pair_key'].nunique()),'combined_checkpoint_count':len(cps),'n_intrinsic':int(cfg['surface']['n_intrinsic']),'n_component':int(cfg['surface']['n_component']),'new_surface_modes':['kappa_only','J_only'],'combined_surface_reused':True,'planned_new_state_rows':int(len(scenarios)*int(cfg['surface']['n_intrinsic'])*int(cfg['surface']['n_component'])*2)}
    src.close()
    if run.get('version')!='1.2.0': raise ValueError('v1.2.0 combined results required')
    if len(cps)!=988: raise ValueError('all 988 combined checkpoints required')
    if mode=='all_support' and not maxn and len(scenarios)!=988: raise ValueError('expected 988 support scenarios')
    return out


def _make_figures(out,bounds,effects,line):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f=out/'figures'; f.mkdir(exist_ok=True)
    b=bounds[(bounds.subset=='core_secure_pairs')&(bounds.projection=='isi')]
    for st in ['WT_exit','balance','SCA3_entry']:
        fig,ax=plt.subplots(figsize=(8.5,5.6))
        for mode in MODES:
            g=b[(b.stage==st)&(b.mode==mode)].sort_values('p_intrinsic')
            ax.plot(g.p_intrinsic,g.median_majority_support,label=mode.replace('_',' '))
        ax.set(xlabel='intrinsic progress',ylabel='required component progress',ylim=(0,1),title=f'ISI {st}: combined vs kappa-only vs J-only'); ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(f/f'01_boundary_{st}.png',dpi=220); plt.close(fig)
    e=effects[effects.subset=='core_secure_pairs']; pis=np.sort(e.p_intrinsic.unique()); qs=np.sort(e.p_component.unique()); M=e.pivot(index='p_component',columns='p_intrinsic',values='interaction_isi_median').reindex(index=qs,columns=pis).to_numpy(float)
    fig,ax=plt.subplots(figsize=(8.5,6)); mesh=ax.pcolormesh(pis,qs,M,shading='auto'); fig.colorbar(mesh,ax=ax,label='ISI interaction contrast'); ax.set(xlabel='intrinsic progress',ylabel='component progress',title='Non-additive kappa x J interaction (core-secure)'); fig.tight_layout(); fig.savefig(f/'02_interaction_surface_ISI.png',dpi=220); plt.close(fig)
    l=line[(line.subset=='core_secure_pairs')&(line.projection=='isi')].sort_values('p_intrinsic'); fig,ax=plt.subplots(figsize=(8.5,5.6))
    ax.plot(l.p_intrinsic,np.abs(l.dA_dkappa_median),label='|dA/d kappa progress|'); ax.plot(l.p_intrinsic,np.abs(l.dA_dJ_median),label='|dA/d J progress|'); ax.plot(l.p_intrinsic,np.abs(l.interaction_gradient_residual_median),label='|interaction residual|'); ax.set(xlabel='intrinsic progress on coupled line',ylabel='local sensitivity magnitude',title='Drive decomposition along coupled WT-SCA3 path'); ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(f/'03_coupled_line_component_sensitivity_ISI.png',dpi=220); plt.close(fig)


def run_all(cfg,resume=True):
    src=_resolve(cfg); run=src.read_json('RUN_SUMMARY.json'); thresholds=_thresholds(run); allsc=src.read_csv('transition_pair_scenarios_v1_2.csv'); old=src.read_csv('biological_pair_surface_summary.csv.gz',compression='gzip',usecols=['biological_pair_key','both_core_secure']).drop_duplicates(); secure_map=dict(zip(old.biological_pair_key,old.both_core_secure.astype(bool)))
    scenarios=allsc.copy(); mode=cfg['surface'].get('scenario_mode','all_support')
    if mode=='best_only': scenarios=scenarios[(scenarios.wt_source=='best')&(scenarios.sca_source=='best')].copy()
    maxn=cfg['surface'].get('max_scenarios');
    if maxn: scenarios=scenarios.head(int(maxn)).copy()
    out=Path(cfg['output']['dir']).resolve(); out.mkdir(parents=True,exist_ok=True); geom=_geom(Path(__file__).resolve().parent.parent)
    workers=int(cfg.get('parallel',{}).get('workers',1)); pending=[]
    for _,sc in scenarios.iterrows():
        if resume and _checkpoint_path(out,'kappa_only',sc.scenario_id).exists() and _checkpoint_path(out,'J_only',sc.scenario_id).exists(): continue
        pending.append(sc)
    if workers>1 and pending:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            fut={ex.submit(_scenario_job,sc,cfg,geom):int(sc.scenario_id) for sc in pending}
            for f in as_completed(fut):
                sid,k,j=f.result(); _write_checkpoint(_checkpoint_path(out,'kappa_only',sid),k); _write_checkpoint(_checkpoint_path(out,'J_only',sid),j)
    else:
        for sc in pending:
            sid,k,j=_scenario_job(sc,cfg,geom); _write_checkpoint(_checkpoint_path(out,'kappa_only',sid),k); _write_checkpoint(_checkpoint_path(out,'J_only',sid),j)
    pis=np.linspace(0,1,int(cfg['surface']['n_intrinsic'])); qs=np.linspace(0,1,int(cfg['surface']['n_component'])); persistence=int(cfg['staging']['persistence_points']); band=float(cfg['decomposition']['boundary_A_band'])
    pair_cross=[]; pair_effect_rows=[]; stage_rows=[]; line_rows=[]
    for pair,gsc in scenarios.groupby('biological_pair_key',sort=True):
        secure=secure_map.get(pair,False); cframes=[]; eframes=[]; weights=[]
        for _,sc in gsc.iterrows():
            kdf=_read_checkpoint(_checkpoint_path(out,'kappa_only',sc.scenario_id)); jdf=_read_checkpoint(_checkpoint_path(out,'J_only',sc.scenario_id)); cd=_sample_combined(src.read_checkpoint(int(sc.scenario_id)),pis,qs)
            for m,d in [('combined',cd),('kappa_only',kdf),('J_only',jdf)]: cframes.append((m,_scenario_crossings(d,sc,secure,thresholds,persistence,m)))
            ef,sr=_scenario_effects(cd,kdf,jdf,sc,secure,thresholds,band); eframes.append(ef); stage_rows.extend(sr); line_rows.extend(_coupled_line_sensitivity(cd,kdf,jdf,sc,secure)); weights.append(float(sc.within_pair_support_weight))
        # aggregate crossings mode by mode preserving scenario ordering
        for m in MODES:
            frames=[f for mm,f in cframes if mm==m]; pair_cross.append(_aggregate_pair_crossings(frames,np.asarray(weights,float)))
        pair_effect_rows.extend(_pair_effect_surface(pair,eframes,np.asarray(weights,float)))
    src.close()
    pair_cross=pd.concat(pair_cross,ignore_index=True); ens_cross=_ensemble_crossings(pair_cross); pair_cross.to_csv(out/'biological_pair_component_boundaries_v1_3.csv',index=False); ens_cross.to_csv(out/'component_boundary_summary_v1_3.csv',index=False)
    pair_effect=pd.DataFrame(pair_effect_rows); ens_effect=_ensemble_effect(pair_effect); pair_effect.to_csv(out/'biological_pair_component_effect_surface_v1_3.csv.gz',index=False,compression='gzip'); ens_effect.to_csv(out/'component_effect_surface_summary_v1_3.csv',index=False)
    pair_stage,ens_stage=_aggregate_stage_effects(pd.DataFrame(stage_rows)); pair_stage.to_csv(out/'biological_pair_interaction_at_stage_boundaries_v1_3.csv',index=False); ens_stage.to_csv(out/'interaction_at_stage_boundaries_v1_3.csv',index=False)
    pair_line,ens_line=_aggregate_line(pd.DataFrame(line_rows)); pair_line.to_csv(out/'biological_pair_coupled_line_component_sensitivity_v1_3.csv',index=False); ens_line.to_csv(out/'coupled_line_component_sensitivity_v1_3.csv',index=False)
    scenarios.to_csv(out/'transition_pair_scenarios_v1_3.csv',index=False)
    _make_figures(out,ens_cross,ens_effect,ens_line)
    # invariants on newly simulated maps at q=pi: both isolated maps must reproduce coupled map.
    diag=[]
    dsrc=_resolve(cfg)
    for sid in scenarios.scenario_id.head(min(20,len(scenarios))):
        sc=scenarios[scenarios.scenario_id==sid].iloc[0]; kdf=_read_checkpoint(_checkpoint_path(out,'kappa_only',sid)); jdf=_read_checkpoint(_checkpoint_path(out,'J_only',sid)); cd=_sample_combined(dsrc.read_checkpoint(int(sid)),pis,qs)
        for pi in pis:
            jj=int(np.argmin(np.abs(qs-pi))); q=qs[jj]
            for proj in ['A_isi','A_active']:
                a=float(cd[(np.isclose(cd.p_intrinsic,pi))&(np.isclose(cd.p_component,q))][proj].iloc[0]); b=float(kdf[(np.isclose(kdf.p_intrinsic,pi))&(np.isclose(kdf.p_component,q))][proj].iloc[0]); c=float(jdf[(np.isclose(jdf.p_intrinsic,pi))&(np.isclose(jdf.p_component,q))][proj].iloc[0]); diag.append(abs(a-b) if np.isfinite(a) and np.isfinite(b) else np.nan); diag.append(abs(a-c) if np.isfinite(a) and np.isfinite(c) else np.nan)
    dsrc.close()
    maxdiag=float(np.nanmax(diag)) if len(diag) else np.nan
    summary={'version':'1.3.0','analysis':'factorized_drive_decomposition','combined_surface_reused_from':'1.2.0','new_surface_modes':['kappa_only','J_only'],'n_selected_scenarios':int(len(scenarios)),'n_biological_pairs':int(scenarios.biological_pair_key.nunique()),'core_secure_pairs':int(sum(1 for p in scenarios.biological_pair_key.unique() if secure_map.get(p,False))),'n_intrinsic':len(pis),'n_component':len(qs),'new_state_rows':int(len(scenarios)*len(pis)*len(qs)*2),'primary_ISI_boundaries':thresholds['isi'],'secondary_active_boundaries':thresholds['active'],'max_coupled_identity_abs_A_difference_diagnostic_first20':maxdiag}
    (out/'RUN_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True),encoding='utf-8')
    core=ens_stage[(ens_stage.subset=='core_secure_pairs')&(ens_stage.projection=='isi')]
    lines=['# Transition ensemble v1.3 — factorized drive decomposition','','Combined drive is decomposed into fitted HR input scaling kappa_I and experimental current protocol J. The frozen v1.2 combined surface is reused; only kappa-only and J-only maps are newly simulated.','','Factorial contrasts at fixed intrinsic progress p_i and matched component progress q:','','- K = A(kappa=q, J=p_i) - A(kappa=p_i, J=p_i)','- J = A(kappa=p_i, J=q) - A(kappa=p_i, J=p_i)','- Combined = A(kappa=q, J=q) - A(kappa=p_i, J=p_i)','- Interaction = Combined - K - J','','The J-only surface is protocol sensitivity, not a disease parameter trajectory. Raw kappa_I is a fitted model coordinate and remains confounded with capacitance; it is not interpreted as an independent biological phenotype.','','## Core-secure ISI interaction near frozen stage boundaries','']
    for _,r in core.iterrows(): lines.append(f"- {r.stage}: |K|={r.median_abs_delta_kappa_median:.4f}, |J|={r.median_abs_delta_J_median:.4f}, |interaction|={r.median_abs_interaction_median:.4f}, signed interaction={r.median_signed_interaction_median:.4f}.")
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return summary
