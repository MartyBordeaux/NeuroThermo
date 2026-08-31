from __future__ import annotations
from pathlib import Path
import json, traceback
import numpy as np
import pandas as pd
from .data import build_cells,validate_bundle,find_seed_row,find_baseline_cell_row,find_threshold_bracket
from .optimize import fit_cell
from .diagnostics import classify_sweep,classify_cell
from .plotting import plot_cell,make_audit_pdf
from .params import params_to_z
from .objective import evaluate_cell
from .identifiability import profile_identifiability

PARAMS=('b','r','s','kappa_I')

def _ensure(out):Path(out).mkdir(parents=True,exist_ok=True)
def _jsonable(x):
    if isinstance(x,(np.floating,np.integer)):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,dict):return {k:_jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [_jsonable(v) for v in x]
    return x

def _load_fit_rows(path):
    return pd.read_csv(path) if Path(path).exists() else pd.DataFrame()

def _threshold_to_dict(row):
    return None if row is None else {k:(v.item() if hasattr(v,'item') else v) for k,v in row.to_dict().items()}

def _cell_result_rows(cell,fit,cfg,threshold_bracket):
    ev=fit['best_eval'];swrows=[];abs_shifts=[];large=0
    warn=float(cfg['loss']['latency_alignment']['large_shift_warning_ms'])
    for s,rr in zip(cell['sweeps'],ev['sweeps']):
        aligned=np.asarray(rr['model_spikes'],float);raw=np.asarray(rr['raw_model_spikes'],float);q=classify_sweep(s['exp_spike_times_ms'],aligned,cfg)
        tau=float(rr.get('latency_shift_ms',np.nan));abs_shifts.append(abs(tau) if np.isfinite(tau) else np.nan)
        if np.isfinite(tau) and abs(tau)>=warn:large+=1
        row={'group':cell['group'],'cell_id':cell['cell_id'],'file':s['file'],'abf_sweep':s['abf_sweep'],'current_pA':s['current_pA'],'capacitance_pF':s['capacitance_pF'],'J':s['J'],
             'fit_end_ms':s['fit_end_ms'],'vp_loss':rr['vp_loss'],'raw_vp_loss':rr['raw_vp_loss'],'count_penalty':rr['count_penalty'],'sweep_loss':rr['loss'],'latency_shift_ms':tau,
             'raw_model_latency_ms':float(raw[0]) if len(raw) else np.nan,'aligned_model_latency_ms':float(aligned[0]) if len(aligned) else np.nan,
             'exp_latency_ms':float(s['exp_spike_times_ms'][0]),'raw_n_model_spikes':len(raw),'aligned_last_spike_ms':float(aligned[-1]) if len(aligned) else np.nan,
             'exp_last_spike_ms':float(s['exp_spike_times_ms'][-1]),'aligned_last_spike_error_ms':rr['aligned_last_spike_error_ms'],
             'raw_model_spike_times_ms':json.dumps(raw.tolist()),'model_spike_times_ms':json.dumps(aligned.tolist()),'exp_spike_times_ms':json.dumps(np.asarray(s['exp_spike_times_ms']).tolist())}
        row.update(q);swrows.append(row)
    cq=classify_cell(swrows,cfg);p=ev['params'];med_abs=float(np.nanmedian(abs_shifts)) if abs_shifts else np.nan
    trow={'group':cell['group'],'cell_id':cell['cell_id'],'capacitance_pF':cell['capacitance_pF'],'spike_train_loss':ev['spike_train_loss'],'threshold_loss':ev['threshold_loss'],'cell_loss':ev['loss'],
          'fit_source':fit['source'],'elapsed_s':fit['elapsed_s'],'n_spiking_sweeps':len(cell['sweeps']),'median_abs_latency_shift_ms':med_abs,'max_abs_latency_shift_ms':float(np.nanmax(abs_shifts)) if abs_shifts else np.nan,
          'n_large_abs_latency_shifts':large,'large_abs_latency_shift_fraction':large/max(len(abs_shifts),1),'threshold_pass':ev['threshold_pass'],'model_spikes_at_nonspiking_current':ev['model_spikes_at_nonspiking_current'],
          'model_spikes_at_first_spiking_current':ev['model_spikes_at_first_spiking_current'],'b':p['b'],'r':p['r'],'s':p['s'],'kappa_I':p['kappa_I'],'kappa_over_Cm':p['kappa_I']/cell['capacitance_pF']}
    trow.update(cq);return trow,swrows

def validate(cfg):return validate_bundle(cfg)

def run(cfg):
    out=Path(cfg['output']['dir']);_ensure(out);_ensure(out/'plots');_ensure(out/'fits');cells,bundle=build_cells(cfg);summary=[];swall=[];trows=[]
    existing=_load_fit_rows(out/'cell_fit_summary.csv');exmap={r.cell_id:r for _,r in existing.iterrows()} if not existing.empty else {}
    cells_by_id={c['cell_id']:c for c in cells};sw_by={};th_by={}
    for cell in cells:
        cid=cell['cell_id'];fitjson=out/'fits'/f'{cid}.json'
        try:
            thr=find_threshold_bracket(bundle,cid);thd=_threshold_to_dict(thr)
            if cfg['output']['resume'] and fitjson.exists():
                saved=json.loads(fitjson.read_text());z=np.asarray(saved['best_z'],float);ev=evaluate_cell(cell,z,cfg,dt_ms=cfg['optimization']['dt_fine_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=thd);fit={'best_z':z,'best_eval':ev,'source':saved.get('source','resume'),'elapsed_s':saved.get('elapsed_s',0.0)}
            else:
                seed=find_seed_row(bundle,cid)
                if seed is None:raise RuntimeError('No v3.8 seed row')
                fit=fit_cell(cell,seed,cfg,threshold_bracket=thd);fitjson.write_text(json.dumps(_jsonable({'best_z':fit['best_z'],'params':fit['best_eval']['params'],'loss':fit['best_eval']['loss'],'source':fit['source'],'elapsed_s':fit['elapsed_s']}),indent=2))
            crow,srows=_cell_result_rows(cell,fit,cfg,thd);base=find_baseline_cell_row(bundle,cid)
            if base is not None:
                for p in PARAMS:
                    crow['v3_8_'+p]=float(base[p])
                    crow['delta_'+p]=crow[p]-float(base[p])
                crow['v3_8_loss']=float(base['cell_loss']) if 'cell_loss' in base.index else np.nan
            summary.append(crow);swall.extend(srows);sw_by[cid]=srows;trows.append({'cell_id':cid,'nonspiking_current_pA':float(thr['nonspiking_current_pA']),'first_spiking_current_pA':float(thr['first_spiking_current_pA']),'model_spikes_at_nonspiking_current':crow['model_spikes_at_nonspiking_current'],'model_spikes_at_first_spiking_current':crow['model_spikes_at_first_spiking_current'],'threshold_pass':crow['threshold_pass']})
            plot_cell(cell,srows,trows[-1],crow,out/'plots'/f'{cid}.png')
        except Exception as e:
            summary.append({'group':cell['group'],'cell_id':cid,'auto_cell_decision':'BAD','error':repr(e),'traceback':traceback.format_exc()})
    sdf=pd.DataFrame(summary);swd=pd.DataFrame(swall);tdf=pd.DataFrame(trows)
    sdf['final_v3_9_decision']=sdf.get('post_v3_9_auto_cell_decision',sdf.get('auto_cell_decision','')).fillna(sdf.get('auto_cell_decision',''))
    sdf['review_comment']=''
    sdf.to_csv(out/'cell_fit_summary.csv',index=False);swd.to_csv(out/'sweep_fit_summary.csv',index=False);tdf.to_csv(out/'threshold_constraint_summary.csv',index=False)
    stress=[]
    for _,r in sdf.iterrows():
        if 's' not in r or not np.isfinite(r.get('s',np.nan)):continue
        stress.append({'cell_id':r['cell_id'],'group':r['group'],'v3_8_s':r.get('v3_8_s',np.nan),'v3_9_s':r['s'],'delta_s':r.get('delta_s',np.nan),'moved_below_old_s_min_0p25':bool(r['s']<0.25-1e-9),'distance_from_new_s_min':float(r['s']-0.05),'fraction_of_new_s_range_above_min':float((r['s']-0.05)/(15.0-0.05))})
    pd.DataFrame(stress).to_csv(out/'s_boundary_stress_summary.csv',index=False)
    make_audit_pdf(cells_by_id,sw_by,{r['cell_id']:r for r in trows},summary,out/'joint_fit_visual_audit_v3_9.pdf')
    return {'n_cells':len(sdf),'output':str(out)}

def identify_final(cfg):
    out=Path(cfg['output']['dir']);p=out/'cell_fit_summary.csv';df=pd.read_csv(p);cells,bundle=build_cells(cfg);idrows=[]
    for cell in cells:
        cid=cell['cell_id'];row=df[df.cell_id==cid]
        if len(row)!=1 or row.iloc[0]['final_v3_9_decision']!='ACCEPT' or int(row.iloc[0].get('n_spiking_sweeps',0))<cfg['identifiability']['min_spiking_sweeps']:continue
        params={k:float(row.iloc[0][k]) for k in PARAMS};z=params_to_z(params,cfg['bounds']);thr=find_threshold_bracket(bundle,cid);thd=_threshold_to_dict(thr)
        ev=evaluate_cell(cell,z,cfg,dt_ms=cfg['identifiability']['dt_ms'],search_tau_ms=cfg['loss']['vp_tau_ms'],threshold_bracket=thd,identifiability_mode=True);res=profile_identifiability(cell,z,ev,cfg,threshold_bracket=thd)
        idrows.extend(res.get('alternative_rows',[]))
    id_df=pd.DataFrame(idrows);id_df.to_csv(out/'final_identifiability.csv',index=False)
    if not id_df.empty:
        up={}
        for cid,g in id_df.groupby('cell_id'):
            per={r.parameter_tested:r.parameter_status for _,r in g.iterrows()};up[cid]={'identifiability':'IDENTIFIABLE' if all(per.get(p)=='IDENTIFIABLE' for p in PARAMS) else 'NON_IDENTIFIABLE'}
            for par in PARAMS:up[cid]['id_'+par]=per.get(par,'NOT_TESTED')
        for i,r in df.iterrows():
            cid=r.cell_id
            if cid in up:
                for k,v in up[cid].items():df.loc[i,k]=v
            elif r['final_v3_9_decision']=='ACCEPT' and r['primary_support']=='SINGLE_SWEEP_ONLY':df.loc[i,'identifiability']='INSUFFICIENT_SPIKING_SWEEPS'
            elif r['final_v3_9_decision']!='ACCEPT':df.loc[i,'identifiability']='NOT_RUN'
        df.to_csv(p,index=False)
    return {'n_identified_cells':int(id_df.cell_id.nunique()) if not id_df.empty else 0}
