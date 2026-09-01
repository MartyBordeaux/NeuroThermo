from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .data import load_inputs, theta_from_row, selected_spikes, near_optimal_sets
from .model import simulate_trace, has_spike
from .metrics import spike_metrics, symmetric_relative_difference, align_first_spike, restrict_spikes
from .phase import supported_phase_profile, phase_descriptors, profile_nrmse, PHASE_KEYS

SCALAR_METRICS=['firing_rate_hz','active_rate_hz','spike_count','train_duration_ms','mean_isi_ms','median_isi_ms','cv_isi',
                'early_isi_ms','late_isi_ms','adaptation_index','cycle_period_ms','x_range','y_range','z_range',
                'mean_speed','peak_speed','mean_divergence','min_divergence','max_divergence',
                'fraction_positive_divergence','cycle_arc_length','z_peak_phase','speed_peak_phase']


def refine_rheobase(theta,tr,cfg):
    lo=float(tr['nonspiking_J']); hi=float(tr['first_spiking_J']); cm=float(tr['capacitance_pF'])
    dt=float(cfg['rheobase']['dt_ms'])
    lo_sp=has_spike(theta,lo,cfg,dt_ms=dt); hi_sp=has_spike(theta,hi,cfg,dt_ms=dt)
    if lo_sp or not hi_sp:
        return {'rheobase_J':np.nan,'rheobase_pA':np.nan,'lo_J':lo,'hi_J':hi,'iterations':0,'status':'BRACKET_VIOLATION'}
    tol_pa=float(cfg['rheobase']['tolerance_pA']); it=0
    while (hi-lo)*cm>tol_pa and it<int(cfg['rheobase']['max_iterations']):
        mid=.5*(lo+hi)
        if has_spike(theta,mid,cfg,dt_ms=dt): hi=mid
        else: lo=mid
        it+=1
    return {'rheobase_J':hi,'rheobase_pA':hi*cm,'lo_J':lo,'hi_J':hi,'iterations':it,'status':'BRACKET_OK'}


def simulate_at_observed_current(theta,sweep,exp_spikes,cfg):
    obs=max(float(cfg['analysis']['observation_end_ms']),float(sweep['fit_end_ms']))
    stim=float(cfg['analysis']['stimulus_duration_ms'])
    t,x,y,z,raw,ok=simulate_trace(theta,float(sweep['J']),cfg,duration_ms=stim,observation_end_ms=obs,dt_ms=float(cfg['analysis']['dt_ms']))
    aligned,tau=align_first_spike(exp_spikes,raw)
    support_end=float(sweep['fit_end_ms'])
    supported=restrict_spikes(aligned,0.0,support_end)
    expmet=spike_metrics(exp_spikes,support_end)
    modmet=spike_metrics(supported,support_end)
    prof=supported_phase_profile(t,x,y,z,raw,aligned,exp_spikes,float(sweep['J']),theta,cfg) if ok and np.isfinite(tau) else None
    pdesc=phase_descriptors(prof)
    out={**modmet,**pdesc,'simulation_ok':bool(ok),'latency_shift_ms':tau,'support_end_ms':support_end,
         'raw_model_spike_count_full_observation':int(len(raw)),'aligned_model_spike_count_supported':int(len(supported)),
         'model_spikes_excluded_after_support':int(max(0,len(raw)-len(supported)))}
    expout={**expmet,'support_end_ms':support_end}
    return out,expout,prof


def interp_scalar(x,y,target):
    x=np.asarray(x,float); y=np.asarray(y,float)
    m=np.isfinite(x)&np.isfinite(y)
    x=x[m]; y=y[m]
    if len(x)==0: return np.nan
    order=np.argsort(x); x=x[order]; y=y[order]
    # combine duplicate x by median
    ux=np.unique(x); uy=np.array([np.median(y[x==v]) for v in ux],float)
    if target<ux[0]-1e-12 or target>ux[-1]+1e-12: return np.nan
    exact=np.where(np.isclose(ux,target,rtol=0,atol=1e-10))[0]
    if len(exact): return float(uy[exact[0]])
    if len(ux)<2: return np.nan
    return float(np.interp(target,ux,uy))


def interp_profile(records,target,n_bins):
    good=[r for r in records if r['profile'] is not None and np.isfinite(r['q'])]
    if not good: return None
    qs=np.array([r['q'] for r in good],float)
    if target<np.min(qs)-1e-12 or target>np.max(qs)+1e-12: return None
    result={'phase':np.linspace(0.0,1.0,n_bins),'n_cycles':0,'period_ms':np.nan}
    for key in PHASE_KEYS:
        mat=np.stack([r['profile'][key] for r in good])
        vals=[]
        for j in range(mat.shape[1]): vals.append(interp_scalar(qs,mat[:,j],target))
        arr=np.asarray(vals,float)
        if not np.all(np.isfinite(arr)): return None
        result[key]=arr
    periods=np.array([r['profile']['period_ms'] for r in good],float)
    result['period_ms']=interp_scalar(qs,periods,target)
    result['n_cycles']=int(np.median([r['profile']['n_cycles'] for r in good]))
    return result


def build_solution_rows(primary,manifest,events,thresholds,alts,cfg):
    trmap={str(r.cell_id):r for _,r in thresholds.iterrows()}
    man=manifest[manifest.cell_id.astype(str).isin(set(primary.cell_id.astype(str)))].copy()
    actual_rows=[]; best_phase_rows=[]; actual_rob=[]; phase_rob=[]; rheo_rows=[]
    q_rows=[]; q_exp_rows=[]; q_rob=[]; q_phase_best=[]; q_phase_rob=[]; coverage=[]
    q_targets=[float(q) for q in cfg['analysis']['q_targets']]
    n_bins=int(cfg['phase']['n_bins'])

    for _,cell in primary.iterrows():
        cid=str(cell.cell_id); group=str(cell.group); tr=trmap[cid]; best_theta=theta_from_row(cell)
        rb=refine_rheobase(best_theta,tr,cfg)
        rheo_rows.append({'group':group,'cell_id':cid,'solution':'best','source':'best',**rb})
        cellman=man[man.cell_id.astype(str).eq(cid)].sort_values('J').copy()
        if not np.isfinite(rb['rheobase_J']) or cellman.empty: continue
        jmax=float(cellman.J.max()); den=jmax-float(rb['rheobase_J'])
        if den<=0: continue
        cellman['q_ref']=(cellman.J.astype(float)-float(rb['rheobase_J']))/den
        qmin=float(cellman.q_ref.min()); qmax=float(cellman.q_ref.max())
        cov={'group':group,'cell_id':cid,'best_rheobase_J':rb['rheobase_J'],'J_min_observed':float(cellman.J.min()),'J_max_observed':jmax,
             'q_min_observed':qmin,'q_max_observed':qmax,'n_observed_spiking_currents':int(len(cellman))}
        for q in q_targets: cov['q_%g_supported'%q]=bool(q>=qmin-1e-12 and q<=qmax+1e-12)
        coverage.append(cov)

        solutions=[{'solution':'best','source':'best','theta':best_theta,'loss':float(cell.get('cell_loss',np.nan))}]
        for alt in near_optimal_sets(alts,cid):
            arb=refine_rheobase(alt['theta'],tr,cfg)
            rheo_rows.append({'group':group,'cell_id':cid,'solution':'alternative','source':alt['source'],'alt_loss':alt['loss'],**arb})
            solutions.append({'solution':'alternative','source':alt['source'],'theta':alt['theta'],'loss':alt['loss']})

        by_solution={}
        phase_by_solution={}
        exp_actual=[]
        for sol in solutions:
            srows=[]; sprofs=[]
            for _,sw in cellman.iterrows():
                exp=selected_spikes(events,sw.sweep_id)
                met,expmet,prof=simulate_at_observed_current(sol['theta'],sw,exp,cfg)
                base={'group':group,'cell_id':cid,'solution':sol['solution'],'source':sol['source'],'solution_loss':sol['loss'],
                      'sweep_id':sw.sweep_id,'current_pA':float(sw.current_pA),'J':float(sw.J),'q_ref':float(sw.q_ref),
                      'fit_end_ms':float(sw.fit_end_ms),'n_exp_spikes':int(len(exp))}
                row={**base,**met}
                for k,v in expmet.items(): row['exp_'+k]=v
                actual_rows.append(row); srows.append(row); sprofs.append({'q':float(sw.q_ref),'sweep_id':sw.sweep_id,'profile':prof})
                if sol['solution']=='best':
                    exp_actual.append({**base,**{'exp_'+k:v for k,v in expmet.items()}})
                    if prof is not None:
                        for i,p in enumerate(prof['phase']):
                            best_phase_rows.append({**base,'phase':float(p),**{k:float(prof[k][i]) for k in PHASE_KEYS}})
            by_solution[sol['source']]=pd.DataFrame(srows); phase_by_solution[sol['source']]=sprofs

        bestdf=by_solution['best']; bestprofmap={r['sweep_id']:r['profile'] for r in phase_by_solution['best']}
        # Same-current robustness: alternatives are compared at exactly the same observed J.
        for _,br in bestdf.iterrows():
            sid=br.sweep_id
            for metric in SCALAR_METRICS:
                vals=[]
                for source,adf in by_solution.items():
                    if source=='best': continue
                    ar=adf[adf.sweep_id.eq(sid)]
                    if ar.empty: continue
                    d=symmetric_relative_difference(float(br.get(metric,np.nan)),float(ar.iloc[0].get(metric,np.nan)))
                    if np.isfinite(d): vals.append(d)
                actual_rob.append({'group':group,'cell_id':cid,'sweep_id':sid,'current_pA':br.current_pA,'J':br.J,'q_ref':br.q_ref,'metric':metric,
                                   'n_near_opt_alternatives':len(vals),'median_srd':float(np.median(vals)) if vals else np.nan,
                                   'max_srd':float(np.max(vals)) if vals else np.nan,
                                   'stable_20pct':bool(vals and np.max(vals)<=float(cfg['robustness']['scalar_srd_threshold'])),
                                   'status':'EVALUATED' if vals else 'NO_NEAR_OPT_ALTERNATIVE_FOUND'})
            bp=bestprofmap.get(sid)
            for key in PHASE_KEYS:
                vals=[]
                for source,recs in phase_by_solution.items():
                    if source=='best': continue
                    ap=next((r['profile'] for r in recs if r['sweep_id']==sid),None)
                    d=profile_nrmse(bp,ap,key)
                    if np.isfinite(d): vals.append(d)
                phase_rob.append({'group':group,'cell_id':cid,'sweep_id':sid,'current_pA':br.current_pA,'J':br.J,'q_ref':br.q_ref,'profile':key,
                                  'n_near_opt_alternatives':len(vals),'median_nrmse':float(np.median(vals)) if vals else np.nan,
                                  'max_nrmse':float(np.max(vals)) if vals else np.nan,
                                  'stable_20pct':bool(vals and np.max(vals)<=float(cfg['robustness']['phase_nrmse_threshold'])),
                                  'status':'EVALUATED' if vals else 'NO_NEAR_OPT_ALTERNATIVE_FOUND'})

        # q-space interpolation, never extrapolating beyond the actually observed spiking-current support.
        expdf=pd.DataFrame(exp_actual)
        for q in q_targets:
            exprow={'group':group,'cell_id':cid,'q':q,'supported_by_observed_current_range':bool(q>=qmin-1e-12 and q<=qmax+1e-12)}
            for metric in ['spike_count','firing_rate_hz','active_rate_hz','train_duration_ms','mean_isi_ms','median_isi_ms','cv_isi','early_isi_ms','late_isi_ms','adaptation_index']:
                exprow[metric]=interp_scalar(expdf.q_ref,expdf['exp_'+metric],q)
            q_exp_rows.append(exprow)
            qsol={}
            qprof={}
            for source,sdf in by_solution.items():
                row={'group':group,'cell_id':cid,'solution':'best' if source=='best' else 'alternative','source':source,'q':q,
                     'supported_by_observed_current_range':bool(q>=qmin-1e-12 and q<=qmax+1e-12)}
                for metric in SCALAR_METRICS:
                    row[metric]=interp_scalar(sdf.q_ref,sdf[metric],q)
                q_rows.append(row); qsol[source]=row
                prof=interp_profile(phase_by_solution[source],q,n_bins); qprof[source]=prof
                if source=='best' and prof is not None:
                    for i,p in enumerate(prof['phase']):
                        q_phase_best.append({'group':group,'cell_id':cid,'q':q,'phase':float(p),**{k:float(prof[k][i]) for k in PHASE_KEYS}})
            bestq=qsol['best']
            for metric in SCALAR_METRICS:
                vals=[]
                for source,row in qsol.items():
                    if source=='best': continue
                    d=symmetric_relative_difference(float(bestq.get(metric,np.nan)),float(row.get(metric,np.nan)))
                    if np.isfinite(d): vals.append(d)
                q_rob.append({'group':group,'cell_id':cid,'q':q,'metric':metric,'n_near_opt_alternatives':len(vals),
                              'median_srd':float(np.median(vals)) if vals else np.nan,'max_srd':float(np.max(vals)) if vals else np.nan,
                              'stable_20pct':bool(vals and np.max(vals)<=float(cfg['robustness']['scalar_srd_threshold'])),
                              'status':'EVALUATED' if vals else 'NO_SUPPORTED_NEAR_OPT_ALTERNATIVE'})
            for key in PHASE_KEYS:
                vals=[]
                for source,prof in qprof.items():
                    if source=='best': continue
                    d=profile_nrmse(qprof['best'],prof,key)
                    if np.isfinite(d): vals.append(d)
                q_phase_rob.append({'group':group,'cell_id':cid,'q':q,'profile':key,'n_near_opt_alternatives':len(vals),
                                    'median_nrmse':float(np.median(vals)) if vals else np.nan,'max_nrmse':float(np.max(vals)) if vals else np.nan,
                                    'stable_20pct':bool(vals and np.max(vals)<=float(cfg['robustness']['phase_nrmse_threshold'])),
                                    'status':'EVALUATED' if vals else 'NO_SUPPORTED_NEAR_OPT_ALTERNATIVE'})

    return tuple(pd.DataFrame(x) for x in [actual_rows,best_phase_rows,actual_rob,phase_rob,rheo_rows,q_rows,q_exp_rows,q_rob,q_phase_best,q_phase_rob,coverage])


def _safe_median(series):
    a=pd.to_numeric(series,errors='coerce').to_numpy(float)
    a=a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else np.nan


def summarize_q(qmodel,qexp,qphase,primary):
    idmap=primary[['cell_id','animal_id','animal_resolved']].copy()
    best=qmodel[qmodel.solution.eq('best')].merge(idmap,on='cell_id',how='left')
    exp=qexp.merge(idmap,on='cell_id',how='left')
    core=['firing_rate_hz','active_rate_hz','spike_count','train_duration_ms','mean_isi_ms','median_isi_ms','adaptation_index']
    model_extra=['cycle_period_ms','x_range','y_range','z_range','mean_speed','mean_divergence','fraction_positive_divergence','cycle_arc_length']
    group=[]; animal=[]
    for source,df,metrics in [('experiment',exp,core),('model_best',best,core+model_extra)]:
        for (g,q),x in df.groupby(['group','q']):
            row={'source':source,'group':g,'q':q,'n_cells_current_support':int(x.supported_by_observed_current_range.fillna(False).sum())}
            for k in metrics:
                row[k+'_median']=_safe_median(x[k]) if k in x else np.nan
                row[k+'_n']=int(x[k].notna().sum()) if k in x else 0
            group.append(row)
        res=df[df.animal_resolved.fillna(False).astype(bool)]
        for (g,a,q),x in res.groupby(['group','animal_id','q']):
            row={'source':source,'group':g,'animal_id':a,'q':q,'n_cells':len(x)}
            for k in metrics: row[k+'_median']=_safe_median(x[k]) if k in x else np.nan
            animal.append(row)
    phase_group=[]
    if not qphase.empty:
        for (g,q,ph),x in qphase.groupby(['group','q','phase']):
            phase_group.append({'group':g,'q':q,'phase':ph,'n_cells':x.cell_id.nunique(),**{k:float(x[k].median()) for k in PHASE_KEYS}})
    return pd.DataFrame(group),pd.DataFrame(animal),pd.DataFrame(phase_group)


def run_all(cfg):
    out=Path(cfg['output']['dir']); out.mkdir(parents=True,exist_ok=True)
    primary,manifest,events,thresholds,alts,ids=load_inputs(cfg)
    (actual,phase_actual,rob_actual,prob_actual,rheo,qmodel,qexp,qrob,qphase,qprob,coverage)=build_solution_rows(primary,manifest,events,thresholds,alts,cfg)
    group,animal,phase_group=summarize_q(qmodel,qexp,qphase,primary)
    files={'observed_current_dynamics_all_solutions.csv':actual,'observed_current_phase_profiles_best.csv':phase_actual,
           'observed_current_scalar_robustness.csv':rob_actual,'observed_current_phase_robustness.csv':prob_actual,
           'rheobase_refinement_all_solutions.csv':rheo,'q_interpolated_model_all_solutions.csv':qmodel,
           'q_interpolated_experiment.csv':qexp,'q_scalar_robustness_near_optimal.csv':qrob,
           'q_phase_profiles_best.csv':qphase,'q_phase_robustness_near_optimal.csv':qprob,
           'q_support_by_cell.csv':coverage,'group_q_medians.csv':group,'animal_q_medians.csv':animal,
           'group_q_phase_median_profiles.csv':phase_group}
    for name,df in files.items(): df.to_csv(out/name,index=False)
    best_rheo=rheo[rheo.solution.eq('best')]
    qtargets=[float(q) for q in cfg['analysis']['q_targets']]
    cov_summary={str(q):int(coverage['q_%g_supported'%q].sum()) for q in qtargets} if not coverage.empty else {}
    evalrob=rob_actual[rob_actual.status.eq('EVALUATED')]
    qeval=qrob[qrob.status.eq('EVALUATED')]
    summary={'version':'2.1.0','analysis_mode':'experimental-support-restricted','primary_cells':int(len(primary)),
             'WT_cells':int((primary.group=='WT').sum()),'SCA3_cells':int((primary.group=='SCA3').sum()),
             'observed_primary_spiking_sweeps':int(actual[actual.solution.eq('best')].sweep_id.nunique()),
             'selected_experimental_spikes_primary':int(events.cell_id.astype(str).isin(set(primary.cell_id.astype(str))).sum()),
             'best_rheobase_bracket_pass':int((best_rheo.status=='BRACKET_OK').sum()),'best_rheobase_total':int(len(best_rheo)),
             'q_targets':qtargets,'q_target_cell_support_counts':cov_summary,
             'same_current_scalar_evaluations':int(len(evalrob)),'same_current_scalar_stable_20pct_fraction':float(evalrob.stable_20pct.mean()) if len(evalrob) else np.nan,
             'q_scalar_evaluations':int(len(qeval)),'q_scalar_stable_20pct_fraction':float(qeval.stable_20pct.mean()) if len(qeval) else np.nan,
             'formal_animal_level_p_values':False,'extrapolation_allowed':False,'last_spike_anchor':False,'time_rescaling':False}
    (out/'RUN_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
