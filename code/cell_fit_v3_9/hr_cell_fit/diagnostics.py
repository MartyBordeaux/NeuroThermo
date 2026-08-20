from __future__ import annotations
import numpy as np
from scipy.optimize import linear_sum_assignment


def one_to_one_circle_match(exp_times, model_times, tol_ms=10.0):
    e=np.asarray(exp_times,float); m=np.asarray(model_times,float)
    if e.size==0 and m.size==0:return {'matched':0,'precision':1.0,'recall':1.0,'f1':1.0,'median_abs_error_ms':0.0}
    if e.size==0:return {'matched':0,'precision':0.0,'recall':1.0,'f1':0.0,'median_abs_error_ms':np.nan}
    if m.size==0:return {'matched':0,'precision':1.0,'recall':0.0,'f1':0.0,'median_abs_error_ms':np.nan}
    C=np.abs(e[:,None]-m[None,:]); big=1e6;C2=np.where(C<=tol_ms,C,big)
    r,c=linear_sum_assignment(C2);good=C2[r,c]<big;rr=r[good];cc=c[good];matched=len(rr)
    precision=matched/len(m);recall=matched/len(e);f1=2*precision*recall/(precision+recall) if precision+recall else 0.0
    mae=float(np.median(np.abs(e[rr]-m[cc]))) if matched else np.nan
    return {'matched':matched,'precision':precision,'recall':recall,'f1':f1,'median_abs_error_ms':mae}


def leading_extra_spikes(exp_times,model_times,tol_ms=5.0):
    e=np.asarray(exp_times,float);m=np.asarray(model_times,float)
    if m.size==0:return 0
    if e.size==0:return int(m.size)
    return int(np.sum(m<e[0]-tol_ms))


def classify_sweep(exp_times,model_times,cfg):
    qc=cfg['sweep_qc']; circ=one_to_one_circle_match(exp_times,model_times,qc['circle_match_tolerance_ms'])
    ne=len(exp_times);nm=len(model_times);count_err=abs(nm-ne)/max(ne,1);lead=leading_extra_spikes(exp_times,model_times,qc['leading_extra_tolerance_ms'])
    if count_err<=qc['accept_max_count_error_fraction'] and circ['f1']>=qc['accept_min_circle_f1'] and lead<=qc['accept_max_leading_extra_spikes']:
        decision='ACCEPT'
    elif count_err>=qc['bad_count_error_fraction'] or circ['f1']<qc['bad_min_circle_f1']:
        decision='BAD'
    else:decision='UNCERTAIN'
    return {'n_exp_spikes':ne,'n_model_spikes':nm,'count_error_fraction':count_err,'leading_extra_spikes':lead,'circle_f1':circ['f1'],'circle_precision':circ['precision'],'circle_recall':circ['recall'],'circle_median_abs_error_ms':circ['median_abs_error_ms'],'joint_sweep_decision':decision}


def classify_cell(sweep_rows,cfg):
    q=cfg['cell_qc'];ds=[r['joint_sweep_decision'] for r in sweep_rows];n=max(len(ds),1)
    frac_acc=sum(x=='ACCEPT' for x in ds)/n;frac_bad=sum(x=='BAD' for x in ds)/n;frac_nonbad=1-frac_bad
    med_f1=float(np.median([r['circle_f1'] for r in sweep_rows])) if sweep_rows else 0.0
    if frac_bad>=q['bad_min_fraction_sweeps_bad'] or med_f1<=q['bad_max_median_circle_f1']:decision='BAD'
    elif frac_acc>=q['accept_min_fraction_sweeps_accept'] and frac_nonbad>=q['accept_min_fraction_sweeps_nonbad'] and med_f1>=q['accept_min_median_circle_f1']:decision='ACCEPT'
    else:decision='UNCERTAIN'
    return {'fraction_sweeps_accept':frac_acc,'fraction_sweeps_bad':frac_bad,'median_circle_f1':med_f1,'auto_cell_decision':decision}
