from __future__ import annotations
import numpy as np


def circle_match_metrics(exp_spikes, model_spikes, tolerance_ms: float = 10.0) -> dict:
    exp=np.sort(np.asarray(exp_spikes,dtype=float)); mod=np.sort(np.asarray(model_spikes,dtype=float))
    tol=float(tolerance_ms); i=j=0; pairs=[]
    while i<len(exp) and j<len(mod):
        delta=mod[j]-exp[i]
        if abs(delta)<=tol:
            pairs.append((i,j,float(delta))); i+=1; j+=1
        elif mod[j] < exp[i]-tol:
            j+=1
        else:
            i+=1
    n=len(pairs); rec=n/max(len(exp),1); prec=n/max(len(mod),1)
    f1=2*rec*prec/(rec+prec) if rec+prec>0 else 0.0
    ae=np.asarray([abs(p[2]) for p in pairs],dtype=float)
    return {
        'circle_matched':int(n),'circle_recall':float(rec),'circle_precision':float(prec),'circle_f1':float(f1),
        'circle_median_abs_timing_error_ms':float(np.median(ae)) if ae.size else np.nan,
        'circle_max_abs_timing_error_ms':float(np.max(ae)) if ae.size else np.nan,
    }


def classify_sweep(vp_loss, exp_spikes, model_spikes, cfg):
    exp=np.sort(np.asarray(exp_spikes,dtype=float)); mod=np.sort(np.asarray(model_spikes,dtype=float))
    q=cfg['sweep_qc']; n_exp=len(exp); n_mod=len(mod)
    count_abs=abs(n_mod-n_exp); count_frac=count_abs/max(n_exp,1)
    cm=circle_match_metrics(exp,mod,float(q['circle_match_tolerance_ms']))
    lead_tol=float(q.get('leading_extra_tolerance_ms',5.0))
    leading=int(np.sum(mod < exp[0]-lead_tol)) if n_exp else n_mod
    trailing=int(np.sum(mod > exp[-1]+lead_tol)) if n_exp else n_mod
    accept=(
        n_mod>0 and count_frac<=float(q['accept_max_count_error_fraction'])
        and cm['circle_f1']>=float(q['accept_min_circle_f1'])
        and leading<=int(q['accept_max_leading_extra_spikes'])
    )
    bad=(
        n_mod==0 or count_frac>float(q['bad_count_error_fraction'])
        or cm['circle_f1']<float(q['bad_min_circle_f1'])
    )
    decision='ACCEPT' if accept else ('BAD' if bad else 'REVIEW')
    reason=[]
    if n_mod==0: reason.append('NO_MODEL_SPIKES')
    if count_frac>float(q['accept_max_count_error_fraction']): reason.append('COUNT')
    if cm['circle_f1']<float(q['accept_min_circle_f1']): reason.append('CIRCLE')
    if leading>int(q['accept_max_leading_extra_spikes']): reason.append('LEADING_EXTRA')
    latency_exp=float(exp[0]) if n_exp else np.nan; latency_mod=float(mod[0]) if n_mod else np.nan
    latency_err=abs(latency_mod-latency_exp) if np.isfinite(latency_exp) and np.isfinite(latency_mod) else np.nan
    return decision, '+'.join(reason) if reason else 'OK', {
        'vp_loss':float(vp_loss),'n_exp_spikes':int(n_exp),'n_model_spikes':int(n_mod),
        'count_error_abs':int(count_abs),'count_error_fraction':float(count_frac),
        'leading_extra_spikes':int(leading),'trailing_extra_spikes':int(trailing),
        'exp_latency_ms':latency_exp,'model_latency_ms':latency_mod,'latency_error_ms':latency_err,
        **cm,
    }


def classify_cell(sweep_rows, cfg):
    if not sweep_rows:
        return 'BAD', {'reason':'NO_SWEEPS'}
    decisions=[r['joint_sweep_decision'] for r in sweep_rows]
    n=len(decisions); n_acc=decisions.count('ACCEPT'); n_bad=decisions.count('BAD'); n_rev=decisions.count('REVIEW')
    frac_acc=n_acc/n; frac_bad=n_bad/n; frac_nonbad=(n_acc+n_rev)/n
    med_f1=float(np.nanmedian([r['circle_f1'] for r in sweep_rows]))
    med_count=float(np.nanmedian([r['count_error_fraction'] for r in sweep_rows]))
    mean_vp=float(np.nanmean([r['vp_loss'] for r in sweep_rows]))
    q=cfg['cell_qc']
    if (
        frac_acc>=float(q['accept_min_fraction_sweeps_accept'])
        and frac_nonbad>=float(q['accept_min_fraction_sweeps_nonbad'])
        and med_f1>=float(q['accept_min_median_circle_f1'])
    ):
        decision='ACCEPT'
    elif frac_bad>=float(q['bad_min_fraction_sweeps_bad']) or med_f1<float(q['bad_max_median_circle_f1']):
        decision='BAD'
    else:
        decision='REVIEW'
    return decision, {
        'n_sweeps':n,'n_sweeps_accept':n_acc,'n_sweeps_review':n_rev,'n_sweeps_bad':n_bad,
        'fraction_sweeps_accept':frac_acc,'fraction_sweeps_nonbad':frac_nonbad,'fraction_sweeps_bad':frac_bad,
        'median_circle_f1':med_f1,'median_count_error_fraction':med_count,'mean_vp_loss':mean_vp,
        'primary_support': 'MULTI_SWEEP' if n>=int(q.get('primary_min_sweeps',2)) else 'SINGLE_SWEEP_ONLY',
    }
