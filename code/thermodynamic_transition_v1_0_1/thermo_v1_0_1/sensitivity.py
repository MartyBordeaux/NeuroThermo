from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def _safe_rel(a,b,eps=1e-12):return np.abs(a-b)/np.maximum(np.maximum(np.abs(a),np.abs(b)),eps)
def compare_noise(reference,half,double,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);names={'reference':Path(reference),'half':Path(half),'double':Path(double)};frames={k:pd.read_csv(v/'biological_pair_thermodynamic_markers.csv') for k,v in names.items()};key='biological_pair_key';cols=['fisher_interior_peak_p_median','entropy_peak_p_median','kl_balance_p_median','epr_onset_p_median','epr_peak_p_median'];m=frames['reference'][[key]+cols].copy().rename(columns={c:c+'__reference' for c in cols})
    for tag in ['half','double']:m=m.merge(frames[tag][[key]+cols].rename(columns={c:c+f'__{tag}' for c in cols}),on=key,how='inner')
    for c in cols:m[c+'__shift_half']=m[c+'__half']-m[c+'__reference'];m[c+'__shift_double']=m[c+'__double']-m[c+'__reference']
    m.to_csv(out/'noise_sensitivity_by_pair.csv',index=False);rows=[]
    for c in cols:
        for tag in ['half','double']:
            d=m[c+f'__shift_{tag}'].dropna().to_numpy(float);rows.append({'marker':c,'noise_condition':tag,'n_pairs':len(d),'median_shift_p':np.median(d) if len(d) else np.nan,'median_abs_shift_p':np.median(np.abs(d)) if len(d) else np.nan,'q90_abs_shift_p':np.quantile(np.abs(d),.9) if len(d) else np.nan})
    pd.DataFrame(rows).to_csv(out/'noise_sensitivity_marker_summary.csv',index=False);return out
def compare_dt(reference,halfdt,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);r=Path(reference);h=Path(halfdt);mr=pd.read_csv(r/'scenario_thermodynamic_markers.csv');mh=pd.read_csv(h/'scenario_thermodynamic_markers.csv');cols=['fisher_interior_peak_p','entropy_peak_p','kl_balance_p','epr_onset_p','epr_peak_p'];m=mr[['scenario_id','biological_pair_key']+cols].merge(mh[['scenario_id']+cols],on='scenario_id',suffixes=('__dt005','__dt0025'))
    for c in cols:m[c+'__delta']=m[c+'__dt0025']-m[c+'__dt005']
    m.to_csv(out/'dt_convergence_markers_by_scenario.csv',index=False);rows=[]
    for c in cols:
        d=m[c+'__delta'].dropna().to_numpy(float);rows.append({'marker':c,'n_scenarios':len(d),'median_shift_p':np.median(d) if len(d) else np.nan,'median_abs_shift_p':np.median(np.abs(d)) if len(d) else np.nan,'max_abs_shift_p':np.max(np.abs(d)) if len(d) else np.nan})
    pd.DataFrame(rows).to_csv(out/'dt_convergence_marker_summary.csv',index=False);sr=pd.read_csv(r/'stationary_metrics_scenarios.csv');sh=pd.read_csv(h/'stationary_metrics_scenarios.csv');keep=set(m.scenario_id.astype(int));q=sr[sr.scenario_id.astype(int).isin(keep)].merge(sh[sh.scenario_id.astype(int).isin(keep)],on=['scenario_id','p'],suffixes=('__dt005','__dt0025'));curve_rows=[]
    for sid,g in q.groupby('scenario_id'):
        rec={'scenario_id':sid}
        for c in ['kl_balance','fisher','epr','entropy_delta_wt']:
            a=g[c+'__dt005'].to_numpy(float);b=g[c+'__dt0025'].to_numpy(float);rec[c+'_median_relative_difference']=float(np.nanmedian(_safe_rel(a,b)));rec[c+'_max_relative_difference']=float(np.nanmax(_safe_rel(a,b)))
        curve_rows.append(rec)
    pd.DataFrame(curve_rows).to_csv(out/'dt_convergence_curve_diagnostics.csv',index=False);return out
def main():
    ap=argparse.ArgumentParser(prog='thermo_v1_0_1.sensitivity');sub=ap.add_subparsers(dest='cmd',required=True);p=sub.add_parser('noise');p.add_argument('--reference',required=True);p.add_argument('--half',required=True);p.add_argument('--double',required=True);p.add_argument('--out',required=True);p=sub.add_parser('dt');p.add_argument('--reference',required=True);p.add_argument('--halfdt',required=True);p.add_argument('--out',required=True);a=ap.parse_args();out=compare_noise(a.reference,a.half,a.double,a.out) if a.cmd=='noise' else compare_dt(a.reference,a.halfdt,a.out);print(out)
if __name__=='__main__':main()
