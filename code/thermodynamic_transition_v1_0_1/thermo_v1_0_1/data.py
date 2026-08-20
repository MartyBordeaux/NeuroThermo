from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np


def package_root():
    return Path(__file__).resolve().parent.parent


def load_frozen():
    root=package_root()/"frozen"
    scenarios=pd.read_csv(root/"transition_pair_scenarios.csv")
    pair_stage=pd.read_csv(root/"biological_pair_stage_summary_v1_1.csv")
    primary=pd.read_csv(root/"PRIMARY_ISI_STAGING.csv")
    boundaries=pd.read_csv(root/"staging_boundary_definitions_v1_1.csv")
    return scenarios,pair_stage,primary,boundaries


def select_scenarios(cfg):
    scenarios,pair_stage,primary,boundaries=load_frozen()
    mode=str(cfg['cohort']['mode'])
    coupled=pair_stage[pair_stage['path_family'].eq('coupled')].copy()
    if mode=='core_secure_all_support':
        pairs=set(coupled.loc[coupled['both_core_secure'].astype(bool),'biological_pair_key'])
        out=scenarios[scenarios['biological_pair_key'].isin(pairs)].copy()
    elif mode=='core_secure_best_only':
        pairs=set(coupled.loc[coupled['both_core_secure'].astype(bool),'biological_pair_key'])
        out=scenarios[scenarios['biological_pair_key'].isin(pairs) & scenarios['wt_source'].eq('best') & scenarios['sca_source'].eq('best')].copy()
    elif mode=='all_pairs_best_only':
        out=scenarios[scenarios['wt_source'].eq('best') & scenarios['sca_source'].eq('best')].copy()
    elif mode=='all_support':
        out=scenarios.copy()
    else:
        raise ValueError(f"Unknown cohort.mode={mode}")
    scenario_ids=cfg['cohort'].get('scenario_ids')
    if scenario_ids is not None:
        keep=set(int(x) for x in scenario_ids)
        out=out[out['scenario_id'].astype(int).isin(keep)].copy()
        missing=sorted(keep-set(out['scenario_id'].astype(int)))
        if missing: raise ValueError(f"Requested scenario_ids not present in selected cohort: {missing}")
    max_s=cfg['cohort'].get('max_scenarios')
    if max_s is not None: out=out.sort_values('scenario_id').head(int(max_s)).copy()
    sums=out.groupby('biological_pair_key')['within_pair_support_weight'].transform('sum')
    out['analysis_within_pair_weight']=out['within_pair_support_weight']/sums
    pairs=out['biological_pair_key'].nunique(); out['analysis_pair_weight']=1.0/max(pairs,1)
    out['analysis_scenario_weight']=out['analysis_pair_weight']*out['analysis_within_pair_weight']
    return out.reset_index(drop=True),coupled,primary,boundaries


def scenario_arrays(row,p_grid):
    p=np.asarray(p_grid,float)
    def lin(a,b): return (1-p)*float(a)+p*float(b)
    def loginterp(a,b):
        a=max(float(a),1e-15); b=max(float(b),1e-15)
        return np.exp((1-p)*np.log(a)+p*np.log(b))
    pars=np.empty((len(p),4),float)
    pars[:,0]=lin(row.wt_b,row.sca_b); pars[:,1]=loginterp(row.wt_r,row.sca_r)
    pars[:,2]=lin(row.wt_s,row.sca_s); pars[:,3]=loginterp(row.wt_kappa_I,row.sca_kappa_I)
    J=lin(row.wt_J_q75,row.sca_J_q75)
    return pars,J
