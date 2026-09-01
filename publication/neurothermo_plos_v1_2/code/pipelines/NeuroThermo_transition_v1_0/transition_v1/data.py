from __future__ import annotations
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

REQUIRED_SUPPORT_COLS={
    'group','cell_id','animal_id','solution_key','source','b','r','s','kappa_I',
    'rheobase_J','q75_firing_rate_hz','q75_mean_isi_ms',
    'cell_weight_within_group','within_cell_support_weight','group_support_weight'
}


def load_csv(root,name):
    p=Path(root)/name
    if not p.exists(): raise FileNotFoundError(p)
    return pd.read_csv(p)


def load_inputs(cfg):
    root=Path(cfg['data']['root'])
    support=load_csv(root,cfg['data']['endpoint_support'])
    cells=load_csv(root,cfg['data']['endpoint_cells'])
    anchors=load_csv(root,cfg['data']['q75_anchors'])
    transform=load_csv(root,cfg['data']['experimental_transform'])
    geometry=load_csv(root,cfg['data']['experimental_geometry'])
    missing=REQUIRED_SUPPORT_COLS-set(support.columns)
    if missing: raise ValueError('endpoint support missing columns: '+','.join(sorted(missing)))
    return support,cells,anchors,transform,geometry


def theta_from_row(r):
    return {k:float(r[k]) for k in ['b','r','s','kappa_I']}


def build_support_with_anchors(support,anchors):
    a=anchors[['group','cell_id','best_rheobase_J','J_max_observed','J_q75','active_support_q75_ms','q75_supported']].copy()
    out=support.merge(a,on=['group','cell_id'],how='left',validate='many_to_one')
    if out[['J_q75','active_support_q75_ms']].isna().any().any():
        bad=out[out[['J_q75','active_support_q75_ms']].isna().any(axis=1)][['group','cell_id']].drop_duplicates()
        raise ValueError('missing q75 anchors for '+bad.to_dict('records').__repr__())
    return out


def enumerate_scenarios(support):
    wt=support[support.group.eq('WT')].copy()
    sc=support[support.group.eq('SCA3')].copy()
    rows=[]; sid=0
    for _,w in wt.iterrows():
        for _,s in sc.iterrows():
            rows.append({
                'scenario_id':sid,
                'wt_cell_id':w.cell_id,'sca_cell_id':s.cell_id,
                'wt_animal_id':w.animal_id,'sca_animal_id':s.animal_id,
                'wt_solution_key':w.solution_key,'sca_solution_key':s.solution_key,
                'wt_source':w.source,'sca_source':s.source,
                'wt_b':w.b,'wt_r':w.r,'wt_s':w.s,'wt_kappa_I':w.kappa_I,
                'sca_b':s.b,'sca_r':s.r,'sca_s':s.s,'sca_kappa_I':s.kappa_I,
                'wt_rheobase_J_endpoint':w.rheobase_J,'sca_rheobase_J_endpoint':s.rheobase_J,
                'wt_J_q75':w.J_q75,'sca_J_q75':s.J_q75,
                'wt_active_support_ms':w.active_support_q75_ms,'sca_active_support_ms':s.active_support_q75_ms,
                'wt_within_cell_support_weight':w.within_cell_support_weight,
                'sca_within_cell_support_weight':s.within_cell_support_weight,
                'wt_group_support_weight':w.group_support_weight,
                'sca_group_support_weight':s.group_support_weight,
                'biological_pair_weight':1.0/(wt.cell_id.nunique()*sc.cell_id.nunique()),
                'within_pair_support_weight':float(w.within_cell_support_weight)*float(s.within_cell_support_weight),
                'scenario_weight':float(w.group_support_weight)*float(s.group_support_weight),
                'biological_pair_key':f'{w.cell_id}__TO__{s.cell_id}',
            })
            sid+=1
    df=pd.DataFrame(rows)
    return df


def sha256_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''): h.update(chunk)
    return h.hexdigest()


def frozen_manifest(root):
    root=Path(root)
    return {p.name:sha256_file(p) for p in sorted(root.glob('*')) if p.is_file()}
