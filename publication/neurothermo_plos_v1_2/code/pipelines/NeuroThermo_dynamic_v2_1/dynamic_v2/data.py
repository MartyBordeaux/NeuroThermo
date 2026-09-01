from __future__ import annotations
from pathlib import Path
import pandas as pd, numpy as np

def load_inputs(cfg):
    root=Path(cfg['data']['root'])
    primary=pd.read_csv(root/cfg['data']['primary_cell_master'])
    manifest=pd.read_csv(root/cfg['data']['spiking_manifest'])
    events=pd.read_csv(root/cfg['data']['selected_spike_events'])
    thresholds=pd.read_csv(root/cfg['data']['threshold_brackets'])
    alts=pd.read_csv(root/cfg['data']['identifiability_alternatives'])
    ids=pd.read_csv(root/cfg['data']['animal_id_map'])
    primary=primary[primary['analysis_set'].eq('PRIMARY_MULTI_SWEEP')].copy()
    return primary,manifest,events,thresholds,alts,ids

def theta_from_row(r):
    return {k:float(r[k]) for k in ['b','r','s','kappa_I']}

def selected_spikes(events,sweep_id):
    x=events[events.sweep_id.astype(str).eq(str(sweep_id))].sort_values('spike_number')
    return x.time_from_onset_ms.to_numpy(float)

def near_optimal_sets(alts, cell_id):
    x=alts[(alts.cell_id.astype(str)==str(cell_id)) & (alts.near_optimal==True)].copy()
    seen=set(); out=[]
    for _,r in x.iterrows():
        theta=theta_from_row(r); key=tuple(round(theta[k],10) for k in ['b','r','s','kappa_I'])
        if key in seen: continue
        seen.add(key)
        out.append({'source':'%s_%s'%(r['parameter'],r['side']),'loss':float(r['loss']),'theta':theta})
    return out
