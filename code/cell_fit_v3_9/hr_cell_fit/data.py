from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any, Dict, List
import numpy as np
import pandas as pd
from .fit_window import fit_window_end_ms

REQUIRED_EVENT_COLUMNS = {
    'group','cell_id','file','abf_sweep','current_pA','capacitance_pF','stimulus_onset_ms','stimulus_offset_ms','spike_time_ms','spike_source'
}


def _norm_file(x):
    if pd.isna(x): return ''
    return Path(str(x)).name


def _load_csv(path):
    return pd.read_csv(Path(path))


def _load_json(path):
    with open(path,'r',encoding='utf-8') as f:return json.load(f)


def load_calibration_bundle(cfg):
    dcfg=cfg['data']; include_col=dcfg['events_include_column']
    events=_load_csv(dcfg['events_file'])
    missing=REQUIRED_EVENT_COLUMNS-set(events.columns)
    if missing:raise ValueError(f'Events file missing required columns: {sorted(missing)}')
    if include_col not in events.columns:raise ValueError(f'Events file does not contain include column {include_col!r}')
    events=events[events[include_col].astype(bool)].copy()
    manifest=_load_csv(dcfg['frozen_sweeps_manifest'])
    overrides=_load_csv(dcfg['peak_overrides_file'])
    baseline_cells=_load_csv(dcfg['baseline_cell_summary_file'])
    baseline_sweeps=_load_csv(dcfg['baseline_sweep_summary_file'])
    baseline_ident=_load_csv(dcfg['baseline_identifiability_file']) if dcfg.get('baseline_identifiability_file') else pd.DataFrame()
    seeds=_load_csv(dcfg['seed_cell_summary_file']) if dcfg.get('seed_cell_summary_file') else baseline_cells.copy()
    threshold=_load_csv(dcfg['threshold_brackets_file'])
    return {'events':events,'manifest':manifest,'overrides':overrides,'baseline_cells':baseline_cells,'baseline_sweeps':baseline_sweeps,'baseline_ident':baseline_ident,'seeds':seeds,'threshold_brackets':threshold}


def _manifest_row(manifest,group,cell_id,file,sweep,current):
    m=manifest[(manifest.group==group)&(manifest.cell_id==cell_id)&(manifest.abf_sweep.astype(int)==int(sweep))]
    if 'file' in manifest.columns:
        m=m[m.file.map(_norm_file)==_norm_file(file)]
    if 'current_pA' in manifest.columns:
        m=m[np.isclose(m.current_pA.astype(float),float(current),atol=1e-6)]
    if len(m)!=1:raise ValueError(f'Frozen manifest expected one row for {cell_id} sweep {sweep}; got {len(m)}')
    return m.iloc[0]


def _metadata_for(events, group, cell_id, file, sweep, current, tol_ms):
    g=events[(events.group==group)&(events.cell_id==cell_id)&(events.abf_sweep.astype(int)==int(sweep))]
    g=g[g.file.map(_norm_file)==_norm_file(file)]
    g=g[np.isclose(g.current_pA.astype(float),float(current),atol=1e-6)]
    if g.empty:raise ValueError(f'No selected events for {cell_id} sweep {sweep}')
    vals={k:float(g[k].iloc[0]) for k in ['current_pA','capacitance_pF','stimulus_onset_ms','stimulus_offset_ms']}
    for k in ['capacitance_pF','stimulus_onset_ms','stimulus_offset_ms']:
        if float(g[k].max()-g[k].min())>tol_ms:raise ValueError(f'Inconsistent {k} in events for {cell_id} sweep {sweep}')
    return vals,g


def _apply_overrides(events_sweep,overrides,group,cell_id,file,sweep):
    times=events_sweep.spike_time_ms.astype(float).to_numpy()
    ov=overrides[(overrides.group==group)&(overrides.cell_id==cell_id)&(overrides.abf_sweep.astype(int)==int(sweep))]
    if 'file' in overrides.columns:
        ov=ov[ov.file.map(_norm_file)==_norm_file(file)]
    drop=[]; add=[]
    for _,r in ov.iterrows():
        action=str(r.get('action','')).strip().lower(); t=float(r['spike_time_ms'])
        if action=='drop': drop.append(t)
        elif action=='add': add.append(t)
        else: raise ValueError(f'Unsupported override action {action!r}')
    keep=np.ones(times.size,dtype=bool)
    for t in drop:
        if times.size==0:raise ValueError(f'Drop override {t} cannot match empty sweep')
        j=int(np.argmin(np.abs(times-t)))
        if abs(times[j]-t)>0.11:raise ValueError(f'Drop override {t} has no matching selected event; nearest={times[j]}')
        keep[j]=False
    final=np.sort(np.concatenate([times[keep],np.asarray(add,dtype=float)]))
    return final


def build_cells(cfg):
    b=load_calibration_bundle(cfg); events=b['events'];manifest=b['manifest'];over=b['overrides'];tol=cfg['data']['metadata_tolerance_ms']
    accepted=manifest[manifest.accepted.astype(bool)].copy() if 'accepted' in manifest.columns else manifest.copy()
    cells=[]
    for cell_id,mg in accepted.groupby('cell_id',sort=True):
        groups=mg.group.dropna().unique()
        if len(groups)!=1:raise ValueError(f'{cell_id}: ambiguous group {groups}')
        group=str(groups[0]); sweeps=[]
        for _,mr in mg.sort_values('current_pA').iterrows():
            file=str(mr['file']);sweep=int(mr['abf_sweep']);current=float(mr['current_pA'])
            md,eg=_metadata_for(events,group,cell_id,file,sweep,current,tol)
            _manifest_row(manifest,group,cell_id,file,sweep,current)
            spikes_abs=_apply_overrides(eg,over,group,cell_id,file,sweep)
            onset=float(md['stimulus_onset_ms']);off=float(md['stimulus_offset_ms'])
            rel=spikes_abs-onset
            rel=np.sort(rel[(rel>=-tol)&(rel<=off-onset+tol)])
            fit_end=fit_window_end_ms(rel,off-onset,cfg['loss']['fit_window'])
            exp_fit=rel[(rel>=0)&(rel<=fit_end+1e-9)]
            if exp_fit.size==0:raise ValueError(f'{cell_id} sweep {sweep}: frozen accepted spiking sweep has no fitted spikes')
            sweeps.append({'group':group,'cell_id':cell_id,'file':_norm_file(file),'abf_sweep':sweep,'current_pA':current,
                           'capacitance_pF':float(md['capacitance_pF']),'J':current/float(md['capacitance_pF']),
                           'stimulus_onset_ms':onset,'stimulus_offset_ms':off,'stimulus_duration_ms':off-onset,
                           'exp_spike_times_ms':exp_fit.astype(float),'all_exp_spike_times_ms':rel.astype(float),'fit_end_ms':float(fit_end)})
        cms={round(s['capacitance_pF'],8) for s in sweeps}
        if len(cms)!=1:raise ValueError(f'{cell_id}: capacitance differs across frozen accepted sweeps')
        cells.append({'group':group,'cell_id':cell_id,'capacitance_pF':sweeps[0]['capacitance_pF'],'sweeps':sweeps})
    return cells,b


def find_seed_row(bundle,cell_id):
    df=bundle['seeds']; row=df[df.cell_id==cell_id]
    return row.iloc[0] if len(row)==1 else None


def find_baseline_cell_row(bundle,cell_id):
    df=bundle['baseline_cells']; row=df[df.cell_id==cell_id]
    return row.iloc[0] if len(row)==1 else None


def find_threshold_bracket(bundle,cell_id):
    df=bundle['threshold_brackets']; row=df[df.cell_id==cell_id]
    if len(row)!=1:return None
    return row.iloc[0]


def validate_bundle(cfg):
    cells,b=build_cells(cfg); errors=[];warnings=[]
    n_sweeps=sum(len(c['sweeps']) for c in cells); n_spikes=sum(len(s['exp_spike_times_ms']) for c in cells for s in c['sweeps'])
    expected_cells=20;expected_sweeps=113;expected_spikes=4884
    if len(cells)!=expected_cells:errors.append(f'Expected {expected_cells} frozen accepted cells, got {len(cells)}')
    if n_sweeps!=expected_sweeps:errors.append(f'Expected {expected_sweeps} frozen accepted spiking sweeps, got {n_sweeps}')
    if n_spikes!=expected_spikes:errors.append(f'Expected {expected_spikes} fitted spikes, got {n_spikes}')
    for c in cells:
        if find_seed_row(b,c['cell_id']) is None:errors.append(f"{c['cell_id']}: missing v3.8 seed row")
        if find_threshold_bracket(b,c['cell_id']) is None:errors.append(f"{c['cell_id']}: missing threshold bracket")
    bb=cfg['bounds']
    expected={'b':(0.5,7.0,'linear'),'r':(0.0001,0.10,'log'),'s':(0.05,15.0,'linear'),'kappa_I':(0.0002,2.0,'log')}
    for p,(lo,hi,sc) in expected.items():
        x=bb[p]
        if not (abs(float(x['min'])-lo)<1e-12 and abs(float(x['max'])-hi)<1e-12 and str(x['scale'])==sc):errors.append(f'{p}: unexpected v3.9 bound {x}')
    return {'status':'PASS' if not errors else 'FAIL','errors':errors,'warnings':warnings,'n_cells':len(cells),'n_spiking_sweeps':n_sweeps,'n_spikes':n_spikes}
