from __future__ import annotations
import argparse, json
from pathlib import Path
import yaml
from .analysis import run_all
from .plotting import make_plots


def load_config(path):
    cfg=yaml.safe_load(Path(path).read_text())
    base=Path(path).resolve().parent.parent
    root=Path(cfg['data']['root'])
    if not root.is_absolute(): cfg['data']['root']=str((base/root).resolve())
    out=Path(cfg['output']['dir'])
    if not out.is_absolute(): cfg['output']['dir']=str((base/out).resolve())
    return cfg


def validate(cfg):
    root=Path(cfg['data']['root'])
    req=[cfg['data'][k] for k in ['primary_cell_master','spiking_manifest','selected_spike_events','threshold_brackets','identifiability_alternatives','animal_id_map']]
    missing=[x for x in req if not (root/x).exists()]
    if missing: raise SystemExit('Missing frozen inputs: '+', '.join(missing))
    import pandas as pd
    p=pd.read_csv(root/cfg['data']['primary_cell_master']); m=pd.read_csv(root/cfg['data']['spiking_manifest']); e=pd.read_csv(root/cfg['data']['selected_spike_events'])
    pp=p[p.analysis_set.eq('PRIMARY_MULTI_SWEEP')]
    assert len(pp)==18 and (pp.group=='WT').sum()==12 and (pp.group=='SCA3').sum()==6
    assert len(m)==113 and len(e)==4884
    primary_ids=set(pp.cell_id.astype(str)); primary_sweeps=int(m.cell_id.astype(str).isin(primary_ids).sum()); primary_spikes=int(e.cell_id.astype(str).isin(primary_ids).sum())
    print(json.dumps({'version':'2.1.0','mode':'experimental-support-restricted','primary_cells':18,'WT':12,'SCA3':6,
                      'primary_spiking_sweeps':primary_sweeps,'primary_selected_spikes':primary_spikes,
                      'all_accepted_spiking_sweeps':113,'all_accepted_selected_spikes':4884,
                      'q_targets':cfg['analysis']['q_targets'],'extrapolation_allowed':False,
                      'same_current_alternative_comparison':True},indent=2))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('command',choices=['validate','run','plot']); ap.add_argument('--config',required=True); a=ap.parse_args(); cfg=load_config(a.config)
    if a.command=='validate': validate(cfg)
    elif a.command=='run':
        validate(cfg); s=run_all(cfg); make_plots(cfg['output']['dir']); print(json.dumps(s,indent=2))
    else: make_plots(cfg['output']['dir'])
if __name__=='__main__': main()
