from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional
from .fit_window import compute_fit_window

REQUIRED_FROZEN = {
    'sweep_id','group','cell_id','sweep_index','current_pA','capacitance_pF','J',
    'stim_onset_ms','stim_offset_ms','final_audit_decision'
}
REQUIRED_EVENTS = {'group','cell_id','sweep_index','current_pA','onset_ms','offset_ms','time_ms'}
REQUIRED_BASELINE = {
    'group','cell_id','capacitance_pF','final_review_decision','primary_support',
    'b','r','s','kappa_I','cell_loss'
}
REQUIRED_THRESHOLD = {
    'group','cell_id','capacitance_pF',
    'nonspiking_sweep_index','nonspiking_current_pA','nonspiking_J',
    'nonspiking_onset_ms','nonspiking_offset_ms','nonspiking_fixed_qc_spike_count',
    'first_spiking_sweep_index','first_spiking_current_pA','first_spiking_J',
    'first_spiking_onset_ms','first_spiking_offset_ms','first_spiking_fixed_qc_spike_count',
    'threshold_bracket_width_pA'
}


def _truthy_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).astype(float) != 0
    vals = s.fillna('').astype(str).str.strip().str.lower()
    return vals.isin({'true','1','yes','y','pass','include','included'})


def _load_overrides(path: Optional[str]):
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError('Peak override file not found: %s' % p)
    df = pd.read_csv(p)
    if df.empty:
        return {}
    if not {'sweep_id','action'}.issubset(df.columns):
        raise ValueError('Peak override file requires sweep_id, action')
    out = {}
    for _, r in df.iterrows():
        sid = str(r['sweep_id']).strip()
        if sid in out:
            raise ValueError('Duplicate peak override for %s' % sid)
        out[sid] = str(r['action']).strip().upper()
    return out


def _apply_override(spikes: np.ndarray, action: Optional[str]):
    x = np.sort(np.asarray(spikes, dtype=float))
    if not action:
        return x, []
    applied = []
    for part in str(action).split('+'):
        part = part.strip().upper()
        if part == 'DROP_FIRST':
            if x.size < 2:
                raise ValueError('DROP_FIRST would leave no spikes')
            x = x[1:]
            applied.append('DROP_FIRST')
        elif part == 'DROP_LAST':
            if x.size < 2:
                raise ValueError('DROP_LAST would leave no spikes')
            x = x[:-1]
            applied.append('DROP_LAST')
        elif part in {'','NONE'}:
            pass
        else:
            raise ValueError('Unknown peak override action %r' % part)
    return x, applied


def _threshold_sweep(row, prefix):
    onset = float(row[prefix + '_onset_ms'])
    offset = float(row[prefix + '_offset_ms'])
    duration = offset - onset
    return {
        'sweep_id': '%s__threshold_%s' % (str(row['cell_id']), prefix),
        'group': str(row['group']),
        'cell_id': str(row['cell_id']),
        'sweep_index': int(row[prefix + '_sweep_index']),
        'current_pA': float(row[prefix + '_current_pA']),
        'capacitance_pF': float(row['capacitance_pF']),
        'J': float(row[prefix + '_J']),
        'stim_onset_ms': onset,
        'stim_offset_ms': offset,
        'stimulus_duration_ms': float(duration),
        # This is only a simulation horizon for binary spike presence, not a voltage-fit window.
        'fit_start_ms': 0.0,
        'fit_end_ms': float(duration),
        'exp_spike_times_ms': np.empty(0, dtype=float),
        'n_exp_spikes': int(row[prefix + '_fixed_qc_spike_count']),
        'threshold_only': True,
    }


def load_v3_6_cells(cfg):
    frozen_path = Path(cfg['data']['frozen_sweeps_manifest'])
    events_path = Path(cfg['data']['events_file'])
    baseline_path = Path(cfg['data']['baseline_cell_summary_file'])
    threshold_path = Path(cfg['data']['threshold_brackets_file'])
    for p in (frozen_path, events_path, baseline_path, threshold_path):
        if not p.exists():
            raise FileNotFoundError('Required v3.6 input not found: %s' % p)

    frozen = pd.read_csv(frozen_path)
    missing = sorted(REQUIRED_FROZEN - set(frozen.columns))
    if missing:
        raise ValueError('Frozen manifest missing columns: %s' % missing)
    if not (frozen['final_audit_decision'].astype(str).str.upper() == 'ACCEPT').all():
        raise ValueError('v3.6 frozen spiking manifest must contain ACCEPT sweeps only')
    if frozen['sweep_id'].duplicated().any():
        raise ValueError('Duplicate sweep_id in frozen v3.6 manifest')

    baseline = pd.read_csv(baseline_path)
    missing = sorted(REQUIRED_BASELINE - set(baseline.columns))
    if missing:
        raise ValueError('v3.1 baseline cell summary missing columns: %s' % missing)
    baseline = baseline[baseline['final_review_decision'].astype(str).str.upper() == 'ACCEPT'].copy()
    if baseline['cell_id'].duplicated().any():
        raise ValueError('Duplicate cell_id in v3.1 baseline cell summary')
    base_by = baseline.set_index('cell_id').to_dict('index')
    if set(frozen['cell_id'].astype(str)) != set(base_by):
        raise ValueError('Frozen v3.6 spiking sweep cells and accepted v3.1 cell summary do not match exactly')

    seed_by = {}
    seed_path_value = cfg['data'].get('seed_cell_summary_file')
    if seed_path_value:
        seed_path = Path(seed_path_value)
        if not seed_path.exists():
            raise FileNotFoundError('v3.8 seed cell summary not found: %s' % seed_path)
        seed_df = pd.read_csv(seed_path)
        required_seed = {'cell_id','b','r','s','kappa_I'}
        miss_seed = sorted(required_seed - set(seed_df.columns))
        if miss_seed:
            raise ValueError('v3.8 seed cell summary missing columns: %s' % miss_seed)
        if seed_df['cell_id'].duplicated().any():
            raise ValueError('Duplicate cell_id in v3.8 seed cell summary')
        seed_by = seed_df.set_index(seed_df['cell_id'].astype(str)).to_dict('index')
        if set(seed_by) != set(base_by):
            raise ValueError('v3.8 seed cells and frozen accepted cells do not match exactly')

    thresholds = pd.read_csv(threshold_path)
    missing = sorted(REQUIRED_THRESHOLD - set(thresholds.columns))
    if missing:
        raise ValueError('Threshold bracket manifest missing columns: %s' % missing)
    if thresholds['cell_id'].duplicated().any():
        raise ValueError('Threshold bracket manifest must have exactly one row per cell')
    if set(thresholds['cell_id'].astype(str)) != set(base_by):
        raise ValueError('Threshold bracket cells and accepted v3.1 cells do not match exactly')
    if not (pd.to_numeric(thresholds['nonspiking_fixed_qc_spike_count']) == 0).all():
        raise ValueError('Every nonspiking threshold sweep must have zero fixed-QC spikes')
    if not (pd.to_numeric(thresholds['first_spiking_fixed_qc_spike_count']) > 0).all():
        raise ValueError('Every first-spiking threshold sweep must have at least one fixed-QC spike')
    if not (pd.to_numeric(thresholds['first_spiking_current_pA']) > pd.to_numeric(thresholds['nonspiking_current_pA'])).all():
        raise ValueError('Invalid threshold bracket ordering')
    threshold_by = thresholds.set_index('cell_id').to_dict('index')

    events = pd.read_csv(events_path)
    missing = sorted(REQUIRED_EVENTS - set(events.columns))
    if missing:
        raise ValueError('Events file missing columns: %s' % missing)
    include_col = cfg['data'].get('events_include_column','fixed_qc_detected')
    if include_col not in events.columns:
        raise ValueError('Events file missing include column %r' % include_col)
    selected = events.loc[_truthy_series(events[include_col])].copy()
    for c in ('sweep_index','current_pA','onset_ms','offset_ms','time_ms'):
        selected[c] = pd.to_numeric(selected[c], errors='coerce')
    if selected[['sweep_index','current_pA','onset_ms','offset_ms','time_ms']].isna().any(axis=None):
        raise ValueError('Selected Stage-1 events contain nonnumeric required values')

    event_groups = {}
    for key, g in selected.groupby(['group','cell_id','sweep_index'], sort=False):
        event_groups[(str(key[0]), str(key[1]), int(key[2]))] = g
    overrides = _load_overrides(cfg['data'].get('peak_overrides_file'))
    tol = float(cfg['data'].get('metadata_tolerance_ms',0.05))
    strict = bool(cfg['data'].get('strict_metadata_match',True))

    sweeps = []
    audit = []
    for _, r in frozen.sort_values(['group','cell_id','sweep_index']).iterrows():
        key = (str(r['group']), str(r['cell_id']), int(r['sweep_index']))
        if key not in event_groups:
            raise ValueError('Frozen accepted spiking sweep missing selected Stage-1 events: %s' % (key,))
        g = event_groups[key]
        current = float(g['current_pA'].iloc[0]); onset = float(g['onset_ms'].iloc[0]); offset = float(g['offset_ms'].iloc[0])
        if strict:
            checks = [
                ('current_pA', current, float(r['current_pA']), 1e-6),
                ('onset_ms', onset, float(r['stim_onset_ms']), tol),
                ('offset_ms', offset, float(r['stim_offset_ms']), tol),
            ]
            for name,a,b,t in checks:
                if abs(a-b) > t:
                    raise ValueError('%s: frozen %s=%s disagrees with Stage-1=%s' % (r['sweep_id'],name,b,a))
        raw_rel = np.sort(g['time_ms'].to_numpy(dtype=float) - onset)
        action = overrides.get(str(r['sweep_id']))
        spikes, applied = _apply_override(raw_rel, action)
        stimulus_duration = offset - onset
        window = compute_fit_window(spikes, stimulus_duration, cfg)
        abf_path = ''
        if 'abf_path' in g.columns and g['abf_path'].notna().any():
            abf_path = str(g.loc[g['abf_path'].notna(),'abf_path'].iloc[0])
        sweep = {
            'sweep_id':str(r['sweep_id']), 'group':str(r['group']), 'cell_id':str(r['cell_id']),
            'sweep_index':int(r['sweep_index']), 'current_pA':float(r['current_pA']),
            'capacitance_pF':float(r['capacitance_pF']), 'J':float(r['J']),
            'stim_onset_ms':onset, 'stim_offset_ms':offset, 'stimulus_duration_ms':float(stimulus_duration),
            'fit_start_ms':float(window['fit_start_ms']), 'fit_end_ms':float(window['fit_end_ms']),
            'fit_guard_ms':float(window['guard_ms']), 'fit_local_isi_ms':window['local_isi_ms'],
            'last_exp_spike_ms':float(window['last_exp_spike_ms']), 'excluded_plateau_ms':float(window['excluded_plateau_ms']),
            'exp_spike_times_ms':np.asarray(spikes,dtype=float), 'n_exp_spikes':int(len(spikes)),
            'raw_reviewed_spike_count':int(len(raw_rel)), 'peak_override':'+'.join(applied) if applied else '',
            'abf_path':abf_path, 'threshold_only':False,
        }
        sweeps.append(sweep)
        audit.append({
            'sweep_id':sweep['sweep_id'],'group':sweep['group'],'cell_id':sweep['cell_id'],
            'sweep_index':sweep['sweep_index'],'current_pA':sweep['current_pA'],
            'raw_reviewed_spikes':sweep['raw_reviewed_spike_count'],'v3_6_spikes':sweep['n_exp_spikes'],
            'peak_override':sweep['peak_override'],'fit_end_ms':sweep['fit_end_ms'],
        })

    cells = []
    index = pd.DataFrame([{'i':i,'group':s['group'],'cell_id':s['cell_id']} for i,s in enumerate(sweeps)])
    for (group, cell_id), xs in index.groupby(['group','cell_id'], sort=True):
        sw = [sweeps[i] for i in xs['i'].tolist()]
        caps = {round(float(s['capacitance_pF']),9) for s in sw}
        if len(caps) != 1:
            raise ValueError('%s: capacitance not constant across frozen sweeps' % cell_id)
        b = base_by[str(cell_id)]
        tr = dict(threshold_by[str(cell_id)])
        tr['cell_id'] = str(cell_id)
        tr['group'] = str(group)
        if abs(float(tr['capacitance_pF']) - float(sw[0]['capacitance_pF'])) > 1e-6:
            raise ValueError('%s: threshold and spiking capacitance disagree' % cell_id)
        bracket = {
            'nonspiking_sweep': _threshold_sweep(tr, 'nonspiking'),
            'first_spiking_sweep': _threshold_sweep(tr, 'first_spiking'),
            'bracket_width_pA': float(tr['threshold_bracket_width_pA']),
            'selection_rule': str(tr.get('selection_rule','')),
        }
        seed = seed_by.get(str(cell_id), b)
        cells.append({
            'group':str(group), 'cell_id':str(cell_id), 'capacitance_pF':float(sw[0]['capacitance_pF']),
            'sweeps':sw, 'n_sweeps':len(sw), 'primary_support':str(b['primary_support']),
            'baseline_theta':{'b':float(seed['b']),'r':float(seed['r']),'s':float(seed['s']),'kappa_I':float(seed['kappa_I'])},
            'v3_8_seed':dict(seed), 'v3_1':dict(b), 'threshold_bracket':bracket,
        })
    return cells, pd.DataFrame(audit), baseline, thresholds
