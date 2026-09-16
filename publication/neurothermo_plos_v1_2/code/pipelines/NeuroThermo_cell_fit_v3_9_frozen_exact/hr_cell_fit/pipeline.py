from __future__ import annotations
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import yaml
from .data import load_v3_6_cells
from .optimize import fit_cell
from .diagnostics import classify_sweep, classify_cell
from .identifiability import practical_identifiability
from .objective import evaluate_theta
from .plotting import plot_cell, make_audit_pdf


def _jsonable(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, dict): return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [_jsonable(v) for v in x]
    if hasattr(x, '__dict__'): return _jsonable(x.__dict__)
    return x


def _file_sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(1024 * 1024)
            if not b: break
            h.update(b)
    return h.hexdigest()


def run_signature(cfg):
    core = {k: v for k, v in cfg.items() if k != 'output'}
    files = {}
    for key in (
        'events_file', 'frozen_sweeps_manifest', 'peak_overrides_file',
        'baseline_cell_summary_file', 'baseline_sweep_summary_file',
        'baseline_identifiability_file', 'seed_cell_summary_file', 'threshold_brackets_file',
    ):
        p = cfg['data'].get(key)
        files[key] = _file_sha(p) if p and Path(p).exists() else None
    payload = json.dumps({'config': core, 'files': files}, sort_keys=True, default=str).encode()
    return hashlib.sha256(payload).hexdigest()


def validate(cfg):
    cells, audit, baseline, thresholds = load_v3_6_cells(cfg)
    primary = [c for c in cells if c['primary_support'] == 'MULTI_SWEEP']
    widths = pd.to_numeric(thresholds['threshold_bracket_width_pA'])
    la = cfg['loss']['latency_alignment']
    return {
        'status': 'validated',
        'version': '3.9.0',
        'boundary_stress_parameter': 's',
        'v3_8_s_lower_bound': 0.25,
        'v3_9_s_lower_bound': float(cfg['bounds']['s']['min']),
        'v3_9_s_upper_bound': float(cfg['bounds']['s']['max']),
        'frozen_v3_1_accepted_cells': len(cells),
        'optimizer_seed_cells': int(sum('v3_8_seed' in c for c in cells)),
        'primary_multisweep_cells': len(primary),
        'single_sweep_cells': [c['cell_id'] for c in cells if c['primary_support'] == 'SINGLE_SWEEP_ONLY'],
        'n_spiking_fit_sweeps': int(sum(len(c['sweeps']) for c in cells)),
        'n_spikes_after_overrides': int(audit['v3_6_spikes'].sum()),
        'n_peak_overrides': int((audit['peak_override'].fillna('') != '').sum()),
        'n_threshold_brackets': int(len(thresholds)),
        'threshold_bracket_widths_pA': sorted(set(float(x) for x in widths)),
        'threshold_rule': 'binary presence only: zero spikes at highest nonspiking step; >=1 spike at first spiking step',
        'latency_alignment_enabled': bool(la.get('enabled', True)),
        'latency_alignment_method': str(la.get('method', 'exact_first_spike')),
        'latency_alignment_rule': 'tau = first experimental spike - first raw model spike; no tau optimization; no post-shift clipping; recomputed for every theta including identifiability alternatives',
        'alignment_preserves_model_spike_count': True,
        'explicit_latency_penalty': False,
        'explicit_isi_loss': False,
        'voltage_plateau_fit': False,
        'groups': {
            g: {
                'n_cells': int(sum(c['group'] == g for c in cells)),
                'n_spiking_fit_sweeps': int(sum(len(c['sweeps']) for c in cells if c['group'] == g)),
            }
            for g in sorted({c['group'] for c in cells})
        },
    }


def _result_rows(cell, fit, cfg):
    rows = []
    large_thr = float(cfg['loss']['latency_alignment'].get('large_shift_warning_ms', 200.0))
    for sweep, ev in zip(cell['sweeps'], fit['sweep_evals']):
        exp = np.asarray(sweep['exp_spike_times_ms'], dtype=float)
        decision, reason, d = classify_sweep(ev.vp_loss, exp, ev.model_spikes, cfg)
        _, _, raw_d = classify_sweep(ev.raw_vp_loss, exp, ev.raw_model_spikes, cfg)
        aligned_lat = float(ev.model_spikes[0]) if len(ev.model_spikes) else np.nan
        raw_lat = float(ev.raw_model_spikes[0]) if len(ev.raw_model_spikes) else np.nan
        exp_lat = float(exp[0]) if len(exp) else np.nan
        exp_last = float(exp[-1]) if len(exp) else np.nan
        raw_last = float(ev.raw_model_spikes[-1]) if len(ev.raw_model_spikes) else np.nan
        aligned_last = float(ev.model_spikes[-1]) if len(ev.model_spikes) else np.nan
        exp_duration = float(exp[-1]-exp[0]) if len(exp) >= 2 else np.nan
        model_duration = float(ev.raw_model_spikes[-1]-ev.raw_model_spikes[0]) if len(ev.raw_model_spikes) >= 2 else np.nan
        rows.append({
            'sweep_id': sweep['sweep_id'], 'group': cell['group'], 'cell_id': cell['cell_id'],
            'sweep_index': sweep['sweep_index'], 'current_pA': sweep['current_pA'], 'J': sweep['J'],
            'peak_override': sweep['peak_override'], 'fit_end_ms': sweep['fit_end_ms'],
            'joint_sweep_decision': decision, 'joint_sweep_reason': reason,
            **d,
            'latency_shift_ms': float(ev.latency_shift_ms),
            'abs_latency_shift_ms': abs(float(ev.latency_shift_ms)),
            'large_latency_shift': bool(abs(float(ev.latency_shift_ms)) >= large_thr),
            'latency_alignment_applied': bool(ev.latency_alignment_applied),
            'count_preserved_by_alignment': bool(ev.count_preserved_by_alignment),
            'raw_vp_loss': float(ev.raw_vp_loss),
            'raw_count_error_fraction': float(ev.raw_count_error_fraction),
            'raw_n_model_spikes': int(len(ev.raw_model_spikes)),
            'aligned_spikes_beyond_fit_end': int(np.sum(np.asarray(ev.model_spikes, dtype=float) > float(sweep['fit_end_ms']) + 1e-9)),
            'aligned_spikes_before_zero': int(np.sum(np.asarray(ev.model_spikes, dtype=float) < -1e-9)),
            'raw_model_latency_ms': raw_lat,
            'aligned_model_latency_ms': aligned_lat,
            'aligned_latency_error_ms': abs(aligned_lat-exp_lat) if np.isfinite(aligned_lat) and np.isfinite(exp_lat) else np.nan,
            'raw_latency_error_ms': abs(raw_lat-exp_lat) if np.isfinite(raw_lat) and np.isfinite(exp_lat) else np.nan,
            'exp_last_spike_ms': exp_last,
            'raw_model_last_spike_ms': raw_last,
            'aligned_model_last_spike_ms': aligned_last,
            'aligned_last_spike_error_ms': abs(aligned_last-exp_last) if np.isfinite(aligned_last) and np.isfinite(exp_last) else np.nan,
            'signed_aligned_last_spike_error_ms': aligned_last-exp_last if np.isfinite(aligned_last) and np.isfinite(exp_last) else np.nan,
            'exp_train_duration_ms': exp_duration,
            'model_train_duration_ms': model_duration,
            'train_duration_error_ms': model_duration-exp_duration if np.isfinite(model_duration) and np.isfinite(exp_duration) else np.nan,
            'last_spike_used_as_anchor': False,
            'time_rescaling_applied': False,
            'raw_circle_f1': float(raw_d['circle_f1']),
            'raw_circle_recall': float(raw_d['circle_recall']),
            'raw_circle_precision': float(raw_d['circle_precision']),
            'model_spike_times_ms': np.asarray(ev.model_spikes, dtype=float),
            'raw_model_spike_times_ms': np.asarray(ev.raw_model_spikes, dtype=float),
        })
    return rows


def _threshold_row(cell, threshold_eval):
    br = cell['threshold_bracket']
    lo = br['nonspiking_sweep']; hi = br['first_spiking_sweep']
    return {
        'group': cell['group'], 'cell_id': cell['cell_id'],
        'nonspiking_sweep_index': lo['sweep_index'], 'nonspiking_current_pA': lo['current_pA'],
        'first_spiking_sweep_index': hi['sweep_index'], 'first_spiking_current_pA': hi['current_pA'],
        'threshold_bracket_width_pA': br['bracket_width_pA'],
        'model_spikes_at_nonspiking_current': int(len(threshold_eval.nonspiking_model_spikes)),
        'model_spikes_at_first_spiking_current': int(len(threshold_eval.first_spiking_model_spikes)),
        'nonspiking_constraint_pass': bool(not threshold_eval.nonspiking_violation),
        'first_spiking_constraint_pass': bool(not threshold_eval.first_spiking_violation),
        'threshold_pass': bool(threshold_eval.pass_constraint),
        'threshold_loss': float(threshold_eval.total_penalty),
        'nonspiking_model_spike_times_ms': np.asarray(threshold_eval.nonspiking_model_spikes, dtype=float),
        'first_spiking_model_spike_times_ms': np.asarray(threshold_eval.first_spiking_model_spikes, dtype=float),
    }


def _fit_one(payload):
    cell, cfg, seed_offset, sig = payload
    fit = fit_cell(cell, cfg, seed_offset)
    sweep_rows = _result_rows(cell, fit, cfg)
    spike_decision, cdiag = classify_cell(sweep_rows, cfg)
    te = fit['threshold_eval']
    threshold_pass = bool(te.pass_constraint)
    cdec = spike_decision
    if bool(cfg.get('threshold_constraint', {}).get('require_pass_for_auto_accept', True)) and cdec == 'ACCEPT' and not threshold_pass:
        cdec = 'REVIEW'

    theta = fit['theta']
    old31 = cell['v3_1']
    old38 = cell['v3_8_seed']
    seed_eval = evaluate_theta(cell['baseline_theta'], cell, cfg, float(cfg['optimization']['dt_fine_ms']), float(cfg['loss']['vp_tau_ms']), 'fine')
    shifts = np.asarray([float(x['latency_shift_ms']) for x in sweep_rows], dtype=float)
    raw_vps = np.asarray([float(x['raw_vp_loss']) for x in sweep_rows], dtype=float)
    aligned_vps = np.asarray([float(x['vp_loss']) for x in sweep_rows], dtype=float)
    row = {
        'group': cell['group'], 'cell_id': cell['cell_id'], 'capacitance_pF': cell['capacitance_pF'],
        'n_frozen_spiking_sweeps': len(cell['sweeps']), 'frozen_v3_1_decision': 'ACCEPT',
        'post_v3_9_auto_cell_decision': cdec, 'spike_only_auto_cell_decision': spike_decision,
        'review': '', 'comment': '', 'primary_support': cell['primary_support'],
        'cell_loss': float(fit['loss']), 'spike_train_loss': float(fit['spike_train_loss']),
        'threshold_loss': float(fit['threshold_loss']), **cdiag,
        'mean_raw_vp_loss': float(np.nanmean(raw_vps)),
        'mean_aligned_vp_loss': float(np.nanmean(aligned_vps)),
        'median_abs_latency_shift_ms': float(np.nanmedian(np.abs(shifts))),
        'max_abs_latency_shift_ms': float(np.nanmax(np.abs(shifts))),
        'n_large_latency_shift_sweeps': int(sum(bool(x['large_latency_shift']) for x in sweep_rows)),
        'n_alignment_count_changes': int(sum(not bool(x['count_preserved_by_alignment']) for x in sweep_rows)),
        'threshold_nonspiking_current_pA': float(cell['threshold_bracket']['nonspiking_sweep']['current_pA']),
        'threshold_first_spiking_current_pA': float(cell['threshold_bracket']['first_spiking_sweep']['current_pA']),
        'threshold_bracket_width_pA': float(cell['threshold_bracket']['bracket_width_pA']),
        'threshold_nonspiking_model_spikes': int(len(te.nonspiking_model_spikes)),
        'threshold_first_spiking_model_spikes': int(len(te.first_spiking_model_spikes)),
        'threshold_nonspiking_pass': bool(not te.nonspiking_violation),
        'threshold_first_spiking_pass': bool(not te.first_spiking_violation),
        'threshold_pass': threshold_pass,
        'b': float(theta['b']), 'r': float(theta['r']), 's': float(theta['s']), 'kappa_I': float(theta['kappa_I']),
        'kappa_over_Cm': float(theta['kappa_I']) / float(cell['capacitance_pF']),
        's_bound_min': float(cfg['bounds']['s']['min']), 's_bound_max': float(cfg['bounds']['s']['max']),
        's_bound_position_fraction': (float(theta['s']) - float(cfg['bounds']['s']['min'])) / (float(cfg['bounds']['s']['max']) - float(cfg['bounds']['s']['min'])),
        's_below_v3_8_lower_bound': bool(float(theta['s']) < 0.25),
        's_near_lower_1pct': bool((float(theta['s']) - float(cfg['bounds']['s']['min'])) / (float(cfg['bounds']['s']['max']) - float(cfg['bounds']['s']['min'])) <= 0.01),
        's_near_lower_2pct': bool((float(theta['s']) - float(cfg['bounds']['s']['min'])) / (float(cfg['bounds']['s']['max']) - float(cfg['bounds']['s']['min'])) <= 0.02),
        's_near_lower_5pct': bool((float(theta['s']) - float(cfg['bounds']['s']['min'])) / (float(cfg['bounds']['s']['max']) - float(cfg['bounds']['s']['min'])) <= 0.05),
        'best_label': fit['best_label'], 'global_rescue_used': bool(fit['rescue_used']), 'elapsed_s': float(fit['elapsed_s']),
        'seed_source': str(old38.get('seed_source', 'frozen_seed')), 'seed_cell_loss': float(old38.get('cell_loss', np.nan)), 'seed_b': float(old38['b']), 'seed_r': float(old38['r']),
        'seed_s': float(old38['s']), 'seed_kappa_I': float(old38['kappa_I']),
        'seed_on_v3_9_objective': float(seed_eval.loss),
        'seed_on_v3_9_spike_train_loss': float(seed_eval.spike_train_loss),
        'seed_on_v3_9_threshold_loss': float(seed_eval.threshold_loss),
        'delta_total_loss_vs_seed_same_objective': float(fit['loss']) - float(seed_eval.loss),
        'delta_spike_train_loss_vs_seed_same_objective': float(fit['spike_train_loss']) - float(seed_eval.spike_train_loss),
        'v3_1_cell_loss': float(old31['cell_loss']), 'v3_1_b': float(old31['b']), 'v3_1_r': float(old31['r']),
        'v3_1_s': float(old31['s']), 'v3_1_kappa_I': float(old31['kappa_I']),
        'v3_8_cell_loss': float(old38.get('cell_loss', np.nan)), 'v3_8_b': float(old38['b']), 'v3_8_r': float(old38['r']),
        'v3_8_s': float(old38['s']), 'v3_8_kappa_I': float(old38['kappa_I']),
        'delta_b_vs_seed': float(theta['b']) - float(old38['b']),
        'delta_r_vs_seed': float(theta['r']) - float(old38['r']),
        'delta_s_vs_seed': float(theta['s']) - float(old38['s']),
        'delta_kappa_I_vs_seed': float(theta['kappa_I']) - float(old38['kappa_I']),
        'identifiability': 'PENDING_REVIEW' if cell['primary_support'] == 'MULTI_SWEEP' else 'INSUFFICIENT_SPIKING_SWEEPS',
    }
    return {
        'run_signature': sig,
        'fit': {k: v for k, v in fit.items() if k not in ('sweep_evals','threshold_eval')},
        'cell_row': row, 'sweep_rows': sweep_rows, 'threshold_row': _threshold_row(cell, te),
    }


def _load_checkpoint(path, sig):
    try:
        obj = json.loads(Path(path).read_text())
    except Exception:
        return None
    return obj if obj.get('run_signature') == sig else None


def run(cfg):
    cells, audit, baseline, thresholds = load_v3_6_cells(cfg)
    sig = run_signature(cfg)
    out = Path(cfg['output']['dir']); out.mkdir(parents=True, exist_ok=True)
    per = out / 'per_cell'; per.mkdir(exist_ok=True)
    plots = out / 'plots'; plots.mkdir(exist_ok=True)

    audit.to_csv(out / 'frozen_spiking_input_audit.csv', index=False)
    thresholds.to_csv(out / 'frozen_threshold_brackets_v3_9.csv', index=False)
    (out / 'run_signature.txt').write_text(sig + '\n')
    with (out / 'resolved_config.yaml').open('w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    for key in (
        'frozen_sweeps_manifest', 'peak_overrides_file', 'baseline_cell_summary_file',
        'baseline_sweep_summary_file', 'baseline_identifiability_file', 'seed_cell_summary_file', 'threshold_brackets_file',
    ):
        p = cfg['data'].get(key)
        if p: shutil.copy2(p, out / Path(p).name)

    results = []; pending = []
    for i, cell in enumerate(cells):
        cp = per / ('%s.json' % cell['cell_id'])
        obj = _load_checkpoint(cp, sig) if bool(cfg['output'].get('resume', True)) and cp.exists() else None
        if obj is not None: results.append(obj)
        else: pending.append((i, cell, cp))

    n_jobs = max(1, int(cfg['optimization'].get('n_jobs', 1)))
    if n_jobs == 1:
        for i, cell, cp in pending:
            obj = _jsonable(_fit_one((cell, cfg, i, sig)))
            cp.write_text(json.dumps(obj, indent=2))
            results.append(obj)
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            fut = {ex.submit(_fit_one, (cell, cfg, i, sig)): (cell, cp) for i, cell, cp in pending}
            for f in as_completed(fut):
                cell, cp = fut[f]
                obj = _jsonable(f.result())
                cp.write_text(json.dumps(obj, indent=2))
                results.append(obj)

    results.sort(key=lambda o: (o['cell_row']['group'], o['cell_row']['cell_id']))
    cell_df = pd.DataFrame([o['cell_row'] for o in results])
    sweep_rows = []
    for o in results: sweep_rows.extend(o['sweep_rows'])
    sweep_df = pd.DataFrame(sweep_rows)
    threshold_df = pd.DataFrame([o['threshold_row'] for o in results])
    for col in ('model_spike_times_ms','raw_model_spike_times_ms'):
        if col in sweep_df:
            sweep_df[col] = sweep_df[col].apply(lambda x: json.dumps(_jsonable(x)))
    for col in ('nonspiking_model_spike_times_ms','first_spiking_model_spike_times_ms'):
        if col in threshold_df:
            threshold_df[col] = threshold_df[col].apply(lambda x: json.dumps(_jsonable(x)))
    cell_df.to_csv(out / 'cell_fit_summary.csv', index=False)
    sweep_df.to_csv(out / 'sweep_fit_summary.csv', index=False)
    threshold_df.to_csv(out / 'threshold_constraint_summary.csv', index=False)

    # Dedicated s-boundary stress table versus the frozen v3.8 solution.
    stress_cols = ['group','cell_id','primary_support','cell_loss','v3_8_cell_loss','s','v3_8_s','delta_s_vs_seed',
                   's_bound_min','s_bound_max','s_bound_position_fraction','s_below_v3_8_lower_bound',
                   's_near_lower_1pct','s_near_lower_2pct','s_near_lower_5pct','post_v3_9_auto_cell_decision','review','comment']
    stress_df = cell_df[[c for c in stress_cols if c in cell_df.columns]].copy()
    if 'cell_loss' in stress_df.columns and 'v3_8_cell_loss' in stress_df.columns:
        stress_df['delta_cell_loss_vs_v3_8'] = pd.to_numeric(stress_df['cell_loss'], errors='coerce') - pd.to_numeric(stress_df['v3_8_cell_loss'], errors='coerce')
    stress_df.to_csv(out / 's_boundary_stress_summary.csv', index=False)

    latency_cols = [
        'sweep_id','group','cell_id','sweep_index','current_pA','latency_shift_ms','abs_latency_shift_ms',
        'large_latency_shift','latency_alignment_applied','count_preserved_by_alignment','exp_latency_ms','raw_model_latency_ms',
        'aligned_model_latency_ms','raw_latency_error_ms','aligned_latency_error_ms','raw_vp_loss','vp_loss',
        'raw_circle_f1','circle_f1','raw_n_model_spikes','n_model_spikes','exp_last_spike_ms','raw_model_last_spike_ms','aligned_model_last_spike_ms','aligned_last_spike_error_ms','signed_aligned_last_spike_error_ms','exp_train_duration_ms','model_train_duration_ms','train_duration_error_ms','last_spike_used_as_anchor','time_rescaling_applied','aligned_spikes_beyond_fit_end','aligned_spikes_before_zero',
    ]
    sweep_df[[c for c in latency_cols if c in sweep_df.columns]].to_csv(out / 'latency_alignment_summary.csv', index=False)

    sol = []
    for o in results:
        cr = o['cell_row']
        for s in o['fit'].get('solutions', []):
            sol.append({'group': cr['group'], 'cell_id': cr['cell_id'], **s})
    pd.DataFrame(sol).to_csv(out / 'optimizer_solutions.csv', index=False)

    primary = cell_df[cell_df['primary_support'] == 'MULTI_SWEEP']
    gs = []
    for group, g in primary.groupby('group'):
        rr = {'group': group, 'n_cells': len(g), 'threshold_pass_fraction': float(g['threshold_pass'].mean())}
        for p in ('b','r','s','kappa_I','kappa_over_Cm','spike_train_loss','cell_loss','median_abs_latency_shift_ms'):
            rr['median_' + p] = g[p].median(); rr['q25_' + p] = g[p].quantile(.25); rr['q75_' + p] = g[p].quantile(.75)
        gs.append(rr)
    pd.DataFrame(gs).to_csv(out / 'group_parameter_summary_frozen_cohort.csv', index=False)

    cells_by = {c['cell_id']: c for c in cells}
    sr_by = {o['cell_row']['cell_id']: o['sweep_rows'] for o in results}
    tr_by = {o['cell_row']['cell_id']: o['threshold_row'] for o in results}
    if bool(cfg['output'].get('make_plots', True)):
        for row in cell_df.to_dict('records'):
            plot_cell(cells_by[row['cell_id']], sr_by[row['cell_id']], tr_by[row['cell_id']], row,
                      plots / ('%s_joint_fit_v3_9.png' % row['cell_id']))
    if bool(cfg['output'].get('make_audit_pdf', True)):
        make_audit_pdf(cells_by, sr_by, tr_by, cell_df.to_dict('records'), out / 'joint_fit_visual_audit_v3_9.pdf')

    summary = {
        'n_cells': len(cell_df), 'n_spiking_fit_sweeps': len(sweep_df), 'n_threshold_brackets': len(threshold_df),
        'post_v3_9_auto_decisions': cell_df['post_v3_9_auto_cell_decision'].value_counts().to_dict(),
        'threshold_pass': cell_df['threshold_pass'].value_counts().to_dict(),
        'primary_multisweep_cells': int((cell_df.primary_support == 'MULTI_SWEEP').sum()),
        'median_abs_latency_shift_ms': float(sweep_df['abs_latency_shift_ms'].median()),
        'large_latency_shift_sweeps': int(sweep_df['large_latency_shift'].sum()),
        'alignment_count_changes': int((~sweep_df['count_preserved_by_alignment'].astype(bool)).sum()),
        'max_aligned_latency_error_ms': float(pd.to_numeric(sweep_df['aligned_latency_error_ms'], errors='coerce').max()),
        'aligned_spikes_beyond_fit_end': int(pd.to_numeric(sweep_df['aligned_spikes_beyond_fit_end'], errors='coerce').fillna(0).sum()),
        'aligned_spikes_before_zero': int(pd.to_numeric(sweep_df['aligned_spikes_before_zero'], errors='coerce').fillna(0).sum()),
        's_bound_min': float(cfg['bounds']['s']['min']),
        's_bound_max': float(cfg['bounds']['s']['max']),
        'cells_with_s_below_v3_8_lower_bound': int(cell_df['s_below_v3_8_lower_bound'].astype(bool).sum()),
        'cells_near_new_s_lower_1pct': int(cell_df['s_near_lower_1pct'].astype(bool).sum()),
        'cells_near_new_s_lower_2pct': int(cell_df['s_near_lower_2pct'].astype(bool).sum()),
        'cells_near_new_s_lower_5pct': int(cell_df['s_near_lower_5pct'].astype(bool).sum()),
        'identifiability_status': 'PENDING_VISUAL_REVIEW',
    }
    (out / 'RUN_SUMMARY.json').write_text(json.dumps(_jsonable(summary), indent=2) + '\n')
    return summary


def _review_decision(row):
    rev = str(row.get('review', '')).strip().lower()
    if rev == 'accept': return 'ACCEPT'
    if rev == 'reject': return 'BAD'
    return str(row['post_v3_9_auto_cell_decision']).upper()


def _identify_one(payload):
    idx, cell, theta, cfg = payload
    ev = evaluate_theta(theta, cell, cfg, float(cfg['identifiability']['dt_ms']), float(cfg['loss']['vp_tau_ms']), 'identifiability')
    ident = practical_identifiability(theta, ev.loss, cell, cfg, idx)
    shifts = np.asarray([e.latency_shift_ms for e in ev.sweep_evals], dtype=float)
    return cell['cell_id'], ident, ev.loss, ev.spike_train_loss, ev.threshold_loss, shifts


def identify_final(cfg):
    cells, audit, baseline, thresholds = load_v3_6_cells(cfg)
    out = Path(cfg['output']['dir']); p = out / 'cell_fit_summary.csv'
    if not p.exists():
        raise FileNotFoundError('Run v3.9 fitting first: %s' % p)
    df = pd.read_csv(p)
    df['final_v3_9_decision'] = df.apply(_review_decision, axis=1)
    final = df[(df.final_v3_9_decision == 'ACCEPT') & (df.primary_support == 'MULTI_SWEEP')].copy()
    cells_by = {c['cell_id']: c for c in cells}
    theta_by = final.set_index('cell_id')[['b','r','s','kappa_I']].to_dict('index')
    payloads = [(i, cells_by[cid], theta_by[cid], cfg) for i, cid in enumerate(final.cell_id.astype(str))]
    outputs = []
    n_jobs = max(1, int(cfg['identifiability'].get('n_jobs', 1)))
    if n_jobs == 1:
        for x in payloads: outputs.append(_identify_one(x))
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            fut = [ex.submit(_identify_one, x) for x in payloads]
            for f in as_completed(fut): outputs.append(f.result())

    id_rows = []; alts = []
    for cid, ident, base_loss, spike_loss, threshold_loss, shifts in outputs:
        row = {
            'cell_id': cid, 'identifiability': ident.get('overall', 'UNKNOWN'),
            'identifiability_base_loss': float(base_loss),
            'identifiability_spike_train_loss': float(spike_loss),
            'identifiability_threshold_loss': float(threshold_loss),
            'identifiability_median_abs_latency_shift_ms': float(np.median(np.abs(shifts))) if shifts.size else np.nan,
            'identifiability_max_abs_latency_shift_ms': float(np.max(np.abs(shifts))) if shifts.size else np.nan,
            'latency_realigned_for_every_alternative': bool(ident.get('latency_realigned_for_every_alternative', False)),
        }
        for par in ('b','r','s','kappa_I'):
            ps = ident.get('parameter_status', {}).get(par, {})
            row['id_' + par] = ps.get('status', 'NA')
            row['id_alt_loss_' + par] = ps.get('best_separated_loss', np.nan)
        id_rows.append(row)
        for alt in ident.get('alternatives', []):
            alts.append({'cell_id': cid, **alt})
    id_df = pd.DataFrame(id_rows)
    id_df.to_csv(out / 'final_identifiability.csv', index=False)
    pd.DataFrame(alts).to_csv(out / 'final_identifiability_alternatives.csv', index=False)

    up = id_df.set_index('cell_id').to_dict('index') if not id_df.empty else {}
    for i, r in df.iterrows():
        cid = str(r['cell_id'])
        if cid in up:
            for k, v in up[cid].items(): df.loc[i, k] = v
        elif r['final_v3_9_decision'] == 'ACCEPT' and r['primary_support'] == 'SINGLE_SWEEP_ONLY':
            df.loc[i, 'identifiability'] = 'INSUFFICIENT_SPIKING_SWEEPS'
        elif r['final_v3_9_decision'] != 'ACCEPT':
            df.loc[i, 'identifiability'] = 'NOT_RUN'
    df.to_csv(p, index=False)

    primary = df[(df.final_v3_9_decision == 'ACCEPT') & (df.primary_support == 'MULTI_SWEEP')]
    gs = []
    for group, g in primary.groupby('group'):
        rr = {'group': group, 'n_cells': len(g), 'threshold_pass_fraction': float(g['threshold_pass'].mean())}
        for par in ('b','r','s','kappa_I','kappa_over_Cm'):
            rr['median_' + par] = g[par].median(); rr['q25_' + par] = g[par].quantile(.25); rr['q75_' + par] = g[par].quantile(.75)
        gs.append(rr)
    pd.DataFrame(gs).to_csv(out / 'group_parameter_summary_final_review.csv', index=False)

    result = {
        'n_identified_cells': len(id_df),
        'final_decisions': df['final_v3_9_decision'].value_counts().to_dict(),
        'threshold_pass_among_final_accept': df.loc[df.final_v3_9_decision == 'ACCEPT','threshold_pass'].value_counts().to_dict(),
        'identifiability': id_df['identifiability'].value_counts().to_dict() if not id_df.empty else {},
        'latency_realigned_for_every_alternative': bool(id_df['latency_realigned_for_every_alternative'].all()) if not id_df.empty else True,
    }
    (out / 'IDENTIFIABILITY_SUMMARY.json').write_text(json.dumps(_jsonable(result), indent=2) + '\n')
    return result
