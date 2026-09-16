from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def _cell_figure(cell, sweep_rows, threshold_row, cell_row, figsize):
    n = len(sweep_rows)
    raster_rows = max(1, math.ceil(n / 2))
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(raster_rows + 1, 2, height_ratios=[1.15] + [1.0] * raster_rows)
    currents = np.asarray([s['current_pA'] for s in cell['sweeps']], dtype=float)
    ne = np.asarray([r['n_exp_spikes'] for r in sweep_rows], dtype=float)
    nm = np.asarray([r['n_model_spikes'] for r in sweep_rows], dtype=float)
    nraw = np.asarray([r['raw_n_model_spikes'] for r in sweep_rows], dtype=float)
    le = np.asarray([r['exp_latency_ms'] for r in sweep_rows], dtype=float)
    lm = np.asarray([r['aligned_model_latency_ms'] for r in sweep_rows], dtype=float)
    lraw = np.asarray([r['raw_model_latency_ms'] for r in sweep_rows], dtype=float)
    i0 = float(threshold_row['nonspiking_current_pA'])
    i1 = float(threshold_row['first_spiking_current_pA'])
    n0 = float(threshold_row['model_spikes_at_nonspiking_current'])
    n1 = float(threshold_row['model_spikes_at_first_spiking_current'])

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(currents, ne, 'o-', label='experimental selected spikes')
    ax.plot(currents, nm, 'o--', label='model (count preserved by alignment)')
    ax.axvspan(i0, i1, alpha=0.08, label='experimental rheobase bracket')
    ax.scatter([i0, i1], [n0, n1], marker='x', s=50, label='model: threshold probes')
    ax.set_xlabel('current (pA)')
    ax.set_ylabel('spike count')
    ax.set_title('Spike trains + binary threshold bracket')
    ax.legend(fontsize=6.5)
    ax.grid(alpha=0.2)

    ax = fig.add_subplot(gs[0, 1])
    me = np.isfinite(le); ma = np.isfinite(lm); mr = np.isfinite(lraw)
    if me.any(): ax.plot(currents[me], le[me], 'o-', label='experimental')
    if mr.any(): ax.plot(currents[mr], lraw[mr], 'o:', alpha=0.7, label='raw model')
    if ma.any(): ax.plot(currents[ma], lm[ma], 'o--', label='first-spike-aligned model')
    ax.axvspan(i0, i1, alpha=0.08)
    ax.set_xlabel('current (pA)')
    ax.set_ylabel('first spike in comparison window (ms)')
    ax.set_title('Exact first-spike alignment; no latency penalty')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2)

    for k, (s, r) in enumerate(zip(cell['sweeps'], sweep_rows)):
        rr = 1 + k // 2; cc = k % 2
        ax = fig.add_subplot(gs[rr, cc])
        exp = np.asarray(s['exp_spike_times_ms'], dtype=float)
        aligned = np.asarray(r['model_spike_times_ms'], dtype=float)
        raw = np.asarray(r['raw_model_spike_times_ms'], dtype=float)
        fit_end = float(s['fit_end_ms'])
        ax.axvspan(0, fit_end, alpha=0.08)
        if exp.size: ax.scatter(exp, np.full_like(exp, 1.0), s=12, label='exp')
        if aligned.size: ax.vlines(aligned, 0.25, 0.70, linewidth=0.8, label='first-spike-aligned model')
        if raw.size: ax.vlines(raw, 0.02, 0.18, linewidth=0.55, alpha=0.55, label='raw model')
        ax.axvline(fit_end, linestyle='--', linewidth=0.8)
        xmax = max(80.0, fit_end + 10.0, float(np.max(aligned)) + 10.0 if aligned.size else 0.0)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([0.10, 0.48, 1.0])
        ax.set_yticklabels(['raw', 'aligned', 'exp'])
        ax.set_title('%g pA | %s | tau=%+.1f ms | F1=%.2f | VP %.2f -> %.2f' % (
            s['current_pA'], r['joint_sweep_decision'], r['latency_shift_ms'], r['circle_f1'],
            r['raw_vp_loss'], r['vp_loss']), fontsize=8.2)
        ax.tick_params(labelsize=7.3)
    if n % 2 == 1:
        fig.add_subplot(gs[raster_rows, 1]).axis('off')

    fig.suptitle(
        '%s (%s) | auto=%s | threshold=%s [%g, %g] pA; model counts=(%d, %d)\n'
        'b=%.5g, r=%.5g, s=%.5g, kappa_I=%.5g | first-spike-aligned loss=%.4f, threshold loss=%.2f, total=%.4f | median |tau|=%.1f ms' % (
            cell['cell_id'], cell['group'], cell_row['post_v3_9_auto_cell_decision'],
            'PASS' if bool(cell_row['threshold_pass']) else 'FAIL', i0, i1, int(n0), int(n1),
            cell_row['b'], cell_row['r'], cell_row['s'], cell_row['kappa_I'],
            cell_row['spike_train_loss'], cell_row['threshold_loss'], cell_row['cell_loss'],
            cell_row['median_abs_latency_shift_ms']
        ), fontsize=10.2,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


def plot_cell(cell, sweep_rows, threshold_row, cell_row, out_path):
    n = len(sweep_rows); raster_rows = max(1, math.ceil(n / 2)); h = max(6.0, 2.4 * (raster_rows + 1))
    fig = _cell_figure(cell, sweep_rows, threshold_row, cell_row, (12, h))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def make_audit_pdf(cells_by_id, sweep_rows_by_cell, threshold_rows_by_cell, cell_rows, out_path):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out_path) as pdf:
        for row in cell_rows:
            cid = row['cell_id']; cell = cells_by_id[cid]; sr = sweep_rows_by_cell[cid]; tr = threshold_rows_by_cell[cid]
            n = len(sr); raster_rows = max(1, math.ceil(n / 2)); h = max(7.0, 2.35 * (raster_rows + 1))
            fig = _cell_figure(cell, sr, tr, row, (11.7, h))
            pdf.savefig(fig, dpi=110)
            plt.close(fig)
