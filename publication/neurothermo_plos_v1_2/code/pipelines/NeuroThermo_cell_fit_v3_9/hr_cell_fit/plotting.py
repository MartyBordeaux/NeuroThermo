from __future__ import annotations
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def _cell_figure(cell, sweep_rows, threshold_row, cell_row, figsize):
    n=len(sweep_rows);raster_rows=max(1,math.ceil(n/2));fig=plt.figure(figsize=figsize);gs=fig.add_gridspec(raster_rows+1,2,height_ratios=[1.15]+[1.0]*raster_rows)
    currents=np.asarray([s['current_pA'] for s in cell['sweeps']],float);ne=np.asarray([r['n_exp_spikes'] for r in sweep_rows],float);nm=np.asarray([r['n_model_spikes'] for r in sweep_rows],float)
    le=np.asarray([r['exp_latency_ms'] for r in sweep_rows],float);lm=np.asarray([r['aligned_model_latency_ms'] for r in sweep_rows],float);lraw=np.asarray([r['raw_model_latency_ms'] for r in sweep_rows],float)
    i0=float(threshold_row['nonspiking_current_pA']);i1=float(threshold_row['first_spiking_current_pA']);n0=float(threshold_row['model_spikes_at_nonspiking_current']);n1=float(threshold_row['model_spikes_at_first_spiking_current'])
    ax=fig.add_subplot(gs[0,0]);ax.plot(currents,ne,'o-',label='experimental selected spikes');ax.plot(currents,nm,'o--',label='model');ax.axvspan(i0,i1,alpha=.08);ax.scatter([i0,i1],[n0,n1],marker='x',s=50);ax.set_xlabel('current (pA)');ax.set_ylabel('spike count');ax.legend(fontsize=6.5);ax.grid(alpha=.2)
    ax=fig.add_subplot(gs[0,1]);me=np.isfinite(le);ma=np.isfinite(lm);mr=np.isfinite(lraw)
    if me.any():ax.plot(currents[me],le[me],'o-',label='experimental')
    if mr.any():ax.plot(currents[mr],lraw[mr],'o:',alpha=.7,label='raw model')
    if ma.any():ax.plot(currents[ma],lm[ma],'o--',label='aligned model')
    ax.axvspan(i0,i1,alpha=.08);ax.set_xlabel('current (pA)');ax.set_ylabel('first spike (ms)');ax.legend(fontsize=7);ax.grid(alpha=.2)
    for k,(s,r) in enumerate(zip(cell['sweeps'],sweep_rows)):
        rr=1+k//2;cc=k%2;ax=fig.add_subplot(gs[rr,cc]);exp=np.asarray(s['exp_spike_times_ms'],float);aligned=np.asarray(r['model_spike_times_ms'],float);raw=np.asarray(r['raw_model_spike_times_ms'],float);fit_end=float(s['fit_end_ms']);ax.axvspan(0,fit_end,alpha=.08)
        if exp.size:ax.scatter(exp,np.full_like(exp,1.0),s=12)
        if aligned.size:ax.vlines(aligned,.25,.70,linewidth=.8)
        if raw.size:ax.vlines(raw,.02,.18,linewidth=.55,alpha=.55)
        ax.axvline(fit_end,linestyle='--',linewidth=.8);ax.set_ylim(0,1.2);ax.set_yticks([.10,.48,1.0]);ax.set_yticklabels(['raw','aligned','exp']);ax.set_title('%g pA | %s | tau=%+.1f ms | F1=%.2f' %(s['current_pA'],r['joint_sweep_decision'],r['latency_shift_ms'],r['circle_f1']),fontsize=8.2)
    if n%2==1:fig.add_subplot(gs[raster_rows,1]).axis('off')
    fig.tight_layout(rect=[0,0,1,.93]);return fig

def plot_cell(cell,sweep_rows,threshold_row,cell_row,out_path):
    n=len(sweep_rows);h=max(6.0,2.4*(max(1,math.ceil(n/2))+1));fig=_cell_figure(cell,sweep_rows,threshold_row,cell_row,(12,h));Path(out_path).parent.mkdir(parents=True,exist_ok=True);fig.savefig(out_path,dpi=150,bbox_inches='tight');plt.close(fig)

def make_audit_pdf(cells_by_id,sweep_rows_by_cell,threshold_rows_by_cell,cell_rows,out_path):
    Path(out_path).parent.mkdir(parents=True,exist_ok=True)
    with PdfPages(out_path) as pdf:
        for row in cell_rows:
            cid=row['cell_id'];fig=_cell_figure(cells_by_id[cid],sweep_rows_by_cell[cid],threshold_rows_by_cell[cid],row,(11.7,8));pdf.savefig(fig,dpi=110);plt.close(fig)
