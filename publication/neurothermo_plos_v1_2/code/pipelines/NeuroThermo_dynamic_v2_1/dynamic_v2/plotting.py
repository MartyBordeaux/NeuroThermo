from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def make_plots(outdir):
    out=Path(outdir); figdir=out/'figures'; figdir.mkdir(exist_ok=True)
    group=pd.read_csv(out/'group_q_medians.csv')
    for metric in ['firing_rate_hz','active_rate_hz','mean_isi_ms','adaptation_index']:
        fig,ax=plt.subplots(figsize=(6,4))
        for (source,g),x in group.groupby(['source','group']):
            col=metric+'_median'
            if col not in x or x[col].notna().sum()==0: continue
            label='%s: %s'%(g,source)
            ax.plot(x.q,x[col],marker='o',label=label)
        ax.set_xlabel('q within observed suprathreshold support'); ax.set_ylabel(metric); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(figdir/('%s_vs_q.png'%metric),dpi=180); plt.close(fig)
    cov=pd.read_csv(out/'q_support_by_cell.csv')
    qcols=[c for c in cov.columns if c.startswith('q_') and c.endswith('_supported')]
    if qcols:
        counts=[int(cov[c].sum()) for c in qcols]
        labels=[c[len('q_'):-len('_supported')] for c in qcols]
        fig,ax=plt.subplots(figsize=(5,4)); ax.bar(labels,counts); ax.set_ylim(0,len(cov)); ax.set_xlabel('q'); ax.set_ylabel('cells with observed-current support')
        fig.tight_layout(); fig.savefig(figdir/'q_support_coverage.png',dpi=180); plt.close(fig)
    ph_path=out/'group_q_phase_median_profiles.csv'
    if ph_path.exists():
        ph=pd.read_csv(ph_path)
        if not ph.empty:
            for metric in ['x','y','z','speed','divergence']:
                for q in sorted(ph.q.unique()):
                    fig,ax=plt.subplots(figsize=(6,4)); z=ph[ph.q==q]
                    for group_name,g in z.groupby('group'):
                        ax.plot(g.phase,g[metric],label=group_name)
                    ax.set_xlabel('phase'); ax.set_ylabel(metric); ax.set_title('q=%g'%q); ax.legend()
                    fig.tight_layout(); fig.savefig(figdir/('phase_%s_q_%g.png'%(metric,q)),dpi=180); plt.close(fig)
