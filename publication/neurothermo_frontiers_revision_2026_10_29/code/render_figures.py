#!/usr/bin/env python3
"""Rebuild revision figures from archived data; no fit or cohort reselection."""
from pathlib import Path
import json,re,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1]; D=ROOT/'source_data';O=ROOT/'figures';O.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'code'))
from nt_pm.hr import simulate,interpolate
from nt_pm.active import active_observables
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':11,'legend.fontsize':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
def start(size=(5.1,3.5)):return plt.subplots(figsize=size,layout='constrained')
def save(fig,name):
 m=re.match(r'Fig[1-4]([A-D])_',name)
 if m:fig.text(.006,.99,m.group(1),va='top',fontweight='bold',fontsize=13)
 fig.savefig(O/name,bbox_inches='tight');plt.close(fig)
for group,lab,fn in [('WT_14','A','Fig1A_native_WT.png'),('SCA3_04','B','Fig1B_native_SCA3.png')]:
 im=Image.open(D/(group+'_native_400pA_audit.png'));crop=im.crop((0,0,im.width,round(im.height*0.263)));cv=Image.new('RGB',(im.width,crop.height+24),'white');cv.paste(crop,(0,24));ImageDraw.Draw(cv).text((8,2),lab,fill='black',font_size=20);cv.save(O/fn)
e=pd.read_csv(D/'experimental_endpoints_archival.csv')
f,a=start((5,3.3))
for i,g in enumerate(['WT','SCA3']):
 z=e[e.group==g].sort_values('cell_id');a.scatter(i+np.linspace(-.13,.13,len(z)),z.exp_isi,s=33,marker='o' if i==0 else '^');a.plot([i-.22,i+.22],[z.exp_isi.median()]*2,lw=2)
a.set(xticks=[0,1],xticklabels=['WT (12)','SCA3 (6)'],ylabel='Experimental mean ISI (ms)',title='Archival endpoints at q = 0.75',xlim=(-.5,1.5));save(f,'Fig1C_ISI_endpoints.pdf')
f,a=start((5,3.3))
for g,m in [('WT','o'),('SCA3','^')]:
 z=e[e.group==g];a.scatter(z.exp_isi,z.model_isi,s=35,marker=m,label=g)
mx=max(e.exp_isi.max(),e.model_isi.max())*1.08;a.plot([0,mx],[0,mx],ls='--',label='Identity');a.set(xlabel='Experimental mean ISI (ms)',ylabel='Model mean ISI (ms)',xlim=(0,mx),ylim=(0,mx),title='Endpoint reproduction');a.legend(frameon=False);save(f,'Fig1D_ISI_reproduction.pdf')
cfg=json.load(open(D/'illustration_config.json'));pars=pd.read_csv(D/'frozen_core_support.csv');dr=pd.read_csv(D/'frozen_drives_from_v1_0_1.csv').set_index('cell_id')
pa=pars[(pars.cell_id=='WT_14')&(pars.solution_key=='best')].iloc[0].to_dict();pb=pars[(pars.cell_id=='SCA3_04')&(pars.solution_key=='best')].iloc[0].to_dict();ref=json.load(open(D/'v110/GLOBAL_STATE_SCALE.json'));mu=np.asarray(ref['center_xyz']);sc=np.asarray(ref['scale_xyz']);replays=[]
for p in [0.,.5,1.]:
 par,J=interpolate(pa,pb,p,dr.loc['WT_14','q75_J_native'],dr.loc['SCA3_04','q75_J_native']);sim=simulate(par,J,cfg);ob=active_observables(sim,cfg)
 if not sim['ok'] or ob['status']!='PASS':raise RuntimeError('Illustration replay failed')
 replays.append((p,sim,ob));pd.DataFrame(ob['state_active'],columns=['x','y','z']).to_csv(D/f'illustration_active_state_p{p:.2f}.csv',index=False);pd.DataFrame({'isi_ms':ob['isi_ms']}).to_csv(D/f'illustration_isi_p{p:.2f}.csv',index=False)
f,a=start((5.2,3.7))
for i,(p,s,o) in enumerate(replays):
 t=np.arange(1,len(s['states'])+1)*.5;k=t<=250;a.plot(t[k],s['states'][k,0]+5*i,lw=.8,label=f'p = {p:.2f}')
a.set(xlabel='Time after current onset (ms)',ylabel='HR x (offset by 5 units per trace)',title='WT_14 to SCA3_04, q = 0.75');a.legend(frameon=False);save(f,'Fig2A_model_traces.pdf')
f,a=start((5.2,3.7))
for p,s,o in replays:
 z=(o['state_active']-mu)/sc;a.scatter(z[:,0],z[:,2],s=2,alpha=.45,label=f'p = {p:.2f}',rasterized=True)
a.set(xlabel='Standardized x',ylabel='Standardized z',title='Active-cycle samples');a.legend(markerscale=4,frameon=False);save(f,'Fig2B_state_projection.pdf')
f,a=start((5.2,3.7))
for p,s,o in replays:
 x=np.sort(o['log_isi']);a.step(x,np.arange(1,len(x)+1)/len(x),where='post',label=f'p = {p:.2f}')
a.set(xlabel='log(ISI / ms)',ylabel='Empirical cumulative probability',ylim=(0,1.02),title='Same accepted ISIs');a.legend(frameon=False);save(f,'Fig2C_ISI_cdf.pdf')
f,a=start((5.2,3.7));no=pd.read_csv(D/'illustrative_pair_nodes.csv').sort_values('p')
for col,lab in [('state_affinity','Active state'),('isi_dist_affinity','log(ISI)')]:a.plot(no.p,no[col],marker='o',ms=3,label=lab)
a.axhline(0,ls=':',lw=.9);a.set(xlabel='Interpolation coordinate p',ylabel='Endpoint affinity A(p)',ylim=(-1.12,1.12),title='Archived affinity curves');a.legend(frameon=False);save(f,'Fig2D_affinities.pdf')
qp=pd.read_csv(D/'v130/Q_PAIR_SUMMARY.csv');qp=qp[(qp.tier=='full')&(qp.observable=='isi_dist')&(qp.q>=.75)];xs=sorted(qp.q.unique())
f,a=start((5.2,3.6))
for met,lab in [('SCENARIO_MARKER_FIRST','Scenario-first'),('AFFINITY_FIELD_FIRST','Field-first')]:
 g=qp[qp.method==met].groupby('q').delta;a.plot(xs,g.median().reindex(xs),marker='o',label=lab);a.fill_between(xs,g.quantile(.25).reindex(xs),g.quantile(.75).reindex(xs),alpha=.14)
a.axhline(0,ls=':',lw=.9);a.set(xlabel='Relative current q',ylabel='Pair median Delta',xticks=xs,title='Median and pair IQR');a.legend(frameon=False);save(f,'Fig3A_delta_q.pdf')
piv=qp[qp.method=='SCENARIO_MARKER_FIRST'].pivot(index='pair',columns='q',values='delta');complete=piv.dropna();assert len(complete)==29
f,a=start((5.2,3.6))
for _,row in complete.iterrows():a.plot(complete.columns,row,alpha=.4,lw=.8)
a.plot(complete.columns,complete.median(),lw=2.8,label='Median of 29 complete pairs');a.axhline(0,ls=':',lw=.9);a.set(xlabel='Relative current q',ylabel='Pair Delta',xticks=xs,title='Same 29 pairs at all six q levels');a.legend(frameon=False);save(f,'Fig3B_complete_pairs.pdf')
c=pd.read_csv(D/'v140/Q_CROSS_VS_NULL_COHORT.csv');c=c[(c.tier=='full')&(c.method=='SCENARIO_MARKER_FIRST')&(c.observable=='isi_dist')];f,a=start((5.1,3.6))
for gr in ['WT-SCA3','WT-WT','SCA3-SCA3']:
 z=c[c.group==gr].sort_values('q');a.plot(z.q,z.fraction_negative_among_evaluable,marker='o',label=gr)
a.set(xlabel='Relative current q',ylabel='Fraction with Delta < 0',ylim=(-.03,1.03),xticks=xs,title='Oriented within-genotype controls');a.legend(frameon=False);save(f,'Fig4A_oriented_null.pdf')
ep=pd.read_csv(D/'v140/CROSS_AND_NULL_ENDPOINT_CONTRASTS.csv');f,a=start((5.1,3.6))
for gr,m in [('WT-SCA3','o'),('WT-WT','s'),('SCA3-SCA3','^')]:
 z=ep[(ep.q==.75)&(ep.group==gr)];a.scatter(z.state_endpoint_w1_weighted_median,z.isi_endpoint_w1_weighted_median,s=25,alpha=.8,marker=m,label=gr)
a.set(xlabel='Endpoint active-state distance',ylabel='Endpoint log(ISI) distance',title='Endpoint geometry at q = 0.75');a.legend(frameon=False);save(f,'Fig4B_endpoint_contrasts.pdf')
ma=pd.read_csv(D/'v140/Q_MATCHED_CONTROL_SUMMARY.csv');ma=ma[ma.method=='SCENARIO_MARKER_FIRST'].sort_values('q');f,a=start((5.1,3.6));a.plot(ma.q,ma.n_distinct_matched_cross_pairs,marker='o',label='Cross pairs with a control');a.plot(ma.q,ma.n_distinct_matched_null_pairs,marker='s',label='Unique within-genotype pairs');a.set(xlabel='Relative current q',ylabel='Number of distinct pairs',ylim=(0,7),xticks=xs,title='Both endpoint distances matched within 25%');a.legend(frameon=False);save(f,'Fig4C_matching_coverage.pdf')
fa=pd.read_csv(D/'v120/FAMILY_SUMMARY.csv');fa=fa[(fa.tier=='fixed_full')&(fa.observable=='isi_dist')&(fa.method=='SCENARIO_MARKER_FIRST')];qa=qp[(qp.q==.75)&(qp.method=='SCENARIO_MARKER_FIRST')];jo=qa[['pair','delta']].merge(fa[['family','delta']],left_on='pair',right_on='family',suffixes=('_q75','_fixed')).dropna();assert len(jo)==23;jo.to_csv(D/'fixed_vs_q75_matched_pairs.csv',index=False);f,a=start((5.1,3.6));a.scatter(jo.delta_q75,jo.delta_fixed,s=30);a.axhline(0,ls=':',lw=.9);a.axvline(0,ls=':',lw=.9);a.set(xlabel='Delta: q = 0.75',ylabel='Delta: fixed physical current',title='Same 23 evaluable pairs');save(f,'Fig4D_fixed_current.pdf')
f,a=start((7.1,8.1));piv=piv.sort_index();im=a.imshow(np.ma.masked_invalid(piv.to_numpy()),aspect='auto');a.set(xticks=np.arange(6),xticklabels=[f'{q:.2f}' for q in xs],yticks=np.arange(len(piv)),yticklabels=[p.replace('__',' / ') for p in piv.index],xlabel='Relative current q',title='All 32 pairs; missing values remain blank');a.tick_params(axis='y',labelsize=7);f.colorbar(im,ax=a,label='Delta');save(f,'FigS1_pair_heatmap.pdf')
f,a=start((6.2,4.4));se=pd.read_csv(D/'v110/SENSITIVITY_PAIR_AUDIT.csv')
for v in ['gap_strict','gap_loose','fine_dt','long_pre']:
 z=se[(se.method=='SCENARIO_MARKER_FIRST')&(se.variant==v)&se.both_defined];a.scatter(z.baseline_delta,z.variant_delta,s=22,alpha=.7,label=f'{v} (n={len(z)})')
a.plot([-.7,.7],[-.7,.7],ls='--',lw=.8);a.axhline(0,ls=':',lw=.7);a.axvline(0,ls=':',lw=.7);a.set(xlabel='Baseline pair Delta',ylabel='Sensitivity pair Delta',title='q = 0.75; support differs between checks');a.legend(frameon=False);save(f,'FigS2_sensitivity.pdf')
print('Rendered',len(list(O.glob('*'))),'panels')
