#!/usr/bin/env python3
"""Verify manuscript numerators, medians and paired protocol comparisons.

This script does not fit models or alter source tables. It writes an audit JSON.
"""
from pathlib import Path
import json,hashlib,sys,platform,re
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1];D=R/'source_data';P=R/'provenance';P.mkdir(exist_ok=True)
a=pd.read_csv(D/'v130/Q_PAIR_SUMMARY.csv')
a=a[(a.tier=='full')&(a.observable=='isi_dist')&(a.q>=.75)]
expected={'SCENARIO_MARKER_FIRST':[(26,29),(28,30),(28,31),(29,32),(30,32),(30,32)],'AFFINITY_FIELD_FIRST':[(26,29),(27,31),(27,31),(28,32),(28,32),(28,32)]}
rows=[]
for met,counts in expected.items():
 g=a[a.method==met]
 for i,(q,z) in enumerate(g.groupby('q')):
  d=z.delta.dropna();got=(int((d<0).sum()),len(d));assert got==counts[i],(met,q,got)
  rows.append({'method':met,'q':float(q),'attempted_pairs':len(z),'evaluable_pairs':len(d),'negative_pairs':got[0],'median_delta':float(d.median()),'q1':float(d.quantile(.25)),'q3':float(d.quantile(.75))})
 pivot=g.pivot(index='pair',columns='q',values='delta').dropna()
 assert len(pivot)==29
 assert int((pivot<0).all(axis=1).sum())==(26 if met=='SCENARIO_MARKER_FIRST' else 25)
pd.DataFrame(rows).to_csv(D/'MANUSCRIPT_Q_SUMMARY.csv',index=False)
null=pd.read_csv(D/'v140/Q_CROSS_VS_NULL_COHORT.csv')
null=null[(null.tier=='full')&(null.observable=='isi_dist')&(null.method=='SCENARIO_MARKER_FIRST')]
expected_null={'WT-WT':[(10,26),(7,25),(8,26),(9,27),(10,27),(9,27)],'SCA3-SCA3':[(1,3),(1,3),(1,3),(0,3),(0,3),(0,4)]}
for gr,ex in expected_null.items():
 got=[(int(x.pairs_negative),int(x.pairs_evaluable)) for _,x in null[null.group==gr].sort_values('q').iterrows()];assert got==ex,(gr,got)
null[['group','q','pairs_total','pairs_evaluable','pairs_negative','median_delta','median_abs_delta']].to_csv(D/'MANUSCRIPT_NULL_SUMMARY.csv',index=False)
f=pd.read_csv(D/'v120/FAMILY_SUMMARY.csv');f=f[(f.tier=='fixed_full')&(f.observable=='isi_dist')&(f.method=='SCENARIO_MARKER_FIRST')]
q=a[(a.q==.75)&(a.method=='SCENARIO_MARKER_FIRST')]
j=q[['pair','delta']].merge(f[['family','delta']],left_on='pair',right_on='family',suffixes=('_q75','_fixed')).dropna()
assert len(j)==23 and (j.delta_q75<0).sum()==21 and (j.delta_fixed<0).sum()==8
npchg=int(((j.delta_q75<0)&(j.delta_fixed>0)).sum());pnchg=int(((j.delta_q75>0)&(j.delta_fixed<0)).sum());assert (npchg,pnchg)==(13,0)
j.to_csv(D/'fixed_vs_q75_matched_pairs.csv',index=False)
fixed={'n_pairs':len(j),'q75_negative':int((j.delta_q75<0).sum()),'fixed_negative':int((j.delta_fixed<0).sum()),'q75_median':float(j.delta_q75.median()),'fixed_median':float(j.delta_fixed.median()),'negative_to_positive':npchg,'positive_to_negative':pnchg}
sup=pd.read_csv(D/'frozen_core_support.csv');assert len(sup)==34 and sup.cell_id.nunique()==12
st=pd.read_csv(D/'v130/Q_SCENARIO_MARKERS.csv');ss=st[(st.tier=='full')&(st.q_key=='q075')&(st.observable=='isi_dist')];assert len(ss)==264 and ss.delta.notna().sum()==233
review=pd.read_csv(P/'REVIEWER_POINT_MAP.csv');assert len(review)==23
tex='\n'.join((R/n).read_text() for n in ['manuscript_body.tex','supplementary.tex','response_to_reviewers.tex'])
keys=set(x.strip() for b in re.findall(r'\\cite\w*\{([^}]*)\}',tex) for x in b.split(','))
bib=set(re.findall(r'@\w+\{\s*([^,]+),',(R/'references.bib').read_text()))
assert keys<=bib,sorted(keys-bib)
report={'status':'PASS','q_summary':rows,'fixed_current_paired':fixed,'retained_sets':34,'endpoint_cells':12,'q075_scenarios':264,'q075_finite_markers':233,'reviewer_points':23,'citation_keys_checked':len(keys),'interpretation':'Recomputed descriptive model summaries; not a new biological test.'}
(P/'NUMERICAL_RECHECK.json').write_text(json.dumps(report,indent=2)+'\n')
(P/'RENDER_ENVIRONMENT.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scope':'Figure rendering and summary verification; original scientific packages target Python 3.9.'},indent=2)+'\n')
print(json.dumps({'status':'PASS','reviewer_points':23,'q075_finite':233,'fixed_current':fixed},indent=2))
