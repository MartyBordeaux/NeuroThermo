from __future__ import annotations
import numpy as np

ISI_WT_EXIT=0.1358293233470019
ISI_BALANCE=0.5
ISI_SCA3_ENTRY=0.7978563373093712
ACTIVE_WT_EXIT=0.1483287106321978
ACTIVE_BALANCE=0.5
ACTIVE_SCA3_ENTRY=0.6832820264389761

TRANSFORMS = {
    'isi': {
        'columns':['log10_J_rheo','log10_mean_isi_ms'],
        'center': np.array([0.2392899006871207,1.2791714223626993],float),
        'scale': np.array([0.3160852097840443,0.1682174857544345],float),
        'wt_centroid': np.array([-0.2857035878224934,-0.4228321906245955],float),
        'sca3_centroid': np.array([2.2389510784632347,2.397345062344755],float),
        'thresholds':(ISI_WT_EXIT,ISI_BALANCE,ISI_SCA3_ENTRY),
    },
    'active_rate': {
        'columns':['log10_J_rheo','log10_active_support_rate_hz'],
        'center': np.array([0.2392899006871207,1.856048818911317],float),
        'scale': np.array([0.3160852097840443,0.1402900089308298],float),
        'wt_centroid': np.array([-0.2857035878224934,-0.3183184835867736],float),
        'sca3_centroid': np.array([2.2389510784632347,-1.7755275101606447],float),
        'thresholds':(ACTIVE_WT_EXIT,ACTIVE_BALANCE,ACTIVE_SCA3_ENTRY),
    }
}


def project(J_rheo, value, projection):
    if not (np.isfinite(J_rheo) and J_rheo>0 and np.isfinite(value) and value>0):
        return np.nan
    t=TRANSFORMS[projection]
    obs=np.array([np.log10(J_rheo),np.log10(value)],float)
    z=(obs-t['center'])/t['scale']
    w=t['sca3_centroid']-t['wt_centroid']
    return float(np.dot(z-t['wt_centroid'],w)/np.dot(w,w))
