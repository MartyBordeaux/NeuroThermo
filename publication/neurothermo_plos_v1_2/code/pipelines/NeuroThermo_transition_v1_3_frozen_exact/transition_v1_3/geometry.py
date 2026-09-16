from __future__ import annotations
import numpy as np
import pandas as pd

PRIMARY = 'isi_primary_v1_0_frozen'
SECONDARY = 'active_rate_experimental_v2_1'


def _projection(refs, transforms, name):
    r = refs[refs.projection.eq(name)]
    if len(r) != 1:
        raise ValueError(f'projection reference {name!r} not unique')
    r = r.iloc[0]
    t = transforms[transforms.projection.eq(name)].set_index('coordinate')
    if name == PRIMARY:
        coords = ['log10_rheobase', 'log10_isi']
    elif name == SECONDARY:
        coords = ['log10_rheobase', 'log10_active_rate']
    else:
        raise ValueError(name)
    center = np.array([float(t.loc[c, 'center']) for c in coords])
    scale = np.array([float(t.loc[c, 'scale']) for c in coords])
    cwt = np.array([float(r.wt_centroid_0), float(r.wt_centroid_1)])
    csc = np.array([float(r.sca3_centroid_0), float(r.sca3_centroid_1)])
    delta = csc - cwt
    den = float(delta @ delta)
    return {
        'name': name, 'coords': coords, 'center': center, 'scale': scale,
        'cwt': cwt, 'csc': csc, 'delta': delta, 'den': den,
        'wt_exit': float(r.wt_exit_A_threshold),
        'sca3_entry': float(r.sca3_entry_A_threshold),
        'corridor_radius_q90': float(r.corridor_radius_q90),
        'cloud_overlap': bool(r.cloud_overlap),
    }


def load_geometry(refs, transforms):
    return {
        'isi': _projection(refs, transforms, PRIMARY),
        'active': _projection(refs, transforms, SECONDARY),
    }


def project_scalar(log_rheobase, log_second, ref):
    if not (np.isfinite(log_rheobase) and np.isfinite(log_second)):
        return np.nan, np.nan
    arr = np.array([float(log_rheobase), float(log_second)])
    z = (arr - ref['center']) / ref['scale']
    A = float(((z - ref['cwt']) @ ref['delta']) / ref['den'])
    foot = ref['cwt'] + A * ref['delta']
    orth = float(np.linalg.norm(z - foot))
    return A, orth


def stage_from_A(A, ref):
    if not np.isfinite(A):
        return 'INVALID'
    if A <= ref['wt_exit']:
        return 'WT_like'
    if A >= ref['sca3_entry']:
        return 'SCA3_like'
    return 'TRANSITION'


def weighted_quantile(values, weights, q):
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[m], w[m]
    if len(v) == 0:
        return np.nan
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    c = c / c[-1]
    return float(np.interp(q, c, v))


def persistent_crossing(x, y, thr, persistence=2):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    for i in range(n):
        if not np.isfinite(y[i]) or y[i] <= thr:
            continue
        j = min(n, i + int(persistence))
        if np.all(np.isfinite(y[i:j])) and np.all(y[i:j] > thr):
            if i == 0:
                return float(x[0])
            if np.isfinite(y[i-1]) and y[i-1] <= thr:
                if y[i] == y[i-1]:
                    return float(x[i])
                return float(x[i-1] + (thr-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]))
            return float(x[i])
    return np.nan
