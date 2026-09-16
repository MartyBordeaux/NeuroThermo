from __future__ import annotations
import numpy as np

ACTIVE_PARAMS = ('b', 'r', 's', 'kappa_I')


def active_params(cfg):
    return ACTIVE_PARAMS


def unit_to_theta(u, cfg):
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    if len(u) != len(ACTIVE_PARAMS):
        raise ValueError('Expected %d active parameters, got %d' % (len(ACTIVE_PARAMS), len(u)))
    theta = {}
    for i, name in enumerate(ACTIVE_PARAMS):
        spec = cfg['bounds'][name]
        lo, hi = float(spec['min']), float(spec['max'])
        if spec.get('scale', 'linear') == 'log':
            if lo <= 0 or hi <= 0:
                raise ValueError('Log-scaled bound %s must be positive' % name)
            theta[name] = float(np.exp(np.log(lo) + u[i] * (np.log(hi) - np.log(lo))))
        else:
            theta[name] = float(lo + u[i] * (hi - lo))
    return theta


def theta_to_unit(theta, cfg):
    out = np.empty(len(ACTIVE_PARAMS), dtype=float)
    for i, name in enumerate(ACTIVE_PARAMS):
        spec = cfg['bounds'][name]
        lo, hi = float(spec['min']), float(spec['max'])
        x = float(theta[name])
        if spec.get('scale', 'linear') == 'log':
            out[i] = (np.log(x) - np.log(lo)) / (np.log(hi) - np.log(lo))
        else:
            out[i] = (x - lo) / (hi - lo)
    return np.clip(out, 0.0, 1.0)
