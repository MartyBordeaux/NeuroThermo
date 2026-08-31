from __future__ import annotations
import math
import numpy as np
PARAM_NAMES=('b','r','s','kappa_I')

def _scale(z,lo,hi,scale):
    z=float(np.clip(z,0,1))
    return math.exp(math.log(lo)+z*(math.log(hi)-math.log(lo))) if scale=='log' else lo+z*(hi-lo)

def _inv(v,lo,hi,scale):
    return (math.log(v)-math.log(lo))/(math.log(hi)-math.log(lo)) if scale=='log' else (v-lo)/(hi-lo)

def z_to_params(z,bounds):return {n:_scale(z[i],bounds[n]['min'],bounds[n]['max'],bounds[n]['scale']) for i,n in enumerate(PARAM_NAMES)}
def params_to_z(params,bounds):return np.asarray([_inv(params[n],bounds[n]['min'],bounds[n]['max'],bounds[n]['scale']) for n in PARAM_NAMES],float)
def pack_params(params):return np.asarray([params[n] for n in PARAM_NAMES],float)
