import math
import numpy as np
from transition_v1_3.decomposition import _interp_theta, _interp_J
from transition_v1_3.mathutils import persistent_crossing

SC={
 'wt_b':2.0,'sca_b':4.0,'wt_s':1.0,'sca_s':3.0,'wt_r':0.01,'sca_r':0.04,
 'wt_kappa_I':0.1,'sca_kappa_I':0.4,'wt_J_q75':1.0,'sca_J_q75':5.0,
 'wt_active_support_ms':900.0,'sca_active_support_ms':800.0,
}

def test_endpoint_interpolation():
    a=_interp_theta(SC,0,0); b=_interp_theta(SC,1,1)
    assert a['b']==2.0 and b['b']==4.0
    assert abs(a['r']-0.01)<1e-12 and abs(b['r']-0.04)<1e-12
    assert abs(a['kappa_I']-0.1)<1e-12 and abs(b['kappa_I']-0.4)<1e-12
    assert _interp_J(SC,0)==1.0 and _interp_J(SC,1)==5.0

def test_factorial_coupled_identity_parameters():
    for p in np.linspace(0,1,7):
        t_comb=_interp_theta(SC,p,p)
        t_k=_interp_theta(SC,p,p)
        assert all(abs(t_comb[k]-t_k[k])<1e-12 for k in t_comb)
        assert abs(_interp_J(SC,p)-_interp_J(SC,p))<1e-12

def test_persistent_crossing():
    x=np.array([0,.25,.5,.75,1.0]); y=np.array([0,.1,.4,.7,.9])
    z=persistent_crossing(x,y,.5,persistence=2)
    assert abs(z-(.5+(.5-.4)*(.75-.5)/(.7-.4)))<1e-12
