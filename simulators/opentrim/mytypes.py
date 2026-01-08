import os
from numba.core.types import UniTuple
from numba.experimental import jitclass
import numpy as np
from numba import jit, float64

PROJ_DTYPE = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("ispec", np.int32),
    ("is_inside", np.bool_)
], align=True)

# Preserve compatibility with vanilla NumPy (with numba disabled)
if os.environ.get("NUMBA_DISABLE_JIT", "") == "1":
    def Projectile(e, pos, dir, ispec=0, is_inside=True):
        rec = np.recarray(1, dtype=PROJ_DTYPE)[0]
        rec['e'] = e
        rec['pos'] = pos    # copied
        rec['dir'] = dir    # copied
        rec['ispec'] = ispec
        rec['is_inside'] = is_inside
        return rec
else:
    @jit(inline = 'always')
    def Projectile(e, pos, dir, ispec=0, is_inside=True):
        rec = np.empty(1, dtype=PROJ_DTYPE)[0]
        rec['e'] = e
        rec['pos'] = pos    # copied
        rec['dir'] = dir    # copied
        rec['ispec'] = ispec
        rec['is_inside'] = is_inside
        return rec
    
@jitclass([("limits", UniTuple(float64, 2))])  # pyright: ignore[reportCallIssue]
class SimParams:
    nspec: int
    nbin: int
    limits: tuple[float, float]
    
    def __init__(self, nspec, nbin, limits):
        self.nspec = nspec
        self.nbin = nbin
        self.limits = limits
        
    def to_tuple(self):
        return self.nspec, self.nbin, self.limits