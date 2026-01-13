import os
from typing import Optional
import numpy as np
from numba import jit
from numba.experimental import jitclass
from numba.core.types import float64, UniTuple

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

@jitclass
class ScatterParams:
    enorm: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    rnorm: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    dirfrac: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    denfrac: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    
    def __init__(self, enorm, rnorm, dirfrac, denfrac):
        self.enorm = (enorm[0], enorm[1])
        self.rnorm = (rnorm[0], rnorm[1])
        self.dirfrac = (dirfrac[0], dirfrac[1])
        self.denfrac = (denfrac[0], denfrac[1])
        
    def to_tuple(self):
        return self.enorm, self.rnorm, self.dirfrac, self.denfrac

@jitclass
class RecoilParams:
    pmax: float
    mean_free_path: float
    
    def __init__(self, pmax, mean_free_path):
        self.pmax = pmax
        self.mean_free_path = mean_free_path
        
    def to_tuple(self):
        return self.pmax, self.mean_free_path

@jitclass
class EstopParams:
    fac_linhard: UniTuple(float64, 2)    # pyright: ignore[reportInvalidTypeForm]
    density: float
    
    def __init__(self, fac_linhard, density):
        self.fac_linhard = fac_linhard
        self.density = density
        
    def to_tuple(self):
        return self.fac_linhard, self.density

@jitclass
class GeometryParams:
    zmin: float
    zmax: float
    
    def __init__(self, zmin, zmax):
        self.zmin = zmin
        self.zmax = zmax
        
    def to_tuple(self):
        return self.zmin, self.zmax

@jitclass
class CascadeParams:
    emin: float
    ed: float
    
    def __init__(self, emin, ed):
        self.emin = emin
        self.ed = ed
        
    def to_tuple(self):
        return self.emin, self.ed

@jitclass
class StatParams:
    nspec: int
    nbin: int
    limits: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    
    def __init__(self, nspec, nbin, limits):
        self.nspec = nspec
        self.nbin = nbin
        self.limits = limits
        
    def to_tuple(self):
        return self.nspec, self.nbin, self.limits

@jitclass
class SimParams:
    stat_params: Optional[StatParams]
    cascade_params: Optional[CascadeParams]
    recoil_params: Optional[RecoilParams]
    geometry_params: Optional[GeometryParams]
    estop_params: Optional[EstopParams]
    scatter_params: Optional[ScatterParams]
    
    def __init__(self, stat_params_tup=None, cascade_params_tup=None, recoil_params_tup=None, geometry_params_tup=None, estop_params_tup=None, scatter_params_tup=None):
        self.stat_params = None
        self.cascade_params = None
        self.recoil_params = None
        self.geometry_params = None
        self.estop_params = None
        self.scatter_params = None
        if stat_params_tup is not None:
            self.stat_params = StatParams(*stat_params_tup)
        if cascade_params_tup is not None:
            self.cascade_params = CascadeParams(*cascade_params_tup)
        if recoil_params_tup is not None:
            self.recoil_params = RecoilParams(*recoil_params_tup)
        if geometry_params_tup is not None:
            self.geometry_params = GeometryParams(*geometry_params_tup)
        if estop_params_tup is not None:
            self.estop_params = EstopParams(*estop_params_tup)
        if scatter_params_tup is not None:
            self.scatter_params = ScatterParams(*scatter_params_tup)
        
    def to_tuple(self):
        return self.stat_params.to_tuple() if self.stat_params is not None else None, \
               self.cascade_params.to_tuple() if self.cascade_params is not None else None, \
               self.recoil_params.to_tuple() if self.recoil_params is not None else None,   \
               self.geometry_params.to_tuple() if self.geometry_params is not None else None, \
               self.estop_params.to_tuple() if self.estop_params is not None else None, \
               self.scatter_params.to_tuple() if self.scatter_params is not None else None