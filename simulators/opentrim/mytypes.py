import os
import numpy as np
from numba import jit
from numba.core.types import float64, UniTuple

PROJ_DTYPE = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("ispec", np.int32),
    ("is_inside", np.bool_)
], align=True)

SIM_PARAMS_DTYPE = np.dtype([
    ("nspec", np.int32),
    ("nbin", np.int32),
    ("limits", np.float64, (2,)),
    ("emin", np.float64),
    ("ed", np.float64),
    ("pmax", np.float64),
    ("mean_free_path", np.float64),
    ("zmin", np.float64),
    ("zmax", np.float64),
    ("fac_linhard", np.float64, (2,)),
    ("density", np.float64),
    ("enorm", np.float64, (2,)),
    ("rnorm", np.float64, (2,)),
    ("dirfrac", np.float64, (2,)),
    ("denfrac", np.float64, (2,)),
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

class SimParams:
    # StatParams
    nspec: int
    nbin: int
    limits: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    
    # CascadeParams
    emin: float
    ed: float
    
    # RecoilParams
    pmax: float
    mean_free_path: float
    
    # GeometryParams
    zmin: float
    zmax: float
    
    # EstopParams
    fac_linhard: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    density: float
    
    # ScatterParams
    enorm: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    rnorm: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    dirfrac: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    denfrac: UniTuple(float64, 2)  # pyright: ignore[reportInvalidTypeForm]
    
    def __init__(self, stat_params_tup=None, cascade_params_tup=None, recoil_params_tup=None, geometry_params_tup=None, estop_params_tup=None, scatter_params_tup=None):
        self.nspec = 0
        self.nbin = 0
        self.limits = (0.0, 0.0)
        self.emin = 0.0
        self.ed = 0.0
        self.pmax = 0.0
        self.mean_free_path = 0.0
        self.zmin = 0.0
        self.zmax = 0.0
        self.fac_linhard = (0.0, 0.0)
        self.density = 0.0
        self.enorm = (0.0, 0.0)
        self.rnorm = (0.0, 0.0)
        self.dirfrac = (0.0, 0.0)
        self.denfrac = (0.0, 0.0)
        
        if stat_params_tup is not None:
            self.nspec, self.nbin, self.limits = stat_params_tup
        if cascade_params_tup is not None:
            self.emin, self.ed = cascade_params_tup 
        if recoil_params_tup is not None:
            self.pmax, self.mean_free_path = recoil_params_tup
        if geometry_params_tup is not None:
            self.zmin, self.zmax = geometry_params_tup     
        if estop_params_tup is not None:
            self.fac_linhard, self.density = estop_params_tup
        if scatter_params_tup is not None:
            self.enorm, self.rnorm, self.dirfrac, self.denfrac = scatter_params_tup

    def to_record(self):
        # TODO Vanilla numpy compatabilit
        rec = np.empty(1, dtype=SIM_PARAMS_DTYPE)[0]
        rec['nspec'] = self.nspec
        rec['nbin'] = self.nbin
        rec['limits'] = np.array(self.limits)
        rec['emin'] = self.emin
        rec['ed'] = self.ed
        rec['pmax'] = self.pmax
        rec['mean_free_path'] = self.mean_free_path
        rec['zmin'] = self.zmin
        rec['zmax'] = self.zmax
        rec['fac_linhard'] = self.fac_linhard
        rec['density'] = self.density
        rec['enorm'] = self.enorm
        rec['rnorm'] = self.rnorm
        rec['dirfrac'] = self.dirfrac
        rec['denfrac'] = self.denfrac
        return rec