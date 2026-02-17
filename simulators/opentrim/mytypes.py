import os
import numpy as np
from numba import jit

PROJ_DTYPE = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("ispec", np.int32),
    ("is_inside", np.bool_)
], align=True)

NLHLIN_COEFS_DTYPE = np.dtype([
    ("z1", np.uint32),
    ("z2", np.uint32),
    ("a1", np.float64),
    ("b1", np.float64),
    ("a2", np.float64),
    ("b2", np.float64),
    ("a3", np.float64),
    ("b3", np.float64),
    ("rmax", np.float64),
], align=True)

STAT_PARAMS_DTYPE = np.dtype([
    ("nspec", np.int32),
    ("nbin", np.int32),
    ("limits", np.float64, (2,)),
], align=True)

CASCADE_PARAMS_DTYPE = np.dtype([
    ("emin", np.float64),
    ("ed", np.float64),
], align=True)

RECOIL_PARAMS_DTYPE = np.dtype([
    ("pmax", np.float64),
    ("mean_free_path", np.float64),
], align=True)

GEOMETRY_PARAMS_DTYPE = np.dtype([
    ("zmin", np.float64),
    ("zmax", np.float64),
], align=True)

ESTOP_PARAMS_DTYPE = np.dtype([
    ("fac_linhard", np.float64, (2,)),
    ("density", np.float64),
], align=True)

SCATTER_PARAMS_DTYPE = np.dtype([
    ("pot_model", "<U16"),
    ("z1", np.uint32),
    ("z2", np.uint32),
    ("enorm", np.float64, (2,)),
    ("rnorm", np.float64, (2,)),
    ("dirfrac", np.float64, (2,)),
    ("denfrac", np.float64, (2,)),
], align=True)

# Main simulation parameters with nested data types
SIM_PARAMS_DTYPE = np.dtype([
    ("stat_params", STAT_PARAMS_DTYPE),
    ("cascade_params", CASCADE_PARAMS_DTYPE),
    ("recoil_params", RECOIL_PARAMS_DTYPE),
    ("geometry_params", GEOMETRY_PARAMS_DTYPE),
    ("estop_params", ESTOP_PARAMS_DTYPE),
    ("scatter_params", SCATTER_PARAMS_DTYPE),
    ("rng_seed", np.uint32),
], align=True)

# Preserve compatibility with vanilla NumPy (with numba disabled)
if os.environ.get("NUMBA_DISABLE_JIT", "") == "1":
    def Projectile(e, pos, dir, ispec=0, is_inside=True):
        """Create a single numpy record with initial properties of a Projectile
        
        This implementation is used when Numba is disabled to preserve compatability
        with vanilla NumPy and allow record field access via its attributes.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            ispec (int): Species index. Defaults to 0
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
        Returns:
            (numpy.record): A single record containing properties of a Projectile
        """
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
        """Create a single numpy record with initial properties of a Projectile
        
        This implementation is supported by Numba and creates a
        compatible instance of type np.void.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            ispec (int): Species index. Defaults to 0
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
        Returns:
            (numpy.void): A single record containing properties of a Projectile
        """
        rec = np.empty(1, dtype=PROJ_DTYPE)[0]
        rec['e'] = e
        rec['pos'] = pos    # copied
        rec['dir'] = dir    # copied
        rec['ispec'] = ispec
        rec['is_inside'] = is_inside
        return rec

class SimParams:
    def __init__(self, stat_params_tup=None, cascade_params_tup=None, recoil_params_tup=None, geometry_params_tup=None, estop_params_tup=None, scatter_params_tup=None, rng_seed=None):
        self.nspec = 0
        self.nbin = 0
        self.limits = (0.0, 0.0)
        self.emin = 0.0
        self.ed = 0.0
        self.pmax = 0.0
        self.mean_free_path = 0.0
        self.zmin = 0.0
        self.zmax = 0.0
        self.fac_linhard = np.zeros(2)
        self.density = 0.0
        self.pot_model = ""
        self.z1 = 0.0
        self.z2 = 0.0
        self.enorm = np.zeros(2)
        self.rnorm = np.zeros(2)
        self.dirfrac = np.zeros(2)
        self.denfrac = np.zeros(2)
        self.rng_seed = np.random.randint(2**31, dtype=np.uint32)
        
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
            self.pot_model, self.z1, self.z2, self.enorm, self.rnorm, self.dirfrac, self.denfrac = scatter_params_tup
        if rng_seed is not None:
            self.rng_seed = np.uint32(rng_seed)

    def to_record(self):
        rec = np.recarray(1, dtype=SIM_PARAMS_DTYPE)[0]
        rec['stat_params']['nspec'] = self.nspec
        rec['stat_params']['nbin'] = self.nbin
        rec['stat_params']['limits'] = np.array(self.limits)
        rec['cascade_params']['emin'] = self.emin
        rec['cascade_params']['ed'] = self.ed
        rec['recoil_params']['pmax'] = self.pmax
        rec['recoil_params']['mean_free_path'] = self.mean_free_path
        rec['geometry_params']['zmin'] = self.zmin
        rec['geometry_params']['zmax'] = self.zmax
        rec['estop_params']['fac_linhard'] = self.fac_linhard
        rec['estop_params']['density'] = self.density
        rec['scatter_params']['pot_model'] = self.pot_model
        rec['scatter_params']['z1'] = self.z1
        rec['scatter_params']['z2'] = self.z2
        rec['scatter_params']['enorm'] = self.enorm
        rec['scatter_params']['rnorm'] = self.rnorm
        rec['scatter_params']['dirfrac'] = self.dirfrac
        rec['scatter_params']['denfrac'] = self.denfrac
        rec['rng_seed'] = self.rng_seed
        return rec