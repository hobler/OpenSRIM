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
    ("fac_lindhard", np.float64, (2,)),
    ("density", np.float64),
], align=True)

SCATTER_PARAMS_DTYPE = np.dtype([
    ("pot_model", "<U16"),
    ("z1", np.uint32),
    ("z2", np.uint32),
    ("enorm", np.float64, (2,)),
    ("rnorm", np.float64, (2,)),
    ("dirfac", np.float64, (2,)),
    ("denfac", np.float64, (2,)),
    ("nlhlin_coefs", NLHLIN_COEFS_DTYPE, (4278,)),
], align=True)

# Main simulation parameters with nested data types
PARAMS_DTYPE = np.dtype([
    ("rng_seed", np.uint32),
    ("stat", STAT_PARAMS_DTYPE),
    ("cascade", CASCADE_PARAMS_DTYPE),
    ("recoil", RECOIL_PARAMS_DTYPE),
    ("geometry", GEOMETRY_PARAMS_DTYPE),
    ("estop", ESTOP_PARAMS_DTYPE),
    ("scatter", SCATTER_PARAMS_DTYPE)
], align=True)

# Preserve compatibility with vanilla NumPy (with nuparams.scattermba disabled)
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
    def __init__(self, stat_params=None, cascade_params=None, recoil_params=None, geometry_params=None, estop_params=None, scatter_params=None, rng_seed=None):
        self.nspec = 0
        self.nbin = 0
        self.limits = (0.0, 0.0)
        self.emin = 0.0
        self.ed = 0.0
        self.pmax = 0.0
        self.mean_free_path = 0.0
        self.zmin = 0.0
        self.zmax = 0.0
        self.fac_lindhard = np.zeros(2)
        self.density = 0.0
        self.pot_model = ""
        self.z1 = 0.0
        self.z2 = 0.0
        self.enorm = np.zeros(2)
        self.rnorm = np.zeros(2)
        self.dirfac = np.zeros(2)
        self.denfac = np.zeros(2)
        self.rng_seed = np.random.randint(2**31, dtype=np.uint32)
        self.nlhlin_coefs = np.zeros(4278, dtype=NLHLIN_COEFS_DTYPE)
        
        if stat_params is not None:
            self.nspec = stat_params.nspec
            self.nbin = stat_params.nbin
            self.limits = stat_params.limits
        if cascade_params is not None:
            self.emin = cascade_params.emin
            self.ed = cascade_params.ed
        if recoil_params is not None:
            self.pmax = recoil_params.pmax
            self.mean_free_path = recoil_params.mean_free_path
        if geometry_params is not None:
            self.zmin = geometry_params.zmin
            self.zmax = geometry_params.zmax     
        if estop_params is not None:
            self.fac_lindhard = estop_params.fac_lindhard
            self.density = estop_params.density
        if scatter_params is not None:
            self.pot_model = scatter_params.pot_model
            self.z1 = scatter_params.z1
            self.z2 = scatter_params.z2
            self.enorm = scatter_params.enorm
            self.rnorm = scatter_params.rnorm
            self.dirfac = scatter_params.dirfac
            self.denfac = scatter_params.denfac
            self.nlhlin_coefs = scatter_params.nlhlin_coefs
        if rng_seed is not None:
            self.rng_seed = np.uint32(rng_seed)

    def to_record(self):
        rec = np.recarray(1, dtype=PARAMS_DTYPE)[0]
        rec['stat']['nspec'] = self.nspec
        rec['stat']['nbin'] = self.nbin
        rec['stat']['limits'] = np.array(self.limits)
        rec['cascade']['emin'] = self.emin
        rec['cascade']['ed'] = self.ed
        rec['recoil']['pmax'] = self.pmax
        rec['recoil']['mean_free_path'] = self.mean_free_path
        rec['geometry']['zmin'] = self.zmin
        rec['geometry']['zmax'] = self.zmax
        rec['estop']['fac_lindhard'] = self.fac_lindhard
        rec['estop']['density'] = self.density
        rec['scatter']['pot_model'] = self.pot_model
        rec['scatter']['z1'] = self.z1
        rec['scatter']['z2'] = self.z2
        rec['scatter']['enorm'] = self.enorm
        rec['scatter']['rnorm'] = self.rnorm
        rec['scatter']['dirfac'] = self.dirfac
        rec['scatter']['denfac'] = self.denfac
        rec['scatter']['nlhlin_coefs'] = self.nlhlin_coefs
        rec['rng_seed'] = self.rng_seed
        return rec