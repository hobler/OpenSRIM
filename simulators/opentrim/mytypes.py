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
        
        This implementation is used when Numba is disabled to preserve 
        compatability with vanilla NumPy and allow record field access via its 
        attributes.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            ispec (int): Species index. Defaults to 0
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
        Returns:
            (PROJ_TYPE): A single record containing properties of a Projectile
        """
        rec = np.recarray(1, dtype=PROJ_DTYPE)[0]
        rec["e"] = e
        rec["pos"] = pos    # copied
        rec["dir"] = dir    # copied
        rec["ispec"] = ispec
        rec["is_inside"] = is_inside
        return rec
else:
    @jit(inline = "always")
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
            (PROJ_TYPE): A single record containing properties of a Projectile
        """
        rec = np.empty(1, dtype=PROJ_DTYPE)[0]
        rec["e"] = e
        rec["pos"] = pos    # copied
        rec["dir"] = dir    # copied
        rec["ispec"] = ispec
        rec["is_inside"] = is_inside
        return rec

