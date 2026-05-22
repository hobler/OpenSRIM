import os
import numpy as np
from numba import jit, from_dtype
from . import config


PROJ_DTYPE = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("dffp_old", np.float64),
    ("dffp_new", np.float64),
    ("ielem", np.int32),
    ("ilayer", np.int32),
    ("is_inside", np.bool_),
    ("first_ffp", np.bool_)
], align=True)
PROJ_NUMBA_DTYPE = from_dtype(PROJ_DTYPE)


# Preserve compatibility with vanilla NumPy (with nuparams.scattermba disabled)
if os.environ.get("NUMBA_DISABLE_JIT", "") == "1":
    def Projectile(e, pos, dir, dffp_old=0.0, dffp_new=0.0, ielem=0, ilayer=0, 
                   is_inside=True, first_ffp=True):
        """Create a single numpy record with initial properties of a Projectile
        
        This implementation is used when Numba is disabled to preserve 
        compatability with vanilla NumPy and allow record field access via its 
        attributes.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            dffp_old (float): Length that has been added to the free flight
                path in the previous collision (often negative)
            dffp_new (float): Length that has to be added to the next free 
                flight path
            ielem (int): Species index
            ilayer (int): Layer index
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
            first_ffp (bool): If the first flight path is to be chosen randomly 
                between 0 and the mean value
        Returns:
            (PROJ_TYPE): A single record containing properties of a Projectile
        """
        rec = np.recarray(1, dtype=PROJ_DTYPE)[0]
        rec["e"] = e
        rec["pos"] = pos    # copied
        rec["dir"] = dir    # copied
        rec["dffp_old"] = dffp_old
        rec["dffp_new"] = dffp_new
        rec["ielem"] = ielem
        rec["ilayer"] = ilayer
        rec["is_inside"] = is_inside
        rec["first_ffp"] = first_ffp
        return rec
else:
    @jit(inline = "always", debug=config.DEBUG)
    def Projectile(e, pos, dir, dffp_old=0.0, dffp_new=0.0, ielem=0, ilayer=0, 
                   is_inside=True, first_ffp=True):
        """Create a single numpy record with initial properties of a Projectile
        
        This implementation is supported by Numba and creates a
        compatible instance of type np.void.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            dffp_old (float): Length that has been added to the free flight
                path in the previous collision (often negative)
            dffp_new (float): Length that has to be added to the next free 
                flight path
            ielem (int): Species index. Defaults to 0
            ilayer (int): Layer index. Defaults to 0
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
            first_ffp (bool): If the first flight path is to be chosen randomly 
                between 0 and the mean value. Defaults to True
        Returns:
            (PROJ_TYPE): A single record containing properties of a Projectile
        """
        rec = np.empty(1, dtype=PROJ_NUMBA_DTYPE)[0]
        rec["e"] = e
        rec["pos"] = pos    # copied
        rec["dir"] = dir    # copied
        rec["dffp_old"] = dffp_old
        rec["dffp_new"] = dffp_new
        rec["ielem"] = ielem
        rec["ilayer"] = ilayer
        rec["is_inside"] = is_inside
        rec["first_ffp"] = first_ffp
        return rec

