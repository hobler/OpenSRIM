import os
import numpy as np
from numba import jit, from_dtype


PROJ_DTYPE = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("ielem", np.int32),
    ("ilayer", np.int32),
    ("is_inside", np.bool_)
], align=True)
PROJ_NUMBA_DTYPE = from_dtype(PROJ_DTYPE)


# Preserve compatibility with vanilla NumPy (with nuparams.scattermba disabled)
if os.environ.get("NUMBA_DISABLE_JIT", "") == "1":
    def Projectile(e, pos, dir, ielem=0, ilayer=0, is_inside=True):
        """Create a single numpy record with initial properties of a Projectile
        
        This implementation is used when Numba is disabled to preserve 
        compatability with vanilla NumPy and allow record field access via its 
        attributes.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            ielem (int): Species index. Defaults to 0
            ilayer (int): Layer index. Defaults to 0
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
        Returns:
            (PROJ_TYPE): A single record containing properties of a Projectile
        """
        rec = np.recarray(1, dtype=PROJ_DTYPE)[0]
        rec["e"] = e
        rec["pos"] = pos    # copied
        rec["dir"] = dir    # copied
        rec["ielem"] = ielem
        rec["ilayer"] = ilayer
        rec["is_inside"] = is_inside
        return rec
else:
    @jit(inline = "always")
    def Projectile(e, pos, dir, ielem=0, ilayer=0, is_inside=True):
        """Create a single numpy record with initial properties of a Projectile
        
        This implementation is supported by Numba and creates a
        compatible instance of type np.void.
    
        Parameters:
            e (float): Energy of the projectile
            pos (np.ndarray[float]): Current position vector
            dir (np.ndarray[float]): Current direction vector
            ielem (int): Species index. Defaults to 0
            ilayer (int): Layer index. Defaults to 0
            is_inside (bool): If the projectile is within the simulation area.
                Defaults to True
        
        Returns:
            (PROJ_TYPE): A single record containing properties of a Projectile
        """
        rec = np.empty(1, dtype=PROJ_DTYPE)[0]
        rec["e"] = e
        rec["pos"] = pos    # copied
        rec["dir"] = dir    # copied
        rec["ielem"] = ielem
        rec["ilayer"] = ilayer
        rec["is_inside"] = is_inside
        return rec

