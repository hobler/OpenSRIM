import os
import numpy as np
from numba import jit


PROJ_DTYPE = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("ielem", np.int32),
    ("ilayer", np.int32),
    ("is_inside", np.bool_)
], align=True)


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


HIST_CONFIG_DTYPE = np.dtype([
    ("nbin", np.int32),
    ("limits_min", np.float64),
    ("limits_max", np.float64),
    ("bin_width", np.float64),
    ("offset", np.int32),
    ("counts_size", np.int32)
], align=True)


@jit(nopython=True)
def create_histogram_configs(nbin_list, limits_min_list, limits_max_list, nvar):
    """
    Create structured array of histogram configurations.
    
    Parameters:
        nbin_list: Array of integers, number of bins per histogram
        limits_min_list: Array of minimum limits per histogram
        limits_max_list: Array of maximum limits per histogram
        nvar: Number of variables (species)
    
    Returns:
        tuple: (configs, total_size)
            configs: Structured array of configurations
            total_size: Total size of flattened counts array
    """
    n_hist = len(nbin_list)
    configs = np.zeros(n_hist, dtype=HIST_CONFIG_DTYPE)
    
    offset = 0
    for i in range(n_hist):
        nbin = nbin_list[i]
        lmin = limits_min_list[i]
        lmax = limits_max_list[i]
        
        configs[i]["nbin"] = nbin
        configs[i]["limits_min"] = lmin
        configs[i]["limits_max"] = lmax
        configs[i]["bin_width"] = (lmax - lmin) / nbin
        configs[i]["offset"] = offset
        configs[i]["counts_size"] = nvar * (nbin + 2)
        offset += configs[i]["counts_size"]
    
    return configs, offset

