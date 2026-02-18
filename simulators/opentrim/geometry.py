"""Target-geometry related operations.

Currently, only a planar target geometry is supported.

Available functions:
    setup: setup module variables.
    is_inside_target: check if a given position is inside the target
"""
import numpy as np
from numba import jit
from mytypes import GEOMETRY_PARAMS_DTYPE


def setup(zmin, zmax):
    """Define the geometry of the target.
    
    Parameters:
        zmin (float): minimum z coordinate of the target (A)
        zmax (float): maximum z coordinate of the target (A)
    Returns:
        (GEOMETRY_PARAMS_DTYPE): geometry parameters
    """
    geometry_params = np.recarray(1, dtype=GEOMETRY_PARAMS_DTYPE)[0]
    geometry_params["zmin"] = zmin
    geometry_params["zmax"] = zmax

    return geometry_params


@jit(inline = 'always')
def is_inside_target(pos, geometry):
    """Check if a given position is inside the target.

    Parameters:
        pos (ndarray): position to check (size 3)
        geometry (np.recarray): Geometry parameters

    Returns:
        (bool): whether the position is inside the target
    """
    return geometry.zmin <= pos[2] <= geometry.zmax