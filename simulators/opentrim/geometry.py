"""Target-geometry related operations.

Currently, only a planar target geometry is supported.

Available functions:
    setup: setup module variables.
    is_inside_target: check if a given position is inside the target
"""
import numpy as np
from numba import jit


def setup(zmin, zmax):
    """Define the geometry of the target.
    
    Parameters:
        zmin (float): minimum z coordinate of the target (A)
        zmax (float): maximum z coordinate of the target (A)
    Returns:
        (GEOMETRY_PARAMS_DTYPE): geometry parameters
    """
    GEOMETRY_PARAMS_DTYPE = np.dtype([
        ("zmin", np.float64),
        ("zmax", np.float64),
    ], align=True)

    geometry_params = np.recarray(1, dtype=GEOMETRY_PARAMS_DTYPE)[0]
    geometry_params["zmin"] = zmin
    geometry_params["zmax"] = zmax

    return geometry_params


@jit(inline = "always")
def is_inside_target(pos, params):
    """Check if a given position is inside the target.

    Parameters:
        pos (ndarray): position to check (size 3)
        params (GEOMETRY_PARAMS_DTYPE): Geometry parameters

    Returns:
        (bool): whether the position is inside the target
    """
    return params.zmin <= pos[2] <= params.zmax