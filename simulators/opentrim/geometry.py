"""Target-geometry related operations.

Currently, only a planar target geometry is supported.

Available functions:
    setup: setup module variables.
    is_inside_target: check if a given position is inside the target
"""

from numba import jit

def setup(zmin, zmax):
    """Define the geometry of the target.
    
    Parameters:
        zmin (float): minimum z coordinate of the target (A)
        zmax (float): maximum z coordinate of the target (A)
    Returns:
        (float): zmin
        (float): zmax
    """
    return zmin, zmax

@jit(inline = 'always')
def is_inside_target(pos, geometry_params):
    """Check if a given position is inside the target.

    Parameters:
        pos (ndarray): position to check (size 3)
        geometry_params (np.recarray): Geometry parameters

    Returns:
        (bool): whether the position is inside the target
    """
    return geometry_params.zmin <= pos[2] <= geometry_params.zmax