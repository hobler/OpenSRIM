"""Target-geometry related operations.

Currently, only a planar target geometry is supported.

Available functions:
    setup: setup module variables.
    is_inside_target: check if a given position is inside the target
"""


def setup(zmin, zmax):
    """Define the geometry of the target.
    
    Parameters:
        zmin (float): minimum z coordinate of the target (A)
        zmax (float): maximum z coordinate of the target (A)
    """
    global ZMIN, ZMAX

    ZMIN = zmin
    ZMAX = zmax


def is_inside_target(pos):
    """Check if a given position is inside the target.

    Parameters:
        pos (NDArray, Shape(3)): position to check

    Returns:
        (bool): whether the position is inside the target
    """
    if ZMIN <= pos[2] <= ZMAX:
        return True
    else:
        return False