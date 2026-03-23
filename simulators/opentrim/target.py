"""Target related definitions and operations.

Currently, only planar targets are supported.

Available functions:
    is_inside_target: check if a given position is inside the target
    get_layer_index: get the layer index for a given position
    get_element_index: randomly select an element index for a given projectile's 
        layer index
"""
import numpy as np
from numba import jit


@jit(inline = "always")
def is_inside_target(pos, geometry_params):
    """Check if a given position is inside the target.

    Parameters:
        pos (ndarray): Position to check (size 3)
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters

    Returns:
        (bool): True if the position is inside the target, False otherwise
    """
    return (geometry_params.x_intf[0] <= pos[0] 
            <= geometry_params.x_intf[geometry_params.nlayers])


@jit(inline = "always")
def get_layer_index(pos, geometry_params):
    """Get the layer index for a given position.

    For pos[0] < geometry_params.x_intf[0] return 0.
    For pos[0] >= geometry_params.x_intf[-1], return the last layer index.

    Note that material index = layer index.

    Parameters:
        pos (ndarray): Position to check (size 3)
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters

    Returns:
        (int): Layer index
    """
    for i in range(1, geometry_params.nlayers):
        if pos[0] < geometry_params.x_intf[i]:
            return i - 1    # return layer index to the left of the interface
        
    return geometry_params.nlayers - 1  # If not found in any layer, 
                                        # return last layer index


@jit(inline = "always")
def get_element_index(proj, materials_params):
    """Randomly select an element index for a given projectile's layer index.

    Parameters:
        proj (Projectile): Projectile to get and set the element index for
        materials_params (MATERIALS_PARAMS_DTYPE): materials parameters

    Returns:
        (int): element index
    """
    imat = proj["ilayer"]
    nelem_mat = materials_params.nelem[imat]

    r = np.random.rand() * sum(
        materials_params.atomic_fraction[imat, :nelem_mat])
    for ielem in range(nelem_mat - 1):
        if r < materials_params.cumulative_fraction[imat, ielem]:
            return materials_params.ielem[imat, ielem]

    return materials_params.ielem[imat, nelem_mat - 1]