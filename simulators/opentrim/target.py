"""Target related definitions and operations.

Currently, only planar targets are supported.

Available functions:
    is_inside_target: check if a given position is inside the target
"""
import numpy as np
from numba import jit


@jit(inline = "always")
def is_inside_target(pos, geometry_params):
    """Check if a given position is inside the target.

    Parameters:
        pos (ndarray): position to check (size 3)
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters

    Returns:
        (bool): whether the position is inside the target
    """
    return (geometry_params.z_intf[0] <= pos[2] 
            <= geometry_params.z_intf[geometry_params.nlayers])


@jit(inline = "always")
def get_layer_index(pos, geometry_params):
    """Get the layer index for a given position.

    For pos[2] < geometry_params.z_intf[0] return 0.
    For pos[2] >= geometry_params.z_intf[-1], return the last layer index.

    Note that material index = layer index.

    Parameters:
        pos (ndarray): position to check (size 3)
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters

    Returns:
        (int): layer index
    """
    for i in range(1, geometry_params.nlayers):
        if pos[2] < geometry_params.z_intf[i]:
            return i - 1    # return layer index to the left of the interface
        
    return geometry_params.nlayers - 1  # If not found in any layer, 
                                        # return last layer index


@jit(inline = "always")
def get_element_index(ilayer, materials_params):
    """Randomly select an element index for a given layer index.

    Parameters:
        ilayer (int): layer index
        materials_params (MATERIALS_PARAMS_DTYPE): materials parameters

    Returns:
        (int): element index
    """
    imat = ilayer

    r = np.random.rand() * sum(
        materials_params.atomic_fraction[imat, :materials_params.nelem[imat]])
    for ielem in range(materials_params.nelem[imat] - 1):
        cumulative_fraction = materials_params.cumulative_fraction[imat, ielem]
        if r < cumulative_fraction:
            return materials_params.ielem[imat, ielem]
    
    return materials_params.ielem[imat, materials_params.nelem[imat] - 1]