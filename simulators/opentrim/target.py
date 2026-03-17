"""Target related definitions and operations.

Currently, only planar targets are supported.

Available functions:
    is_inside_target: check if a given position is inside the target
"""
import numpy as np
from numba import jit


@jit(inline = "always")
def set_is_inside_target(proj, geometry_params):
    """Check if a given position is inside the target.

    Parameters:
        proj (Projectile): Projectile to check against the geometry
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters
    """
    proj["is_inside"] = (geometry_params.x_intf[0] <= proj["pos"][0] 
            <= geometry_params.x_intf[geometry_params.nlayers])


@jit(inline = "always")
def set_layer_index(proj, geometry_params):
    """Get the layer index for a given projectile and apply it.

    For pos[0] < geometry_params.x_intf[0] return 0.
    For pos[0] >= geometry_params.x_intf[-1], return the last layer index.

    Note that material index = layer index.

    Parameters:
        proj (Projectile): Projectile to get and set the index for
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters
    """
    for i in range(1, geometry_params.nlayers):
        if proj["pos"][0] < geometry_params.x_intf[i]:
            proj["ilayer"] = i - 1    # return layer index to the left of the interface
        
    proj["ilayer"] = geometry_params.nlayers - 1  # If not found in any layer, 
                                        # return last layer index


@jit(inline = "always")
def set_element_index(proj, materials_params):
    """Randomly select an element index for a given projectile's layer index.

    Parameters:
        proj (Projectile): Projectile to get and set the element index for
        materials_params (MATERIALS_PARAMS_DTYPE): materials parameters

    Returns:
        (int): element index
    """
    imat = proj["ilayer"]

    r = np.random.rand() * sum(
        materials_params.atomic_fraction[imat, :materials_params.nelem[imat]])
    for ielem in range(materials_params.nelem[imat] - 1):
        cumulative_fraction = materials_params.cumulative_fraction[imat, ielem]
        if r < cumulative_fraction:
            proj["ielem"] = materials_params.ielem[imat, ielem]
    
    proj["ielem"] = materials_params.ielem[imat, materials_params.nelem[imat] - 1]