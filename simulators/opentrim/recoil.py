"""Create the recoil position for the next collision.

Currently, only amorphous targets are supported. The free path length to
the next collision is assumed to be constant and equal to the atomic
density to the power -1/3.

Available functions:
    select_recoil: select the recoil position and atom species.
"""
from math import sqrt, sin, cos
import numpy as np
from numba import jit
from .mytypes import PROJ_DTYPE
from .target import get_layer_index, get_element_index, is_inside_target
from . import config


@jit(debug=config.DEBUG)
def select_recoil(proj, recoil, params):
    """Get the position of the recoil hit after the next free flight path.

    The recoil more precisely is a recoil candidate, since it is not guaranteed 
    that the recoil has enough energy to leave its position.

    The recoil position is determined by a deterministic free path length 
    and sampling a random impact parameter.

    Due to a limitation of Numba, we cannot return a recoil structured array 
    if we create it here. Instead, we receive a recoil structured array as an 
    argument and modify it in-place.
    
    Parameters:
        proj (Projectile): state of the projectile
        recoil (Projectile): the recoil with position and is_inside, element 
            and layer index, set (modified in-place)
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (float): free path length to the next collision (A)
        (float): impact parameter = distance between collision point and 
            recoil (A)
        (ndarray): direction vector from collision point to recoil (size 3)
    """
    e = proj["e"]
    pos = proj["pos"][:]
    dir = proj["dir"][:]
    ielem1 = proj["ielem"]
    ilayer = proj["ilayer"]
    imat = ilayer
    
    # free flight path and impact parameter
    #pmax = params.cascade.pmax[ilayer]
    if proj["is_inside"]:
        pmax = np.interp(e, params.cascade.pmax_energies[ielem1, imat], 
                         params.cascade.pmax_vals)
    else:
        pmax = np.interp(e, params.cascade.pmax_energies_surface[ielem1, imat], 
                         params.cascade.pmax_vals)

    #free_path = params.cascade.mean_free_path[ilayer]
    free_path = 1 / (params.materials[imat].density * np.pi * pmax**2)
    if False:  # determistic free flight path except first one
        if proj["first_ffp"]:
            free_path *= np.random.rand()
        proj["first_ffp"] = False
    else:      # statistical free flight path
        free_path *= -np.log(np.random.rand())

    p = pmax * sqrt(np.random.rand())
    collision_pos = pos[:] + free_path * dir[:]

    # Azimuthal angle fi
    fi = 2 * np.pi * np.random.rand()
    cos_fi = cos(fi)
    sin_fi = sin(fi)

    # Convert direction vector to polar angles
    # make k point to the smallest dir[:] so sin_alpha > sqrt(2/3)
    k = np.argmin(np.abs(dir[:]))
    i = (k + 1) % 3
    j = (i + 1) % 3
    cos_alpha = dir[k]
    sin_alpha = sqrt( dir[i]**2 + dir[j]**2 )
    cos_phi = dir[i] / sin_alpha
    sin_phi = dir[j] / sin_alpha

    # direction vector from collision point to recoil
    dirp = np.empty(3)
    dirp[i] = cos_fi*cos_alpha*cos_phi - sin_fi*sin_phi
    dirp[j] = cos_fi*cos_alpha*sin_phi + sin_fi*cos_phi
    dirp[k] = - cos_fi*sin_alpha
    norm = np.linalg.norm(dirp)
    dirp /= norm

    # recoil position
    recoil["pos"] = collision_pos[:] + p * dirp[:]
    recoil["ilayer"] = get_layer_index(recoil["pos"], params)
    recoil["ielem"] = get_element_index(recoil, params)
    recoil["is_inside"] = is_inside_target(recoil["pos"], params)

    return free_path, p, dirp[:]