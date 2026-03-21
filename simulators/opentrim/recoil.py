"""Create the recoil position for the next collision.

Currently, only amorphous targets are supported. The free path length to
the next collision is assumed to be constant and equal to the atomic
density to the power -1/3.

Available functions:
    get_recoil_position: get the recoil position.
"""
from math import sqrt, sin, cos
import numpy as np
from numba import jit
from .mytypes import PROJ_DTYPE
from .target import get_layer_index, get_element_index, is_inside_target

@jit
def get_recoil(proj, recoil, params):
    """Get the position of the recoil hit after the next free flight path.

    The recoil more precisely is a recoil candidate, since it is not guaranteed 
    that the recoil has enough energy to leave its position.

    The recoil position is determined by a deterministic free path length 
    and sampling a random impact parameter.

    We cannot return a recoil structured array here, since Numba apparently 
    does not allow returning structured arrays from jit functions. Instead, we 
    return the recoil position as a separate array, and the caller can 
    construct the recoil structured array if needed.

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
    pos = proj["pos"][:]
    dir = proj["dir"][:]
    ilayer = proj["ilayer"]
    
    # free flight path and impact parameter
    free_path = params.recoil.mean_free_path[ilayer]
    collision_pos = pos[:] + free_path * dir[:]
    p = params.recoil.pmax[ilayer] * sqrt(np.random.rand())

    # Azimuthal angle fi
    fi = 2 * np.pi * np.random.rand()
    cos_fi = cos(fi)
    sin_fi = sin(fi)

    # Convert direction vector to polar angles
    # make k point to the smallest dir(:) so sinalf > sqrt(2/3)
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
    recoil["ilayer"] = get_layer_index(recoil["pos"], params.geometry)
    recoil["ielem"] = get_element_index(recoil, params.materials)
    recoil["is_inside"] = is_inside_target(recoil["pos"], params.geometry)

    return free_path, p, dirp[:]