"""Create the recoil position for the next collision.

Currently, only amorphous targets are supported. The free path length to
the next collision is assumed to be constant and equal to the atomic
density to the power -1/3.

Available functions:
    setup: setup module variables.
    get_recoil_position_position: get the recoil position.
"""
from math import sqrt, sin, cos
import numpy as np
from numba import jit


def setup(input_params):
    """Setup module variables depending on target density.

    Parameters:
        input_params (dict): input parameters dictionary

    Returns:
        (RECOIL_PARAMS_DTYPE): Recoil parameters
    """
    densities = np.array(input_params["layers"]["density"])
    nlayers = len(densities)
    
    RECOIL_PARAMS_DTYPE = np.dtype([
        ("pmax", np.float64, (nlayers,)),
        ("mean_free_path", np.float64, (nlayers,)),
    ], align=True)

    recoil_params = np.recarray(1, dtype=RECOIL_PARAMS_DTYPE)[0]
    recoil_params["mean_free_path"] = densities**(-1/3)
    recoil_params["pmax"] = recoil_params["mean_free_path"] / sqrt(np.pi)
    
    return recoil_params


@jit
def get_recoil_position(pos, dir, params):
    """Get the recoil hit after the next free flight path.

    The recoil more precisely is a recoil candidate, since it is not guaranteed 
    that the recoil has enough energy to leave its position.

    The recoil position is determined by a deterministicrandom free path length 
    and sampling a random impact parameter,

    Parameters:
        pos (ndarray): position of the projectile (size 3)
        dir (ndarray): direction vector of the projectile (size 3)
        params (RECOIL_PARAMS_DTYPE): Recoil parameters

    Returns:
        (float): free path length to the next collision (A)
        (float): impact parameter = distance between collision point and 
            recoil (A)
        (ndarray): direction vector from collision point to recoil (size 3)
        (Projectile): the recoil
    """
    free_path = params.recoil.mean_free_path[0]  # TODO: use the correct layer index
    collision_pos = pos[:] + free_path * dir[:]

    p = params.recoil.pmax[0] * sqrt(np.random.rand())  # TODO: use the correct layer index
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

    # recoil properties (except energy and direction)
    recoil_pos = collision_pos[:] + p * dirp[:]

    return free_path, p, dirp[:], recoil_pos