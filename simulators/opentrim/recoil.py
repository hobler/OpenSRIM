"""Create the recoil position for the next collision.

Currently, only amorphous targets are supported. The free path length to
the next collision is assumed to be constant and equal to the atomic
density to the power -1/3.

Available functions:
    setup: setup module variables.
    get_recoil_position_position: get the recoil position.
"""
from math import sqrt, sin, cos
from collections import namedtuple
import numpy as np
from numba import jit
from .mytypes import PROJ_DTYPE


def setup(input_params):
    """Setup module variables depending on target density.

    Parameters:
        input_params (dict): input parameters dictionary

    Returns:
        (RECOIL_PARAMS_DTYPE): Recoil parameters
    """
    densities = np.array(input_params["layers"]["density"])
    nlayers = len(densities)
    
    #RECOIL_PARAMS_DTYPE = np.dtype([
    #    ("pmax", np.float64, (nlayers,)),
    #    ("mean_free_path", np.float64, (nlayers,)),
    #], align=True)

    #recoil_params = np.recarray(1, dtype=RECOIL_PARAMS_DTYPE)[0]
    #recoil_params["mean_free_path"] = densities**(-1/3)
    #recoil_params["pmax"] = recoil_params["mean_free_path"] / sqrt(np.pi)

    RecoilParams = namedtuple("RecoilParams", ["pmax", "mean_free_path"])
    recoil_params = RecoilParams(
        pmax = densities**(-1/3) / sqrt(np.pi),
        mean_free_path = densities**(-1/3),
    )

    return recoil_params


@jit
def get_recoil_position(proj, recoil_params):
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
        recoil_params (RECOIL_PARAMS_DTYPE): Recoil parameters

    Returns:
        (float): free path length to the next collision (A)
        (float): impact parameter = distance between collision point and 
            recoil (A)
        (ndarray): direction vector from collision point to recoil (size 3)
        (ndarray): the recoil position (size 3)
    """
    pos = proj["pos"][:]
    dir = proj["dir"][:]
    ilayer = proj["ilayer"]
    
    # free flight path and impact parameter
    free_path = recoil_params.mean_free_path[ilayer]
    collision_pos = pos[:] + free_path * dir[:]
    p = recoil_params.pmax[ilayer] * sqrt(np.random.rand())

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
    recoil_pos = collision_pos[:] + p * dirp[:]

    return free_path, p, dirp[:], recoil_pos