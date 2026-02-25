"""Treat the scattering of a projectile on a target atom.

Currently, the ZBL potential (Ziegler, Biersack, Littmark, The Stopping and 
Range of Ions in Matter, Pergamon Press, 1985) and the NLHlin potential
(to be published) are implemented. The scattering integrals are evaluated
numerically, or Biersack"s "magic formula" is usedfor the scattering angle.

Available functions:
    scatter: treat a scattering event.
"""
import os
import math
from collections import namedtuple
import numpy as np
from numba import jit
from .zbl import magic
from .cm_scatter import scatter_integrals


@jit(inline = "always")
def normalize_if_needed(vec, fallback):
    """Fast normalization with fallback – in‑place, no new array.
    
    Parameters:
        vec (np.ndarray): Vector to be normalized (size 3)
        fallback (np.ndarray): Vector to replace the original with if 
            len(norm) == 0
        
    Returns:
        np.ndarray: The normalized vector (size 3)
    """
    norm_sq = vec[0]**2 + vec[1]**2 + vec[2]**2
    if norm_sq == 0.0:
        # fallback is already a unit vector, just copy it
        vec[0], vec[1], vec[2] = fallback[0], fallback[1], fallback[2]
        return vec
    norm = math.sqrt(norm_sq)
    vec[0] /= norm
    vec[1] /= norm
    vec[2] /= norm
    return vec


@jit
def scatter(proj, p, dirp, ielem2, scatter_params):
    """Treat a scattering event.

    The atomic numbers and masses of the ion and the target atom enter the
    calculation via the ENORM, PNORM, DIRFAC, and DENFAC parameters.

    The direction vectors proj["dir"] and dirp[:] are assumed to be normalized 
    to unit length.

    Parameters:
        proj (Projectile): state of the projectile (will be modified in-place)
        p (float): impact parameter (A)
        dirp (ndarray): direction vector of the impact parameter
            (= from the collision point to the recoil position before 
            the collision) (unit vector, size 3)
        ielem2 (int): index of the target element
        scatter_params (np.recarray): Scatter parameters
        is_magic (bool): If magic function should be used (otherwise 
            scatter_integrals)
    
    Returns:
        (ndarray): direction vector of the recoil after the collision 
            (size 3)
        (float): energy of the projectile after the collision
    """
    # scattering angle theta in the center-of-mass system
    ielem1 = proj["ielem"]
    proj_e = proj["e"]
    proj_dir = proj["dir"][:]

    enorm = scatter_params.enorm[ielem1, ielem2]
    rnorm = scatter_params.rnorm[ielem1, ielem2]
    dirfac = scatter_params.dirfac[ielem1, ielem2]
    denfac = scatter_params.denfac[ielem1, ielem2]
    pot_model = scatter_params.pot_model
    
    if pot_model == "ZBL_magic":
        cos_half_theta = magic(proj_e/enorm, p/rnorm)
        sin_half_theta = math.sqrt(1 - cos_half_theta**2)
    elif pot_model.endswith("magic"):
        raise ValueError(f"Unknown potential model {pot_model}")
    else:
        theta, _ = scatter_integrals(proj_e/enorm, p/rnorm, pot_model)
        sin_half_theta = math.sin(0.5 * theta)
        cos_half_theta = math.cos(0.5 * theta)

    # directions of the recoil and the projectile after the collision
    recoil_dir = dirfac * sin_half_theta * (
        sin_half_theta*proj_dir[:] + cos_half_theta*dirp[:])
    dir_new = proj_dir[:] - recoil_dir[:]
    dir_new = normalize_if_needed(dir_new, proj_dir[:])
    recoil_dir = normalize_if_needed(recoil_dir, proj_dir[:])

    # Copy dir_new buffer content into proj["dir"] buffer
    proj["dir"][:] = dir_new

    # energy after scattering
    recoil_e = denfac * proj_e * sin_half_theta**2
    proj["e"] -= recoil_e

    return recoil_dir[:], recoil_e