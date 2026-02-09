"""Treat the scattering of a projectile on a target atom.

Currently, only the ZBL potential (Ziegler, Biersack, Littmark,
The Stopping and Range of Ions in Matter, Pergamon Press, 1985) is 
implemented, along with Biersack's "magic formula" for the scattering 
angle.

Available functions:
    setup: setup module variables.
    scatter: treat a scattering event.
"""

import math
from zbl import magic
from cm_scatter import scatter_integrals
import numpy as np
from numba import jit

@jit(inline = 'always')
def normalize_if_needed(vec, fallback):
    """Fast normalization with fallback – in‑place, no new array."""
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
def scatter(proj, p, dirp, screen_fun, scatter_params, is_magic):
    """Treat a scattering event.

    The atomic numbers and masses of the ion and the target atom enter the
    calculation via the module variables ENORM, PNORM, DIRFAC, and DENFAC.

    The direction vectors proj.dir and dirp are assumed to be normalized to 
    unit length.

    Parameters:
        proj (Projectile): state of the projectile before the collision
        p (float): impact parameter (A)
        dirp (ndarray): direction vector of the impact parameter
            (= from the collision point to the recoil position before 
            the collision) (unit vector, size 3)
        screen_fun (object): Screening function
        scatter_params (np.recarray): Scatter parameters
        is_magic (bool): If magic function should be used (otherwise scatter_integrals)
    
    Returns:
        (Projectile): state of the projectile after the collision 
        (ndarray): direction vector of the recoil after the collision 
            (size 3)
        (float): energy of the projectile after the collision
    """
    # scattering angle theta in the center-of-mass system
    enorm = scatter_params.enorm
    rnorm = scatter_params.rnorm
    dirfrac = scatter_params.dirfrac
    denfrac = scatter_params.denfrac
    
    ispec = proj.ispec
    proj_e = proj.e
    if is_magic:
        cos_half_theta = magic(proj_e/enorm[ispec], 
                               p/rnorm[ispec],
                               screen_fun[ispec])
        sin_half_theta = math.sqrt(1 - cos_half_theta**2)
    else:
        theta, _ = scatter_integrals(proj_e/enorm[ispec], 
                                     p/rnorm[ispec], 
                                     screen_fun[ispec])
        sin_half_theta = math.sin(0.5 * theta)
        cos_half_theta = math.cos(0.5 * theta)

    # directions of the recoil and the projectile after the collision
    recoil_dir = dirfrac[ispec] * sin_half_theta * (sin_half_theta*proj.dir[:] 
                                                 + cos_half_theta*dirp[:])
    dir_new = proj.dir[:] - recoil_dir[:]
    dir_new = normalize_if_needed(dir_new, proj['dir'][:])
    recoil_dir = normalize_if_needed(recoil_dir, proj['dir'][:])

    # Copy dir_new buffer content into proj.dir buffer
    proj.dir[:] = dir_new

    # energy after scattering
    recoil_e = denfrac[ispec] * proj_e * sin_half_theta**2
    proj.e -= recoil_e

    return recoil_dir[:], recoil_e

# Excluded from being JIT-Compiled
def setup(z1, m1, z2, m2, pot_model):
    """Setup module variables depending on projectile and target species.

    Each of the module variables ENORM, RNORM, DIRFAC, and DENFAC is a tuple
    with two entries: one for the ion species 0 and one for moving atom 
    species 1. Currently we assume there is only on target atoms species.

    Parameters:
        z1 (int): atomic number of projectile
        m1 (float): mass of projectile (amu)
        z2 (int): atomic number of target
        m2 (float): mass of target (amu)
        pot_model (str): potential model for scattering
        
    Returns:
        (str): Model identifier (name)
        (int): Z1
        (int): Z2
        (np.ndarray): ENORM
        (np.ndarray): RNORM
        (np.ndarray): DIRFAC
        (np.ndarray): DENFAC
    """
    m1_m2 = m1 / m2
    if pot_model.startswith('ZBL'):
        rnorm = (0.4685 / (z1**0.23 + z2**0.23),
                 0.4685 / (z2**0.23 + z2**0.23))                  # A
    else:
        rnorm = (0.4685 / math.sqrt(math.sqrt(z1) + math.sqrt(z2)),
                 0.4685 / math.sqrt(math.sqrt(z2) + math.sqrt(z2)))     # A
    enorm = np.array((14.39979 * z1 * z2 / rnorm[0] * (1 + m1_m2),
                14.39979 * z2 * z2 / rnorm[1] * (1 + 1)))            # eV
    dirfac = np.array((2 / (1 + m1_m2),
                1))
    denfac = np.array((4 * m1_m2 / (1 + m1_m2)**2,
                1))
              
    return pot_model, z1, z2, enorm, rnorm, dirfac, denfac