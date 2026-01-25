"""Treat the scattering of a projectile on a target atom.

Currently, only the ZBL potential (Ziegler, Biersack, Littmark,
The Stopping and Range of Ions in Matter, Pergamon Press, 1985) is 
implemented, along with Biersack's "magic formula" for the scattering 
angle.

Available functions:
    setup: setup module variables.
    scatter: treat a scattering event.
"""

from math import sqrt, exp
from zblpot import magic
import numpy as np


def setup(z1, m1, z2, m2):
    """Setup module variables depending on projectile and target species.

    Each of the module variables ENORM, RNORM, DIRFAC, and DENFAC is a tuple
    with two entries: one for the ion species 0 and one for moving atom 
    species 1. Currently we assume there is only on target atoms species.

    Parameters:
        z1 (int): atomic number of projectile
        m1 (float): mass of projectile (amu)
        z2 (int): atomic number of target
        m2 (float): mass of target (amu)
    """
    global ENORM, RNORM, DIRFAC, DENFAC

    m1_m2 = m1 / m2
    RNORM = (0.4685 / (z1**0.23 + z2**0.23),
             0.4685 / (z2**0.23 + z2**0.23))                  # A
    ENORM = (14.39979 * z1 * z2 / RNORM[0] * (1 + m1_m2),
             14.39979 * z2 * z2 / RNORM[1] * (1 + 1))            # eV
    DIRFAC = (2 / (1 + m1_m2),
              1)
    DENFAC = (4 * m1_m2 / (1 + m1_m2)**2,
              1)
        
        
def scatter(proj, p, dirp):
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
    
    Returns:
        (Projectile): state of the projectile after the collision 
        (ndarray): direction vector of the recoil after the collision 
            (size 3)
        (float): energy of the projectile after the collision
    """
    # scattering angle theta in the center-of-mass system
    sin_half_theta, cos_half_theta = magic(proj.e/ENORM[proj.ispec], 
                                           p/RNORM[proj.ispec])

    # directions of the recoil and the projectile after the collision
    sin_psi = cos_half_theta
    cos_psi = sin_half_theta
    recoil_dir = DIRFAC[proj.ispec] * cos_psi * (cos_psi*proj.dir[:] 
                                                 + sin_psi*dirp[:])
    dir_new = proj.dir[:] - recoil_dir[:]
    norm = np.linalg.norm(dir_new[:])
    if norm == 0:
        dir_new = proj.dir[:]
    else:
        dir_new /= norm
    norm = np.linalg.norm(recoil_dir[:])
    if norm == 0:
        recoil_dir = proj.dir[:]
    else:
        recoil_dir /= norm
    proj.dir = dir_new[:]

    # energy after scattering
    recoil_e = DENFAC[proj.ispec] * proj.e * (1 - cos_half_theta**2)
    proj.e -= recoil_e

    return proj, recoil_dir[:], recoil_e