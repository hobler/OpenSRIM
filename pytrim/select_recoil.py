"""Create the recoil position for the next collision.

Currently, only amorphous targets are supported. The free path length to
the next collision is assumed to be constant and equal to the atomic
density to the power -1/3.

Available functions:
    setup: setup module variables.
    get_recoil_position: get the recoil position.
"""
from math import sqrt, sin, cos
import numpy as np


def setup(density):
    """Setup module variables depending on target density.

    Parameters:
        density (float): target density (atoms/A^3)

    Returns:
        None    
    """
    global PMAX, MEAN_FREE_PATH

    MEAN_FREE_PATH = density**(-1/3)
    PMAX = MEAN_FREE_PATH / sqrt(np.pi)


def get_recoil_position(pos, dir):
    """Get the recoil position based on the projectile position and direction.

    Parameters:
        pos (ndarray): position of the projectile (size 3)
        dir (ndarray): direction vector of the projectile (size 3)

    Returns:
        float: free path length to the next collision (A)
        float: impact parameter = distance between collision point and 
            recoil (A)
        ndarray: direction vector from collision point to recoil (size 3)
        ndarray: position of the recoil (size 3)
    """
    free_path = MEAN_FREE_PATH
    pos_collision = pos[:] + free_path * dir[:]

    p = PMAX * sqrt(np.random.rand())
    # Azimuthal angle fi
    fi = 2 * np.pi * np.random.rand()
    cos_fi = cos(fi)
    sin_fi = sin(fi)

    # Convert direction vector to polar angles
    k = np.argmin( np.abs(dir[:]) )   # make k point to the smallest dir(:) so sinalf > sqrt(2/3)
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

    # position of the recoil
    pos_recoil = pos_collision[:] + p * dirp[:]

    return free_path, p, dirp, pos_recoil
