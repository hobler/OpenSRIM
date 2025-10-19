"""Simulate projectile trajectories.

Available functions:
    setup: setup module variables.
    trajectory: simulate one trajectory."""
from select_recoil import get_recoil_position
from scatter import scatter
from estop import eloss
from geometry import is_inside_target

def setup():
    """Setup module variables.

    Parameters:
        None

    Returns:
        None    
    """
    global EMIN

    EMIN = 5.0  # eV


def trajectory(pos_init, dir_init, e_init):
    """Simulate one trajectory.
    
    Parameters:
        pos_init (ndarray): initial position of the projectile (size 3)
        dir_init (ndarray): initial direction of the projectile (size 3)
        e_init (float): initial energy of the projectile (eV)

    Returns:
        ndarray: final position of the projectile (size 3)
        ndarray: final direction of the projectile (size 3)
        float: final energy of the projectile (eV)
        bool: True if projectile is stopped inside the target, 
            False otherwise
    """
    pos = pos_init.copy()
    dir = dir_init.copy()
    e = e_init
    is_inside = True

    while e > EMIN:
        free_path, p, dirp, _ = get_recoil_position(pos, dir)
        dee = eloss(e, free_path)
        e -= dee
        pos += free_path * dir
        if not is_inside_target(pos):
            is_inside = False
            break
        dir, e, _, _ = scatter(e, dir, p, dirp)

    return pos, dir, e, is_inside