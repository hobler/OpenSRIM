"""Simulate projectile trajectories.

Available functions:
    setup: setup module variables.
    trajectory: simulate one trajectory.
"""
from select_recoil import get_recoil_position
from scatter import scatter
from estop import eloss
from geometry import is_inside_target


def setup():
    """Setup module variables."""
    global EMIN

    EMIN = 5.0  # eV


def trajectory(proj):
    """Simulate one projectile trajectory.
    
    Parameters:
        proj: (Projectile) the initial state of the projectile

    Returns:
        (list[Projectile]) list of final projectile states
        (bool) whether the projectile stopped inside the target
    """
    is_inside = True
    proj_lst = []

    while proj.e > EMIN:
        free_path, p, dirp, _ = get_recoil_position(proj.pos, proj.dir)
        dee = eloss(proj.e, free_path)
        proj.e -= dee
        proj.pos += free_path * proj.dir
        if not is_inside_target(proj.pos):
            is_inside = False
            break
        proj.dir, proj.e, _, _ = scatter(proj.e, proj.dir, p, dirp)

    proj_lst.append(proj)

    return proj_lst, is_inside