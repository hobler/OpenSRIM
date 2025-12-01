"""Simulate projectile trajectories.

The trajectory function may call itself recursively to follow recoil
trajectories.

Available functions:
    setup: setup module variables.
    trajectory: simulate one trajectory.
"""
from select_recoil import get_recoil_position
from scatter import scatter
from estop import eloss
from geometry import is_inside_target
from mytypes import Projectile


def setup():
    """Setup module variables."""
    global EMIN, ED

    EMIN = 5.0  # eV
    ED = 15.0   # eV


def trajectory(proj, follow_recoils=False):
    """Simulate one projectile trajectory.
    
    Parameters:
        proj: (Projectile) the initial state of the projectile
        follow_recoils: (bool) whether to follow recoil trajectories
        
    Returns:
        (list[Projectile]) list of final projectile states
        (bool) whether the projectile stopped inside the target
    """
    proj_lst = []

    while proj.e > EMIN:
        free_path, p, dirp, recoil_pos = get_recoil_position(proj.pos[:], 
                                                             proj.dir[:])
        dee = eloss(proj, free_path)
        proj.e -= dee
        proj.pos += free_path * proj.dir[:]
        if not is_inside_target(proj.pos[:]):
            proj.is_inside = False
            break
        proj, recoil_dir, recoil_e = scatter(proj, p, dirp[:])
        if follow_recoils and recoil_e > ED:
            recoil = Projectile(
                e = recoil_e,
                pos = recoil_pos[:],
                dir = recoil_dir[:],
                ispec = 1               # assuming recoil is of type 1
            )
            recoil_lst = trajectory(recoil)
            proj_lst.extend(recoil_lst)

    proj_lst.append(proj)

    return proj_lst