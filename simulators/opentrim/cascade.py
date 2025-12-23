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
import numpy as np
from numba import jit


def setup():
    """Setup module variables."""
    global EMIN, ED

    EMIN = 5.0  # eV
    ED = 15.0   # eV


@jit(fastmath=True)
def trajectory(initial_proj, follow_recoils=False, prealloc=100):
    """Simulate one projectile trajectory.
    
    Parameters:
        initial_proj: (Projectile) the initial state of the first projectile
        follow_recoils: (bool) whether to follow recoil trajectories
        prealloc: (int) number of recoil projectiles to pre-allocate space for (for better performance)
        
    Returns:
        (numpy.ndarray[Projectile]) list of final projectile states
    """
    proj_lst = np.full(1 if not follow_recoils else prealloc, initial_proj)
    tail = 0
    head = 1

    while tail < head:
        proj = proj_lst[tail]
        while proj.e > EMIN:
            free_path, p, dirp, recoil_pos = get_recoil_position(proj.pos[:], proj.dir[:])
            
            dee = eloss(proj, free_path)
            proj.e -= dee
            proj.pos += free_path * proj.dir[:]
            
            if not is_inside_target(proj.pos[:]):
                proj.is_inside = False
                break
            
            recoil_dir, recoil_e = scatter(proj, p, dirp[:])        
            if follow_recoils and recoil_e > ED:
                if head == proj_lst.size:
                    proj_lst = np.append(proj_lst, np.full(int(1.5 * proj_lst.size), initial_proj))
                
                proj_head = proj_lst[head]
                proj_head.e = recoil_e
                proj_head.pos[:] = recoil_pos
                proj_head.dir[:] = recoil_dir
                proj_head.ispec = 1
                proj_head.is_inside = True
                head += 1
        tail+=1

    return proj_lst[:head]