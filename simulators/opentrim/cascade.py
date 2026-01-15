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
    """Setup module variables.
    
    Returns:
        (float): EMIN
        (float): ED
    """
    emin = 5.0  # eV
    ed = 15.0   # eV
    return emin, ed


@jit(fastmath=True)
def trajectory(initial_proj, sim_params_arr, follow_recoils=False, prealloc=100):
    """Simulate one projectile trajectory.
    
    Parameters:
        initial_proj: (Projectile) the initial state of the first projectile
        sim_params: (SimParams) Simulation parameters
        follow_recoils: (bool) whether to follow recoil trajectories
        prealloc: (int) number of recoil projectiles to pre-allocate space for (for better performance)
        
    Returns:
        (numpy.ndarray[Projectile]) list of final projectile states
    """
    sim_params = sim_params_arr[0]
    proj_lst = np.full(1 if not follow_recoils else prealloc, initial_proj)
    tail = 0
    head = 1

    while tail < head:
        proj = proj_lst[tail]
        # scatter_params = sim_params.scatter_params # SimParams is now flattened
        while proj.e > sim_params.cascade_params.emin:
            free_path, p, dirp, recoil_pos = get_recoil_position(proj.pos[:], proj.dir[:], sim_params.recoil_params)
            
            dee = eloss(proj, free_path, sim_params.estop_params)
            proj.e -= dee
            proj.pos += free_path * proj.dir[:]
            
            if not is_inside_target(proj.pos[:], sim_params.geometry_params):
                proj.is_inside = False
                break
            
            recoil_dir, recoil_e = scatter(proj, p, dirp[:], sim_params.scatter_params)        
            if follow_recoils and recoil_e > sim_params.cascade_params.ed:
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