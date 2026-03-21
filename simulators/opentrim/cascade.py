"""Simulate projectile trajectories.

Available functions:
    cascade: simulate one cascade.
"""
import os
from collections import namedtuple

import numpy as np
from numba import jit, typed
from .mytypes import PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .recoil import get_recoil
from .scatter import scatter
from .estop import eloss
from .target import get_layer_index, is_inside_target
from .stats import score


@jit
def cascade(initial_proj, params, stats):
    """Simulate one projectile trajectory.
    
    Parameters:
        initial_proj: (Projectile) the initial state of the first projectile
        params: (PARAMS_DTYPE) Simulation parameters
        stats: (STATS_DTYPE) statistical data container
        
    Returns:
        ndarray[Projectile]: list of final projectile states
    """
    GROWTH_FACTOR = 1.5
    
    emin = params.cascade.emin
    ed = params.cascade.ed
    
    # Fully simulated projectiles
    proj_lst = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    
    # Projectiles to be simulated
    proj_stack = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_stack.append(initial_proj)
    
    # Recoils of the currently simulated projectile
    recoils = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    
    # A recoil to be added to the simulation list
    recoil = np.empty(1, dtype=PROJ_DTYPE)[0]

    while len(proj_stack) > 0:
        proj = proj_stack[-1]
    
        while proj["e"] > emin:
            # set recoil parameters (modified in-place)
            free_path, p, dirp = get_recoil(proj, recoil, params)
            
            # step projectile forward and update energy
            dee = eloss(proj, free_path, params.estop, params.materials)
            proj["e"] -= dee
            proj["pos"] += free_path * proj["dir"]
            proj["ilayer"] = get_layer_index(proj["pos"], params.geometry)
            proj["is_inside"] = is_inside_target(proj["pos"], params.geometry)
            
            if not proj["is_inside"]:
                break
            
            # scattering event, modifying proj and recoil in-place
            scatter(proj, p, dirp, recoil, params.scatter)
            # TODO: We may want to score the recoil energy even when it is 
            # below ed
            if params.cascade.follow_recoils and recoil["e"] > ed:
                recoils.append(recoil)    # Proj copied to the list (not a reference)
        
        proj_lst.append(proj)
        score(stats, proj)

        # Remove currently processed projectile
        proj_stack.pop()
        
        for i in range(len(recoils) - 1, -1, -1):
            proj_stack.append(recoils[i])
            recoils.pop()

    # Return fully simulated projectiles in the correct order
    return proj_lst[::-1]
