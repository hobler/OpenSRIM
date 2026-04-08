"""Simulate projectile trajectories.

Available functions:
    cascade: simulate one cascade.
"""
import os
from collections import namedtuple

import numpy as np
from numba import jit, typed
from .mytypes import PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .recoil import select_recoil
from .scatter import scatter
from .estop import eloss
from .target import get_layer_index, is_inside_target
from .stats import score_eed, score_ned, score_start, score_end


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
    emin = params.cascade.emin
    ed = params.cascade.ed
    
    # Fully simulated projectiles
    proj_lst = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    
    # Projectiles to be simulated
    proj_stack = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_stack.append(initial_proj)
    
    # Due to a limitation of Numba, we cannot create a recoil structured array
    # in select_recoil and return it from there, so we create it here and 
    # modify it in select_recoil in-place
    recoil = np.empty(1, dtype=PROJ_DTYPE)[0]

    # Loop over collision events until there are no more projectiles to 
    # simulate
    while len(proj_stack) > 0:
        proj = proj_stack[-1]
    
        # set recoil parameters (recoil modified in-place)
        free_path, p, dirp = select_recoil(proj, recoil, params)

        # step projectile forward considering electronic energy loss
        dee = eloss(proj, free_path, params.estop, params.materials)
        proj["e"] -= dee
        proj["pos"] += free_path * proj["dir"]
        proj["ilayer"] = get_layer_index(proj["pos"], params.geometry)
        proj["is_inside"] = is_inside_target(proj["pos"], params.geometry)
        score_eed(stats, proj, dee) 

        # terminate trajectory if the projectile is outside the target, or if 
        # it has no more energy
        if not proj["is_inside"] or proj["e"] <= emin:
            proj_lst.append(proj)
            score_end(stats, proj)
            proj_stack.pop()
            continue

        # scattering event, modifying proj and recoil in-place
        scatter(proj, p, dirp, recoil, params.scatter)

        # terminate trajectory if the projectile has no more energy
        if proj["e"] <= emin:
            proj_lst.append(proj)
            score_end(stats, proj)
            proj_stack.pop()

        # start a new cascade if the recoil has enough energy to leave its 
        # position
        if True:
            imat = recoil["ilayer"]
            ielem = params.materials[imat].ielem[recoil["ielem"]]
            ed = params.materials[imat].displacement_energy[ielem]
        if (params.cascade.follow_recoils and recoil["e"] > ed):
            proj_stack.append(recoil)
            score_start(stats, recoil, params.nelem_target)
        else:
            score_ned(stats, recoil)

    # Return fully simulated projectiles in the correct order
    return proj_lst[::-1]
