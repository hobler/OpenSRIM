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
from .target import check_exit_and_move
from .stats import score_eed, score_ned, score_start, score_stop, score_exit
from . import config


@jit(debug=config.DEBUG)
def cascade(initial_proj, params, stats):
    """Simulate trajectories of one ion and its recoils.
    
    Parameters:
        initial_proj: (Projectile) the initial state of the first projectile
        params: (PARAMS_DTYPE) Simulation parameters
        stats: (STATS_DTYPE) statistics data container
        
    Returns:
        ndarray[Projectile]: list of final projectile states
    """
    emin = params.cascade.emin
    ed = params.cascade.ed  # ignored (see below)
    
    # Fully simulated projectiles
    final_proj_lst = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    
    # Projectiles to be simulated
    proj_stack = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_stack.append(initial_proj)
    
    # Due to a limitation of Numba, we cannot create a recoil structured array
    # in select_recoil and return it from there, so we create it here and 
    # modify it in-place in select_recoil
    recoil = np.empty(1, dtype=PROJ_DTYPE)[0]

    # Loop over collision events until there are no more projectiles to 
    # simulate
    while len(proj_stack) > 0:

        proj = proj_stack[-1]
    
        # set recoil parameters before interaction(recoil is modified in-place)
        free_path, p, dirp = select_recoil(proj, recoil, params)
        free_path += proj["dffp_new"]

        # consider electronic energy loss
        dee = eloss(proj, free_path + proj["dffp_old"], params)
        proj["e"] -= dee
        score_eed(stats, proj, dee) 

        # move projectile forward and terminate trajectory if projectile exits 
        # the target (proj is modified in-place)
        # note that an exiting projectile cannot recoil a target atom
        exiting = check_exit_and_move(proj, free_path, params)
        if exiting:
            final_proj_lst.append(proj)
            score_exit(stats, proj)
            proj_stack.pop()
            continue

        # terminate trajectory if the projectile has no more energy
        if proj["e"] <= emin:
            final_proj_lst.append(proj)
            score_stop(stats, proj)
            proj_stack.pop()
            continue

        # treat scattering event and recoil
        if recoil["is_inside"]:
            # scattering event, modifying proj and recoil in-place
            #print("got to scatter")
            scatter(proj, p, dirp, recoil, params)
            #print("got past scatter")

            # terminate trajectory if the projectile has lost too much energy
            if proj["e"] <= emin:
                final_proj_lst.append(proj)
                score_stop(stats, proj)
                proj_stack.pop()

            # start a new cascade if the recoil has enough energy to leave its 
            # position
            if True:
                imat = recoil["ilayer"]
                ielem = params.materials[imat].ielem[recoil["ielem"]]
                ed = params.materials[imat].edisp[ielem]
                eb = params.materials[imat].ebulk[ielem]
            if (params.cascade.follow_recoils and recoil["e"] > ed):
                recoil["e"] -= eb
                proj_stack.append(recoil)
                score_start(stats, recoil, params.nelem_target)
            else:
                score_ned(stats, recoil)

    # Return fully simulated projectiles in the correct order
    return final_proj_lst[::-1]
