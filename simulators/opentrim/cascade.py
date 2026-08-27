"""Simulate projectile trajectories.

Available functions:
    cascade: simulate one cascade.
"""
import os
NUMBA_DISABLE_JIT = os.environ.get("NUMBA_DISABLE_JIT", "") == "1"
from copy import deepcopy

import numpy as np
from numba import jit, typed
from .mytypes import PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .recoil import select_recoil
from .scatter import scatter
from .estop import eloss
from .target import check_exit_and_move
from .stats import (score_eed, score_ned, score_start, score_stop, 
                    score_backscattered, score_transmitted, 
                    score_yield_back, score_yield_in, score_yield_trans)
from . import config


@jit(debug=config.DEBUG)
def _check_replacement_collision(proj, recoil, params):
    """Check if the recoil is a replacement collision and modify the projectile 
    accordingly.

    Parameters:
        proj (Projectile): state of the projectile (modified in-place)
        recoil (Projectile): state of the recoil (modified in-place)
        params (PARAMS_DTYPE): Simulation parameters
    """
    if proj["ielem"] != recoil["ielem"]:
        return

    imat = proj["ilayer"]
    ielem_mat = params.materials[imat].ielem_mat[recoil["ielem"]]
    e_disp = params.materials[imat].edisp[ielem_mat]

    if proj["e"] < e_disp and recoil["e"] > e_disp and recoil["e"] > proj["e"]:
        proj, recoil = recoil, proj


@jit(debug=config.DEBUG)
def _get_damage_kp(recoil, params):
    """Get the damage produced by a recoil using the Kinchin-Pease model.

    Parameters:
        recoil (Projectile): state of the recoil
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (float): energy deposited by the recoil into nuclear collisions (NED)
        (float): number of displacements produced by the recoil
            """
    kd_KP = params.elements[recoil["ielem"]].kd_KP
    fd_KP = params.elements[recoil["ielem"]].fd_KP

    ed = fd_KP * recoil["e"]
    ned = (recoil["e"] 
           / (1.0 + kd_KP * (ed + 0.40244*ed**(3/4) + 3.4008*ed**(1/6))))

    imat = recoil["ilayer"]
    ielem_mat = params.materials[imat].ielem_mat[recoil["ielem"]]
    e_disp = params.materials[imat].edisp[ielem_mat]
    if ned < e_disp:
        nvac = 0.0
    elif ned < 2.5 * e_disp:
        nvac = 1.0
    else:
        nvac = ned / (2.5 * e_disp)

    return ned, nvac


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
    #print("Entered cascade...")
    # In Python mode, we need to create a deep copy of the initial projectile 
    # to avoid modifying it in-place when it is appended to the projectile 
    # stack. In Numba mode, this is not necessary, since the initial projectile 
    # is passed by value to the cascade function.
    if NUMBA_DISABLE_JIT:
        initial_proj = deepcopy(initial_proj)

    emin = params.cascade.emin
    
    # Fully simulated projectiles
    final_proj_lst = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    
    # Projectiles to be simulated
    proj_stack = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_stack.append(initial_proj)
    
    # Due to a limitation of Numba, we cannot create a recoil structured array
    # in select_recoil and return it from there, so we create it here and 
    # modify it in-place in select_recoil
    recoil = np.empty(1, dtype=PROJ_DTYPE)[0]

    # Counters for yields
    nin = np.zeros(stats["i"]["nvar"], dtype=np.float64)
    nback = np.zeros(stats["b"]["nvar"], dtype=np.float64)
    ntrans = np.zeros(stats["t"]["nvar"], dtype=np.float64)

    # TODO: Apply the surface binding energy to the projectile when it enters 
    # the target. For that, we need a way to specify the surface binding energy.

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
            if proj["dir"][0] < 0:
                score_backscattered(stats, proj)
                nback[proj["ielem"]] += 1
            else:
                score_transmitted(stats, proj)
                ntrans[proj["ielem"]] += 1
            proj_stack.pop()
            continue

        # terminate trajectory if the projectile has no more energy
        if proj["e"] <= emin:
            final_proj_lst.append(proj)
            score_ned(stats, proj)
            score_stop(stats, proj)
            nin[proj["ielem"]] += 1
            proj_stack.pop()
            continue

        # treat scattering event and recoil
        if recoil["is_inside"]:
            # scattering event, modifying proj and recoil in-place
            scatter(proj, p, dirp, recoil, params)
            if params.cascade.replacement_collisions:
                _check_replacement_collision(proj, recoil, params)

            # terminate trajectory if the projectile has lost too much energy
            if proj["e"] <= emin:
                final_proj_lst.append(proj)
                score_ned(stats, proj)
                score_stop(stats, proj)
                nin[proj["ielem"]] += 1
                proj_stack.pop()

            # start a new sub-cascade if the recoil has enough energy to leave 
            # its position
            imat = recoil["ilayer"]
            ielem_mat = params.materials[imat].ielem_mat[recoil["ielem"]]
            if recoil["dist_surf"] > 10.0:
                e_disp = params.materials[imat].edisp[ielem_mat]
            else:  # if close to the surface, use surface binding energy
                e_disp = params.materials[imat].esurf[ielem_mat]
            e_bulk = params.materials[imat].ebulk[ielem_mat]
            if params.cascade.follow_recoils:
                if recoil["e"] > e_disp:
                    recoil["e"] -= e_bulk
                    proj_stack.append(recoil)  # stores a copy of recoil
                    score_start(stats, recoil, params.nelem_target)
                    nin[recoil["ielem"] + params.nelem_target] += 1
                    recoil["e"] = e_bulk  # score the binding energy as NED
                    score_ned(stats, recoil)
                else:
                    score_ned(stats, recoil)
            else:
                ned, nvac = _get_damage_kp(recoil, params)
                eed = recoil["e"] - ned
                score_stop(stats, recoil, weight=nvac)
                nin[recoil["ielem"]] += nvac
                recoil["e"] = ned
                score_ned(stats, recoil)
                recoil["e"] = eed
                score_eed(stats, recoil, dee=eed)

    # Score the yields
    score_yield_in(stats, nin)
    score_yield_back(stats, nback)
    score_yield_trans(stats, ntrans)                

    # Return fully simulated projectiles in the correct order
    return final_proj_lst[::-1]
