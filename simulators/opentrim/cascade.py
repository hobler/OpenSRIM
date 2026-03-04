"""Simulate projectile trajectories.

Available functions:
    cascade: simulate one cascade.
"""
import os
from collections import namedtuple

import numpy as np
from numba import jit
from .mytypes import PROJ_DTYPE
from .recoil import get_recoil_position
from .scatter import scatter
from .estop import eloss
from .target import get_layer_index, get_element_index, is_inside_target
from . import stats as statistics


# TODO: Make follow_recoils an input parameter, passed via cascade_params
@jit
def cascade(initial_proj, params, follow_recoils=False, prealloc=400):
    """Simulate one projectile trajectory.
    
    Parameters:
        initial_proj: (Projectile) the initial state of the first projectile
        params: (PARAMS_DTYPE) Simulation parameters
        follow_recoils: (bool) whether to follow recoil trajectories
        prealloc: (int) number of recoil projectiles to pre-allocate space for (for better performance)
        
    Returns:
        tuple[ndarray[Projectile], ndarray[int32], ndarray[float64]]:
            list of final projectile states,
            results buffer of Histogram_1d class,
            results buffer of Moment_1d class
    """
    GROWTH_FACTOR = 1.5
    INITIAL_STACK_SIZE = 100
    INITIAL_RECOILS_SIZE = 5
    
    emin = params.cascade.emin
    ed = params.cascade.ed
    #emin = 5.0
    #ed = 15.0
    
    stat = statistics.Statistics(params.stat)

    # NOTE: We cannot create record arrays within numba-jitted functions, so we 
    # use regular structured arrays and access fields by name (rather than
    # by attributes).

    # Fully simulated projectiles
    proj_lst = np.empty(1 if not follow_recoils else prealloc, dtype=PROJ_DTYPE)
    lst_tail = 0
    
    # Projectiles to be simulated
    stack = np.empty(INITIAL_STACK_SIZE, dtype=PROJ_DTYPE)
    stack[0] = initial_proj
    stack_tail = 1
    
    # Recoils of the currently simulated projectile
    recoils = np.empty(INITIAL_RECOILS_SIZE, dtype=PROJ_DTYPE)

    while stack_tail > 0:
        stack_tail -= 1
        proj = stack[stack_tail]
        recoils_tail = 0
    
        while proj["e"] > emin:
            free_path, p, dirp, recoil_pos = get_recoil_position(
                proj, params.recoil)
            #print(f"Trajectory start: proj_e={proj['e']:.2f} eV")
            
            # step projectile forward and update energy
            dee = eloss(proj, free_path, params.estop, params.materials)
            proj["e"] -= dee
            proj["pos"] += free_path * proj["dir"]
            proj["ilayer"] = get_layer_index(proj["pos"], params.geometry)
            proj["is_inside"] = is_inside_target(proj["pos"], params.geometry)
            
            #print(f"geometry_params={params.geometry}")
            #print(f"proj_e={proj['e']:.2f} eV, proj_pos={proj['pos']}, is_inside={proj['is_inside']}")

            if not proj["is_inside"]:
                break
            
            # get chemical element of recoil
            recoil_ilayer = get_layer_index(recoil_pos, 
                                            params.geometry)
            recoil_ielem = get_element_index(recoil_ilayer, 
                                             params.materials)
            
            # scattering event
            recoil_dir, recoil_e = scatter(proj, p, dirp, recoil_ielem,
                                           params.scatter)
            # Create recoil projectile and add to recoils list
            # TODO: We may want to score the recoil energy even when it is 
            # below ed
            if follow_recoils and recoil_e > ed:
                recoil_is_inside = is_inside_target(recoil_pos, 
                                                    params.geometry)

                if recoils_tail >= recoils.size:
                    #print("Growing recoils array from size", recoils.size)
                    recoils = np.append(
                        recoils, 
                        np.empty(int((GROWTH_FACTOR - 1.0)* recoils.size), 
                                 dtype=PROJ_DTYPE))
                    #print("to size", recoils.size)
                recoils[recoils_tail]["e"] = recoil_e
                recoils[recoils_tail]["pos"] = recoil_pos
                recoils[recoils_tail]["dir"] = recoil_dir
                recoils[recoils_tail]["ielem"] = recoil_ielem
                recoils[recoils_tail]["ilayer"] = recoil_ilayer
                recoils[recoils_tail]["is_inside"] = recoil_is_inside
                recoils_tail += 1
        
        if lst_tail >= proj_lst.size:
            proj_lst = np.append(
                proj_lst, np.empty(int((GROWTH_FACTOR - 1.0) * proj_lst.size), 
                                   dtype=PROJ_DTYPE))
        proj_lst[lst_tail] = proj
        lst_tail += 1
        
        stat.score(proj)
        
        for i in range(recoils_tail - 1, -1, -1):
            if stack_tail >= stack.size:
                stack = np.append(stack, np.empty(int(GROWTH_FACTOR * stack.size), dtype=PROJ_DTYPE))
            stack[stack_tail] = recoils[i]
            stack_tail += 1

    # Return continuous arrays
    hist_results, mom_results = stat.results
    return proj_lst[:lst_tail][::-1].copy(), hist_results, mom_results.copy()
