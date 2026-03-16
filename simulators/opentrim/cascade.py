"""Simulate projectile trajectories.

Available functions:
    cascade: simulate one cascade.
"""
import os
from collections import namedtuple

import numpy as np
from numba import jit, typed
from .mytypes import PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .recoil import get_recoil_position
from .scatter import scatter
from .estop import eloss
from .target import set_layer_index, set_element_index, set_is_inside_target
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
    
    emin = params.cascade.emin
    ed = params.cascade.ed
    
    stat = statistics.Statistics(params.stat)

    # NOTE: Record arrays cannot be created within numba-jitted functions, so we 
    # use regular structured arrays and access fields by name (rather than
    # by attributes).

    # Fully simulated projectiles
    proj_lst = np.empty(1 if not follow_recoils else prealloc, dtype=PROJ_DTYPE)
    lst_tail = 0
    
    # Projectiles to be simulated
    stack = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    stack.append(initial_proj)
    
    # Recoils of the currently simulated projectile
    recoils = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    
    # A recoil to be added to the simulation list
    recoil_proj = np.empty(1, dtype=PROJ_DTYPE)[0]

    while len(stack) > 0:
        proj = stack[-1]
    
        while proj["e"] > emin:
            free_path, p, dirp = get_recoil_position(
                proj, recoil_proj, params.recoil)
            
            # step projectile forward and update energy
            dee = eloss(proj, free_path, params.estop, params.materials)
            proj["e"] -= dee
            proj["pos"] += free_path * proj["dir"]
            set_layer_index(proj, params.geometry)
            set_is_inside_target(proj, params.geometry)
            
            if not proj["is_inside"]:
                break
            
            # get chemical element of recoil
            set_layer_index(recoil_proj, params.geometry)
            set_element_index(recoil_proj, params.materials)
            
            # scattering event
            scatter(proj, p, dirp, recoil_proj, params.scatter)
            set_is_inside_target(recoil_proj, params.geometry)
            # TODO: We may want to score the recoil energy even when it is 
            # below ed
            if follow_recoils and recoil_proj["e"] > ed:
                recoils.append(recoil_proj)    # Proj copied to the list (not a reference)
        
        if lst_tail >= proj_lst.size:
            proj_lst = np.append(
                proj_lst, np.empty(int((GROWTH_FACTOR - 1.0) * proj_lst.size), 
                                   dtype=PROJ_DTYPE))
        proj_lst[lst_tail] = proj   # TODO copied?
        lst_tail += 1
        
        stat.score(proj)
        stack.pop() # Remove currently processed projectile
        
        for i in range(len(recoils) - 1, -1, -1):
            stack.append(recoils[i])
            recoils.pop()

    # Return continuous arrays
    hist_results, mom_results = stat.results
    return proj_lst[:lst_tail][::-1].copy(), hist_results, mom_results.copy()
