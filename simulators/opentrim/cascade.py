"""Simulate projectile trajectories.

Available functions:
    setup: setup module variables.
    cascade: simulate one cascade.
"""
import numpy as np
from numba import jit
from .select_recoil import get_recoil_position
from .scatter import scatter
from .estop import eloss
from .geometry import is_inside_target
from . import stats as statistics


def setup():
    """Setup module variables.
    
    Returns:
        (CASCADE_PARAMS_DTYPE): The cascade parameters
    """

    CASCADE_PARAMS_DTYPE = np.dtype([
        ("emin", np.float64),
        ("ed", np.float64),
    ], align=True)

    cascade_params = np.recarray(1, dtype=CASCADE_PARAMS_DTYPE)[0]

    cascade_params.emin = 5.0  # eV
    cascade_params.ed = 15.0   # eV
    
    return cascade_params


@jit
def cascade(initial_proj, params_arr, screen_fun, follow_recoils=False, prealloc=400):
    """Simulate one projectile trajectory.
    
    Parameters:
        initial_proj: (Projectile) the initial state of the first projectile
        params: (PARAMS_DTYPE) Simulation parameters
        screen_fun (object): Screening function
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
    
    params = params_arr[0]
    emin = params.cascade.emin
    ed = params.cascade.ed
    is_magic = (params.scatter.pot_model == "ZBL_magic")
    stat = params.stat
    hist = statistics.Histogram_1d(stat.nspec, stat.nbin, (stat.limits[0], stat.limits[1]))
    mom = statistics.Moment_1d(stat.nspec, 4)

    #proj_lst = np.empty(1 if not follow_recoils else prealloc, dtype=initial_proj.dtype)
    #proj_lst[0] = initial_proj

    # proj_lst[1:] will be overwritten with recoils
    proj_lst = np.full(1 if not follow_recoils else prealloc, initial_proj)
    lst_tail = 0
    
    stack = np.full(INITIAL_STACK_SIZE, initial_proj)
    stack_tail = 1
    
    recoils = np.full(INITIAL_RECOILS_SIZE, initial_proj)

    while stack_tail > 0:
        stack_tail -= 1
        proj = stack[stack_tail]
        recoils_tail = 0
        
        while proj.e > emin:
            free_path, p, dirp, recoil_pos = get_recoil_position(proj.pos[:], proj.dir[:], params.recoil)
            
            dee = eloss(proj, free_path, params.estop)
            proj.e -= dee
            proj.pos += free_path * proj.dir[:]
            
            if not is_inside_target(proj.pos[:], params.geometry):
                proj.is_inside = False
                break
            
            recoil_dir, recoil_e = scatter(proj, p, dirp[:], screen_fun, params.scatter, is_magic)        
            if follow_recoils and recoil_e > ed:
                if recoils_tail >= recoils.size:
                    recoils = np.append(recoils, np.full(int(GROWTH_FACTOR * recoils.size), initial_proj))
                
                recoils[recoils_tail].e = recoil_e
                recoils[recoils_tail].pos[:] = recoil_pos
                recoils[recoils_tail].dir[:] = recoil_dir
                recoils[recoils_tail].ispec = 1
                recoils[recoils_tail].is_inside = True
                recoils_tail += 1
        
        if lst_tail >= proj_lst.size:
            proj_lst = np.append(proj_lst, np.full(int(GROWTH_FACTOR * proj_lst.size), initial_proj))
        proj_lst[lst_tail] = proj
        lst_tail+=1
        
        if proj.is_inside:
            hist.score(proj.ispec, proj.pos[2])
            mom.score(proj.ispec, proj.pos[2])
        
        for i in range(recoils_tail - 1, -1, -1):
            if stack_tail >= stack.size:
                stack = np.append(stack, np.full(int(GROWTH_FACTOR * stack.size), initial_proj))
            stack[stack_tail] = recoils[i]
            stack_tail += 1

    # Return continuous arrays
    return proj_lst[:lst_tail][::-1].copy(), hist.results.copy(), mom.results.copy()