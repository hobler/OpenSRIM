"""Calculate the electronic stopping power.

Currently, only the Lindhard model (Phys. Rev. 124, (1961) 128) with 
a correction factor is implemented.

Available functions:
    eloss: calculate the electronic energy loss.
"""
from math import dist, sqrt
import numpy as np
from numba import jit
from . import config
from .target import get_distance_from_surface


@jit(debug=config.DEBUG)
def estop_lindhard(e, fac_lindhard):
    """Calculate the Lindhard electronic stopping power for a given energy.

    Parameters:
        e (float): projectile energy (eV)
        fac_lindhard (float): pre-factor for the Lindhard stopping power
    
    Returns:
        (float): electronic stopping cross section (eV*A^2)
    """
    return fac_lindhard * sqrt(e)


@jit(inline = "always", debug=config.DEBUG)
def eloss(proj, free_path, params):
    """Calculate the electronic energy loss over a given free path length.

    No electronic stopping is applied outside the target.
    TODO: Consider different stopping powers when layer boundaries are crossed.

    Parameters:
        proj (Projectile): state of the projectile before the free flight path
        free_path (float): free path length (A)
        params (PARAMS_DTYPE): Simulation parameters

    Returns:
        (float): energy loss (eV)
    """
    e = proj["e"]
    imat = proj["ilayer"]
    ielem1 = proj["ielem"]

    # No electronic stopping outside the target
    dist_beg, beamside_beg = get_distance_from_surface(proj["pos"], params)
    pos_end = proj["pos"] + free_path * proj["dir"]
    dist_end, beamside_end = get_distance_from_surface(pos_end, params)

    if (dist_beg > 0.0) and (dist_end > 0.0):
        pass
    elif dist_beg > 0.0:
        free_path *= dist_beg / (dist_beg - dist_end)
    elif dist_end > 0.0:
        free_path *= dist_end / (dist_end - dist_beg)
    elif beamside_end == beamside_beg:
        return 0.0
    else:
        dist = abs(proj["pos"][0] - pos_end[0])
        free_path *= (dist + dist_beg + dist_end) / dist

    weighted_se = 0.0
    for i in range(params.materials.nelem[imat]):
        ielem2 = params.materials.ielem[imat, i]
        if params.estop.model == "Lindhard":
            se = estop_lindhard(e, params.estop.fac_lindhard[ielem1, ielem2])
        else:
            se = np.interp(e, params.estop.srim_energies, 
                           params.estop.srim_table[ielem1, ielem2])
        atomic_fraction = params.materials.atomic_fraction[imat, i]
        weighted_se += atomic_fraction * se
    
    dee = weighted_se * params.materials.density[imat] * free_path

    if dee > proj["e"]:
        dee = proj["e"]
    
    return dee