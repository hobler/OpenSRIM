"""Calculate the electronic stopping power.

Currently, only the Lindhard model (Phys. Rev. 124, (1961) 128) with 
a correction factor is implemented.

Available functions:
    eloss: calculate the electronic energy loss.
"""
from math import sqrt
from collections import namedtuple
import numpy as np
from numba import jit


@jit
def estop_lindhard(e, fac_lindhard):
    """Calculate the Lindhard electronic stopping power for a given energy.

    Parameters:
        e (float): projectile energy (eV)
        fac_lindhard (float): pre-factor for the Lindhard stopping power
    
    Returns:
        (float): electronic stopping cross section (eV*A^2)
    """
    return fac_lindhard * sqrt(e)


@jit(inline = "always")
def eloss(proj, free_path, estop_params, materials_params):
    """Calculate the electronic energy loss over a given free path length.

    Parameters:
        proj (Projectile): state of the projectile before the free flight path
        free_path (float): free path length (A)
        estop_params (ESTOP_PARAMS_DTYPE): Electronic stopping parameters
        materials_params (MATERIALS_PARAMS_DTYPE): Material parameters

    Returns:
        (float): energy loss (eV)
    """
    e = proj["e"]
    imat = proj["ilayer"]
    ielem1 = proj["ielem"]

    weighted_se = 0.0
    for i in range(materials_params.nelem[imat]):
        ielem2 = materials_params.ielem[imat, i]
        if estop_params.model == "Lindhard":
            se = estop_lindhard(e, estop_params.fac_lindhard[ielem1, ielem2])
        else:
            se = np.interp(e, estop_params.srim_energies, 
                           estop_params.srim_table[ielem1, ielem2])
        atomic_fraction = materials_params.atomic_fraction[imat, i]
        weighted_se += atomic_fraction * se
    
    dee = weighted_se * materials_params.density[imat] * free_path

    if dee > proj["e"]:
        dee = proj["e"]
    
    return dee