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
    imat = proj["ilayer"]
    ielem1 = proj["ielem"]

    weighted_fac_lindhard = 0.0
    for i in range(materials_params.nelem[imat]):
        ielem2 = materials_params.ielem[imat, i]
        atomic_fraction = materials_params.atomic_fraction[imat, i]
        weighted_fac_lindhard += (
            atomic_fraction * estop_params.fac_lindhard[ielem1, ielem2])
    
    dee = (weighted_fac_lindhard * sqrt(proj["e"]) 
           * materials_params.density[imat] * free_path)

    if dee > proj["e"]:
        dee = proj["e"]
    
    return dee