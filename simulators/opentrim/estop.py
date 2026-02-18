"""Calculate the electronic stopping power.

Currently, only the Lindhard model (Phys. Rev. 124, (1961) 128) with 
a correction factor is implemented.

Available functions:
    setup: setup module variables.
    eloss: calculate the electronic energy loss.
"""
from math import sqrt
import numpy as np
from numba import jit
from .mytypes import ESTOP_PARAMS_DTYPE


def setup(corr_lindhard1, z1, m1, corr_lindhard2, z2, m2, density):
    """Setup module variables for electronic stopping.

    Parameters:
        corr_lindhard (float): Correction factor to Lindhard stopping power
        z1 (int): atomic number of projectile
        m1 (float): mass of projectile (amu)
        z2 (int): atomic number of the target atom
        m2 (float): mass of the target atom (amu)
        density (float): target density (atoms/A^3)

    Returns:
        (ESTOP_PARAMS_DTYPE): estop parameters
    """
    fac_lindhard = np.array([corr_lindhard1 * 1.212 * z1**(7/6) * z2 / (
        (z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1) ),
        corr_lindhard2 * 1.212 * z2**(7/6) * z2 / (
        (z2**(2/3) + z2**(2/3))**(3/2) * sqrt(m2) )])         # eV/A

    estop_params = np.recarray(1, dtype=ESTOP_PARAMS_DTYPE)[0]
    estop_params["fac_lindhard"] = fac_lindhard
    estop_params["density"] = density

    return estop_params


@jit(inline = "always")
def eloss(proj, free_path, params):
    """Calculate the electronic energy loss over a given free path length.

    Parameters:
        proj (Projectile): state of the projectile before the free flight path
        free_path (float): free path length (A)
        params (ESTOP_PARAMS_DTYPE): Estop parameters

    Returns:
        (float): energy loss (eV)
    """
    dee = params.fac_lindhard[proj.ispec] * params.density * sqrt(proj.e) * free_path
    if dee > proj.e:
        dee = proj.e

    return dee