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


def setup(input_params):
    """Setup module variables for electronic stopping.

    Parameters:
        input_params (dict): input parameters

    Returns:
        (ESTOP_PARAMS_DTYPE): estop parameters
    """
    corr_lindhard1 = input_params["models"]["Lindhard correction"]["B->Si"]
    corr_lindhard2 = input_params["models"]["Lindhard correction"]["Si->Si"]
    z1 = input_params["beam"]["Z"]
    m1 = input_params["beam"]["M"]
    z2 = input_params["layers"]["material"][0]["Z"][0]
    m2 = input_params["layers"]["material"][0]["M"][0]
    density = input_params["layers"]["density"][0]

    fac_lindhard = np.array([corr_lindhard1 * 1.212 * z1**(7/6) * z2 / (
        (z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1) ),
        corr_lindhard2 * 1.212 * z2**(7/6) * z2 / (
        (z2**(2/3) + z2**(2/3))**(3/2) * sqrt(m2) )])         # eV/A

    ESTOP_PARAMS_DTYPE = np.dtype([
        ("fac_lindhard", np.float64, (2,)),
        ("density", np.float64),
    ], align=True)

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