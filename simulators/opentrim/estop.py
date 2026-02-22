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


def setup(input_params, elements_params):
    """Setup module variables for electronic stopping.

    Parameters:
        input_params (dict): input parameters
        elements_params (ELEMENTS_PARAMS_DTYPE): element parameters

    Returns:
        (ESTOP_PARAMS_DTYPE): estop parameters
    """
    model = input_params["models"]["electronic stopping"]
    if model != "Lindhard":
        raise ValueError(f"Unsupported electronic stopping model: {model}")

    # correction factors to Lindhard stopping power
    nelem = len(elements_params)
    corr_lindhard = np.ones((nelem, nelem))  # default = 1.0
    for corr in input_params["models"]["Lindhard correction"]:
        elem1, elem2 = corr.split("->")
        for ielem1 in range(nelem):
            for ielem2 in range(nelem):
                if (elements_params[ielem1]["symbol"].strip() == elem1 and 
                    elements_params[ielem2]["symbol"].strip() == elem2):
                    corr_lindhard[ielem1, ielem2] = (
                        input_params["models"]["Lindhard correction"][corr])

    # Prefactor for Lindhard stopping (sqrt(eV)*A^2)
    fac_lindhard = np.empty((nelem, nelem))
    for ielem1 in range(nelem):
        z1 = elements_params[ielem1]["Z"]
        m1 = elements_params[ielem1]["M"]
        for ielem2 in range(nelem):
            z2 = elements_params[ielem2]["Z"]
            fac_lindhard[ielem1, ielem2] = (
                corr_lindhard[ielem1, ielem2] * 1.212 * z1**(7/6) * z2 
                / ((z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1)) )

    ESTOP_PARAMS_DTYPE = np.dtype([
        ("fac_lindhard", np.float64, (nelem, nelem)),
    ], align=True)

    estop_params = np.recarray(1, dtype=ESTOP_PARAMS_DTYPE)[0]
    estop_params["fac_lindhard"] = fac_lindhard

    return estop_params


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
    ilayer = proj["ilayer"]
    mat = materials_params[ilayer]

    ielem1 = proj["ielem"]
    weighted_fac_lindhard = 0.0
    for i in range(mat.nelem):
        ielem2 = mat.ielem[i]
        atomic_fraction = mat.atomic_fraction[i]
        weighted_fac_lindhard += (
            atomic_fraction * estop_params.fac_lindhard[ielem1, ielem2])
    
    dee = weighted_fac_lindhard * sqrt(proj["e"]) * mat.density * free_path

    if dee > proj["e"]:
        dee = proj["e"]
    
    return dee