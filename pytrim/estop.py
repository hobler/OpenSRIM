"""Calculate the electronic stopping power.

Currently, only the Lindhard model (Phys. Rev. 124, (1961) 128) with 
a correction factor is implemented.

Available functions:
    setup: setup module variables.
    eloss: calculate the electronic energy loss.
"""
from math import sqrt


def setup(corr_lindhard, z1, m1, z2, density):
    """Setup module variables depending on target density.

    Parameters:
        corr_lindhard (float): Correction factor to Lindhard stopping power
        z1 (int): atomic number of projectile
        m1 (float): mass of projectile (amu)
        z2 (int): atomic number of 
        density (float): target density (atoms/A^3)

    Returns:
        None    
    """
    global FAC_LINDHARD, DENSITY

    FAC_LINDHARD = corr_lindhard * 1.212 * z1**(7/6) * z2 / (
        (z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1) )
    DENSITY = density


def eloss(e, free_path):
    """Calculate the electronic energy loss over a given free path length.

    Parameters:
        e (float): energy of the projectile (eV)
        free_path (float): free path length (A)

    Returns:
        float: energy loss (eV)
    """
    dee = FAC_LINDHARD * DENSITY * sqrt(e) * free_path
    if dee > e:
        dee = e

    return dee