"""Calculate the electronic stopping power.

Currently, only the Lindhard model (Phys. Rev. 124, (1961) 128) with 
a correction factor is implemented.

Available functions:
    setup: setup module variables.
    eloss: calculate the electronic energy loss.
"""
from math import sqrt


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
        None    
    """
    global FAC_LINDHARD, DENSITY

    FAC_LINDHARD = (corr_lindhard1 * 1.212 * z1**(7/6) * z2 / (
        (z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1) ),
        corr_lindhard2 * 1.212 * z2**(7/6) * z2 / (
        (z2**(2/3) + z2**(2/3))**(3/2) * sqrt(m2) ))         # eV/A
    DENSITY = density


def eloss(proj, free_path):
    """Calculate the electronic energy loss over a given free path length.

    Parameters:
        proj (Projectile): state of the projectile before the free flight path
        free_path (float): free path length (A)

    Returns:
        (float): energy loss (eV)
    """
    dee = FAC_LINDHARD[proj.ispec] * DENSITY * sqrt(proj.e) * free_path
    if dee > proj.e:
        dee = proj.e

    return dee