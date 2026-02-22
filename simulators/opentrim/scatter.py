"""Treat the scattering of a projectile on a target atom.

Currently, the ZBL potential (Ziegler, Biersack, Littmark, The Stopping and 
Range of Ions in Matter, Pergamon Press, 1985) and the NLHlin potential
(to be published) are implemented. The scattering integrals are evaluated
numerically, or Biersack"s "magic formula" is usedfor the scattering angle.

Available functions:
    setup: setup module variables.
    scatter: treat a scattering event.
"""
import os
import math
import numpy as np
from numba import jit
from .zbl import magic
from .cm_scatter import scatter_integrals


def _get_nlhlin_coefs(z1, z2):
    """Get the coefficients for the NLHlin screening function.

    Parameters:
        z1: (int) atomic number of atom 1
        z2: (int) atomic number of atom 2
        
    Returns:
        a1, a2, a3: prefactors
        b1, b2, b3: 1/screening lengths (1/A)
        rmax: maximum range of the potential (A)
    """
    fname = os.path.join(os.path.dirname(__file__), "dmol_coeffs_rmax.dat")
    if not os.path.exists(fname):
        raise OSError(f"get_nlhlin_coefs: Coefficients file {fname} not found")
    
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            items = line.split()
            if int(items[0]) == min(z1, z2) and int(items[1]) == max(z1, z2):
                return [float(c) for c in items[2:-1]]
    
    raise ValueError(f"get_nlhlin_coefs: Coefficients for Z1={z1}, Z2={z2} "
                     f"not found in {fname}")


def setup(input_params, elements_params):
    """Setup module variables depending on projectile and target species.

    Each of the module variables ENORM, RNORM, DIRFAC, and DENFAC is a tuple
    with two entries: one for the ion species 0 and one for moving atom 
    species 1. Currently we assume there is only on target atoms species.

    Parameters:
        input_params (dict): input parameters dictionary
        elements_params (ELEMENT_PARAMS_DTYPE): parameters of chemical elements 
            in the target
            
    Returns:
        (SCATTER_PARAMS_DTYPE): Scatter parameters
    """
    pot_model = input_params["models"]["potential"]
    if input_params["models"]["scattering integrals"]["algorithm"] == "magic":
        pot_model += "_magic"

    nelem = len(elements_params)
    rnorm = np.empty((nelem, nelem), dtype=np.float64)
    enorm = np.empty((nelem, nelem), dtype=np.float64)
    dirfac = np.empty((nelem, nelem), dtype=np.float64)
    denfac = np.empty((nelem, nelem), dtype=np.float64)

    if pot_model == "NLHlin":
        NLHLIN_COEFS_DTYPE = np.dtype([
            ("a1", np.float64),
            ("b1", np.float64),
            ("a2", np.float64),
            ("b2", np.float64),
            ("a3", np.float64),
            ("b3", np.float64),
            ("c", np.float64),
            ("d", np.float64),
            ("rmax", np.float64),
        ], align=True)
        nlhlin_coefs = np.empty((nelem, nelem), dtype=NLHLIN_COEFS_DTYPE)

    for ielem1 in range(nelem):
        for ielem2 in range(nelem):
            z1 = elements_params[ielem1]["Z"]
            z2 = elements_params[ielem2]["Z"]
            m1 = elements_params[ielem1]["M"]
            m2 = elements_params[ielem2]["M"]
            m1_m2 = m1 / m2

            if pot_model.startswith("ZBL"):
                rnorm[ielem1, ielem2] = 0.4685 / (z1**0.23 + z2**0.23)
            elif pot_model.startswith("NLHlin"):
                rnorm[ielem1, ielem2] = 0.4685 / (
                    math.sqrt(math.sqrt(z1) + math.sqrt(z2)))
            enorm[ielem1, ielem2] = (14.39979 * z1 * z2 / rnorm[ielem1, ielem2] 
                                     * (1 + m1_m2))
            dirfac[ielem1, ielem2] = 2 / (1 + m1_m2)
            denfac[ielem1, ielem2] = 4 * m1_m2 / (1 + m1_m2)**2
            
            if pot_model == "NLHlin":
                a1, b1, a2, b2, a3, b3, rmax = _get_nlhlin_coefs(z1, z2)
                b1 *= rnorm[ielem1, ielem2]
                b2 *= rnorm[ielem1, ielem2]
                b3 *= rnorm[ielem1, ielem2]
                c = 1.0 - (a1 + a2 + a3)
                rmax /= rnorm[ielem1, ielem2]
                d = (a1*np.exp(-b1*rmax) + a2*np.exp(-b2*rmax) 
                     + a3*np.exp(-b3*rmax) + c)
                nlhlin_coefs[ielem1, ielem2]["a1"] = a1
                nlhlin_coefs[ielem1, ielem2]["b1"] = b1
                nlhlin_coefs[ielem1, ielem2]["a2"] = a2
                nlhlin_coefs[ielem1, ielem2]["b2"] = b2
                nlhlin_coefs[ielem1, ielem2]["a3"] = a3
                nlhlin_coefs[ielem1, ielem2]["b3"] = b3
                nlhlin_coefs[ielem1, ielem2]["c"] = c
                nlhlin_coefs[ielem1, ielem2]["d"] = d
                nlhlin_coefs[ielem1, ielem2]["rmax"] = rmax 

    SCATTER_PARAMS_DTYPE = np.dtype([
        ("pot_model", "<U16"),
        ("enorm", np.float64, (nelem, nelem)),
        ("rnorm", np.float64, (nelem, nelem)),
        ("dirfac", np.float64, (nelem, nelem)),
        ("denfac", np.float64, (nelem, nelem)),
    ], align=True)
    if pot_model == "NLHlin":
        SCATTER_PARAMS_DTYPE = np.dtype(SCATTER_PARAMS_DTYPE.descr + [
            ("nlhlin_coefs", NLHLIN_COEFS_DTYPE, (nelem, nelem)),
        ], align=True)

    scatter_params = np.recarray(1, dtype=SCATTER_PARAMS_DTYPE)[0]
    scatter_params["pot_model"] = pot_model
    scatter_params["enorm"] = enorm
    scatter_params["rnorm"] = rnorm
    scatter_params["dirfac"] = dirfac
    scatter_params["denfac"] = denfac
    if pot_model == "NLHlin":
        scatter_params["nlhlin_coefs"] = nlhlin_coefs

    return scatter_params


@jit(inline = "always")
def normalize_if_needed(vec, fallback):
    """Fast normalization with fallback – in‑place, no new array.
    
    Parameters:
        vec (np.ndarray): Vector to be normalized (size 3)
        fallback (np.ndarray): Vector to replace the original with if 
            len(norm) == 0
        
    Returns:
        np.ndarray: The normalized vector (size 3)
    """
    norm_sq = vec[0]**2 + vec[1]**2 + vec[2]**2
    if norm_sq == 0.0:
        # fallback is already a unit vector, just copy it
        vec[0], vec[1], vec[2] = fallback[0], fallback[1], fallback[2]
        return vec
    norm = math.sqrt(norm_sq)
    vec[0] /= norm
    vec[1] /= norm
    vec[2] /= norm
    return vec


@jit
def scatter(proj, p, dirp, ielem2, scatter_params):
    """Treat a scattering event.

    The atomic numbers and masses of the ion and the target atom enter the
    calculation via the ENORM, PNORM, DIRFAC, and DENFAC parameters.

    The direction vectors proj["dir"] and dirp[:] are assumed to be normalized 
    to unit length.

    Parameters:
        proj (Projectile): state of the projectile (will be modified in-place)
        p (float): impact parameter (A)
        dirp (ndarray): direction vector of the impact parameter
            (= from the collision point to the recoil position before 
            the collision) (unit vector, size 3)
        ielem2 (int): index of the target element
        scatter_params (np.recarray): Scatter parameters
        is_magic (bool): If magic function should be used (otherwise 
            scatter_integrals)
    
    Returns:
        (ndarray): direction vector of the recoil after the collision 
            (size 3)
        (float): energy of the projectile after the collision
    """
    # scattering angle theta in the center-of-mass system
    ielem1 = proj["ielem"]
    proj_e = proj["e"]
    proj_dir = proj["dir"][:]

    enorm = scatter_params.enorm[ielem1, ielem2]
    rnorm = scatter_params.rnorm[ielem1, ielem2]
    dirfac = scatter_params.dirfac[ielem1, ielem2]
    denfac = scatter_params.denfac[ielem1, ielem2]
    pot_model = scatter_params.pot_model
    
    if pot_model == "ZBL_magic":
        cos_half_theta = magic(proj_e/enorm, p/rnorm)
        sin_half_theta = math.sqrt(1 - cos_half_theta**2)
    elif pot_model.endswith("magic"):
        raise ValueError(f"Unknown potential model {pot_model}")
    else:
        theta, _ = scatter_integrals(proj_e/enorm, p/rnorm, pot_model)
        sin_half_theta = math.sin(0.5 * theta)
        cos_half_theta = math.cos(0.5 * theta)

    # directions of the recoil and the projectile after the collision
    recoil_dir = dirfac * sin_half_theta * (
        sin_half_theta*proj_dir[:] + cos_half_theta*dirp[:])
    dir_new = proj_dir[:] - recoil_dir[:]
    dir_new = normalize_if_needed(dir_new, proj_dir[:])
    recoil_dir = normalize_if_needed(recoil_dir, proj_dir[:])

    # Copy dir_new buffer content into proj["dir"] buffer
    proj["dir"][:] = dir_new

    # energy after scattering
    recoil_e = denfac * proj_e * sin_half_theta**2
    proj["e"] -= recoil_e

    return recoil_dir[:], recoil_e