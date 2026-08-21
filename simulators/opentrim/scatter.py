"""Treat the scattering of a projectile on a target atom.

Currently, the ZBL potential (Ziegler, Biersack, Littmark, The Stopping and 
Range of Ions in Matter, Pergamon Press, 1985) and the NLHlin potential
(to be published) are implemented. The scattering integrals are evaluated
numerically, or Biersack"s "magic formula" is usedfor the scattering angle.

Available functions:
    scatter: treat a scattering event.
"""
import math
import numpy as np
from numba import jit, from_dtype
from numba.core.extending import register_jitable
from . import nlhlin
from . import zbl
from . import config
from .cm_scatter import scatter_integrals


@register_jitable
def screen_fun_wrapper(r, pot_model, pot_coefs):
    """Wrapper to call the appropriate screening function based on the 
    potential model.

    Parameters:
        r (float): Distance between the ion and the target atom (RNORM)
        pot_model (str): The potential model to use ("ZBL" or "NLHlin")
        pot_coefs: Parameters for the potential model, including the model type.
    
    Returns:
        (float): Screening function value
        (float): Derivative of the screening function
    """
    if pot_model == "ZBL":
        return zbl.screen_fun(r, pot_coefs)
    elif pot_model == "NLHlin":
        return nlhlin.screen_fun(r, pot_coefs)
    else:
        raise ValueError(f"Potential model {pot_model} not recognized")


@jit(inline = "always", debug=config.DEBUG)
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


@jit(debug=config.DEBUG)
def scatter(proj, p, dirp, recoil, params):
    """Treat a scattering event.

    The atomic numbers and masses of the ion and the target atom enter the
    calculation via the ENORM, PNORM, DIRFAC, and DENFAC parameters.

    The direction vectors proj["dir"] and dirp[:] are assumed to be normalized 
    to unit length.

    Note: Attempts to assign zbl.screen_fun or nlhlin.screen_fun to a variable 
    and to pass it to scatter_integrals failed.

    Parameters:
        proj (Projectile): state of the projectile (modified in-place)
        p (float): impact parameter (A)
        dirp (ndarray): direction vector of the impact parameter
            (= from the collision point to the recoil position before 
            the collision) (unit vector, size 3)
        recoil (Projectile): the recoil projectile (modified in-place)
        params (PARAMS_DTYPE): Simulation parameters
    """
    p = max(p, 1e-10)  # Avoid 0/0 in p*tan(theta/2)

    # Some abbreviations
    ielem1 = proj["ielem"]
    ielem2 = recoil["ielem"]
    enorm = params.scatter.enorm[ielem1, ielem2]
    rnorm = params.scatter.rnorm[ielem1, ielem2]
    dirfac = params.scatter.dirfac[ielem1, ielem2]
    denfac = params.scatter.denfac[ielem1, ielem2]
    integrate_algorithm = params.scatter.integrate_algorithm
    pot_model = params.scatter.pot_model
    pot_coefs = params.scatter.pot_coefs[ielem1, ielem2]
    
    # Scattering angle in the CM system and time integral
    if integrate_algorithm == "magic":
        if pot_model == "ZBL":
            cos_half_theta = zbl.magic(proj["e"]/enorm, p/rnorm, pot_coefs)
            sin_half_theta = math.sqrt(1 - cos_half_theta**2)
        else:
            raise ValueError(f"Potential model {pot_model} deactivated for now")
    elif integrate_algorithm == "Legendre":
        pi_minus_theta, tau = scatter_integrals(proj["e"]/enorm, p/rnorm, 
                                                pot_model, pot_coefs)
        sin_half_theta = math.cos(0.5 * pi_minus_theta)
        cos_half_theta = math.sin(0.5 * pi_minus_theta)
        tau *= rnorm
    else:
        raise ValueError(f"Unknown scattering integrals algorithm "
                         f"{integrate_algorithm}")

    # Directions of the recoil and the projectile after the collision
    recoil_dir = dirfac * sin_half_theta * (
        sin_half_theta * proj["dir"][:] + cos_half_theta * dirp[:])
    proj_dir = proj["dir"][:] - recoil_dir[:]
    proj_dir = normalize_if_needed(proj_dir[:], proj["dir"][:])
    recoil_dir = normalize_if_needed(recoil_dir[:], proj["dir"][:])

    # Determine turning points ("tp") of trajectries
    x12 = p * sin_half_theta / cos_half_theta
    #x12 = 0.0
    if integrate_algorithm == "magic":  # tau undefined
        x2 = 0.0
    else:
        x2 = (2.0 - dirfac) * (x12 - tau)
    x1 = x2 - x12  # x1 is Eckstein's x1p
    #x1 = 0.0  # for testing
    #x2 = 0.0  # for testing
    proj_tp = proj["pos"][:] + x1 * proj["dir"][:]
    recoil_tp = recoil["pos"][:] + x2 * proj["dir"][:]

    # Energy transfer to the recoil
    recoil_e = denfac * proj["e"] * sin_half_theta**2

    # New projectile properties
    proj["e"] -= recoil_e
    proj["pos"] = proj_tp[:]
    proj["dir"] = proj_dir[:]
    proj["dffp_old"] = x1
    proj["dffp_new"] = np.dot((recoil["pos"][:] - proj_tp[:]), proj_dir[:])
    #proj["dffp_new"] = 0.0  # for testing

    # New recoil properties
    recoil["e"] = recoil_e
    recoil["pos"] = recoil_tp[:]
    recoil["dir"] = recoil_dir[:]
    recoil["dffp_old"] = 0.0
    recoil["dffp_new"] = 0.0
