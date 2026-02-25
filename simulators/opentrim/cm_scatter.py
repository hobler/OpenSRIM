"""Characterize a scattering event in the center-of-mass system.

Calculate scattering angle and time integral for given interatomic potential,
energy, and impact parameter, using Gauss-Legendre quadrature.
"""
# import os
# os.environ["NUMBA_DISABLE_JIT"] = "1"
from numba import jit
import numpy as np
from scipy.special import roots_legendre
from .zbl import zbl_screen, estimate_apsis

# TODO: Put ROOTS_LEGENDRE into params, since global variables might not
# work with cached jit functions.
ROOTS_LEGENDRE = roots_legendre(4)


@jit
def calc_phi_chi(u, r0, screen_fun):
    """Calculate the screening function and chi at the given u values.

    phi is evaluated at r0/(1-u**2).
    chi is the difference of phi and phi0 (=phi(r0)) divided by u^2, with a 
    Taylor expansion used for small u to avoid numerical issues.

    Parameters:
        u (array-like): Integration variable.
        r0 (float): Distance of closest approach (RNORM).
        screen_fun (callable object): Function to calculate the screening 
            function for given distance r (RNORM).
    Returns:
        (float): Value of chi function.
    """
    # u = np.asarray(u)
    phi0, dphi0 = screen_fun(r0)
    phi, _ = screen_fun(r0/(1-u**2))
    chi = np.where(u < 3e-4, phi0 - r0*dphi0, (phi0 - phi*(1-u**2)) / u**2)
    return phi, chi


@jit
def scatter_integrals(e, p, pot_model):
    """Calculate scattering angle and time integral.

    The calculation uses Gauss-Legendre quadrature with a fixed number of
    abscissae given by ROOTS_LEGENDRE. The abscrissae and weights must have 
    been set up before calling this function.

    Parameters:
        e (float): Reduced energy.
        p (float): Reduced impact parameter.
        pot_model (str): Name of the potential model to be used.
    
    Returns:
        (float): Scattering angle (rad)
        (float): Time integral (RNORM).
    """
    if pot_model.startswith("ZBL"):
        screen_fun = zbl_screen
        rmax = np.inf
    elif pot_model.startswith("NLHlin"):
        raise NotImplementedError("NLHlin potential not implemented yet")
    else:
        raise ValueError(f"Unknown potential model {pot_model}")
    
    if p >= rmax:
        return 0.0, 0.0
    elif p == 0.0:
        return np.pi, 0.0
    
    # TODO: Use more general apsis calculation ffrom the apsis module
    if pot_model.startswith("ZBL"):
        r0 = estimate_apsis(e, p)
    else:
        raise NotImplementedError("Apsis estimation not implemented for this "
                                  "potential model")

    def integrands(u):
        phi, chi = calc_phi_chi(u, r0, screen_fun)
        rho = r0 / (e*p**2)
        g = np.sqrt(rho*chi + (2-u**2))
        integrand_theta = 1 / g
        integrand_tau = (1 + rho * phi/(1-u**2)) / (g * (1 + u*p/r0*g))
        return integrand_theta, integrand_tau
    
    def integrand_two_arccos(u):
        return 4 / np.sqrt(2 - u**2)

    umax = np.sqrt(1 - r0 / rmax)
    u_vals, weights = ROOTS_LEGENDRE
    u_vals = 0.5 * umax * (u_vals + 1)
    weights = 0.5 * umax * weights
    integrand_theta_vals, integrand_tau_vals = integrands(u_vals)
    two_arccos_num = np.sum(weights * integrand_two_arccos(u_vals))
    if rmax == np.inf:
        theta = np.pi * (1 
                - 4/two_arccos_num * np.sum(weights * integrand_theta_vals))
        tau = (r0 - 2 * p * np.sum(weights * integrand_tau_vals))
    else:    
        theta = 2*np.arccos(p/rmax) * (1 
                - 4/two_arccos_num * np.sum(weights * integrand_theta_vals))
        tau = (r0 - (rmax - np.sqrt(rmax**2 - p**2)) 
            - 2 * p * np.sum(weights * integrand_tau_vals))

    return theta, tau
