"""Characterize a scattering event in the center-of-mass system.

Calculate scattering angle and time integral for given interatomic potential,
energy, and impact parameter, using Gauss-Legendre quadrature.
"""
# import os
# os.environ["NUMBA_DISABLE_JIT"] = "1"
from numba import jit
import numpy as np
from scipy.special import roots_legendre
from . import config
from . import nlhlin
from . import zbl

# TODO: Put ROOTS_LEGENDRE into params, since global variables might not
# work with cached jit functions.
ROOTS_LEGENDRE = roots_legendre(4)


@jit(debug=config.DEBUG)
def scatter_integrals(e, p, pot_model, pot_coefs):
    """Calculate scattering angle and time integral.

    The calculation uses Gauss-Legendre quadrature with a fixed number of
    abscissae given by ROOTS_LEGENDRE. The abscrissae and weights must have 
    been set up before calling this function.

    Parameters:
        e (float): Reduced energy.
        p (float): Reduced impact parameter.
        pot_model (str): The potential model to use ("ZBL" or "NLHlin").
        pot_coefs: Parameters for the potential model.
    
    Returns:
        (float): Pi minus scattering angle (rad)
        (float): Time integral (RNORM).
    """
    if p >= pot_coefs.rmax:
        return np.pi, 0.0
    
    r0, _ = get_apsis(e, p, pot_model, pot_coefs)

    def calc_phi_chi(u):
        if pot_model == "ZBL":
            phi0, dphi0 = zbl.screen_fun(r0, pot_coefs)
            phi, _ = zbl.screen_fun(r0/(1-u**2), pot_coefs)
        else:
            phi0, dphi0 = nlhlin.screen_fun(r0, pot_coefs)
            phi, _ = nlhlin.screen_fun(r0/(1-u**2), pot_coefs)
        chi = np.where(u < 3e-4, phi0 - r0*dphi0, (phi0 - phi*(1-u**2)) / u**2)
        return phi, chi

    def integrands(u):
        phi, chi = calc_phi_chi(u)
        h = np.sqrt(r0/e*chi + p**2*(2-u**2))
        integrand_theta = p / h
        integrand_tau = (p**2 + r0/e * phi/(1-u**2)) / (h * (1 + u/r0*h))
        return integrand_theta, integrand_tau
    
    def integrand_arccos_half(u):
        return 1 / np.sqrt(2 - u**2)

    rmax = pot_coefs.rmax
    umax = np.sqrt(1 - r0 / rmax)
    u_vals, weights = ROOTS_LEGENDRE
    u_vals = 0.5 * umax * (u_vals + 1)
    weights = 0.5 * umax * weights
    integrand_theta_vals, integrand_tau_vals = integrands(u_vals)
    arccos_num_half = np.sum(weights * integrand_arccos_half(u_vals))
    if rmax == np.inf:
        pi_minus_theta = (np.pi / arccos_num_half 
                          * np.sum(weights * integrand_theta_vals))
        tau = (r0 - 2 * np.sum(weights * integrand_tau_vals))
    else:    
        pi_minus_theta = (2 * np.arccos(r0/rmax) / arccos_num_half 
                          * np.sum(weights * integrand_theta_vals)
                          + 2 * np.arcsin(p/rmax))
        tau = (r0 - (rmax - np.sqrt(rmax**2 - p**2)) 
            - 2 * np.sum(weights * integrand_tau_vals))

    return pi_minus_theta, tau


@jit(debug=config.DEBUG)
def get_apsis(e, p, pot_model, pot_coefs):
    """Calculate the distance of closest approach (apsis) in a colllision.

    As initial condition for Newton's method is calculated from a 
    piecewise approximation to the screening function.

    Parameters:
        e (float): energy of projectile before the collision (ENORM)
        p (float): impact parameter (RNORM)
        pot_model (str): The potential model to use ("ZBL" or "NLHlin").
        pot_coefs: Parameters for the potential model.

    Returns:
        r0 (float): Estimated apsis of the collision (RNORM)
        count(int): Number of iterations used to converge the apsis
    """
    psq = p**2
    k1, k2, k3, k4 = pot_coefs.k[:4]
    r34 = pot_coefs.r34
    rmax = pot_coefs.rmax
    
    # Initial condition: Assume r0 > r34
    if rmax is np.inf:  # Use TRIM85 algorithm
        r0 = max(1e-10, p)
        r0_try = -2.7 * np.log(e*r0)
        if r0_try > p:
            r0_try = -2.7 * np.log(e*r0_try)
            if r0_try > p:
                r0 = r0_try
        done = r0 > r34
    else:
        if psq > r34**2 - k3/(e*r34**2):
            a = e + k4
            b = - k4 * rmax
            c = - e * psq
            r0 = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)
            done = True
        else:
            done = False

    # Initial condition: Use piecewise approximation if r0 <= r34
    if not done:
        r0sq = psq + k2/e
        if r0sq > k3 / k2:
            r0 = np.sqrt(psq/2 + np.sqrt(psq**2/4 + k3/e))
        elif r0sq >= k2 / k1:
            r0 = np.sqrt(r0sq)
        else:
            r0 = (1 + np.sqrt(1 + 4*e*(e+k1)*psq)) / (2*(e+k1))

    # Newton iteration
    delta_r0 = np.inf

    def fun(r):
        if pot_model == "ZBL":
            screen, dscreen = zbl.screen_fun(r, pot_coefs)
        else:
            screen, dscreen = nlhlin.screen_fun(r, pot_coefs)
        return r - screen[0]/e - p**2/r, 1 - dscreen[0]/e + p**2/r**2
    
    count = 0
    while abs(delta_r0) > 1e-3 * r0:
        f, df = fun(r0)
        delta_r0 = - f / df
        r0 += delta_r0
        count += 1

    return r0, count
