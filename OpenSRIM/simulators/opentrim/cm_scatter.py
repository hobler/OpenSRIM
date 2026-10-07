"""Characterize a scattering event in the center-of-mass system.

Calculate scattering angle and time integral for given interatomic potential,
energy, and impact parameter, using Gauss-Legendre quadrature.
"""

from numba import jit
import numpy as np
from scipy.special import roots_legendre
from . import config
from . import nlhlin
from . import zbl


# Immutable rules prepared outside JIT. Select by argument, never by mutating
# globals after compilation. Existing callers retain the four-point rule.
LEGENDRE_ORDERS = (4, 16, 32, 64, 128, 256)
LEGENDRE_RULES = tuple(roots_legendre(n) for n in LEGENDRE_ORDERS)
ROOTS_LEGENDRE = LEGENDRE_RULES[0]


@jit(debug=config.DEBUG)
def _integrands(u, e, p, r0, pot_model, pot_coefs):
    if pot_model == "ZBL":
        phi0, dphi0 = zbl.screen_fun(r0, pot_coefs)
        phi, _ = zbl.screen_fun(r0 / (1 - u**2), pot_coefs)
    else:
        phi0, dphi0 = nlhlin.screen_fun(r0, pot_coefs)
        phi, _ = nlhlin.screen_fun(r0 / (1 - u**2), pot_coefs)

        r = r0 / (1 - u**2)
        phi = np.empty_like(r)
        for i in range(r.size):
            if r[i] >= pot_coefs.rmax:
                phi[i] = 0.0
            else:
                phi[i], _ = nlhlin.screen_fun(r[i], pot_coefs)

    chi = np.where(
        u < 3e-4,
        phi0 - r0 * dphi0,
        (phi0 - phi * (1 - u**2)) / u**2
    )

    h = np.sqrt(r0 / e * chi + p**2 * (2 - u**2))

    integrand_theta = p / h

    integrand_tau = (
        (p**2 + r0 / e * phi / (1 - u**2))
        / (h * (1 + u / r0 * h))
    )

    # Rationalized difference between the free and interacting integrands.
    base = np.sqrt(2-u**2)
    correction = (r0/e*chi) / (base*h*(h+p*base))
    return integrand_theta, integrand_tau, correction


@jit(debug=config.DEBUG)
def _integrand_arccos_half(u):
    return 1 / np.sqrt(2 - u**2)


@jit(debug=config.DEBUG)
def scatter_integrals(e, p, pot_model, pot_coefs, n_absc=4):
    """Calculate scattering angle and time integral.

    The calculation uses Gauss-Legendre quadrature with an explicit order.
    Omitting n_absc preserves the four-point configuration.

    Parameters:
        e (float): Reduced energy.
        p (float): Reduced impact parameter.
        pot_model (str): The potential model to use ("ZBL" or "NLHlin").
        pot_coefs: Parameters for the potential model.
        n_absc (int): Quadrature order: 4, 16, 32, 64, 128 or 256 (default 4).

    Returns:
        (float): Pi minus scattering angle (rad)
        (float): Time integral (RNORM).
    """

    rule_index = -1
    for i in range(len(LEGENDRE_ORDERS)):
        if n_absc == LEGENDRE_ORDERS[i]:
            rule_index = i
            break
    if rule_index < 0:
        raise ValueError("Unsupported Legendre order; use 4, 16, 32, 64, 128 or 256")

    rmax = pot_coefs.rmax

    if p >= rmax:
        return np.pi, 0.0

    r0, _ = get_apsis(e, p, pot_model, pot_coefs)

    if rmax == np.inf:
        umax = 1.0
    else:
        umax = np.sqrt((rmax-r0)/rmax)
        sp = np.sqrt((rmax-p)*(rmax+p))
        s0 = np.sqrt((rmax-r0)*(rmax+r0))
        if umax == 0.0:
            return np.pi-2*np.arctan2(sp, p), r0-p**2/(rmax+sp)

    u_vals, weights = LEGENDRE_RULES[rule_index]
    u_vals = 0.5 * umax * (u_vals + 1)
    weights = 0.5 * umax * weights

    integrand_theta_vals, integrand_tau_vals, correction = _integrands(
        u_vals,
        e,
        p,
        r0,
        pot_model,
        pot_coefs,
    )

    arccos_num_half = np.sum(
        weights * _integrand_arccos_half(u_vals)
    )

    if rmax == np.inf:
        pi_minus_theta = (
            np.pi
            / arccos_num_half
            * np.sum(weights * integrand_theta_vals)
        )

        tau = (
            r0
            - 2 * np.sum(weights * integrand_tau_vals)
        )

    else:
        geometric = 2*np.arcsin((r0-p)*(r0+p)/(r0*sp+p*s0))
        theta = (geometric + 2*np.arctan2(s0, r0) / arccos_num_half
                 * np.sum(weights * correction))
        # Preserve the public API's complementary-angle convention.
        pi_minus_theta = np.pi - theta

        tau = (
            r0
            - (rmax - np.sqrt(rmax**2 - p**2))
            - 2 * np.sum(weights * integrand_tau_vals)
        )

    return pi_minus_theta, tau


@jit(debug=config.DEBUG)
def get_apsis(e, p, pot_model, pot_coefs):
    """Calculate the distance of closest approach (apsis).

    The initial condition for Newton's method is calculated from a
    piecewise approximation to the screening function.

    Parameters:
        e (float): Reduced energy.
        p (float): Reduced impact parameter.
        pot_model (str): Repulsive potential model ("ZBL" or "NLHlin").
        pot_coefs: Parameters returned by nlhlin.get_coefs().

    Returns:
        r0 (float): Estimated apsis of the collision (RNORM)
        count (int): Number of iterations used to converge the apsis
    """

    psq = p**2

    # NLHlin potential coefficients:
    k = pot_coefs.k
    r34 = pot_coefs.r34
    rmax = pot_coefs.rmax

    if p >= rmax:
        return p, 0

    k1, k2, k3 = k[:3]

    # Initial condition for infinite-range potentials.
    #
    # This is the TRIM85-style initial estimate.
    if rmax == np.inf:
        r0 = max(1e-10, p)

        r0_try = -2.7 * np.log(e * r0)

        if r0_try > p:
            r0_try = -2.7 * np.log(e * r0_try)

            if r0_try > p:
                r0 = r0_try

        done = r0 > r34

    else:
        # Finite-range potential.
        if psq > r34**2 - k3 / (e * r34**2):
            k4 = k[3]
            a_quad = e + k4
            b_quad = -k4 * rmax
            c_quad = -e * psq

            r0 = (
                -b_quad
                + np.sqrt(b_quad**2 - 4 * a_quad * c_quad)
            ) / (2 * a_quad)

            done = True
        else:
            done = False

    # Initial condition from the piecewise approximation.
    if not done:
        r0sq = psq + k2 / e

        if r0sq > k3 / k2:
            r0 = np.sqrt(
                psq / 2
                + np.sqrt(psq**2 / 4 + k3 / e)
            )

        elif r0sq >= k2 / k1:
            r0 = np.sqrt(r0sq)

        else:
            r0 = (
                1
                + np.sqrt(1 + 4 * e * (e + k1) * psq)
            ) / (2 * (e + k1))

    # The repulsive-potential root cannot lie below the impact parameter.
    # Safeguard iterates, not the returned scattering angle.
    r0 = max(p, r0)
    if np.isfinite(rmax):
        r0 = min(r0, rmax)

    def fun(r):
        if pot_model == "ZBL":
            screen, dscreen = zbl.screen_fun(r, pot_coefs)
        else:
            screen, dscreen = nlhlin.screen_fun(r, pot_coefs)
        return (r - p) * (1 + p/r) - screen/e, 1 - dscreen/e + p**2/r**2

    for count in range(1, 101):
        f, df = fun(r0)
        if not np.isfinite(f) or not np.isfinite(df) or df <= 0:
            raise ValueError("Invalid Newton step in get_apsis")
        delta_r0 = -f / df
        candidate = max(p, r0 + delta_r0)
        if np.isfinite(rmax):
            candidate = min(candidate, rmax)
        step = candidate - r0
        r0 = candidate
        if abs(step) <= 1e-3 * r0:
            return r0, count

    raise RuntimeError("get_apsis did not converge within 100 iterations")
