"""Define the NLHlin screening function.

Available functions:
    NLHlin_screen: Callable object for the NLHlin screening function.
    
Moreover, there are plotting functions for testing and visualization:
    plot_screen: plot NLHlin screening function for various atomic numbers.
    plot_ZBLscreen: plot ZBL screening function for various atomic numbers.
"""
import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from numba import jit
from numba.core.extending import register_jitable
from . import config


def get_coefs(Z1, Z2):
    """Read NLHlin coefficients and calculate derived parameters.

    Parameters:
        Z1 (int): Atomic number of atom 1
        Z2 (int): Atomic number of atom 2

    Returns:
        rnorm (float): Normalization length (RNORM)
        a (ndarray): Coefficients for the screening function
        b (ndarray): Screening lengths (1/RNORM)
        c (float): Constant term
        d (float): Coefficient for the linear term
        rmax (float): Maximum distance for the screening function
        k (ndarray): Array of coefficients for the piecewise approximation
    """
    if Z1 > Z2:
        Z1, Z2 = Z2, Z1

    fname = os.path.join(os.path.dirname(__file__), "dmol_coeffs_rmax.dat")
    if not os.path.exists(fname):
        print(f"get_nlhlin_coefs: Coefficients file {fname} not found")
        sys.exit()
    
    # Search for line with the correct Z1 and Z2, and read the coefficients
    found = False
    with open(fname, "r") as f:
        for line in f:
            if line[0] == "#":
                continue
            coefs = line.split()[:-1]   # exclude "error" column
            Z1_ = int(coefs[0])
            Z2_ = int(coefs[1])
            if (Z1_, Z2_) == (Z1, Z2):
                found = True
                break
    
    if not found:
        print(f"get_nlhlin_coefs: Coefficients for Z1={Z1} and Z2={Z2} "
              "not found \n" f"in {fname}")
        sys.exit()

    # Extract a, b, and rmax
    a = np.array([float(coefs[i]) for i in range(2, 7, 2)])
    b = np.array([float(coefs[i]) for i in range(3, 8, 2)])
    rmax = float(coefs[8])

    # Calculate derived parameters for the NLHlin function
    rnorm = 0.4685 / np.sqrt(np.sqrt(Z1) + np.sqrt(Z2))
    b *= rnorm
    rmax /= rnorm
    c = 1 - np.sum(a[:])
    d = np.sum(a[:] * np.exp(-b[:] * rmax)) + c

    # Calculate parameters for the piecewise approximation
    def screen_fun(r):
        exp = np.exp(-b*r)
        screen = np.sum(a*exp) + c - d*r/rmax
        dscreen = - np.sum(a*b*exp) - d/rmax
        return screen, dscreen
    
    # k2/R part
    def fun2(r):
        screen, dscreen = screen_fun(r)
        return screen + r * dscreen

    r_touch2 = brentq(fun2, 0.0, rmax)
    k2 = r_touch2 * screen_fun(r_touch2)[0]

    # k3/R^3 part
    def fun3(r):
        screen, dscreen = screen_fun(r)
        return screen + 1/3 * r * dscreen

    r_touch3 = brentq(fun3, 0.0, rmax)
    k3 = r_touch3**3 * screen_fun(r_touch3)[0]

    # 1 - k1*R part
    k1 = 1 / (4*k2)

    # k4*(Rmax-R) part
    r34 = 0.75 * rmax
    k4 = 3*k3 / r34**4

    k = np.array([k1, k2, k3, k4])

    return rnorm, a, b, c, d, r34, rmax, k


def impulse_integral(p, pot_coefs):
    """Evaluate the integral that appears in the impulse approximation.
    
    This integral is defined as one half of the integral of 
    
        Phi(r) - r*Phi'(r)
        ------------------ * p
                r^3

    along a straight line which passes by the center of the potential at 
    a distance p.

    The calculation uses quad from SciPy for numerical integration.
    
    Parameters:
        p (float): impact parameter (RNORM)

    Returns:
        (float): value of the integral
    """
    if p >= pot_coefs.rmax:
        return 0.0
    
    def integrand(x, p):
        r = np.sqrt(x**2 + p**2)
        screen, dscreen = screen_fun(r, pot_coefs)
        return (screen[0] - r*dscreen[0]) * p / r**3
    
    xmax = np.sqrt(pot_coefs.rmax**2 - p**2) if p < pot_coefs.rmax else 0
    integral, abserr = quad(integrand, 0, xmax, args=(p,))
    #print(p, xmax, pot_coefs.rmax, integral, abserr)
    
    return integral


@register_jitable(debug=config.DEBUG)
def screen_fun(r, pot_coefs):
    """Calculate the NLHlin screening function and its derivative.
    
    Parameters:
        r (float): Distance (RNORM)
        pot_coefs: Parameters needed for the evaluation of the
            screening function
            
    Returns:
        (float): NLHlin screening function at distance r
        (float): derivative of NLHlin screening function at distance r
            (1/RNORM)
    """
    r = np.asarray(r)
    a = pot_coefs.a[:3]
    b = pot_coefs.b[:3]
    c = pot_coefs.c
    d = pot_coefs.d
    rmax = pot_coefs.rmax

    exp = np.exp(-np.outer(r, b))
    screen = np.sum(a*exp, axis=1) + c - d*r/rmax
    dscreen = - np.sum(a*b*exp, axis=1) - d/rmax

    return screen, dscreen

