"""Models specific to the universal ZBL potential.

The ZBL potential is described in Ziegler, Biersack, Littmark,
The Stopping and Range of Ions in Matter, Pergamon Press, 1985).

The functions in this module use normalized units:
- Distances are given in units of RNORM = 0.4685 / (Z1**0.23 + Z2**0.23) A
- Energies are given in units of 
      ENORM = 14.39979 * Z1 * Z2 / RNORM * (1 + m1/m2) eV

Available functions:
    zbl_screen: Callable object for the ZBL screening function.
    estimate_apsis: estimate distance of closest approach in a collision.
    magic: calculate scattering angle using Biersack's magic formula.
"""
from math import sqrt
import numpy as np
from scipy import special
from scipy.optimize import brentq
from numba import jit
from numba.core.extending import register_jitable
from . import config


# Constants for ZBL screening function
A0 = 0.18175
A1 = 0.50986
A2 = 0.28022
A3 = 0.02817

B0 = 3.1998
B1 = 0.94229
B2 = 0.4029
B3 = 0.20162

def get_coefs(Z1, Z2):
    """Read ZBL coefficients and calculate derived parameters.

    Parameters:
        Z1 (int): Atomic number of atom 1
        Z2 (int): Atomic number of atom 2

    Returns:
        rnorm (float): Normalization length (RNORM)
        a (ndarray): Coefficients for the screening function
        b (ndarray): Screening lengths (1/RNORM)
        r34 (float): Maximum range for piecewise approximation of the 
            potential (RNORM)
        k (ndarray): Array of coefficients for the piecewise approximation
    """
    rnorm = 0.4685 / (Z1**0.23 + Z2**0.23)
    a = np.array([A0, A1, A2, A3])
    b = np.array([B0, B1, B2, B3])

   # Calculate parameters for the piecewise approximation
    def screen_fun(r):
        exp = np.exp(-b*r)
        screen = np.sum(a*exp)
        dscreen = - np.sum(a*b*exp)
        return screen, dscreen
    
    # k2/R part
    def fun2(r):
        screen, dscreen = screen_fun(r)
        return screen + r * dscreen

    r_touch2 = brentq(fun2, 0.0, 40.0)
    k2 = r_touch2 * screen_fun(r_touch2)[0]

    # k3/R^3 part
    def fun3(r):
        screen, dscreen = screen_fun(r)
        return screen + 1/3 * r * dscreen

    r_touch3 = brentq(fun3, 0.0, 40.0)
    k3 = r_touch3**3 * screen_fun(r_touch3)[0]

    # 1 - k1*R part
    k1 = 1 / (4*k2)

    k = np.array([k1, k2, k3])
    
    # Radius beyond which the piecewise approximation is not used anymore
    r34 = r_touch3

    return rnorm, a, b, r34, k


def impulse_integral(p, pot_coefs):
    """Evaluate the integral that appears in the impulse approximation.
    
    This integral is defined as one half of the integral of 

         d  Phi(r)
        --  ------
        dp    r
    
    along a straight line which passes by the center of the potential at 
    a distance p.

    The calculation uses the modified Bessel function of the second kind 
    and order 1 (scipy.special.kn).
    
    Parameters:
        p (float): impact parameter (RNORM)

    Returns:
        (float): value of the integral
    """
    k0 = special.kn(1, pot_coefs.b[0] * p)
    k1 = special.kn(1, pot_coefs.b[1] * p)
    k2 = special.kn(1, pot_coefs.b[2] * p)
    k3 = special.kn(1, pot_coefs.b[3] * p)

    integral = (pot_coefs.a[0]*pot_coefs.b[0]*k0 
                + pot_coefs.a[1]*pot_coefs.b[1]*k1 
                + pot_coefs.a[2]*pot_coefs.b[2]*k2 
                + pot_coefs.a[3]*pot_coefs.b[3]*k3)

    return integral


@register_jitable(debug=config.DEBUG)
def screen_fun(r, pot_coefs):
    """Calculate the ZBL screening function and its derivative.

    Parameters:
        r (float): Distance (RNORM)
        pot_coefs: Parameters needed for the evaluation of the
            screening function
            
    Returns:
        (float): ZBL screening function at distance r
        (float): derivative of ZBL screening function at distance r
            (1/RNORM)
    """
    r = np.asarray(r)
    a = pot_coefs.a[:4]
    b = pot_coefs.b[:4]

    exp = np.exp(-np.outer(r, b))
    screen = np.sum(a*exp, axis=1)
    dscreen = - np.sum(a*b*exp, axis=1)

    return screen, dscreen




# Constants for apsis estimation for the ZBL potential
K2 = 0.38           # factor of the 1/R part
K3 = 7.2            # factor of the 1/R^3 part
K1 = 1/(4*K2)
R12sq = (2*K2)**2
R23sq = K3 / K2
NITER = 1           # number of Newton-Raphson iterations

@jit(debug=config.DEBUG)
def estimate_apsis(e, p, pot_coefs):
    """Estimate the distance of closest approach (apsis) in a colllision.

    Parameters:
        e (float): energy of projectile before the collision (ENORM)
        p (float): impact parameter (RNORM)
        pot_coefs: Parameters needed for the evaluation of the
            screening function

    Returns:
        (float): Estimated apsis of the collision (RNORM)
    """
    psq = p**2
    r0sq = 0.5 * (psq + sqrt(psq**2 + 4*K3/e))

    if r0sq < R23sq:
        r0sq = psq + K2/e
        if r0sq < R12sq:
            r0 = (1 + sqrt(1 + 4*e*(e+K1)*psq)) / (2*(e+K1))
        else:
            r0 = sqrt(r0sq)
    else:
        r0 = sqrt(r0sq)
    
    # Do Newton-Raphson iterations to improve the estimate
    r0 = np.asarray(r0)
    for _ in range(NITER):
        screen, dscreen = screen_fun(r0, pot_coefs)
        numerator = r0*(r0-screen[0]/e) - p**2
        denominator = 2*r0 - (screen[0]+r0*dscreen[0])/e
        r0 -= numerator/denominator

        residuum = 1 - screen[0]/(e*r0) - p**2/r0**2
        if abs(residuum) < 1e-4:
            break

    return r0

# Constants for Biersack's magic formula
C1 = 0.99229
C2 = 0.011615
C3 = 0.007122
C4 = 14.813
C5 = 9.3066

@jit(debug=config.DEBUG)
def magic(e, p, pot_coefs):
    """Calculate CM scattering angle using Biersack's magic formula.

    Parameters:
        e (float): energy of projectile before the collision (ENORM)
        p (float): impact parameter (RNORM)
        pot_coefs: Parameters needed for the evaluation of the
            screening function

    Returns:
        (float): cosine of half the scattering angle in the center-of-mass 
            system
    """
    r0 = estimate_apsis(e, p, pot_coefs)
    screen, dscreen = screen_fun(r0, pot_coefs)

    rho = 2*(e*r0-screen[0]) / (screen[0]/r0-dscreen[0])
    sqrte = sqrt(e)
    alpha = 1 + C1/sqrte
    beta = (C2+sqrte) / (C3+sqrte)
    gamma = (C4+e) / (C5+e)
    a = 2 * alpha * e * p**beta
    g = gamma / (sqrt(1+a**2)-a)
    delta = a * (r0-p) / (1+g)

    cos_half_theta = (p + rho + delta) / (r0 + rho)
    if cos_half_theta > 1:
        print("Warning: cos_half_theta > 1:", cos_half_theta)
        print("  e =", e, "p =", p, "r0 =", r0, "rho =", rho, "delta =", delta)

    return cos_half_theta


