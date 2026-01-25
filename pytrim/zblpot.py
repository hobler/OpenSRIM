"""Functions specific to the ZBL interatomic potential.

The ZBL potential (Ziegler, Biersack, Littmark is decribed in "The Stopping and 
Range of Ions in Matter, Pergamon Press, 1985". The functions in this module 
use normalized units:
- Distances are given in units of 
      RNORM = 0.4685 / (z1**0.23 + z2**0.23) A
- Energies are given in units of 
      ENORM = 14.39979 * z1 * z2 / RNORM * (1 + m1/m2) eV

Available functions:
    ZBLscreen: calculate ZBL screening function.
    estimate_apsis: estimate distance of closest approach in a collision.
    magic: calculate scattering angle using Biersack's magic formula.
"""
import numpy as np


# Constants for ZBL screening function
A1 = 0.18175
A2 = 0.50986
A3 = 0.28022
A4 = 0.02817

B1 = 3.1998
B2 = 0.94229
B3 = 0.4029
B4 = 0.20162

A1B1 = A1 * B1
A2B2 = A2 * B2
A3B3 = A3 * B3
A4B4 = A4 * B4

def ZBLscreen(r):
    """Calculate the ZBL screening function and its derivative.

    Parameters:
        r (float): Distance (RNORM)

    Returns:
        (float): ZBL potential at distance r (ENORM)
        (float): derivative of ZBL potential at distance r (ENORM/RNORM)
    """
    exp1 = np.exp(-B1 * r)
    exp2 = np.exp(-B2 * r)
    exp3 = np.exp(-B3 * r)
    exp4 = np.exp(-B4 * r)
    screen = A1*exp1 + A2*exp2 + A3*exp3 + A4*exp4
    dscreen = - (A1B1*exp1 + A2B2*exp2 + A3B3*exp3 + A4B4*exp4)
    
    return screen, dscreen


# Constants for apsis estimation for the ZBL potential
K2 = 0.38           # factor of the 1/R part
K3 = 7.2            # factor of the 1/R^3 part
K1 = 1/(4*K2)
R12sq = (2*K2)**2
R23sq = K3 / K2
NITER = 1           # number of Newton-Raphson iterations

def estimate_apsis(e, p):
    """Estimate the distance of closest approach (apsis) in a colllision.

    Parameters:
        e (float): energy of projectile before the collision (ENORM)
        p (float): impact parameter (RNORM)

    Returns:
        (float): Estimated apsis of the collision (RNORM)
    """
    psq = p**2
    r0sq = 0.5 * (psq + np.sqrt(psq**2 + 4*K3/e))

    if r0sq < R23sq:
        r0sq = psq + K2/e
        if r0sq < R12sq:
            r0 = (1 + np.sqrt(1 + 4*e*(e+K1)*psq)) / (2*(e+K1))
        else:
            r0 = np.sqrt(r0sq)
    else:
        r0 = np.sqrt(r0sq)
    
    # Do Newton-Raphson iterations to improve the estimate
    for _ in range(NITER):
        screen, dscreen = ZBLscreen(r0)
        numerator = r0*(r0-screen/e) - p**2
        denominator = 2*r0 - (screen+r0*dscreen)/e
        r0 -= numerator/denominator

        residuum = 1 - screen/(e*r0) - p**2/r0**2
        if abs(residuum) < 1e-4:
            break

    return r0


C1 = 0.99229
C2 = 0.011615
C3 = 0.007122
C4 = 14.813
C5 = 9.3066

def magic(e, p):
    """Calculate CM scattering angle using Biersack's magic formula.

    Parameters:
        e (float): energy of projectile before the collision (ENORM)
        p (float): impact parameter (RNORM)
    
    Returns:
        (float): sine of half the scattering angle in the center-of-mass 
            system
        (float): cosine of half the scattering angle in the center-of-mass 
            system
    """
    r0 = estimate_apsis(e, p)
    screen, dscreen = ZBLscreen(r0)

    rho = 2*(e*r0-screen) / (screen/r0-dscreen)
    sqrte = np.sqrt(e)
    alpha = 1 + C1/sqrte
    beta = (C2+sqrte) / (C3+sqrte)
    gamma = (C4+e) / (C5+e)
    a = 2 * alpha * e * p**beta
    g = gamma / (np.sqrt(1+a**2)-a)
    delta = a * (r0-p) / (1+g)

    cos_half_theta = (p + rho + delta) / (r0 + rho)
    if cos_half_theta > 1:
        print("Warning: cos_half_theta > 1:", cos_half_theta)
        print("  e =", e, "p =", p, "r0 =", r0, "rho =", rho, "delta =", delta)

    sin_half_theta = np.sqrt(1 - cos_half_theta**2)

    return sin_half_theta, cos_half_theta


