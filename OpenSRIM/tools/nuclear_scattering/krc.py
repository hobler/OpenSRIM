"""Defines the KrC screening function.
"""
import numpy as np

from apsis import Apsis


# Constants for KrC screening function
A0 = 0.335381
A1 = 0.473674
A2 = 0.190945

B0 = 1.919249
B1 = 0.637174
B2 = 0.278544

A0B0 = A0 * B0
A1B1 = A1 * B1
A2B2 = A2 * B2

# Constants for Biersack's magic formula (after Eckstein, 1991)
C1 = 1.0144
C2 = 0.235809
C3 = 0.126
C4 = 69350
C5 = 83550


class KrC_screen:
    """Defines the KrC screening function.
    """
    def __init__(self, Z1=None, Z2=None, rnorm=None, magic=False):
        """Setup KrC screening function for given screening length.

        Define global variables for the KrC screening function.
        
        Parameters:
            Z1: atomic number of atom 1
            Z2: atomic number of atom 2
            rnorm: screening length (A)
        """
        self.Z1 = Z1
        self.Z2 = Z2

        if rnorm is None:
            self.a = (A0, A1, A2)
            self.b = (B0, B1, B2)
            self.ab = (A0B0, A1B1, A2B2)
        else:
            if Z1 is None or Z2 is None:
                print("KrCsetup: Z1 and Z2 must be given if rnorm is given")
                return
            rnormKrC = 0.4685 / (np.sqrt(Z1) + np.sqrt(Z2))**(2/3)
            factor = rnorm / rnormKrC
            self.a = (A0, A1, A2)
            self.b = (B0*factor, B1*factor, B2*factor)   
            self.ab = (A0B0*factor, A1B1*factor, A2B2*factor)
        
        self.rmax = np.inf  # KrC screening function is defined for all r
        self.apsis = Apsis(self)

    def __call__(self, r):
        """Calculate the KrC screening function and its derivative.

        Parameters:
            r (float): Distance (RNORM)

        Returns:
            (float): KrC screening function at distance r
            (float): derivative of KrC screening function at distance r 
                (1/RNORM)
        """
        r = np.asarray(r, dtype=float)

        exp0 = np.exp(-self.b[0] * r)
        exp1 = np.exp(-self.b[1] * r)
        exp2 = np.exp(-self.b[2] * r)

        screen = self.a[0]*exp0 + self.a[1]*exp1 + self.a[2]*exp2
        dscreen = - self.ab[0]*exp0 - self.ab[1]*exp1 - self.ab[2]*exp2
        
        return screen, dscreen

    def magic(self, e, p):
        """Calculate CM scattering angle using Biersack's magic formula.

        Note: This function has not been tested.

        Parameters:
            e (float): energy of projectile before the collision (ENORM)
            p (float): impact parameter (RNORM)
        
        Returns:
            (float): cosine of half the scattering angle in the center-of-mass 
                system
        """
        r0 = self.apsis(e, p)
        screen, dscreen = self(r0)

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

        return cos_half_theta


