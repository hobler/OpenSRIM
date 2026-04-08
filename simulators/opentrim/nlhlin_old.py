"""Define the NLHlin screening function.

Available functions:
    NLHlin_screen: Callable object for the NLHlin screening function.
    
Moreover, there are plotting functions for testing and visualization:
    plot_screen: plot NLHlin screening function for various atomic numbers.
    plot_ZBLscreen: plot ZBL screening function for various atomic numbers.
"""
import os
import sys
from numba.core.types import UniTuple, float64
from numba.experimental import jitclass
import numpy as np
from .apsis import Apsis


@jitclass
class NLHlin_screen:
    """Defines the NLHlin screening function.
    """
    Z1: int
    Z2: int
    a: UniTuple(float64, 3)  # ty:ignore[invalid-type-form]
    b: UniTuple(float64, 3) # ty:ignore[invalid-type-form]
    ab: UniTuple(float64, 3)    # ty:ignore[invalid-type-form]
    c: float
    d: float
    rmax: float
    _apsis: Apsis
    def __init__(self, Z1, Z2, rnorm=None, coefs=None):
        """Setup NLHlin screening function for given atomic numbers.

        Parameters:
            Z1: atomic number of atom 1
            Z2: atomic number of atom 2
            rnorm: screening length (A), None for the default value of
                0.4685 / sqrt(sqrt(Z1) + sqrt(Z2))
            coefs: Table of coefficients (for Z1 and Z2)
        """
        if coefs is None:
            print("Required coefficients table missing for NLHlin_screen")
            return
        self.Z1 = Z1
        self.Z2 = Z2
        mask = ((coefs.z1 == Z1) & (coefs.z2 == Z2)) | \
               ((coefs.z1 == Z2) & (coefs.z2 == Z1))
        rec = coefs[mask][0]

        self.a = (rec.a1, rec.a2, rec.a3)
        if rnorm is None:
            rnorm = 0.4685 / np.sqrt(np.sqrt(Z1) + np.sqrt(Z2))
        self.b = (rec.b1*rnorm, rec.b2*rnorm, rec.b3*rnorm)
        self.ab = (self.a[0]*self.b[0], self.a[1]*self.b[1], 
                   self.a[2]*self.b[2])

        self.rmax = rec.rmax / rnorm
        self.c = 1 - self.a[0] - self.a[1] - self.a[2]
        self.d = (self.a[0]*np.exp(-self.b[0]*self.rmax) 
                  + self.a[1]*np.exp(-self.b[1]*self.rmax) 
                  + self.a[2]*np.exp(-self.b[2]*self.rmax)
                  + self.c)
        self._apsis = Apsis(self)

    def apsis(self, e, p):
        """Calculate the distance of closest approach (apsis) in a colllision.

        As initial condition, the larger of the apsis estimate from the table 
        and the impact parameter is used.

        Parameters:
            e (float): energy of projectile before the collision (ENORM)
            p (float): impact parameter (RNORM)

        Returns:
            (float): Estimated apsis of the collision (RNORM)
            (int): Number of iterations used to converge the apsis
        """
        return self._apsis.call(e, p, self)

    def call(self, r):
        """Calculate the NLHlin screening function and its derivative.

        Parameters:
            r (float or ndarray): Distance (RNORM)

        Returns:
            (ndarray): NLHlin screening function at distance r
            (ndarray): Derivative of NLHlin screening function at distance r
                (1/RNORM)
        """
        exp0 = np.exp(-self.b[0]*r)
        exp1 = np.exp(-self.b[1]*r)
        exp2 = np.exp(-self.b[2]*r)

        screen = (self.a[0]*exp0 + self.a[1]*exp1 
                        + self.a[2]*exp2 + self.c - self.d*r/self.rmax)
        dscreen = (- self.ab[0]*exp0 - self.ab[1]*exp1 
                        - self.ab[2]*exp2 - self.d/self.rmax)
        
        mask = r < self.rmax    # should always be true except for testing
        screen = np.where(mask, screen, 0.0)
        dscreen = np.where(mask, dscreen, 0.0)
        return screen, dscreen


def read_coefs():
    """Read NLHlin screening coefficients from the data file.

    Returns:
        (recarray): A record array containing the coefficients.
    """
    fname = os.path.join(os.path.dirname(__file__), "dmol_coeffs_rmax.dat")
    if not os.path.exists(fname):
        print(f"NLHlin_screen: Coefficients file {fname} not found")
        sys.exit()
    
    coef_rows = []
    with open(fname, "r") as f:
        for line in f:
            if line[0] == "#":
                continue
            coefs = line.split()[:-1]   # exclude "error" column
            coef_rows.append(tuple([float(c) for c in coefs]))

    # TODO: Remove duplicate definition of SCATTER_PARAMS_DTYPE
    NLHLIN_COEFS_DTYPE = np.dtype([
        ("z1", np.uint32),
        ("z2", np.uint32),
        ("a1", np.float64),
        ("b1", np.float64),
        ("a2", np.float64),
        ("b2", np.float64),
        ("a3", np.float64),
        ("b3", np.float64),
        ("rmax", np.float64),
    ], align=True)

    return np.array(coef_rows, dtype=NLHLIN_COEFS_DTYPE).view(np.recarray)


def post_plot(p1, p2, Z2):
    """Do post-plot setup for NLHlin screening function plots.
    
    Parameters:
        p1 (float): exponent for atomic number scaling
        p2 (float): exponent for atomic number scaling
        Z2 (int or None): atomic number of second atom, or None for Z2=Z1
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                label=r"atomic number Z$_1$", 
                ax=plt.gca(),
                ticks=ticks)
    plt.yscale("log")
    if (p1, p2) == (0, 0):
        plt.xlim(0.0, 3.0)
    else:
        plt.xlim(0.0, 30.0)
    plt.ylim(1e-4, 1.0)
    if (p1, p2) == (0, 0):
        text = ""
    else:
        text = r"a$_\mathrm{I}$=0.4685$\rm\AA$/"
        if p1 == 1:
            text += fr"(Z$_1$+Z$_2$)"
        else:
            text += fr"(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)"
        if p2 != 1:
            text += fr"$^{{{p2:.2f}}}$"
        text += "\n"
    if Z2 is None:
        text += f" Z$_2$=Z$_1$"
    else:
        text += f" Z$_2$={Z2}"
    plt.text(0.95, 0.95, text, 
             horizontalalignment="right", verticalalignment="top",
             transform=plt.gca().transAxes, fontsize="medium")
    if (p1, p2) == (0, 0):
        plt.xlabel(r"distance ($\rm\AA$)")
    else:
        plt.xlabel("reduced distance")
    plt.ylabel("NLHlin screening function")
    plt.tight_layout()


def plot_screen(p1, p2, z2=None):
    """Plot the NLHlin screening function for testing purposes.
    
    Parameters:
        p1 (float): exponent for atomic number scaling
        p2 (float): exponent for atomic number scaling
        z2 (int or None): atomic number of second atom, or None for Z2=Z1
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    rmax_A  = 3.0

    plt.figure()
    #cmap = mpl.cm.get_cmap("jet", 92)
    cmap = mpl.cm.get_cmap("viridis", 92)
    plt.rcParams.update({"font.size": 14})
    #plt.gca().set_facecolor("darkgray")

    coefs = read_coefs()
    for Z1 in range(1, 93):
        if z2 is None:
            Z2 = Z1
        else:
            Z2 = z2
        if (p1, p2) == (0, 0):
            rnorm = 1.0
        else:
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
        rmax = rmax_A / rnorm


        r = np.linspace(0.0, rmax, 101)
        screen, _ = NLHlin_screen(Z1, Z2, coefs, rnorm).call(r)
        plt.plot(r, screen, color=cmap((Z1-1)/92), zorder=Z1)

    if (p1, p2) == (0.23, 1):
        screen, _ = ZBL_screen().call(r)
        plt.plot(r, screen, "k--", label="ZBL", zorder=100)
        plt.legend(loc="right")
    # TODO uncomment
    # elif (p1, p2) == (1/2, 2/3):
    #     screen, _ = KrC_screen(r)
    #     plt.plot(r, screen, "k--", label="KrC", zorder=100)
    #     plt.legend(loc="right")

    post_plot(p1, p2, Z2=z2)

    plt.show()


def plot_ZBLscreen(p1, p2, z2=None):
    """Plot the ZBL screening function for comparison.
    
    Parameters:
        p1 (float): exponent for atomic number scaling
        p2 (float): exponent for atomic number scaling
        z2 (int or None): atomic number of second atom, or None for Z2=Z1
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    rmax_A  = 3.0

    plt.figure()
    #cmap = mpl.cm.get_cmap("jet", 92)
    cmap = mpl.cm.get_cmap("viridis", 92)
    plt.rcParams.update({"font.size": 14})
    #plt.gca().set_facecolor("darkgray")

    for Z1 in range(1, 93):
        if z2 is None:
            Z2 = Z1
        else:
            Z2 = z2
        a_ZBL = 0.4685 / (Z1**0.23 + Z2**0.23)
        if (p1, p2) == (0, 0):
            rnorm = 1.0
        else:
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2

        r = np.linspace(0.0, rmax_A, 101)
        screen, _ = ZBL_screen().call(r/a_ZBL)
        plt.plot(r/rnorm, screen, color=cmap((Z1-1)/92), zorder=Z1)

    post_plot(p1, p2, Z2=r"Z$_1$")
    plt.ylabel("ZBL screening function")

    plt.show()


if __name__ == "__main__":
    from zbl import ZBL_screen
    # from krc import KrC_screen    # TODO uncomment
    Z2 = 29
    #plot_screen(p1=0, p2=0)      # unscaled
    #plot_screen(p1=1/2, p2=2/3)  # Firsov
    #plot_screen(p1=2/3, p2=1/2)  # Lindhard
    #plot_screen(p1=0.23, p2=1)   # ZBL
    #plot_screen(p1=1/4, p2=1, z2=Z2)    # Suggested by M. Hou (AI generated, true?)
    #plot_screen(p1=1/2, p2=1/2)  # New suggestion
    #plot_screen(p1=1, p2=1/4, z2=Z2)    # Alternate suggestion   

    #plot_ZBLscreen(p1=0, p2=0)   # unscaled ZBL
    plot_ZBLscreen(p1=0.23, p2=1) # ZBL
