"""Define the NLH screening function.

Available functions:
    NLH_screen: Callable object for the NLH screening function. Also 
        defines a method impulse integral for calculating the integral that 
        appears in the impulse approximation.
    
Moreover, there are plotting functions for testing and visualization:
    plot_screen: plot NLHlin screening function for various atomic numbers.
    plot_ZBLscreen: plot ZBL screening function for various atomic numbers.
"""
import os, sys
import numpy as np
from scipy import special
from apsis import Apsis
from utils import atom, ask_if_save


class NLH_screen:
    """Defines the NLH screening function.
    """
    def __init__(self, Z1, Z2, rnorm=None):
        """Setup NLH screening function for given atomic numbers.

        Parameters:
            Z1: atomic number of atom 1
            Z2: atomic number of atom 2
            rnorm: screening length (A), None for the default value of
                0.4685 / sqrt(sqrt(Z1) + sqrt(Z2))
        """
        self.name = "NLH"
        self.Z1 = Z1
        self.Z2 = Z2
        
        fname = os.path.join(os.path.dirname(__file__), 
                             '../../data/NLH/nlh_coeffs.dat')
        if not os.path.exists(fname):
            print(f'NLH_screen: Coefficients file {fname} not found')
            sys.exit()

        with open(fname) as f:
            for line in f:
                if line[0] == '#':
                    continue
                z1, z2, a1, b1, a2, b2, a3, b3, *error = line.split()
                z1 = int(z1)
                z2 = int(z2)
                if min(z1,z2) == min(Z1, Z2) and max(z1, z2) == max(Z1, Z2):
                    break

        self.a = (float(a1), float(a2), float(a3))
        if rnorm is None:
            rnorm = 0.4685 / np.sqrt(np.sqrt(Z1) + np.sqrt(Z2))
        self.b = (float(b1)*rnorm, float(b2)*rnorm, float(b3)*rnorm)
        self.ab = (self.a[0]*self.b[0], self.a[1]*self.b[1], 
                   self.a[2]*self.b[2])
        self.rnorm = rnorm

        self.rmax = np.inf
        
        self.apsis = Apsis(self)

    def __call__(self, r):
        """Calculate the NLH screening function and its derivative.

        Parameters:
            r (float or ndarray): Distance (RNORM)

        Returns:
            (ndarray): NLH screening function at distance r
            (ndarray): Derivative of NLH screening function at distance r
                (1/RNORM)
        """
        r = np.asarray(r, dtype=float)
        
        exp0 = np.exp(-self.b[0]*r)
        exp1 = np.exp(-self.b[1]*r)
        exp2 = np.exp(-self.b[2]*r)

        screen = self.a[0]*exp0 + self.a[1]*exp1 + self.a[2]*exp2
        dscreen = - self.ab[0]*exp0 - self.ab[1]*exp1 - self.ab[2]*exp2
            
        return screen, dscreen

    def impulse_integral(self, p):
        """Evaluate the integral that appears in the impulse approximation.
        
        This integral is defined as one half of the integral of 

             d   Phi(r)
            -- * ------
            dp     r
        
        along a straight line which passes by the center of the potential at 
        a distance p.

        The calculation uses the modified Bessel function of the second kind 
        and order 1 (scipy.special.kn).
        
        Parameters:
            p (float): impact parameter (RNORM)

        Returns:
            (float): value of the integral
        """
        k0 = special.kn(1, self.b[0] * p)
        k1 = special.kn(1, self.b[1] * p)
        k2 = special.kn(1, self.b[2] * p)
        k3 = special.kn(1, self.b[3] * p)

        integral = self.ab[0]*k0 + self.ab[1]*k1 + self.ab[2]*k2 + self.ab[3]*k3

        return integral


def post_plot(p1, p2, Z2, xmax=None, ymin=None):
    """Do post-plot setup for NLH screening function plots.
    
    Parameters:
        p1 (float): exponent for atomic number scaling
        p2 (float): exponent for atomic number scaling
        Z2 (int or None): atomic number of second atom, or None for Z2=Z1
        xmax (float or None): maximum x value for plot
        ymin (float or None): minimum y value for plot
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    plt.plot(0, 1, 'k:', label='V<1eV')
    plt.legend(loc=(0.635, 0.6),frameon=False)

    plt.yscale('log')
    if xmax is None:
         xmax = 3.0 if (p1, p2) == (0, 0) else 30.0
    if ymin is None:
        ymin = 1e-4
    plt.xlim(0.0, xmax)
    plt.ylim(ymin, 1.0)

    if (p1, p2) == (0, 0):
        text = ''
    else:
        text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
        if p1 == 1:
            text += fr'(Z$_1$+Z$_2$)'
        else:
            text += fr'(Z$_1^{{{round(p1, 2)}}}$+Z$_2^{{{round(p1, 2)}}}$)'
        if p2 != 1:
            text += fr'$^{{{round(p2, 2)}}}$'
        text += '\n'
    if Z2 is None:
        text += f' Z$_2$=Z$_1$'
    else:
        text += f' Z$_2$={Z2}'
    plt.text(0.95, 0.95, text, 
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes, fontsize='medium')

    if (p1, p2) == (0, 0):
        plt.xlabel(r'distance ($\rm\AA$)')
    else:
        plt.xlabel('reduced distance')
    plt.ylabel('NLH screening function')

    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    #bounds = np.linspace(1, 10, 10)
    #cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                 ax=plt.gca(),
                 label=r'atomic number Z$_1$', 
                 ticks=ticks)

    plt.tight_layout()


def plot_screen(p1, p2, z2=None, xmax=None, ymin=None):
    """Plot the NLH screening function for testing purposes.
    
    Parameters:
        p1 (float): exponent for atomic number scaling
        p2 (float): exponent for atomic number scaling
        z2 (int or None): atomic number of second atom, or None for Z2=Z1
        xmax (float or None): maximum x value for plot
        ymin (float or None): minimum y value for plot
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    rmax_A  = 4.0

    fig = plt.figure()
    #cmap = plt.get_cmap('jet', 92)
    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})
    #plt.gca().set_facecolor('darkgray')

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
        screen, _ = NLH_screen(Z1, Z2, rnorm)(r)
        plt.plot(r, screen, ':', color=cmap((Z1-1)/92), zorder=Z1)
        mask = (14.4 * Z1 * Z2 / np.maximum(1e-10, r*rnorm) * screen > 1)
        plt.plot(r[mask], screen[mask], color=cmap((Z1-1)/92), zorder=Z1)

    if (p1, p2) == (0.23, 1):
        screen, _ = ZBL_screen()(r)
        plt.plot(r, screen, 'k--', label='ZBL', zorder=100)
        #plt.legend(loc='right')
    elif (p1, p2) == (1/2, 2/3):
        screen, _ = KrC_screen()(r)
        plt.plot(r, screen, 'k--', label='KrC', zorder=100)
        #plt.legend(loc='right')

    post_plot(p1, p2, Z2=z2, xmax=xmax, ymin=ymin)
    plt.show()

    if (p1, p2) == (0, 0):
        fname = f"figs/nlh_unscaled.pdf"
    elif z2 is None:
        fname = f"figs/nlh_p{round(p1, 2)}_p{round(p2, 2)}.pdf"
    else:
        fname = f"figs/nlh_p{round(p1, 2)}_p{round(p2, 2)}_{atom[Z2]}.pdf"
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))


def plot_ZBLscreen(p1, p2, z2=None, xmax=None, ymin=None):
    """Plot the ZBL screening function for comparison.
    
    Parameters:
        p1 (float): exponent for atomic number scaling
        p2 (float): exponent for atomic number scaling
        z2 (int or None): atomic number of second atom, or None for Z2=Z1
        xmax (float or None): maximum x value for plot
        ymin (float or None): minimum y value for plot
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    rmax_A  = 3.0

    fig = plt.figure()
    #cmap = plt.get_cmap('jet', 10)
    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})
    #plt.gca().set_facecolor('darkgray')

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
        screen, _ = ZBL_screen()(r/a_ZBL)
        plt.plot(r/rnorm, screen, ':', color=cmap((Z1-1)/92), zorder=Z1)
        mask = (14.4 * Z1 * Z2 / np.maximum(1e-10, r*rnorm) * screen > 1)
        plt.plot(r[mask]/rnorm, screen[mask], color=cmap((Z1-1)/92), zorder=Z1)

    post_plot(p1, p2, Z2=r'Z$_1$', xmax=xmax, ymin=ymin)
    plt.ylabel('ZBL screening function')
    plt.show()

    if (p1, p2) == (0, 0):
        fname = f"figs/zbl_unscaled.pdf"
    else:
        fname = f"figs/zbl_p{round(p1, 2)}_p{round(p2, 2)}.pdf"
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))


if __name__ == "__main__":
    from zbl import ZBL_screen
    from krc import KrC_screen
    Z2 = 2
    #plot_screen(p1=0, p2=0)      # unscaled
    #plot_screen(p1=1/2, p2=2/3, xmax=14.0, ymin=0.01)  # Firsov
    #plot_screen(p1=2/3, p2=1/2)  # Lindhard
    #plot_screen(p1=0.23, p2=1, xmax=14.0, ymin=0.01)   # ZBL
    #plot_screen(p1=1/4, p2=1)    # Suggested by M. Hou (AI generated, true?)
    #plot_screen(p1=1/2, p2=1/2, xmax=10.0, ymin=0.01)  # New suggestion
    #plot_screen(p1=1, p2=1/4, z2=Z2, xmax=10.0, ymin=0.01)    # Alternative suggestion   
    for Z2 in range(1, 93):
        plot_screen(p1=1/2, p2=1/2, z2=Z2, xmax=30.0, ymin=1e-5)  # New suggestion
    #plot_screen(p1=1/4, p2=1, z2=Z2, xmax=10.0, ymin=0.01)    # Alternative suggestion   

    #plot_ZBLscreen(p1=0, p2=0)   # unscaled ZBL
    #plot_ZBLscreen(p1=0.23, p2=1) # ZBL
