"""Define the NLHlin screening function.

Available functions:
    NLHlin_screen: Callable object for the NLHlin screening function.
    
Moreover, there are plotting functions for testing and visualization:
    plot_screen: plot NLHlin screening function for various atomic numbers.
    plot_ZBLscreen: plot ZBL screening function for various atomic numbers.
"""
import os, sys
import numpy as np
from apsis import Apsis


class NLHlin_screen:
    """Defines the NLHlin screening function.
    """
    def __init__(self, Z1, Z2, rnorm=None):
        """Setup NLHlin screening function for given atomic numbers.

        Parameters:
            Z1: atomic number of atom 1
            Z2: atomic number of atom 2
            rnorm: screening length (A), None for the default value of
                0.4685 / sqrt(sqrt(Z1) + sqrt(Z2))
        """
        self.Z1 = Z1
        self.Z2 = Z2
        
        fname = os.path.join(os.path.dirname(__file__), 'dmol_coeffs_rmax.dat')
        if not os.path.exists(fname):
            print(f'NLHlin_screen: Coefficients file {fname} not found')
            sys.exit()

        with open(fname) as f:
            for line in f:
                if line[0] == '#':
                    continue
                z1, z2, a1, b1, a2, b2, a3, b3, rmax, error = line.split()
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

        self.rmax = float(rmax) / rnorm
        self.c = 1 - self.a[0] - self.a[1] - self.a[2]
        self.d = (self.a[0]*np.exp(-self.b[0]*self.rmax) 
                  + self.a[1]*np.exp(-self.b[1]*self.rmax) 
                  + self.a[2]*np.exp(-self.b[2]*self.rmax)
                  + self.c)
        self.apsis = Apsis(self)


    def __call__(self, r):
        """Calculate the NLHlin screening function and its derivative.

        Parameters:
            r (float or ndarray): Distance (RNORM)

        Returns:
            (ndarray): NLHlin screening function at distance r
            (ndarray): Derivative of NLHlin screening function at distance r
                (1/RNORM)
        """
        r = np.asarray(r, dtype=float)
        
        if np.all(r < self.rmax):   # should always be true except for testing
            exp0 = np.exp(-self.b[0]*r)
            exp1 = np.exp(-self.b[1]*r)
            exp2 = np.exp(-self.b[2]*r)

            screen = (self.a[0]*exp0 + self.a[1]*exp1 
                            + self.a[2]*exp2 + self.c - self.d*r/self.rmax)
            dscreen = (- self.ab[0]*exp0 - self.ab[1]*exp1 
                            - self.ab[2]*exp2 - self.d/self.rmax)
        else:
            mask = np.asarray(r < self.rmax)

            exp0 = np.exp(-self.b[0]*r[mask])
            exp1 = np.exp(-self.b[1]*r[mask])
            exp2 = np.exp(-self.b[2]*r[mask])

            screen = np.zeros_like(r)
            screen[mask] = (self.a[0]*exp0 + self.a[1]*exp1 
                            + self.a[2]*exp2 + self.c 
                            - self.d*r[mask]/self.rmax)
            dscreen = np.zeros_like(r)
            dscreen[mask] = (- self.ab[0]*exp0 - self.ab[1]*exp1 
                            - self.ab[2]*exp2 - self.d/self.rmax)
            
        return screen, dscreen

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
                label=r'atomic number Z$_1$', 
                ticks=ticks)
    plt.yscale('log')
    if (p1, p2) == (0, 0):
        plt.xlim(0.0, 3.0)
    else:
        plt.xlim(0.0, 30.0)
    plt.ylim(1e-4, 1.0)
    if (p1, p2) == (0, 0):
        text = ''
    else:
        text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
        if p1 == 1:
            text += fr'(Z$_1$+Z$_2$)'
        else:
            text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
        if p2 != 1:
            text += fr'$^{{{p2:.2f}}}$'
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
    plt.ylabel('NLHlin screening function')
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
    #cmap = mpl.cm.get_cmap('jet', 92)
    cmap = mpl.cm.get_cmap('viridis', 92)
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
        screen, _ = NLHlin_screen(Z1, Z2, rnorm)(r)
        plt.plot(r, screen, color=cmap((Z1-1)/92), zorder=Z1)

    if (p1, p2) == (0.23, 1):
        screen, _ = ZBL_screen()(r)
        plt.plot(r, screen, 'k--', label='ZBL', zorder=100)
        plt.legend(loc='right')
    elif (p1, p2) == (1/2, 2/3):
        screen, _ = KrC_screen(r)
        plt.plot(r, screen, 'k--', label='KrC', zorder=100)
        plt.legend(loc='right')

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
    #cmap = mpl.cm.get_cmap('jet', 92)
    cmap = mpl.cm.get_cmap('viridis', 92)
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
        plt.plot(r/rnorm, screen, color=cmap((Z1-1)/92), zorder=Z1)

    post_plot(p1, p2, Z2=r'Z$_1$')
    plt.ylabel('ZBL screening function')

    plt.show()


if __name__ == "__main__":
    from zbl import ZBL_screen
    from krc import KrC_screen
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
