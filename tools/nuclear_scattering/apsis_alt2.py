"""Calculate the apsis of collision for arbitrary screening function.

Uses initial conditions based on the tangent points at r=0 and r=rmax.

Available classes:
- Apsis: Main class for calculating apsis of collision.
"""
import numpy as np
from scipy.optimize import brentq


class Apsis:
    """Calculate apsis of collision for arbitrary screening function.
    """
    def __init__(self, screen_fun):
        """Initialize apsis calculation with given screening function.

        In particular, the coefficients for the piecewise approximation to the
        screening function are calculated.

        Parameters:
            screen_fun (callable): Screening function
        """
        self.screen_fun = screen_fun

        _, self.dscreen_0 = screen_fun(0.0)
        _, self.dscreen_rmax = screen_fun(screen_fun.rmax)

    def __call__(self, e, p):
        """Calculate the distance of closest approach (apsis) in a colllision.

        The initial condition is calculated using the piecewise approximation
        to the screening function.

        Parameters:
            e (float): energy of projectile before the collision (ENORM)
            p (float): impact parameter (RNORM)

        Returns:
            (float): Estimated apsis of the collision (RNORM)
            (int): Number of iterations used to converge the apsis
        """
        if p >= self.screen_fun.rmax:
            return p, 0

        # Initial condition
        a = 1.0 - self.dscreen_0 / e
        b = - 1 / e
        c = - p**2
        r1 = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)

        if self.screen_fun.rmax == np.inf:
            r2 = p
        else:
            a = 1.0 - self.dscreen_rmax / e
            b = self.screen_fun.rmax / e * self.dscreen_rmax
            c = - p**2
            r2 = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)

        r0 = max(r1, r2)

        # Newton iteration
        delta_r0 = np.inf

        def fun(r):
            screen, dscreen = self.screen_fun(r)
            return r - screen/e - p**2/r, 1 - dscreen/e + p**2/r**2
        
        count = 0
        while abs(delta_r0) > 1e-3 * r0:
            f, df = fun(r0)
            delta_r0 = - f / df
            r0 += delta_r0
            count += 1

        return r0, count


def plot_iteration_counts(screen_fun, Z1=None, Z2=None):
    """Plot the number of iterations required for apsis calculation.

    The plot shows the number of iterations needed to converge the apsis
    as a function of energy and impact parameter to a relative accuracy of
    1e-3 (If quadratic convergence has been reached, the relative accuracy
    actually is 1e-6).

    Parameters:
        screen_fun (callable): Screening function taking distance r as argument
        Z1 (int or None): atomic number of first atom, or None if not needed
        Z2 (int or None): atomic number of second atom, or None if not needed
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    apsis = Apsis(screen_fun)

    e_values = np.logspace(-6, 2, 9)
    p_values = 10.0 * np.linspace(0, 2, 11)**2  * 0.75
    e_idx = np.arange(len(e_values)+1) - 0.5
    p_idx = np.arange(len(p_values)+1) - 0.5

    n_iter = np.empty((len(p_values), len(e_values)), dtype=int)
    for i, e in enumerate(e_values):
        for j, p in enumerate(p_values):
            r0, niter = apsis(e, p)
            n_iter[j, i] = niter

    plt.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(8,6))
    ax = plt.gca()
    nmax = 5
    bounds = np.linspace(-0.5, nmax+0.5, nmax+2)
    ticks = np.linspace(0, nmax, nmax+1)
    cmap = mpl.cm.viridis   # viridis, plasma, jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    c = ax.pcolormesh(e_idx - 0.5, p_idx - 0.5, n_iter,
                      shading='auto', cmap=cmap, norm=norm)
    cbar = fig.colorbar(c, ax=ax, ticks=ticks)
    cbar.set_label(r"number of iterations $n$")

    ax.set_xticks(e_idx[:-1])
    ax.set_xticklabels(np.asarray(np.log10(e_values), dtype=int))
    ax.set_xlabel(r"$^{10}$log(energy $\varepsilon$)")
    ax.set_yticks(p_idx[:-1])
    ax.set_yticklabels(np.asarray(10*p_values, dtype=int) / 10)
    ax.set_ylabel(r"impact parameter $P$")

    if type(screen_fun) is ZBL_screen:
        ax.set_title("Universal ZBL potential", fontsize="small")
    elif type(screen_fun) is NLHlin_screen:
        ax.set_title(fr"NLHlin potential, Z$_1$={Z1}, Z$_2$={Z2}", 
                     fontsize="small")
    plt.tight_layout()
    
    plt.show()

if __name__ == "__main__":
    from zbl import ZBL_screen
    from nlhlin import NLHlin_screen

    #screen_fun = ZBL_screen()

    Z1 = 74
    Z2 = 74
    rnorm = 0.4685 / (np.sqrt(np.sqrt(Z1)) + np.sqrt(np.sqrt(Z2)))
    screen_fun = NLHlin_screen(Z1, Z2, rnorm)
    
    plot_iteration_counts(screen_fun, Z1, Z2)

