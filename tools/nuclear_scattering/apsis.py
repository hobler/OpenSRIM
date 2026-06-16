"""Calculate the apsis of collision for arbitrary screening function.

Available functions:
- apsis_setup: Setup apsis table for head-on collisions.
- calc_apsis: Calculate the distance of closest approach (apsis) in a collision.
"""
import os
import numpy as np
from table1d import Table1D
from utils import atom, ask_if_save


class Apsis:
    """Calculate apsis of collision for arbitrary screening function.
    """
    def __init__(self, screen_fun):
        """Initialize apsis calculation with given screening function.

        In particular, this sets up the apsis table for head-on collisions.
        In the table, the apsis is stored for a range of reduced energies
        (lowest energy/largest apsis first, highest energy/smallest apsis last).

        Parameters:
            screen_fun (callable): Screening function
        """
        self.screen_fun = screen_fun

        emax = 1e4
        emin = 1e-8

        energies = []
        apses = []
        dapses_de = []

        # Initial conditions for the maximum reduced energy
        e = emax
        r0 = 1 / emax
        delta_r0 = np.inf

        # Iterating down the reduced energy by factors of 2
        while e >= emin/2:
            while abs(delta_r0) > 1e-6 * r0:
                screen, dscreen = screen_fun(r0)
                f = r0 - screen/e
                df = 1 - dscreen/e
                delta_r0 = - f / df
                r0 += delta_r0
            dr0_de = - r0**2 / (screen - r0 * dscreen)
            energies.append(e)
            apses.append(r0)
            dapses_de.append(dr0_de)

            # Initial conditions for the next lower reduced energy
            # (can use old r0 as initial condition for new r0)
            e /= 2
            delta_r0 = np.inf

        energies = np.array(energies)
        apses = np.ravel(np.array(apses))
        dapses_de = np.ravel(np.array(dapses_de))

        self.apsis_headon_table = Table1D(energies[::-1],
                                          apses[::-1],
                                          dapses_de[::-1],
                                          powerof2=True)

    def __call__(self, e, p, relerr=1e-3):
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
        if p >= self.screen_fun.rmax:
            return p, 0
        
        if e < self.apsis_headon_table.x[0]:
            r0 = self.apsis_headon_table.y[0]
        elif e > self.apsis_headon_table.x[-1]:
            r0 = 1 / e
        else:
            r0 = self.apsis_headon_table.interpolate(e)
        
        r0 = max(r0, p)
        delta_r0 = np.inf

        def fun(r):
            screen, dscreen = self.screen_fun(r)
            return r - screen/e - p**2/r, 1 - dscreen/e + p**2/r**2
        
        count = 0
        while abs(delta_r0) > relerr * r0:
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
    p_values = 10.0 * np.linspace(0, 2, 11)**2
    if type(screen_fun) is NLHlin_screen:
        p_values *= 0.75
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
    bounds = np.linspace(-0.5, 4.5, 6)
    ticks = np.linspace(0, 4, 5)
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
        fname = "figs/iter_apsis_zbl.pdf"
    elif type(screen_fun) is NLHlin_screen:
        ax.set_title(fr"NLHlin potential, Z$_1$={Z1}, Z$_2$={Z2}", 
                     fontsize="small")
        fname = f"figs/iter_apsis_nlhlin_{atom[Z1]}_{atom[Z2]}.pdf"
    plt.tight_layout()
    plt.show()

    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))


if __name__ == "__main__":
    from zbl import ZBL_screen
    from nlhlin import NLHlin_screen

    #screen_fun = ZBL_screen()

    Z1 = 33
    Z2 = 14
    rnorm = 0.4685 / (np.sqrt(np.sqrt(Z1)) + np.sqrt(np.sqrt(Z2)))
    screen_fun = NLHlin_screen(Z1, Z2, rnorm)
    
    plot_iteration_counts(screen_fun, Z1, Z2)

