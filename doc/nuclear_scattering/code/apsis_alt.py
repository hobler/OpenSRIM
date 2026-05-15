"""Calculate the apsis of collision for arbitrary screening function.

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

        # determine rmax so that the screening function and its derivative
        # are non-zero (not guaranteed to succeed, but ok for ZBL)
        rmax = min(100.0, screen_fun.rmax)

        # k2/R part
        def fun2(r):
            screen, dscreen = screen_fun(r)
            return screen + r * dscreen
        print("fun2 at R=0:", fun2(0.0))
        print("fun2 at Rmax:", fun2(rmax))
        r_touch2 = brentq(fun2, 0.0, rmax)
        self.k2 = r_touch2 * screen_fun(r_touch2)[0]

        # k3/R^3 part
        def fun3(r):
            screen, dscreen = screen_fun(r)
            return screen + 1/3 * r * dscreen
        print("fun3 at R=0:", fun3(0.0))
        print("fun3 at Rmax:", fun3(rmax))
        r_touch3 = brentq(fun3, 0.0, rmax)
        self.k3 = r_touch3**3 * screen_fun(r_touch3)[0]
        
        # 1 - k1*R part
        self.k1 = 1 / (4*self.k2)

        # k4*(Rmax-R) part
        if screen_fun.rmax == np.inf:
            self.r34 = r_touch3
        else:
            self.r34 = 0.75 * rmax
            self.k4 = 3*self.k3 / self.r34**4

        print(f"Touching point for 1/R part: r={r_touch2:.3e}")
        print(f"Touching point for 1/R^3 part: r={r_touch3:.3e}")
        print(f"Coefficients: k1={self.k1:.3e}, k2={self.k2:.3e}, k3={self.k3:.3e}")

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
        
        psq = p**2
        # Initial condition
        if self.screen_fun.rmax is np.inf:  # Use TRIM85 algorithm
            r0 = max(1e-10, p)
            r0_try = -2.7 * np.log(e*r0)
            if r0_try > p:
                r0_try = -2.7 * np.log(e*r0_try)
                if r0_try > p:
                    r0 = r0_try
            done = r0 > self.r34
#            done = True
        else:
            if psq > self.r34**2 - self.k3/(e*self.r34**2):
                a = e + self.k4
                b = - self.k4 * self.screen_fun.rmax
                c = - e * psq
                r0 = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)
                done = True
            else:
                done = False

        if not done:
            r0sq = psq + self.k2/e
            if r0sq > self.k3 / self.k2:
                r0 = np.sqrt(psq/2 + np.sqrt(psq**2/4 + self.k3/e))
            elif r0sq >= self.k2 / self.k1:
                r0 = np.sqrt(r0sq)
            else:
                r0 = (1 + np.sqrt(1 + 4*e*(e+self.k1)*psq)) / (2*(e+self.k1))

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

    Z1 = 33
    Z2 = 14
    rnorm = 0.4685 / (np.sqrt(np.sqrt(Z1)) + np.sqrt(np.sqrt(Z2)))
    screen_fun = NLHlin_screen(Z1, Z2, rnorm)
    
    plot_iteration_counts(screen_fun, Z1, Z2)

