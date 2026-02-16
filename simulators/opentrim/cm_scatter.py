"""Characterize a scattering event in the center-of-mass system.

Calculate scattering angle and time integral for given interatomic potential,
energy, and impact parameter, using Gauss-Legendre quadrature.
"""
# import os
# os.environ["NUMBA_DISABLE_JIT"] = "1"
from numba import jit
import numpy as np
from scipy.special import roots_legendre


def setup(n_absc):
    """Setup Gauss-Legendre abscissae and weights.

    Parameters:
        n_absc (int): Number of abscissae.
    """
    global ROOTS_LEGENDRE

    ROOTS_LEGENDRE = roots_legendre(n_absc)

@jit
def calc_phi_chi(u, r0, screen_fun):
    """Calculate the screening function and chi at the given u values.

    phi is evaluated at r0/(1-u**2).
    chi is the difference of phi and phi0 (=phi(r0)) divided by u^2, with a 
    Taylor expansion used for small u to avoid numerical issues.

    Parameters:
        u (array-like): Integration variable.
        r0 (float): Distance of closest approach (RNORM).
        screen_fun (object): Object with `.call()` to calculate the screening function
            for given distance r (RNORM).
    Returns:
        (float): Value of chi function.
    """
    # u = np.asarray(u)
    phi0, dphi0 = screen_fun.call(r0)
    phi, _ = screen_fun.call(r0/(1-u**2))
    chi = np.where(u < 3e-4, phi0 - r0*dphi0, (phi0 - phi*(1-u**2)) / u**2)
    return phi, chi

@jit
def scatter_integrals(e, p, screen_fun):
    """Calculate scattering angle and time integral.

    The calculation uses Gauss-Legendre quadrature with a fixed number of
    abscissae given by ROOTS_LEGENDRE. The abscrissae and weights must have 
    been set up before calling this function using setup().

    Parameters:
        e (float): Reduced energy.
        p (float): Reduced impact parameter.
        screen_fun (object): Object with `.call()` to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Scattering angle (rad)
        (float): Time integral (RNORM).
    """
    if p >= screen_fun.rmax:
        return 0.0, 0.0
    elif p == 0.0:
        return np.pi, 0.0
    
    r0, _ = screen_fun.apsis(e, p)

    def integrands(u):
        phi, chi = calc_phi_chi(u, r0, screen_fun)
        rho = r0 / (e*p**2)
        g = np.sqrt(rho*chi + (2-u**2))
        integrand_theta = 1 / g
        integrand_tau = (1 + rho * phi/(1-u**2)) / (g * (1 + u*p/r0*g))
        return integrand_theta, integrand_tau
    
    def integrand_two_arccos(u):
        return 4 / np.sqrt(2 - u**2)

    rmax = screen_fun.rmax
    umax = np.sqrt(1 - r0 / rmax)
    u_vals, weights = ROOTS_LEGENDRE
    u_vals = 0.5 * umax * (u_vals + 1)
    weights = 0.5 * umax * weights
    integrand_theta_vals, integrand_tau_vals = integrands(u_vals)
    two_arccos_num = np.sum(weights * integrand_two_arccos(u_vals))
    if rmax == np.inf:
        theta = np.pi * (1 
                - 4/two_arccos_num * np.sum(weights * integrand_theta_vals))
        tau = (r0 - 2 * p * np.sum(weights * integrand_tau_vals))
    else:    
        theta = 2*np.arccos(p/rmax) * (1 
                - 4/two_arccos_num * np.sum(weights * integrand_theta_vals))
        tau = (r0 - (rmax - np.sqrt(rmax**2 - p**2)) 
            - 2 * p * np.sum(weights * integrand_tau_vals))

    return theta, tau



def plot_chi(r0_vals, screen_fun):
    """Plot the chi function for scattering integrals.

    chi is part of the integrand for scattering angle and time integral 
    calculations. The code here repeats the calculation of chi from within
    scatter_integrals for plotting purposes.

    Parameters:
        u (float): Integration variable.
        r0_vals (array-like): Distance of closest approach (RNORM).
        screen_fun (object): Object with `.call()` to calculate the screening function
            for given distance r (RNORM).
    Returns:
        (float): Value of chi function.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 14})

    u = np.linspace(0.001, 0.999, 999)

    for r0 in r0_vals:
        phi, chi = calc_phi_chi(u, r0, screen_fun)
        plt.plot(u, chi, label=fr'$R_0$={r0}')

    plt.xlabel('u')
    plt.ylabel(r'$\chi$(u)')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()
    title = "Potential"
    if isinstance(screen_fun, ZBL_screen):
        title = 'ZBL potential'
    elif isinstance(screen_fun, NLHlin_screen):
        title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                    fr'Z$_2$={screen_fun.Z2}')
    plt.title(title, fontsize='medium')
    plt.show()


def plot_chi_near_zero(r0_vals, screen_fun):
    """Plot the chi function for scattering integrals.

    chi is part of the integrand for scattering angle and time integral 
    calculations. The code here repeats the calculation of chi from within
    scatter_integrals for plotting purposes.

    Parameters:
        u (float): Integration variable.
        r0_vals (array-like): Distance of closest approach (RNORM).
        screen_fun (object): Object with `.call()` to calculate the screening function
            for given distance r (RNORM).
    Returns:
        (float): Value of chi function.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 14})

    u = np.logspace(-10, 0, 101)

    for r0 in r0_vals:
        phi0, dphi0 = screen_fun.call(r0)
        phi, _ = screen_fun.call(r0/(1-u**2))
        #chi = (phi0 - phi*(1-u**2)) / u**2
        _, chi = calc_phi_chi(u, r0, screen_fun)
        plt.plot(u, np.full_like(u, phi0) - r0*dphi0, 'k--')
        plt.plot(u, chi, label=fr'$R_0$={r0}')

    plt.xscale('log')
    plt.xlabel('u')
    plt.ylabel(r'$\chi$(u)')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()
    title = "Potential"
    if isinstance(screen_fun, ZBL_screen):
        title = 'ZBL potential'
    elif isinstance(screen_fun, NLHlin_screen):
        title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                    fr'Z$_2$={screen_fun.Z2}')
    plt.title(title, fontsize='medium')
    plt.show()


def plot_chi_near_one(r0_vals, screen_fun):
    """Plot the chi function for scattering integrals.

    chi is part of the integrand for scattering angle and time integral 
    calculations. The code here repeats the calculation of chi from within
    scatter_integrals for plotting purposes.

    Parameters:
        u (float): Integration variable.
        r0_vals (array-like): Distance of closest approach (RNORM).
        screen_fun (object): Object with `.call()` to calculate the screening function
            for given distance r (RNORM).
    Returns:
        (float): Value of chi function.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 14})

    u1 = np.logspace(-10, 0, 101)
    u = 1 - u1

    for r0 in r0_vals:
        phi0, dphi0 = screen_fun.call(r0)
        phi, _ = screen_fun.call(r0/(1-u**2))
        #chi = (phi0 - phi*(1-u**2)) / u**2
        _, chi = calc_phi_chi(u, r0, screen_fun)
        plt.plot(u1, np.full_like(u, phi0), 'k--')
        plt.plot(u1, chi, label=fr'$R_0$={r0}')

    plt.xscale('log')
    plt.xlabel('1-u')
    plt.ylabel(r'$\chi$(u)')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()
    title = "Potential"
    if isinstance(screen_fun, ZBL_screen):
        title = 'ZBL potential'
    elif isinstance(screen_fun, NLHlin_screen):
        title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                    fr'Z$_2$={screen_fun.Z2}')
    plt.title(title, fontsize='medium')
    plt.show()


def plot_theta_error(screen_fun):
    """Plot the error in the scattering angle.
    
    Parameters:
        screen_fun (object): Screening function instance with `.call()` method
    """
    import matplotlib.pyplot as plt
    from zbl import magic
    plt.rcParams.update({'font.size': 13})

    p_vals = (np.linspace(0.02, 30, 101), 
              np.linspace(0.02, 12, 101), 
              np.linspace(0.02, 4, 101))

    n_absc_vals = (1, 2, 3, 4, 5, 10)

    for ie, e in enumerate((1e-4, 0.1, 100)):
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), layout='constrained')
        theta_ref_vals = []
        theta_magic_vals = []
        setup(n_absc=32)
        for p in p_vals[ie]:
            theta_ref, _ = scatter_integrals(e, p, screen_fun)
            theta_ref_vals.append(theta_ref)
            if type(screen_fun) is ZBL_screen:
                theta_magic = 2 * np.arccos(magic(e, p, screen_fun))
                theta_magic_vals.append(theta_magic)
        theta_ref_vals = np.array(theta_ref_vals)

        for i_absc, n_absc in enumerate(n_absc_vals):
            setup(n_absc)
            errors = []
            for ip, p in enumerate(p_vals[ie]):
                theta, _ = scatter_integrals(e, p, screen_fun)
                theta_ref = theta_ref_vals[ip]
                err = np.abs(theta - theta_ref) / np.abs(theta_ref)
                errors.append(err)
            plt.semilogy(p_vals[ie], errors, f'C{i_absc}')

        if type(screen_fun) is ZBL_screen:
            theta_magic_vals = np.array(theta_magic_vals)
            errors = (np.abs(theta_magic_vals - theta_ref_vals) 
                      / np.abs(theta_ref_vals))
            plt.semilogy(p_vals[ie], errors, 'k--')

        lines = [f'n={n}' for n in n_absc_vals]
        if type(screen_fun) is ZBL_screen:
            lines.append('magic')
        plt.legend(lines, loc=(1.03, 0.2))

        plt.xlim(0, p_vals[ie][-1])
        plt.ylim(1e-7, 1)
        plt.yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
        plt.xlabel('impact parameter P')
        plt.ylabel(r'relative error in scattering angle $\theta$')
        title = "Potential"
        if type(screen_fun) is ZBL_screen:
            title = 'ZBL potential'
        elif type(screen_fun) is NLHlin_screen:
            title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                     fr'Z$_2$={screen_fun.Z2}')
        title += fr', $\varepsilon$={e}'
        plt.title(title, fontsize='medium')
        plt.show()


def plot_tau_error(screen_fun):
    """Plot the error in the time integral.
    
    Parameters:
        screen_fun (object): Screening function instance with `.call()` method
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 13})

    p_vals = (np.linspace(0.02, 30, 101), 
              np.linspace(0.02, 12, 101), 
              np.linspace(0.02, 4, 101))

    n_absc_vals = (1, 2, 3, 4, 5, 10)

    for ie, e in enumerate((1e-4, 0.1, 100)):
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), layout='constrained')
        tau_ref_vals = []
        setup(n_absc=32)
        for p in p_vals[ie]:
            _, tau_ref = scatter_integrals(e, p, screen_fun)
            tau_ref_vals.append(tau_ref)
        tau_ref_vals = np.array(tau_ref_vals)

        for i_absc, n_absc in enumerate(n_absc_vals):
            setup(n_absc)
            errors = []
            for ip, p in enumerate(p_vals[ie]):
                _, tau = scatter_integrals(e, p, screen_fun)
                tau_ref = tau_ref_vals[ip]
                err = np.abs(tau - tau_ref)# / abs(tau_ref)
                errors.append(err)
            plt.semilogy(p_vals[ie], errors, f'C{i_absc}')

        lines = [f'n={n}' for n in n_absc_vals]
        plt.legend(lines, loc=(1.03, 0.25))

        plt.xlim(0, p_vals[ie][-1])
        plt.ylim(1e-7, 1)
        plt.yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
        plt.xlabel('impact parameter P')
        plt.ylabel(r'absolute error in time integral T')
        #plt.ylabel(r'relative error in time integral T')
        title = "Potential"
        if isinstance(screen_fun, ZBL_screen):
            title = 'ZBL potential'
        elif isinstance(screen_fun, NLHlin_screen):
            title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                     fr'Z$_2$={screen_fun.Z2}')
        title += fr', $\varepsilon$={e}'
        plt.title(title, fontsize='medium')
        plt.show()


def plot_tau_over_theta(screen_fun):
    """Plot the time integral as a function of scattering angle.
    
    Parameters:
        screen_fun (object): Screening function instance with `.call()` method
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 13})

    e_vals = np.logspace(-4, 1, 6)
    p_vals = np.linspace(0.02, 30, 201)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4), layout='constrained')

    setup(n_absc=32)
    for ie, e in enumerate(e_vals):
        theta_vals = []
        tau_vals = []
        for p in p_vals:
            theta, tau = scatter_integrals(e, p, screen_fun)
            theta_vals.append(theta)
            tau_vals.append(tau)
        theta_vals = np.array(theta_vals)
        tau_vals = np.array(tau_vals)

        plt.plot(theta_vals*180/np.pi, tau_vals, f'C{ie}', 
                 label=fr'$\varepsilon$={e:.0e}')

    plt.xlim(0, 180)
    #plt.ylim(1e-7, 1)
    #plt.yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    plt.xlabel('scattering angle θ (rad)')
    plt.ylabel(r'time integral τ (RNORM)')
    title = "Potential"
    if isinstance(screen_fun, ZBL_screen):
        title = 'ZBL potential'
    elif isinstance(screen_fun, NLHlin_screen):
        title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                    fr'Z$_2$={screen_fun.Z2}')
    plt.title(title, fontsize='medium')
    plt.legend()
    plt.show()


def plot_sinhalftheta_over_p(screen_fun):
    """Plot the scattering angle as a function of impact parameter.
    
    Parameters:
        screen_fun (object): Screening function instance with `.call()` method
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 13})

    e_vals = np.logspace(-4, 2, 7)
    p_vals = np.linspace(0.01, 30, 201)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4), layout='constrained')

    setup(n_absc=32)
    for ie, e in enumerate(e_vals):
        theta_vals = []
        for p in p_vals:
            theta, _ = scatter_integrals(e, p, screen_fun)
            theta_vals.append(theta)
        theta_vals = np.array(theta_vals)

        plt.plot(p_vals, np.sin(theta_vals/2)**2, f'C{ie}', 
                 label=fr'$\varepsilon$={e:.0e}')
        plt.axvline(2/e, color=f'C{ie}', linestyle='dashed')

    plt.xlim(0, p_vals[-1])
    plt.ylim(0, 1)
    plt.xlabel('impact parameter P')
    plt.ylabel(r'sin$^2$(θ/2)')
    title = "Potential"
    if isinstance(screen_fun, ZBL_screen):
        title = 'ZBL potential'
    elif isinstance(screen_fun, NLHlin_screen):
        title = (fr'NLHlin potential, Z$_1$={screen_fun.Z1}, '
                    fr'Z$_2$={screen_fun.Z2}')
    plt.title(title, fontsize='medium')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    from zbl import ZBL_screen
    from nlhlin import NLHlin_screen, read_coefs

    #screen_fun = ZBL_screen()

    Z1 = 33
    Z2 = 14
    rnorm = 0.4685 / (np.sqrt(np.sqrt(Z1)) + np.sqrt(np.sqrt(Z2)))
    coefs = read_coefs()
    screen_fun = NLHlin_screen(Z1, Z2, rnorm, coefs)
    
    r0_vals = [0.01, 0.1, 1, 10]
    #print(screen_fun.call(1.0))
    #plot_chi(r0_vals, screen_fun)
    #plot_chi_near_zero(r0_vals, screen_fun)
    #plot_chi_near_one(r0_vals, screen_fun)

    #plot_theta_error(screen_fun)
    #plot_tau_error(screen_fun)

    #plot_tau_over_theta(screen_fun)
    plot_sinhalftheta_over_p(screen_fun)
    