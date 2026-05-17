"""Calculate the nuclear stopping cross section.
"""
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
from cm_scatter import setup, scatter_integrals

def calc_sn(e, screen_fun):
    """Calculate nuclear stopping S_n.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear stopping S_n.
    """
    def func(p, e, screen_fun):
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {theta=}')
        return 2*e*p*np.cos(pi_minus_theta / 2)**2
    setup(n_absc=4)
    sn, err, *infodict = quad(func, 0, screen_fun.rmax, args=(e, screen_fun),
                             limit=100, epsabs=0, epsrel=1e-3, full_output=True)
    relerr = abs(err / sn)
    if relerr > 1e-2:
        print(f'Warning: High relative error in S_n calculation:')
        print(f'         {e=}Z1={screen_fun.Z1}, Z2={screen_fun.Z2}, '
              f'{relerr=:.2e}')
        sn = np.nan
    #print(f'{e=}, {sn=}, {err=}, neval={infodict[0]["neval"]}')
    return sn


def calc_qn(e, screen_fun):
    """Calculate nuclear straggling Q_n.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear stopping S_n.
    """
    def func(p, e, screen_fun):
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {pi_minus_theta=}')
        return 2 * e**2 * p * (np.cos(pi_minus_theta / 2))**4
    setup(n_absc=4)
    qn, err, *infodict = quad(func, 0, screen_fun.rmax, args=(e, screen_fun),
                              epsabs=0, epsrel=1e-3, full_output=True)
    relerr = abs(err / qn)
    if relerr > 1e-2:
        print(f'Warning: High relative error in Q_n calculation:')
        print(f'         {e=}, Z1={screen_fun.Z1}, Z2={screen_fun.Z2}, '
              f'{relerr=:.2e}')
        qn = np.nan
    #print(f'{e=}, {qn=}, {err=}, neval={infodict[0]["neval"]}')
    return qn


def calc_sn_split(e, screen_fun):
    """Calculate nuclear stopping S_n splitting integral in two parts.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear stopping S_n.
    """
    def func(p, e, screen_fun):
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {pi_minus_theta=}')
        return 2*e*p*np.cos(pi_minus_theta / 2)**2
    setup(n_absc=4)
    pmax = 2 / e
    if pmax > screen_fun.rmax:
        sn1, err1, infodict1 = quad(func, 0, screen_fun.rmax, args=(e, screen_fun),
                                    epsabs=0, epsrel=1e-2, full_output=True)
        sn2, err2, infodict2 = 0, 0, {"neval": 0}
    else:
        sn1, err1, infodict1 = quad(func, 0, pmax, args=(e, screen_fun),
                                    epsabs=0, epsrel=1e-2, full_output=True)
        sn2, err2, infodict2, *rest = quad(func, pmax, screen_fun.rmax, 
                                           args=(e, screen_fun),
                                           epsabs=sn1*1e-2, epsrel=1e-2, 
                                           full_output=True)
    sn = sn1 + sn2
    err = err1 + err2
    neval1 = infodict1["neval"]
    neval2 = infodict2["neval"]
    neval = neval1 + neval2
    print(f'*** {e=}, {sn1=}, {err1=}, {neval1=} ***')
    print(f'*** {e=}, {sn2=}, {err2=}, {neval2=} ***')
    print(f'*** {e=}, {sn=}, {err=}, {neval=} ***')
    return sn


def calc_sn_p2(e, screen_fun):
    """Calculate nuclear stopping S_n.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear stopping S_n.
    """
    def func(p2, e, screen_fun):
        p = np.sqrt(p2)
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {pi_minus_theta=}')
        return e*np.cos(pi_minus_theta / 2)**2
    setup(n_absc=4)
    sn, err, infodict = quad(func, 0, screen_fun.rmax**2, args=(e, screen_fun), 
                             limit=100, epsabs=0, epsrel=1e-2, full_output=True)
    print(f'{e=}, {sn=}, {err=}, neval={infodict["neval"]}')
    return sn


def calc_sn_laguerre(e, screen_fun):
    """Calculate nuclear stopping S_n using Gauss-Laguerre integration.

    DOES NOT WORK!

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear stopping S_n.
    """
    def func(p, e, screen_fun):
        pi_minus_theta = np.array([scatter_integrals(e, p_val, screen_fun)[0] 
                          for p_val in p])
        #print(f'{p=}, {pi_minus_theta=}')
        return 2*e*p*np.cos(pi_minus_theta / 2)**2
    setup(n_absc=4)
    n_iter = 100
    u, w = roots_legendre(n_iter)
    sn = np.dot(w, func(u, e, screen_fun)*np.exp(u))
    print(f'{e=}, {sn=}, neval={n_iter}')
    return sn

def tab_sn_qn(e_vals, screen_fun, p1, p2):
    """Tabulate S_n and Q_n for given reduced energies.

    Parameters:
        e_vals (array-like): Reduced energy values.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
    """
    Z1 = screen_fun.Z1
    Z2 = screen_fun.Z2
    p1 = round(p1, 2)
    p2 = round(p2, 2)
    if type(screen_fun) is ZBL_screen:
        fname = f'../data/p{p1}_{p2}/sn_qn_table_zbl.txt'
    elif type(screen_fun) is KrC_screen:
        fname = f'../data/p{p1}_{p2}/sn_qn_table_krc.txt'
    elif type(screen_fun) is NLHlin_screen:
        fname = f'../data/p{p1}_{p2}/sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt'
    with open(fname, 'w') as f:
        f.write('#  epsilon       S_n          Q_n\n')
        for e in e_vals:
            sn = calc_sn(e, screen_fun)
            qn = calc_qn(e, screen_fun)
            f.write(f'{e:.5e}  {sn:.5e}  {qn:.5e}\n')

def tab_all_sn_qn(p1, p2):
    """Tabulate S_n and Q_n for all Z1-Z2 combinations using NLHlin.
    """
    e_vals_10 = np.logspace(-4, 4, 9)
    e_vals = np.sort(np.concatenate((e_vals_10, 0.3*e_vals_10)))

    for Z1 in range(1, 93):
        for Z2 in range(Z1, 93):
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
            screen_fun = NLHlin_screen(Z1, Z2, rnorm)
            tab_sn_qn(e_vals, screen_fun, p1, p2)

def plot_sn(p1, p2):
    """Plot S_n for homonuclear pairs.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = mpl.cm.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})

    p1 = round(p1, 2)
    p2 = round(p2, 2)
    
    Z1_vals = range(1, 93)
    sn_vals_e0p01 = []

    for Z1 in Z1_vals:
            Z2 = Z1
            fname = f'../data/p{p1}_{p2}/sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt'
            data = np.loadtxt(fname, comments='#')
            e_vals = data[:, 0]
            sn_vals = data[:, 1]
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
            mask = e_vals *2*14.4*Z1*Z2/rnorm > 10
            plt.loglog(e_vals[mask], sn_vals[mask], color=cmap((Z1-1)/92), 
                       zorder=Z1)
            plt.loglog(e_vals, sn_vals, ':', color=cmap((Z1-1)/92), zorder=-Z1)
            idx = np.argmin(np.abs(e_vals - 0.01))
            sn_vals_e0p01.append(sn_vals[idx])

    if (p1, p2) == (0.23, 1):
        fname = f'../data/p{p1}_{p2}/sn_qn_table_zbl.txt'
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        sn_vals = data[:, 1]
        plt.loglog(e_vals, sn_vals, 'k--', label='ZBL', zorder=100)
        plt.legend(loc='upper right')
        idx = np.argmin(np.abs(e_vals - 0.01))
        sn_val_zbl_e0p01 = sn_vals[idx]

    if (p1, p2) == (0.5, 0.67):
        fname = f'../data/p{p1}_{p2}/sn_qn_table_krc.txt'
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        sn_vals = data[:, 1]
        plt.loglog(e_vals, sn_vals, 'k--', label='KrC', zorder=100)
        plt.legend(loc='upper right')
        idx = np.argmin(np.abs(e_vals - 0.01))
        sn_val_krc_e0p01 = sn_vals[idx]

    text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
    if p1 == 1:
        text += fr'(Z$_1$+Z$_2$)'
    else:
        text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2:.2f}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')

    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                label=r'atomic number Z$_1$=Z$_2$', 
                ticks=ticks)
    plt.xlim(1e-4, 1e3)
    plt.ylim(1e-4, 1)
    plt.xlabel('Reduced energy ε')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$')
    plt.tight_layout()
    plt.show()

    plt.plot(Z1_vals, sn_vals_e0p01, 'o-')
    if (p1, p2) == (0.23, 1):
        plt.axhline(sn_val_zbl_e0p01, color='k', linestyle='--', 
                    label='ZBL', zorder=100)
        plt.legend(loc='center right')
    if (p1, p2) == (0.5, 0.67):
        plt.axhline(sn_val_krc_e0p01, color='k', linestyle='--', 
                    label='KrC', zorder=100)
        plt.legend(loc='center right')
    text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
    if p1 == 1:
        text += fr'(Z$_1$+Z$_2$)'
    else:
        text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2:.2f}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')
    plt.xlim(0, 93)
    plt.ylim(0, None)
    plt.xlabel('Atomic number Z$_1$=Z$_2$')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$ at ε=0.01')
    plt.tight_layout()
    plt.show()


def plot_sn_at_e(e, p1, p2):
    """Plot S_n for all atom pairs.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = mpl.cm.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})

    p1 = round(p1, 2)
    p2 = round(p2, 2)
    
    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                label=r'atomic number Z$_2$', 
                ticks=ticks)

    Z1_vals = range(1, 93)
    sn_vals_e0p01 = []

    if (p1, p2) == (0.23, 1):
        fname = f'../data/p{p1}_{p2}/sn_qn_table_zbl.txt'
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        sn_vals = data[:, 1]
        idx = np.argmin(np.abs(e_vals - e))
        sn_zbl = sn_vals[idx]

    if (p1, p2) == (0.5, 0.67):
        fname = f'../data/p{p1}_{p2}/sn_qn_table_krc.txt'
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        sn_vals = data[:, 1]
        idx = np.argmin(np.abs(e_vals - e))
        sn_krc = sn_vals[idx]


    for Z1 in Z1_vals:
#        for Z2 in range(Z1, 93):
        for Z2 in Z1_vals:
            Z1_ = min(Z1, Z2)
            Z2_ = max(Z1, Z2)
            fname = f'../data/p{p1}_{p2}/sn_qn_table_nlhlin_{Z1_:02d}_{Z2_:02d}.txt'
            data = np.loadtxt(fname, comments='#')
            e_vals = data[:, 0]
            sn_vals = data[:, 1]
            idx = np.argmin(np.abs(e_vals - e))
            sn = sn_vals[idx]
            plt.plot(Z1, sn, '.', color=cmap((Z2-1)/92), zorder=Z2)

    if (p1, p2) == (0.23, 1):
        plt.axhline(sn_zbl, color='k', linestyle='--', 
                    label='ZBL', zorder=100)
        plt.legend(loc='center right')
    if (p1, p2) == (0.5, 0.67):
        plt.axhline(sn_krc, color='k', linestyle='--', 
                    label='KrC', zorder=100)
        plt.legend(loc='center right')
    text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
    if p1 == 1:
        text += fr'(Z$_1$+Z$_2$)'
    else:
        text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2:.2f}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')
    plt.xlim(0, 93)
    plt.ylim(0, None)
    plt.xlabel('Atomic number Z$_1$')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$ at' f' ε={e}')
    plt.tight_layout()
    plt.show()


def plot_qn(p1, p2):
    """Plot Q_n for homonuclear pairs.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = mpl.cm.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})

    p1 = round(p1, 2)
    p2 = round(p2, 2)
    
    Z1_vals = range(1, 93)
    qn_vals_e0p01 = []

    for Z1 in Z1_vals:
            Z2 = Z1
            fname = f'../data/p{p1}_{p2}/sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt'
            data = np.loadtxt(fname, comments='#')
            e_vals = data[:, 0]
            qn_vals = data[:, 2]
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
            mask = e_vals *2*14.4*Z1*Z2/rnorm > 10
            plt.loglog(e_vals[mask], qn_vals[mask], color=cmap((Z1-1)/92), 
                       zorder=Z1)
            plt.loglog(e_vals, qn_vals, ':',color=cmap((Z1-1)/92), zorder=Z1)
            idx = np.argmin(np.abs(e_vals - 0.01))
            qn_vals_e0p01.append(qn_vals[idx])

    if (p1, p2) == (0.23, 1):
        fname = f'../data/p{p1}_{p2}/sn_qn_table_zbl.txt'
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        qn_vals = data[:, 2]
        plt.plot(e_vals, qn_vals, 'k--', label='ZBL', zorder=100)
        plt.legend(loc='upper left')
        idx = np.argmin(np.abs(e_vals - 0.01))
        qn_val_zbl_e0p01 = qn_vals[idx]

    if (p1, p2) == (0.5, 0.67):
        fname = f'../data/p{p1}_{p2}/sn_qn_table_krc.txt'
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        qn_vals = data[:, 2]
        plt.plot(e_vals, qn_vals, 'k--', label='KrC', zorder=100)
        plt.legend(loc='upper left')
        idx = np.argmin(np.abs(e_vals - 0.01))
        qn_val_krc_e0p01 = qn_vals[idx]

    text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
    if p1 == 1:
        text += fr'(Z$_1$+Z$_2$)'
    else:
        text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2:.2f}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')

    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                label=r'atomic number Z$_1$=Z$_2$', 
                ticks=ticks)
    plt.xlim(1e-4, 1e3)
    plt.ylim(1e-6, 1)
    plt.xlabel('Reduced energy ε')
    plt.ylabel(r'Reduced nuclear straggling q$_\mathrm{n}$')
    plt.tight_layout()
    plt.show()

    plt.plot(Z1_vals, qn_vals_e0p01, 'o-')
    if (p1, p2) == (0.23, 1):
        plt.axhline(qn_val_zbl_e0p01, color='k', linestyle='--', 
                    label='ZBL', zorder=100)
        plt.legend(loc='center right')
    if (p1, p2) == (0.5, 0.67):
        plt.axhline(qn_val_krc_e0p01, color='k', linestyle='--', 
                    label='KrC', zorder=100)
        plt.legend(loc='center right')
    text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
    if p1 == 1:
        text += fr'(Z$_1$+Z$_2$)'
    else:
        text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2:.2f}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')
    plt.xlim(0, 93)
    plt.ylim(0, None)
    plt.xlabel('Atomic number Z$_1$=Z$_2$')
    plt.ylabel(r'Reduced nuclear straggling q$_\mathrm{n}$ at ε=0.01')
    plt.tight_layout()
    plt.show()


def test_zbl():
    """Test ZBL stopping against tabulated values.
    """
    import matplotlib.pyplot as plt

    fname = '../data/p0.23_1/sn_qn_table_zbl.txt'
    data = np.loadtxt(fname, comments='#')
    energies = data[:, 0]
    sn_vals = data[:, 1]
    qn_vals = data[:, 2]

    def sn_fit_func(e, a, b, c, d):
         return np.where(e < 30,
                         0.5*np.log(1 + a*e) / (e + b*e**c + d*e**0.5),
                         0.5*np.log(e) / e)

    sn_vals_zbl = sn_fit_func(energies, 1.1383, 0.01321, 0.21226, 0.19593)

    plt.loglog(energies, sn_vals, 'b.', label='Tabulated ZBL')
    plt.loglog(energies, sn_vals_zbl, 'r-', label='ZBL fit')
    plt.xlabel('Reduced energy ε')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$')
    plt.legend()
    plt.show()

    def qn_fit_func(e, a, b, c, d):
         return 1 / (4 + a*e**b + c*e**d)
    
    qn_vals_zbl = qn_fit_func(energies, 0.19676, -1.6991, 6.5841, -1.0494)

    plt.loglog(energies, qn_vals, 'b.', label='Tabulated ZBL')
    plt.loglog(energies, qn_vals_zbl, 'r-', label='ZBL fit')
    plt.xlabel('Reduced energy ε')
    plt.ylabel(r'Reduced nuclear straggling q$_\mathrm{n}$')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    from zbl import ZBL_screen
    from krc import KrC_screen
    from nlhlin import NLHlin_screen
    import time

    #test_zbl()
    
    #screen_fun = ZBL_screen()
    screen_fun = KrC_screen()

    Z1 = 7
    Z2 = 7
    #p1, p2 = 0.23, 1    # ZBL
    #p1, p2 = 1/2, 1/2   # NLHlin
    p1, p2 = 1/2, 2/3   # Firsov
    rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
    #screen_fun = NLHlin_screen(Z1, Z2, rnorm)
    
    sum_time = 0.0
    e_vals_10 = np.logspace(-4, 4, 9)
    e_vals = np.sort(np.concatenate((e_vals_10, 0.3*e_vals_10)))
    #e_vals = np.logspace(-4, 4, 9)
    #for e in e_vals:
        #start_time = time.time()
        #sn = calc_sn(e, screen_fun)         # better convergence
        #sn = calc_sn_p2(e, screen_fun)
        #sn = calc_sn_laguerre(e, screen_fun)
        #sn = calc_sn_split(e, screen_fun)         # better convergence
        #print(f'{e=}, {sn=}, time taken: {time.time() - start_time:.3f} s\n')
        #qn = calc_qn(e, screen_fun)
        #print(f'{e=}, {qn=}, time taken: {time.time() - start_time:.3f} s\n')
        #sum_time += time.time() - start_time
    #print(f'Total time: {sum_time:.3f} s')

    #start_time = time.time()
    tab_sn_qn(e_vals, screen_fun, p1, p2)
    #tab_all_sn_qn(p1, p2)
    #print(f'Total time: {time.time() - start_time:.3f} s')

    #plot_sn(p1, p2)
    #plot_qn(p1, p2)
    #plot_sn_at_e(0.01, p1, p2)
