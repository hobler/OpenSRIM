"""Calculate the nuclear stopping cross section.
"""
import os
import time
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad
from zbl import ZBL_screen
from krc import KrC_screen
from nlh import NLH_screen
from nlhlin import NLHlin_screen
from cm_scatter import setup, scatter_integrals
from utils import atom, ask_if_save, get_mass, get_density


def calc_sn(e, screen_fun, pmax=None):
    """Calculate nuclear stopping cross section S_n.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
        (float): Maximum impact parameter used in the integration (RNORM).
    
    Returns:
        (float): Nuclear stopping cross section S_n.
        (int): Number of function evaluations in the integration.
    """
    pmax = screen_fun.rmax if pmax is None else pmax

    def func(p, e, screen_fun):
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {theta=}')
        return 2*e*p*np.cos(pi_minus_theta / 2)**2
    
    setup(n_absc=4)
    sn, err, infodict, *rest = quad(func, 0, pmax, 
                                    args=(e, screen_fun),
                                    limit=100, epsabs=0, epsrel=1e-3, 
                                    full_output=True)
    relerr = abs(err / sn)
    neval = infodict["neval"]
    if relerr > 1e-2:
        print(f'Warning: High relative error in S_n calculation:')
        print(f'         {e=}Z1={screen_fun.Z1}, Z2={screen_fun.Z2}, '
              f'{relerr=:.2e}')
        sn = np.nan
    #print(f'{e=}, {sn=}, {err=}, neval={infodict["neval"]}')
    return sn, neval


def calc_qn(e, screen_fun):
    """Calculate nuclear straggling Q_n.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear straggling Q_n.
        (int): Number of function evaluations in the integration.
    """
    def func(p, e, screen_fun):
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {pi_minus_theta=}')
        return 2 * e**2 * p * (np.cos(pi_minus_theta / 2))**4
    setup(n_absc=4)
    qn, err, infodict, *rest = quad(func, 0, screen_fun.rmax, 
                                    args=(e, screen_fun),
                                    epsabs=0, epsrel=1e-3, full_output=True)
    relerr = abs(err / qn)
    neval = infodict["neval"]
    if relerr > 1e-2:
        print(f'Warning: High relative error in Q_n calculation:')
        print(f'         {e=}, Z1={screen_fun.Z1}, Z2={screen_fun.Z2}, '
              f'{relerr=:.2e}')
        qn = np.nan
    #print(f'{e=}, {qn=}, {err=}, neval={infodict["neval"]}')
    return qn, neval


def calc_sn_split(e, screen_fun):
    """Calculate nuclear stopping S_n splitting integral in two parts.

    Parameters:
        e (float): Reduced energy.
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
    
    Returns:
        (float): Nuclear stopping S_n.
        (int): Number of function evaluations in the integration.
    """
    def func(p, e, screen_fun):
        pi_minus_theta, _ = scatter_integrals(e, p, screen_fun)
        #print(f'{p=}, {pi_minus_theta=}')
        return 2*e*p*np.cos(pi_minus_theta / 2)**2
    setup(n_absc=4)
    pmax = 2 / e
    if pmax > screen_fun.rmax:
        sn1, err1, infodict1, *rest = quad(func, 0, screen_fun.rmax, 
                                           args=(e, screen_fun),
                                           epsabs=0, epsrel=1e-2, 
                                           full_output=True)
        sn2, err2, infodict2 = 0, 0, {"neval": 0}
    else:
        sn1, err1, infodict1, *rest = quad(func, 0, pmax, args=(e, screen_fun),
                                           epsabs=0, epsrel=1e-2, 
                                           full_output=True)
        sn2, err2, infodict2, *rest = quad(func, pmax, screen_fun.rmax, 
                                           args=(e, screen_fun),
                                           epsabs=sn1*1e-2, epsrel=1e-2, 
                                           full_output=True)
    sn = sn1 + sn2
    err = err1 + err2
    neval1 = infodict1["neval"]
    neval2 = infodict2["neval"]
    neval = neval1 + neval2
    #print(f'*** {e=}, {sn1=}, {err1=}, {neval1=} ***')
    #print(f'*** {e=}, {sn2=}, {err2=}, {neval2=} ***')
    #print(f'*** {e=}, {sn=}, {err=}, {neval=} ***')
    return sn, neval


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
    sn, err, infodict, *rest = quad(func, 0, screen_fun.rmax**2, 
                                    args=(e, screen_fun), 
                                    limit=100, epsabs=0, epsrel=1e-2, 
                                    full_output=True)
    neval = infodict["neval"]
    #print(f'{e=}, {sn=}, {err=}, neval={infodict["neval"]}')
    return sn, neval


def tab_sn_qn(screen_fun):
    """Tabulate S_n and Q_n for given reduced energies.

    Parameters:
        screen_fun (callable): Function to calculate the screening function
            for given distance r (RNORM).
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
    """
    Z1 = screen_fun.Z1
    Z2 = screen_fun.Z2
    e_vals_10 = np.logspace(-4, 4, 9)
    e_vals = np.sort(np.concatenate((e_vals_10, 0.3*e_vals_10)))
    if type(screen_fun) is ZBL_screen:
        fname = (f'../../data/nuclear_scattering/zbl/sn_qn_tables/'
                 f'sn_qn_table_zbl.txt')
    elif type(screen_fun) is KrC_screen:
        fname = (f'../../data/nuclear_scattering/krc/sn_qn_tables/'
                 f'sn_qn_table_krc.txt')
    elif type(screen_fun) is NLHlin_screen:
        fname = (f'../../data/nuclear_scattering/nlhlin/sn_qn_tables/'
                 f'sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt')
    elif type(screen_fun) is NLH_screen:
        fname = (f'../../data/nuclear_scattering/nlh/sn_qn_tables/'
                 f'sn_qn_table_nlh_{Z1:02d}_{Z2:02d}.txt')
    fname = os.path.join(os.path.dirname(__file__), fname)
    with open(fname, 'w') as f:
        f.write('#  epsilon       S_n          Q_n\n')
        for e in e_vals:
            sn, _ = calc_sn(e, screen_fun)
            qn, _ = calc_qn(e, screen_fun)
            f.write(f'{e:.5e}  {sn:.5e}  {qn:.5e}\n')


def tab_all_sn_qn():
    """Tabulate S_n and Q_n for all Z1-Z2 combinations using NLHlin.
    """
    for Z1 in range(1, 93):
        for Z2 in range(Z1, 93):
            screen_fun = NLHlin_screen(Z1, Z2)
            tab_sn_qn(screen_fun)

def tab_all_sn_qn_nlh():
    """Tabulate S_n and Q_n for all Z1-Z2 combinations using NLH.
    """
    for Z1 in range(1, 93):
        for Z2 in range(Z1, 93):
            print(f"Calculating S_n and Q_n for Z1={Z1}, Z2={Z2} using NLH...")
            screen_fun = NLH_screen(Z1, Z2)
            tab_sn_qn(screen_fun)


def plot_sn(p1, p2):
    """Plot reducedS_n for homonuclear pairs.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})

    Z1_vals = range(1, 93)
    sn_vals_e0p01 = []

    fig = plt.figure()

    for Z1 in Z1_vals:
            Z2 = Z1
            fname = os.path.join(
                os.path.dirname(__file__), 
                '../../data/nuclear_scattering/nlhlin/sn_qn_tables', 
                f'sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt'
            )
            data = np.loadtxt(fname, comments='#')
            e_vals = data[:, 0]
            sn_vals = data[:, 1]
            # convert from NLHlin to desired scaling
            fac = (Z1**p1 + Z2**p1)**p2 / (Z1**0.5 + Z2**0.5)**0.5
            e_vals /= fac
            sn_vals *= fac
            # mask values below 10 eV
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
            mask = e_vals *2*14.4*Z1*Z2/rnorm > 10
            plt.loglog(e_vals[mask], sn_vals[mask], color=cmap((Z1-1)/92), 
                       zorder=Z1)
            plt.loglog(e_vals, sn_vals, ':', color=cmap((Z1-1)/92), zorder=-Z1)
            
            sn_val_e0p01 = np.exp(
                np.interp(np.log(0.01), np.log(e_vals), np.log(sn_vals))
            )
            sn_vals_e0p01.append(sn_val_e0p01)

    p1 = round(p1, 2)
    p2 = round(p2, 2)
    
    if (p1, p2) == (0.23, 1):
        fname = os.path.join(
            os.path.dirname(__file__), 
            '../../data/nuclear_scattering/zbl/sn_qn_tables',
            'sn_qn_table_zbl.txt'
        )
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        sn_vals = data[:, 1]
        plt.loglog(e_vals, sn_vals, 'k--', label='ZBL', zorder=100)
        plt.legend(loc='upper right')
        idx = np.argmin(np.abs(e_vals - 0.01))
        sn_val_zbl_e0p01 = sn_vals[idx]

    if (p1, p2) == (0.5, 0.67):
        fname = os.path.join(
            os.path.dirname(__file__), 
            '../../data/nuclear_scattering/krc/sn_qn_tables',
            'sn_qn_table_krc.txt'
        )
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
        text += fr'(Z$_1^{{{p1}}}$+Z$_2^{{{p1}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')

    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                 ax=plt.gca(),
                 label=r'atomic number Z$_1$=Z$_2$', 
                 ticks=ticks)
    plt.xlim(1e-4, 1e3)
    plt.ylim(1e-4, 1)
    plt.xlabel('Reduced energy ε')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$')
    plt.tight_layout()
    plt.show()

    fname = os.path.join(
        os.path.dirname(__file__), 
        f'figs/sn_p{p1}_p{p2}_homo.pdf'
    )
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")

    fig = plt.figure()

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
        text += fr'(Z$_1^{{{p1}}}$+Z$_2^{{{p1}}}$)'
    if p2 != 1:
        text += fr'$^{{{p2}}}$'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')
    plt.xlim(0, 93)
    plt.ylim(0, None)
    plt.xlabel('Atomic number Z$_1$=Z$_2$')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$ at ε=0.01')
    plt.tight_layout()
    plt.show()

    fname = os.path.join(
        os.path.dirname(__file__), 
        f'figs/sn_p{p1}_p{p2}_homo_e0.01.pdf'
    )
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


def plot_sn_at_e(e, p1, p2):
    """Plot S_n for all atom pairs.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})
    fig = plt.figure()

    p1 = round(p1, 2)
    p2 = round(p2, 2)
    
    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 ax=plt.gca(),
                 label=r'atomic number Z$_2$', 
                 ticks=ticks)

    Z1_vals = range(1, 93)

    if (p1, p2) == (0.23, 1):
        fname = os.path.join(
            os.path.dirname(__file__), 
            '../../data/nuclear_scattering/zbl/sn_qn_tables',
            'sn_qn_table_zbl.txt'
        )
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        sn_vals = data[:, 1]
        idx = np.argmin(np.abs(e_vals - e))
        sn_zbl = sn_vals[idx]

    if (p1, p2) == (0.5, 0.67):
        fname = os.path.join(
            os.path.dirname(__file__), 
            '../../data/nuclear_scattering/krc/sn_qn_tables',
            'sn_qn_table_krc.txt'
        )
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
            fname = os.path.join(
                os.path.dirname(__file__), 
                '../../data/nuclear_scattering/nlhlin/sn_qn_tables',
                f'sn_qn_table_nlhlin_{Z1_:02d}_{Z2_:02d}.txt'
            )
            data = np.loadtxt(fname, comments='#')
            e_vals = data[:, 0]
            sn_vals = data[:, 1]
            # convert from NLHlin to desired scaling
            fac = (Z1**p1 + Z2**p1)**p2 / (Z1**0.5 + Z2**0.5)**0.5
            e_vals /= fac
            sn_vals *= fac
            sn = np.exp(np.interp(np.log(e), np.log(e_vals), np.log(sn_vals)))
            plt.plot(Z1, sn, '.', color=cmap((Z2-1)/92), zorder=Z2)

    if (p1, p2) == (0.23, 1):
        potname = "zbl"
        plt.axhline(sn_zbl, color='k', linestyle='--', 
                    label='ZBL', zorder=100)
        plt.legend(loc='center right')
    if (p1, p2) == (0.5, 0.67):
        potname = "krc"
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
    #plt.text(0.95, 0.05, text, 
    #         horizontalalignment='right', verticalalignment='bottom',
    #         transform=plt.gca().transAxes, fontsize='medium')
    plt.xlim(0, 93)
    plt.ylim(0, None)
    plt.xlabel('Atomic number Z$_1$')
    plt.ylabel(r'Reduced nuclear stopping s$_\mathrm{n}$ at' f' ε={e}')
    plt.tight_layout()
    plt.show()

    fname = os.path.join(
        os.path.dirname(__file__), 
        f'figs/sn_{potname}_all_e0.01.pdf'
    )
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


def plot_sn_cutoff(energies, Z1_list, Z2_list):
    """Calculate S_n for ZBL with and without cutoff impact parameter.
    
    The cutoff impact parameter is calculated from the target density.

    Parameters:
        energies (list): List of reduced energies to calculate S_n.
        Z1_list (int or list of int): Atomic numbers of projectiles.
        Z2_list (int or list of int): Atomic numbers of targets.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 14})

    # convert Z1 and Z2 to lists if they are not already
    if not isinstance(Z1_list, list):
            Z1_list = [Z1_list]
    if not isinstance(Z2_list, list):
            Z2_list = [Z2_list]

    lines = []
    for i, (Z1, Z2) in enumerate(zip(Z1_list, Z2_list)):
        # screening function
        screen_fun = ZBL_screen(Z1, Z2)
        print('rnorm=', screen_fun.rnorm)

        # Convert energies and pmax to reduced units
        M1 = get_mass(Z1)
        M2 = get_mass(Z2)
        enorm = (M1 + M2) / M2 * 14.4 * Z1 * Z2 / screen_fun.rnorm
        eps_vals = np.array(energies) / enorm
        print('eps=', eps_vals)

        if Z2 == 6:
            density = 0.176
        else:
            density = get_density(Z2)
        pmax_A = np.pi**(-1/2) * density**(-1/3)  # in Angstroms
        pmax = pmax_A / screen_fun.rnorm
        print('pmax= ', pmax)

        # Calculate S_n without cutoff
        sn_vals = []
        for eps in eps_vals:
            sn, _ = calc_sn(eps, screen_fun)
            sn_vals.append(sn)

        # Calculate S_n with cutoff
        sn_cutoff_vals = []
        for eps in eps_vals:
            sn_cutoff, _ = calc_sn(eps, screen_fun, pmax=pmax)
            sn_cutoff_vals.append(sn_cutoff)

        # Convert S_n back to physical units
        unscale = np.pi * 4*M1*M2/(M1+M2)**2 * enorm * screen_fun.rnorm**2
        sn_vals = np.array(sn_vals)
        sn_vals *= unscale
        sn_cutoff_vals = np.array(sn_cutoff_vals)
        sn_cutoff_vals *= unscale

        line, = plt.loglog(energies, sn_vals, f'C{i}', 
                   label=f'{atom[Z1]} in {atom[Z2]}')
        lines.append(line)
        plt.loglog(energies, sn_cutoff_vals, f'C{i}--')

    first_legend = plt.gca().legend(handles=lines, loc='upper left')

    lines = []
    line, = plt.plot(energies[0], sn_vals[0], 'k-', 
                     label=r'p$_\mathrm{max} \rightarrow \infty$', zorder=0)
    lines.append(line)
    line, = plt.plot(energies[0], sn_cutoff_vals[0], 'k--', 
                     label=r'p$_\mathrm{max}=\pi^{-1/2}N^{-1/3}$', zorder=0)
    lines.append(line)
    plt.legend(handles=lines, loc='lower right')

    plt.gca().add_artist(first_legend)

    plt.ylim(10, 1e4)
    plt.xlabel('Energy (eV)')
    plt.ylabel(r'Nuclear stopping cross section (eVÅ$^2$)')
    plt.title(f'{screen_fun.name} potential', fontsize='medium')
    plt.tight_layout()
    plt.show()


def plot_qn(p1, p2):
    """Plot Q_n for homonuclear pairs.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})

    p1 = round(p1, 2)
    p2 = round(p2, 2)
    
    fig = plt.figure()

    Z1_vals = range(1, 93)
    qn_vals_e0p01 = []

    for Z1 in Z1_vals:
            Z2 = Z1
            fname = os.path.join(
                os.path.dirname(__file__), 
                '../../data/nuclear_scattering/nlhlin/sn_qn_tables',
                f'sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt'
            )
            data = np.loadtxt(fname, comments='#')
            e_vals = data[:, 0]
            qn_vals = data[:, 2]
            # convert from NLHlin to desired scaling
            fac = (Z1**p1 + Z2**p1)**p2 / (Z1**0.5 + Z2**0.5)**0.5
            e_vals /= fac
            # mask values below 10 eV
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2  # small error due to rounding
            mask = e_vals *2*14.4*Z1*Z2/rnorm > 10
            plt.loglog(e_vals[mask], qn_vals[mask], color=cmap((Z1-1)/92), 
                       zorder=Z1)
            plt.loglog(e_vals, qn_vals, ':',color=cmap((Z1-1)/92), zorder=Z1)
            qn_val_e0p01 = np.exp(
                np.interp(np.log(0.01), np.log(e_vals), np.log(qn_vals))
            )
            qn_vals_e0p01.append(qn_val_e0p01)

    if (p1, p2) == (0.23, 1):
        fname = os.path.join(
            os.path.dirname(__file__), 
            '../../data/nuclear_scattering/zbl/sn_qn_tables',
            'sn_qn_table_zbl.txt'
        )
        data = np.loadtxt(fname, comments='#')
        e_vals = data[:, 0]
        qn_vals = data[:, 2]
        plt.plot(e_vals, qn_vals, 'k--', label='ZBL', zorder=100)
        plt.legend(loc='upper left')
        idx = np.argmin(np.abs(e_vals - 0.01))
        qn_val_zbl_e0p01 = qn_vals[idx]

    if (p1, p2) == (0.5, 0.67):
        fname = os.path.join(
            os.path.dirname(__file__), 
            '../../data/nuclear_scattering/krc/sn_qn_tables',
            'sn_qn_table_krc.txt'
        )
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
                 ax=plt.gca(), 
                label=r'atomic number Z$_1$=Z$_2$', 
                ticks=ticks)
    plt.xlim(1e-4, 1e3)
    plt.ylim(1e-6, 1)
    plt.xlabel('Reduced energy ε')
    plt.ylabel(r'Reduced nuclear straggling q$_\mathrm{n}$')
    plt.tight_layout()
    plt.show()

    fname = os.path.join(
        os.path.dirname(__file__), 
        f'figs/qn_p{p1}_p{p2}_homo.pdf'
    )
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")

    fig = plt.figure()

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

    fname = os.path.join(
        os.path.dirname(__file__), 
        f'figs/qn_p{p1}_p{p2}_homo_e0.01.pdf'
    )
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


def test_zbl():
    """Test ZBL stopping against tabulated values.
    """
    import matplotlib.pyplot as plt

    fname = os.path.join(
        os.path.dirname(__file__), 
        '../../data/nuclear_scattering/zbl/sn_qn_tables',
        'sn_qn_table_zbl.txt'
    )
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


def compare_calc_sn_qn(screen_fun):
    """Compare values and function evaluations for Sn calculation."""""

    e_vals = np.logspace(-4, 4, 9)
    
    time_sn = 0.0
    time_sn_p2 = 0.0
    time_sn_split = 0.0
    time_qn = 0.0

    for e in e_vals:
        start_time = time.time()
        sn, neval = calc_sn(e, screen_fun)
        end_time = time.time()
        time_sn += end_time - start_time
        print(f'calc_sn      : {e=:5.0f}, {sn=:.6f}, neval={neval}')
        
        start_time = time.time()
        sn, neval = calc_sn_p2(e, screen_fun)
        end_time = time.time()
        time_sn_p2 += end_time - start_time
        print(f'calc_sn_p2   : {e=:5.0f}, {sn=:.6f}, neval={neval}')

        start_time = time.time()
        sn, neval = calc_sn_split(e, screen_fun)
        end_time = time.time()
        time_sn_split += end_time - start_time
        print(f'calc_sn_split: {e=:5.0f}, {sn=:.6f}, neval={neval}')

        start_time = time.time()
        qn, neval = calc_qn(e, screen_fun)
        end_time = time.time()
        time_qn += end_time - start_time
        print(f'calc_qn      : {e=:5.0f}, {qn=:.6f}, neval={neval}')

    print(f'Total time for calc_sn: {time_sn:.3f} s')
    print(f'Total time for calc_sn_p2: {time_sn_p2:.3f} s')
    print(f'Total time for calc_sn_split: {time_sn_split:.3f} s')
    print(f'Total time for calc_qn: {time_qn:.3f} s')


if __name__ == "__main__":

    #test_zbl()
    
    #screen_fun = ZBL_screen()
    #screen_fun = KrC_screen()

    #Z1 = 4
    #Z2 = 70
    #aZBL = 0.4685 / (Z1**0.23 + Z2**0.23)**1
    #aNLHlin = 0.4685 / (Z1**0.5 + Z2**0.5)**0.5
    #a = aZBL
    #screen_fun = NLHlin_screen(Z1, Z2)#, a)
    #print(calc_sn(0.3*aNLHlin/a, screen_fun))
    
    #print(calc_sn(1e-5, screen_fun))

    #compare_calc_sn_qn(screen_fun)

    #start_time = time.time()
    #tab_sn_qn(screen_fun)
    #tab_all_sn_qn()
    #tab_all_sn_qn_nlh()
    #print(f'Total time: {time.time() - start_time:.3f} s')

    #plot_sn(p1=0.5, p2=0.5)
    #plot_sn(p1=0.23, p2=1)
    #plot_sn(p1=0.5, p2=2/3)
    #plot_qn(p1=0.5, p2=0.5)
    #plot_qn(p1=0.23, p2=1)
    #plot_qn(p1=0.5, p2=2/3)
    #plot_sn_at_e(0.01, p1=0.5, p2=0.5)
    #plot_sn_at_e(0.01, p1=0.23, p2=1)
    #plot_sn_at_e(0.01, p1=0.5, p2=2/3)

    energies = np.logspace(1, 4, 7)
    Z1 = [7, 79, 74]
    Z2 = [6, 14, 74]
    plot_sn_cutoff(energies, Z1, Z2)