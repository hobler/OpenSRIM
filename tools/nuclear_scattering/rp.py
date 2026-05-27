import os
import numpy as np
from matplotlib.ticker import MultipleLocator
from table1d import Table1D
from sn_fit import sn_fit_func
from scipy.integrate import solve_ivp
from utils import atom, ask_if_save


class Se:
    """Defines the electronic stopping cross section.
    """
    def __init__(self, Z1, Z2):
        """Setup electronic stopping cross section for given atomic numbers.

        Parameters:
            Z1: atomic number of atom 1
            Z2: atomic number of atom 2
        """
        self.Z1 = Z1
        self.Z2 = Z2

        fname = os.path.join(os.path.dirname(__file__), 
                             f'../../data/SRIM_setab/SRIM2013-{Z1:02d}.dat')
        self.data = np.loadtxt(fname, comments='#')
        e_vals = self.data[:, 0]
        se_vals = self.data[:, Z2]

        self.se_table = Table1D(e_vals, se_vals)
    
    def __call__(self, e):
        """Calculate the electronic stopping cross section.

        Parameters:
            e (float): energy of projectile (eV)

        Returns:
            (float): electronic stopping cross section (eV A^2)
        """
        if e < self.se_table.x[0]:
            return self.se_table.y[0]
        elif e > self.se_table.x[-1]:
            return self.se_table.y[-1]
        else:
            return self.se_table.interpolate(e)

class Sn:
    """Defines the nuclear stopping cross section.
    """
    def __init__(self, Z1, Z2, M1, M2, zbl=False, krc=False):
        """Setup nuclear stopping cross section for given atomic numbers.

        Parameters:
            Z1: atomic number of atom 1
            Z2: atomic number of atom 2
            M1: mass of atom 1 (amu)
            M2: mass of atom 2 (amu)
            zbl: whether to use ZBL screening function for nuclear stopping.
            krc: whether to use the KRC screening function for nuclear stopping.
        """
        self.Z1 = min(Z1, Z2)
        self.Z2 = max(Z1, Z2)
        self.M1 = M1
        self.M2 = M2

        directory = os.path.join(os.path.dirname(__file__), 
                                 '../../data/nuclear_scattering')
        if zbl:
            fname = os.path.join(directory, "zbl/sn_fit_params_zbl.txt")
            a, b, c, d, *_ = np.loadtxt(fname, unpack=True)
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            p1 = 0.23
            p2 = 1
        elif krc:
            fname = os.path.join(directory, f"krc/sn_fit_params_krc.txt")
            a, b, c, d, *_ = np.loadtxt(fname, unpack=True)
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            p1 = 0.5
            p2 = 2/3
        else:
            fname = os.path.join(directory, f"nlhlin/sn_fit_params_nlhlin.txt")
            Z1_, Z2_, a, b, c, d, *_ = np.loadtxt(fname, unpack=True)
            idx = np.where((Z1_ == self.Z1) & (Z2_ == self.Z2))[0]
            self.a = a[idx]
            self.b = b[idx]
            self.c = c[idx]
            self.d = d[idx]
            p1 = 0.5
            p2 = 0.5
        print(f'{Z1=}, {Z2=}:', self.a, self.b, self.c, self.d)
        self.rnorm = 0.4685 / (self.Z1**p1 + self.Z2**p1)**p2
        self.enorm = (M1+M2)/M2 * 14.4 * self.Z1 * self.Z2 / self.rnorm
    
    def __call__(self, e):
        """Calculate the nuclear stopping cross section.

        Parameters:
            e (float): energy of projectile (eV)

        Returns:
            (float): nuclear stopping cross section (eV A^2)
        """
        e_norm = e / self.enorm
        sn_norm = sn_fit_func(e_norm, self.a, self.b, self.c, self.d)
        sn = (4 * np.pi * self.M1 / (self.M1 + self.M2)
              * self.Z1 * self.Z2 * 14.4  * self.rnorm) * sn_norm
        return sn


def sn_zbl(e, Z1, Z2, M1, M2):
    """Nuclear stopping cross section from the ZBL screening function.

    Parameters:
        e (float): energy of projectile (eV)
        Z1: atomic number of atom 1
        Z2: atomic number of atom 2
        M1: mass of atom 1 (amu)
        M2: mass of atom 2 (amu)
    
    Returns:
        (float): nuclear stopping cross section (eV A^2)
    """
    rnorm = 0.4685 / (Z1**0.23 + Z2**0.23)
    enorm = (M1+M2)/M2 * 14.4 * Z1 * Z2 / rnorm
    e_norm = e / enorm
    sn_zbl_norm = (0.5 * np.log(1 + 1.1383*e_norm) 
                   / (e_norm + 0.01321*e_norm**0.21226 + 0.19593*e_norm**0.5))
    #print(f"ZBL: {e=}, {e_norm=} {sn_zbl_norm=}")
    return (4*M1/(M1+M2) * np.pi * Z1 * Z2 * 14.4 * rnorm) * sn_zbl_norm
        

def get_mass(Z):
    """Get the mass of an atom given its atomic number.

    Parameters:
        Z (int): atomic number of the atom

    Returns:
        (float): mass of the atom in amu
    """
    fname = os.path.join(os.path.dirname(__file__), 
                         '../../data/atom_data/ATOMDATA')
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.split()
            Z_val = int(items[0])
            M_val = float(items[4])
            if Z_val == Z:
                return M_val
            
    raise ValueError(f"Atomic number {Z} not found in atomic masses data.")


def get_density(Z):
    """Get the density of an atom given its atomic number.

    Parameters:
        Z (int): atomic number of the atom

    Returns:
        (float): density of the atom in 1/Å^3
    """
    fname = os.path.join(os.path.dirname(__file__), 
                         '../../data/atom_data/ATOMDATA')
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.split()
            Z_val = int(items[0])
            dens_val = float(items[7]) * 1e-24
            if Z_val == Z:
                return dens_val
            
    raise ValueError(f"Atomic number {Z} not found in atomic densities data.")


def projected_range(emax, Z1, Z2, e_vals=None, zbl=False, krc=False):
    """Calculate the projected range of a projectile in a target.

    Projected ranges are calculated between 1 eV and emax at energies either 
    specified by the user or selected by the solver.

    Parameters:
        emax (float): maximum energy of projectile (eV)
        Z1: atomic number of atom 1
        Z2: atomic number of atom 2
        e_vals (ndarray or None): energies at which to calculate projected 
            range, None for energies selected by the solver.
        zbl (bool): whether to use ZBL screening function for nuclear stopping.
        krc (bool): whether to use KrC screening function for nuclear stopping.    
    
    Returns:
        (ndarray): energies of projectile (eV)
        (ndarray): projected range (A)
    """
    M1 = get_mass(Z1)
    M2 = get_mass(Z2)
    density = get_density(Z2)

    sn = Sn(Z1, Z2, M1, M2, zbl=zbl, krc=krc)
    se = Se(Z1, Z2)

    def func(e, rp, M1, M2, density):
        sn_force = density * sn(e)
        se_force = density * se(e)
        return (1 - M2/M1 * sn_force / (2*e) * rp) / (sn_force + se_force)

    sol = solve_ivp(func, [1, emax], [0], args=(M1, M2, density), max_step=1e3,
                    t_eval=e_vals)
    return sol.t, sol.y[0, :]


def plot_sn_se(Z1, Z2, M1, M2, p1, p2):
    """Plot nuclear and electronic stopping cross sections.

    Parameters:
        Z1 (int): atomic number of first atom
        Z2 (int): atomic number of second atom
        M1 (float): mass of first atom (amu)
        M2 (float): mass of second atom (amu)
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
    """
    se = Se(Z1, Z2)
    sn = Sn(Z1, Z2, M1, M2)

    energies = np.linspace(1, 300, 300)
    se_vals = np.array([se(e*1000) for e in energies])
    sn_vals = np.array([sn(e*1000) for e in energies])
    sn_zbl_vals = np.array([sn_zbl(e*1000, Z1, Z2, M1, M2) for e in energies])

    import matplotlib.pyplot as plt
    plt.plot(energies, se_vals, 'b-', label='Se')
    plt.plot(energies, sn_vals, 'r-', label='Sn')
    plt.plot(energies, sn_zbl_vals, 'r--', label='Sn (ZBL)')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Stopping cross section (eV A$^2$)')
    plt.title(f'Stopping cross sections for Z1={Z1}, Z2={Z2}')
    plt.legend()
    plt.grid(True, which='both', ls='--')
    plt.show()        


def tab_rp(Z1, Z2, zbl=False, krc=False):
    """Tabulate projected range as a function of energy.

    Parameters:
        Z1 (int): atomic number of first atom
        Z2 (int): atomic number of second atom
        zbl (bool): whether to use ZBL potential.
        krc (bool): whether to use KRC potential.
    """
    e_vals_10 = np.logspace(1, 7, 7)
    e_vals = np.sort(np.concatenate((e_vals_10, 0.3*e_vals_10)))
    emax = e_vals[-1]
    directory = os.path.join(os.path.dirname(__file__), 
                             f"../../data/nuclear_scattering")
    if zbl:
        fname = f"rp_table_zbl_{Z1:02d}_{Z2:02d}.txt"
        energies, rp_vals = projected_range(emax, Z1, Z2,
                                            e_vals=e_vals, zbl=True)
        potname = "zbl"
    elif krc:
        fname = f"rp_table_krc_{Z1:02d}_{Z2:02d}.txt"
        energies, rp_vals = projected_range(emax, Z1, Z2,
                                            e_vals=e_vals, krc=True)
        potname = "krc"
    else:
        fname = f"rp_table_nlhlin_{Z1:02d}_{Z2:02d}.txt"
        energies, rp_vals = projected_range(emax, Z1, Z2,e_vals=e_vals)
        potname = "nlhlin"
    fname = os.path.join(directory, 
                         f"{potname}/rp_tables",
                         f"rp_table_{potname}_{Z1:02d}_{Z2:02d}.txt")
    np.savetxt(fname, np.column_stack((energies, rp_vals)), fmt='%.3e',
               header='Energy(eV) ProjectedRange(A)')


def tab_all_rp(zbl=False, krc=False):
    """Tabulate projected range for all pairs of atomic numbers.

    Parameters:
        zbl (bool): whether to use ZBL potential.
        krc (bool): whether to use KrC potential.
"""
    for Z1 in range(1, 93):
        for Z2 in range(1, 93):
            tab_rp(Z1, Z2, zbl=zbl, krc=krc)


def plot_rp(Z1, Z2, coarse=False):
    """Plot projected range as a function of energy.

    Parameters:
        Z1 (int): atomic number of first atom
        Z2 (int): atomic number of second atom
    """
    fig = plt.figure()
    
    if coarse:
        e_vals_10 = np.logspace(1, 7, 7)
        e_vals = np.sort(np.concatenate((e_vals_10, 0.3*e_vals_10)))
        emax = e_vals[-1]
    else:
        e_vals = None
        emax = 1e7
    energies, rp_vals = projected_range(emax, Z1, Z2, e_vals=e_vals)
    plt.loglog(energies, rp_vals, 'b-', label='NLHlin')

    energies, rp_vals = projected_range(emax, Z1, Z2, e_vals=e_vals, 
                                        zbl=True)
    plt.loglog(energies, rp_vals, 'r:', label='ZBL')

    #energies, rp_vals = projected_range(emax, Z1, Z2, e_vals=e_vals, 
    #                                    krc=True)
    #plt.loglog(energies, rp_vals, 'k-.', label='KrC', zorder=-1)

    if Z1 == 33 and Z2 == 14:
        fname = os.path.join(os.path.dirname(__file__),
                            'SRIM_output/Arsenic_in_Silicon.txt')
        with open(fname, 'r') as f:
            lines = f.readlines()
            energies_srim = []
            rp_srim = []
            for line in lines[23:]:
                if line.startswith('#'):
                    continue
                if line.startswith('----'):
                    break
                line = line.replace(',', '.')
                items = line.split()
                if items[1] == 'keV':
                    energies_srim.append(float(items[0]) * 1e3)
                elif items[1] == 'MeV':
                    energies_srim.append(float(items[0]) * 1e6)
                else:
                    energies_srim.append(float(items[0]))
                if items[5] == 'um':
                    rp_srim.append(float(items[4]) * 1e4)
                else:
                    rp_srim.append(float(items[4]))

        plt.loglog(energies_srim, rp_srim, 'g--', label='SRIM')

    plt.xlabel('Energy (eV)')
    plt.ylabel('Projected range (A)')
    text = fr'Z$_1$={Z1}, Z$_2$={Z2}'
    plt.text(0.95, 0.05, text, 
             horizontalalignment='right', verticalalignment='bottom', 
             transform=plt.gca().transAxes)
    plt.legend()
    #plt.grid(True, which='both', ls='--')  
    plt.xlim(1, 1e7)
    plt.ylim(1, 1e5) 
    #plt.xticks([1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    plt.minorticks_on()
    plt.tight_layout()
    plt.show()

    fname = os.path.join(os.path.dirname(__file__), 
                        f"figs/rp_{atom[Z1]}_{atom[Z2]}.pdf")
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")



def plot_rp_ratio(Z2=None, zbl=True, krc=False):
    """Plot the ratio of projected ranges between NLHlin and ZBL or KrC.

    Parameters:
        Z2 (int): atomic number of second atom, None for Z1=Z2
        zbl (bool): whether to use ZBL potential for comparison
        krc (bool): whether to use KrC potential for comparison instead of ZBL
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.ticker import MaxNLocator

    if zbl == krc:
        raise ValueError("plot_rp_ratio: exactly one of zbl and krc must "
                         "be True")
    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})
    fig = plt.figure()

    if krc:
        potname = "krc"
    else:
        potname = "zbl"

    if Z2 is None:
        Z2_text = "homo"
    else:
        Z2_text = f"{atom[Z2]}"

    Z1_vals = range(1, 93)
    rp_ratios = []

    directory = os.path.join(os.path.dirname(__file__), 
                            f"../../data/nuclear_scattering")

    for Z1 in Z1_vals:
        if Z2 is None:
            Z2_ = Z1
        else:
            Z2_ = Z2
        fname = f'nlhlin/rp_tables/rp_table_nlhlin_{Z1:02d}_{Z2_:02d}.txt'
        fname = os.path.join(directory, fname)
        data_nlhlin = np.loadtxt(fname)
        energies_nlhlin = data_nlhlin[:, 0]
        rp_nlhlin = data_nlhlin[:, 1]
        fname = f'{potname}/rp_tables/rp_table_{potname}_{Z1:02d}_{Z2_:02d}.txt'
        fname = os.path.join(directory, fname)
        data_ref = np.loadtxt(fname)
        energies_ref = data_ref[:, 0]
        rp_ref = data_ref[:, 1]
        if energies_nlhlin.shape != energies_ref.shape:
            raise ValueError(f"Energy arrays for Z1={Z1} do not match between"
                            f" NLHlin and ZBL/KrC.")
        rp_ratio = rp_nlhlin / rp_ref
        plt.semilogx(energies_nlhlin, rp_ratio, color=cmap((Z1-1)/92), 
                     zorder=Z1)
        rp_ratios.append(rp_ratio)

    rp_ratios = np.array(rp_ratios)

    ticks = range(0, 100, 10)
    if False:
        plt.axhline(1, color='k', ls='--', zorder=100)
        bounds = np.linspace(1, 92, 92)
        cmap = mpl.cm.viridis
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=plt.gca(), 
                    label=r'atomic number Z$_1$', 
                    ticks=ticks)
        plt.text(0.95, 0.95, r'Z$_2$=' f'{Z2 if Z2 is not None else "Z$_1$"}', 
                horizontalalignment='right', verticalalignment='top',
                transform=plt.gca().transAxes, fontsize='medium')
        plt.xlim(10, 1e7)
        plt.xticks([1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
        plt.ylim(0, 5)
        plt.xlabel('Energy E (eV)')
        if krc:
            plt.ylabel(r'R$_p$(NLHlin)/R$_p$(KrC)')
        else:
            plt.ylabel(r'R$_p$(NLHlin)/R$_p$(ZBL)')
        plt.tight_layout()
        plt.show()

        fname = os.path.join(os.path.dirname(__file__), 
                            f"figs/rp_ratio_{Z2_text}_{potname}.pdf")
        fname = ask_if_save(fname)
        if fname is not None:
            fig.savefig(os.path.join(os.path.dirname(__file__), fname))
            print(f"Saved figure to {fname}")
    else:
        plt.clf()

    fig = plt.figure()

    for e in [10, 100, 1000, 10000]:
        idx = np.argmin(np.abs(energies_nlhlin - e))
        rp_ratios_e = rp_ratios[:, idx]
        plt.plot(Z1_vals, rp_ratios_e, 'o-', 
                 label=fr'E=10$^{{{int(np.log10(e))}}}$ eV')
    plt.axhline(1, color='k', ls='--')
    #plt.text(0.95, 0.5, text, 
    #         horizontalalignment='right', verticalalignment='bottom',
    #         transform=plt.gca().transAxes, fontsize='medium')
    plt.xlim(0, 93)
    plt.ylim(0, 6)
    plt.xlabel('Atomic number Z$_1$')
    plt.text(0.95, 0.05, fr'Z$_2$={Z2 if Z2 is not None else "Z$_1$"}', 
             horizontalalignment='right', verticalalignment='bottom',
             transform=plt.gca().transAxes, fontsize='medium')
    if krc:
        plt.ylabel(r'R$_p$(NLHlin)/R$_p$(KrC) at E')
    else:
        plt.ylabel(r'R$_p$(NLHlin)/R$_p$(ZBL) at E')
    plt.xticks(ticks)
    plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
    plt.grid(True, which='major', ls='--')
    plt.legend()
    plt.tight_layout()
    plt.show()

    fname = os.path.join(os.path.dirname(__file__), 
                         f"figs/rp_ratio_{Z2_text}_{potname}_at_e.pdf")
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


def plot_rp_ratio_at_e(e_vals, Z2=None, zbl=True, krc=False):
    """Plot R_p(NLHlin)/R_p(ZBL or KrC) for all atom pairs.

    Parameters:
        e_vals (list): List of energies (eV)
        Z2 (int): Atomic number of second atom
        zbl (bool): Whether to use ZBL potential
        krc (bool): Whether to use KrC potential
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    if zbl == krc:
        raise ValueError("plot_rp_ratio_at_e: exactly one of zbl and krc must "
                         "be True")

    cmap = plt.get_cmap('viridis', 92)
    plt.rcParams.update({'font.size': 14})
    fig = plt.figure()

    if krc:
        potname = "krc"
    else:
        potname = "zbl"

    ticks = range(0, 100, 10)
    bounds = np.linspace(1, 92, 92)
    cmap = mpl.cm.viridis
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                 ax=plt.gca(),
                 label=r'atomic number Z$_2$', 
                 ticks=ticks)

    Z1_vals = range(1, 93)
    if Z2 is None:
        Z2_vals = range(1, 93)
    else:
        Z2_vals = [Z2]

    directory = os.path.join(os.path.dirname(__file__), 
                             f"../../data/nuclear_scattering")

    for e in e_vals:
        for Z1 in Z1_vals:
            for Z2 in Z2_vals:
                fname = f'nlhlin/rp_tables/rp_table_nlhlin_{Z1:02d}_{Z2:02d}.txt'
                fname = os.path.join(directory, fname)
                data = np.loadtxt(fname, comments='#')
                e_vals = data[:, 0]
                rp_vals = data[:, 1]
                idx = np.argmin(np.abs(e_vals - e))
                rp = rp_vals[idx]
                fname = f'{potname}/rp_tables/rp_table_{potname}_{Z1:02d}_{Z2:02d}.txt'
                fname = os.path.join(directory, fname)
                data_ref = np.loadtxt(fname)
                e_vals_ref = data_ref[:, 0]
                rp_vals_ref = data_ref[:, 1]
                idx_ref = np.argmin(np.abs(e_vals_ref - e))
                rp_ref = rp_vals_ref[idx_ref]
                rp_ratio = rp / rp_ref
                plt.plot(Z1, rp_ratio, '.', color=cmap((Z2-1)/92), zorder=Z2)
                if rp_ratio > 2:
                    print(f"Large ratio for Z1={Z1}, Z2={Z2}: {rp_ratio:.2f}")

    plt.axhline(1, color='k', linestyle='--', zorder=100)
    plt.xlim(0, 93)
    plt.ylim(0, 9)
    plt.xlabel('Atomic number Z$_1$')
    if krc:
        plt.ylabel(r'R$_p$(NLHlin)/R$_p$(KrC) at E=' f'{e} eV')
    else:
        plt.ylabel(r'R$_p$(NLHlin)/R$_p$(ZBL) at E=' f'{e} eV')
    plt.xticks(ticks)
    plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
    plt.grid(True, which='major', ls='--')
    plt.tight_layout()
    plt.show()

    fname = os.path.join(os.path.dirname(__file__), 
                         f"figs/rp_ratio_all_{potname}_at_e.pdf")
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 14})
    import time

    #p1, p2 = 0.23, 1    # ZBL
    #p1, p2 = 1/2, 1/2   # NLHlin
    p1, p2 = 1/2, 2/3   # Firsov
    Z1, Z2 = 33, 14
    M1 = get_mass(Z1)
    #print(f"Mass of Z={Z1} is {M1} amu")
    M2 = get_mass(Z2)
    #plot_sn_se(Z1, Z2, M1, M2, p1, p2)

    #plot_rp(Z1, Z2)
    #plot_rp(Z1, Z2, coarse=True)
    #start_time = time.time()
    #tab_rp(Z1, Z2)
    #tab_rp(Z1, Z2, zbl=True)
    #tab_rp(Z1, Z2, krc=True)
    #print(f'Time taken: {time.time() - start_time:.3f} s\n')
    #tab_all_rp()
    #tab_all_rp(zbl=True)
    #tab_all_rp(krc=True)

    plot_rp_ratio(zbl=True)
    #plot_rp_ratio(Z2=74, zbl=True)
    #plot_rp_ratio(krc=True)
    #plot_rp_ratio_at_e(100, krc=False)
    #plot_rp_ratio_at_e([100], Z2=6, krc=False)
    #plot_rp_ratio_at_e(100, krc=True)
