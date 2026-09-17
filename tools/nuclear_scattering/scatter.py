"""Tests for the scattering functions."""
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
from krc import KrC_screen
from zbl import ZBL_screen
from nlhlin import NLHlin_screen
from cm_scatter import setup, scatter_integrals
from rp import Sn
from utils import atom, ask_if_save


def md_scatter(e, p, screen_fun, dx=0.001):
    """Calculate the position and velocity vectors after a collision.

    The scattering event is treated within the plane defined by the incoming
    projectile and the recoil. The recoil's initial position is at the origin
    of the 2D coordinate system, and its initial velocity is zero. The incoming
    projectile is initially moving along the x-axis, and the impact parameter
    is along the y-axis.

    The units of length and time are Angstroms and femtoseconds, respectively.
    The screening length of the screening function must be 1 Angstrom.

    Parameters:
        e (float): kinetic energy of the incoming particle (eV)
        p (float): impact parameter (A)
        screen_fun (callable): screening function object
        dx (float): initial step size (A)

    Returns:
        (2darray): position of ion
        (2darray): position of recoil
        (2darray): velocity of ion
        (2darray): velocity of recoil
    """
    assert screen_fun.Z1 == screen_fun.Z2,  \
           "md_scatter: only implemented for homo-nuclear case"
    assert screen_fun.rnorm == 1.0,  \
           "md_scatter: screening length must be 1 A"

    Z = screen_fun.Z1
    with open(os.path.join(os.path.dirname(__file__), 
                           '../../data/atom_data/ATOMDATA'), 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            items = line.split()
            Z_val = int(items[0])
            M = float(items[4])
            if Z_val == Z:
                break

    # Initial positions in A
    x1 = np.array((-np.sqrt(5.0**2 - p**2), p))
    x2 = np.array((0.0, 0.0))

    # Initial velocities in A/fs
    v1 = np.array((0.1389*np.sqrt(e/M), 0.0))
    v2 = np.array((0.0, 0.0))

    # Further initializations
    dt = dx / v1[0]
    x1_old = x1 - v1 * dt
    x2_old = x2 - v2 * dt
    r12 = x2 - x1
    r = np.linalg.norm(r12)
    rmax = r

    # Loop over time steps
    while r <= rmax:
        screen, dscreen = screen_fun(r)
        a = Z**2 / M * 0.1389 * (screen - dscreen*r) / r**2
        a2 = a * r12 / r
        a1 = -a2
        x1_new = 2 * x1 - x1_old + a1 * dt**2
        x2_new = 2 * x2 - x2_old + a2 * dt**2
        v1 = (x1_new - x1_old) / (2*dt)
        v2 = (x2_new - x2_old) / (2*dt)

        x1_old, x1 = x1, x1_new
        x2_old, x2 = x2, x2_new
        r12 = x2 - x1
        r = np.linalg.norm(r12)

    return x1, x2, v1, v2


def plot_psi_and_de(e, Z1, Z2, screen_fun):
    """Plot scattering angle and energy transfer vs impact parameter.
    
    Exact results are compared to the impulse approximation.

    Parameters:
        e: float
            Kinetic energy of the incoming particle (eV).
        Z1: int
            Atomic number of the incoming particle.
        Z2: int
            Atomic number of the target particle.
        screen_fun: callable
            Screening function object
    """
    data = np.genfromtxt(
        os.path.join(os.path.dirname(__file__), 
                     '../../data/atom_data/ATOMDATA'),
        dtype='i8, U2, U16, i8, f8, f8, f8, f8, f8, f8, f8, f8',
        names=['Z', 'symbol', 'name', 'mass', 'M1', 'M2', 'density', 'N',
            'vF', 'Esurf', 'density_gas', 'Ngas'],
        max_rows=92
    )
    print(data[Z1-1])
    M1 = data[Z1-1]['M1']
    M2 = data[Z2-1]['M2']
    rnorm = screen_fun.rnorm
    print(rnorm)
    vcoul = 14.4 * Z1 * Z2 / rnorm
    enorm = vcoul * (M1 + M2) / M2

    p_vals = np.linspace(0.01, 5.0, 500)  # impact parameters (A)

    psi_vals = []
    de_vals = []
    for p in p_vals:
        pi_minus_theta, _ = scatter_integrals(e/enorm, p/rnorm, screen_fun)
        psi = np.arctan(
            np.sin(pi_minus_theta) / (-np.cos(pi_minus_theta) + M1/M2)
        )
        psi_vals.append(psi)
        de = 4*M1*M2/(M1+M2)**2 * e * np.cos(pi_minus_theta/2)**2
        de_vals.append(de)
    psi_vals = np.array(psi_vals)
    de_vals = np.array(de_vals)

    psi_impulse_vals = []
    de_impulse_vals = []
    for p in p_vals:
        integral = screen_fun.impulse_integral(p/rnorm)
        psi_impulse = vcoul * integral / e
        psi_impulse_vals.append(psi_impulse)
        de_impulse = M1/M2 * vcoul**2 * integral**2 / e
        de_impulse_vals.append(de_impulse)
    psi_impulse_vals = np.array(psi_impulse_vals)
    de_impulse_vals = np.array(de_impulse_vals)

    fig = plt.figure()
    plt.plot(p_vals[psi_vals>0], np.degrees(psi_vals[psi_vals>0]), 
            'C0',label='Exact scattering angle')
    plt.plot(p_vals, np.degrees(psi_impulse_vals), 
            'C0--', label='Impulse approximation')
    plt.plot(p_vals, de_vals, 
            'C1', label='Exact energy transfer')
    plt.plot(p_vals, de_impulse_vals, 
            'C1--', label='Impulse approximation')
    plt.xlabel(r'Impact parameter ($\rm\AA$)')
    plt.ylabel('Scattering angle (deg), Energy transfer (eV)')
    plt.ylim(0.1, 100)
    plt.yscale('log')
    plt.title(rf'{screen_fun.name} potential, Z1={Z1}, Z2={Z2}, E={e} eV',
              fontsize='medium')
    plt.legend()
    plt.xlim(0, 2)
    plt.tight_layout()
    plt.show()

    if screen_fun.name == "ZBL":
        potname = f"zbl_{atom[screen_fun.Z1]}_{atom[screen_fun.Z2]}"
    elif screen_fun.name == "NLHlin":
        potname = f"nlhlin_{atom[screen_fun.Z1]}_{atom[screen_fun.Z2]}"
    else:
        potname = screen_fun.name
    fname = f"figs/impulse_{potname}.pdf"
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


def plot_avg_and_min_psi(psi_avg, Z1, Z2, screen_fun):
    """Plot average cumulated and minimum scattering angle vs reduced energy.
    
    Parameters:
        psi_avg: float
            Average cumulated scattering angle (deg)    
        e: float
            Kinetic energy of the incoming particle (eV).
        Z1: int
            Atomic number of the incoming particle.
        Z2: int
            Atomic number of the target particle.
        screen_fun: callable
            Screening function object
    """
    data = np.genfromtxt(
        os.path.join(os.path.dirname(__file__), 
                     '../../data/atom_data/ATOMDATA'),
        dtype='i8, U2, U16, i8, f8, f8, f8, f8, f8, f8, f8, f8',
        names=['Z', 'symbol', 'name', 'mass', 'M1', 'M2', 'density', 'N',
            'vF', 'Esurf', 'density_gas', 'Ngas'],
        max_rows=92
    )
    print(data[Z1-1])
    M1 = data[Z1-1]['M1']
    M2 = data[Z2-1]['M2']
    rnorm = screen_fun.rnorm
    print(rnorm)
    vcoul = 14.4 * Z1 * Z2 / rnorm
    enorm = vcoul * (M1 + M2) / M2
    print(enorm)

    eps_vals = np.logspace(-4, 2, 61)  # reduced energies
    psi_min_vals = []
    for eps in np.logspace(-4, 2, 61):
        sn = Sn(Z1, Z2, M1, M2, 
                zbl=screen_fun.name == "ZBL", 
                krc=screen_fun.name == "KrC")
        e = eps * enorm
        pmax = np.sqrt(M2/M1 * sn(e) / (np.pi*e)) / np.radians(psi_avg)
        pmax_ZBL = 2*rnorm * np.sqrt(np.log(1+eps) / (0.02*(1+M1/M2)**2*(eps**2 + 0.1*eps**1.38)))
        print(f'ε={eps:.4e}, pmax={pmax:.4e} A, pmax_ZBL={pmax_ZBL:.4e} A')
        integral = screen_fun.impulse_integral(pmax/rnorm)        
        psi_min = np.degrees(M2 / (M1 + M2) * integral / eps)
        psi_min_vals.append(psi_min)
    psi_min_vals = np.array(psi_min_vals)

    plt.axhline(psi_avg, color='k', ls='--', 
                label='Average cumulated scattering angle')
    plt.plot(eps_vals, psi_min_vals, 'C1', label='Minimum scattering angle')
    plt.xscale('log')
    plt.xlabel(r'Reduced energy $\epsilon$')
    plt.ylabel('Scattering angle (deg)')
    plt.title(rf'{screen_fun.name} potential, Z1={Z1}, Z2={Z2}, '
              rf'$\psi_{{avg}}$={psi_avg:.2f} deg',
              fontsize='medium')
    plt.legend()
    plt.tight_layout()
    plt.show()


def compare_pmax(t_min, psi_min, psi_avg, Z1, Z2, screen_fun):
    """Plot maximum impact parameters vs reduced energy.

    This function plots the maximum impact parameters as a function of reduced 
    energy, applying three different criteria: 
    
    - the minimum scattering angle, 
    - the average cumulated scattering angle, and 
    - the minimum energy transfer.

    Parameters:
        psi_min: float
            Minimum scattering angle (deg)
        psi_avg: float
            Average cumulated scattering angle (deg)    
        e: float
            Kinetic energy of the incoming particle (eV).
        Z1: int
            Atomic number of the incoming particle.
        Z2: int
            Atomic number of the target particle.
        screen_fun: callable
            Screening function object
    """
    data = np.genfromtxt(
        os.path.join(os.path.dirname(__file__), 
                     '../../data/atom_data/ATOMDATA'),
        dtype='i8, U2, U16, i8, f8, f8, f8, f8, f8, f8, f8, f8',
        names=['Z', 'symbol', 'name', 'mass', 'M1', 'M2', 'density', 'N',
            'vF', 'Esurf', 'density_gas', 'Ngas'],
        max_rows=92
        )
    print(data[Z1-1])
    M1 = data[Z1-1]['M1']
    M2 = data[Z2-1]['M2']
    N = data[Z2-1]['N'] / 1e24  # convert to atoms/A^3
    gamma = 4*M1*M2 / (M1 + M2)**2
    rnorm = screen_fun.rnorm
    print(rnorm)
    vcoul = 14.4 * Z1 * Z2 / rnorm
    enorm = vcoul * (M1 + M2) / M2
    print(enorm)

    sn = Sn(Z1, Z2, M1, M2, 
            zbl=screen_fun.name == "ZBL", 
            krc=screen_fun.name == "KrC")

    eps_vals = np.logspace(-6, 4, 101)  # reduced energies
    pmax_psi_avg_vals = []
    pmax_psi_avg_ZBL_vals = []
    pmax_de_ZBL_vals = []
    for eps in eps_vals:
        e = eps * enorm
        pmax_psi_avg_vals.append(
            np.sqrt(M2/M1 * sn(e) / (np.pi*e)) / np.radians(psi_avg)
            )
        pmax_psi_avg_ZBL_vals.append(
            2*rnorm * np.sqrt(
                np.log(1+eps) / (0.02 * (1+M1/M2)**2 * (eps**2 + 0.1*eps**1.38))
                )
            )
        xi = np.sqrt(eps * t_min/enorm / gamma)
        pmax_de_ZBL_vals.append(
            rnorm / (xi + xi**0.5 + 0.125*xi**0.1)
            )

    pmax_vals = np.linspace(0.02, 4, 200)
    eps_psi_min_vals = []
    eps_de_min_vals = []
    for pmax in pmax_vals:
        integral = screen_fun.impulse_integral(pmax/rnorm)
        print(f'pmax={pmax:.3f}, integral={integral:.3f}')
        eps_psi_min_vals.append(
            vcoul * integral / np.radians(psi_min) / enorm
            )
        eps_de_min_vals.append(
            M1/M2 * vcoul**2 * integral**2 / t_min / enorm
            )
    eps_psi_min_vals = np.array(eps_psi_min_vals)
    eps_de_min_vals = np.array(eps_de_min_vals)

    # plot over reduced energy
    plt.plot(eps_vals, pmax_psi_avg_vals, 'C0', 
             label=r'$\psi_\mathrm{avg}$=' 
                   + fr'{psi_avg:.2f}$\degree$ (B-W)')  # Bohr-Williams
    plt.plot(eps_vals, pmax_psi_avg_ZBL_vals, 'C0--', 
             label=r'$\psi_\mathrm{avg}$=' + r'5.78$\degree$ (SRIM)')
    plt.plot(eps_psi_min_vals, pmax_vals, 'C1',
             label=r'$\psi_\mathrm{min}$=' + fr'{psi_min:.2f}$\degree$')
    plt.plot(eps_de_min_vals, pmax_vals, 'C2',
             label=r'$T_\mathrm{min}$=' + f'{t_min:.2f} eV')
    plt.plot(eps_vals, pmax_de_ZBL_vals, 'C2--',
             label=r'$T_\mathrm{min}$=' + f'{t_min:.2f} eV (SRIM)')
    plt.xlim(1e-4, 1e2)
    plt.xscale('log')
    plt.xlabel(r'Reduced energy $\epsilon$')
    plt.ylim(0, 4)
    plt.ylabel(r'Maximum impact parameter $p_\mathrm{max}$ ($\rm\AA$)')
    plt.title(rf'{screen_fun.name} potential, Z$_1$={Z1}, Z$_2$={Z2}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # plot over energy in eV
    e_psi_min_vals = eps_psi_min_vals * enorm
    e_de_min_vals = eps_de_min_vals * enorm
    e_vals = eps_vals * enorm

    # with and without considering t_min
    for flag in (False, True):
        # SRIM range calculations
        pmaxmax_psi_avg_vals = np.full_like(e_vals, 0.5642*N**(-1/3))
        plt.plot(e_vals, np.minimum(pmaxmax_psi_avg_vals, pmax_psi_avg_ZBL_vals), 
                'C1', label=r'$\psi_{avg}$=' + r'5.78$\degree$ (B-W/SRIM)')
        plt.plot(e_vals, pmaxmax_psi_avg_vals, 'C1:')
        plt.plot(e_vals, pmax_psi_avg_ZBL_vals, 'C1:') 

        if flag:
            # SRIM defect calculations
            pmaxmax_de_vals = np.full_like(e_vals, 1.2407*N**(-1/3))
            plt.plot(e_vals, np.minimum(pmaxmax_de_vals, pmax_de_ZBL_vals), 
                     'C0', label=r'$T_\mathrm{min}$=' + f'{t_min} eV (SRIM)')
            plt.plot(e_vals, pmaxmax_de_vals, 'C0:')
            plt.plot(e_vals, pmax_de_ZBL_vals, 'C0:')

            # SRIM maximum
            plt.plot(e_vals, np.maximum(
                np.minimum(pmaxmax_psi_avg_vals, pmax_psi_avg_ZBL_vals), 
                np.minimum(pmaxmax_de_vals, pmax_de_ZBL_vals)
                ), 
                'k--', lw=2, label='SRIM maximum') 

        # OpenSRIM angle and energy criteria
        plt.plot(e_psi_min_vals, pmax_vals, 'C3',
                label=r'$\psi_\mathrm{min}$=' + fr'{psi_min}$\degree$ (OpenSRIM)')
        if flag:
            plt.plot(e_de_min_vals, pmax_vals, 'C2',
                    label=r'$T_\mathrm{min}$=' + f'{t_min} eV (OpenSRIM)')
            plt.plot(np.maximum(e_psi_min_vals, e_de_min_vals), pmax_vals, 'k', 
                    lw=2, label='OpenSRIM maximum') 

        plt.xlim(1e1, 1e7)
        plt.xscale('log')
        plt.xlabel(r'Energy $E$ (eV)')
        plt.ylim(0, 3)
        plt.ylabel(r'Maximum impact parameter $p_\mathrm{max}$ ($\rm\AA$)')
        plt.title(rf'{screen_fun.name} potential, {atom[Z1]} in {atom[Z2]}',
                fontsize='medium')
        plt.legend(loc='upper right', fontsize='small')
        plt.tight_layout()
        plt.show()

    # Same plot without ZBL results
    # psi_avg
    pmaxmax_psi_avg_vals = np.full_like(e_vals, 0.5642*N**(-1/3))
    plt.plot(e_vals, np.minimum(pmaxmax_psi_avg_vals, pmax_psi_avg_vals), 
            'C1', label=r'$\psi_{avg}$=' + r'5.78$\degree$ (+cutoff)')
    plt.plot(e_vals, pmaxmax_psi_avg_vals, 'C1:')
    plt.plot(e_vals, pmax_psi_avg_vals, 'C1:') 

    # psi_min
    plt.plot(e_psi_min_vals, pmax_vals, 'C3',
            label=r'$\psi_\mathrm{min}$=' + fr'{psi_min}$\degree$')

    # de_min
    plt.plot(e_de_min_vals, pmax_vals, 'C2',
            label=r'$T_\mathrm{min}$=' + f'{t_min} eV')

    # psi_avg + de_min
    pmaxmax_de_vals = np.full_like(e_vals, 1.2407*N**(-1/3))
    pmax_de_min_vals_interp = np.interp(np.log(e_vals), 
                                        np.log(e_de_min_vals[::-1]), 
                                        pmax_vals[::-1])
    plt.plot(e_vals, np.maximum(
        np.minimum(pmaxmax_psi_avg_vals, pmax_psi_avg_vals), 
        pmax_de_min_vals_interp), 
        'k--', lw=2, 
        label=r'combining $\psi_\mathrm{avg}$ and $T_\mathrm{min}$')

    # psi_min + de_min
    plt.plot(np.maximum(e_psi_min_vals, e_de_min_vals), pmax_vals, 'k', 
            lw=2, label=r'combining $\psi_\mathrm{min}$ and $T_\mathrm{min}$') 

    plt.xlim(1e1, 1e7)
    plt.xscale('log')
    plt.xlabel(r'Energy $E$ (eV)')
    plt.ylim(0, 3)
    plt.ylabel(r'Maximum impact parameter $p_\mathrm{max}$ ($\rm\AA$)')
    plt.title(rf'{screen_fun.name} potential, {atom[Z1]} in {atom[Z2]}',
            fontsize='medium')
    plt.legend(loc='upper right', fontsize='small')
    plt.tight_layout()
    plt.show()


def compare_theta_and_tau(screen_fun):
    """Compare scattering angle and time integral between MD and numerical."""
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    p_vals = np.linspace(0.0, 3.0, 300)
    #p_vals = np.linspace(0.0, 1e-7, 10)
    for i, eps in enumerate((1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)):
    #for i, eps in enumerate((1e-4, 100.0)):
        e = eps * 2 * 14.4 * screen_fun.Z1**2
        print(f'ε={eps}, e={e} eV')
        
        # MD simulations
        theta_vals = []
        tau_vals = []
        mask = np.ones_like(p_vals, dtype=bool)
        for p in p_vals:
            x1, x2, v1, v2 = md_scatter(e, p, screen_fun)
            if v1[1] == 0.0:
                mask[p_vals == p] = False
                continue
            theta = np.pi - 2.0 * np.arctan2(-v2[1], v2[0])
            theta_vals.append(np.degrees(theta))
            x1_hat = x1[0] + (p - x1[1]) * v1[0] / v1[1]
            tau = -x1_hat  # this is valid for the homo-nuclear case only
            tau_vals.append(tau)
        ax1.plot(p_vals[mask], theta_vals, label=f'ε={eps}')
        ax2.plot(p_vals[mask], tau_vals, label=f'ε={eps}')
        ax3.plot(theta_vals, tau_vals, label=f'ε={eps}')
        ax3.plot(theta_vals, p_vals[mask]*np.tan(np.radians(theta_vals)/2), 
                 'k:')

        # numeric values
        theta_vals = []
        tau_vals = []
        for p in p_vals[mask]:
            pi_minus_theta, tau = scatter_integrals(eps, p, screen_fun)
            theta = np.pi - pi_minus_theta
            theta_vals.append(np.degrees(theta))
            tau_vals.append(tau)
        ax1.plot(p_vals[mask], theta_vals, f'C{i}--')
        ax2.plot(p_vals[mask], tau_vals, f'C{i}--')
        ax2.plot(p_vals[mask], p_vals[mask]*np.tan(np.radians(theta_vals)/2), 
                 f'C{i}:')
        ax3.plot(theta_vals, tau_vals, f'C{i}--')
    ax1.set_xlabel(r'Impact parameter ($\rm\AA$)')
    ax1.set_ylabel('Scattering angle (deg)')
    ax2.set_xlabel(r'Impact parameter ($\rm\AA$)')
    ax2.set_ylabel('Time integral (A)')
    ax3.set_xlabel('Scattering angle (deg)')
    ax3.set_ylabel('Time integral (A)')
    ax1.legend()
    fig.suptitle(rf'{screen_fun.name} potential, Z1=Z2={Z1}', fontsize='medium')
    plt.tight_layout()
    plt.show()

    if screen_fun.name == "ZBL":
        potname = "zbl"
    elif screen_fun.name == "NLHlin":
        potname = f"nlhlin_{atom[screen_fun.Z1]}_{atom[screen_fun.Z2]}"
    else:
        potname = screen_fun.name
    fname = f"figs/theta_tau_{potname}.pdf"
    fname = ask_if_save(fname)
    if fname is not None:
        fig.savefig(os.path.join(os.path.dirname(__file__), fname))
        print(f"Saved figure to {fname}")


if __name__ == "__main__":
    setup(n_absc=4)
    Z1 = 15
    Z2 = 14
    #screen_fun = ZBL_screen(Z1, Z2)
    screen_fun = NLHlin_screen(Z1, Z2) #, rnorm=1.0)
    
    #plot_psi_and_de(e=3000, Z1=Z1, Z2=Z2, screen_fun=screen_fun)
    #plot_avg_and_min_psi(psi_avg=5.0, Z1=Z1, Z2=Z2, screen_fun=screen_fun)
    t_min = 15
    psi_min = 1
    psi_avg = 5.78
    compare_pmax(t_min, psi_min, psi_avg, Z1, Z2, screen_fun)

    #compare_theta_and_tau(screen_fun)
