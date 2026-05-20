"""Tests for the scattering functions."""
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
from zbl import ZBL_screen
from nlhlin import NLHlin_screen
from cm_scatter import setup, scatter_integrals
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
    Z1 = 2
    Z2 = 74
    #screen_fun = ZBL_screen(Z1, Z2, rnorm=1.0)
    screen_fun = NLHlin_screen(Z1, Z2, rnorm=1.0)
    plot_psi_and_de(e=3000, Z1=Z1, Z2=Z2, screen_fun=screen_fun)

    #compare_theta_and_tau(screen_fun)
