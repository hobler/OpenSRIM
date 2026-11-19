"""
Fitting a sum-of-exponentials plus linear function to DMol data

phi(r) = sum(a*exp(-b*r)) + c - d*r/r0

where phi(0)=1 and phi(r0)=0. 

Results are plotted if plot=True. Data are saved to file dmol_coeffs.dat
if save=True. The variables plot and save are hard-coded.
"""
import os, sys
import numpy as np
from numpy import exp
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit, newton
import matplotlib.pyplot as plt


# The following flags control whether the results are plotted and/or saved to
# the file NLHlin_{vmin_relweight}.dat
plot = False
save = True

# Minimum potential in eV with relative weight for fitting (weight=1/v)
# Below vmin_relweight the weight is constant (weight=1/vmin_weight)
vmin_relweight = 3

atom = [ None,
          'H', 'He', 'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne',
         'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',  'K', 'Ca',
         'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
         'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',  'Y', 'Zr',
         'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
         'Sb', 'Te',  'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
         'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
         'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
         'Pa',  'U' ]

b1 = None
b2 = None
b3 = None
narg = len(sys.argv) - 1
if narg in (3, 5, 7) or narg >= 9:
    sys.exit(f'USAGE: python3 {os.path.basename(__file__)} '
             '[Z1 [Z2 [a1 b1 [a2 b2 [a3 b3]]]]]')
if narg >= 1:
    z1 = int(sys.argv[1])
    z1_range = range(z1, z1+1)
else:
    z1_range = range(1, 93)
if narg >= 2:
    z2 = int(sys.argv[2])
    z2_range = range(z2, z2+1)
else:
    z2_range = range(1, 93)
    
if narg >= 4:
    a1 = float(sys.argv[3])
    b1 = float(sys.argv[4])
if narg >= 6:
    a2 = float(sys.argv[5])
    b2 = float(sys.argv[6])
if narg >= 8:
    a3 = float(sys.argv[7])
    b3 = float(sys.argv[8])
    
if save:
    with open(f'NLHlin_{vmin_relweight}eV.dat', 'w') as fobj:
        fobj.write('# Z1 Z2     a1        b1       a2       b2       '
                   'a3       b3  rmax[A]    error\n')


def screen_fun_1(r, a1, b1):
    global r0, screen_0, c, d
    c = 1 - a1
    d = a1*exp(-b1*r0) + c - screen_0
    return a1*exp(-b1*r) + c - d*r/r0

def screen_fun_2(r, a1, a2, b1, b2):
    global r0, screen_0, c, d
    c = 1 - a1 - a2
    d = a1*exp(-b1*r0) + a2*exp(-b2*r0) + c - screen_0
    return a1*exp(-b1*r) + a2*exp(-b2*r) + c - d*r/r0

def screen_fun_3(r, a1, a2, a3, b1, b2, b3):
    global r0, screen_0, c, d
    c = 1 - a1 - a2 - a3
    d = a1*exp(-b1*r0) + a2*exp(-b2*r0) + a3*exp(-b3*r0) + c - screen_0
    return a1*exp(-b1*r) + a2*exp(-b2*r) + a3*exp(-b3*r) + c - d*r/r0

def screen_fun_1_prime(r, a1, b1):
    global r0, screen_0, c, d
    c = 1 - a1
    d = a1*exp(-b1*r0) + c - screen_0
    return -a1*b1*exp(-b1*r) - d/r0

def screen_fun_2_prime(r, a1, a2, b1, b2):
    global r0, screen_0, c, d
    c = 1 - a1 - a2
    d = a1*exp(-b1*r0) + a2*exp(-b2*r0) + c - screen_0
    return -a1*b1*exp(-b1*r) - a2*b2*exp(-b2*r) - d/r0

def screen_fun_3_prime(r, a1, a2, a3, b1, b2, b3):
    global r0, screen_0, c, d
    c = 1 - a1 - a2 - a3
    d = a1*exp(-b1*r0) + a2*exp(-b2*r0) + a3*exp(-b3*r0) + c - screen_0
    return -a1*b1*exp(-b1*r) - a2*b2*exp(-b2*r) - a3*b3*exp(-b3*r) - d/r0

def init(r0, screen_1, screen_2, screen_0):
    """Approximate fit through r1=r0/3, r2=2*r0/3, r0."""
    aa = 1 - 2*screen_1 + screen_2
    bb = 1 + screen_0 - screen_1 - screen_2
    cc = screen_0 + screen_1 - 2*screen_2
    B = (bb - np.sqrt(bb**2 - 4*aa*cc)) / (2*aa)
    b1 = -3/r0 * np.log(B)
    a1 = (2 + screen_0 - 3*screen_1) / (2-3*B+B**3)
    rmse = np.sqrt((((screen_fun_1(rs[mask], a1, b1)-screens[mask])
                    / sigmas[mask])**2).mean())
    nexp = 1
    if plot:
        print('Initial condition:')
        print(f'{a1:.4f} {c:.4f} | {b1:.4f} {d:.4f} | e={rmse:.4f}', '\n')
    plt.plot(rs_fine, screen_fun_1(rs_fine, a1, b1), color='C1',
             label=f'initial, error={rmse:.2%}')
    return a1, b1, rmse

def fit_1(p0):
    """Fit to a1*exp(-b1*r) + c - d*r/r0."""
    try:
        popt, pcov = curve_fit(screen_fun_1, rs, screens, p0=p0,
                               bounds=(0., (9.9999, np.inf)),
                               sigma=sigmas)
    except Exception as e:
        print(f'One-exponential fit:\n{e}\n')
        popt = None
        rmse1 = np.inf
    else:
        rmse1 = np.sqrt((((screen_fun_1(rs[mask], *popt)-screens[mask])
                          / sigmas[mask])**2).mean())
        if plot:
            a1, b1 = popt
            print('Fit using one exponential:')
            print(f'{a1:.4f} {c:.4f} | {b1:.4f} {d:.4f} | e={rmse1:.4f}')
            print(np.diag(pcov), '\n')
    return popt, rmse1
    
def fit_2(p0):
    """Fit to sum_i=1^2(ai*exp(-bi*r)) + c - d*r/r0."""
    try:
        popt, pcov = curve_fit(screen_fun_2, rs, screens, p0=p0,
                               bounds=(0., (9.9999, 9.9999, 
                                            np.inf, np.inf)),
                               sigma=sigmas)
    except Exception as e:
        print(f'Two-exponentials fit:\n{e}\n')
        popt = None
        rmse2 = np.inf
    else:
        rmse2 = np.sqrt((((screen_fun_2(rs[mask], *popt)-screens[mask])
                          / sigmas[mask])**2).mean())
        if plot:
            a1, a2, b1, b2 = popt
            print('Fit using two exponentials:')
            print(f'{a1:.4f} {a2:.4f} {c:.4f} ' +
                  f'| {b1:.4f} {b2:.4f} {d:.4f} ' +
                  f'| e={rmse2:.4f}')
            print(np.diag(pcov), '\n')
    return popt, rmse2

def fit_3(p0):
    """Fit to sum_i=1^3(ai*exp(-bi*r)) + c - d*r/r0."""
    try:
        popt, pcov = curve_fit(screen_fun_3, rs, screens, p0=p0,
                               bounds=(0., (9.9999, 9.9999, 9.9999, 
                                            np.inf, np.inf, np.inf)), 
                               sigma=sigmas)
    except Exception as e:
        print(f'Three-exponential fit:\n{e}\n')
        popt = None
        rmse3 = np.inf
    else:
        rmse3 = np.sqrt((((screen_fun_3(rs[mask], *popt)-screens[mask])
                          / sigmas[mask])**2).mean())
        if plot:
            a1, a2, a3, b1, b2, b3 = popt
            print('Fit using three exponentials:')
            print(f'{a1:.4f} {a2:.4f} {a3:.4f} {c:.4f} ' +
                  f'| {b1:.4f} {b2:.4f} {b3:.4f} {d:.4f} ' +
                  f'| e={rmse3:.4f}')
            print(np.diag(pcov), '\n')
    return popt, rmse3


for z1 in z1_range:
    for z2 in z2_range:
        if save and z2 < z1:
            continue
        # read MP2 data
        fname = os.path.join('mp2data', 
                             f'{atom[z2]}_{atom[z1]}_unc-pc-1.dat')
        if os.path.exists(fname):
            rs, screens = np.loadtxt(fname, unpack=True, usecols=(0,4))
            plt.plot(rs, screens, 'C7:', label='MP2')
        
        # read screening function
        fname = os.path.join('dmol_repulsive_screen', 
                             f'screen_{min(z1, z2)}_{max(z1, z2)}.dat')
        rs, screens = np.loadtxt(fname, unpack=True)
        if rs[-2] == rs[-1]:
            rs = np.delete(rs, -2)
            screens = np.delete(screens, -2)
        vnns = int(z1)*int(z2)*14.400/rs[1:]
        vs = vnns * screens[1:]
        mask = np.insert((vs > 1), 0, True)
        not_mask = np.logical_not(mask)
        #sigmas_min = screens[mask][-1]
        sigmas_min = screens[np.insert((vs > vmin_relweight), 0, True)][-1]
        sigmas = np.maximum(screens, sigmas_min)
        #sigmas = screens.copy()
        #sigmas[not_mask] = np.sqrt(screens[not_mask]**3/(screens[mask][-1]))
        rs_fine = np.linspace(0, 1.2*rs[-1], 100)
        plt.plot(rs[mask], screens[mask], 'C0.', label='DMol (V>1eV)')
        plt.plot(rs[not_mask], screens[not_mask], 'C0+', label='DMol (V<1eV)')
        
        # radius at which v=0
        try:
            screen_spline = CubicSpline(rs, screens)
        except Exception as e:
            print(f'Z1={z1}, Z2={z2}: Erroneous screening function data:')
            print(e, '\n')
            answer = input('Press any key to continue')
            continue
        r0 = rs[-1]
        screen_0 = 0.            

        rmse = np.inf

        if narg <= 2:
            # initial condition
            screen_1 = screen_spline(r0/3)
            screen_2 = screen_spline(2*r0/3)
            a1, b1, rmse = init(r0, screen_1, screen_2, screen_0)
            nexp = 1

        if narg <= 4:
            # fit using one exponential
            p0 = (a1, b1)
            popt, rmse1 = fit_1(p0)

            if rmse1 < 0.99*rmse:
                a1, b1 = popt
                rmse = rmse1
                nexp = 1
                plt.plot(rs_fine, screen_fun_1(rs_fine, *popt), color='C2',
                         label=f'1 term, error={rmse:.2%}')

        if narg <= 4:
            p0 = (0.5*a1, 0.5*a1, 1.1*b1, 0.9*b1)
        else:
            p0 = (a1, a2, b1, b2)

        if narg <= 6:
            # fit using two exponentials
            popt, rmse2 = fit_2(p0)

            if rmse2 < 0.99*rmse:
                a1, a2, b1, b2 = popt
                rmse = rmse2
                nexp = 2
                plt.plot(rs_fine, screen_fun_2(rs_fine, *popt),
                          color='C3', label=f'2 terms, error={rmse:.2%}')
            else:
                a2 = 0
                b2 = 1

        if narg <= 6:
            p0 = (0.5*a1, 0.5*a1, a2, 1.1*b1, 0.9*b1, b2)
        else:
            p0 = (a1, a2, a3, b1, b2, b3)
        
        if narg <= 8:
            # fit using three exponentials, splitting first exponential
            popt, rmse3 = fit_3(p0)
                
            p0 = (a1, 0.5*a2, 0.5*a2, b1, 1.1*b2, 0.9*b2)

            if rmse3 < 0.99*rmse:
                a1, a2, a3, b1, b2, b3 = popt
                rmse = rmse3
                nexp = 3
                plt.plot(rs_fine, screen_fun_3(rs_fine, *popt),
                          color='C4', label=f'3 terms, error={rmse:.2%}')
            else:
                a3 = 0
                b3 = 1

        if narg <= 6:
            # fit using three exponentials, splitting second exponential
            # (p0 has been set above)
            popt, rmse3 = fit_3(p0)

            if rmse3 < 0.99*rmse:
                a1, a2, a3, b1, b2, b3 = popt
                rmse = rmse3
                nexp = 3
                plt.plot(rs_fine, screen_fun_3(rs_fine, *popt),
                          color='C5', label=f'3 terms, error={rmse:.2%}')
        
        # if the derivative is positive, find the other root
        screen_0_prime = screen_fun_3_prime(r0, a1, a2, a3, b1, b2, b3)
        if screen_0_prime > 0:
            r0 = newton(screen_fun_3, 0., fprime=screen_fun_3_prime, 
                        args=(a1, a2, a3, b1, b2, b3))
            if plot:
                print(f'r0 corrected to {r0}')
        plt.plot(rs_fine, screen_fun_3(rs_fine, a1, a2, a3, b1, b2, b3),
                 'k--', label='final result')

        plt.legend()
        plt.axhline(color='k')
        plt.title(rf'$_{{{z1}}}${atom[z1]} : $_{{{z2}}}${atom[z2]}')
        plt.gca().set_yscale('log')        
        plt.xlabel(fr'r [$\rm\AA$]')
        plt.ylabel('screening function')
        plt.xlim(0, rs_fine[-1])
        ymin = 10.**(np.floor(np.log10(screens[-2])))
        plt.ylim(ymin, 1)
        line = (f'{z1:2d} {z2:2d} {a1:8.5f} {b1:9.5f} '
                f'{a2:8.5f} {b2:8.5f} {a3:8.5f} {b3:8.5f} '
                f'{r0:8.5f} {rmse:8.5f}')
        print(line)
        if plot:
            plt.show()
        if save:
            with open(f'NLHlin_{vmin_relweight}eV.dat', 'a') as fobj:
                fobj.write(line + '\n')
