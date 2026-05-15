"""Fit the Sn and Qn values by an analytical function."""
import numpy as np
from scipy.optimize import curve_fit


def sn_fit_func(eps, a, b, c, d):
    """Analytical function to fit Sn values.

    Parameters:
        eps (float or np.ndarray): Energy value(s)
        a, b, c, d (float): Fit parameters
    Returns:
        (float or ndarray): Fitted Sn value(s)
    """
    return np.log(1 + a*eps) / (2*(eps + b*eps**c + d*eps**0.5))


def sn_fit_func_corrb(eps, sn, a, c, d):
    """Adjust b so sn_fit_func matches the sn at eps."""
    b = ((np.log(1 + a*eps) / (2*sn)) - eps - d*eps**0.5) / eps**c
    return b


def fit_sn(p1, p2, fname, plot=True, Z1=None, Z2=None):
    """Fit Sn values from a file using the analytical function.

    Parameters:
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
        fname (str): Path to the file containing energy and Sn values.

    Returns:
        popt (tuple): Optimal values for the fit parameters.
        rms_relerr (float): Root mean square relative error (%) of the fit.
        max_relerr (float): Maximum relative error (%) of the fit.
        Z1, Z2 (int): atomic numbers of the two atoms (if provided)
    """
    fname = f"data/p{p1}_{round(p2, 2)}/" + fname
    data = np.loadtxt(fname)
    energies = data[:, 0]
    sn_data = data[:, 1]
    #fit_energies = energies[energies < 31]
    fit_energies = energies
    fit_sn_data = sn_data[:len(fit_energies)]

    # Initial guess for the parameters a, b, c, d (ZBL parameters)
    initial_guess = [1.1383, 0.01321, 0.21226, 0.19593]
    #initial_guess = [0.94127, 0.09252, 0.40137, 0.15247]  # for 2 38
    #initial_guess = [0.74490,  0.58785,  0.44206, -0.54043]  # for 2 52
    initial_guess = [0.99600,  0.00007, -0.19004,  0.38612]  # for 4 70
    # Adjust initial guess for b based on first data point
    #initial_guess[1] = (initial_guess[0] / (2*sn_data[0]) 
    #                    * energies[0]**(1-initial_guess[2]))
    #initial_guess[1] = (((np.log(1 + initial_guess[0]*energies[0]) 
    #                      / (2*sn_data[0]))
    #                     - energies[0] - initial_guess[3]*energies[0]**0.5) 
    #                    / energies[0]**initial_guess[2])
    e_val = energies[0]
    sn_val = sn_data[0]
    a, b, c, d = initial_guess
    initial_guess[1] = sn_fit_func_corrb(e_val, sn_val, a, c, d)

    initial_sn_values = sn_fit_func(energies, *initial_guess)
    relerr_initial = 100 * (initial_sn_values - sn_data) / sn_data
    rms_relerr_initial = np.sqrt(np.mean(relerr_initial**2))
    max_relerr_initial = np.max(np.abs(relerr_initial))

    # Fit the data
    try:
        popt, pcov = curve_fit(sn_fit_func, fit_energies, fit_sn_data, 
                            sigma=fit_sn_data, p0=initial_guess)
    except RuntimeError as e:
        print(fname, ':')
        print(f"Fit failed: {e}")  # Return initial guess if fit fails
        #return initial_guess, rms_relerr_initial, max_relerr_initial
        popt = (a, b, c, d)
    fitted_sn_values = sn_fit_func(energies, *popt)

    relerr = 100 * (fitted_sn_values - sn_data) / sn_data
    rms_relerr = np.sqrt(np.mean(relerr**2))
    max_relerr = np.max(np.abs(relerr))
    
    if plot:
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 14})
        print("Fitted Sn parameters:")
        print(f"a = {popt[0]}")
        print(f"b = {popt[1]}")
        print(f"c = {popt[2]}")
        print(f"d = {popt[3]}")
        print(f"Initial guess RMS/maximum relerr(%)="
              f"{rms_relerr_initial:5.2f}/{max_relerr_initial:5.2f}")
        #print(f"relerr(%)={relerr}")
        print(f"RMS/maximum error(%): {rms_relerr:5.2f}/{max_relerr:5.2f}")
        plt.loglog(energies, sn_data, 'b.', label='Data')
        plt.loglog(energies, fitted_sn_values, 'r-', label='Fit')
        #plt.loglog(energies, initial_sn_values, 'g--', label='Initial Guess')
        plt.loglog(energies, np.log(1+energies)/(2*energies), 'k:', 
                   label=r'$\ln \varepsilon / 2 \varepsilon$')
        plt.xlim(1e-5, 1e4)
        plt.ylim(1e-3, 1)
        plt.xlabel(r'Reduced energy $\varepsilon$')
        plt.ylabel(r'Reduced nuclear stopping $s_n$')
        plt.legend()

        if Z1 is not None and Z2 is not None:
            text = fr"Z$_1$={Z1}, Z$_2$={Z2}"
            plt.text(0.05, 0.8, text, 
                    horizontalalignment='left', verticalalignment='top',
                    transform=plt.gca().transAxes, fontsize='medium')

        text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
        if p1 == 1:
            text += fr'(Z$_1$+Z$_2$)'
        else:
            text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
        if p2 != 1:
            text += fr'$^{{{p2:.2f}}}$'
        plt.text(0.5, 0.05, text, 
                horizontalalignment='center', verticalalignment='bottom',
                transform=plt.gca().transAxes, fontsize='medium')

        plt.tight_layout()
        plt.show()

    return popt, rms_relerr, max_relerr


def qn_fit_func(eps, a, b, c, d):
    """Analytical function to fit Sn values.

    Parameters:
        eps (float or np.ndarray): Energy value(s)
        a, b, c, d (float): Fit parameters
    Returns:
        (float or ndarray): Fitted Sn value(s)
    """
    return 1 / (4 + a*eps**b + c*eps**d)


def fit_qn(p1, p2, fname, plot=True, Z1=None, Z2=None):
    """Fit Qn values from a file using the analytical function.

    Parameters:
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
        fname (str): Path to the file containing energy and Qn values.
        Z1, Z2 (int): atomic numbers of the two atoms (if provided)

    Returns:
        popt (tuple): Optimal values for the fit parameters.
        rms_relerr (float): Root mean square relative error (%) of the fit.
        max_relerr (float): Maximum relative error (%) of the fit.
    """
    fname = f"data/p{p1}_{round(p2, 2)}/" + fname
    data = np.loadtxt(fname)
    energies = data[:, 0]
    qn_data = data[:, 2]
    mask = np.logical_not(np.isnan(qn_data))
    energies = energies[mask]  # Excludes nan values
    qn_data = qn_data[mask]

    # Initial guess for the parameters a, b, c, d (ZBL parameters)
    initial_guess = [0.19676, -1.6991, 6.5841, -1.0494]
    # Adjust initial guess for b based on first data point
    initial_guess[0] = 1 / (qn_data[0] * energies[0]**initial_guess[1])
    initial_qn_values = qn_fit_func(energies, *initial_guess)
    relerr_initial = 100 * (initial_qn_values - qn_data) / qn_data
    rms_relerr_initial = np.sqrt(np.mean(relerr_initial**2))
    max_relerr_initial = np.max(np.abs(relerr_initial))

    # Fit the data
    try:
        popt, pcov = curve_fit(qn_fit_func, energies, qn_data, 
                            sigma=qn_data, p0=initial_guess)
    except (RuntimeError, ValueError) as e:
        print(fname, ':')
        print(f"Fit failed: {e}")  # Return initial guess if fit fails
        return initial_guess, rms_relerr_initial, max_relerr_initial
    fitted_qn_values = qn_fit_func(energies, *popt)

    relerr = 100 * (fitted_qn_values - qn_data) / qn_data
    rms_relerr = np.sqrt(np.mean(relerr**2))
    max_relerr = np.max(np.abs(relerr))
    
    if plot:
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 14})
        print("Fitted Qn parameters:")
        print(f"a = {popt[0]}")
        print(f"b = {popt[1]}")
        print(f"c = {popt[2]}")
        print(f"d = {popt[3]}")
        print(f"Initial guess RMS/maximum relerr(%)="
              f"{rms_relerr_initial:5.2f}/{max_relerr_initial:5.2f}")
        #print(f"relerr(%)={relerr}")
        print(f"RMS/maximum error(%): {rms_relerr:5.2f}/{max_relerr:5.2f}")
        plt.loglog(energies, qn_data, 'b.', label='Data')
        plt.loglog(energies, fitted_qn_values, 'r-', label='Fit')
        #plt.loglog(energies, initial_qn_values, 'g--', label='Initial Guess')
        plt.xlim(1e-5, 1e4)
        plt.ylim(1e-8, 1)
        plt.xlabel(r'Reduced energy $\varepsilon$')
        plt.ylabel(r'Reduced nuclear straggling $q_n$')
        plt.legend(loc="center right")

        if Z1 is not None and Z2 is not None:
            text = fr"Z$_1$={Z1}, Z$_2$={Z2}"
            plt.text(0.05, 0.95, text, 
                    horizontalalignment='left', verticalalignment='top',
                    transform=plt.gca().transAxes, fontsize='medium')

        text = r'a$_\mathrm{I}$=0.4685$\rm\AA$/'
        if p1 == 1:
            text += fr'(Z$_1$+Z$_2$)'
        else:
            text += fr'(Z$_1^{{{p1:.2f}}}$+Z$_2^{{{p1:.2f}}}$)'
        if p2 != 1:
            text += fr'$^{{{p2:.2f}}}$'
        plt.text(0.5, 0.05, text, 
                horizontalalignment='center', verticalalignment='bottom',
                transform=plt.gca().transAxes, fontsize='medium')

        plt.tight_layout()
        plt.show()

    return popt, rms_relerr, max_relerr


def tab_sn_fit(p1, p2):
    """Tabulate parameters for Sn fit of all Z1-Z2 combinations.
    
    Parameters:
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
    """
    p1 = round(p1, 2)
    p2 = round(p2, 2)
    directory = f"data/p{p1}_{p2}"
    with open(f"{directory}/sn_fit_params.txt", "w") as f:
        f.write("# Z1 Z2 a b c d rms_err(%) max_err(%)\n")
        for Z1 in range(1, 93):
            for Z2 in range(Z1, 93):
                fname = f"sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt"
                #print(f"Fitting Sn for Z1={Z1}, Z2={Z2} from {fname}")
                popt, rms_relerr, max_relerr = fit_sn(p1, p2, fname, plot=False)
                #print(f"Z1={Z1}, Z2={Z2}, a={popt[0]}, b={popt[1]}, "
                #      f"c={popt[2]}, d={popt[3]}")
                f.write(f"{Z1:2d} {Z2:2d} {popt[0]:8.5f} {popt[1]:8.5f}"
                        f" {popt[2]:8.5f} {popt[3]:8.5f} "
                        f"{rms_relerr:5.2f} {max_relerr:5.2f}\n")


def tab_qn_fit(p1, p2):
    """Tabulate parameters for Qn fit of all Z1-Z2 combinations.
    
    Parameters:
        p1 (float): power in RNORM formula.
        p2 (float): power in RNORM formula.
    """
    p1 = round(p1, 2)
    p2 = round(p2, 2)
    directory = f"data/p{p1}_{p2}"
    with open(f"{directory}/qn_fit_params.txt", "w") as f:
        f.write("# Z1 Z2 a b c d rms_err(%) max_err(%)\n")
        for Z1 in range(1, 93):
            for Z2 in range(Z1, 93):
                fname = f"sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt"
                #print(f"Fitting Qn for Z1={Z1}, Z2={Z2} from {fname}")
                popt, rms_relerr, max_relerr = fit_qn(p1, p2, fname, plot=False)
                #print(f"Z1={Z1}, Z2={Z2}, a={popt[0]}, b={popt[1]}, "
                #      f"c={popt[2]}, d={popt[3]}")
                f.write(f"{Z1:2d} {Z2:2d} {popt[0]:8.5f} {popt[1]:8.5f}"
                        f" {popt[2]:8.5f} {popt[3]:8.5f} "
                        f"{rms_relerr:5.2f} {max_relerr:5.2f}\n")


def plot_sn_ratio_at_e(e, p1, p2, Z1=None, Z2=None, krc=False):
    """Plot Sn ratios as a function of Z at specific energies.
    
    The ratio Sn(NLHlin)/Sn(ZBL) or Sn(NLHlin)/Sn(KrC) is plotted.
    Sn values are taken from the fitted parameters for all Z1-Z2 combinations.
    Currently only implemented for Z2=Z1.
    
    Parameters:
        e (array-like): Energy value to plot Sn at.
        p1 (float): power in RNORM formula. 
        p2 (float): power in RNORM formula.
        Z1 (int or None): atomic number of first atom, None for all.
        Z2 (int or None): atomic number of second atom, None for Z1.
        krc (bool): whether to plot ratio to KrC instead of ZBL.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    plt.rcParams.update({'font.size': 14})

    energies = np.asarray(e)

    if Z1 is None:
        Z1_vals = range(1, 93)
    else:
        Z1_vals = [Z1]
    assert Z2 is None, "Currently only implemented for Z2=Z1"

    for energy in energies:
        sn_ratios = []
        for Z1 in Z1_vals:
            Z2 = Z1
            # NLHlin potential
            rnorm = 0.4685 / (Z1**p1 + Z2**p1)**p2
            M1_M2 = 1
            enorm = (M1_M2+1) * 14.4 * Z1 * Z2 / rnorm

            fname = f"data/p{round(p1, 2)}_{round(p2, 2)}/sn_fit_params.txt"
            Z1_, Z2_, a, b, c, d, *_ = np.loadtxt(fname, unpack=True)
            idx = np.where((Z1_ == Z1) & (Z2_ == Z2))[0]

            eps = energy / enorm
            sn = sn_fit_func(eps, a[idx], b[idx], c[idx], d[idx])
            sn *= rnorm

            # KrC potential
            if krc:
                rnorm = 0.4685 / (Z1**(1/2) + Z2**(1/2))**(2/3)
                M1_M2 = 1
                enorm = (M1_M2+1) * 14.4 * Z1 * Z2 / rnorm

                fname = f"data/p0.5_0.67/sn_fit_params_krc.txt"
                a, b, c, d, *_ = np.loadtxt(fname, unpack=True)

                eps = energy / enorm
                sn_ref = sn_fit_func(eps, a, b, c, d)
                sn_ref *= rnorm
            # ZBL potential
            else:
                rnorm = 0.4685 / (Z1**0.23 + Z2**0.23)
                M1_M2 = 1
                enorm = (M1_M2+1) * 14.4 * Z1 * Z2 / rnorm

                fname = f"data/p0.23_1/sn_fit_params_zbl.txt"
                a, b, c, d, *_ = np.loadtxt(fname, unpack=True)

                eps = energy / enorm
                sn_ref = sn_fit_func(eps, a, b, c, d)
                sn_ref *= rnorm

            sn_ratio = sn / sn_ref
            sn_ratios.append(sn_ratio)

        plt.plot(Z1_vals, sn_ratios, '.-', label=f'E={energy} eV')

    plt.axhline(1, color='k', linestyle='--', zorder=100)
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
    plt.ylim(0, 2)
    plt.xlabel('Atomic number Z$_1$=Z$_2$')
    if krc:
        plt.ylabel(r'S$_\mathrm{n}$(NLHlin)/S$_\mathrm{n}$(KrC) at E')
    else:
        plt.ylabel(r'S$_\mathrm{n}$(NLHlin)/S$_\mathrm{n}$(ZBL) at E')
    ticks = range(0, 100, 10)
    plt.xticks(ticks)
    plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
    plt.grid(True, which='major', ls='--')
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    #p1, p2 = 0.23, 1    # ZBL
    p1, p2 = 1/2, 1/2   # NLHlin
    #p1, p2 = 1/2, 2/3   # Firsov (KrC)
    Z1 = 1
    Z2 = 67
    
    #fit_sn(p1, p2, f"sn_qn_table_zbl.txt")
    #fit_sn(p1, p2, f"sn_qn_table_krc.txt")
    #fit_sn(p1, p2, f"sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt", Z1=Z1, Z2=Z2)
    #tab_sn_fit(p1, p2)

    #fit_qn(p1, p2,f"sn_qn_table_zbl.txt")
    #fit_qn(p1, p2, f"sn_qn_table_krc.txt")
    #fit_qn(p1, p2, f"sn_qn_table_nlhlin_{Z1:02d}_{Z2:02d}.txt", Z1=Z1, Z2=Z2)
    #tab_qn_fit(p1, p2)

    energies = [10, 100, 1000, 10000]
    plot_sn_ratio_at_e(energies, p1, p2)
    #plot_sn_ratio_at_e(energies, p1, p2, krc=True)
