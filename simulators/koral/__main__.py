from collections.abc import Callable
import time
import tomllib
import os

import numpy as np
from scipy import interpolate

from koral_settings import KORALSettings
import misc, stopping_powers, solver

KORAL_VERSION = 'alpha'

# TODO: Abklären wie diese funktion aus __main__.py mit Parametern aufgerufen werden kann
def KORAL(working_directory_path: str) -> int:
    dir_path = os.path.normpath(working_directory_path)

    files = os.listdir(dir_path)
    for file in files:
        if file.split('.')[-1] == 'toml':
            toml_file = os.path.join(dir_path, file)
            with open(toml_file, "rb") as f:
                toml = tomllib.load(f)

    settings = KORALSettings(
        nr_iterations=toml['settings']['nr_iterations'],
        integration_method=toml['settings']['integration_method'],
        rtol=toml['settings']['rtol'],
        atol=toml['settings']['atol']
    )
    
    f_targets = np.divide(toml['params']['c_target'], np.sum(toml['params']['c_target']))

    # start_time = time.time()
    E = np.geomspace(toml['params']['start_energy'],
                     toml['params']['stop_energy'],
                     toml['params']['nr_values'])
    
    s_e_targets:    list[Callable[[list[float]], list[float]]] = []
    s_n_targets:    list[Callable[[list[float]], list[float]]] = []
    mu_targets:     list[Callable[[list[float]], list[float]]] = []
    s_t_targets:    list[Callable[[list[float]], list[float]]] = []
    q_n_targets:    list[Callable[[list[float]], list[float]]] = []
    q_t_targets:    list[Callable[[list[float]], list[float]]] = []
    w_targets:      list[Callable[[list[float]], list[float]]] = []
    k_targets:      list[Callable[[list[float]], list[float]]] = []
    l_targets:      list[Callable[[list[float]], list[float]]] = []

    for i, z_target in enumerate(toml['params']['z_target']):
        m_target = toml['params']['m_target'][i]
        d_target = toml['params']['d_target'][i]
        f_target = f_targets[i]

        mu_target = misc.mu(toml['params']['m_ion'], m_target)
        mu_targets.append(mu_target)
        gamma_target = misc.gamma(mu_target)

        s_e_target = stopping_powers.S_e_SRIM(toml['params']['z_ion'],
                                        z_target,
                                        d_target,
                                        f_target)
        s_e_targets.append(s_e_target)

        if toml['params']['method'] == 'ZBL':
            a_u_target = misc.a_ZBL(toml['params']['z_ion'], z_target)
        elif toml['params']['method'] == 'NLH':
            a_u_target = misc.a_NLH(toml['params']['z_ion'], z_target)
            sn_nlh_params = misc.get_sn_nlh_params(toml['params']['z_ion'], z_target)
            qn_nlh_params = misc.get_qn_nlh_params(toml['params']['z_ion'], z_target)
        else:
            return -1 # TODO: Handle unknown method

        e_u_target = misc.E_u(mu_target, toml['params']['z_ion'], z_target, a_u_target)
        s_u_target = stopping_powers.S_u(gamma_target, e_u_target, a_u_target, d_target)
        q_u_target = stopping_powers.Q_u(s_u_target, gamma_target, e_u_target)

        if toml['params']['method'] == 'ZBL':
            def s_n_target(e: list[float]) -> list[float]:
                epsilon = misc.epsilon(e, e_u_target)
                return(
                    stopping_powers.S_n_ZBL(epsilon, s_u_target) * f_target
                )
            s_n_targets.append(s_n_target)

            def q_n_target(e: list[float]) -> list[float]:
                epsilon = misc.epsilon(e, e_u_target)
                return(
                    stopping_powers.Q_n_ZBL(epsilon, q_u_target) * f_target
                )
            q_n_targets.append(q_n_target)
            
        elif toml['params']['method'] == 'NLH':
            def s_n_target(e: list[float]) -> list[float]:
                epsilon = misc.epsilon(e, e_u_target)
                return(
                    stopping_powers.S_n_NLH(epsilon, sn_nlh_params, s_u_target) * f_target
                )
            s_n_targets.append(s_n_target)

            def q_n_target(e: list[float]) -> list[float]:
                epsilon = misc.epsilon(e, e_u_target)
                return(
                    stopping_powers.Q_n_NLH(epsilon, qn_nlh_params, q_u_target) * f_target
                )
            q_n_targets.append(q_n_target)

        else:
            return -1 # TODO: Handle unknown method

        def q_t_target(e: list[float]) -> list[float]:
            # TODO: Sobald Q_e und Q_ne implementiert wurde
            # muss q_t(e) = q_n(e) + 2*q_ne(e) + q_e(e)
            # implementiert werden
            return q_n_target(e)
        q_t_targets.append(q_t_target)
    
        def s_t_target(e: list[float]) -> list[float]:
            return (
                s_n_target(e) + s_e_target(e)
            )
        s_t_targets.append(s_t_target)
        
        def w_target(e: list[float]) -> list[float]:
            return (
                np.divide((1-2*mu_target)*q_n_target(e), 8*np.pow(e,2)) -
                np.divide(mu_target * s_n_target(e), 2*e)
            )
        w_targets.append(w_target)
        
        def k_target(e: list[float]) -> list[float]:
            return (
                np.divide(np.pow(mu_target,2)*q_n_target(e), 24*np.pow(e,2))
            )
        k_targets.append(k_target)
        
        def l_target(e: list[float]) -> list[float]:
            return (
                np.divide(-1*mu_target*q_n_target(e), 2*e)
            )
        l_targets.append(l_target)

    def s_e(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in s_e_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def s_n(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in s_n_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def mu_s_n(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for i, f in enumerate(s_n_targets):
            y = np.multiply(f(e), mu_targets[i])
            sum = sum + y
        return sum
    
    def s_t(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in s_t_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def q_n(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in q_n_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def q_t(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in q_t_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def w(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in w_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def k(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in k_targets:
            y = f(e)
            sum = sum + y
        return sum
    
    def l(e: list[float]) -> list[float]:
        sum = np.zeros_like(e)
        for f in l_targets:
            y = f(e)
            sum = sum + y
        return sum

    r_p = solver.solve_rp(
        E, s_t, mu_s_n, q_t, l, w, k, settings
    )
    r_p_interpolated = interpolate.PchipInterpolator(E, r_p)

    r_c = solver.solve_rc(
        E, s_t, mu_s_n, q_t, l, w, k, r_p_interpolated, settings
    )

    r_r = solver.solve_rr(
        E, s_t, mu_s_n, q_t, l, w, k, r_p_interpolated, settings
    )

    sigma_x = np.sqrt((r_c + 2*r_r)/3 - np.pow(r_p, 2))
    sigma_z = np.sqrt((r_c - r_r)/3)

    # end_time = time.time()
    # print(f'Computation time: {end_time-start_time:.3f} s')

    np.savetxt(
        os.path.join(dir_path, 'koral.csv'),
        header=f'KORAL Version: {KORAL_VERSION}\nE/eV; R_p/angstrom; sigma_x/angstrom; sigma_z/angstrom; S_e/eV/angstrom; S_n/eV/angstrom; Q_n/eV^2/angstrom',
        X=np.array([E, r_p, sigma_x, sigma_z, s_e(E), s_n(E), q_n(E)]).T,
        delimiter=';'  
    )
    return 1

# if __name__ == '__main__':
#     from matplotlib import pyplot as plt

#     result = KORAL('./simulators/koral/kt')

#     # plt.loglog(result[0,:], result[1,:], label=r'$S_e$')
#     # plt.loglog(result[0,:], result[2,:], label=r'$S_n$')
#     # plt.loglog(result[0,:], result[3,:], label=r'$Q_n$')
#     plt.loglog(result[0,:], result[4,:], label=r'$\overline{R_p}$')
#     # plt.loglog(result[0,:], result[5,:], label=r'$\overline{R_c^2}$')
#     # plt.loglog(result[0,:], result[6,:], label=r'$\overline{R_r^2}$')
#     plt.loglog(result[0,:], result[7,:], label=r'$\sigma_x$')
#     plt.loglog(result[0,:], result[8,:], label=r'$\sigma_z$')

#     plt.legend(loc="upper left", ncol=2)
#     plt.xlabel('Energy/eV')
#     plt.ylabel(r'Stopping power/eV$/\AA$, Straggling/eV$^2/\AA$, Range/$\AA$')
#     plt.grid(True)
#     plt.show()