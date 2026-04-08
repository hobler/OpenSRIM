import numpy as np

from koral_input import KORALInput
from koral_settings import KORALSettings
from base.validator import valid
from model import misc, stopping_powers

def KORAL(input_params: KORALInput, settings: KORALSettings) -> list[list[float]]:
    if not valid(input_params, settings):
        print('Invalid function parameters!')
        return 0
    
    E = np.geomspace(input_params.start_energy,
                     input_params.stop_energy,
                     input_params.nr_values)

    s_e = stopping_powers.S_e_SRIM(input_params.z_ion,
             input_params.z_target,
             input_params.d_target,
             input_params.s_e_f)
    
    # TODO: wie berechnet sich mu für mehrere ziel-elemente?
    mu = misc.mu(input_params.m_ion, input_params.m_target)
    gamma = misc.gamma(mu)
    
    if input_params.method == 'ZBL':
        a_u = misc.a_ZBL(input_params.z_ion, input_params.z_target)
        E_u = misc.E_u(mu, input_params.z_ion, input_params.z_target, a_u)
        s_u = stopping_powers.S_u(gamma, E_u, a_u, input_params.d_target)
        q_u = stopping_powers.Q_u(s_u, gamma, E_u, input_params.d_target)
        
        def s_n(e: list[float]) -> list[float]:
            epsilon = misc.epsilon(e, E_u)
            return(stopping_powers.S_n_ZBL(epsilon, s_u))

        def q_n(e: list[float]) -> list[float]:
            epsilon = misc.epsilon(e, E_u)
            return(stopping_powers.Q_n_ZBL(epsilon, q_u))
        
    def s_t(e: list[float]) -> list[float]:
        return (
            s_n(e) + s_e(e)
        )
    
    def w(e: list[float]) -> list[float]:
        return (
            (1-2*mu)*q_n(e)
        )

    S_e = s_e(E)
    S_n = s_n(E)
    Q_n = q_n(E)
    S_tot = np.maximum(S_e + S_n, 1e-30)

    # CSDA projected range Rp(E) = ∫₀ᴱ dE' / S_tot(E')
    inv_s = 1.0 / S_tot
    dE = np.diff(E)
    avg_inv_s = (inv_s[:-1] + inv_s[1:]) * 0.5
    Rp = np.concatenate(([0.0], np.cumsum(dE * avg_inv_s)))

    # Longitudinal straggling σ_x²(E) = ∫₀ᴱ Q_n / S_tot² dE'
    integrand_x = Q_n / S_tot**2
    avg_x = (integrand_x[:-1] + integrand_x[1:]) * 0.5
    sigma_x2 = np.concatenate(([0.0], np.cumsum(dE * avg_x)))
    sigma_x = np.sqrt(np.maximum(sigma_x2, 0.0))

    # Lateral straggling σ_z²(E) = ∫₀ᴱ Q_n · Rp / S_tot dE'
    integrand_z = Q_n * Rp / S_tot
    avg_z = (integrand_z[:-1] + integrand_z[1:]) * 0.5
    sigma_z2 = np.concatenate(([0.0], np.cumsum(dE * avg_z)))
    sigma_z = np.sqrt(np.maximum(sigma_z2, 0.0))

    # result rows: E, S_e, S_n, Q_n, Rp, sigma_x, sigma_z
    return np.array([E, S_e, S_n, Q_n, Rp, sigma_x, sigma_z])

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    input_params = KORALInput(
        method='NLH',
        z_ion=33,
        m_ion=74.992,
        z_target=[14],
        m_target=[28.085],
        d_target=0.04996,
        s_e_f=[1],
        start_energy=1,
        stop_energy=10e6,
        nr_values=100)

    settings = KORALSettings()

    result = KORAL(input_params, settings)
    E       = result[0, :]
    S_e     = result[1, :]
    S_n     = result[2, :]
    Q_n     = result[3, :]
    Rp      = result[4, :]
    sigma_x = result[5, :]
    sigma_z = result[6, :]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.loglog(E, S_e,  label='S_e')
    ax1.loglog(E, S_n,  label='S_n')
    ax1.loglog(E, Q_n,  label='Q_n')
    ax1.set_xlabel('Energy (eV)')
    ax1.set_ylabel('Stopping (eV/Å)')
    ax1.legend()
    ax1.grid(True)
    ax1.set_title('Stopping Powers')

    ax2.loglog(E, Rp,      label='Rp (proj. range)')
    ax2.loglog(E, sigma_x, label='σ_x (long. straggling)')
    ax2.loglog(E, sigma_z, label='σ_z (lat. straggling)')
    ax2.set_xlabel('Energy (eV)')
    ax2.set_ylabel('Distance (Å)')
    ax2.legend()
    ax2.grid(True)
    ax2.set_title('Range & Straggling')

    plt.tight_layout()
    plt.show()