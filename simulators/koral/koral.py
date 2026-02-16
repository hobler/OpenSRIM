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

    return np.array([E,s_e(E),s_n(E),q_n(E)])

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    input_params = KORALInput(
        method='ZBL',
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
    plt.loglog(result[0,:], result[1,:])
    plt.loglog(result[0,:], result[2,:])
    plt.loglog(result[0,:], result[3,:])
    plt.grid(True)
    plt.show()