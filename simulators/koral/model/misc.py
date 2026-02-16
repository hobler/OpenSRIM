import scipy.constants as constants
import numpy as np

def g_over_cm3_to_1_over_A3(g_over_cm3: float, m_element: float) -> float:
    """Converts element densitiy from g/cm^3 to 1/A^3

    Args:
        g_over_cm3 (float): Element density in g/cm^3
        m_element (float): Relative mass of element

    Returns:
        float: Element density in 1/A^3
    """    
    return(
        g_over_cm3 / m_element * constants.Avogadro / 1e24
    )

def mu(m_ion: int, m_target: list[int]) -> list[float]:
    return(
        np.divide(m_target,m_ion)
    )

def gamma(mu: list[float]) -> list[float]:
    numerator = 4*mu
    denominator = ((1+mu)**2)
    return(
        np.divide(numerator, denominator)
    )

def a_ZBL(z_ion: int, z_target: list[int]) -> list[float]:
    denominator = np.pow(z_ion, 0.23) + np.pow(z_target, 0.23)
    return(
        np.divide(0.46850, denominator)
    )

def a_F(z_ion: int, z_target: list[int]) -> list[float]:
    denominator = np.pow(np.pow(z_ion, 0.5) + np.pow(z_target, 0.5), 2/3)
    return(
        np.divide(0.46850, denominator)
    )

def a_L(z_ion: int, z_target: list[int]) -> list[float]:
    denominator = np.pow(np.pow(z_ion, 2/3) + np.pow(z_target, 2/3), 0.5)
    return(
        np.divide(0.46850, denominator)
    )

def a_NLH(z_ion: int, z_target: list[int]) -> list[float]:
    denominator = np.pow(np.pow(z_ion, 0.5) + np.pow(z_target, 0.5), 0.5)
    return(
        np.divide(0.46850, denominator)
    )

def E_u(mu: list[float], z_ion: int, z_target: list[int], a_u: list[float]) -> list[float]:
    numerator = np.multiply((1+mu), z_target) * z_ion * constants.elementary_charge
    denominator = np.multiply(mu, a_u) * 4 * constants.pi * constants.epsilon_0 * 1e-10
    return (
        np.divide(numerator, denominator)
    )

def epsilon(e: list[float], e_u: list[float]) -> list[float]:
    """Calculates $\\varepsilon(E)$

    $$
    \\varepsilon(E) = \\frac{E}{E_\\text{U}}
    $$

    Returns:
        float: result of $\\varepsilon(E)$
    """

    return (
        np.divide(e, e_u)
    )