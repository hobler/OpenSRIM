import os
from collections.abc import Callable

import numpy as np
import scipy.constants as constants
from scipy.interpolate import PchipInterpolator

def S_u(gamma: float,
        e_u: float,
        a_u: float,
        d_target: float) -> float:
    return (
        np.multiply(np.multiply(gamma, e_u), np.pow(a_u, 2)) * constants.pi * d_target
    )

def Q_u(s_u: float,
        gamma: float,
        e_u: float) -> float:
    return (
        np.multiply(s_u, np.multiply(gamma, e_u))
    )

def S_e_SRIM(z_ion: int,
            z_target: int,
            d_target: float,
            f_target: float) -> Callable[[list[float]], list[float]]:
    # TODO: define file path globally
    srim_setab_dir = './data/SRIM_setab/'
    filename = f'SRIM2013-{z_ion:02d}.dat'

    data = np.loadtxt(os.path.join(srim_setab_dir, filename), skiprows=6, dtype=float)
    e = data[:,0]
    s_e = np.multiply(np.multiply(data[:, z_target], d_target), f_target)
    return PchipInterpolator(e, s_e)


###############################################################################
# Universal ZBL
###############################################################################
def S_n_ZBL(epsilon: list[float], s_u: float) -> list[float]:
    """Calculates the universal nuclear stopping power $S_n(\epsilon)$

    $$
    S_\\text{n}(E) = 
    \\frac{\\text{ln}(1+1.1383\\varepsilon)}
    {2(\\varepsilon+0.01321\\varepsilon^{0.21226}+0.19593\\varepsilon^{0.5})}
    S_\\text{U}
    $$
    """
    # TODO: Check if changes for multi element target are necessary
    numerator = s_u * np.log(1+1.1383*epsilon)
    denominator = (2*(epsilon + 0.01321*epsilon**0.21226 + 0.19593*epsilon**0.5))

    return (
        np.divide(numerator, denominator)
    )

def Q_n_ZBL(epsilon: list[float], q_u: float) -> list[float]:
    """Calculates the nuclear energy loss $Q_n(E)$

    $$
    Q_n(E) = \\frac{1}
    {4+0.197\\varepsilon^{-1.6991}+6.584\\varepsilon^{-1.0494}}Q_\\text{U}
    $$
    """

    denominator = (4+0.197*epsilon**-1.6991 + 6.584*epsilon**-1.0494)

    return (
        np.divide(q_u, denominator)
    )

###############################################################################
# NLH
###############################################################################
def S_n_NLH(epsilon: list[float], params: list[float], s_u: float) -> list[float]:
    a, b, c, d = params

    numerator = s_u * np.log(1 + a*epsilon)
    denominator = 2 * (epsilon + b*epsilon**c + d*epsilon**0.5)

    return (
        np.divide(numerator, denominator)
    )

def Q_n_NLH(epsilon: list[float], params: list[float], q_u: float) -> list[float]:
    a, b, c, d = params

    denominator = (4 + a*epsilon**b + c*epsilon**d)

    return (
        np.divide(q_u, denominator)
    )