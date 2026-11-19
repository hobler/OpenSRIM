import os
from typing import Callable

import numpy as np
import scipy.constants as constants
from scipy.interpolate import PchipInterpolator

def S_u(gamma: list[float],
        e_u: list[float],
        a_u: list[float],
        d_target: float) -> float:
    return (
        np.multiply(np.multiply(gamma, e_u), np.pow(a_u, 2)) * constants.pi * d_target
    )

def Q_u(s_u: list[float],
        gamma: list[float],
        e_u: list[float],
        d_target: float) -> list[float]:
    return (
        np.multiply(s_u, np.multiply(gamma, e_u)) * d_target
    )

def S_e_SRIM(z_ion: int,
            z_target: list[int],
            d_target: float,
            s_e_f: list[float]) -> Callable[[list[float]], list[float]]:
    # TODO: define file path globally
    srim_setab_dir = './data/SRIM_setab/'
    filename = f'SRIM2013-{z_ion:02d}.dat'

    se_vals = []
    e = None
    for z in z_target:
        data = np.loadtxt(os.path.join(srim_setab_dir, filename), skiprows=6, dtype=float)
        if e is None:
            e = data[:, 0]
        # TODO: check if implementation is correct
        se_vals.append(data[:, z] * d_target)
    s_e = np.inner(np.atleast_2d(se_vals).T, np.atleast_2d(s_e_f)).T[0] * d_target
    print(s_e)
    return PchipInterpolator(e, s_e)


###############################################################################
# Universal ZBL
###############################################################################
def S_n_ZBL(epsilon: list[float], s_u: list[float]):
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

def Q_n_ZBL(epsilon: list[float], q_u: list[float]):
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
def S_n_NLH(epsilon: list[float], params: list[float], s_u: list[float]):
    a, b, c, d = params

    numerator = s_u * np.log(1 + a*epsilon)
    denominator = 2 * (epsilon + b*epsilon**c + d*epsilon**0.5)

    return (
        np.divide(numerator, denominator)
    )

def Q_n_NLH(epsilon: list[float], params: list[float], q_u: list[float]):
    a, b, c, d = params

    denominator = (4 + a*epsilon**b + c*epsilon**d)

    return (
        np.divide(q_u, denominator)
    )