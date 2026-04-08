from collections.abc import Callable

import numpy as np
from scipy import integrate, interpolate, differentiate

from koral_settings import KORALSettings
import differential_equations

def solve_koral(
        e: list[float],
        s_t: Callable[[list[float]], list[float]],
        mu_s_n: Callable[[list[float]], list[float]],
        q_t: Callable[[list[float]], list[float]],
        l: Callable[[list[float]], list[float]],
        w: Callable[[list[float]], list[float]],
        k: Callable[[list[float]], list[float]],
        settings: KORALSettings
) -> list[list[float]]:
    pass

def solve_rp(
        E: list[float],
        s_t: Callable[[list[float]], list[float]],
        mu_s_n: Callable[[list[float]], list[float]],
        q_t: Callable[[list[float]], list[float]],
        l: Callable[[list[float]], list[float]],
        w: Callable[[list[float]], list[float]],
        k: Callable[[list[float]], list[float]],
        settings: KORALSettings
) -> list[list[float]]:
    E_0 = E[0]

    def a(e: list[float]) -> list[float]:
        return s_t(e)
    
    def b(e: list[float]) -> list[float]:
        return (
            np.divide(mu_s_n(e), 2*e)
        )
    
    def c(e: list[float]) -> list[float]:
        return np.ones_like(e)
    
    def d(e: list[float]) -> list[float]:
        return (
            np.divide(q_t(e), 2)
        )
    
    def e(e: list[float]) -> list[float]:
        return (
            -1*(s_t(e) + l(e))
        )
    
    def f(e: list[float]) -> list[float]:
        return w(e)
    
    init_0 = (
        np.divide(2*E_0, s_t(E_0) + mu_s_n(E_0))
    )
    def init_1(init_e: float) -> float:
        # TODO: Implement
        # (1-Q_t(E_start)*y0(E_start)/(8*E_start^2)-
        # (S_t(E_start)+L(E_start))*y0(E_start)/(2*E_start)+
        # W(E_start)*y0(E_start))*2*E_start/(S_t(E_start)+materials.mu*S_n(E_start))

        # return 0.0 as long as Q_t = Q_n
        return 0.0

    y = solve_iterative_refinement(
        E, a, b, c, d, e, f, init_0, init_1, settings
    )
    
    return y

def solve_rc(
        E: list[float],
        s_t: Callable[[list[float]], list[float]],
        mu_s_n: Callable[[list[float]], list[float]],
        q_t: Callable[[list[float]], list[float]],
        l: Callable[[list[float]], list[float]],
        w: Callable[[list[float]], list[float]],
        k: Callable[[list[float]], list[float]],
        r_p: Callable[[list[float]], list[float]],
        settings: KORALSettings
) -> list[list[float]]:
    E_0 = E[0]

    def a(e: list[float]) -> list[float]:
        return s_t(e)
    
    def b(e: list[float]) -> list[float]:
        return (
            np.zeros_like(e)
        )
    
    def c(e: list[float]) -> list[float]:
        return (
            np.multiply(r_p(e), 2)
        )
    
    def d(e: list[float]) -> list[float]:
        return (
            np.divide(q_t(e), 2)
        )
    
    def e(e: list[float]) -> list[float]:
        return (
            np.multiply(-1,s_t(e))
        )
    
    def f(e: list[float]) -> list[float]:
        return np.zeros_like(e)
    
    init_0 = (
        np.divide(2*r_p(E_0)*E_0, s_t(E_0))
    )
    def init_1(i_0: float) -> float:
        return 0.0

    y = solve_iterative_refinement(
        E, a, b, c, d, e, f, init_0, init_1, settings
    )
    
    return y

def solve_rr(
        E: list[float],
        s_t: Callable[[list[float]], list[float]],
        mu_s_n: Callable[[list[float]], list[float]],
        q_t: Callable[[list[float]], list[float]],
        l: Callable[[list[float]], list[float]],
        w: Callable[[list[float]], list[float]],
        k: Callable[[list[float]], list[float]],
        r_p: Callable[[list[float]], list[float]],
        settings: KORALSettings
) -> list[list[float]]:
    E_0 = E[0]

    def a(e: list[float]) -> list[float]:
        return s_t(e)
    
    def b(e: list[float]) -> list[float]:
        return (
            np.divide(3*mu_s_n(e), 2*e)
        )
    
    def c(e: list[float]) -> list[float]:
        return (
            np.multiply(2, r_p(e))
        )
    
    def d(e: list[float]) -> list[float]:
        return (
            np.divide(q_t(e), 2)
        )
    
    def e(e: list[float]) -> list[float]:
        return (
            -1*(s_t(e) + 3*l(e))
        )
    
    def f(e: list[float]) -> list[float]:
        return (
            3*(w(e) + 3*k(e))
        )
    
    init_0 = (
        np.divide(2*r_p(E_0)*E_0, s_t(E_0)+3/2*mu_s_n(E_0))
    )
    def init_1(i_0: float) -> float:
        return (
            E_0 * (2*r_p(E_0) - (s_t(E_0)+3*l(E_0)) * i_0/E_0 + 3*(w(E_0)+3*k(E_0))*i_0) /
            (s_t(E_0)+3/2*mu_s_n(E_0)) 
        )

    y = solve_iterative_refinement(
        E, a, b, c, d, e, f, init_0, init_1, settings
    )
    
    return y

def solve_iterative_refinement(
        x: list[float],
        a: Callable[[list[float]], list[float]],
        b: Callable[[list[float]], list[float]],
        c: Callable[[list[float]], list[float]],
        d: Callable[[list[float]], list[float]],
        e: Callable[[list[float]], list[float]],
        f: Callable[[list[float]], list[float]],
        init_0 : float,
        init_1 : Callable[[list[float]], list[float]],
        settings: KORALSettings
) -> list[float]:
    ivp_settings = {
        'method': settings.integration_method,
        'rtol': settings.rtol,
        'atol': settings.atol
        }
    x_0 = x[0]

    df_0 = differential_equations.Bowyer_17(a, b, c)
    y_0 = solve_initial_value_problem(df_0, x, [init_0], ivp_settings)
    dy_0 = df_0.equation(x, [y_0])

    for i in range(settings.nr_interations):
        f_0_interpolated = interpolate.PchipInterpolator(x, y_0)
        df_0_interpolated = interpolate.PchipInterpolator(x, dy_0)

        # Precalculate to enhance calculation speed
        ddy_0 = differentiate.derivative(df_0_interpolated, x).df
        ddf_0_interpolated = interpolate.PchipInterpolator(x, ddy_0)

        def D(x: list[float]):
            return (
                np.multiply(d(x), ddf_0_interpolated(x))
            )
        
        def E(x: list[float]):
            return (
                np.multiply(e(x), df_0_interpolated(x))
            )
        
        def F(x: list[float]):
            return (
                np.multiply(f(x), f_0_interpolated(x))
            )
        
        df_1 = differential_equations.Bowyer_18(a, b, c, D, E, F)
        y_1 = solve_initial_value_problem(df_1, x, [init_1(f_0_interpolated(x_0))], ivp_settings)
        dy_1 = df_1.equation(x, [y_1])

        y_0 = y_0 + y_1
        dy_0 = dy_0 + dy_1

    return y_0

def solve_initial_value_problem(
        differential_equation: differential_equations.DifferentialEquation,
        range: list[float],
        f_0: float,
        settings # TODO: Welchen Datentyp verwenden die Settings?
) -> list[float]:
    result = integrate.solve_ivp(
            differential_equation.equation,
            (np.min(range), np.max(range)),
            f_0,
            t_eval=range,
            **settings)
    
    return result.y[0]