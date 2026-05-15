"""Calculate an integral using Gauss-Legendre quadrature.
"""
import numpy as np
from scipy.special import roots_legendre


def quad_legendre(func, a, b, n):
    """Calculate an integral using Gauss-Legendre quadrature.

    Parameters:
        func (callable): Function to integrate. Must take a single argument.
        a (float): Lower limit of integration.
        b (float): Upper limit of integration.
        n (int): Number of quadrature points.

    Returns:
        float: Approximation of the integral of func from a to b.
    """
    # Get the roots and weights for the standard interval [-1, 1]
    roots, weights = roots_legendre(n)
    # Transform roots to the interval [a, b]
    transformed_roots = 0.5 * (b - a) * roots + 0.5 * (b + a)
    # Evaluate the function at the transformed roots
    func_values = func(transformed_roots)
    # Compute the integral approximation
    integral = 0.5 * (b - a) * np.dot(weights, func_values)

    return integral
