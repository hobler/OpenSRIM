"""Table class that allows 1D linear or cubic interpolation."""
# import os
# os.environ["NUMBA_DISABLE_JIT"] = "1"
import math
from typing import Optional
from numba import jit
from numba.core.types import float64
from numba.experimental import jitclass
import numpy as np

@jitclass
class Table1D:
    """1D lookup table with linear or cubic interpolation.

    Cubic interpolation is done if derivatives at the table points are 
    provided during initialization.

    Attributes:
        linear (bool): whether to use linear interpolation.
        regular (bool): whether the x values are regularly spaced.
        powerof2 (bool): whether the x values are a constant times a 
            power of 2.
        x (np.ndarray): x values of the table (length N).
        y (np.ndarray): y values of the table (length N).
        a, b, c, d (np.ndarray): coefficients used to interpolate y for a 
            given value of x (length N - 1):
                
                y = a + b * t + c * t^2 + d * t^3

            where t = (x - x0) / (x1 - x0) and [x0, x1] is the interval 
            containing x.
            For linear interpolation c and d are not defined.
    """
    
    regular: bool
    powerof2: bool
    x: float64[:]   # ty:ignore[invalid-type-form]
    y: float64[:]   # ty:ignore[invalid-type-form]
    a: float64[:]   # ty:ignore[invalid-type-form]
    b: float64[:]   # ty:ignore[invalid-type-form]
    c: Optional[float64[:]] # ty:ignore[invalid-type-form]
    d: Optional[float64[:]] # ty:ignore[invalid-type-form]
    def __init__(self, x, y, dydx=None, regular=False, powerof2=False):
        self.regular = regular
        self.powerof2 = powerof2

        if regular:
            dx = np.diff(x)
            if not np.allclose(dx, dx[0]):
                raise ValueError("x values are not regularly spaced.")
        if powerof2:
            if not np.allclose(x[1:] / x[:-1], 2.0):
                raise ValueError("x values are not a constant times a "
                                 "power of 2.")

        self.x = x
        self.y = y
        self.a = y[:-1]
        self.c = None
        self.d = None
        if dydx is None:
            self.b = np.diff(y)
        else:
            dx = np.diff(x)
            self.b = dydx[:-1] * dx
            self.c = 3 * np.diff(y) - 2 * dydx[:-1] * dx - dydx[1:] * dx
            self.d = -2 * np.diff(y) + dydx[:-1] * dx + dydx[1:] * dx

    # TODO inline?
    def get_index(self, x_val): 
        """Get the index of the table for a given x value.

        Parameters:
            x_val (float): x value to find index for.
        
        Returns:
            (int): Index of the table.
        """
        if self.regular:
            idx = math.floor((x_val - self.x[0]) / (self.x[1] - self.x[0]))
        elif self.powerof2:
            idx = math.floor(math.log2(x_val / self.x[0]))
        else:
            idx = np.searchsorted(self.x, x_val) - 1

        return max(0, min(len(self.a) - 1, idx))

    def interpolate(self, x_val, extrapolate=True):
        """Interpolate the table at a given x value.

        Parameters:
            x_val (float): x value to interpolate.
            extrapolate (bool): whether to allow extrapolation outside 
                the table range.

        Returns:
            (tuple[bool, float]): Success, Interpolated y value.
        """
        if not extrapolate and (x_val <= self.x[0] or x_val >= self.x[-1]):
            return False, 0.0
        idx = self.get_index(x_val)

        t = (x_val - self.x[idx]) / (self.x[idx + 1] - self.x[idx])
        if self.c is not None and self.d is not None:
            return True, (self.a[idx] + t * (self.b[idx] + 
                    t * (self.c[idx] + t * self.d[idx])))
        else:
            return True, self.a[idx] + self.b[idx] * t

@jit
def test():
    def f(x):
        return x**3 - 2*x + 1
    def df(x):
        return 3*x**2 - 2
    # Test of the Table1D class
    x = np.array([-1.5, 0.0, 2.0, 5.0])
    y = f(x)
    dydx = df(x)

    table_cubic = Table1D(x, y, dydx, False, False)
    x_vals = np.linspace(-2.0, 7.0, 101)
    for x_val in x_vals:
        y_true = f(x_val)
        succ, y_interp = table_cubic.interpolate(x_val, True)
        assert np.isclose(y_true, y_interp)
    print("Cubic interpolation test passed.")

    for x_val in x_vals:
        y_true = f(x_val)
        succ, y_interp = table_cubic.interpolate(x_val, False)
        if not succ and x_val < x[0] or x_val > x[-1]:
            print("Caught expected out-of-bounds error for x =", x_val)
        else:
            assert np.isclose(y_true, y_interp)
    print("Cubic interpolation with bounds test passed.")

if __name__ == "__main__":
    test()