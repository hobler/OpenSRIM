from abc import ABC, abstractmethod
import numpy as np

class DifferentialEquation(ABC):
    def __init__(self):
        super().__init__()

    @abstractmethod
    def equation(self, x, f_x):
        pass

class Bowyer_17(DifferentialEquation):
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
    
    def equation(self, x, y):
        E = x
        y0 = y[0]

        a = self.a
        b = self.b
        c = self.c

        df_dx = np.divide((c(E) - b(E)*y0), (a(E)))
        return df_dx

class Bowyer_18(DifferentialEquation):
    def __init__(self, a, b, c, d, e, f):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
    
    def equation(self, x, y):
        E = x
        y0 = y[0]

        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e
        f = self.f

        df_dx = np.divide(c(E) + d(E) + e(E) + f(E) - b(E)*y0, a(E))
        return df_dx