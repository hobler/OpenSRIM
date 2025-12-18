import numpy as np
from numba.experimental import jitclass
from numba import float64, int32, bool

@jitclass(spec = [
    ("e", float64),
    ("pos", float64[:]),
    ("dir", float64[:]),
    ("ispec", int32),
    ("is_inside", bool)
])
class Projectile:
    """Data class holding projectile properties.
    
    Attributes:
        e (float): energy (eV)
        pos (ndarray): position (A, size 3)
        dir (ndarray): direction (unit vector, size 3)
        ispec (int): atom species index
        is_inside (bool): whether the projectile is inside the target"""
    # TODO find out how to include default values (overload? None check?)
    def __init__(self, e, pos, dir, ispec, is_inside):
        self.e = e
        self.pos = pos
        self.dir = dir
        self.ispec = ispec
        self.is_inside = is_inside
    
    def copy(self):
        return Projectile(
            self.e, 
            self.pos.copy(), 
            self.dir.copy(), 
            self.ispec,
            self.is_inside
        )
