import numpy as np


class Projectile:
    """Data class holding projectile properties.
    
    Attributes:
        e (float): energy (eV)
        pos (ndarray): position (A, size 3)
        dir (ndarray): direction (unit vector, size 3)
        ispec (int): atom species index
        is_inside (bool): whether the projectile is inside the target"""
    def __init__(self, e, pos, dir, ispec=0, is_inside=True):
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
