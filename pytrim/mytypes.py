import numpy as np


class Projectile:
    """Data class holding projectile properties.
    
    Attributes:
        ispec (int): atom species index
        e (float): energy (eV)
        pos (NDAarray, Shape(3)): position (A)
        dir (NDAarray, Shape(3)): direction (unit vector)
    """
    def __init__(self, e, pos, dir, ispec=0):
         self.e = e
         self.pos = pos
         self.dir = dir
         self.ispec = ispec
    
    def copy(self):
        return Projectile(self.e, self.pos.copy(), self.dir.copy(), self.ispec)
