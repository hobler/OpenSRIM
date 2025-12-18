import numpy as np
from numba import jit

proj_dtype = np.dtype([
    ("e", np.float64),
    ("pos", np.float64, (3,)),
    ("dir", np.float64, (3,)),
    ("ispec", np.int32),
    ("is_inside", np.bool_)
], align=True)

@jit(inline = 'always')
def Projectile(e, pos, dir, ispec, is_inside):
    rec = np.empty(1, dtype=proj_dtype)[0]
    rec['e'] = e
    rec['pos'] = pos    # copied
    rec['dir'] = dir    # copied
    rec['ispec'] = ispec
    rec['is_inside'] = is_inside
    return rec