"""Remove the negative values of the screening function."""
import os
import numpy as np
from scipy.interpolate import CubicSpline

for Z1 in range(1, 93):
    for Z2 in range(Z1, 93):
        fname = f'screening_zbl_{Z1}_{Z2}.dat'
        fname_save = fname + '_save'
        if not os.path.exists(fname):
            #print(fname, ' does not exist.')
            continue
        r, phi = np.loadtxt(fname, unpack=True)
        if not np.any(phi < 0):
            #print(fname, ': no negative values.')
            continue
        os.rename(fname, fname_save)
        phi_spline = CubicSpline(r, phi)
        r0 = phi_spline.roots(extrapolate=False)[0]
        mask = (r < r0)
        r = np.concatenate((r[mask], (r0,)))
        phi = np.concatenate((phi[mask], (0.,)))
        #print(fname, ':')
        with open(fname, 'w') as f:
            for i in range(len(r)):
                r[i] = round(r[i], 4)
                if round(r[i], 2) == r[i]:
                    #print(f'{r[i]:4.2f}{phi[i]:11.3e}')
                    f.write(f'{r[i]:4.2f}{phi[i]:11.3e}\n')
                else:
                    #print(f'{r[i]}{phi[i]:10.3e}')
                    f.write(f'{r[i]}{phi[i]:10.3e}\n')

