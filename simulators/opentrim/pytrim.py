"""PyTRIM aims to be a Python implementation of TRIM.

TRIM (Transport of Ions in Matter) is a widely used software package
for simulating the interaction of ions with matter, particularly for
ion implantation in semiconductors. It comes as part of SRIM, see
www.srim.org.

PyTRIM seeks to replicate the core functionalities of TRIM using
Python, making it more accessible and easier to integrate with other
Python-based tools and workflows.

Currently, the input parameters are hardcoded in this script, but future
versions may include a more user-friendly interface for specifying
simulation parameters. Also, recoils are not yet followed, and only the
mean and the straggling of the penetration depth of the primary ions are
recorded.
"""
# import os
# os.environ["NUMBA_DISABLE_JIT"] = "1"

import time
import numpy as np
import select_recoil
import scatter
import estop
import geometry
import cascade
import pytrim_stats as statistics
from mytypes import Projectile, SimParams
from numba import jit, prange

nion = 1000             # number of projectiles to simulate

zmin = 0.0              # minimum z coordinate of the target (A)
zmax = 4000.0           # maximum z coordinate of the target (A)
z1 = 5                  # atomic number of projectile
m1 = 11.009             # mass of projectile (amu)
z2 = 14                 # atomic number of target
m2 = 28.086             # mass of target atom (amu)
density = 0.04994       # target density (atoms/A^3)
corr_lindhard1 = 1.5    # Correction factor to Lindhard stopping power (B->Si)
corr_lindhard2 = 1.0    # Correction factor to Lindhard stopping power (Si->Si)

# Setup modules
select_recoil.setup(density)
scatter.setup(z1, m1, z2, m2)
estop.setup(corr_lindhard1, z1, m1, corr_lindhard1, z2, m2, density)
geometry.setup(zmin, zmax)
cascade.setup()
sim_params = SimParams( nspec = 2, 
                        nbin = 40, 
                        limits = (0.0, 4000.0))
statistics.setup(nspec=sim_params.nspec, nbin=sim_params.nbin, limits=sim_params.limits)

@jit(fastmath=True, cache=False, parallel=True, nogil=True)
def simulate(nion, sim_params_tup, follow_recoils=False):
    # Initial conditions of the projectile
    proj_init = Projectile(
        50000.0,                         # energy (eV)
        np.array([0.0, 0.0, 0.0]),     # position (A)
        np.array([0.0, 0.0, 1.0]),     # direction (unit vector)
        0,
        True
    )
    proj_dummy = np.full(1, proj_init)
    # Linked list with pointers to dummy array
    # Pointers will be overwritten during simulation
    proj_sim = [proj_dummy for _ in range(nion)]

    # Simulate the trajectories
    for i in prange(nion):
        proj_sim[i] = cascade.trajectory(proj_dummy[0], follow_recoils)
    
    # TODO alternatives???
    proj_count = 0
    sim_params = SimParams(*sim_params_tup)
    hist = statistics.Histogram_1d(sim_params.nspec, sim_params.nbin, sim_params.limits)
    mom = statistics.Moment_1d(sim_params.nspec, 4)
    for proj_lst in proj_sim:
        proj_count += proj_lst.size
        for proj in proj_lst:
            if proj.is_inside:
                hist.score(proj.ispec, proj.pos[2])
                mom.score(proj.ispec, proj.pos[2])
    return proj_count, hist.counts, mom._mom

if __name__ == "__main__":
    times = []
    proj_counts = []
    counts = [1000, 10000]
    simulate(10, sim_params.to_tuple(), follow_recoils=True)
    for c in counts:
        start_time = time.time()
        proj_count, hist_buf, mom_buf = simulate(c, sim_params.to_tuple(), follow_recoils=True)
        times.append(time.time() - start_time)
        proj_counts.append(proj_count)
        
        statistics.hist.counts = hist_buf
        statistics.mom._mom = mom_buf
    print(times, proj_counts)
    
    # Output the results
    statistics.print_results()
    statistics.plot_results(log=True)
