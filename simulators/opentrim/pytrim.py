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
import time
import os
import config
import numpy as np
import select_recoil
import scatter
import cm_scatter
import estop
import geometry
import cascade
import pytrim_stats as statistics
from mytypes import SimParams
from nlhlin import read_coefs
from simulator import simulate, simulate_chunked, simulate_adaptive

zmin = 0.0              # minimum z coordinate of the target (A)
zmax = 4000.0           # maximum z coordinate of the target (A)
pot_model = 'ZBL_magic'  # potential model for scattering
z1 = 5                  # atomic number of projectile
m1 = 11.009             # mass of projectile (amu)
z2 = 14                 # atomic number of target
m2 = 28.086             # mass of target atom (amu)
density = 0.04994       # target density (atoms/A^3)
corr_lindhard1 = 1.5    # Correction factor to Lindhard stopping power (B->Si)
corr_lindhard2 = 1.0    # Correction factor to Lindhard stopping power (Si->Si)

start = time.time()
# Setup modules
recoil_params_tup = select_recoil.setup(density)
scatter_params_tup = scatter.setup(z1, m1, z2, m2, pot_model)
cm_scatter.setup(n_absc=4)
estop_params_tup = estop.setup(corr_lindhard1, z1, m1, corr_lindhard2, z2, m2, density)
geometry_params_tup = geometry.setup(zmin, zmax)
cascade_params_tup = cascade.setup()
sim_params = SimParams( rng_seed = np.random.randint(2**31, dtype=np.uint32),
                        # Seed can be specified manually or generated automatically (default)
                        stat_params_tup = (2, 40, (0.0, 4000.0)),
                        cascade_params_tup = cascade_params_tup,
                        recoil_params_tup = recoil_params_tup,
                        geometry_params_tup = geometry_params_tup,
                        estop_params_tup = estop_params_tup,
                        scatter_params_tup = scatter_params_tup)
statistics.setup(nspec=sim_params.nspec, nbin=sim_params.nbin, limits=sim_params.limits)

if __name__ == "__main__":
    if not config.ENABLE_CACHING:
        print("##### CACHING DISABLED #####")
    if os.environ.get("NUMBA_DISABLE_JIT", '') == "1":
        print("##### NUMBA DISABLED #####")
    
    print("Startup time:", time.time() - start)
    iter_cnt = 1
    counts = [100, 200]
    chunk_size = 100
    # avg_chunk_time = 0.1    # seconds
    
    proj_counts = [[] for _ in range(len(counts))]
    times = [[] for _ in range(len(counts))]
    coefs = read_coefs()
    
    for _ in range(iter_cnt):
        for i, c in enumerate(counts):
            # empty stats for each nion count
            statistics.setup(nspec=sim_params.nspec, nbin=sim_params.nbin, limits=sim_params.limits)
            
            start_time = time.time()
            # proj_count, hist_buf, mom_buf = simulate_adaptive(avg_chunk_time, c, sim_params.to_record(), coefs, follow_recoils=True)
            proj_count, hist_buf, mom_buf = simulate_chunked(chunk_size, c, sim_params.to_record(), coefs, follow_recoils=True)
            # proj_count, hist_buf, mom_buf = simulate(c, sim_params.to_record(), coefs, follow_recoils=True, sim_idx=0)
            statistics.hist.results = hist_buf
            statistics.mom.results = mom_buf
            times[i].append(time.time() - start_time)
            proj_counts[i].append(proj_count)
    print([sum(t)/iter_cnt for t in times], [sum(t)/iter_cnt for t in proj_counts])
    print(times, proj_counts)
    print("--------------------")
    
    # Output the results
    start = time.time()
    statistics.print_results()
    print("Stats time:", time.time() - start)
    statistics.plot_results(log=True)