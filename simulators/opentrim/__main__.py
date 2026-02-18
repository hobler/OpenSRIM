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
import sys
import os
from pathlib import Path

if __package__ is None:
    # Running as a script: add parent dir to sys.path
    project_root = str(Path(__file__).parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    __package__ = str(Path(__file__).parent.name)

from . import config
from . import select_recoil
from . import scatter
from . import cm_scatter
from . import estop
from . import geometry
from . import cascade
from . import stats as statistics
from .mytypes import SimParams
from .nlhlin import read_coefs
from .simulator import simulate, simulate_chunked, simulate_adaptive  # noqa: F401

import numpy as np

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
nlhlin_coefs = read_coefs()
recoil_params_tup = select_recoil.setup(density)
scatter_params_tup = scatter.setup(z1, m1, z2, m2, pot_model, nlhlin_coefs)
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
    else:
        print("----- CACHING ENABLED -----")
    if os.environ.get("NUMBA_DISABLE_JIT", '') == "1":
        print("##### NUMBA DISABLED #####")
    else:
        print("----- NUMBA ENABLED -----")
    
    print("Startup time:", time.time() - start)
    iter_cnt = 1
    counts = [10000]
    chunk_size = 100
    # avg_chunk_time = 0.1    # seconds
    
    proj_counts = [[] for _ in range(len(counts))]
    times = [[] for _ in range(len(counts))]
    
    for _ in range(iter_cnt):
        for i, c in enumerate(counts):
            # empty stats for each nion count
            statistics.setup(nspec=sim_params.nspec, nbin=sim_params.nbin, limits=sim_params.limits)
            
            start_time = time.time()
            # proj_count, hist_buf, mom_buf = simulate_adaptive(avg_chunk_time, c, sim_params.to_record(), follow_recoils=True)
            proj_count, hist_buf, mom_buf = simulate_chunked(chunk_size, c, sim_params.to_record(), follow_recoils=True)
            # proj_count, hist_buf, mom_buf = simulate(c, sim_params.to_record(), follow_recoils=True, sim_idx=0)
            if statistics.hist is not None:
                statistics.hist.results = hist_buf
            if statistics.mom is not None:
                statistics.mom.results = mom_buf
            times[i].append(time.time() - start_time)
            proj_counts[i].append(proj_count)
    
    # Output the results
    start = time.time()
    statistics.print_results()
    end = time.time() - start
    print("--------------------")
    for in_count, times_per_count, out_counts in zip(counts, times, proj_counts):
        print(f"Statistics for {in_count} initial projectiles:")
        print(f"- Average time ({iter_cnt} iterations) [s]:", sum(times_per_count) / iter_cnt)
        print(f"- All simulation times ({iter_cnt} iterations) [s]:", [round(t, 3) for t in times_per_count])
        print("- Average interactions (output projectiles):", sum(out_counts) / iter_cnt)
    print("Stats calculation time [s]:", end)
    print("--------------------")
    statistics.plot_results(log=True)