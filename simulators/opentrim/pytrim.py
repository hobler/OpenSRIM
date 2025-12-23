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
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        sim_params_tup: (tuple) Simulation parameters (provided by `SimParams.to_tuple()`)
        follow_recoils: (bool) If the simulation should be performed for recoils aswell
        
    Returns:
        tuple[int, np.ndarray, np.ndarray]:
            Total number of simulated projectiles,
            Result buffer for `Histogram_1d` class,
            Result buffer for `Moment_1d` class
    """
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
    return proj_count, hist.results, mom.results

def simulate_adaptive(avg_chunk_time, nion, *args, **kwargs):
    """Adaptive, chunked simulation with each chunk taking around avg_sim_time seconds
    
    Parameters:
        avg_sim_time: (int) Desired simulation time per in seconds
        nion: (int) Total number of projectiles to simulate
        *args, **kwargs: As in `simulate()`
        
    Returns:
        tuple[int, np.ndarray, np.ndarray]:
            Total number of simulated projectiles,
            Result buffer for `Histogram_1d` class,
            Result buffer for `Moment_1d` class
    """
    min_chunk_size = 100
    chunk_size = min_chunk_size
    processed_count = 0
    
    total_proj_count = 0
    total_hist_buf = None
    total_mom_buf = None

    while processed_count < nion:
        current_batch = min(chunk_size, nion - processed_count)
        
        start_time = time.time()
        proj_count, hist_buf, mom_buf = simulate(current_batch, *args, **kwargs)
        # NOTE: Saving can be performed here
        
        total_proj_count += proj_count
        if total_hist_buf is None:
            total_hist_buf = hist_buf.copy()
            total_mom_buf = mom_buf.copy()
        else:
            total_hist_buf += hist_buf
            total_mom_buf += mom_buf  # pyright: ignore[reportOperatorIssue]
        duration = time.time() - start_time
        
        processed_count += current_batch
        # Calculate optimal chunk size
        new_chunk = int((current_batch / duration) * avg_chunk_time)
        chunk_size = max(min_chunk_size, new_chunk)
    
    return total_proj_count, total_hist_buf, total_mom_buf

def simulate_chunked(chunk_size, nion, *args, **kwargs):
    """Chunked simulation for nion projectiles
    
    Parameters:
        chunk_size: (int) Size to split total count into
        nion: (int) Total number of projectiles to simulate
        *args, **kwargs: As in `simulate()`
        
    Returns:
        tuple[int, np.ndarray, np.ndarray]:
            Total number of simulated projectiles,
            Result buffer for `Histogram_1d` class,
            Result buffer for `Moment_1d` class
    """
    assert nion % chunk_size == 0, "Total projectile count must be a multiple of chunk size"
    
    total_proj_count = 0
    total_hist_buf = None
    total_mom_buf = None
    for _ in range(0, nion, chunk_size):
        proj_count, hist_buf, mom_buf = simulate(chunk_size, *args, **kwargs)
        # NOTE: Saving can be performed here
        
        if total_hist_buf is None:
            total_hist_buf = hist_buf.copy()
            total_mom_buf = mom_buf.copy()
        else:
            total_hist_buf += hist_buf
            total_mom_buf += mom_buf  # pyright: ignore[reportOperatorIssue]
        total_proj_count += proj_count
    return total_proj_count, total_hist_buf, total_mom_buf

if __name__ == "__main__":
    times = []
    proj_counts = []
    counts = [1000, 10000]
    # chunk_size = 100
    avg_chunk_time = 0.1    # seconds
    simulate(10, sim_params.to_tuple(), follow_recoils=True)    # pre-compile
    for i, c in enumerate(counts):
        # empty stats for each nion count
        statistics.setup(nspec=sim_params.nspec, nbin=sim_params.nbin, limits=sim_params.limits)
        
        start_time = time.time()
        proj_count, hist_buf, mom_buf = simulate_adaptive(avg_chunk_time, c, sim_params.to_tuple(), follow_recoils=True)
        # proj_count, hist_buf, mom_buf = simulate_chunked(chunk_size, c, sim_params.to_tuple(), follow_recoils=True)
        statistics.hist.results = hist_buf
        statistics.mom.results = mom_buf
        times.append(time.time() - start_time)
        proj_counts.append(proj_count)
    print(times, proj_counts)
    
    # Output the results
    statistics.print_results()
    statistics.plot_results(log=True)
