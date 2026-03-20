"""OpenTRIM aims to be a Python implementation of TRIM.

TRIM (Transport of Ions in Matter) is a widely used software package
for simulating the interaction of ions with matter, particularly for
ion implantation in semiconductors. It comes as part of SRIM, see
www.srim.org.

OpenTRIM seeks to replicate the core functionalities of TRIM using
Python, making it more accessible and easier to integrate with other
Python-based tools and workflows. Computational efficiency is achieved 
through the use of Numba for JIT compilation.

Currently, the input parameters are hardcoded in this script, but future
versions may include a more user-friendly interface for specifying
simulation parameters. Also, only depth histograms and moments of the
penetration depth distribution are recorded.
"""
import time
import sys
import os
from pathlib import Path
from numba import jit, prange

if __package__ is None:
    # Running as a script: add parent dir to sys.path
    project_root = str(Path(__file__).parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    __package__ = str(Path(__file__).parent.name)

from . import config  # import config early to set up caching and parallel settings
start = time.time()
from .read_params import read_params  # get the input parameters as a dictionary
from .init_params import get_params  # defines the params structured array
from .stats import init_stats, zero_stats, print_moments, plot_histograms  # defines the stats structured array
from .simulator import simulate, simulate_chunked, simulate_adaptive  # noqa: F401


input_params = read_params()
params, stats_params = get_params(input_params)  # TODO. Remove stats_params
stats = init_stats(params[0].nelem, input_params, stats_params)  # remove stats_params

print("params is in globals():", "params" in globals())
print("stats:", stats)

if __name__ == "__main__":
    if not config.ENABLE_CACHING:
        print("##### CACHING DISABLED #####")
    else:
        print("----- CACHING ENABLED -----")
    if os.environ.get("NUMBA_DISABLE_JIT", "") == "1":
        print("##### NUMBA DISABLED #####")
    else:
        print("----- NUMBA ENABLED -----")
    
    print("Startup time:", time.time() - start)
    iter_cnt = 1
    counts = [10000, 10000] #[10000]
    chunk_size = 100
    # avg_chunk_time = 0.1    # seconds
    
    proj_counts = [[] for _ in range(len(counts))]
    times = [[] for _ in range(len(counts))]
    
    for _ in range(iter_cnt):
        for i, c in enumerate(counts):
            
            zero_stats(stats)

            start_time = time.time()
            # simulate_adaptive(avg_chunk_time, c, params, follow_recoils=True)
            simulate_chunked(chunk_size, c, 
                             params, stats,
                             follow_recoils=True)
            # simulate(c, params, follow_recoils=True, sim_idx=0)
            times[i].append(time.time() - start_time)
            
            proj_count = stats[0]['x']['power_sums'][1,0]  # this is only approximate
            proj_counts[i].append(proj_count)

    # Output the results
    start = time.time()
    print_moments(stats[0])
    end = time.time() - start
    print("--------------------")
    for in_count, times_per_count, out_counts in zip(counts, times, proj_counts):
        print(f"Statistics for {in_count} initial projectiles:")
        if iter_cnt > 1:    # Exclude 1st longer run
            avg_time = sum(times_per_count[1:]) / (iter_cnt - 1)
            avg_iter = iter_cnt - 1
        else:
            avg_time = sum(times_per_count) / iter_cnt
            avg_iter = iter_cnt
        print(f"- Average time ({avg_iter} iterations) [s]:", avg_time)
        print(f"- All simulation times ({iter_cnt} iterations) [s]:", [round(t, 3) for t in times_per_count])
        print("- Average interactions (output projectiles):", sum(out_counts) / iter_cnt)
    print("Stats calculation time [s]:", end)
    print("--------------------")
    plot_histograms(stats[0], log=True)    