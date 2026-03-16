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
from . import stats_old as statistics
start = time.time()
from .init_params import get_params  # defines the params structured array
from .stats import init_stats, plot_results  # defines the stats structured array
from .simulator import simulate, simulate_chunked, simulate_adaptive  # noqa: F401


params = get_params()
stats = init_stats(params[0].nelem, params[0].stats[0])

print("params is in globals():", "params" in globals())
print("stats:", stats)

@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, nogil=config.PARALLEL)
def test_params(params):
    print("Testing params:")
#    print("params:", params)
    for _ in prange(2):
        print("params[0].rng_seed:", params[0].rng_seed)
        print("params[0].stats:", params[0].stats)
        print("params[0].cascade:", params[0].cascade)
        print("params[0].recoil:", params[0].recoil)
        print("params[0].geometry:", params[0].geometry)
        print("params[0].elements:", params[0].elements)
        print("params[0].materials:", params[0].materials)
        print("params[0].estop:", params[0].estop)
        print("params[0].scatter:", params[0].scatter)
        print("params[0].nelem:", params[0].nelem)
        print("Finished!")

#test_params(params); exit()

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
    counts = [10000] #[10000]
    chunk_size = 100
    # avg_chunk_time = 0.1    # seconds
    
    proj_counts = [[] for _ in range(len(counts))]
    times = [[] for _ in range(len(counts))]
    
    for _ in range(iter_cnt):
        for i, c in enumerate(counts):
            # empty stats for each nion count
            statistics.setup(params[0].stats)
            
            start_time = time.time()
            # proj_count, hist_buf, mom_buf = simulate_adaptive(avg_chunk_time, c, params, follow_recoils=True)
            proj_count, hist_buf, mom_buf = simulate_chunked(chunk_size, c, 
                                                             params, stats,
                                                             follow_recoils=True)
            # proj_count, hist_buf, mom_buf = simulate(c, params, follow_recoils=True, sim_idx=0)
            if statistics.stat is not None:
                statistics.stat.results = (hist_buf, mom_buf)
            times[i].append(time.time() - start_time)
            proj_counts[i].append(proj_count)

    print(f'stats["inside"]["histx"]["counts"]=', stats[0]["inside"]["histx"]["counts"][:, 1:-1])

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
    plot_results(stats[0], log=True)    