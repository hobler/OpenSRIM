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

if __package__ is None:
    # Running as a script: add parent dir to sys.path
    project_root = str(Path(__file__).parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    __package__ = str(Path(__file__).parent.name)

from . import config  # import config early to set up caching and parallel settings
from . import stats as statistics
from .init import init_params
from .simulator import simulate, simulate_chunked  # noqa: F401

start = time.time()
params = init_params()

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
    counts = [10000]
    chunk_size = 100
    # avg_chunk_time = 0.1    # seconds
    
    proj_counts = [[] for _ in range(len(counts))]
    times = [[] for _ in range(len(counts))]
    
    for _ in range(iter_cnt):
        for i, c in enumerate(counts):
            # empty stats for each nion count
            statistics.setup(nspec=params.stat.nspec, nbin=params.stat.nbin, limits=tuple(params.stat.limits))
            
            start_time = time.time()
            # proj_count, hist_buf, mom_buf = simulate_adaptive(avg_chunk_time, c, params, follow_recoils=True)
            proj_count, hist_buf, mom_buf = simulate_chunked(chunk_size, c, params, follow_recoils=True)
            # proj_count, hist_buf, mom_buf = simulate(c, params, follow_recoils=True, sim_idx=0)
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