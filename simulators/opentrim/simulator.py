import time

from . import config
import numpy as np
from numba import jit, prange, typed, int32
import numba as nb
from . import cascade
from .mytypes import Projectile, PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .stats import STATS_DTYPE, merge_stats, zero_stats


empty_stats = None


def simulate(nion, params, stats, sim_idx=0):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        stats: (STATS_DTYPE) Statistical data container to store results in
        sim_idx: (int) Simulation index (for chunked simulations)
    """
    global empty_stats
    
    if empty_stats is None:
        empty_stats = stats[0].copy()
        zero_stats(empty_stats)

    # Construct an array of stats for each ion, since lists cannot be used in 
    # Numba-jitted functions
    stats_per_ion = np.array([empty_stats.copy() for _ in range(nion)], 
                             dtype=STATS_DTYPE)

    _simulate(nion, params, stats_per_ion, sim_idx)

    # Merge stats from each ion into the total stats
    for i in range(len(stats_per_ion)):
        merge_stats(stats, stats_per_ion[i])

    return


@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, nogil=config.PARALLEL)
def _simulate(nion, params, stats_per_ion, sim_idx):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        stats_per_ion: (ndarray[STATS_DTYPE]) Array of stats for each ion
        sim_idx: (int) Simulation index (for chunked simulations)
    """
    # Initial conditions of the projectile
    proj_init = Projectile(
        params[0].beam.energy,  # energy (eV)
        np.array([0.0, 0.0, 0.0]),  # position (A)
        np.array([np.cos(np.radians(params[0].beam.tilt)), 
                  np.sin(np.radians(params[0].beam.tilt)), 0.0]),
                  # direction (unit vector)
        0,
        0,
        True
    )
    proj_dummy_list = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_sim = [proj_dummy_list for _ in range(nion)]
    
    proj_dummy = np.full(1, proj_init)

    # Parallel loop over collision cascades
    for i in prange(nion):  # ty:ignore[not-iterable]
        np.random.seed(params[0].rng_seed + sim_idx + i)
        proj_sim[i] = cascade.cascade(
            proj_dummy[0], params[0], stats_per_ion[i])
    
    proj_count = 0
    for proj_lst in proj_sim:
        proj_count += len(proj_lst)
    
    return


def simulate_adaptive(avg_chunk_time, nion, *args, **kwargs):
    """Adaptive, chunked simulation with each chunk taking around avg_sim_time seconds
    
    Parameters:
        avg_sim_time: (int) Desired simulation time per in seconds
        nion: (int) Total number of projectiles to simulate
        *args, **kwargs: As in `simulate()`
    """
    # TODO Doesn't work with fixed seed (due to varying chunk sizes)
    min_chunk_size = 100
    chunk_size = min_chunk_size
    processed_count = 0
    
    while processed_count < nion:
        current_batch = min(chunk_size, nion - processed_count)
        
        start_time = time.time()
        simulate(
            current_batch,
            *args,
            sim_idx=processed_count,
            **kwargs
        )
        # NOTE: Saving or adding data to queue can be performed here
        
        duration = time.time() - start_time
        
        processed_count += current_batch
        # Calculate optimal chunk size
        new_chunk = int((current_batch / duration) * avg_chunk_time)
        chunk_size = max(min_chunk_size, new_chunk)
    
    return


def simulate_chunked(chunk_size, nion, *args, **kwargs):
    """Chunked simulation for nion projectiles
    
    Parameters:
        chunk_size: (int) Size to split total count into
        nion: (int) Total number of projectiles to simulate
        *args, **kwargs: As in `simulate()`
    """    
    
    def _process_chunks(chunk_size, sim_idx):
        if chunk_size == 0:
            return
        
        simulate(
            chunk_size,
            *args,
            sim_idx=sim_idx,
            **kwargs
        )
        # NOTE: Saving can be performed here
        
    for processed_count in range(0, nion, chunk_size):
        _process_chunks(chunk_size, processed_count)
    remainder = nion % chunk_size
    _process_chunks(remainder, nion - remainder) # Process remainder

    return
