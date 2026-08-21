import time

from . import config
import numpy as np
from numba import jit, prange, typed, get_thread_id
import numba as nb
from . import cascade
from .mytypes import Projectile, PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .stats import merge_stats, zero_stats
from .save_output import write_stats, save_progress


empty_stats = None


def simulate(nion, params, stats, nion_processed=0):
    """Perform simulation on given number of projectiles

    Wrapper function for the Numba-compiled `_simulate()` function, which 
    performs the actual simulation. This wrapper handles the creation of 
    thread-local statistics buffers and merges them after the simulation.
    The wrapper is necessary because Numba does not support global variables.
        
    Parameters:
        nion: (int) Total number of projectiles to simulate (of this chunk
            in chunked simulations)
        params: (PARAMS_DTYPE) Simulation parameters
        stats: (STATS_DTYPE) Statistical data container to store results in
        nion_processed: (int) Already processed projectiles (used for chunked 
            simulations)
    """
    global empty_stats
    
    if empty_stats is None:
        empty_stats = stats[0].copy()
        zero_stats(empty_stats)

    # Construct an array of stats for each Numba worker thread. Each parallel
    # iteration selects its buffer via get_thread_id(), avoiding concurrent
    # writes to the same statistics record.
    nthreads = nb.get_num_threads()
    stats_per_thread = np.array([empty_stats.copy() for _ in range(nthreads)],
                                dtype=stats.dtype)

    _simulate(nion, params, stats_per_thread, nion_processed)

    #print("Chunk processed")

    # Merge stats from each thread into the total stats
    for i in range(len(stats_per_thread)):
        merge_stats(stats, stats_per_thread[i])

    return


@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, 
     nogil=config.PARALLEL, debug=config.DEBUG)
def _simulate(nion, params, stats_per_thread, nion_processed):
    """Perform simulation on given number of projectiles.
    
    Parameters:
        nion: (int) Total number of projectiles to simulate (in this chunk
            for chunked simulations)
        params: (PARAMS_DTYPE) Simulation parameters
        stats_per_thread: (ndarray[STATS_DTYPE]) Array of stats for each thread
        nion_processed: (int) Already processed projectiles (used for chunked 
            simulations)
    """
    # Initial conditions of the projectile, starting outside target
    dirx = np.cos(np.radians(params[0].beam.tilt))
    diry = np.sin(np.radians(params[0].beam.tilt))
    dirz = 0.0
    xinit = params[0].geometry.x_intf[0] - params[0].cascade.pmax_max
    yinit = xinit * diry / dirx
    zinit = xinit * dirz / dirx
    proj_init = Projectile(
        params[0].beam.energy,  # energy (eV)
        np.array([xinit, yinit, zinit]),  # position (A)
        np.array([dirx, diry, dirz])  # direction (unit vector)
    )
    proj_init_array = np.full(1, proj_init)

    # Container for the returned projectile states
    proj_dummy_list = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_sim = [proj_dummy_list for _ in range(nion)]
    
    # Parallel loop over collision cascades
    for i in prange(nion):  # ty:ignore[not-iterable]
        np.random.seed(params[0].rng_seed + nion_processed + i)
        tid = get_thread_id()
        proj_sim[i] = cascade.cascade(
            proj_init_array[0], params[0], stats_per_thread[tid])

    # TODO: Use returned projectile states
    proj_count = 0
    for proj_lst in proj_sim:
        proj_count += len(proj_lst)
    
    return


def simulate_adaptive(avg_chunk_time, nion, params, stats, input_params=None, 
                      upd_callback=None):
    """Adaptive, chunked simulation, each chunk taking ~avg_chunk_time seconds.
    
    Parameters:
        avg_chunk_time: (int) Desired simulation time per chunk in seconds
        nion: (int) Total number of projectiles to simulate
        params, stats: As in `simulate()`
        input_params (dict): Simulation configuration (for data saving)
        upd_callback (callable): A function to call on simulation data update
    """
    # TODO Doesn't work with fixed seed (due to varying chunk sizes)
    min_chunk_size = 100
    nion_chunksize = min_chunk_size
    nion_processed = 0
    
    while nion_processed < nion:
        nion_chunk = min(nion_chunksize, nion - nion_processed)
        
        start_time = time.time()
        simulate(nion_chunk, params, stats, nion_processed)
        duration = time.time() - start_time
        
        nion_processed += nion_chunk
        if input_params:
            workdir = input_params["simulation"]["workdir"]
            write_stats(params[0], stats, workdir)
            save_progress(workdir, nion_processed, nion)
        if upd_callback:
            upd_callback(nion_processed, nion, stats)
            # TODO: make use of callback to return user stop request; break
        # Calculate optimal chunk size
        new_chunk = int((nion_chunk / duration) * avg_chunk_time)
        nion_chunksize = max(min_chunk_size, new_chunk)
    
    return


def simulate_chunked(nion_chunksize, nion, params, stats, input_params=None, 
                     upd_callback=None):
    """Chunked simulation for nion projectiles
    
    Parameters:
        nion_chunksize: (int) Desired size to split total count into
        nion: (int) Total number of projectiles to simulate
        params, stats: As in `simulate()`
        input_params (dict): Simulation configuration (for data saving)
        upd_callback (callable): A function to call on simulation data update
    """    
    
    def _process_chunks(nion_chunk, nion_processed):
        if nion_chunksize == 0:
            return
        
        simulate(nion_chunk, params, stats, nion_processed)

        if input_params:
            done = nion_processed + nion_chunk
            workdir = input_params["simulation"]["workdir"]
            write_stats(params[0], stats, workdir)
            save_progress(workdir, done, nion)
        if upd_callback:
            upd_callback(nion_processed+nion_chunk, nion, stats)
            # TODO: make use of callback to return user stop request; break
        
    nion_processed = 0
    while nion_processed < nion:
        nion_chunk = min(nion_chunksize, nion - nion_processed)
        _process_chunks(nion_chunk, nion_processed)
        nion_processed += nion_chunk

    return
