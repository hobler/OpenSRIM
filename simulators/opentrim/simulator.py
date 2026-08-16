import time

from . import config
import numpy as np
from numba import jit, prange, typed, get_thread_id
import numba as nb
from . import cascade
from .mytypes import Projectile, PROJ_DTYPE, PROJ_NUMBA_DTYPE
from .stats import merge_stats, zero_stats
from .process_data import write_stats, save_progress


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

    # Construct an array of stats for each Numba worker thread. Each parallel
    # iteration selects its buffer via get_thread_id(), avoiding concurrent
    # writes to the same statistics record.
    nthreads = nb.get_num_threads()
    stats_per_thread = np.array([empty_stats.copy() for _ in range(nthreads)],
                                dtype=stats.dtype)

    _simulate(nion, params, stats_per_thread, sim_idx)

    #print("Chunk processed")

    # Merge stats from each thread into the total stats
    for i in range(len(stats_per_thread)):
        merge_stats(stats, stats_per_thread[i])

    return


@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, 
     nogil=config.PARALLEL, debug=config.DEBUG)
def _simulate(nion, params, stats_per_thread, sim_idx):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        stats_per_thread: (ndarray[STATS_DTYPE]) Array of stats for each thread
        sim_idx: (int) Simulation index (for chunked simulations)
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
    proj_dummy_list = typed.List.empty_list(PROJ_NUMBA_DTYPE)
    proj_sim = [proj_dummy_list for _ in range(nion)]
    
    proj_dummy = np.full(1, proj_init)

    # Parallel loop over collision cascades
    for i in prange(nion):  # ty:ignore[not-iterable]
        np.random.seed(params[0].rng_seed + sim_idx + i)
        tid = get_thread_id()
        proj_sim[i] = cascade.cascade(
            proj_dummy[0], params[0], stats_per_thread[tid])
    
    proj_count = 0
    for proj_lst in proj_sim:
        proj_count += len(proj_lst)
    
    return


def simulate_adaptive(avg_chunk_time, nion, params, stats, input_params=None, upd_callback=None):
    """Adaptive, chunked simulation with each chunk taking around avg_sim_time seconds
    
    Parameters:
        avg_sim_time: (int) Desired simulation time per in seconds
        nion: (int) Total number of projectiles to simulate
        params, stats: As in `simulate()`
        input_params (dict): Simulation configuration (for data saving)
        upd_callback (callable): A function to call on simulation data update
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
            params,
            stats,
            processed_count
        )
        duration = time.time() - start_time
        
        processed_count += current_batch
        if input_params:
            write_stats(input_params, stats)
            save_progress(input_params, processed_count, nion)
        if upd_callback:
            upd_callback(processed_count, nion, stats)
            # TODO: make use of callback to return user stop request; break
        # Calculate optimal chunk size
        new_chunk = int((current_batch / duration) * avg_chunk_time)
        chunk_size = max(min_chunk_size, new_chunk)
    
    return


def simulate_chunked(chunk_size, nion, params, stats, input_params=None, 
                     upd_callback=None):
    """Chunked simulation for nion projectiles
    
    Parameters:
        chunk_size: (int) Size to split total count into
        nion: (int) Total number of projectiles to simulate
        params, stats: As in `simulate()`
        input_params (dict): Simulation configuration (for data saving)
        upd_callback (callable): A function to call on simulation data update
    """    
    
    def _process_chunks(chunk_size, sim_idx):
        if chunk_size == 0:
            return
        
        simulate(
            chunk_size,
            params,
            stats,
            sim_idx
        )
        if input_params:
            done = sim_idx + chunk_size
            write_stats(input_params, stats)
            save_progress(input_params, done, nion)
        if upd_callback:
            upd_callback(sim_idx+chunk_size, nion, stats)
            # TODO: make use of callback to return user stop request; break
        
    processed_count = 0
    while processed_count < nion:
        current_batch = min(chunk_size, nion - processed_count)
        _process_chunks(current_batch, processed_count)
        processed_count += current_batch

    return
