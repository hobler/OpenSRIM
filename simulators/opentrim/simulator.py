import time
from . import config
import numpy as np
from numba import jit, prange
from . import cascade
from .mytypes import Projectile, PROJ_DTYPE


def simulate(nion, params, follow_recoils=False, sim_idx=0):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        follow_recoils: (bool) If the simulation should be performed for recoils aswell
        sim_idx: (int) Simulation index (for chunked simulations)
        
    Returns:
        tuple[int, np.ndarray, np.ndarray]:
            Total number of simulated projectiles,
            Result buffers for `Histogram_1d` class,
            Result buffers for `Moments_1d` class
    """
    #print(f"params is C contiguous = {params.flags.c_contiguous}")
    #print(f"Size of params: {params.nbytes/1024:.3f} kB")

    proj_count, hist_buf, mom_buf = _simulate(nion, params, follow_recoils, sim_idx)
    return proj_count, np.sum(hist_buf, axis=0, dtype=np.int32), np.sum(mom_buf, axis=0, dtype=np.float64)


@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, nogil=config.PARALLEL)
def _simulate(nion, params, follow_recoils, sim_idx):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        sim_params_tup: (tuple) Simulation parameters (provided by `SimParams.to_tuple()`)
        follow_recoils: (bool) If the simulation should be performed for recoils aswell
        sim_idx: (int) Simulation index (for chunked simulations)
        
    Returns:
        tuple[int, list[np.ndarray], list[np.ndarray]]:
            Total number of simulated projectiles,
            List of result buffers for `Histogram_1d` class (for each `nion`),
            List of result buffers for `Moments_1d` class (for each `nion`)
    """
    print(f"index={sim_idx}, type(params)={type(params)}")

    # Initial conditions of the projectile
    proj_init = Projectile(
        50000.0,                         # energy (eV)
        np.array([0.0, 0.0, 0.0]),     # position (A)
        np.array([0.0, 0.0, 1.0]),     # direction (unit vector)
        0,
        0,
        True
    )
    proj_dummy_arr = np.empty(1, dtype=PROJ_DTYPE)
    proj_sim = [proj_dummy_arr for _ in range(nion)]
    proj_dummy_arr[0] = proj_init

    # Fixes weird Numba error by passing array instead of single record
    #params_arr = np.full(1, params)
    
    hist_dummy = np.empty((1, 1), dtype=np.int32)
    mom_dummy = np.empty((1, 1), dtype=np.float64)
    hist_results = [hist_dummy for _ in range(nion)]
    mom_results = [mom_dummy for _ in range(nion)]
    
    # Simulate the collision cascades in parallel
    for i in prange(nion):
        np.random.seed(params.rng_seed + sim_idx + i)
        proj_sim[i], hist_results[i], mom_results[i] = cascade.cascade(
            proj_dummy_arr[0], params, follow_recoils)
    
    proj_count = 0
    for proj_lst in proj_sim:
        proj_count += proj_lst.size
    
    return proj_count, hist_results, mom_results


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
    # TODO Doesn't work with fixed seed (due to varying chunk sizes)
    min_chunk_size = 100
    chunk_size = min_chunk_size
    processed_count = 0
    
    total_proj_count = 0
    total_hist_buf = None
    total_mom_buf = None

    while processed_count < nion:
        current_batch = min(chunk_size, nion - processed_count)
        
        start_time = time.time()
        proj_count, hist_buf, mom_buf = simulate(
            current_batch,
            *args,
            sim_idx=processed_count,
            **kwargs
        )
        # NOTE: Saving or adding data to queue can be performed here
        
        total_proj_count += proj_count
        if total_hist_buf is None:
            total_hist_buf = hist_buf.copy()
            total_mom_buf = mom_buf.copy()
        else:
            total_hist_buf += hist_buf
            total_mom_buf += mom_buf
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
    total_proj_count = 0
    total_hist_buf = None
    total_mom_buf = None
    
    def _process_chunks(chunk_size):
        nonlocal total_hist_buf, total_mom_buf, total_proj_count
        if chunk_size == 0:
            return
        
        proj_count, hist_buf, mom_buf = simulate(
            chunk_size,
            *args,
            sim_idx=processed_count,
            **kwargs
        )
        # NOTE: Saving can be performed here
        
        if total_hist_buf is None:
            total_hist_buf = hist_buf.copy()
            total_mom_buf = mom_buf.copy()
        else:
            total_hist_buf += hist_buf
            total_mom_buf += mom_buf
        total_proj_count += proj_count
    
    for processed_count in range(0, nion, chunk_size):
        _process_chunks(chunk_size)
    _process_chunks(nion // chunk_size) # Process remainder
    return total_proj_count, total_hist_buf, total_mom_buf