import time
from . import config
import numpy as np
from numba import jit, prange, typed, int32
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
        tuple[int, list[np.ndarray], np.ndarray]:
            Total number of simulated projectiles,
            A list of buffers (for each hist) for `Histogram_1d` class,
            Result buffers for `Moments_1d` class
    """
    #print(f"params is C contiguous = {params.flags.c_contiguous}")
    #print(f"Size of params: {params.nbytes/1024:.3f} kB")

    proj_count, hist_buf, mom_buf = _simulate(nion, params, follow_recoils, sim_idx)
    #_simulate.inspect_types()  # For debugging Numba type inference issues
    for i, lst in enumerate(hist_buf):
        if i == 0:
            continue
        for j, arr in enumerate(lst):
            hist_buf[0][j] += arr
    return proj_count, hist_buf[0], np.sum(mom_buf, axis=0, dtype=np.float64)


@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, nogil=config.PARALLEL)
def _simulate(nion, params, follow_recoils, sim_idx):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        follow_recoils: (bool) If the simulation should be performed for recoils as well
        sim_idx: (int) Simulation index (for chunked simulations)
        
    Returns:
        tuple[int, list[np.ndarray], list[np.ndarray]]:
            Total number of simulated projectiles,
            List of result buffers for `Histogram_1d` class (for each `nion`),
            List of result buffers for `Moments_1d` class (for each `nion`)
    """
    # Initial conditions of the projectile
    proj_init = Projectile(
        params[0].beam.energy,  # energy (eV)
        np.array([0.0, 0.0, 0.0]),  # position (A)
        np.array([np.sin(np.radians(params[0].beam.tilt)), 0.0, 
                  np.cos(np.radians(params[0].beam.tilt))]), # direction (unit vector)
        0,
        0,
        True
    )
    proj_dummy = np.empty(1, dtype=PROJ_DTYPE)
    proj_sim = [proj_dummy for _ in range(nion)]
    proj_dummy[0] = proj_init

    hist_dummy = np.empty((1, 1), dtype=np.int32)
    mom_dummy = np.empty((1, 1), dtype=np.float64)
    # hist_results: [[hist1[:, :], hist2[:, :], ...], ...]
    # where len(hist_results) == nion
    hist_results = typed.List.empty_list(typed.List.empty_list(int32[:,:]))
    for _ in range(nion):
        hist_results.append(typed.List.empty_list(int32[:,:]))
    mom_results = [mom_dummy for _ in range(nion)]
    
    # Simulate the collision cascades in parallel
    for i in prange(nion):  # ty:ignore[not-iterable]
        np.random.seed(params[0].rng_seed + sim_idx + i)
        proj_sim[i], hist_results[i], mom_results[i] = cascade.cascade(
            proj_dummy[0], params[0], follow_recoils)
    
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
        tuple[int, list[np.ndarray], np.ndarray]:
            Total number of simulated projectiles,
            A list of buffers (for each hist) for `Histogram_1d` class,
            Result buffers for `Moments_1d` class
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
            total_hist_buf = hist_buf
            total_mom_buf = mom_buf.copy()
        else:
            for i, arr in enumerate(hist_buf):
                total_hist_buf[i] += arr
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
        tuple[int, list[np.ndarray], np.ndarray]:
            Total number of simulated projectiles,
            A list of buffers (for each hist) for `Histogram_1d` class,
            Result buffers for `Moments_1d` class
    """    
    total_proj_count = 0
    total_hist_buf = None
    total_mom_buf = None
    
    def _process_chunks(chunk_size, sim_idx):
        nonlocal total_hist_buf, total_mom_buf, total_proj_count
        if chunk_size == 0:
            return
        
        proj_count, hist_buf, mom_buf = simulate(
            chunk_size,
            *args,
            sim_idx=sim_idx,
            **kwargs
        )
        # NOTE: Saving can be performed here
        
        if total_hist_buf is None:
            total_hist_buf = hist_buf
            total_mom_buf = mom_buf.copy()
        else:
            for i, arr in enumerate(hist_buf):
                total_hist_buf[i] += arr
            total_mom_buf += mom_buf
        total_proj_count += proj_count
    
    for processed_count in range(0, nion, chunk_size):
        _process_chunks(chunk_size, processed_count)
    remainder = nion % chunk_size
    _process_chunks(remainder, nion - remainder) # Process remainder
    return total_proj_count, total_hist_buf, total_mom_buf
