import time
from . import config
import numpy as np
from numba import jit, prange
from . import cascade
from .mytypes import Projectile, PROJ_DTYPE, HIST_CONFIG_DTYPE, create_histogram_configs


def simulate(nion, params, hist_configs=None, follow_recoils=False, sim_idx=0):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        hist_configs: (ndarray) Structured array of histogram configurations.
                     If None, creates single histogram from params.stat
        follow_recoils: (bool) If the simulation should be performed for recoils aswell
        sim_idx: (int) Simulation index (for chunked simulations)
        
    Returns:
        tuple[int, np.ndarray, np.ndarray]:
            Total number of simulated projectiles,
            Result buffers for `Histogram_1d` class,
            Result buffers for `Moments_1d` class
    """
    # Create default histogram configuration if not provided
    hist_config_was_none = hist_configs is None
    if hist_config_was_none:
        hist_configs, _ = create_histogram_configs(
            np.array([params.stat.nbin], dtype=np.int32),
            np.array([params.stat.limits[0]], dtype=np.float64),
            np.array([params.stat.limits[1]], dtype=np.float64),
            params.stat.nspec
        )
    assert hist_configs is not None
    
    proj_count, flat_hist_buf_list, mom_buf_list = _simulate(
        nion, params, hist_configs, follow_recoils, sim_idx
    )
    
    # Aggregate histogram buffers across all projectiles
    flat_hist_aggregated = np.sum(flat_hist_buf_list, axis=0, dtype=np.int32)
    mom_aggregated = np.sum(mom_buf_list, axis=0, dtype=np.float64)
    
    return proj_count, flat_hist_aggregated, mom_aggregated.reshape(-1, 9)  # TODO get rid of reshape


@jit(cache=config.ENABLE_CACHING, parallel=config.PARALLEL, nogil=config.PARALLEL)
def _simulate(nion, params, hist_configs, follow_recoils, sim_idx):
    """Perform simulation on given number of projectiles
    
    Parameters:
        nion: (int) Total number of projectiles to simulate
        params: (PARAMS_DTYPE) Simulation parameters
        hist_configs: (ndarray) Structured array of histogram configurations
        follow_recoils: (bool) If the simulation should be performed for recoils aswell
        sim_idx: (int) Simulation index (for chunked simulations)
        
    Returns:
         tuple[int, np.ndarray, np.ndarray]:
             Total number of simulated projectiles,
             2D array of flattened histogram result buffers (shape: nion x flat_counts_size),
             2D array of moment result buffers (shape: nion x mom_size)
    """
    # Calculate total buffer size for histograms
    flat_counts_size = hist_configs["counts_size"].sum()
    
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
    params_arr = np.full(1, params)
    
    # Pre-allocate 2D arrays for histogram and moment results (fixes Numba type inference)
    mom_size = params.stat.nspec * 9  # nvar * (2*nmax+1)
    hist_results_2d = np.zeros((nion, flat_counts_size), dtype=np.int32)
    mom_results_2d = np.zeros((nion, mom_size), dtype=np.float64)
    
    for i in prange(nion):  # ty:ignore[not-iterable]
        np.random.seed(params_arr[0].rng_seed + sim_idx + i)
        proj_sim[i], hist_flat, mom_flat = cascade.cascade(
            proj_dummy_arr[0], params_arr[0], hist_configs, follow_recoils)
        
        # Store histogram and moment results in pre-allocated 2D arrays
        hist_results_2d[i, :] = hist_flat.reshape(flat_counts_size)
        mom_results_2d[i, :] = mom_flat.reshape(mom_size)
    
    proj_count = 0
    for proj_lst in proj_sim:
        proj_count += proj_lst.size
    
    return proj_count, hist_results_2d, mom_results_2d


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
            total_hist_buf = hist_buf.copy()
            total_mom_buf = mom_buf.copy()
        else:
            total_hist_buf += hist_buf
            total_mom_buf += mom_buf
        total_proj_count += proj_count
    
    for processed_count in range(0, nion, chunk_size):
        _process_chunks(chunk_size, processed_count)
    remainder = nion % chunk_size
    _process_chunks(remainder, nion - remainder) # Process remainder
    return total_proj_count, total_hist_buf, total_mom_buf