"""Initialize the array for statistics.

"""
import numpy as np


HISTOGRAM_1D_DTYPE = np.dtype([
    ("nbins", np.int32),
    ("limits", np.float64, (2,)),
    ("bin_width", np.float64),
    ("bin_counts", np.float64, (122,))
], align=True)

INSIDE_DTYPE = np.dtype([
    ("histx", HISTOGRAM_1D_DTYPE),
], align=True)

STATS_DTYPE = np.dtype([
    ("inside", INSIDE_DTYPE),
], align=True)


def init_stats(stat_params):
    """Initialize the stats structured array."""
    stats = np.recarray(1, dtype=STATS_DTYPE)
    
    # Initialize the histogram parameters
    stats["inside"]["histx"]["nbins"] = stat_params.nbin
    stats["inside"]["histx"]["limits"] = stat_params.limits
    stats["inside"]["histx"]["bin_width"] = (
        (stats[0].inside.histx.limits[1] - stats[0].inside.histx.limits[0]) 
         / stats[0].inside.histx.nbins)
#        (stat_params.limits[1] - stat_params.limits[0]) / stat_params.nbin)
    stats["inside"]["histx"]["bin_counts"] = 0.0  # Initialize bin counts to zero
    
    return stats