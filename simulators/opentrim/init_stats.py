"""Initialize the array for statistics.

"""
import numpy as np


HISTOGRAM_1D_DTYPE = np.dtype([
    ("nbins", np.int32),
    ("limits", np.float64, (2,)),
    ("bin_width", np.float64),
    ("bin_counts", np.float64, (120,))
], align=True)

INSIDE_DTYPE = np.dtype([
    ("histx", HISTOGRAM_1D_DTYPE),
], align=True)

STATS_DTYPE = np.dtype([
    ("inside", INSIDE_DTYPE),
], align=True)
