"""Initialize the array for statistics.

"""
from collections import namedtuple

import numpy as np
from numba import jit, literal_unroll
from . import config


STATS_DTYPE = None
fields = None


def init_stats(nelem, stat_params):
    """Initialize the stats structured array."""
    global STATS_DTYPE, fields

    HISTOGRAM_1D_DTYPE = np.dtype([
        ("nbins", np.int32),
        ("limits", np.float64, (2,)),
        ("bin_width", np.float64),
        ("counts", np.float64, (nelem, stat_params.nbin + 2))
    ], align=True)
        
    INSIDE_DTYPE = np.dtype([
        ("histx", HISTOGRAM_1D_DTYPE),
    ], align=True)

    STATS_DTYPE = np.dtype([
        ("inside", INSIDE_DTYPE),
    ], align=True)

    # For each STATS_DTYPE field, define the subfields
    fields = (
        ("inside", INSIDE_DTYPE.names),  
    )
    for field, subfields in fields:
        print(f"Field '{field}' has subfields: {subfields}")

    # Create the stats structured array
    # Do not create a structured scalar, since this would cause issues with
    # Numba-jitted functions
    stats = np.empty(1, dtype=STATS_DTYPE)
    
    # Initialize the histogram parameters
    stats["inside"]["histx"]["nbins"] = stat_params.nbin
    stats["inside"]["histx"]["limits"] = stat_params.limits
    stats["inside"]["histx"]["bin_width"] = (
        (stat_params.limits[1] - stat_params.limits[0]) / stat_params.nbin)
    stats["inside"]["histx"]["counts"].fill(0.0)
    
    return stats


def zero_stats(stats):
    """Reset the statistics to zero.
    
    We must use direct access to the fields here, since generating a list of 
    field names cannot be done using compile-time constants.
    """
    for field, subfields in fields:
        for subfield in subfields:
            stats[field][subfield]["counts"].fill(0.0)

    return stats


def merge_stats(stats, stat):
    """Merge the statistics from one projectile into the total statistics.
    
    Field names of structured arrays cannot be queried in Numba-jitted 
    functions, so we have to hardcode the field names here.
    """
    for field, subfields in fields:
        for subfield in subfields:
            if subfield.startswith("hist"):
                stats_counts = stats[field][subfield]["counts"]
                stat_counts = stat[field][subfield]["counts"]
                stats_counts += stat_counts


@jit(cache=config.ENABLE_CACHING)
def score(stats, proj):
    """Score a projectile's contribution to the statistics."""
    ivar = proj["ielem"]

    if proj["is_inside"]:
        # Determine the bin index for the projectile's x position
        x = proj["pos"][0]
        hist = stats["inside"]["histx"]
        if x < hist["limits"][0]:
            ibin = 0  # Underflow bin
        elif x < hist["limits"][1]:
            ibin = int((x - hist["limits"][0]) / hist["bin_width"]) + 1
        else:
            ibin = -1  # Overflow bin
        hist["counts"][ivar, ibin] += 1.0


def plot_results(stats, log=False):
    """Plot the histogram using matplotlib."""
    import matplotlib.pyplot as plt

    hist = stats["inside"]["histx"]
    for ivar in range(hist["counts"].shape[0]):
        plt.stairs(hist["counts"][ivar, 1:-1],
                   edges=np.linspace(hist["limits"][0], hist["limits"][1], 
                                     hist["nbins"]+1),
                   label=f"Species {ivar}, Hist 'depth distribution'")
    if log:
        plt.yscale("log")
    plt.xlabel("Penetration depth (A)")
    plt.ylabel("Counts")
    plt.title("Histogram of Penetration Depths")
    plt.legend()
    plt.show()
