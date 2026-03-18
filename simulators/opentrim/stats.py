"""Handle moments and histograms.

In order to allow JIT compilation of the scoring function, we need to define 
a structured array stats the may be passed through the Python-nopython
interface.

Available functions:
    init_stats: Initialize the stats structured array.
    zero_stats: Reset the statistics to zero.
    merge_stats: Merge the statistics from one projectile into the total 
        statistics.
    score: Score a projectile's contribution to the statistics.
    standardize_moments: Calculate the standardized moments from the power sums.
    print_moments: Print the standardized moments.
    plot_histograms: Plot the histogram using matplotlib.
"""
import math
import numpy as np
from numba import jit, literal_unroll
from . import config


STATS_DTYPE = None
fields = None


def init_stats(nelem, stat_params):
    """Initialize the stats structured array.
    
    The stats array contains subarrays for all the moments and histograms.

    Arguments:
        nelem: The number of different atom species
        stat_params: The parameters for the statistics

    Returns:
        The initialized stats structured array.
    """
    global STATS_DTYPE, fields

    max_order = 4    # TODO: Make max_order an input parameter, passed via stat_params
    assert 1 <= max_order <= 4, "max_order must be between 1 and 4"
    MOMENTS_1D_DTYPE = np.dtype([
        ("nvar", np.int32),
        ("max_order", np.int32),
        ("power_sums", np.float64, (nelem, 2*max_order + 1,))
    ], align=True)

    HISTOGRAM_1D_DTYPE = np.dtype([
        ("nvar", np.int32),
        ("nbins", np.int32),
        ("limits", np.float64, (2,)),
        ("bin_width", np.float64),
        ("counts", np.float64, (nelem, stat_params.nbin + 2))
    ], align=True)
        
    INSIDE_DTYPE = np.dtype([
        ("momx", MOMENTS_1D_DTYPE),
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
    
    # Initialize the moments parameters
    stats["inside"]["momx"]["nvar"] = nelem
    stats["inside"]["momx"]["max_order"] = max_order
    stats["inside"]["momx"]["power_sums"].fill(0.0)
    
    # Initialize the histogram parameters
    stats["inside"]["histx"]["nvar"] = nelem
    stats["inside"]["histx"]["nbins"] = stat_params.nbin
    stats["inside"]["histx"]["limits"] = stat_params.limits
    stats["inside"]["histx"]["bin_width"] = (
        (stat_params.limits[1] - stat_params.limits[0]) / stat_params.nbin)
    stats["inside"]["histx"]["counts"].fill(0.0)
    
    return stats


def zero_stats(stats):
    """Reset the statistics to zero.
    
    Loop over the fields and subfields of the stats structured array and reset 
    the power sums and histogram counts to zero.
    
    Arguments:
        stats: The stats structured array to be reset (modified in-place)
    
    Returns: 
        The reset stats structured array.
    """
    for field, subfields in fields:
        for subfield in subfields:
            if subfield.startswith("mom"):
                stats[field][subfield]["power_sums"].fill(0.0)
            elif subfield.startswith("hist"):
                stats[field][subfield]["counts"].fill(0.0)

    return stats


def merge_stats(total_stats, stats):
    """Merge the statistics from one projectile into the total statistics.
    
    Arguments:
        total_stats: The total statistics to be updated (modified in-place)
        stats: The statistics from a single projectile to be merged into the 
            total statistics
    """
    for field, subfields in fields:
        for subfield in subfields:
            if subfield.startswith("mom"):
                total_mom_values = total_stats[field][subfield]["power_sums"]
                mom_values = stats[field][subfield]["power_sums"]
                total_mom_values += mom_values
            elif subfield.startswith("hist"):
                total_hist_counts = total_stats[field][subfield]["counts"]
                hist_counts = stats[field][subfield]["counts"]
                total_hist_counts += hist_counts


@jit
def score(stats, proj):
    """Score a projectile's contribution to the statistics."""
    ivar = proj["ielem"]

    if proj["is_inside"]:
        x = proj["pos"][0]

        mom = stats["inside"]["momx"]
        max_order = mom["max_order"]
        increment = x ** np.arange(2*max_order + 1)
        mom["power_sums"][ivar, :] += increment

        hist = stats["inside"]["histx"]
        if x < hist["limits"][0]:
            ibin = 0  # Underflow bin
        elif x < hist["limits"][1]:
            ibin = int((x - hist["limits"][0]) / hist["bin_width"]) + 1
        else:
            ibin = -1  # Overflow bin
        hist["counts"][ivar, ibin] += 1.0


def standardize_moments(mom, ivar):
    """Calculate the standardized moments from the power sums.
    
    We define the standardized moments (abbreviated as std_moements) here as 
    count, mean, standard deviation, skewness, and kurtosis, although the term 
    is nomally used only for the latter two.
    
    Arguments:
        mom: The moments structured array containing the power sums
        ivar: The index of the variable (atom species) for which to calculate 
            the moments
    
    Returns:
        The standardized moments.
    """
    max_order = mom["max_order"]
    std_moments = np.zeros(mom["max_order"] + 1)
    std_moments_err = np.zeros(mom["max_order"] + 1)

    # Counts
    power_sums = mom["power_sums"][ivar, :]
    count = power_sums[0]
    std_moments[0] = count
    if count == 0:
        return std_moments, std_moments_err

    # Raw moments
    moments = power_sums[:] / count

    # Central moments
    central_moments = np.array(
        [sum(math.comb(i, j) * moments[i-j] * (-moments[1])**j 
             for j in range(i+1))
         for i in range(2*max_order + 1)]
    )

    # Central moments errors
    central_moments_err = np.array(
        [np.sqrt((central_moments[2*i]
                  - 2*i*central_moments[i-1]*central_moments[i+1] 
                  - central_moments[i]**2 
                  + i**2*central_moments[2]*central_moments[i-1]**2)
                 / count) for i in range(max_order + 1)]
    )

    # Mean value
    std_moments[1] = moments[1]
    std_moments_err[1] = np.sqrt(central_moments[2] / count)
    if max_order == 1:
        return std_moments, std_moments_err
    
    # standard deviation
    std_moments[2] = np.sqrt(central_moments[2])
    if std_moments[2] == 0.0:
        return std_moments, std_moments_err
    
    std_moments_err[2] = central_moments_err[2] / (2*std_moments[2]) 
    if max_order == 2 or std_moments[2] == 0.0:
        return std_moments, std_moments_err

    # skewness and kurtosis
    for i in range(3, max_order + 1):
        std_moments[i] = central_moments[i] / std_moments[2]**i
        std_moments_err[i] = central_moments_err[i] / std_moments[2]**i

    return std_moments, std_moments_err


def print_moments(stats):
    """Print the standardized moments."""
    mom = stats["inside"]["momx"]
    max_order = mom["max_order"]
    for ivar in range(mom["nvar"]):
        std_moments, std_moments_err = standardize_moments(mom, ivar)
        
        print(f"Statistics for atom species {ivar}:")

        print(f"   Number of atoms stopped inside the target: "
              f"{std_moments[0]:.0f}")
        print(f"   Mean penetration depth: "
              f"{std_moments[1]:.2f} A +/- {std_moments_err[1]:.2f} A")
        if max_order >= 2:
            print(f"   Standard deviation of penetration depth: "
                  f"{std_moments[2]:.2f} A +/- {std_moments_err[2]:.2f} A")
        if max_order >= 3:
            print(f"   Skewness: "
                  f"{std_moments[3]:.2f} +/- {std_moments_err[3]:.2f}")
        if max_order >= 4:
            print(f"   Kurtosis: "
                  f"{std_moments[4]:.2f} +/- {std_moments_err[4]:.2f}")


def plot_histograms(stats, log=False):
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
