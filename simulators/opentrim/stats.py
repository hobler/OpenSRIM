"""Handle moments and histograms.

In order to allow JIT compilation of the scoring functions, we need to define 
a structured array "stats" the may be passed through the Python-nopython
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
stats_fields = None

max_order = 4    # TODO: Make max_order an input parameter, passed via input_params
                 # e.g., input_params["output"]["depth distribution"].max_order
assert 1 <= max_order <= 4, "max_order must be between 1 and 4"


def init_stats(nelem, input_params):
    """Initialize the stats structured array.
    
    The stats array contains subarrays for all the moments and histograms.

    Arguments:
        nelem: (int) The number of different atom species
        input_params: (dict) The input parameters dictionary, used to extract 
            the histogram configurations from the output configuration section.
    
    Returns:
        The initialized stats structured array.
    """
    global STATS_DTYPE, stats_fields

    short_names = {
        "depth distribution": "x",
        "lateral distribution": "y",
        "backscattered atoms distribution": "b",
        "transmitted atoms distribution": "t",
        "ion/recoils": "",
        "nuclear energy deposition": "n",
        "electronic energy deposition": "e",
        "energy": "e",
        "angle": "a",
    }
    # Build a flattened dictionary "stats_configs" of statistics configurations 
    # from the output configuration. The keys of the dictionary are constructed 
    # from the keys in the output configuration, e.g. 
    #   "depth distribution.ion/recoils" -> "x", 
    #   "lateral distribution.nuclear energy deposition" -> "yn", 
    #   "backscattered atoms distribution.energy" -> "be", etc. 
    # The histogram parameters (number of bins, limits, etc.) are taken from 
    # the corresponding section in the output configuration. 
    output_params = input_params["output"]
    stats_configs = {}
    
    for name in output_params:
        if name == "trajectories":
            continue
        if name not in short_names:
            raise ValueError(f"Unknown output field: {name}")
        if not isinstance(output_params[name], dict):
            raise ValueError(f"Expected a dictionary for output field "
                             f"{name}")
        
        for subname in output_params[name]:
            if subname not in short_names:
                raise ValueError(f"Unknown output type: {name}.{subname}")        
           
            stats_config = output_params[name][subname]
            if not isinstance(stats_config, dict):
                raise ValueError(f"Expected a dictionary for output field "
                                 f"{name}.{subname}")

            short_name = f"{short_names[name]}{short_names[subname]}"
            stats_configs[short_name] = stats_config

    # Build a structured array data type of statistics parameters
    for i, name in enumerate(stats_configs):
        stats_dtype = np.dtype([
            ("score", np.int64),  # whether to score this statistics (use integer for JIT compatibility)
            ("nvar", np.int32),  # number of variables (e.g. atom species) for this statistics
            ("nbins", np.int32),  # number of bins for this statistics
            ("limits", np.float64, (2,)),  # limits for this statistics
            ("bin_width", np.float64),
            ("counts", np.float64, (nelem,  # TODO: may depend on follow_recoils and other factors
                                    stats_configs[name]["nbins"] + 2)),
            ("power_sums", np.float64, (nelem, 2*max_order + 1)),
        ], align=True)

        if i == 0:
            STATS_DTYPE = np.dtype([
                (name, stats_dtype),
            ], align=True)
        else:
            STATS_DTYPE = np.dtype(STATS_DTYPE.descr + [
                (name, stats_dtype),
            ], align=True)
    STATS_DTYPE = np.dtype(STATS_DTYPE.descr, align=True)

    # Create the structured array of statistics
    stats = np.recarray(1, dtype=STATS_DTYPE)
    for name, stats_config in stats_configs.items():
        stats[0][name]["nvar"] = nelem    # TODO: may depend on follow_recoils and other factors
        for field in stats_config:
            if field not in ["score", "nbins", "limits"]:
                raise ValueError(f"Unknown histogram config field: {field} "
                                 f"in {name}")
            stats[0][name][field] = stats_config[field]

        stats[0][name]["bin_width"] = (
            (stats_config["limits"][1] - stats_config["limits"][0]) 
            / stats_config["nbins"])
        stats[0][name]["counts"].fill(0.0)
        stats[0][name]["power_sums"].fill(0.0)

    if False:
        print("-----------")
        print(f"stats: {stats}")
        print(f"stats fields: {stats.dtype.names}")
        for name in stats.dtype.names:
            print(f"   {name}: {stats[name]}")

    stats_fields = STATS_DTYPE.names

#    print(f"stats_fields: '{stats_fields}'")
    
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
    for field in stats_fields:
        stats[field]["counts"].fill(0.0)
        stats[field]["power_sums"].fill(0.0)


def merge_stats(total_stats, stats):
    """Merge the statistics from one projectile into the total statistics.
    
    Arguments:
        total_stats: The total statistics to be updated (modified in-place)
        stats: The statistics from a single projectile to be merged into the 
            total statistics
    """
    for field in stats_fields:
            total_power_sums = total_stats[field]["power_sums"]
            power_sums = stats[field]["power_sums"]
            total_power_sums += power_sums
            total_counts = total_stats[field]["counts"]
            counts = stats[field]["counts"]
            total_counts += counts


@jit
def score(stats, proj):
    """Score a projectile's contribution to the statistics."""
    ivar = proj["ielem"]

    if proj["is_inside"] and stats["x"]["score"]:
        x = proj["pos"][0]

        for i in range(2*max_order + 1):
            stats["x"]["power_sums"][ivar, i] += x ** i

        if x < stats["x"]["limits"][0]:
            ibin = 0  # Underflow bin
        elif x < stats["x"]["limits"][1]:
            ibin = int((x - stats["x"]["limits"][0]) / stats["x"]["bin_width"]) + 1
        else:
            ibin = -1  # Overflow bin
        stats["x"]["counts"][ivar, ibin] += 1.0


def standardize_moments(stats, ivar):
    """Calculate the standardized moments from the power sums.
    
    We define the standardized moments (abbreviated as std_moements) here as 
    count, mean, standard deviation, skewness, and kurtosis, although the term 
    is nomally used only for the latter two.
    
    Arguments:
        stats: The statistics structured array containing the power sums
        ivar: The index of the variable (atom species) for which to calculate 
            the moments
    
    Returns:
        The standardized moments.
    """
    std_moments = np.zeros(max_order + 1)
    std_moments_err = np.zeros(max_order + 1)

    # Counts
    power_sums = stats["power_sums"][ivar, :]
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
    for ivar in range(stats["x"]["nvar"]):
        std_moments, std_moments_err = standardize_moments(stats["x"], ivar)
        
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

    hist = stats["x"]
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
