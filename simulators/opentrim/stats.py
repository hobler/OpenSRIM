"""Handle moments and histograms.

In order to allow JIT compilation of the scoring functions, we need to define 
a structured array "stats" that may be passed through the Python-nopython
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


def init_stats(NELEM_ION, NELEM_TARGET, input_params):
    """Initialize the stats structured array.
    
    The stats array contains subarrays for all the moments and histograms.

    Arguments:
        NELEM_ION: (int) The maximum number of different ion atom species
        NELEM_TARGET: (int) The maximum number of different target atom species
        input_params: (dict) The input parameters dictionary, used to extract 
            the histogram configurations from the output configuration section.
    
    Returns:
        The initialized stats structured array.
    """
    global STATS_DTYPE, stats_fields

    follow_recoils = input_params["simulation"]["follow_recoils"]
    NELEM = NELEM_ION + NELEM_TARGET

    # Define short names for the statistics, which will be used to construct the 
    # keys of the stats structured array 
    short_names = {
        "depth_distribution": "x",
        "lateral_distribution": "y",
        "backscattered_atoms": "b",
        "transmitted_atoms": "t",
        "ion_recoils": "",
        "nuclear_energy_deposition": "n",
        "electronic_energy_deposition": "e",
        "energy": "e",
        "angle": "a",
    }

    # Number of variables to reserve memory for
    nvar = {
        "x": NELEM + NELEM_TARGET if follow_recoils else NELEM,
        "y": NELEM + NELEM_TARGET if follow_recoils else NELEM,
        "xy": NELEM + NELEM_TARGET if follow_recoils else NELEM,
        "xn": NELEM,
        "yn": NELEM,
        "xyn": NELEM,
        "xe": NELEM,
        "ye": NELEM,
        "xye": NELEM,
        "be": NELEM if follow_recoils else NELEM_ION,
        "ba": NELEM if follow_recoils else NELEM_ION,
        "te": NELEM if follow_recoils else NELEM_ION,
        "ta": NELEM if follow_recoils else NELEM_ION,
    }
    # Build a flattened dictionary "stats_configs" of statistics configurations 
    # from the output confiuration defined in the input parameters. The keys of 
    # the dictionary are constructed from the keys in the input parameters, e.g. 
    #   "depth distribution.ion/recoils" -> "x", 
    #   "lateral distribution.nuclear energy deposition" -> "yn", 
    #   "backscattered atoms distribution.energy" -> "be", etc. 
    # The histogram parameters (number of bins, limits, etc.) are taken from 
    # the corresponding section in the input parameters. 
    output_params = input_params["output"]
    stats_configs = {}
    
    for name in output_params:
        if name in ["trajectories", "distribution_2d"]:
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

    paired_configs = {
        "xy": ("x", "y"),
        "xyn": ("xn", "yn"),
        "xye": ("xe", "ye"),
    }
    distribution_2d = output_params["distribution_2d"]
    distribution_2d_names = {
        "xy": "ion_recoils",
        "xyn": "nuclear_energy_deposition",
        "xye": "electronic_energy_deposition",
    }
    for short_name in paired_configs:
        cfg = distribution_2d[distribution_2d_names[short_name]]
        nbins = cfg["nbins"]
        limits = cfg["limits"]
        stats_configs[short_name] = {
            "score": bool(cfg["score"]),
            "x_nbins": int(nbins[0]),
            "y_nbins": int(nbins[1]),
            "x_limits": tuple(limits[0]),
            "y_limits": tuple(limits[1]),
        }

    # Build a structured array data type of statistics parameters
    for i, short_name in enumerate(stats_configs):
        stats_config = stats_configs[short_name]
        if short_name in paired_configs:
            stats_dtype = np.dtype([
                ("score", np.int64),
                ("nvar", np.int32),
                ("x_nbins", np.int32),
                ("y_nbins", np.int32),
                ("_pad", np.int32),
                ("x_limits", np.float64, (2,)),
                ("y_limits", np.float64, (2,)),
                ("x_bin_width", np.float64),
                ("y_bin_width", np.float64),
                ("counts", np.float64, (nvar[short_name],
                                        stats_config["x_nbins"] + 2,
                                        stats_config["y_nbins"] + 2)),
            ], align=True)
        else:
            stats_dtype = np.dtype([
                ("score", np.int64),  # whether to score this statistics (use integer for JIT compatibility)
                ("nvar", np.int32),  # number of variables (e.g. atom species) for this statistics
                ("nbins", np.int32),  # number of bins for this statistics
                ("limits", np.float64, (2,)),  # limits for this statistics
                ("bin_width", np.float64),
                ("counts", np.float64, (nvar[short_name],
                                        stats_config["nbins"] + 2)),
                ("power_sums", np.float64, (nvar[short_name], 2*max_order + 1)),
            ], align=True)

        if i == 0:
            STATS_DTYPE = np.dtype([
                (short_name, stats_dtype),
            ], align=True)
        else:
            STATS_DTYPE = np.dtype(STATS_DTYPE.descr + [
                (short_name, stats_dtype),
            ], align=True)
    STATS_DTYPE = np.dtype(STATS_DTYPE.descr, align=True)

    # Create the structured array of statistics
    stats = np.recarray(1, dtype=STATS_DTYPE)
    for short_name, stats_config in stats_configs.items():
        stats[0][short_name]["nvar"] = nvar[short_name]
        if short_name in paired_configs:
            stats[0][short_name]["score"] = stats_config["score"]
            stats[0][short_name]["x_nbins"] = stats_config["x_nbins"]
            stats[0][short_name]["y_nbins"] = stats_config["y_nbins"]
            stats[0][short_name]["_pad"] = 0
            stats[0][short_name]["x_limits"] = stats_config["x_limits"]
            stats[0][short_name]["y_limits"] = stats_config["y_limits"]
            stats[0][short_name]["x_bin_width"] = (
                (stats_config["x_limits"][1] - stats_config["x_limits"][0])
                / stats_config["x_nbins"])
            stats[0][short_name]["y_bin_width"] = (
                (stats_config["y_limits"][1] - stats_config["y_limits"][0])
                / stats_config["y_nbins"])
            stats[0][short_name]["counts"].fill(0.0)
            continue

        for field in stats_config:
            if field not in ["score", "nbins", "limits"]:
                raise ValueError(f"Unknown histogram config field: {field} "
                                 f"in {short_name}")
            stats[0][short_name][field] = stats_config[field]

        stats[0][short_name]["bin_width"] = (
            (stats_config["limits"][1] - stats_config["limits"][0]) 
            / stats_config["nbins"])
        if "counts" in stats[0][short_name].dtype.names:
            stats[0][short_name]["counts"].fill(0.0)
        stats[0][short_name]["power_sums"].fill(0.0)

    if False:
        print("-----------")
        print(f"stats: {stats}")
        print(f"stats fields: {stats.dtype.names}")
        for short_name in stats.dtype.names:
            print(f"   {short_name}: {stats[short_name]}")

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
        if "counts" in stats[field].dtype.names:
            stats[field]["counts"].fill(0.0)
        if "power_sums" in stats[field].dtype.names:
            stats[field]["power_sums"].fill(0.0)


def merge_stats(total_stats, stats):
    """Merge the statistics from one projectile into the total statistics.
    
    Arguments:
        total_stats: The total statistics to be updated (modified in-place)
        stats: The statistics from a single projectile to be merged into the 
            total statistics
    """
    for field in stats_fields:
            if "power_sums" in total_stats[field].dtype.names:
                total_power_sums = total_stats[field]["power_sums"]
                power_sums = stats[field]["power_sums"]
                total_power_sums += power_sums
            if "counts" in total_stats[field].dtype.names:
                total_counts = total_stats[field]["counts"]
                counts = stats[field]["counts"]
                total_counts += counts


@jit(debug=config.DEBUG)
def _score1d(stats_distribution, value, ivar, weight=1.0):
    """Score a projectile's contribution to a statistics distribution.
    
    Arguments:
        stats: The stats structured array to be updated (modified in-place)
        value: The value to be scored (e.g. penetration depth)
        ivar: The index of the variable (e.g. atom species) for which to score
        weight: The weight of this contribution (default 1.0)
    """
    if stats_distribution["score"]:
        for i in range(2*max_order + 1):
            stats_distribution["power_sums"][ivar, i] += weight * value ** i

        if value < stats_distribution["limits"][0]:
            ibin = 0  # Underflow bin
        elif value < stats_distribution["limits"][1]:
            ibin = int((value - stats_distribution["limits"][0]) / 
                        stats_distribution["bin_width"]) + 1
        else:
            ibin = -1  # Overflow bin
        stats_distribution["counts"][ivar, ibin] += weight


@jit(debug=config.DEBUG)
def _score2d(stats_distribution, x, y, ivar, weight=1.0):
    if stats_distribution["score"]:
        if x < stats_distribution["x_limits"][0]:
            ix = 0
        elif x < stats_distribution["x_limits"][1]:
            ix = int((x - stats_distribution["x_limits"][0]) /
                     stats_distribution["x_bin_width"]) + 1
        else:
            ix = -1

        if y < stats_distribution["y_limits"][0]:
            iy = 0
        elif y < stats_distribution["y_limits"][1]:
            iy = int((y - stats_distribution["y_limits"][0]) /
                     stats_distribution["y_bin_width"]) + 1
        else:
            iy = -1

        stats_distribution["counts"][ivar, ix, iy] += weight


@jit(debug=config.DEBUG)
def _score_stop(stats, proj):
    """Score a projectile that has stopped inside the target."""
    ivar = proj["ielem"]
    x = proj["pos"][0]
    y = proj["pos"][1]

    _score1d(stats["x"], x, ivar)
    _score1d(stats["y"], y, ivar)
    _score2d(stats["xy"], x, y, ivar)


@jit(debug=config.DEBUG)
def _score_backscattered(stats, proj):
    """Score a backscattered projectile."""
    ivar = proj["ielem"]
    energy = proj["e"]
    angle = math.degrees(math.atan2(proj["dir"][1], -proj["dir"][0]))

    _score1d(stats["be"], energy, ivar)
    _score1d(stats["ba"], angle, ivar)


@jit(debug=config.DEBUG)
def _score_transmitted(stats, proj):
    """Score a transmitted projectile."""
    ivar = proj["ielem"]
    energy = proj["e"]
    angle = math.degrees(math.atan2(proj["dir"][1], proj["dir"][0]))

    _score1d(stats["te"], energy, ivar)
    _score1d(stats["ta"], angle, ivar)


@jit(debug=config.DEBUG)
def score_eed(stats, proj, dee):
    """Score the electronic energy deposition for a projectile."""
    ivar = proj["ielem"]
    x = proj["pos"][0]  # TODO: Take center of point and previous point
    y = proj["pos"][1]

    _score1d(stats["xe"], x, ivar, weight=dee)
    _score1d(stats["ye"], y, ivar, weight=dee)
    _score2d(stats["xye"], x, y, ivar, weight=dee)


@jit(debug=config.DEBUG)
def score_ned(stats, proj):
    """Score the nuclear energy deposition for a projectile."""
    ivar = proj["ielem"]
    x = proj["pos"][0]
    y = proj["pos"][1]
    ned = proj["e"]

    _score1d(stats["xn"], x, ivar, weight=ned)
    _score1d(stats["yn"], y, ivar, weight=ned)
    _score2d(stats["xyn"], x, y, ivar, weight=ned)


@jit(debug=config.DEBUG)
def score_start(stats, proj, nelem_target):
    """Score a projectile at its starting position."""
    ivar = proj["ielem"] + nelem_target
    x = proj["pos"][0]
    y = proj["pos"][1]

    _score1d(stats["x"], x, ivar)
    _score1d(stats["y"], y, ivar)
    _score2d(stats["xy"], x, y, ivar)


@jit(debug=config.DEBUG)
def score_stop(stats, proj):
    """Score a projectile that has stopped inside the target."""
    score_ned(stats, proj)
    _score_stop(stats, proj)


@jit(debug=config.DEBUG)
def score_exit(stats, proj):
    """Score a projectile that has exited the target."""
    if proj["dir"][0] < 0:
        _score_backscattered(stats, proj)
    else:
        _score_transmitted(stats, proj)


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

    for field in stats.dtype.names:
        if not stats[field]["score"]:
            continue
        if "nbins" not in stats[field].dtype.names:
            continue
        hist = stats[field]
        for ivar in range(hist["counts"].shape[0]):
            plt.stairs(hist["counts"][ivar, 1:-1],
                    edges=np.linspace(hist["limits"][0], hist["limits"][1], 
                                        hist["nbins"]+1),
                    label=f"Species {ivar}, Hist '{field}'")
        if log:
            plt.yscale("log")
        if field.startswith("x"):
            label = "x (A)"
        elif field.startswith("y"):
            label = "y (A)"
        elif field.endswith("e"):
            label = "Energy (eV)"
        elif field.endswith("a"):
            label = "Angle (degrees)"
        else:
            label = field
        plt.xlabel(label)
        plt.ylabel("Counts")
        plt.title("OpenSRIM")
        plt.legend()
        plt.show()
