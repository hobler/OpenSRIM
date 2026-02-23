"""Module for collecting and reporting statistics of projectile trajectories.

- Classes are defined for 1D moments and histograms. 
- The moments and histograms are module attributes. 
- Funtions are provided to setup, score, and print/plot results.

There are two module-level attributes:
    mom: Moment_1d instance for calculating moments.
    hist: Histogram_1d instance for calculating histograms.
"""
import math
from numba.core.types import UniTuple
import numpy as np
from numba.experimental import jitclass
from numba.extending import overload, register_jitable
from numba import int32, float64, jit, from_dtype
from .mytypes import HIST_CONFIG_DTYPE, create_histogram_configs


mom = None
hist = None

@register_jitable
def fct(n: int):
    """Calculate the factorial of a number (Numba-compatible)
    
    Parameters:
        n: (int) Number to calculate the factorial for
        
    Returns:
        int: The factorial of the provided number
    """
    fact = 1
    for num in range(2, n + 1):
        fact *= num
    return fact

@overload(math.comb)
def impl(n: int, k: int):
    """Calculate k-of-n combinations (Numba-compatible) 
    
    Parameters:
        k: (int) Number of items in the chosen set
        n: (int) Total number of items to choose from
        
    Returns:
        int: The number of possible combinations
    """
    def inner(n: int, k: int):
        return 0 if k > n else fct(n) / (fct(k) * fct(n - k))
    return inner
    

@jitclass(spec = [
    ("nvar", int32),
    ("nmax", int32),
    ("count", float64[:]),
    ("_orders", int32[:]),
    ("_mom", float64[:,:]),
    ("_cenmom", float64[:,:])
])
class Moment_1d:
    """Calculate moments of 1D data.
    
    Moments up to twice the desired order are stored as they are needed 
    to compute the standard deviations of the moments.

    To calculate one of the standard moments (mean, std, skewness, kurtosis),
    create an instance of this class with the desired number of variables and
    maximum order of moments. Score data points using the score() method. After
    all data points have been scored, call the central_moments() method to
    compute the central moments from the raw moments. Finally, call the desired
    moment methods (mean(), std(), skewness(), kurtosis()) to retrieve the
    moments and their standard errors. 

    Note that the standard errors are derived assuming independent data points.
    If they are correlated such as when recoil positions stem from the same
    collision cascade, the standard errors will be underestimated.

    Attributes:
        nvar (int): number of variables for which moments are desired
        nmax (int): maximum order of moments
        count (ndarray[int]): number of scored data points per variable (size 
            nvar)
        _orders (ndarray[int]): orders of all moments to be calculated (size 
            2*nmax+1)
        _mom (ndarray[float]): sum of values raised to the order of the moment 
            to be calculated (shape (nvar, 2*nmax+1))
        results: (ndarray[float]) A public getter / setter for `_mom`
        _cenmom (ndarray[float]): central moments (shape (nvar, 2*nmax+1))
    """
    def __init__(self, nvar, nmax):
        self.nvar = nvar
        if nmax < 1 or nmax > 4:
            raise ValueError("nmax must be between 1 and 4.")
        self.nmax = nmax
        self._orders = np.arange(0, 2*nmax + 1, dtype=np.int32)
        self._mom = np.zeros((nvar, 2*nmax + 1), dtype=np.float64)
        # TODO variable init?
        # self._cenmom = np.zeros((nvar, 2*nmax + 1), dtype=np.float64)
        # self.count = np.zeros(nvar, dtype=np.float64)

    def score(self, ivar, value):
        """Score a new data point for variable ivar."""
        # Original line causing __powidf2 missing error:
        # self._mom[ivar,:] += value**self._orders[:]

        # Workaround
        powers_of_value = np.empty_like(self._orders, dtype=np.float64)
        powers_of_value[0] = 1.0
        if len(self._orders) > 1:
            powers_of_value[1] = value

        # Compute higher powers iteratively using multiplication
        # This avoids the generic float-to-float power function
        for i in range(2, len(self._orders)):
            powers_of_value[i] = powers_of_value[i-1] * value

        self._mom[ivar,:] += powers_of_value
        
    @property
    def results(self):
        return self._mom
        
    @results.setter
    def results(self, new):
        self._mom = new


    def central_moments(self):
        """Compute central moments up to order 2*nmax."""
        self.count = self._mom[:,0]
        mom = self._mom[:,:] / self.count[:,np.newaxis]
        self._cenmom = np.zeros_like(self._mom)
        for ivar in range(self.nvar):
            for i in range(2*self.nmax + 1):
                self._cenmom[ivar,i] = mom[ivar,i]
                for j in range(1, i+1):
                    self._cenmom[ivar,i] += (math.comb(i, j) 
                            * mom[ivar,i-j] * (-mom[ivar,1])**j)

    def mean(self):
        """Return the mean values and their standard errors."""
        mean_ = self._mom[:,1] / self.count[:]
        mean_err = np.sqrt(self._cenmom[:,2] / self.count[:]) 
        return mean_[:], mean_err[:]
    
    def std(self):
        """Return the standard deviations and their standard errors."""
        std_ = np.sqrt(self._cenmom[:,2])
        std_err = self._cenmom_err(2)[:] / (2*std_[:])
        return std_[:], std_err[:]
    
    def skewness(self):
        """Return the skewnesses and their standard errors."""
        skewness_ = self._cenmom[:,3] / self._cenmom[:,2]**1.5
        skewness_err = self._cenmom_err(3)[:] / self._cenmom[:,2]**1.5
        return skewness_[:], skewness_err[:]
    
    def kurtosis(self):
        """Return the kurtoses and their standard errors."""
        kurtosis_ = self._cenmom[:,4] / self._cenmom[:,2]**2
        kurtosis_err = self._cenmom_err(4)[:] / self._cenmom[:,2]**2
        return kurtosis_[:], kurtosis_err[:]

    def _cenmom_err(self, i):
        """Return the standard error of the central moment of order i."""
        cenmom_err = np.sqrt((self._cenmom[:,2*i]
                              - 2*i*self._cenmom[:,i-1]*self._cenmom[:,i+1] 
                              - self._cenmom[:,i]**2 
                              + i**2*self._cenmom[:,2]*self._cenmom[:,i-1]**2)
                             / self.count[:])
        return cenmom_err[:]


@jitclass(spec = [
    ("nvar", int32),
    ("n_hist", int32),
    ("hist_configs", from_dtype(HIST_CONFIG_DTYPE)[:]),
    ("flat_counts", int32[:]),
    ("total_counts_size", int32)
])
class Histogram_1d:
    """Calculate multiple 1D histograms with different parameters.

    To calculate histograms, create an instance of this class with the
    desired number of variables, and an array of histogram configurations.
    Score data points using the score() method, which updates all histogram
    buffers simultaneously. Individual histogram counts can be extracted
    using the get_histogram_counts() method.
    
    This class supports multiple histogram buffers with different binning
    schemes and limits while maintaining Numba JIT cache compatibility.

    Attributes:
        nvar (int): number variables for which histograms are desired
        n_hist (int): number of histogram configurations
        hist_configs (ndarray): structured array of histogram configurations
        flat_counts (ndarray[int]): flattened counts for all histograms
    """
    def __init__(self, nvar, hist_configs):
        self.nvar = nvar
        self.n_hist = len(hist_configs)
        self.hist_configs = hist_configs
        self.total_counts_size = hist_configs["counts_size"].sum()
        self.flat_counts = np.zeros(self.total_counts_size, dtype=np.int32)

    def score(self, ivar, value):
        """Score a new data point into ALL histogram buffers.
        
        Parameters:
            ivar: Variable index (species)
            value: Data point value to score
        """
        for ihist in range(self.n_hist):
            config = self.hist_configs[ihist]
            nbin = config["nbin"]
            limits_min = config["limits_min"]
            limits_max = config["limits_max"]
            bin_width = config["bin_width"]
            offset = config["offset"]
            
            # Calculate bin index
            if value < limits_min:
                ibin = 0
            elif value >= limits_max:
                ibin = -1
            else:
                ibin = int((value - limits_min) / bin_width) + 1
            
            # Calculate index in flattened array
            flat_idx = offset + ivar * (nbin + 2) + ibin
            if ibin == -1:
                flat_idx = offset + ivar * (nbin + 2) + (nbin + 1)
            
            self.flat_counts[flat_idx] += 1
    
    def get_histogram_counts(self, ihist):
        """Extract counts for a specific histogram.
        
        Parameters:
            ihist: Histogram index (0 to n_hist-1)
        
        Returns:
            np.ndarray: Counts of shape (nvar, nbin+2) for the specified histogram
        """
        config = self.hist_configs[ihist]
        offset = config["offset"]
        nbin = config["nbin"]
        counts_size = config["counts_size"]
        
        counts = self.flat_counts[offset:offset+counts_size].copy()
        return counts.reshape(self.nvar, nbin + 2)
    
    @property
    def results(self):
        return self.flat_counts
        
    @results.setter
    def results(self, new):
        self.flat_counts = new


def setup(nspec, nbin, limits):
    """Setup module variables and pre-compile functions

    Parameters:
        nspec(int): number of atom species
        nbin (int): number of bins
        limits (float[2]): [min, max] limits of the histogram (size 2)

    Returns:
        (STAT_PARAMS_DTYPE): Statistics parameters
    """
    global mom, hist

    STAT_PARAMS_DTYPE = np.dtype([
        ("nspec", np.int32),
        ("nbin", np.int32),
        ("limits", np.float64, (2,)),
    ], align=True)

    stat_params = np.recarray(1, dtype=STAT_PARAMS_DTYPE)[0]
    stat_params["nspec"] = nspec
    stat_params["nbin"] = nbin
    stat_params["limits"] = np.array(limits)

    # Create single histogram configuration for backward compatibility
    hist_configs, flat_size = create_histogram_configs(
        np.array([nbin], dtype=np.int32),
        np.array([limits[0]], dtype=np.float64),
        np.array([limits[1]], dtype=np.float64),
        nspec
    )

    mom = Moment_1d(nvar=nspec, nmax=4)
    hist = Histogram_1d(nspec, hist_configs)
    
    mom.central_moments()
    mom.mean()
    mom.std()
    mom.skewness()
    mom.kurtosis()

    return stat_params


def print_results():
    """Print statistics of the scored projectiles."""
    global mom
    assert mom is not None

    mom.central_moments()
    mean, mean_err = mom.mean()
    std, std_err = mom.std()
    skewness, skewness_err = mom.skewness()
    kurtosis, kurtosis_err = mom.kurtosis()

    for ivar in range(mom.nvar):
        print(f"Statistics for atom species {ivar}:")
    
        if mom.count[ivar] == 0:
            print("   No atoms stopped inside the target.")
            continue

        print(f"   Number of atoms stopped inside the target: "
              f"{int(mom.count[ivar])}")
        print(f"   Mean penetration depth: "
              f"{mean[ivar]:.2f} A +/- {mean_err[ivar]:.2f} A")
        print(f"   Standard deviation of penetration depth: "
              f"{std[ivar]:.2f} A +/- {std_err[ivar]:.2f} A")
        print(f"   Skewness: "
              f"{skewness[ivar]:.2f} +/- {skewness_err[ivar]:.2f}")
        print(f"   Kurtosis: "
              f"{kurtosis[ivar]:.2f} +/- {kurtosis_err[ivar]:.2f}")


def plot_results(log=False):
    """Plot the histogram using matplotlib.
    
    For multiple histograms, this will plot all configurations.
    """
    import matplotlib.pyplot as plt
    assert hist is not None

    for ihist in range(hist.n_hist):
        for ivar in range(hist.nvar):
            counts = hist.get_histogram_counts(ihist)
            config = hist.hist_configs[ihist]
            
            plt.stairs(counts[ivar,1:-1],
                      edges=np.linspace(config["limits_min"], config["limits_max"], 
                                        config["nbin"]+1),
                      label=f"Species {ivar}, Hist {ihist}")
    if log:
        plt.yscale("log")
    plt.xlabel("Penetration depth (A)")
    plt.ylabel("Counts")
    plt.title("Histogram of Penetration Depths")
    plt.legend()
    plt.show()

