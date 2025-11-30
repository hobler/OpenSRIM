"""Module for collecting and reporting statistics of projectile trajectories.

- Classes are defined for 1D moments and histograms. 
- The moments and histograms are module attributes. 
- Funtions are provided to setup, score, and print/plot results.

There are two module-level attributes:
    mom: Moment_1d instance for calculating moments.
    hist: Histogram_1d instance for calculating histograms.
"""
import math
import numpy as np


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

    Attributes:
        nvar (int): number variables for which moments are desired
        nmax (int): maximum order of moments
        count (NDArray(int), Shape(nvar)): number of scored data points per 
            variable
        _mom (NDArray(float), Shape(nvar, 2*nmax+1)): sum of values raised to 
            the order of the moment to be calculated
        _cenmom (NDArray(float), Shape(nvar, 2*nmax+1)): central moments
    """
    def __init__(self, nvar, nmax):
        self.nvar = nvar
        if nmax < 1 or nmax > 4:
            raise ValueError("nmax must be between 1 and 4.")
        self.nmax = nmax
        self._mom = np.zeros((nvar, 2*nmax + 1), dtype=float)
        self._orders = np.arange(0, 2*nmax + 1, dtype=int)

    def score(self, ivar, value):
        """Score a new data point for variable ivar."""
        self._mom[ivar,:] += value**self._orders[:]

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


class Histogram_1d:
    """Calculate 1D histograms.

    To calculate histograms, create an instance of this class with the
    desired number of variables, number of bins, and limits. Score data points
    using the score() method. The histogram counts can be accessed via the 
    'counts' attribute.
    
    Attributes:
        nvar (int): number variables for which histograms are desired
        nbin (int): number of bins
        limits (Tuple[float], Shape(2)): (min, max) limits of the histogram
        counts (NDArray[int], Shape(nvar,nbin+2)): counts per bin including 
            underflow and overflow bins
        bin_width (float): width of each bin
    """
    def __init__(self, nvar, nbin, limits):
        self.nvar = nvar
        self.nbin = nbin
        self.limits = limits
        self.bin_width = (self.limits[1] - self.limits[0]) / self.nbin
        self.counts = np.zeros((nvar, nbin+2), dtype=int)

    def score(self, ivar, value):
        """Score a new data point to the histogram of variable ivar."""
        if value < self.limits[0]:
            bin_index = 0                # underflow bin
        elif value >= self.limits[1]:
            bin_index = -1               # overflow bin
        else:
            bin_index = int((value - self.limits[0]) / self.bin_width) + 1
        
        self.counts[ivar,bin_index] += 1


def setup(nspec, nbin, limits):
    """Setup module variables.

    Parameters:
        nspec(int): number of atom species
        nbin (int): number of bins
        limits (Tuple[float], Shape(2)): (min, max) limits of the histogram
    """
    global mom, hist

    mom = Moment_1d(nvar=nspec, nmax=4)
    hist = Histogram_1d(nspec, nbin, limits)


def score(proj):
    """Score the final position of the projectile into the histogram.

    Parameters:
        proj (Projectile): the state of the projectile to be scored
    """
    global mom, hist

    mom.score(proj.ispec, proj.pos[2])
    hist.score(proj.ispec, proj.pos[2])


def print_results():
    """Print statistics of the scored projectiles."""
    global mom

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


def plot_results():
    """Plot the histogram using matplotlib."""
    import matplotlib.pyplot as plt

    for ivar in range(hist.nvar):
        plt.stairs(hist.counts[ivar,1:-1],
                   edges=np.linspace(hist.limits[0], hist.limits[1], 
                                     hist.nbin+1),
                   label=f'Species {ivar}')
    plt.xlabel('Penetration depth (A)')
    plt.ylabel('Counts')
    plt.title(f'Histogram of Penetration Depths')
    plt.legend()
    plt.show()

