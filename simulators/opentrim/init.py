"""Initialize the simulation parameters.

- Read the input parameters from a TOML file (not implemented yet, currently 
  hardcoded).
- Calculate derived parameters.
"""
import numpy as np
from . import select_recoil
from . import scatter
from . import cm_scatter
from . import estop
from . import geometry
from . import cascade
from . import stats as statistics
from .mytypes import SimParams
from .nlhlin import read_coefs



def init_params():
    """Initialize the simulation parameters.

    Returns:
        (SimParams): A structured array containing all simulation parameters.
    """

    # Example hardcoded parameters (to be replaced with file input)
    zmin = 0.0              # minimum z coordinate of the target (A)
    zmax = 4000.0           # maximum z coordinate of the target (A)
    pot_model = 'ZBL_magic'  # potential model for scattering
    z1 = 5                  # atomic number of projectile
    m1 = 11.009             # mass of projectile (amu)
    z2 = 14                 # atomic number of target
    m2 = 28.086             # mass of target atom (amu)
    density = 0.04994       # target density (atoms/A^3)
    corr_lindhard1 = 1.5    # Correction factor to Lindhard stopping power (B->Si)
    corr_lindhard2 = 1.0    # Correction factor to Lindhard stopping power (Si->Si)

    nspec = 2              # number of species to record (e.g. projectile and first recoil)
    nbin = 40              # number of bins for depth distribution
    limits = (0.0, 4000.0) # limits for depth distribution

    nlhlin_coefs = read_coefs()
    
    recoil_params = select_recoil.setup(density)
    scatter_params = scatter.setup(z1, m1, z2, m2, pot_model, nlhlin_coefs)
    cm_scatter.setup(n_absc=4)
    estop_params = estop.setup(corr_lindhard1, z1, m1, corr_lindhard2, z2, m2, density)
    geometry_params = geometry.setup(zmin, zmax)
    cascade_params = cascade.setup()
    stat_params = statistics.setup(nspec, nbin, limits)

    params = SimParams( rng_seed = np.random.randint(2**31, dtype=np.uint32),
                        # Seed can be specified manually or generated automatically (default)
                        stat_params = stat_params,
                        cascade_params = cascade_params,
                        recoil_params = recoil_params,
                        geometry_params = geometry_params,
                        estop_params = estop_params,
                        scatter_params = scatter_params)

    statistics.setup(nspec=params.nspec, nbin=params.nbin, limits=params.limits)

    return params