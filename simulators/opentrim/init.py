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
from .nlhlin import read_coefs


def read_input():
    """Read the input parameters from a file.

    Currently hardcoded, but should read from a TOML file in the future.

    Returns:
        (dict): A dictionary containing the input parameters.
    """
    input_params = {
        "simulation": {
            "nions": 10000,        # number of ions to simulate
            "nions_update": 100,     # update parameters every n ions
            "rng_seed": 12345,      # random seed for reproducibility
            "workdir": "./",          # working directory for output files
        },
        "beam": {
            "symbol": "B",       # chemical symbol of the incoming ions
            "name": "Boron",      # full name of chemical element
            "Z": 5,              # atomic number of the incoming ions
            "M": 11.009,          # mass of the incoming ions (amu)
            "energy": 50.0,       # energy of the incoming ions (keV)
            "tilt": 0.0,         # tilt angle of the beam (degrees)
        },
        "layers" : {
            "name": ["Layer 1"],  # names of the layers
            "width": [4000.0],    # width of each layer (A)
            "density": [0.04994],  # density of each layer (atoms/A^3)
            "compound correction": [1.0],  # correction factor for compound targets
            "gas": [False],        # whether the layer is a gas (True) or solid (False) 
            "material": [
                {
                "symbol": ["Si"],  # chemical symbol of the target atoms
                "name": ["Silicon"], # full name of chemical element
                "Z": [14],          # atomic number of the target atoms
                "M": [28.086],       # mass of the target atoms (amu)
                "stoichiometry": [1],   # stoichiometric ratio of the target atoms in the layer
                "displacement_energy": [15.0],  # displacement energy of the target atoms (eV)
                }
            ],
        },
        "models": {
            "potential": "ZBL",  # potential model for scattering
            "scattering integrals": {
                "algorithm": "magic",  # algorithm for numerical integration of scattering integrals
                                                # "magic" or "Guass-Legendre"
                "n_absc": 4,       # number of abscissas for numerical integration of scattering integrals
            },
            "electronic stopping": "Lindhard",  # model for electronic stopping power
            "Lindhard correction": {
                "B->Si": 1.5,       # Correction factor to Lindhard stopping power for B->Si
                "Si->Si": 1.0,      # Correction factor to Lindhard stopping power for Si->Si
            },
        },
        "output": {
            "trajectories": {
                "start": False,      # whether to record starting points of trajectories
                "end": False,        # whether to record ending points of trajectories
                "collisions": False,  # whether to record collision points of trajectories
            },
            "depth distribution": {
                "nbins": 40,         # number of bins for depth distribution
                "limits": (0.0, 4000.0), # limits for depth distribution (A)
                "ion/recoils": True,  # whether to record depth distribution for both ions and recoils
                "phonons": False,      # whether to record depth distribution for phonons
                "ionization": False,      # whether to record depth distribution for ionization events
            },
            "lateral distribution": {
                "nbins": 40,         # number of bins for lateral distribution
                "limits": (-2000.0, 2000.0), # limits for lateral distribution (A)
                "ion/recoils": True,  # whether to record lateral distribution for both ions and recoils
                "phonons": False,      # whether to record lateral distribution for phonons
                "ionization": False,      # whether to record lateral distribution for ionization events
            },
            "backscattered atoms distribution": {
                "energy": {
                    "nbins": 40,         # number of bins for energy distribution
                    "limits": (0.0, 50.0), # limits for energy distribution (keV)
                },
                "angle": {
                    "nbins": 40,         # number of bins for angle distribution
                    "limits": (-90.0, 90.0), # limits for angle distribution
                },
            },
            "transmitted atoms distribution": {
                "energy": {
                    "nbins": 40,         # number of bins for energy distribution
                    "limits": (0.0, 50.0), # limits for energy distribution (keV)
                },
                "angle": {
                    "nbins": 40,         # number of bins for angle distribution
                    "limits": (-90.0, 90.0), # limits for angle distribution
                },
            },
        }
    }
    
    return input_params


def init_params():
    """Initialize the simulation parameters.

    Returns:
        (PARAMS_DTYPE): A structured array containing all simulation parameters.
    """

    # Example hardcoded parameters (to be replaced with file input)
    zmin = 0.0              # minimum z coordinate of the target (A)
    zmax = 4000.0           # maximum z coordinate of the target (A)
    pot_model = "ZBL_magic"  # potential model for scattering
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

    input_params = read_input()

    nlhlin_coefs = read_coefs()
    
    recoil_params = select_recoil.setup(input_params)
    scatter_params = scatter.setup(input_params, nlhlin_coefs)
    cm_scatter.setup(input_params["models"]["scattering integrals"]["n_absc"])
    estop_params = estop.setup(input_params)
    geometry_params = geometry.setup(input_params)
    cascade_params = cascade.setup()
    stat_params = statistics.setup(nspec, nbin, limits)  # TODO: use input_params as argument

    PARAMS_DTYPE = np.dtype([
        ("rng_seed", np.uint64),    # 64 bits needed for alignment
        ("stat", stat_params.dtype),
        ("cascade", cascade_params.dtype),
        ("recoil", recoil_params.dtype),
        ("geometry", geometry_params.dtype),
        ("estop", estop_params.dtype),
        ("scatter", scatter_params.dtype),
    ], align=True)

    params = np.recarray(1, dtype=PARAMS_DTYPE)[0]
    params["stat"] = stat_params
    params["cascade"] = cascade_params
    params["recoil"] = recoil_params
    params["geometry"] = geometry_params
    params["estop"] = estop_params
    params["scatter"] = scatter_params

    return params