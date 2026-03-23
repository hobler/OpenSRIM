def read_params():
    """Read the input parameters from a file.

    Currently hardcoded, but should read from a TOML file in the future.

    Returns:
        (dict): A dictionary containing the input parameters.
    """
    input_params = {
        "simulation": {
            "follow recoils": True,  # whether to follow recoils in the simulation
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
                },
            ],
        },
#
# The layers part of the TOML file would look something like this:
#
# [layers]
# name = ["Layer 1", ...]
# width = [4000.0, ...]
# density = [0.04994, ...]
# compound_correction = [1.0, ...]
# gas = [false, ...] 
#
# [[layers.material]]
# symbol = ["Si", ...]
# name = ["Silicon", ...]
# Z = [14, ...]
# M = [28.086, ...]
# stoichiometry = [1, ...]
# displacement_energy = [15.0, ...]
#
# [[layers.material]]
# ...
#
# TODO: 
#       "layer": [
#           {
#               "name": "Layer 1",
#               "width": 4000.0,
#               "density": 0.04994,
#               "compound correction": 1.0,
#               "gas": False,
#               "element": [
#                   {
#                       "symbol": "Si",
#                       "name": "Silicon",
#                       "Z": 14,
#                       "M": 28.086,
#                       "stoichiometry": 1,
#                       "displacement energy": 15.0,
#                   },
#                   ...
#               ],
#           },
#           ...
#       ],
#
# would maybe result in better readability of the TOML file:
#
# [[layer]]
# name = "Layer 1"
# width = 4000.0
# density = 0.04994
# compound_correction = 1.0
# gas = false
#
# [[layer.element]]
# symbol = "Si"
# name = "Silicon"
# Z = 14
# M = 28.086
# stoichiometry = 1
# displacement_energy = 15.0
#
# [[layer.element]]
# ...
#
# [[layer]]
# ...
#
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
                "collisions": False,  # whether to record collision points of trajectories
                "stopped": False,        # whether to record ending points of trajectories
                "backscattered": False, # whether to record backscattered projectiles
                "transmitted": False,   # whether to record transmitted projectiles
            },
            "depth distribution": {
                "ion/recoils": {
                    "score": True,       # whether to score depth distribution for both ions and recoils
                    "nbins": 40,         # number of bins for depth distribution
                    "limits": (0.0, 4000.0), # limits for depth distribution (A)
                },
                "nuclear energy deposition": {
                    "score": False,      # whether to record depth distribution for NED events
                    "nbins": 40,         # number of bins for depth distribution
                    "limits": (0.0, 4000.0), # limits for depth distribution (A)
                },
                "electronic energy deposition": {
                    "score": False,      # whether to record depth distribution for EED events
                    "nbins": 40,         # number of bins for depth distribution
                    "limits": (0.0, 4000.0), # limits for depth distribution (A)
                }
            },
            "lateral distribution": {
                "ion/recoils": {
                    "score": True,       # whether to score lateral distribution for both ions and recoils
                    "nbins": 40,         # number of bins for lateral distribution
                    "limits": (-2000.0, 2000.0), # limits for lateral distribution (A)
                },
                "nuclear energy deposition": {
                    "score": False,      # whether to record lateral distribution for NED events
                    "nbins": 40,         # number of bins for lateral distribution
                    "limits": (-2000.0, 2000.0), # limits for lateral distribution (A)
                },
                "electronic energy deposition": {
                    "score": False,      # whether to record lateral distribution for EED events
                    "nbins": 40,         # number of bins for lateral distribution
                    "limits": (-2000.0, 2000.0), # limits for lateral distribution (A)
                }
            },
            "backscattered atoms distribution": {
                "energy": {
                    "score": True,       # whether to score energy distribution for backscattered atoms
                    "nbins": 40,         # number of bins for energy distribution
                    "limits": (0.0, 50.0), # limits for energy distribution (keV)
                },
                "angle": {
                    "score": True,       # whether to score angle distribution for backscattered atoms
                    "nbins": 40,         # number of bins for angle distribution
                    "limits": (-90.0, 90.0), # limits for angle distribution
                },
            },
            "transmitted atoms distribution": {
                "energy": {
                    "score": True,       # whether to score energy distribution for transmitted atoms
                    "nbins": 40,         # number of bins for energy distribution
                    "limits": (0.0, 50.0), # limits for energy distribution (keV)
                },
                "angle": {
                    "score": True,       # whether to score angle distribution for transmitted atoms
                    "nbins": 40,         # number of bins for angle distribution
                    "limits": (-90.0, 90.0), # limits for angle distribution
                },
            },
        }
    }
    
    return input_params
