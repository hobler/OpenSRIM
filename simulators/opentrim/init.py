"""Initialize the simulation parameters.

The params named tuple is defined during import of this module, and contains 
all the parameters needed for the simulation. This code must not be put into a 
function, since the params named tuple needs to be defined at the module level 
for Numba compatibility.

- Read the input parameters from a TOML file (not implemented yet, currently 
  hardcoded).
- Calculate derived parameters.
"""
import os
from collections import namedtuple
from math import sqrt
import numpy as np
#from . import cm_scatter
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


def _get_nlhlin_coefs(z1, z2):
    """Get the coefficients for the NLHlin screening function.

    Parameters:
        z1: (int) atomic number of atom 1
        z2: (int) atomic number of atom 2
        
    Returns:
        a1, a2, a3: prefactors
        b1, b2, b3: 1/screening lengths (1/A)
        rmax: maximum range of the potential (A)
    """
    fname = os.path.join(os.path.dirname(__file__), "dmol_coeffs_rmax.dat")
    if not os.path.exists(fname):
        raise OSError(f"get_nlhlin_coefs: Coefficients file {fname} not found")
    
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            items = line.split()
            if int(items[0]) == min(z1, z2) and int(items[1]) == max(z1, z2):
                return [float(c) for c in items[2:-1]]
    
    raise ValueError(f"get_nlhlin_coefs: Coefficients for Z1={z1}, Z2={z2} "
                     f"not found in {fname}")


input_params = read_input()

# Initialize the geometry
layers_params = input_params["layers"]
nlayers = len(layers_params["name"])

#NLAYERS = 5
NLAYERS = nlayers
z_intf = np.empty(NLAYERS + 1, dtype=np.float64)
z_intf[0] = 0.0
for i in range(nlayers):
    z_intf[i+1] = z_intf[i] + layers_params["width"][i]

GeometryParams = namedtuple("GeometryParams", ["nlayers", "z_intf"])

geometry_params = GeometryParams(
    nlayers = nlayers,
    z_intf = z_intf,
)

# Extract information about the target materials and construct a more 
# convenient data structure.
# In the new data structure, each material has the format
# {
#     "name": str,
#     "density": float,
#     "compound_correction": float,
#     "gas": bool,
#     "elements": [
#         {
#             "symbol": str,
#             "name": str,
#             "Z": int,
#             "M": float,
#         },
#         ...
#     ],
#     "atomic_fractions": [float, ...],
#     "displacement_energies": [float, ...],
# }
# TODO: maybe the UI can already provide the materials in this format
# But this would mean that the input file format would be more complex, so 
# maybe it's better to keep it simple and do the conversion in code for now.

materials = []
for ilayer in range(nlayers):
    mat_elements = []
    atomic_fractions = []
    displacement_energies = []
    mat = layers_params["material"][ilayer]
    for ielem in range(len(mat["symbol"])):
        element = {key: mat[key][ielem] 
                    for key in ["symbol", "name", "Z", "M", 
                                "displacement_energy"]}
        mat_elements.append(element)
        atomic_fraction = (mat["stoichiometry"][ielem] 
                            / sum(mat["stoichiometry"]))
        atomic_fractions.append(atomic_fraction)
        displacement_energies.append(mat["displacement_energy"][ielem])
    
    material = {
        "name": layers_params["name"][ilayer],
        "density": layers_params["density"][ilayer],
        "compound_correction": layers_params["compound correction"][ilayer],
        "gas": layers_params["gas"][ilayer],
        "nelem": len(mat_elements),
        "elements": mat_elements,
        "atomic_fractions": atomic_fractions,
        "displacement_energy": displacement_energies,
    }
    materials.append(material)
nmat = len(materials)

# Determine the number of distinct chemical elements in the simulation
# and replace the element information in the materials with indices of the 
# elements in the elements list.
element = {key: input_params["beam"][key] 
            for key in ["symbol", "name", "Z", "M"]}
element["displacement_energy"] = 0.0
elements = [element]

for mat in materials:
    ielem = []
    for elem in mat["elements"]:
        if elem in elements:
            ielem.append(elements.index(elem))
        else:
            ielem.append(len(elements))
            elements.append(elem)
    mat["ielem"] = ielem
    del mat["elements"]
nelem = len(elements)

# Test:
print(f"Number of distinct chemical elements: {nelem}")
print("Elements:")
for elem in elements:
    print(elem)
print("Materials:")
for mat in materials:
    print(mat)

# Define the elements parameters
ElementsParams = namedtuple("ElementParams", ["symbol", "name", "Z", "M"])

NELEM = 5
#NELEM = nelem
elements_params = ElementsParams(
    symbol = np.array([elem["symbol"] for elem in elements], dtype="<U2"),
    name = np.array([elem["name"] for elem in elements], dtype="<U12"),
    Z = np.array([elem["Z"] for elem in elements], dtype=np.int32),
    M = np.array([elem["M"] for elem in elements], dtype=np.float64),
)
print(f"elements_params={elements_params}")

# Define the materials parameters
NMAT = 5
#NMAT = nmat
mat_names = np.empty(NMAT, dtype="<U12")
mat_densities = np.empty(NMAT, dtype=np.float64)
mat_compound_corrections = np.ones(NMAT, dtype=np.float64)
mat_gas = np.empty(NMAT, dtype=bool)
mat_nelem = np.empty(NMAT, dtype=np.int32)
mat_ielem = np.empty((NMAT, NELEM), dtype=np.int32)
mat_atomic_fractions = np.empty((NMAT, NELEM), dtype=np.float64)
mat_cumulative_fractions = np.empty((NMAT, NELEM), dtype=np.float64)
mat_displacement_energies = np.empty((NMAT, NELEM), dtype=np.float64)

for imat, mat in enumerate(materials):
    mat_names[imat] = mat["name"]
    mat_densities[imat] = mat["density"]
    mat_compound_corrections[imat] = mat["compound_correction"]
    mat_gas[imat] = mat["gas"]
    mat_nelem[imat] = mat["nelem"]
    for ielem in range(mat_nelem[imat]):
        mat_ielem[imat, ielem] = mat["ielem"][ielem]
        mat_atomic_fractions[imat, ielem] = mat["atomic_fractions"][ielem]
        mat_cumulative_fractions[imat, ielem] = (
            np.sum(mat["atomic_fractions"][:ielem+1]))
        mat_displacement_energies[imat, ielem] = (
            mat["displacement_energy"][ielem])

MaterialsParams = namedtuple("MaterialsParams", [
    "name", "density", "compound_correction", "gas", "nelem",
    "ielem", "atomic_fraction", "cumulative_fraction", "displacement_energy",
])

materials_params = MaterialsParams(
    name = mat_names,
    density = mat_densities,
    compound_correction = mat_compound_corrections,
    gas = mat_gas,
    nelem = mat_nelem,
    ielem = mat_ielem,
    atomic_fraction = mat_atomic_fractions,
    cumulative_fraction = mat_cumulative_fractions,
    displacement_energy = mat_displacement_energies,
)

# Define electronic stopping parameters
model = input_params["models"]["electronic stopping"]
if model != "Lindhard":
    raise ValueError(f"Unsupported electronic stopping model: {model}")

# Correction factors to Lindhard stopping power
corr_lindhard = np.ones((NELEM, NELEM))  # default = 1.0
for corr in input_params["models"]["Lindhard correction"]:
    elem1, elem2 = corr.split("->")
    for ielem1 in range(nelem):
        for ielem2 in range(nelem):
            if (elements_params.symbol[ielem1].strip() == elem1 and 
                elements_params.symbol[ielem2].strip() == elem2):
                corr_lindhard[ielem1, ielem2] = (
                    input_params["models"]["Lindhard correction"][corr])

# Prefactor for Lindhard stopping (sqrt(eV)*A^2)
fac_lindhard = np.empty((NELEM, NELEM))
for ielem1 in range(nelem):
    z1 = elements_params.Z[ielem1]
    m1 = elements_params.M[ielem1]
    for ielem2 in range(nelem):
        z2 = elements_params.Z[ielem2]
        fac_lindhard[ielem1, ielem2] = (
            corr_lindhard[ielem1, ielem2] * 1.212 * z1**(7/6) * z2 
            / ((z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1)) )

EstopParams = namedtuple("EstopParams", ["fac_lindhard"])
estop_params = EstopParams(
    fac_lindhard = fac_lindhard,
)

# Initialize the recoiling parameters
densities = np.array(input_params["layers"]["density"])

RecoilParams = namedtuple("RecoilParams", ["pmax", "mean_free_path"])
recoil_params = RecoilParams(
    pmax = densities**(-1/3) / sqrt(np.pi),
    mean_free_path = densities**(-1/3),
)

# Initialize the scattering parameters
pot_model = input_params["models"]["potential"]
if input_params["models"]["scattering integrals"]["algorithm"] == "magic":
    pot_model += "_magic"

rnorm = np.empty((NELEM, NELEM), dtype=np.float64)
enorm = np.empty((NELEM, NELEM), dtype=np.float64)
dirfac = np.empty((NELEM, NELEM), dtype=np.float64)
denfac = np.empty((NELEM, NELEM), dtype=np.float64)

if pot_model == "NLHlin":
    NLHLIN_COEFS_DTYPE = np.dtype([
        ("a1", np.float64),
        ("b1", np.float64),
        ("a2", np.float64),
        ("b2", np.float64),
        ("a3", np.float64),
        ("b3", np.float64),
        ("c", np.float64),
        ("d", np.float64),
        ("rmax", np.float64),
    ], align=True)
    nlhlin_coefs = np.empty((NELEM, NELEM), dtype=NLHLIN_COEFS_DTYPE)

for ielem1 in range(nelem):
    for ielem2 in range(nelem):
        z1 = elements_params.Z[ielem1]
        z2 = elements_params.Z[ielem2]
        m1 = elements_params.M[ielem1]
        m2 = elements_params.M[ielem2]
        m1_m2 = m1 / m2

        if pot_model.startswith("ZBL"):
            rnorm[ielem1, ielem2] = 0.4685 / (z1**0.23 + z2**0.23)
        elif pot_model.startswith("NLHlin"):
            rnorm[ielem1, ielem2] = 0.4685 / (sqrt(sqrt(z1) + sqrt(z2)))
        enorm[ielem1, ielem2] = (14.39979 * z1 * z2 / rnorm[ielem1, ielem2] 
                                 * (1 + m1_m2))
        dirfac[ielem1, ielem2] = 2 / (1 + m1_m2)
        denfac[ielem1, ielem2] = 4 * m1_m2 / (1 + m1_m2)**2
        
        if pot_model == "NLHlin":
            a1, b1, a2, b2, a3, b3, rmax = _get_nlhlin_coefs(z1, z2)
            b1 *= rnorm[ielem1, ielem2]
            b2 *= rnorm[ielem1, ielem2]
            b3 *= rnorm[ielem1, ielem2]
            c = 1.0 - (a1 + a2 + a3)
            rmax /= rnorm[ielem1, ielem2]
            d = (a1*np.exp(-b1*rmax) + a2*np.exp(-b2*rmax) 
                 + a3*np.exp(-b3*rmax) + c)
            nlhlin_coefs[ielem1, ielem2]["a1"] = a1
            nlhlin_coefs[ielem1, ielem2]["b1"] = b1
            nlhlin_coefs[ielem1, ielem2]["a2"] = a2
            nlhlin_coefs[ielem1, ielem2]["b2"] = b2
            nlhlin_coefs[ielem1, ielem2]["a3"] = a3
            nlhlin_coefs[ielem1, ielem2]["b3"] = b3
            nlhlin_coefs[ielem1, ielem2]["c"] = c
            nlhlin_coefs[ielem1, ielem2]["d"] = d
            nlhlin_coefs[ielem1, ielem2]["rmax"] = rmax 

ScatterParams = namedtuple("ScatterParams", [
    "pot_model", "enorm", "rnorm", "dirfac", "denfac", "nlhlin_coefs",
])
scatter_params = ScatterParams(
    pot_model = pot_model,
    enorm = enorm,
    rnorm = rnorm,
    dirfac = dirfac,
    denfac = denfac,
    nlhlin_coefs = nlhlin_coefs if pot_model == "NLHlin" else None,
)

# TODO: include n_absc in params
#cm_scatter.setup(input_params["models"]["scattering integrals"]["n_absc"])

# Initialize the cascade parameters
CascadeParams = namedtuple("CascadeParams", ["emin", "ed"])
cascade_params = CascadeParams(
    emin = 5.0,  # eV
    ed = 15.0,   # eV
)

# Example hardcoded parameters (to be replaced with file input)
#nelem = len(elements_params)  # number of species to record
nbin = 40              # number of bins for depth distribution
limits = (0.0, 4000.0) # limits for depth distribution

StatParams = namedtuple("StatParams", ["nspec", "nbin", "limits"])
stat_params = StatParams(
    nspec = NELEM,
    nbin = nbin,
    limits = limits,
)

#statistics.setup(nelem, nbin, limits)  # TODO: use input_params as argument

Params = namedtuple("Params", ["rng_seed", "nelem", "nmat",
                               "cascade", "recoil", "geometry", 
                               "elements", "materials", "estop", "scatter", 
                               "stat", 
                               ])
params = Params(
    rng_seed = input_params["simulation"]["rng_seed"],
    nelem = nelem,
    nmat = nmat,
    recoil = recoil_params,
    geometry = geometry_params,
    elements = elements_params,
    materials = materials_params,
    estop = estop_params,
    scatter = scatter_params,
    cascade = cascade_params,
    stat = stat_params,
)
