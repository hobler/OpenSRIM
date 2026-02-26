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

NLAYERS = 5
#NLAYERS = nlayers
z_intf = np.empty(NLAYERS + 1, dtype=np.float64)
z_intf[0] = 0.0
for i in range(nlayers):
    z_intf[i+1] = z_intf[i] + layers_params["width"][i]

GEOMETRY_PARAMS_DTYPE = np.dtype([
    ("nlayers", np.int32),
    ("z_intf", np.float64, (NLAYERS+1,)),
], align=True)

geometry_params = np.recarray(1, dtype=GEOMETRY_PARAMS_DTYPE)[0]
geometry_params.nlayers = nlayers
geometry_params.z_intf = z_intf

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
ELEMENT_PARAMS_DTYPE = np.dtype([
    ("symbol", "<U2"),
    ("name", "<U12"),
    ("Z", np.int32),
    ("M", np.float64),
], align=True)

NELEM = 5
#NELEM = nelem

elements_params = np.recarray(NELEM, dtype=ELEMENT_PARAMS_DTYPE)
for ielem, elem in enumerate(elements):
    elements_params.symbol[ielem] = elem["symbol"]
    elements_params.name[ielem] = elem["name"]
    elements_params.Z[ielem] = elem["Z"]
    elements_params.M[ielem] = elem["M"]
print(f"elements_params={elements_params}")

# Define the materials parameters

MATERIALS_PARAMS_DTYPE = np.dtype([
    ("name", "<U12"),
    ("density", np.float64),
    ("compound_correction", np.float64),
    ("gas", bool),
    ("nelem", np.int32),
    ("ielem", np.int32, (NELEM,)),
    ("atomic_fraction", np.float64, (NELEM,)),
    ("cumulative_fraction", np.float64, (NELEM,)),
    ("displacement_energy", np.float64, (NELEM,)),
], align=True)

NMAT = 5
#NMAT = nmat

materials_params = np.recarray(NMAT, dtype=MATERIALS_PARAMS_DTYPE)
for imat, mat in enumerate(materials):
    materials_params.name[imat] = mat["name"]
    materials_params.density[imat] = mat["density"]
    materials_params.compound_correction[imat] = mat["compound_correction"]
    materials_params.gas[imat] = mat["gas"]
    materials_params.nelem[imat] = mat["nelem"]
    for ielem in range(mat["nelem"]):
        materials_params.ielem[imat, ielem] = mat["ielem"][ielem]
        materials_params.atomic_fraction[imat, ielem] = mat["atomic_fractions"][ielem]
        materials_params.cumulative_fraction[imat, ielem] = (
            np.sum(mat["atomic_fractions"][:ielem+1]))
        materials_params.displacement_energy[imat, ielem] = (
            mat["displacement_energy"][ielem])
print(f"materials_params={materials_params}")

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


ESTOP_PARAMS_DTYPE = np.dtype([
    ("fac_lindhard", np.float64, (NELEM, NELEM)),
], align=True)

estop_params = np.recarray(1, dtype=ESTOP_PARAMS_DTYPE)[0]
estop_params.fac_lindhard = fac_lindhard

# Initialize the recoiling parameters
densities = np.array(input_params["layers"]["density"])

RECOIL_PARAMS_DTYPE = np.dtype([
    ("pmax", np.float64, (NMAT,)),
    ("mean_free_path", np.float64, (NMAT)),
], align=True)

recoil_params = np.recarray(1, dtype=RECOIL_PARAMS_DTYPE)[0]
recoil_params.pmax = densities**(-1/3) / sqrt(np.pi)
recoil_params.mean_free_path = densities**(-1/3)

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

SCATTER_PARAMS_DTYPE = np.dtype([
    ("pot_model", "<U12"),
    ("enorm", np.float64, (NELEM, NELEM)),
    ("rnorm", np.float64, (NELEM, NELEM)),
    ("dirfac", np.float64, (NELEM, NELEM)),
    ("denfac", np.float64, (NELEM, NELEM)),
], align=True)

if pot_model == "NLHlin":
    SCATTER_PARAMS_DTYPE = np.dtype(SCATTER_PARAMS_DTYPE.descr + [
        ("nlhlin_coefs", NLHLIN_COEFS_DTYPE, (NELEM, NELEM)),
    ], align=True)

scatter_params = np.recarray(1, dtype=SCATTER_PARAMS_DTYPE)[0]
scatter_params.pot_model = pot_model
scatter_params.enorm = enorm
scatter_params.rnorm = rnorm
scatter_params.dirfac = dirfac
scatter_params.denfac = denfac
if pot_model == "NLHlin":
    scatter_params.nlhlin_coefs = nlhlin_coefs

# TODO: include n_absc in params
#cm_scatter.setup(input_params["models"]["scattering integrals"]["n_absc"])

# Initialize the cascade parameters
CASCADE_PARAMS_DTYPE = np.dtype([
    ("emin", np.float64),
    ("ed", np.float64),
], align=True)
cascade_params = np.recarray(1, dtype=CASCADE_PARAMS_DTYPE)[0]
cascade_params.emin = 5.0
cascade_params.ed = 15.0

# Example hardcoded parameters (TODO: to be replaced with file input)
#nelem = len(elements_params)  # number of species to record
nbin = 40              # number of bins for depth distribution
limits = (0.0, 4000.0) # limits for depth distribution

STAT_PARAMS_DTYPE = np.dtype([
    ("nspec", np.int32),
    ("nbin", np.int32),
    ("limits", np.float64, (2,)),
], align=True)

stat_params = np.recarray(1, dtype=STAT_PARAMS_DTYPE)[0]
stat_params.nspec = NELEM
stat_params.nbin = nbin
stat_params.limits = limits

#statistics.setup(nelem, nbin, limits)  # TODO: use input_params as argument
if False:
    PARAMS_DTYPE = np.dtype([
        ("rng_seed", np.int32),
        ("nelem", np.int32),
        ("nmat", np.int32),
        ("cascade", CASCADE_PARAMS_DTYPE),
        ("recoil", RECOIL_PARAMS_DTYPE),
        ("geometry", GEOMETRY_PARAMS_DTYPE),
        ("elements", ELEMENT_PARAMS_DTYPE, (NELEM,)),
        ("materials", MATERIALS_PARAMS_DTYPE, (NMAT,)),
        ("estop", ESTOP_PARAMS_DTYPE),
        ("scatter", SCATTER_PARAMS_DTYPE),
        ("stat", STAT_PARAMS_DTYPE),
    ], align=True)

    params = np.recarray(1, dtype=PARAMS_DTYPE)[0]
    params.rng_seed = input_params["simulation"]["rng_seed"]
    params.nelem = nelem
    params.nmat = nmat
    params.cascade = cascade_params
    params.recoil = recoil_params
    params.geometry = geometry_params
    params.elements = elements_params
    params.materials = materials_params
    params.estop = estop_params
    params.scatter = scatter_params
    params.stat = stat_params
else:
    Params = namedtuple("Params", [
        "rng_seed",
        "nelem",
        "nmat",
        "cascade",
        "recoil",
        "geometry",
        "elements",
        "materials",
        "estop",
        "scatter",
        "stat",
    ])
    params = Params(
        rng_seed=input_params["simulation"]["rng_seed"],
        nelem=nelem,
        nmat=nmat,
        cascade=cascade_params,
        recoil=recoil_params,
        geometry=geometry_params,
        elements=elements_params,
        materials=materials_params,
        estop=estop_params,
        scatter=scatter_params,
        stat=stat_params,
    )