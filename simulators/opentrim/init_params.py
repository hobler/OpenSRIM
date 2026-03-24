"""Initialize the simulation parameters except the output parameters.

Available functions:
    get_params: Initialize parameters from the input parameters dictionary.
"""
import os
from collections import namedtuple
from math import sqrt
import numpy as np
#from . import cm_scatter
from .nlhlin import read_coefs


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


def _get_beam_params(input_params):
    """Get the beam parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.

    Returns:
        beam_params: (np.recarray) The beam parameters.
    """
    BEAM_PARAMS_DTYPE = np.dtype([
        ("ielem", np.int32),
        ("energy", np.float64),
        ("tilt", np.float64),
    ], align=True)

    beam_params = np.recarray(1, dtype=BEAM_PARAMS_DTYPE)
    beam_params[0].ielem = 0  # only one beam element for now
    beam_params[0].energy = 1000.0 * input_params["beam"]["energy"]
    beam_params[0].tilt = input_params["beam"]["tilt"]

    return beam_params


def _get_geometry_params(input_params):
    """Get the geometry parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.
    
    Returns:    
        geometry_params: (np.recarray) The geometry parameters.
    """
    # Initialize the geometry
    layers_params = input_params["layer"]
    nlayers = len(layers_params)
    NLAYERS = nlayers   # could be set to max(3, nlayers) to avoid frequent
                        # recompilation of Numba functions when nlayers changes

    x_intf = np.empty(NLAYERS + 1, dtype=np.float64)
    x_intf[0] = 0.0
    for i in range(nlayers):
        x_intf[i+1] = x_intf[i] + layers_params[i]["width"]

    GEOMETRY_PARAMS_DTYPE = np.dtype([
        ("nlayers", np.int32),
        ("x_intf", np.float64, (NLAYERS+1,)),
    ], align=True)

    geometry_params = np.recarray(1, dtype=GEOMETRY_PARAMS_DTYPE)
    geometry_params[0].nlayers = nlayers
    geometry_params[0].x_intf = x_intf

    return geometry_params


def _get_elements_and_materials_params(input_params):
    """Get the elements and materials parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.

    Returns:
        elements_params: (np.recarray) The elements parameters.
        materials_params: (np.recarray) The materials parameters.
    """
    global NELEM, NMAT
    
    # Calculate the atomic fractions for each material from the stoichiometry
    materials = input_params["layer"]
    for mat in materials:
        atomic_fractions = []
        for element in mat["element"]:
            atomic_fractions.append(element["stoichiometry"])
            del element["stoichiometry"]
        atomic_fractions = np.array(atomic_fractions, dtype=np.float64)
        atomic_fractions /= np.sum(atomic_fractions)
        mat["atomic_fractions"] = atomic_fractions
        mat["nelem"] = len(mat["element"])
    nmat = len(materials)
    NMAT = nmat     # could be set to max(3, nmat) to avoid frequent 
                    # recompilation of Numba functions when nmat changes

    # Determine the number of distinct chemical elements in the simulation
    # and replace the element information in the materials with indices of the 
    # elements in the elements list.
    element = {key: input_params["beam"][key] 
               for key in ["symbol", "name", "Z", "M"]}
    element["displacement_energy"] = 0.0
    elements = [element]

    for mat in materials:
        ielem = []
        for elem in mat["element"]:
            if elem in elements:
                ielem.append(elements.index(elem))
            else:
                ielem.append(len(elements))
                elements.append(elem)
        mat["ielem"] = ielem
    nelem = len(elements)
    NELEM = nelem   # could be set to max(5, nelem) to avoid frequent 
                    # recompilation of Numba functions when nelem changes

    # Test:
    if False:
        print(f"Number of distinct chemical elements: {nelem}")
        print("Elements:")
        for elem in elements:
            print(elem)
        print("Materials:")
        for mat in materials:
            print(mat)

    ### Define the elements parameters
    ELEMENT_PARAMS_DTYPE = np.dtype([
        ("symbol", "<U2"),
        ("name", "<U12"),
        ("Z", np.int32),
        ("M", np.float64),
    ], align=True)

    elements_params = np.recarray(NELEM, dtype=ELEMENT_PARAMS_DTYPE)
    for ielem, elem in enumerate(elements):
        elements_params[ielem].symbol = elem["symbol"]
        elements_params[ielem].name = elem["name"]
        elements_params[ielem].Z = elem["Z"]
        elements_params[ielem].M = elem["M"]
#    (f"elements_params={elements_params}")

    ### Define the materials parameters
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

    materials_params = np.recarray(NMAT, dtype=MATERIALS_PARAMS_DTYPE)
    for imat, mat in enumerate(materials):
        materials_params[imat].name = mat["name"]
        materials_params[imat].density = mat["density"]
        materials_params[imat].compound_correction = mat["compound correction"]
        materials_params[imat].gas = mat["gas"]
        materials_params[imat].nelem = mat["nelem"]
        for ielem in range(mat["nelem"]):
            materials_params[imat].ielem[ielem] = mat["ielem"][ielem]
            materials_params[imat].atomic_fraction[ielem] = (
                mat["atomic_fractions"][ielem])
            materials_params[imat].cumulative_fraction[ielem] = (
                np.sum(mat["atomic_fractions"][:ielem+1]))
            materials_params[imat].displacement_energy[ielem] = (
                mat["element"][ielem]["displacement energy"])
#    print(f"materials_params={materials_params}")

    return nelem, elements_params, nmat, materials_params


def _get_estop_params(input_params, nelem, elements_params):
    """Get the electronic stopping parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.
        nelem: (int) The number of distinct chemical elements.
        elements_params: (np.recarray) The elements parameters.

    Returns:
        estop_params: (np.recarray) The electronic stopping parameters.
    """
    ### Define electronic stopping parameters
    model = input_params["models"]["electronic stopping"]
    if model != "Lindhard":
        raise ValueError(f"Unsupported electronic stopping model: {model}")

    # Correction factors to Lindhard stopping power
    corr_lindhard = np.ones((NELEM, NELEM))  # default = 1.0
    for corr in input_params["models"]["Lindhard correction"]:
        elem1, elem2 = corr.split("->")
        for ielem1 in range(nelem):
            for ielem2 in range(nelem):
                if (elements_params[ielem1].symbol.strip() == elem1 and 
                    elements_params[ielem2].symbol.strip() == elem2):
                    corr_lindhard[ielem1, ielem2] = (
                        input_params["models"]["Lindhard correction"][corr])

    # Prefactor for Lindhard stopping (sqrt(eV)*A^2)
    fac_lindhard = np.empty((NELEM, NELEM))
    for ielem1 in range(nelem):
        z1 = elements_params[ielem1].Z
        m1 = elements_params[ielem1].M
        for ielem2 in range(nelem):
            z2 = elements_params[ielem2].Z
            fac_lindhard[ielem1, ielem2] = (
                corr_lindhard[ielem1, ielem2] * 1.212 * z1**(7/6) * z2 
                / ((z1**(2/3) + z2**(2/3))**(3/2) * sqrt(m1)) )

    ESTOP_PARAMS_DTYPE = np.dtype([
        ("fac_lindhard", np.float64, (NELEM, NELEM)),
    ], align=True)

    estop_params = np.recarray(1, dtype=ESTOP_PARAMS_DTYPE)
    estop_params[0].fac_lindhard = fac_lindhard

    return estop_params


def _get_recoil_params(input_params):
    """Get the recoil parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.

    Returns:
        recoil_params: (np.recarray) The recoil parameters.
    """
    densities = np.array([layer["density"] for layer in input_params["layer"]])

    RECOIL_PARAMS_DTYPE = np.dtype([
        ("pmax", np.float64, (NMAT,)),
        ("mean_free_path", np.float64, (NMAT)),
    ], align=True)

    recoil_params = np.recarray(1, dtype=RECOIL_PARAMS_DTYPE)
    recoil_params[0].pmax = densities**(-1/3) / sqrt(np.pi)
    recoil_params[0].mean_free_path = densities**(-1/3)

    return recoil_params


def _get_scatter_params(input_params, nelem, elements_params):
    """Get the scattering parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.
        nelem: (int) The number of distinct chemical elements.
        elements_params: (np.recarray) The elements parameters.

    Returns:
        scatter_params: (np.recarray) The scattering parameters.
    """
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
            z1 = elements_params[ielem1].Z
            z2 = elements_params[ielem2].Z
            m1 = elements_params[ielem1].M
            m2 = elements_params[ielem2].M
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

    scatter_params = np.recarray(1, dtype=SCATTER_PARAMS_DTYPE)
    scatter_params[0].pot_model = pot_model
    scatter_params[0].enorm = enorm
    scatter_params[0].rnorm = rnorm
    scatter_params[0].dirfac = dirfac
    scatter_params[0].denfac = denfac
    if pot_model == "NLHlin":
        scatter_params.nlhlin_coefs = nlhlin_coefs
    
    return scatter_params


def _get_cascade_params(input_params):
    """Get the cascade parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.

    Returns:
        cascade_params: (np.recarray) The cascade parameters.
    """
    CASCADE_PARAMS_DTYPE = np.dtype([
        ("follow_recoils", np.int64),  # stored as int for better compatibility with Numba
        ("emin", np.float64),
        ("ed", np.float64),
    ], align=True)

    cascade_params = np.recarray(1, dtype=CASCADE_PARAMS_DTYPE)
    cascade_params[0].follow_recoils = (
        input_params["simulation"]["follow recoils"])
    cascade_params[0].emin = 5.0
    cascade_params[0].ed = 15.0

    return cascade_params


def get_params(input_params):
    """Initialize the simulation parameters except the output parameters.
    
    The parameters are returned as a structured array of subarrays. Do not
    create a structured scalar by saying params = params[0] and similarly for
    the subarrays, since this would cause issues with parallelization in Numba.

    Parameters:
        input_params: (dict) The input parameters dictionary.

    Returns:
        params: (np.recarray) The input parameters structured array.
    """
    global NELEM, NMAT
    
    # rearrange parameters and calculate derived parameters
    beam_params = _get_beam_params(input_params)
    geometry_params = _get_geometry_params(input_params)
    nelem, elements_params, nmat, materials_params = (
        _get_elements_and_materials_params(input_params))
    estop_params = _get_estop_params(input_params, nelem, elements_params)
    recoil_params = _get_recoil_params(input_params)
    scatter_params = _get_scatter_params(input_params, nelem, elements_params)
    cascade_params = _get_cascade_params(input_params)

    # TODO: include n_absc in params
    #cm_scatter.setup(input_params["models"]["scattering integrals"]["n_absc"])

    # collect subarrays into a single structured array
    PARAMS_DTYPE = np.dtype([
        ("rng_seed", np.int64),     # need an even number of int32 for alignment
        ("nelem", np.int32),        # without padding, which causes trouble with Numba when
        ("nmat", np.int32),         # align=True
        ("beam", beam_params.dtype),
        ("cascade", cascade_params.dtype),
        ("recoil", recoil_params.dtype),
        ("geometry", geometry_params.dtype),
        ("elements", elements_params.dtype, (elements_params.size,)),
        ("materials", materials_params.dtype, (materials_params.size,)),
        ("estop", estop_params.dtype),
        ("scatter", scatter_params.dtype),
    ], align=True)

    params = np.recarray(1, dtype=PARAMS_DTYPE)
    # NOTE: Do not do "params = params[0]", since this would create a structured
    # scalar, which causes issues with parallelization in Numba.
    params[0].rng_seed = input_params["simulation"]["rng_seed"]
    params[0].nelem = nelem
    params[0].nmat = nmat
    params[0].beam = beam_params
    params[0].cascade = cascade_params
    params[0].recoil = recoil_params
    params[0].geometry = geometry_params
    params[0].elements = elements_params
    params[0].materials = materials_params
    params[0].estop = estop_params
    params[0].scatter = scatter_params

    return params