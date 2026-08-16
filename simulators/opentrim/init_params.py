"""Initialize the simulation parameters except the output parameters.

Available functions:
    get_params: Initialize parameters from the input parameters dictionary.
"""
import os
from collections import namedtuple
from math import sqrt
import numpy as np
#from . import cm_scatter
from . import nlhlin
from . import zbl


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
    global NELEM, NELEM_ION, NELEM_TARGET, NMAT
    
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
    nelem_ion = 1  # only one ion element for now
    NELEM_ION = nelem_ion

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
    nelem_target = nelem - nelem_ion
    NELEM_TARGET = nelem_target

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
        ("edisp", np.float64, (NELEM,)),
        ("esurf", np.float64, (NELEM,)),
        ("ebulk", np.float64, (NELEM,)),
    ], align=True)

    materials_params = np.recarray(NMAT, dtype=MATERIALS_PARAMS_DTYPE)
    for imat, mat in enumerate(materials):
        materials_params[imat].name = mat["name"]
        materials_params[imat].density = mat["density"]
        materials_params[imat].compound_correction = mat["compound_correction"]
        materials_params[imat].gas = mat["gas"]
        materials_params[imat].nelem = mat["nelem"]
        for ielem in range(mat["nelem"]):
            materials_params[imat].ielem[ielem] = mat["ielem"][ielem]
            materials_params[imat].atomic_fraction[ielem] = (
                mat["atomic_fractions"][ielem])
            materials_params[imat].cumulative_fraction[ielem] = (
                np.sum(mat["atomic_fractions"][:ielem+1]))
            materials_params[imat].edisp[ielem] = (
                mat["element"][ielem]["displacement_energy"])
            materials_params[imat].esurf[ielem] = (
                mat["element"][ielem]["surface_binding_energy"])
            materials_params[imat].ebulk[ielem] = (
                mat["element"][ielem]["lattice_binding_energy"])
#    print(f"materials_params={materials_params}")

    return nelem_target, nelem, elements_params, materials_params


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
    model = input_params["models"]["electronic_stopping"]

    if model == "Lindhard":
        # Correction factors to Lindhard stopping power
        corr_lindhard = np.ones((NELEM, NELEM))  # default = 1.0
        for corr in input_params["models"]["lindhard_correction"]:
            elem1, elem2 = corr.split("->")
            for ielem1 in range(nelem):
                for ielem2 in range(nelem):
                    if (elements_params[ielem1].symbol.strip() == elem1 and 
                        elements_params[ielem2].symbol.strip() == elem2):
                        corr_lindhard[ielem1, ielem2] = (
                            input_params["models"]["lindhard_correction"][corr])

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

    elif model == "SRIM":
        path = os.path.join(os.path.dirname(__file__), "../../data/SRIM_setab/")
        srim_energies = np.empty(219, dtype=np.float64)
        srim_table = np.empty((NELEM, NELEM, 219), dtype=np.float64)
        for ielem1 in range(nelem):
            z1 = elements_params[ielem1].Z
            filename = os.path.join(path, f'SRIM2013-{z1:02d}.dat')
            srim_energies = np.loadtxt(filename, skiprows=6, usecols=0)
            for ielem2 in range(nelem):
                z2 = elements_params[ielem2].Z
                srim_table[ielem1, ielem2] = np.loadtxt(
                    filename, skiprows=6, usecols=z2)

    else:
        raise ValueError(f"Unknown electronic stopping model: {model}")

    ESTOP_PARAMS_DTYPE = np.dtype([
        ("model", "<U8"),
        ("fac_lindhard", np.float64, (NELEM, NELEM)),
        ("srim_energies", np.float64, 219),  # for possible future use
        ("srim_table", np.float64, (NELEM, NELEM, 219)),  # for possible future use
    ], align=True)

    estop_params = np.recarray(1, dtype=ESTOP_PARAMS_DTYPE)
    estop_params[0].model = model
    if model == "Lindhard":
        estop_params[0].fac_lindhard = fac_lindhard
    elif model == "SRIM":
        estop_params[0].srim_energies = srim_energies
        estop_params[0].srim_table = srim_table

    return estop_params


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
    integrate_algorithm = (
        input_params["models"]["scattering_integrals"]["algorithm"])

    rnorm = np.empty((NELEM, NELEM), dtype=np.float64)
    enorm = np.empty((NELEM, NELEM), dtype=np.float64)
    dirfac = np.empty((NELEM, NELEM), dtype=np.float64)
    denfac = np.empty((NELEM, NELEM), dtype=np.float64)

    POT_COEFS_DTYPE = np.dtype([
        ("a", np.float64, (4,)),
        ("b", np.float64, (4,)),
        ("c", np.float64),
        ("d", np.float64),
        ("rmax", np.float64),
        ("r34", np.float64),
        ("k", np.float64, (4,)),
    ], align=True)
    pot_coefs = np.empty((NELEM, NELEM), dtype=POT_COEFS_DTYPE)

    for ielem1 in range(nelem):
        for ielem2 in range(nelem):
            z1 = elements_params[ielem1].Z
            z2 = elements_params[ielem2].Z
            m1 = elements_params[ielem1].M
            m2 = elements_params[ielem2].M
            m1_m2 = m1 / m2
            dirfac[ielem1, ielem2] = 2 / (1 + m1_m2)
            denfac[ielem1, ielem2] = 4 * m1_m2 / (1 + m1_m2)**2

            if pot_model == "NLHlin":
                rnorm_, a, b, c, d, r34, rmax, k = nlhlin.get_coefs(z1, z2)
                pot_coefs[ielem1, ielem2]["c"] = c
                pot_coefs[ielem1, ielem2]["d"] = d
                pot_coefs[ielem1, ielem2]["rmax"] = rmax
            else:
                rnorm_, a, b, r34, k = zbl.get_coefs(z1, z2)
                pot_coefs[ielem1, ielem2]["c"] = 0.0
                pot_coefs[ielem1, ielem2]["d"] = 0.0
                pot_coefs[ielem1, ielem2]["rmax"] = np.inf
            rnorm[ielem1, ielem2] = rnorm_
            pot_coefs[ielem1, ielem2]["a"][:len(a)] = a
            pot_coefs[ielem1, ielem2]["b"][:len(b)] = b
            pot_coefs[ielem1, ielem2]["r34"] = r34
            pot_coefs[ielem1, ielem2]["k"][:len(k)] = k

            enorm[ielem1, ielem2] = (14.39979 * z1 * z2 / rnorm[ielem1, ielem2] 
                                    * (1 + m1_m2))

    SCATTER_PARAMS_DTYPE = np.dtype([
        ("pot_model", "<U16"),
        ("integrate_algorithm", "<U8"),
        ("enorm", np.float64, (NELEM, NELEM)),
        ("rnorm", np.float64, (NELEM, NELEM)),
        ("dirfac", np.float64, (NELEM, NELEM)),
        ("denfac", np.float64, (NELEM, NELEM)),
        ("pot_coefs", POT_COEFS_DTYPE, (NELEM, NELEM)),
    ], align=True)

    scatter_params = np.recarray(1, dtype=SCATTER_PARAMS_DTYPE)
    scatter_params[0].pot_model = pot_model
    scatter_params[0].integrate_algorithm = integrate_algorithm
    scatter_params[0].enorm = enorm
    scatter_params[0].rnorm = rnorm
    scatter_params[0].dirfac = dirfac
    scatter_params[0].denfac = denfac
    scatter_params[0].pot_coefs = pot_coefs

    return scatter_params


def _get_cascade_params(input_params, nelem, elements_params, materials_params,
                        scatter_params):
    """Get the cascade parameters from the input parameters.

    Parameters:
        input_params: (dict) The input parameters dictionary.

    Returns:
        cascade_params: (np.recarray) The cascade parameters.
    """
    NPMAX = 64
    CASCADE_PARAMS_DTYPE = np.dtype([
        ("follow_recoils", np.int64),  # stored as int for better compatibility with Numba
        ("emin", np.float64),
        ("ed", np.float64),
        ("pmax_max", np.float64),
        ("pmax", np.float64, (NMAT,)),
        ("mean_free_path", np.float64, (NMAT,)),
        ("pmax_vals", np.float64, (NPMAX,)),
        ("pmax_energies", np.float64, (NELEM, NMAT, NPMAX)),
    ], align=True)

    cascade_params = np.recarray(1, dtype=CASCADE_PARAMS_DTYPE)
    cascade_params[0].follow_recoils = (
        input_params["simulation"]["follow_recoils"])
    cascade_params[0].emin = 5.0
    # TODO: get ed from input_params
    cascade_params[0].ed = 15.0

    densities = np.array([layer["density"] for layer in input_params["layer"]])
    cascade_params[0].pmax = densities**(-1/3) / sqrt(np.pi)
    cascade_params[0].mean_free_path = densities**(-1/3)

    # TODO: get psimin and demin from input_params
    psimin = np.radians(5.0)
    demin = 15.0

    # TODO: get pmaxmin and pmaxmax from input_params
    pmaxmax = 4.0
    cascade_params[0].pmax_max = pmaxmax  # needed for surface layer

    pmaxmin = 0
    #pmaxmax = 1.53
    #pmaxmin = pmaxmax
    pmaxmin = max(pmaxmin, pmaxmax / NPMAX)
    pmax_vals = np.linspace(pmaxmin, pmaxmax, NPMAX)
    cascade_params[0].pmax_vals = pmax_vals[::-1]


    nmat = len(input_params["layer"])

    for ielem1 in range(nelem):
        z1 = elements_params[ielem1].Z
        m1 = elements_params[ielem1].M
        for imat in range(nmat):
            pmax_energies = np.zeros(NPMAX, dtype=np.float64)
            for ielem in range(materials_params[imat].nelem):
                ielem2 = materials_params[imat].ielem[ielem]
                z2 = elements_params[ielem2].Z
                m2 = elements_params[ielem2].M
                rnorm = scatter_params[0].rnorm[ielem1, ielem2]
                for i, pmax in enumerate(pmax_vals):
                    if scatter_params[0].pot_model == "NLHlin":
                        integral = nlhlin.impulse_integral(
                            pmax/rnorm, 
                            scatter_params[0].pot_coefs[ielem1, ielem2]
                        )
                    else:
                        integral = zbl.impulse_integral(
                            pmax/rnorm, 
                            scatter_params[0].pot_coefs[ielem1, ielem2]
                        )
                    energy_psi = 14.39979 * z1 * z2 / rnorm * integral / psimin
                    energy_de = m1/m2 * (
                        (14.39979 * z1 * z2 / rnorm * integral)**2 / demin)
                    pmax_energies[i] = max(pmax_energies[i], 
                                           energy_psi, energy_de)
            cascade_params[0].pmax_energies[ielem1, imat] = pmax_energies[::-1]
            
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
    global NELEM_ION, NELEM_TARGET
    
    # rearrange parameters and calculate derived parameters
    beam_params = _get_beam_params(input_params)
    geometry_params = _get_geometry_params(input_params)
    nelem_target, nelem, elements_params, materials_params = (
        _get_elements_and_materials_params(input_params))
    estop_params = _get_estop_params(input_params, nelem, elements_params)
    scatter_params = _get_scatter_params(input_params, nelem, elements_params)
    cascade_params = _get_cascade_params(input_params, nelem, elements_params,
                                         materials_params, scatter_params)

    # TODO: include n_absc in params
    #cm_scatter.setup(input_params["models"]["scattering integrals"]["n_absc"])

    # collect subarrays into a single structured array
    PARAMS_DTYPE = np.dtype([
        ("rng_seed", np.int64),      # need an even number of int32 for alignment
        ("nelem_target", np.int32),  # without padding, to avoid trouble with Numba when
        ("nelem", np.int32),         # align=True
        ("beam", beam_params.dtype),
        ("cascade", cascade_params.dtype),
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
    params[0].nelem_target = nelem_target
    params[0].nelem = nelem
    params[0].beam = beam_params
    params[0].cascade = cascade_params
    params[0].geometry = geometry_params
    params[0].elements = elements_params
    params[0].materials = materials_params
    params[0].estop = estop_params
    params[0].scatter = scatter_params

    return NELEM_ION, NELEM_TARGET, params