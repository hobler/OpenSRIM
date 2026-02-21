"""Target related definitions and operations.

Currently, only planar targets are supported.

Available functions:
    setup: setup module variables.
    is_inside_target: check if a given position is inside the target
"""
import numpy as np
from numba import jit


def setup(input_params):
    """Define the target.

    Target properties include geometry and materials. Chemical elements
    (incuding beam elements) are defined in a subarray for convenience.
    
    Parameters:
        input_params (dict): input parameters

    Returns:
        (TARGET_PARAMS_DTYPE): target parameters
    """
    layers_params = input_params["layers"]

    zmin = 0.0
    zmax = layers_params["width"][0]

    GEOMETRY_PARAMS_DTYPE = np.dtype([        
        ("zmin", np.float64),
        ("zmax", np.float64),
    ], align=True)

    geometry_params = np.recarray(1, dtype=GEOMETRY_PARAMS_DTYPE)[0]
    geometry_params.zmin = zmin
    geometry_params.zmax = zmax

    # Construct a list of more convenient material dictionaries.
    # Each material has the format
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
    for ilayer in range(len(layers_params["name"])):
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
            "elements": mat_elements,
            "atomic_fractions": atomic_fractions,
            "displacement_energy": displacement_energies,
        }
        materials.append(material)

    # Determine number of distinct chemical elements in the simulation
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
        mat["nelem"] = len(mat["elements"])
        mat["ielem"] = ielem
        del mat["elements"]  # we only need the indices of the elements in the 
                             # materials

    nelem = len(elements)
    # Test:
    print(f"Number of distinct chemical elements: {nelem}")
    print("Elements:")
    for elem in elements:
        print(elem)
    print("Materials:")
    for mat in materials:
        print(mat)

    ELEMENTS_PARAMS_DTYPE = np.dtype([
        ("symbol", "<U2"),
        ("name", "<U12"),
        ("Z", np.uint32),
        ("M", np.float64),
    ], align=True)

    elements_params = np.recarray(nelem, dtype=ELEMENTS_PARAMS_DTYPE)
    for ielem, elem in enumerate(elements):
        elements_params[ielem].symbol = elem["symbol"]
        elements_params[ielem].name = elem["name"]
        elements_params[ielem].Z = elem["Z"]
        elements_params[ielem].M = elem["M"]

    MATERIALS_PARAMS_DTYPE = np.dtype([
        ("name", "<U16"),
        ("density", np.float64),
        ("compound_correction", np.float64),
        ("gas", np.bool_),
        ("nelem", np.uint32),
        ("ielem", np.uint32, nelem),
        ("atomic_fraction", np.float64, nelem),
        ("displacement_energy", np.float64, nelem),
    ], align=True)

    materials_params = np.recarray(len(materials), dtype=MATERIALS_PARAMS_DTYPE)
    for imat, mat in enumerate(materials):
        materials_params[imat].name = mat["name"]
        materials_params[imat].density = mat["density"]
        materials_params[imat].compound_correction = mat["compound_correction"]
        materials_params[imat].gas = mat["gas"]
        materials_params[imat].nelem = mat["nelem"]
        materials_params[imat].ielem[:mat["nelem"]] = mat["ielem"]
        materials_params[imat].atomic_fraction[:mat["nelem"]] = (
            mat["atomic_fractions"])
        materials_params[imat].displacement_energy[:mat["nelem"]] = (
            mat["displacement_energy"])

    TARGET_PARAMS_DTYPE = np.dtype([
        ("geometry", GEOMETRY_PARAMS_DTYPE),
        ("elements", ELEMENTS_PARAMS_DTYPE, nelem),
        ("materials", MATERIALS_PARAMS_DTYPE, len(materials)),
    ], align=True)

    target_params = np.recarray(1, dtype=TARGET_PARAMS_DTYPE)[0]
    target_params.geometry = geometry_params
    target_params.elements = elements_params
    target_params.materials = materials_params

    print(f"target_params.geometry={target_params.geometry}")
    print(f"target_params.elements={target_params.elements}")
    print(f"target_params.materials={target_params.materials}")
    #exit()

    return target_params


@jit(inline = "always")
def is_inside_target(pos, geometry_params):
    """Check if a given position is inside the target.

    Parameters:
        pos (ndarray): position to check (size 3)
        geometry_params (GEOMETRY_PARAMS_DTYPE): geometry parameters

    Returns:
        (bool): whether the position is inside the target
    """
    return geometry_params.zmin <= pos[2] <= geometry_params.zmax