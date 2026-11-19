"""
Numerical comparison tests for CM scattering.

Implementations compared:
    - Python (live or validated saved references)
    - Numba
    - Fortran (live or validated saved references)
"""

from pathlib import Path
from functools import lru_cache

import numpy as np
import pytest

from numba import njit
from reference_data import PythonScreen, python_scatter, python_chi, PY_SOURCE, reference_filename

from OpenSRIM.simulators.opentrim import init_params
from OpenSRIM.simulators.opentrim.read_params import read_params
from OpenSRIM.simulators.opentrim import cm_scatter as cm_scatter_nb
from OpenSRIM.simulators.opentrim import zbl as zbl_nb, nlhlin as nlhlin_nb


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Boron (5) on Silicon (14) is the production default for OpenTRIM.
Z1 = 5
Z2 = 14

# The two supported screening models have different normalization formulas.
RNORM_ZBL = 0.4685 / (Z1**0.23 + Z2**0.23)
RNORM_NLHlin = 0.4685 / np.sqrt(np.sqrt(Z1) + np.sqrt(Z2))

# define the different inintial radius values to test the chi function for each model
R0_VALUES = [0.01, 0.1, 1.0, 10.0]

# define the energy values used for the error tests, which are the same for both models
E_VALUES_ERROR = (1e-4, 0.1, 100.0)

# Define three different impact parameter grids.
# Each grid contains 101 equally spaced values.
P_VALUES_ERROR = (
    np.linspace(0.02, 30.0, 101),
    np.linspace(0.02, 12.0, 101),
    np.linspace(0.02, 4.0, 101),
)

# Used by plot_tau_over_theta()
# np.logspace creates logarithmically spaced values.
E_VALUES_TAU_THETA = np.logspace(-4, 1, 6)
P_VALUES_TAU_THETA = np.linspace(0.02, 30.0, 201)

# Used by plot_sinhalftheta_over_p()
E_VALUES_SIN = np.logspace(-4, 2, 7)
P_VALUES_SIN = np.linspace(0.01, 30.0, 201)

# The production default selects the four-point rule from immutable cached rules.
NUMBA_N_ABSC = 4

# ---------------------------------------------------------------------------
# Numerical comparison tolerances
# ---------------------------------------------------------------------------

# All backends use Legendre quadrature
# tolerances are set for each quantity and backend pair.
# The tolerances are used in the compare_three_way function to check if the results from different backends are within acceptable limits.
TOLERANCES = {
    # Fallback (rtol, atol) for assert_arrays_equal; current callers override both.
    "default": (1e-3, 1e-6),
    "chi": {"numba": (1e-11, 1e-12), "fortran": (1e-11, 1e-12)},
    "theta": {"numba": (3e-4, 1e-11), "fortran": (1e-8, 2e-11)},
    "sin2": {"numba": (5e-4, 1e-10), "fortran": (1e-8, 1e-11)},
    "tau": {"numba": (1e-3, 1e-7), "fortran": (1e-8, 2e-8)},
    "boundary": {"numba": (0, 1e-14), "fortran": (0, 1e-14)},
}

# fortran reference data directory, relative to this file.
F90_DIR = Path(__file__).resolve().parent / "data" / "f90"


# ---------------------------------------------------------------------------
# Potential setup
# ---------------------------------------------------------------------------

def get_python_screen_fun(model):
    """Return a callable Python screening function for the requested model and the default B/Si normalization.

    Args:
        model (str): The screening model to use, either "ZBL" or "NLHlin".

    Returns:
        PythonScreen: A callable screening function for the specified model.
    """

    rnorm = RNORM_ZBL if model == "ZBL" else RNORM_NLHlin
    return PythonScreen(model, Z1, Z2, rnorm)


def build_numba_potential(model):
    """Build the Numba potential for the requested model.

    Args:
        model (str): The potential model to use, either "ZBL" or "NLHlin".

    Returns:
        tuple: A tuple containing the potential model and its coefficients.
    """

    # Check that the Numba JIT is enabled and the cached Legendre roots are available.
    assert hasattr(cm_scatter_nb.scatter_integrals, "py_func"), "Numba JIT is disabled"
    assert len(cm_scatter_nb.ROOTS_LEGENDRE[0]) == NUMBA_N_ABSC

    # Load the default parameters from the TOML file and set the potential model.
    input_params = read_params(Path(cm_scatter_nb.__file__).with_name("defaults.toml"))
    input_params["models"]["potential"] = model

    # Get the parameters for the specified model and extract the elements.
    _, _, params = init_params.get_params(input_params)

    # Build a mapping from atomic number to index in the elements list.
    elements = params.elements[0]
    element_index = {
        int(element.Z): idx
        for idx, element in enumerate(elements)
    }

    # Get the potential coefficients for the specified element pair (Z1, Z2).
    pot_coefs = params.scatter[0].pot_coefs[
        element_index[Z1],
        element_index[Z2],
    ]

    return input_params["models"]["potential"], pot_coefs


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def all_models():
    """Return the two supported potential models for parameterized tests.
    """
    return ("ZBL", "NLHlin")


def chi_u_values():
    """Return 999 evenly spaced u samples for chi tests.
    """
    return np.linspace(0.001, 0.999, 999)


def chi_near_zero_u_values():
    """Generate 101 logarithmically spaced u samples for the small-radius chi diagnostic.
    The last value is the next float below 1.0 to avoid a singularity in the integrand.
    """
    u = np.logspace(-10, 0, 101)
    u[-1] = np.nextafter(1.0, 0.0)
    return u


def chi_near_one_values():
    """Generate 101 logarithmically spaced u samples for the large-radius chi diagnostic.
    The last value is the next float below 1.0 to avoid a singularity in the integrand.
    """
    u1 = np.logspace(-10, 0, 101)
    u1[-1] = np.nextafter(1.0, 0.0)
    u = 1.0 - u1
    return u1, u


def calculate_scatter(e, p_values, screen, pot_model, pot_coefs, *, grid_name):
    """Evaluate both Python and Numba scattering integrals for the same inputs and return their results.

    Args:
        e (float): The energy value for the scattering calculation.
        p_values (array-like): The impact parameter values for the scattering calculation.
        screen (callable): The Python screening function for the specified model.
        pot_model (str): The potential model to use, either "ZBL" or "NLHlin".
        pot_coefs (array-like): The potential coefficients for the specified element pair

        grid_name (str): Explicit impact-grid label used in reference filenames.

    Returns:
        tuple: A tuple containing the results from the Python and Numba scattering calculations.
    """

    python = calculate_scatter_python(e, p_values, screen, grid_name=grid_name)
    numba = calculate_scatter_numba(e, p_values, pot_model, pot_coefs)
    return python, numba


def calculate_scatter_numba(e, p_values, pot_model, pot_coefs, n_absc=None):
    """Evaluate live-or-cached Numba at an explicit quadrature order.

    The default is four points, matching the production Numba default.

    Args:
        e (float): The energy value for the scattering calculation.
        p_values (array-like): The impact parameter values for the scattering calculation.
        pot_model (str): The potential model to use, either "ZBL" or "NLHlin".
        pot_coefs (array-like): The potential coefficients for the specified element pair
        n_absc (int, optional): The number of abscissas for the quadrature rule. If None, the default of 4 is used.

    Returns:
        np.ndarray: A 2D array containing the results from the Numba scattering calculation, with shape (2, len(p_values)),
        where the first row contains the scattering angles and the second row contains the scattering times. The first row is 
        adjusted to return pi minus theta, as the native API returns pi minus theta.
    """
    if n_absc is None:
        values = [cm_scatter_nb.scatter_integrals(e, p, pot_model, pot_coefs) for p in p_values]
    else:
        values = [cm_scatter_nb.scatter_integrals(e, p, pot_model, pot_coefs, n_absc) for p in p_values]
    result = np.asarray(values).T

    # The native API returns pi minus theta, so we adjust the first row to return theta instead.
    result[0] = np.pi - result[0]
    return result


def calculate_scatter_python(e, p_values, screen, n_absc=NUMBA_N_ABSC, *, grid_name):
    """Evaluate the Python scattering integrals at an explicit quadrature order.

    Args:
        e (float): The energy value for the scattering calculation.
        p_values (array-like): The impact parameter values for the scattering calculation.
        screen (callable): The Python screening function for the specified model.
        n_absc (int, optional): The number of abscissas for the quadrature rule. If None, the default of 4 is used.

        grid_name (str): Explicit impact-grid label used in reference filenames.

    Returns:
        np.ndarray: A 2D array containing the results from the Python scattering calculation
    """
    return python_scatter(e, p_values, screen, n_absc, grid_name=grid_name)


def assert_arrays_equal(reference, result, *,
                        rtol=TOLERANCES["default"][0], atol=TOLERANCES["default"][1]):
    """Checks that two arrays are equal within a given tolerance, and that all values are finite.

    Args:
        reference (np.ndarray): The reference array to compare against.
        result (np.ndarray): The result array to compare.
        rtol (float, optional): Relative tolerance, defaulting to TOLERANCES["default"][0].
        atol (float, optional): Absolute tolerance, defaulting to TOLERANCES["default"][1].

    Raises:
        AssertionError: If the arrays are not equal within the given tolerance, or if any values are non-finite.

    Returns:
        None
    """
    assert np.all(np.isfinite(reference)), "Non-finite reference output"
    assert np.all(np.isfinite(result)), "Non-finite result output"
    np.testing.assert_allclose(
        result, reference, equal_nan=False, rtol=rtol, atol=atol,
    )


def load_fortran_data(filename, grid, columns=2):
    """Load a Fortran reference file and check its shape and grid values.

    Args:
        filename (str): The name of the Fortran reference file to load.
        grid (np.ndarray): The expected grid values to compare against the first column of the loaded data.
        columns (int, optional): The expected number of columns in the loaded data. Default is 2.

    Returns:
        np.ndarray: The loaded Fortran reference data as a 2D array.

    Raises:
        pytest.fail: If the Fortran reference file is missing or cannot be read, or if the loaded data does not have the expected shape or grid values.
    """

    # Check that the Fortran reference file exists and is readable, 
    # and that its shape matches the expected grid and number of columns.
    path = F90_DIR / filename
    if not path.is_file():
        pytest.fail(f"Missing result file: {path}")
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except (OSError, ValueError) as exc:
        pytest.fail(f"Could not read numerical data from {path}: {exc}")
    assert data.shape == (len(grid), columns), f"Unexpected data shape in {path}"

    # Check that the first column of the loaded data matches the expected grid values within a tight tolerance.
    assert_arrays_equal(grid, data[:, 0], rtol=0, atol=1e-14)
    return data


def record_max_error(summaries, name, left, right, coordinates, *,
                     context="", coordinate_name="p", labels=("Python", "Numba"),
                     details=None):
    """Accumulate separate absolute/relative maxima over repeated grid comparisons.

    Relative error uses abs(right) as denominator and excludes zero references
    only from that statistic. Absolute checks still include zeros. Record values,
    coordinates and optional details at the absolute maximum; 
    Count nonfinite pairs for the final report.

    Args:
        summaries (dict): A dictionary to store the accumulated error summaries.
        name (str): The name of the comparison being recorded.
        left (np.ndarray): The left-hand array of values to compare.
        right (np.ndarray): The right-hand array of values to compare.
        coordinates (np.ndarray): The coordinates corresponding to the values being compared.
        context (str, optional): Additional context information to include in the error report. Default is an empty string.
        coordinate_name (str, optional): The name of the coordinate being compared. Default is "p".
        labels (tuple, optional): A tuple of labels for the left and right arrays. Default is ("Python", "Numba").
        details (dict, optional): A dictionary of additional details to include in the error report. Default is None.

    Returns:
        None
    """

    # Convert the left and right arrays to NumPy arrays and create a boolean mask for finite values.
    left, right = np.asarray(left), np.asarray(right)
    finite = np.isfinite(left) & np.isfinite(right)

    # Use the setdefault method to create a new entry in the summaries dictionary if it doesn't already exist.
    entry = summaries.setdefault(name, {
        "error": -np.inf, "location": "no finite pair", "values": {}, "details": {},
        "relative_error": -np.inf, "relative_location": "no nonzero reference",
        "relative_values": {}, "nonfinite": 0, "zero_reference": 0,
    })
    entry["nonfinite"] += int(np.count_nonzero(~finite))
    entry["zero_reference"] += int(np.count_nonzero(finite & (right == 0)))

    def location(index):
        """Return a string describing the location of the maximum error in the comparison.
        """
        point = f"{coordinate_name}={coordinates[index]}"
        return f"{context}, {point}" if context else point

    # Find the indices of the finite pairs and calculate the absolute differences between the left and right arrays at those indices.
    indices = np.flatnonzero(finite)
    if indices.size:
        difference = np.abs(left[indices] - right[indices])
        index = indices[np.argmax(difference)]
        if difference.max() > entry["error"]:
            entry.update(error=difference.max(), location=location(index),
                         values=dict(zip(labels, (left[index], right[index]))),
                         details={key: values[index] for key, values in (details or {}).items()})

    # Find the indices of the finite pairs where the right array is not zero and calculate 
    # the relative differences between the left and right arrays at those indices.
    indices = np.flatnonzero(finite & (right != 0))
    if indices.size:
        relative = np.abs(left[indices] - right[indices]) / np.abs(right[indices])
        index = indices[np.argmax(relative)]
        if relative.max() > entry["relative_error"]:
            entry.update(relative_error=relative.max(), relative_location=location(index),
                         relative_values=dict(zip(labels, (left[index], right[index]))))


def print_results(title, summaries, *, counts=None, note=None):
    """Print a summary of the accumulated error statistics for each comparison.

    Args:
        title (str): The title of the comparison being reported.
        summaries (dict): A dictionary containing the accumulated error summaries for each comparison.
        counts (dict, optional): A dictionary of additional counts to include in the report. Default is None.
        note (str, optional): An additional note to include in the report. Default is None.

    Returns:
        None
    """

    # Print the title of the comparison and iterate over each comparison in the summaries dictionary.
    lines = [f"\n--- {title} ---"]
    for name, summary in summaries.items():
        lines.extend([
            f"{name}:",
            f"  max absolute difference : {summary['error']:.10e}",
            f"  at                      : {summary['location']}",
        ])
        for label, value in summary["values"].items():
            lines.append(f"  {label:<24}: {value:.10e}")
        if summary["relative_error"] >= 0:
            lines.append(f"  max relative difference : {summary['relative_error']:.10e}")
            lines.append(f"  relative maximum at     : {summary['relative_location']}")
            for label, value in summary["relative_values"].items():
                lines.append(f"  relative {label:<15}: {value:.10e}")
        else:
            lines.append("  relative difference     : undefined (all references zero)")
        lines.append(f"  non-finite pairs        : {summary['nonfinite']}")
        lines.append(f"  zero-reference pairs    : {summary['zero_reference']}")
        for label, value in summary["details"].items():
            lines.append(f"  {label:<24}: {value:.16e}")

    # Print any additional counts and notes provided in the arguments.
    for label, value in (counts or {}).items():
        lines.append(f"{label}: {value}")
    if note:
        lines.append(note)
    print("\n".join(lines))


def compare_three_way(summaries, failures, values, coordinates, context, quantity,
                      *, tolerances, prefix="", coordinate_name="p"):
    """Compare the values from three different backends and record any discrepancies.

    Args:
        summaries (dict): A dictionary to store the accumulated error summaries.
        failures (list): A list to store any failures encountered during the comparison.
        values (dict): A dictionary containing the values from each backend to compare.
        coordinates (np.ndarray): The coordinates corresponding to the values being compared.
        context (str): Additional context information to include in the error report.
        quantity (str): The quantity being compared (e.g., "chi", "theta", "sin2", "tau", "boundary").
        tolerances (dict): A dictionary containing the tolerances for each quantity and backend pair.
        prefix (str, optional): A prefix to include in the comparison name. Default is an empty string.
        coordinate_name (str, optional): The name of the coordinate being compared. Default is "p".

    Returns:
        None
    """

    # Iterate over each pair of backends and compare their values, recording any discrepancies in the summaries and failures lists.
    for left_label, right_label in (("Python", "Numba"), ("Python", "Fortran"),
                                    ("Numba", "Fortran")):

        # Get the values for the left and right backends and construct a name for the comparison.
        left, right = values[left_label], values[right_label]
        name = f"{prefix}{left_label} vs {right_label}"

        # Record the maximum error between the two backends and check if the values are equal within the specified tolerances.
        record_max_error(summaries, name, left, right, coordinates, context=context,
                         coordinate_name=coordinate_name, labels=(left_label, right_label))

        # Determine the appropriate tolerance for the comparison based on the quantity and the backends being compared.
        tolerance = tolerances[quantity]["numba" if "Numba" in (left_label, right_label) else "fortran"]

        # Use the assert_arrays_equal function to check if the left and right arrays are equal within the specified tolerances, 
        # and record any failures encountered during the comparison.
        try:
            assert_arrays_equal(right, left, rtol=tolerance[0], atol=tolerance[1])
        except AssertionError as exc:
            failures.append(f"{context}, {name}, rtol={tolerance[0]}, atol={tolerance[1]}: {exc}")


def finish_comparison(title, summaries, failures, **kwargs):
    """Print the results of the comparison and fail the test if any discrepancies were found.

    Args:
        title (str): The title of the comparison being reported.
        summaries (dict): A dictionary containing the accumulated error summaries for each comparison.
        failures (list): A list of any failures encountered during the comparison.
        **kwargs: Additional keyword arguments to pass to the print_results function.

    Returns:
        None

    Raises:
        pytest.fail: If any discrepancies were found during the comparison, with a message summarizing the 
                     number of failures and the first few failure messages.
    """

    # Print the results of the comparison using the print_results function, passing in the title, summaries, and any additional keyword arguments.
    print_results(title, summaries, **kwargs)
    if failures:
        pytest.fail(f"{len(failures)} comparison(s) failed:\n" + "\n".join(failures[:3]))


@njit
def evaluate_numba_screen(radii, model, pot_coefs):
    """ Evaluate the Numba screening function for a given set of radii, model, and potential coefficients.

    Args:
        radii (np.ndarray): An array of radii at which to evaluate the screening function.
        model (str): The screening model to use, either "ZBL" or "NLHlin".
        pot_coefs (np.ndarray): The potential coefficients for the specified element pair.

    Returns:
        tuple: A tuple containing two arrays: the first array contains the screening function values (phi) 
        and the second array contains the derivatives of the screening function with respect to the radii.
    """

    # create empty arrays to hold the screening function values and derivatives for each radius
    phi = np.empty(radii.size)
    derivative = np.empty(radii.size)

    # iterate over each radius and evaluate the screening function and its derivative using the Numba implementation
    for index, radius in enumerate(radii):
        phi[index], derivative[index] = cm_scatter_nb.screen_with_cutoff(radius, model, pot_coefs)
    return phi, derivative


def compare_chi_reference(model, suffix, coordinates, u):
    """ Compare the chi function values from Python, Numba, and Fortran implementations for a given screening model and suffix.

    Args:
        model (str): The screening model to use, either "ZBL" or "NLHlin".
        suffix (str): A suffix to append to the Fortran reference file names.
        coordinates (np.ndarray): The coordinates corresponding to the chi function values.
        u (np.ndarray): The input values for the chi function.

    Returns:
        None
    """

    # Get the Python screening function and build the Numba potential for the specified model.
    screen_py = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)

    # Compare the chi function values from Python, Numba, and Fortran implementations for each initial radius value in R0_VALUES.
    summaries, failures = {}, []
    for r0 in R0_VALUES:
        grid_name = suffix.lstrip('_') or 'regular'
        _, chi_py = python_chi(u, r0, screen_py, grid_name=grid_name)
        _, chi_nb = cm_scatter_nb.calc_phi_chi(u, r0, pot_model, pot_coefs)
        data = load_fortran_data(
            reference_filename('chi', dict(screen_py.parameters, r0=r0),
                               norm=screen_py.norm, grid_name=grid_name, extension='txt'), coordinates,
        )

        # Compare the chi function values from Python, Numba, and Fortran implementations and record any discrepancies in the summaries and failures lists.
        compare_three_way(
            summaries, failures, dict(Python=chi_py, Numba=chi_nb, Fortran=data[:, 1]),
            coordinates, f"r0={r0}", "chi", tolerances=TOLERANCES,
            coordinate_name="1-u" if suffix == "_near_one" else "u",
        )

    # Print the results of the comparison and fail the test if any discrepancies were found.
    finish_comparison(
        f"{model} chi{suffix}: native chi", summaries, failures,
        note=f"Python source: {PY_SOURCE}; compares production chi routines",
    )


# more element pairs for extended tests, including the default B/Si pair and some extreme cases
EXTENDED_PAIRS = ((1, 1), (79, 79), (1, 79), (5, 14))


@lru_cache(maxsize=1)
def coefficient_template():
    """Return a fresh structured record for an element pair, with all fields initialized to zero.

    The template is cached to avoid repeated allocations, but a copy should be made for each use.

    Returns:
        np.ndarray: A structured record with fields for the potential coefficients, all initialized to zero.
    """

    return build_numba_potential("NLHlin")[1]


def numba_coefficients(model, z1, z2):
    """Return the normalization and coefficients for a given element pair and screening model.

    Args:
        model (str): The screening model to use, either "ZBL" or "NLHlin".
        z1 (int): The atomic number of the first element.
        z2 (int): The atomic number of the second element.

    Returns:
        tuple: A tuple containing the normalization factor and a structured record of the potential 
        coefficients for the specified element pair and screening model.
    """

    # Get a fresh copy of the coefficient template and initialize all fields to zero.
    coefs = coefficient_template().copy()
    coefs.a[:] = 0
    coefs.b[:] = 0

    # Depending on the specified screening model, retrieve the normalization factor 
    # and potential coefficients for the given element pair (z1, z2) using the appropriate Numba implementation.
    if model == 'NLHlin':
        norm, a, b, c, d, r34, rmax, k = nlhlin_nb.get_coefs(z1, z2)
    else:
        norm, a, b, r34, k = zbl_nb.get_coefs(z1, z2)
        c = d = 0.
        rmax = np.inf

    # Populate the structured record with the retrieved coefficients and return the normalization factor and the structured record.
    coefs.a[:len(a)] = a
    coefs.b[:len(b)] = b
    coefs.c, coefs.d, coefs.rmax, coefs.r34 = c, d, rmax, r34
    coefs.k[:] = 0
    coefs.k[:len(k)] = k
    return norm, coefs


def pair_inputs(model, z1, z2):
    """Return a PythonScreen instance and the Numba coefficients for a given element pair and screening model.

    Args:
        model (str): The screening model to use, either "ZBL" or "NLHlin".
        z1 (int): The atomic number of the first element.
        z2 (int): The atomic number of the second element.

    Returns:
        tuple: A tuple containing a PythonScreen instance for the specified element pair and screening model, 
        and a structured record of the potential coefficients for the specified element pair and screening model.
    """

    # Get the normalization factor and potential coefficients for the specified element pair and screening model using the numba_coefficients function.
    norm, coefs = numba_coefficients(model, z1, z2)
    return PythonScreen(model, z1, z2, norm), coefs


def extended_grid(kind, model, rmax):
    """Return the energy and impact parameter grids for extended tests based on the specified kind and screening model.

    Args:
        kind (str): The type of grid to generate, either "pairs" or "boundary".
        model (str): The screening model to use, either "ZBL" or "NLHlin".
        rmax (float): The maximum radius for the boundary grid.

    Returns:
        tuple: A tuple containing two arrays: the first array contains the energy values for the grid, 
        and the second array contains the impact parameter values for the grid.

    Raises:
        ValueError: If the specified kind of grid is not supported.
    """

    # Depending on the specified kind of grid, generate the appropriate energy and impact parameter values for the extended tests.
    if kind == 'pairs':
        energies = (1e-4, .1, 100.)
        points = np.array([.02, .1, .5, 1., 2., 3., 4.])
    elif kind == 'boundary':
        energies = (.001, .1, 10.)
        points = np.array([1e-8, 1e-6, 1e-4, .01, .02, .03, .04])
        if model == 'NLHlin':
            points[3:] = [rmax*(1-1e-4), rmax*(1-1e-6), rmax*(1-1e-8), np.nextafter(rmax, -np.inf)]
    else:
        raise ValueError(f'Unsupported extended grid: {kind}')
    return energies, points


def extended_outputs(directory, kind, model, pair, order, energies, points):
    """Load the Fortran reference outputs for a given element pair, screening model, and quadrature order.

    Args:
        directory (Path): The directory containing the Fortran reference files.
        kind (str): The type of grid to load, either "pairs" or "boundary".
        model (str): The screening model to use, either "ZBL" or "NLHlin".
        pair (tuple): A tuple containing the atomic numbers of the two elements (z1, z2).
        order (int): The quadrature order used in the Fortran reference files.
        energies (np.ndarray): The energy values for the grid.
        points (np.ndarray): The impact parameter values for the grid.

    Returns:
        np.ndarray: A 3D array containing the Fortran reference outputs, with shape (len(energies), 2, len(points)), 
        where the first dimension corresponds to the energy values

    Raises:
        pytest.fail: If the Fortran reference file is missing or cannot be read, or if the loaded data does not have the expected shape or grid values.
    """

    z1, z2 = pair
    values = []
    for energy in energies:
        filename = reference_filename(
            'scatter', dict(model=model, z1=z1, z2=z2, e=energy, n_absc=order),
            norm='native', grid_name=kind, extension='txt')
        path = directory / filename
        if not path.is_file():
            pytest.fail(f'Missing Fortran reference: {path}')
        data = np.loadtxt(path, ndmin=2)
        assert data.shape == (len(points), 4), str(path)
        assert_arrays_equal(np.full(len(points), energy), data[:, 0], rtol=1e-14, atol=0)
        grid_tolerance = 0 if kind == 'boundary' else 1e-14
        assert_arrays_equal(points, data[:, 1], rtol=grid_tolerance, atol=grid_tolerance)
        values.append(data[:, 2:].T)
    return np.asarray(values)
