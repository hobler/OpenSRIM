"""
Numerical comparison tests for CM scattering.

The numerical test cases are based on the data ranges used by the
plot_* functions in cm_scatter.py, but no plots are generated.

Implementations compared:
    - Python
    - Numba
    - Fortran (live or validated saved references)
"""

from pathlib import Path
from functools import lru_cache

import numpy as np
import pytest

from numba import njit
from reference_data import PythonScreen, python_scatter, python_chi, diagnostic_chi, PY_SOURCE

# Only the Numba implementation comes from simulators/opentrim.
from OpenSRIM.simulators.opentrim import init_params
from OpenSRIM.simulators.opentrim.read_params import read_params
from OpenSRIM.simulators.opentrim import cm_scatter as cm_scatter_nb
from OpenSRIM.simulators.opentrim import zbl as zbl_nb, nlhlin as nlhlin_nb


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

Z1 = 5
Z2 = 14

RNORM_ZBL = 0.4685 / (Z1**0.23 + Z2**0.23)
RNORM_NLHlin = 0.4685 / np.sqrt(np.sqrt(Z1) + np.sqrt(Z2))

R0_VALUES = [0.01, 0.1, 1.0, 10.0]

E_VALUES_ERROR = (1e-4, 0.1, 100.0)

P_VALUES_ERROR = (
    np.linspace(0.02, 30.0, 101),
    np.linspace(0.02, 12.0, 101),
    np.linspace(0.02, 4.0, 101),
)

# Used by plot_tau_over_theta()
E_VALUES_TAU_THETA = np.logspace(-4, 1, 6)
P_VALUES_TAU_THETA = np.linspace(0.02, 30.0, 201)

# Used by plot_sinhalftheta_over_p()
E_VALUES_SIN = np.logspace(-4, 2, 7)
P_VALUES_SIN = np.linspace(0.01, 30.0, 201)

# The current Numba implementation has ROOTS_LEGENDRE = roots_legendre(4).
NUMBA_N_ABSC = 4

# ---------------------------------------------------------------------------
# Numerical comparison tolerances
# ---------------------------------------------------------------------------

# Python <-> Numba
RTOL = 1e-3
ATOL = 1e-6

# The Fortran tau routines interpolate legacy tables rather than using the
# four-point quadrature. Keep their accuracy budget separate from theta.
ATOL_TAU_FORTRAN = 0.015

# Tolerances are quantity- and backend-specific. They are not fitted per point.
# Fortran Legendre theta agrees closely with Python; Numba uses another apsis
# initial estimate. Native Mehler tau is checked with the previous 0.025 budget.
TOLERANCES = {
    "chi": {"numba": (1e-8, 1e-10), "fortran": (1e-8, 1e-10)},
    "theta": {"numba": (3e-4, 1e-11), "fortran": (1e-8, 2e-11)},
    "sin2": {"numba": (5e-4, 1e-10), "fortran": (1e-8, 1e-11)},
    "tau": {"numba": (1e-3, 1e-7), "fortran": (1e-3, 0.025)},
    "boundary": {"numba": (0, 1e-14), "fortran": (0, 1e-14)},
}

F90_DIR = Path(__file__).resolve().parent / "data" / "f90"


# ---------------------------------------------------------------------------
# Potential setup
# ---------------------------------------------------------------------------

def get_python_screen_fun(model):
    """Return the Python screening function."""
    rnorm = RNORM_ZBL if model == "ZBL" else RNORM_NLHlin
    return PythonScreen(model, Z1, Z2, rnorm)


def build_numba_potential(model):
    """Build the OpenTRIM/Numba potential coefficients."""
    assert hasattr(cm_scatter_nb.scatter_integrals, "py_func"), "Numba JIT is disabled"
    assert len(cm_scatter_nb.ROOTS_LEGENDRE[0]) == NUMBA_N_ABSC
    input_params = read_params(Path(cm_scatter_nb.__file__).with_name("defaults.toml"))
    input_params["models"]["potential"] = model

    _, _, params = init_params.get_params(input_params)

    elements = params.elements[0]

    element_index = {
        int(element.Z): idx
        for idx, element in enumerate(elements)
    }

    pot_coefs = params.scatter[0].pot_coefs[
        element_index[Z1],
        element_index[Z2],
    ]

    return input_params["models"]["potential"], pot_coefs


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def all_models():
    """Return the potential models used by the tests."""
    return ("ZBL", "NLHlin")


def chi_u_values():
    """u values used by plot_chi()."""
    return np.linspace(0.001, 0.999, 999)


def chi_near_zero_u_values():
    """u values used by plot_chi_near_zero()."""
    u = np.logspace(-10, 0, 101)
    u[-1] = np.nextafter(1.0, 0.0)
    return u


def chi_near_one_values():
    """Return (1-u, u) values used by plot_chi_near_one()."""
    u1 = np.logspace(-10, 0, 101)
    u1[-1] = np.nextafter(1.0, 0.0)

    u = 1.0 - u1

    return u1, u


def calculate_scatter(e, p_values, screen, pot_model, pot_coefs):
    """Call both production APIs on one grid, returning (theta, tau) arrays."""
    python = calculate_scatter_python(e, p_values, screen)
    numba = calculate_scatter_numba(e, p_values, pot_model, pot_coefs)
    return python, numba


def calculate_scatter_numba(e, p_values, pot_model, pot_coefs, n_absc=None):
    """Call the native JIT API; None exercises the unchanged default call shape."""
    if n_absc is None:
        values = [cm_scatter_nb.scatter_integrals(e, p, pot_model, pot_coefs) for p in p_values]
    else:
        values = [cm_scatter_nb.scatter_integrals(e, p, pot_model, pot_coefs, n_absc) for p in p_values]
    result = np.asarray(values).T
    # Native API returns pi minus theta. No accuracy corrections are applied.
    result[0] = np.pi - result[0]
    return result


def calculate_scatter_python(e, p_values, screen, n_absc=NUMBA_N_ABSC):
    """Evaluate Python with an explicit quadrature order on every call."""
    return python_scatter(e, p_values, screen, n_absc)


def assert_arrays_equal(reference, result, *, rtol=RTOL, atol=ATOL):
    """Compare finite outputs; the second allclose argument sets relative scale."""
    assert np.all(np.isfinite(reference)), "Non-finite reference output"
    assert np.all(np.isfinite(result)), "Non-finite result output"
    np.testing.assert_allclose(
        result, reference, equal_nan=False, rtol=rtol, atol=atol,
    )


def load_fortran_data(filename, grid, columns=2):
    """Read reference outputs and verify their full shape and input grid."""
    path = F90_DIR / filename
    if not path.is_file():
        pytest.fail(f"Missing result file: {path}")
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except (OSError, ValueError) as exc:
        pytest.fail(f"Could not read numerical data from {path}: {exc}")
    assert data.shape == (len(grid), columns), f"Unexpected data shape in {path}"
    # Output finiteness is checked by compare_three_way after recording diagnostics.
    assert_arrays_equal(grid, data[:, 0], rtol=0, atol=1e-14)
    return data


def record_max_error(summaries, name, left, right, coordinates, *,
                     context="", coordinate_name="p", labels=("Python", "Numba"),
                     details=None):
    """Record absolute/relative maxima separately; non-finite outputs stay visible."""
    left, right = np.asarray(left), np.asarray(right)
    finite = np.isfinite(left) & np.isfinite(right)
    entry = summaries.setdefault(name, {
        "error": -np.inf, "location": "no finite pair", "values": {}, "details": {},
        "relative_error": -np.inf, "relative_location": "no nonzero reference",
        "relative_values": {}, "nonfinite": 0, "zero_reference": 0,
    })
    entry["nonfinite"] += int(np.count_nonzero(~finite))
    entry["zero_reference"] += int(np.count_nonzero(finite & (right == 0)))
    def location(index):
        point = f"{coordinate_name}={coordinates[index]}"
        return f"{context}, {point}" if context else point
    indices = np.flatnonzero(finite)
    if indices.size:
        difference = np.abs(left[indices] - right[indices])
        index = indices[np.argmax(difference)]
        if difference.max() > entry["error"]:
            entry.update(error=difference.max(), location=location(index),
                         values=dict(zip(labels, (left[index], right[index]))),
                         details={key: values[index] for key, values in (details or {}).items()})
    indices = np.flatnonzero(finite & (right != 0))
    if indices.size:
        relative = np.abs(left[indices] - right[indices]) / np.abs(right[indices])
        index = indices[np.argmax(relative)]
        if relative.max() > entry["relative_error"]:
            entry.update(relative_error=relative.max(), relative_location=location(index),
                         relative_values=dict(zip(labels, (left[index], right[index]))))


def print_results(title, summaries, *, counts=None, note=None):
    """Central output routine for all numerical test summaries."""
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
    for label, value in (counts or {}).items():
        lines.append(f"{label}: {value}")
    if note:
        lines.append(note)
    print("\n".join(lines))


def compare_three_way(summaries, failures, values, coordinates, context, quantity,
                      *, tolerances, prefix="", coordinate_name="p", fortran_tolerance=None):
    """Check all three pairs and retain diagnostics even when a pair fails."""
    for left_label, right_label in (("Python", "Numba"), ("Python", "Fortran"),
                                    ("Numba", "Fortran")):
        left, right = values[left_label], values[right_label]
        name = f"{prefix}{left_label} vs {right_label}"
        record_max_error(summaries, name, left, right, coordinates, context=context,
                         coordinate_name=coordinate_name, labels=(left_label, right_label))
        tolerance = tolerances[quantity]["numba" if "Numba" in (left_label, right_label) else "fortran"]
        if right_label == "Fortran" and quantity == "tau":
            tolerance = fortran_tolerance or tolerances[quantity]["fortran"]
        try:
            assert_arrays_equal(right, left, rtol=tolerance[0], atol=tolerance[1])
        except AssertionError as exc:
            failures.append(f"{context}, {name}, rtol={tolerance[0]}, atol={tolerance[1]}: {exc}")


def finish_comparison(title, summaries, failures, **kwargs):
    print_results(title, summaries, **kwargs)
    if failures:
        pytest.fail(f"{len(failures)} comparison(s) failed:\n" + "\n".join(failures[:3]))


@njit
def evaluate_numba_screen(radii, model, pot_coefs):
    """Evaluate native screening with the same finite domain as scattering.

    The NLHlin screening kernel itself extrapolates beyond rmax; its caller
    in production applies the cutoff. This adapter supplies that domain
    convention and vector/scalar conversion, not a fitted output correction.
    """
    phi = np.zeros(radii.size)
    derivative = np.zeros(radii.size)
    for index, radius in enumerate(radii):
        if radius < pot_coefs.rmax:
            if model == "ZBL":
                phi[index], derivative[index] = zbl_nb.screen_fun(radius, pot_coefs)
            else:
                phi[index], derivative[index] = nlhlin_nb.screen_fun(radius, pot_coefs)
    return phi, derivative


def compare_chi_reference(model, suffix, coordinates, u):
    """Compare the plot diagnostic applied to each version's screening.

    Like plot_chi*, this diagnoses screening behaviour through chi. It is
    not a claim to call three independently exposed production chi APIs.
    Python uses the existing tools/pytrim calc_phi_chi outputs (live or saved).
    Numba and Fortran form the same plot diagnostic from native screening.
    No coefficients or computed outputs are changed to force agreement.
    """
    screen_py = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)

    def screen_nb(radius):
        radii = np.asarray(radius, dtype=float)
        phi, derivative = evaluate_numba_screen(radii.reshape(-1), pot_model, pot_coefs)
        return phi.reshape(radii.shape), derivative.reshape(radii.shape)

    summaries, failures = {}, []
    for index, r0 in enumerate(R0_VALUES, start=1):
        _, chi_py = python_chi(u, r0, screen_py)
        _, chi_nb = diagnostic_chi(u, r0, screen_nb)
        data = load_fortran_data(
            f"chi_{model.lower()}{suffix}_r0_{index}.txt", coordinates,
        )
        compare_three_way(
            summaries, failures, dict(Python=chi_py, Numba=chi_nb, Fortran=data[:, 1]),
            coordinates, f"r0={r0}", "chi", tolerances=TOLERANCES,
            coordinate_name="1-u" if suffix == "_near_one" else "u",
        )
    finish_comparison(
        f"{model} chi{suffix}: common diagnostic of native screening", summaries, failures,
        note=f"Python source: {PY_SOURCE}; compares screening-derived chi, not internal scattering chi implementations",
    )




# Independent input coverage: light/light, heavy/heavy, very unequal and B/Si.
EXTENDED_PAIRS = ((1, 1), (79, 79), (1, 79), (5, 14))
MEHLER_ORDERS = (512, 2048, 8192, 32768, 131072, 524288)


@lru_cache(maxsize=1)
def coefficient_template():
    return build_numba_potential("NLHlin")[1]


def numba_coefficients(model, z1, z2):
    """Populate a native record using only the backend's own coefficient loader."""
    coefs = coefficient_template().copy()
    coefs.a[:] = 0
    coefs.b[:] = 0
    if model == 'NLHlin':
        norm, a, b, c, d, r34, rmax, k = nlhlin_nb.get_coefs(z1, z2)
    else:
        norm, a, b, r34, k = zbl_nb.get_coefs(z1, z2)
        c = d = 0.
        rmax = np.inf
    coefs.a[:len(a)] = a
    coefs.b[:len(b)] = b
    coefs.c, coefs.d, coefs.rmax, coefs.r34 = c, d, rmax, r34
    coefs.k[:] = 0
    coefs.k[:len(k)] = k
    return norm, coefs


def pair_inputs(model, z1, z2):
    norm, coefs = numba_coefficients(model, z1, z2)
    return PythonScreen(model, z1, z2, norm), coefs


def extended_grid(kind, model, rmax):
    """Fixed input grids, independent of measured errors or output values."""
    if kind == 'pairs':
        energies = (1e-4, .1, 100.)
        points = np.array([.02, .1, .5, 1., 2., 3., 4.])
    elif kind == 'boundary':
        energies = (.001, .1, 10.)
        points = np.array([1e-8, 1e-6, 1e-4, .01, .02, .03, .04])
        if model == 'NLHlin':
            points[3:] = [rmax*(1-1e-4), rmax*(1-1e-6), rmax*(1-1e-8), np.nextafter(rmax, -np.inf)]
    else:
        energies = (.01, 1., 10.)
        points = np.array([.02, .1, 1., 3., 10., 20., 30.])
        if model == 'NLHlin':
            points = rmax*np.array([.01, .1, .3, .5, .9, .99, .999])
    return energies, points


def extended_outputs(directory, kind, model, pair, order, energies, points):
    """Validate both input columns before exposing the native outputs."""
    z1, z2 = pair
    path = directory / f'{kind}_{model}_{z1}_{z2}_n{order}.txt'
    if not path.is_file():
        pytest.fail(f'Missing Fortran reference: {path}')
    data = np.loadtxt(path, ndmin=2)
    assert data.shape == (len(energies)*len(points), 4), str(path)
    assert_arrays_equal(np.repeat(energies, len(points)), data[:, 0], rtol=1e-14, atol=0)
    # At the immediately adjacent float below rmax, even a one-ULP input
    # difference would invalidate the interpretation of an output mismatch.
    grid_tolerance = 0 if kind == 'boundary' else 1e-14
    assert_arrays_equal(np.tile(points, len(energies)), data[:, 1],
                        rtol=grid_tolerance, atol=grid_tolerance)
    return data[:, 2:].reshape(len(energies), len(points), 2).transpose(0, 2, 1)
