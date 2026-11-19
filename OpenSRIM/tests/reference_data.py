"""Isolate standalone Python imports and persist input-matched reference arrays.

Live Python code is loaded directly from tools/pytrim, whose NLHlin loader
reads the coefficient table relative to its source file.
Missing code requires a saved result for the exact model, parameters and grid;
Numba is never used as a replacement for an unavailable Python implementation.
"""

import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path
import sys
from functools import lru_cache

import numpy as np
import pytest

# sets paths where python source and reference data are expected to be found
TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parents[1]
COEFFICIENTS = TESTS.parent / 'data/NLHlin/NLHlin_3eV.dat'
COEFFICIENT_HASH = hashlib.sha256(COEFFICIENTS.read_bytes()).hexdigest()
PYTHON_SOURCE_DIRECTORY = TESTS.parent / 'tools/pytrim'
# A missing override may exercise cached mode, but live implementations must
# always come from the canonical standalone Python source directory.
PY_SOURCE = Path(os.environ.get('OPEN_SRIM_PYTHON_SOURCE', PYTHON_SOURCE_DIRECTORY))
PY_CACHE = TESTS / 'data/py/reference_v1'


@lru_cache(maxsize=1)
def python_modules():
    """ Load the standalone Python implementation from tools/pytrim, if available, and return a dict of modules.
        if any of the expected source files are missing, return None to indicate that the implementation is unavailable.

    Returns:
        dict or None: A dictionary mapping module names to loaded modules, or None if the implementation is unavailable.

    Raises:
        pytest.fail: If the source directory is not the expected one, to prevent using an alternative implementation.
    """

    # module names must match the source files in tools/pytrim, and the order matters for dependencies
    names = ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter')

    # Check that all expected source files are present in the specified source directory
    if not all((PY_SOURCE / f'{name}.py').is_file() for name in names):
        print(f'Python source unavailable: {PY_SOURCE}; reading saved Python results')
        return None

    # Ensure that the source directory is the expected one, to avoid using an alternative implementation
    if PY_SOURCE.resolve() != PYTHON_SOURCE_DIRECTORY.resolve():
        pytest.fail(f'Python implementation must come from {PYTHON_SOURCE_DIRECTORY}; '
                    f'refusing alternative source directory: {PY_SOURCE}')
    modules = {}

    # Use pytest.MonkeyPatch to temporarily alias the loaded modules into sys.modules, so that dependencies can be imported by name while loading the modules.
    # This resolves imports such as from apsis import Apsis to the intended Python modules without permanently changing sys.modules.
    # The context manager restores the previous sys.modules entries when module loading finishes; imported references remain available in the loaded modules.
    # The loaded modules are stored in a dictionary and returned to the caller, so that they can be used directly without relying on the global import system.
    with pytest.MonkeyPatch.context() as aliases:
        for name in names:
            spec = importlib.util.spec_from_file_location(f'_comparison_{name}', PY_SOURCE / f'{name}.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            modules[name] = module
            aliases.setitem(sys.modules, name, module)
    return modules


def reference_filename(kind, parameters, *, norm, grid_name, extension='npz'):
    """Build a shared readable reference name; exact inputs remain in metadata.

    Numeric labels use 12 decimal places in scientific notation. Metadata checks
    reject collisions, including inputs too close to distinguish in the label.
    """
    def label(value):
        value = str(value).lower()
        if not re.fullmatch(r'[a-z0-9_]+', value):
            raise ValueError(f'Invalid reference label: {value!r}')
        return value

    def number(value):
        value = float(value)
        if not np.isfinite(value):
            raise ValueError('Reference parameters must be finite')
        mantissa, exponent = f'{value:.12e}'.split('e')
        return mantissa.rstrip('0').rstrip('.') + 'e' + exponent

    name = (f'{label(kind)}_{label(parameters["model"])}'
            f'_z{int(parameters["z1"])}_z{int(parameters["z2"])}_norm_{label(norm)}')
    if 'e' in parameters:
        name += f'_e{number(parameters["e"])}'
    if 'r0' in parameters:
        name += f'_r0_{number(parameters["r0"])}'
    if 'n_absc' in parameters:
        name += f'_n{int(parameters["n_absc"])}'
    return f'{name}_grid_{label(grid_name)}.{label(extension)}'


def python_manifest(*, allow_missing=False):
    """ Return the manifest of cached Python reference outputs, validating its schema and hashes.

    Args:
        allow_missing (bool): If True, return an empty manifest if the manifest file is missing. 
                              If False, raise pytest.fail if the manifest file is missing or invalid.

    Returns:
        dict: A dictionary containing the manifest of cached Python reference outputs, including the schema version and the SHA256 hashes of the outputs.

    Raises:
        pytest.fail: If the manifest file is missing or invalid, unless allow_missing is True.
    """

    # The manifest file is expected to be located in the PY_CACHE directory and named 'python_manifest.json'.
    path = PY_CACHE / 'python_manifest.json'

    # If allow_missing is True and the manifest file does not exist, return an empty manifest with schema version 1 and an empty outputs dictionary.
    if allow_missing and not path.exists():
        return {'schema': 1, 'outputs': {}}

    # If the manifest file exists, attempt to read and parse it as JSON. Validate that the schema version is 1, that the 
    # outputs field is a dictionary, and that each output name has a valid SHA256 hash. If any of these checks fail, raise 
    # pytest.fail with an appropriate error message.
    try:
        manifest = json.loads(path.read_text())
        if (manifest.get('schema') != 1 or not isinstance(manifest.get('outputs'), dict)
                or not all(isinstance(name, str) and isinstance(digest, str)
                           and re.fullmatch(r'[0-9a-f]{64}', digest)
                           for name, digest in manifest['outputs'].items())):
            raise ValueError('invalid schema or output hashes')
    except (OSError, ValueError, AttributeError) as exc:
        pytest.fail(f'Missing or invalid Python reference manifest: {path}: {exc}')
    return manifest


def python_result(kind, parameters, grid, compute, *, norm, grid_name):
    """ Return a live or cached result for the given kind, parameters and grid.
    If the Python implementation is available, compute the result and cache it.
    If the implementation is unavailable, load a cached result if it exists.

    Args:
        kind (str): The kind of result, e.g., 'screen', 'scatter', 'chi'.
        parameters (dict): The parameters that define the calculation.
        grid (array-like): The grid of input values for the calculation.
        compute (callable): A function that computes the result if the implementation is available.
        norm (str): Stable normalization/configuration label.
        grid_name (str): Stable test-grid label; distinct inputs require distinct labels.

    Returns:
        np.ndarray: The computed or cached result.

    Raises:
        pytest.fail: If the Python implementation is unavailable and no matching cached result exists, or if
        the cached result's metadata does not match the inputs.
    """

    # Convert the input grid into a NumPy array with 64-bit floating-point precision.
    # This ensures consistent numerical precision and data representation.
    grid = np.asarray(grid, dtype=np.float64)

    # Create metadata describing the exact inputs used for the calculation.
    # The metadata includes the calculation type, parameters, coefficient hash and grid.
    metadata = {'schema': 1, 'kind': kind, 'parameters': parameters,
                'coefficient_sha256': COEFFICIENT_HASH, 'grid': grid.tolist(),
                'norm': norm, 'grid_name': grid_name}
    
    # Convert the metadata dictionary into a sorted JSON string.
    # Sorting ensures that identical metadata always produces the same string.
    encoded = json.dumps(metadata, sort_keys=True)

    path = PY_CACHE / reference_filename(kind, parameters, norm=norm, grid_name=grid_name)

    # A readable name must never silently overwrite a different input grid or
    # parameter set. Coefficient changes may regenerate the same named case.
    if path.is_file() and python_modules() is not None:
        with np.load(path, allow_pickle=False) as saved:
            previous = json.loads(saved['metadata'].item())
        previous['coefficient_sha256'] = COEFFICIENT_HASH
        if previous != metadata:
            pytest.fail(f'Python reference name collision: {path}; use a distinct grid/norm label')

    # If the Python implementation is available, compute the result and cache it.
    if python_modules() is not None:
        manifest = python_manifest(allow_missing=True)
        result = np.asarray(compute(), dtype=float)
        # Ensure that the cache directory exists, creating it if necessary.
        PY_CACHE.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix('.tmp.npz')

        # Calculate and store the SHA256 hashes of all Python source files.
        # These hashes document which source versions produced the results.
        provenance = {name: hashlib.sha256((PY_SOURCE / f'{name}.py').read_bytes()).hexdigest()
                      for name in ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter')}

        # Save the metadata, numerical result and source hashes in an uncompressed NPZ archive.
        # All three pieces of information are stored together.
        np.savez(tmp, metadata=encoded, result=result, source_hashes=json.dumps(provenance))

        # Replace the previous cache file with the newly generated result.
        # This prevents an incomplete temporary file from becoming the final cache.
        tmp.replace(path)

        # Hash the complete archive, including results, metadata and provenance.
        # Update only the freshly computed file; never endorse other old files.
        manifest['outputs'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_path = PY_CACHE / 'python_manifest.json'
        manifest_tmp = manifest_path.with_suffix('.tmp')
        manifest_tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
        manifest_tmp.replace(manifest_path)

        return result

    # If the Python implementation is unavailable, attempt to load a cached result.
    if not path.is_file():
        pytest.fail(f'Python implementation missing and no matching saved result: {path}')
    manifest = python_manifest()
    expected = manifest['outputs'].get(path.name)
    if expected != hashlib.sha256(path.read_bytes()).hexdigest():
        pytest.fail(f'Missing or modified cached Python result: {path}')
    with np.load(path, allow_pickle=False) as saved:
        if saved['metadata'].item() != encoded:
            pytest.fail(f'Python reference metadata does not match inputs: {path}')
        return saved['result'].copy()


@lru_cache(maxsize=1)
def coefficient_rows():
    """ Return a dictionary mapping (Z1, Z2) pairs to their corresponding coefficient rows from the NLHlin coefficient table.
    """

    # Load the NLHlin coefficient table from the specified file and create a dictionary where each key is a tuple of (Z1, Z2) and the value is the corresponding row of coefficients.
    return {(int(row[0]), int(row[1])): row for row in np.loadtxt(COEFFICIENTS)}


class PythonScreen:
    """ A wrapper class for the native Python screening implementation, which can be either ZBL or NLHlin.

        allows test grids to be evaluated with the same interface, regardless of which screening model is used.
    """
    def __init__(self, model, z1, z2, rnorm, *, norm="native"):
        """ Initialize the PythonScreen object with the specified model, element pair, and normalization factor.

        Args:
            model (str): The screening model to use, either 'ZBL' or 'NLHlin'.
            z1 (int): The atomic number of the first element.
            z2 (int): The atomic number of the second element.
            rnorm (float): The normalization length in Angstrom used to express radii in reduced units.
            norm (str): Name of the normalization or screening configuration.

        Returns:
            None
        """
        self.norm = norm
        self.parameters = {'model': model, 'z1': z1, 'z2': z2, 'rnorm': rnorm}

        # If the Python implementation is available for the specified model, create an instance of the corresponding screening class (ZBL_screen or NLHlin_screen) and store it in self.implementation.
        modules = python_modules()
        self.implementation = None

        if modules is not None:

            # Select the appropriate screening class based on the requested model.
            # ZBL uses ZBL_screen; other models use NLHlin_screen.
            cls = modules['zbl'].ZBL_screen if model == 'ZBL' else modules['nlhlin'].NLHlin_screen
            self.implementation = cls(z1, z2, rnorm)

        # Determine the maximum screening radius.
        # ZBL has an infinite range, while NLHlin obtains rmax from the coefficient table.
        # The NLHlin value is converted into reduced length units by dividing by rnorm.
        self.rmax = (np.inf if model == 'ZBL' else
                     coefficient_rows()[tuple(sorted((z1, z2)))][8] / rnorm)

    def __call__(self, radii, *, grid_name):
        """ Calculate the screening potential and its derivative for the given radius using the selected screening model.

        Args:
            radii (array-like): The input radius values for which to calculate the screening potential and its derivative.
            grid_name (str): Name of this radius grid in the reference filename.

        Returns:
            np.ndarray: An array containing the screening potential and its derivative for each input radius.
        """

        # Calculate or load the screening results for the supplied radius values.
        # The lambda function executes the live implementation when available.
        # The result contains the screening potential phi and its radial derivative.
        return python_result('screen', self.parameters, radii,
                             lambda: self.implementation(radii), norm=self.norm, grid_name=grid_name)


def python_scatter(e, p_values, screen, n_absc, *, grid_name):
    """ Return native scattering integrals for the given energy, impact parameters, screening model, and number of abscissas.

    This calls cm_scatter.scatter_integrals and does not reconstruct its output.

    Args:
        e (float): The energy value for which to compute the scattering integrals.
        p_values (array-like): The impact parameter values for which to compute the scattering integrals.
        screen (PythonScreen): An instance of the PythonScreen class representing the screening model to use.
        n_absc (int): The number of abscissas to use in the numerical integration of the scattering integrals.
        grid_name (str): Name of the impact-parameter grid.

    Returns:
        np.ndarray: Shape (2, len(p_values)), with scattering angle theta in the first row and time integral tau in the second.
    """
    def compute():
        """ Compute the scattering integrals using the native Python implementation.
        """

        # Load the cm_scatter module from the Python implementation and set up the number of abscissas for the integration.
        module = python_modules()['cm_scatter']
        module.setup(n_absc)

        # Compute the scattering integrals for each impact parameter using the native implementation and return the results as a transposed NumPy array.
        return np.asarray([module.scatter_integrals(e, p, screen.implementation)
                           for p in p_values]).T

    # Return the computed or cached scattering integrals for the given energy, impact parameters, screening model, and number of abscissas.
    return python_result('scatter', dict(screen.parameters, e=float(e), n_absc=n_absc), p_values, compute,
                         norm=screen.norm, grid_name=grid_name)


def python_chi(u, r0, screen, *, grid_name):
    """ Return the native Python evaluation of (phi, chi) for integration coordinates u, reference radius r0, and screening model.

    Args:
        u (array-like): The dimensionless integration coordinates for which to compute (phi, chi).
        grid_name (str): Name of the integration-coordinate grid.
        r0 (float): The reference radius used in the chi(u) calculation.
        screen (PythonScreen): An instance of the PythonScreen class representing the screening model to use

    Returns:
        np.ndarray: An array containing phi and chi for each input u; shape (2, len(u)) for a one-dimensional grid.
    """

    # Return the computed or cached (phi, chi) values for the integration coordinates, reference radius, and screening model.
    return python_result('chi', dict(screen.parameters, r0=r0), u,
                         lambda: python_modules()['cm_scatter'].calc_phi_chi(u, r0, screen.implementation),
                         norm=screen.norm, grid_name=grid_name)
