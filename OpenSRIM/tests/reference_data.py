"""Isolate standalone Python imports and persist input-matched reference arrays.

Live Python code comes exclusively from tools/pytrim. A temporary installation
supplies its expected coefficient filename without editing production sources.
Missing code requires a saved result for the exact model, parameters and grid;
Numba is never used as a replacement for an unavailable Python implementation.
"""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from functools import lru_cache

import numpy as np
import pytest

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parents[1]
COEFFICIENTS = TESTS.parent / 'data/NLHlin/NLHlin_3eV.dat'
COEFFICIENT_HASH = hashlib.sha256(COEFFICIENTS.read_bytes()).hexdigest()
PYTHON_SOURCE_DIRECTORY = TESTS.parent / 'tools/pytrim'
# A missing override may exercise cached mode, but live implementations must
# always come from the canonical standalone Python source directory.
PY_SOURCE = Path(os.environ.get('OPEN_SRIM_PYTHON_SOURCE', PYTHON_SOURCE_DIRECTORY))
PY_CACHE = TESTS / 'data/py/reference_v1'
# Copy code byte-for-byte into an isolated test installation, selecting the
# canonical input under the filename that the standalone loader expects.
# Keep this object alive while cached modules refer to files in the installation.
# TemporaryDirectory removes the directory when cleaned up, normally at shutdown.
_installation = tempfile.TemporaryDirectory(prefix='opensrim-python-reference-')


@lru_cache(maxsize=1)
def python_modules():
    """Load an isolated, byte-for-byte copy of the standalone Python modules.
    
    Return a name-to-module mapping, or None if any required source is missing.
    A complete live implementation must resolve to tools/pytrim; another tree
    is rejected. The temporary copy lets NLHlin find the canonical coefficient
    table next to its own __file__ as dmol_coeffs_rmax.dat, without changing the
    original source tree or rewriting its file-loading code.
    
    Modules are loaded in dependency order. Temporary sys.modules aliases
    resolve their unqualified imports to these Python copies and are restored
    after loading. The returned module objects retain their imported references.
    The cached mapping and temporary installation are reused within the process;
    import errors in a complete source tree propagate rather than using old data.
    """
    names = ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter')
    if not all((PY_SOURCE / f'{name}.py').is_file() for name in names):
        print(f'Python source unavailable: {PY_SOURCE}; reading saved Python results')
        return None
    if PY_SOURCE.resolve() != PYTHON_SOURCE_DIRECTORY.resolve():
        pytest.fail(f'Python implementation must come from {PYTHON_SOURCE_DIRECTORY}; '
                    f'refusing alternative source directory: {PY_SOURCE}')
    directory = Path(_installation.name)
    for name in names:
        shutil.copy2(PY_SOURCE / f'{name}.py', directory / f'{name}.py')
    # NLHlin resolves this filename relative to its module's __file__, not cwd.
    shutil.copy2(COEFFICIENTS, directory / 'dmol_coeffs_rmax.dat')
    modules = {}
    # Dependency order matters: apsis imports table1d; screening imports apsis.
    # Restore previous aliases afterwards so other imports are not redirected.
    with pytest.MonkeyPatch.context() as aliases:
        for name in names:
            spec = importlib.util.spec_from_file_location(f'_comparison_{name}', directory / f'{name}.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            modules[name] = module
            aliases.setitem(sys.modules, name, module)
    return modules


def python_result(kind, parameters, grid, compute):
    """Compute a live result or load one with exactly matching input metadata.
    
    The cache key includes kind, parameters, the complete float64 grid and the
    coefficient digest. Live code is always executed, even if a cache file exists;
    compute is never called when the implementation is unavailable. Source hashes
    are saved as provenance, not enforced as part of the fallback cache key.
    
    Results are written through a temporary NPZ file. NaN/Inf are retained for
    the numerical assertions to diagnose rather than silently repaired.
    """
    grid = np.asarray(grid, dtype=np.float64)
    metadata = {'schema': 1, 'kind': kind, 'parameters': parameters,
                'coefficient_sha256': COEFFICIENT_HASH, 'grid': grid.tolist()}
    # Hash inputs, not filenames/timestamps: a changed grid must select a new cache.
    encoded = json.dumps(metadata, sort_keys=True)
    key = hashlib.sha256(encoded.encode()).hexdigest()
    path = PY_CACHE / f'{kind}_{key}.npz'
    if python_modules() is not None:
        result = np.asarray(compute(), dtype=float)
        # Preserve non-finite production outputs for diagnostics, never repair them.
        PY_CACHE.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix('.tmp.npz')
        provenance = {name: hashlib.sha256((PY_SOURCE / f'{name}.py').read_bytes()).hexdigest()
                      for name in ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter')}
        np.savez(tmp, metadata=encoded, result=result, source_hashes=json.dumps(provenance))
        tmp.replace(path)
        return result
    if not path.is_file():
        pytest.fail(f'Python implementation missing and no matching saved result: {path}')
    with np.load(path, allow_pickle=False) as saved:
        if saved['metadata'].item() != encoded:
            pytest.fail(f'Python reference metadata does not match inputs: {path}')
        return saved['result'].copy()


@lru_cache(maxsize=1)
def coefficient_rows():
    """Cache canonical coefficient rows indexed by their stored atomic-number pair.
    
    Callers sort a requested pair before lookup; column 8 supplies physical rmax.
    """
    return {(int(row[0]), int(row[1])): row for row in np.loadtxt(COEFFICIENTS)}


class PythonScreen:
    """Describe screening inputs and expose a live-or-cached callable.
    
    Keeping model, pair, normalization and rmax separate from the implementation
    allows test grids to be constructed even when the Python sources are absent.
    """
    def __init__(self, model, z1, z2, rnorm):
        """Select the native screening class if available and retain cache parameters.
        
        NLHlin rmax is converted from the coefficient-table length to reduced units
        using rnorm; ZBL has infinite range. No Numba substitute is loaded.
        """
        self.parameters = {'model': model, 'z1': z1, 'z2': z2, 'rnorm': rnorm}
        modules = python_modules()
        self.implementation = None
        if modules is not None:
            cls = modules['zbl'].ZBL_screen if model == 'ZBL' else modules['nlhlin'].NLHlin_screen
            self.implementation = cls(z1, z2, rnorm)
        self.rmax = (np.inf if model == 'ZBL' else
                     coefficient_rows()[tuple(sorted((z1, z2)))][8] / rnorm)

    def __call__(self, radii):
        """Return phi and dphi/dr for radii in this object's reduced length units.
        
        The complete radius grid participates in cache matching.
        """
        return python_result('screen', self.parameters, radii,
                             lambda: self.implementation(radii))


def python_scatter(e, p_values, screen, n_absc):
    """Return theta and tau as an array of shape (2, len(p_values)).
    
    Energy, model, element pair, normalization and quadrature order identify the
    calculation. Setting the Python rule on each live call prevents a preceding
    high-order diagnostic from changing the next default-order comparison.
    """
    def compute():
        """Set the requested native rule, then evaluate each impact parameter.
        """
        module = python_modules()['cm_scatter']
        module.setup(n_absc)
        return np.asarray([module.scatter_integrals(e, p, screen.implementation)
                           for p in p_values]).T
    return python_result('scatter', dict(screen.parameters, e=float(e), n_absc=n_absc), p_values, compute)


def python_chi(u, r0, screen):
    """Return native (phi, chi) diagnostics on u, with r0 included in the cache key.
    
    This calls tools/pytrim.calc_phi_chi and does not reconstruct its output.
    """
    return python_result('chi', dict(screen.parameters, r0=r0), u,
                         lambda: python_modules()['cm_scatter'].calc_phi_chi(u, r0, screen.implementation))


def diagnostic_chi(u, r0, screening):
    """Evaluate the plot diagnostic using a supplied native screening callable.
    
    For u < 3e-4 use the Taylor limit phi(r0)-r0*phi'(r0); elsewhere use
    [phi(r0)-(1-u**2)*phi(r0/(1-u**2))]/u**2. Return (phi, chi).
    This adapts Numba screening for diagnostics, not its internal scattering API.
    """
    phi0, derivative = screening(r0)
    phi, _ = screening(r0/(1-u**2))
    return phi, np.where(u < 3e-4, phi0-r0*derivative, (phi0-phi*(1-u**2))/u**2)
