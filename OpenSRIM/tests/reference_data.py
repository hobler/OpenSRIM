"""Shared inputs and persistent, input-validated Python comparison results."""
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
_installation = tempfile.TemporaryDirectory(prefix='opensrim-python-reference-')


@lru_cache(maxsize=1)
def python_modules():
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
    shutil.copy2(COEFFICIENTS, directory / 'dmol_coeffs_rmax.dat')
    modules = {}
    with pytest.MonkeyPatch.context() as aliases:
        for name in names:
            spec = importlib.util.spec_from_file_location(f'_comparison_{name}', directory / f'{name}.py')
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            modules[name] = module
            aliases.setitem(sys.modules, name, module)
    return modules


def python_result(kind, parameters, grid, compute):
    """Regenerate live results, or require exactly matching saved inputs."""
    grid = np.asarray(grid, dtype=np.float64)
    metadata = {'schema': 1, 'kind': kind, 'parameters': parameters,
                'coefficient_sha256': COEFFICIENT_HASH, 'grid': grid.tolist()}
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
    return {(int(row[0]), int(row[1])): row for row in np.loadtxt(COEFFICIENTS)}


class PythonScreen:
    """Input description usable even without the standalone Python sources."""
    def __init__(self, model, z1, z2, rnorm):
        self.parameters = {'model': model, 'z1': z1, 'z2': z2, 'rnorm': rnorm}
        modules = python_modules()
        self.implementation = None
        if modules is not None:
            cls = modules['zbl'].ZBL_screen if model == 'ZBL' else modules['nlhlin'].NLHlin_screen
            self.implementation = cls(z1, z2, rnorm)
        self.rmax = (np.inf if model == 'ZBL' else
                     coefficient_rows()[tuple(sorted((z1, z2)))][8] / rnorm)

    def __call__(self, radii):
        return python_result('screen', self.parameters, radii,
                             lambda: self.implementation(radii))


def python_scatter(e, p_values, screen, n_absc):
    def compute():
        module = python_modules()['cm_scatter']
        module.setup(n_absc)
        return np.asarray([module.scatter_integrals(e, p, screen.implementation)
                           for p in p_values]).T
    return python_result('scatter', dict(screen.parameters, e=float(e), n_absc=n_absc), p_values, compute)


def python_chi(u, r0, screen):
    return python_result('chi', dict(screen.parameters, r0=r0), u,
                         lambda: python_modules()['cm_scatter'].calc_phi_chi(u, r0, screen.implementation))


def diagnostic_chi(u, r0, screening):
    """Plot diagnostic for Numba screening; not a new production solver."""
    phi0, derivative = screening(r0)
    phi, _ = screening(r0/(1-u**2))
    return phi, np.where(u < 3e-4, phi0-r0*derivative, (phi0-phi*(1-u**2))/u**2)
