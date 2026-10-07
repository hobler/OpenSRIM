"""Three-way screening diagnostics corresponding to the original plot grids.

All backends read the same NLHlin_3eV coefficients. Radius/derivative conversions
below only express the APIs in common units; no output fitting is performed.
"""
import numpy as np
import pytest

import scattering_support as cm
from reference_data import PythonScreen

# Normalization exponents (p1, p2), with an optional fixed second element Z2=29.
CONFIGURATIONS = ((0., 0., None), (0.23, 1., None), (0.5, 2/3, 29))
SCREEN_TOLERANCES = {'screen': {'numba': (1e-10, 1e-11), 'fortran': (1e-10, 1e-11)}}


def compare_screening(model, config_index, directory):
    """Evaluate native screening APIs and report all pairwise maximum errors."""
    p1, p2, fixed_z2 = CONFIGURATIONS[config_index]
    summaries, failures = {}, []
    derivative_cutoff_points = 0
    for z1 in range(1, 93):
        z2 = fixed_z2 or z1
        rnorm = 1. if config_index == 0 else 0.4685/(z1**p1+z2**p1)**p2
        # Python comes from tools/pytrim; the Numba loader comes from simulators/opentrim.
        native_norm, coefs = cm.numba_coefficients(model, z1, z2)
        if model == 'NLHlin':
            grid = np.linspace(0., 3./rnorm, 101)
            screen = PythonScreen(model, z1, z2, rnorm)
            python = screen(grid)
            factor = rnorm/native_norm
        else:
            grid = np.linspace(0., 3., 101)  # physical radius in Angstrom
            screen = PythonScreen(model, z1, z2, native_norm)
            python = screen(grid/native_norm)
            python[1] /= native_norm  # derivative with respect to Angstrom
            factor = 1./native_norm
        numba = np.asarray(cm.evaluate_numba_screen(grid*factor, model, coefs))
        numba[1] *= factor
        path = directory / f'screen_{model.lower()}_c{config_index+1}_z{z1}.txt'
        if not path.is_file():
            pytest.fail(f'Missing Fortran screening reference: {path}')
        data = np.loadtxt(path, ndmin=2)
        assert data.shape == (len(grid), 3)
        cm.assert_arrays_equal(grid, data[:, 0], rtol=1e-14, atol=1e-14)
        # NLHlin is continuous at rmax, but its derivative jumps there.
        # Unit conversions may place a mathematically identical cutoff one ULP
        # to either side. Compare phi everywhere, derivatives only where defined.
        at_cutoff = (np.isclose(grid, screen.rmax, rtol=1e-14, atol=0)
                     if model == 'NLHlin' else np.zeros(grid.size, dtype=bool))
        derivative_cutoff_points += int(np.count_nonzero(at_cutoff))
        for i, quantity in enumerate(('phi', 'dphi')):
            selected = np.ones(grid.size, dtype=bool) if i == 0 else ~at_cutoff
            cm.compare_three_way(
                summaries, failures,
                dict(Python=python[i, selected], Numba=numba[i, selected], Fortran=data[selected, i+1]),
                grid[selected], f'Z1={z1}, Z2={z2}, config={config_index+1}', 'screen',
                tolerances=SCREEN_TOLERANCES, prefix=f'{quantity}: ', coordinate_name='r',
            )
        # Check the finite-range contract as well as pairwise agreement.
        if model == 'NLHlin':
            tail = (grid > screen.rmax) & ~at_cutoff
            for label, values in (('Python', python), ('Numba', numba), ('Fortran', data[:, 1:].T)):
                if not np.all(values[:, tail] == 0):
                    failures.append(f'{label} nonzero screening beyond cutoff: Z1={z1}, Z2={z2}')
    cm.finish_comparison(f'{model} screening config={config_index+1}', summaries, failures,
                         counts={'undefined derivative points at cutoff': derivative_cutoff_points})


@pytest.mark.parametrize('config_index', range(3))
def test_plot_screen_data_python_numba_fortran(config_index, screening_fortran):
    """Compare NLHlin screening and its derivative across all three backends.

    Each normalization covers 92 element pairs and 101 radii per pair.
    Both outputs must vanish exactly beyond the cutoff."""
    # Only the discontinuous derivative is omitted at the cutoff; phi is checked everywhere.
    compare_screening('NLHlin', config_index, screening_fortran)


@pytest.mark.parametrize('config_index', range(3))
def test_plot_zblscreen_data_python_numba_fortran(config_index, screening_fortran):
    """Compare ZBL screening and its derivative across all three backends.

    Each configuration covers 92 element pairs, including Z=92. Radii
    and derivatives are converted to common physical units."""
    # ZBL has no finite cutoff, so no derivative points are excluded.
    compare_screening('ZBL', config_index, screening_fortran)
