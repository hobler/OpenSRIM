"""Three-way screening diagnostics corresponding to the original plot grids.

    The following tests compare the Python, Numba, and Fortran screening implementations
    for the NLHlin and ZBL models across three normalization configurations. Each configuration
    covers 92 element pairs and 101 radii per pair. The tests validate that the screening
    functions and their derivatives are consistent across all three backends, and that they
    vanish exactly beyond the cutoff for the NLHlin model. The tests also check that the
    derivatives are correctly computed and that the Fortran reference data is available for
    comparison. Any failures encountered during the comparisons are reported at the end of
    each test.
    
    - compare_screening: Helper function that performs the actual comparisons for a given model, 
                         configuration index, and directory containing the Fortran reference data. It checks the screening 
                         functions and their derivatives, validates the Fortran input columns, and reports any failures 
                         encountered during the comparisons.
    
    The following tests are included:

    - test_nlhlin_screening_python_numba_fortran: Compares NLHlin screening and its derivative across all three backends for each normalization configuration.
    - test_zbl_screening_python_numba_fortran: Compares ZBL screening and its derivative across all three backends for each normalization configuration, 
                                               including a configuration with a fixed partner element (Z2=29).

   

"""
import numpy as np
import pytest

import scattering_support as cm
from reference_data import PythonScreen, reference_filename

# Normalization configurations: (p1, p2, fixed_Z2)
CONFIGURATIONS = ((0., 0., None), (0.23, 1., None), (0.5, 2/3, 29))

# Tolerances for screening comparisons: (relative tolerance, absolute tolerance)
SCREEN_TOLERANCES = {'screen': {'numba': (1e-12, 1e-13), 'fortran': (1e-12, 1e-13)}}


def compare_screening(model, config_index, directory):
    """Compare screening functions and their derivatives across Python, Numba, and Fortran backends.

    This function performs a three-way comparison of the screening functions and their derivatives for 
    the specified model and normalization configuration. It checks that the Fortran reference data is available, 
    validates the Fortran input columns, and reports any failures encountered during the comparisons. The function also
    checks that the screening functions vanish exactly beyond the cutoff for the NLHlin model and that the derivatives 
    are correctly computed. The comparisons are performed for all 92 element pairs and 101 radii per pair, and the 
    results are summarized at the end of the function.

    Args:
        model (str): The screening model to compare ('NLHlin' or 'ZBL').
        config_index (int): The index of the normalization configuration to use (0, 1, or 2).
        directory (Path): The directory containing the Fortran reference data for the specified model and configuration.

    Returns:
        None

    Raises:
        pytest.fail: If the Fortran reference data is missing or if any failures are encountered during the comparisons.
    """

    # Unpack the normalization configuration parameters for the specified index
    p1, p2, fixed_z2 = CONFIGURATIONS[config_index]
    summaries, failures = {}, []
    derivative_cutoff_points = 0

    # Iterate over all element pairs (Z1, Z2) and perform the comparisons for the specified model and configuration
    for z1 in range(1, 93):
        z2 = fixed_z2 or z1

        # Compute the normalization factor based on the configuration parameters and element charges
        rnorm = 1. if config_index == 0 else 0.4685/(z1**p1+z2**p1)**p2

        # Compute the native normalization and coefficients for the specified model and element pair
        native_norm, coefs = cm.numba_coefficients(model, z1, z2)
        if model == 'NLHlin':
            grid = np.linspace(0., 3./rnorm, 101)
            screen = PythonScreen(model, z1, z2, rnorm, norm=f"c{config_index+1}")
            python = screen(grid, grid_name="screening")
            factor = rnorm/native_norm
        else:
            grid = np.linspace(0., 3., 101)  # physical radius in Angstrom
            screen = PythonScreen(model, z1, z2, native_norm, norm=f"c{config_index+1}")
            python = screen(grid/native_norm, grid_name="screening")
            python[1] /= native_norm  # derivative with respect to Angstrom
            factor = 1./native_norm

        # Compute the screening function and its derivative using the Numba backend for 
        # the specified model, element pair, and normalization factor
        numba = np.asarray(cm.evaluate_numba_screen(grid*factor, model, coefs))
        numba[1] *= factor

        # Check that the Fortran reference data is available for the specified model, element pair, and configuration
        path = directory / reference_filename('screen', screen.parameters, norm=screen.norm,
                                              grid_name='screening', extension='txt')
        if not path.is_file():
            pytest.fail(f'Missing Fortran screening reference: {path}')
        data = np.loadtxt(path, ndmin=2)

        # Validate that the Fortran reference data has the expected number of columns (3) for the screening function and its derivative
        assert data.shape == (len(grid), 3)
        cm.assert_arrays_equal(grid, data[:, 0], rtol=1e-14, atol=1e-14)

        # Check for points at the cutoff radius where the derivative is undefined for the NLHlin model
        at_cutoff = (np.isclose(grid, screen.rmax, rtol=1e-14, atol=0)
                     if model == 'NLHlin' else np.zeros(grid.size, dtype=bool))
        derivative_cutoff_points += int(np.count_nonzero(at_cutoff))

        # Perform a three-way comparison of the screening function and its derivative across Python, Numba, and Fortran backends
        for i, quantity in enumerate(('phi', 'dphi')):
            selected = np.ones(grid.size, dtype=bool) if i == 0 else ~at_cutoff

            # Compare the screening function and its derivative across Python, Numba, and Fortran 
            # backends for the specified quantity (phi or dphi) and element pair (Z1, Z2)
            cm.compare_three_way(
                summaries, failures,
                dict(Python=python[i, selected], Numba=numba[i, selected], Fortran=data[selected, i+1]),
                grid[selected], f'Z1={z1}, Z2={z2}, config={config_index+1}', 'screen',
                tolerances=SCREEN_TOLERANCES, prefix=f'{quantity}: ', coordinate_name='r',
            )

        # For the NLHlin model, check that the screening function and its derivative vanish exactly beyond 
        # the cutoff radius, and report any failures encountered during the comparisons
        if model == 'NLHlin':
            tail = (grid > screen.rmax) & ~at_cutoff
            for label, values in (('Python', python), ('Numba', numba), ('Fortran', data[:, 1:].T)):
                if not np.all(values[:, tail] == 0):
                    failures.append(f'{label} nonzero screening beyond cutoff: Z1={z1}, Z2={z2}')

    # Finalize the comparison and report any failures encountered during the test for the specified model, configuration, and element pairs
    cm.finish_comparison(f'{model} screening config={config_index+1}', summaries, failures,
                         counts={'undefined derivative points at cutoff': derivative_cutoff_points})


@pytest.mark.parametrize('config_index', range(3))
def test_nlhlin_screening_python_numba_fortran(config_index, screening_fortran):
    """Compare NLHlin screening and its derivative across all three backends.

    Each configuration covers 92 element pairs, including Z=92. Radii
    and derivatives are converted to common physical units.

    config_index:
                    1: p1=0, p2=0, no fixed Z2
                    2: p1=0.23, p2=1, no fixed Z2
                    3: p1=0.5, p2=2/3, fixed Z2=29 (copper)

    Args:
        config_index (int): The index of the normalization configuration to use (0, 1, or 2).
        screening_fortran (Path): The directory containing the Fortran reference data for the NLHlin model and the specified configuration.

    Returns:
        None
    """

    # Compare the NLHlin screening function and its derivative across Python, Numba, and 
    # Fortran backends for the specified normalization configuration
    compare_screening('NLHlin', config_index, screening_fortran)


@pytest.mark.parametrize('config_index', range(3))
def test_zbl_screening_python_numba_fortran(config_index, screening_fortran):
    """Compare ZBL screening and its derivative across all three backends.

    Each configuration covers 92 element pairs, including Z=92. Radii
    and derivatives are converted to common physical units.

    config_index:
                    1: p1=0, p2=0, no fixed Z2
                    2: p1=0.23, p2=1, no fixed Z2
                    3: p1=0.5, p2=2/3, fixed Z2=29 (copper)

    Args:
        config_index (int): The index of the normalization configuration to use (0, 1, or 2).
        screening_fortran (Path): The directory containing the Fortran reference data for the ZBL model and the specified configuration.

    Returns:
        None
    """

    # Compare the ZBL screening function and its derivative across Python, Numba, and 
    # Fortran backends for the specified normalization configuration
    compare_screening('ZBL', config_index, screening_fortran)
