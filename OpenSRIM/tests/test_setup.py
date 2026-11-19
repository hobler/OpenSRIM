""" Tests reference_data.py and its helpers for loading Python modules and caching results.
    
Its only purpose is to verify that the Python loader and cache logic work correctly, 
and that the Fortran reference data is available for comparison. It does not test the 
actual screening functions or their derivatives, which are tested in test_extended_element_pairs.py.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import reference_data as refs
import conftest as generation


def test_python_loader_uses_original_files_from_any_cwd(tmp_path, monkeypatch):
    """Verify that the Python loader uses the original source files from any current working directory.
    
    This test checks that the Python loader correctly loads the original source files from the specified
    source directory, regardless of the current working directory. It verifies that the loaded modules 
    match the expected source files and that the screening coefficients are correctly loaded from the 
    Fortran reference data. The test also checks that the screening function and its coefficients 
    are consistent with the expected values for a specific element pair (Z1=5, Z2=14).

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        monkeypatch (pytest.MonkeyPatch): A pytest fixture for modifying environment variables and attributes.

    Returns:
        None

    Raises:
        AssertionError: If the loaded modules do not match the expected source files or if the screening
                        coefficients do not match the expected values.
    """

    # Set the current working directory to the temporary path and load the Python modules from the original source files
    monkeypatch.chdir(tmp_path)
    names = ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter')
    previous = {name: sys.modules.get(name) for name in names}
    modules = refs.python_modules.__wrapped__()
    if modules is None:
        pytest.skip('Standalone Python sources unavailable')
    for name in names:
        assert Path(modules[name].__file__).resolve() == (
            refs.PYTHON_SOURCE_DIRECTORY / f'{name}.py').resolve()
        assert sys.modules.get(name) is previous[name]

    # Verify that the screening coefficients are correctly loaded from the Fortran reference data for a specific element pair (Z1=5, Z2=14)
    screen = modules['nlhlin'].NLHlin_screen(5, 14, 1.0)
    row = refs.coefficient_rows()[(5, 14)]
    np.testing.assert_array_equal(screen.a, row[[2, 4, 6]])
    np.testing.assert_array_equal(screen.b, row[[3, 5, 7]])
    assert screen.rmax == row[8]


def test_python_loader_rejects_alternative_sources(tmp_path, monkeypatch):
    """Verify that the Python loader rejects alternative source directories and raises an exception.

    This test checks that the Python loader correctly rejects alternative source directories and raises an exception 
    when an attempt is made to load modules from a different source directory. It verifies that the loader raises 
    a pytest.fail.Exception with the expected error message when an alternative source directory is provided.

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        monkeypatch (pytest.MonkeyPatch): A pytest fixture for modifying environment variables and attributes.

    Returns:
        None

    Raises:
        pytest.fail.Exception: If the Python loader does not raise an exception when an alternative source directory is provided.
    """

    # Create a temporary directory for the alternative source and populate it with dummy Python files that raise a RuntimeError when executed
    source = tmp_path / 'other-implementation'
    source.mkdir()
    for name in ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter'):
        (source / f'{name}.py').write_text('raise RuntimeError("wrong source executed")\n')

    # Set the PY_SOURCE attribute to the alternative source directory and verify that the 
    # Python loader raises an exception when attempting to load modules from it
    monkeypatch.setattr(refs, 'PY_SOURCE', source)
    with pytest.raises(pytest.fail.Exception, match='refusing alternative source directory'):
        refs.python_modules.__wrapped__()
    monkeypatch.setattr(refs, 'PY_SOURCE', tmp_path / 'missing')
    assert refs.python_modules.__wrapped__() is None


def test_python_cache_roundtrip_and_input_validation(tmp_path, monkeypatch):
    """Verify that the Python cache correctly saves and retrieves results, and that input validation works as expected.

    This test checks that the Python cache correctly saves and retrieves results for a specific screening function and its coefficients. 
    It verifies that the cached results match the expected values and that the loader raises an exception when an attempt 
    is made to retrieve results with different input parameters. The test also checks that the loader raises an exception
    when an attempt is made to retrieve results with a different coefficient hash.

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        monkeypatch (pytest.MonkeyPatch): A pytest fixture for modifying environment variables and attributes.

    Returns:
        None

    Raises:
        pytest.fail.Exception: If the Python cache does not correctly save and retrieve results, 
                               or if the loader does not raise an exception when input validation fails.
    """

    # Create a temporary directory for the source files and set the PY_SOURCE and PY_CACHE attributes
    # to the appropriate paths. Verify that the Python cache correctly saves and retrieves results for a
    # specific screening function and its coefficients, and that input validation works as expected
    source = tmp_path / 'source'
    source.mkdir()
    for name in ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter'):
        (source / f'{name}.py').write_text('# test source\n')
    monkeypatch.setattr(refs, 'PY_SOURCE', source)
    monkeypatch.setattr(refs, 'PY_CACHE', tmp_path / 'saved')
    monkeypatch.setattr(refs, 'python_modules', lambda: {})

    # Verify that the Python cache correctly saves and retrieves results for a specific screening function and its coefficients
    expected = refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0}, [0., 1.], lambda: [1., 0.5], norm="native", grid_name="test")
    monkeypatch.setattr(refs, 'python_modules', lambda: None)

    def unavailable():
        """Simulate a missing Python implementation by raising an exception when called.

        Raises:
            AssertionError: Always raised to simulate a missing Python implementation.
        """
        raise AssertionError('Cached mode must not call implementation')
    
    # Verify that the cached results match the expected values and that the loader raises an 
    # exception when an attempt is made to retrieve results with different input parameters
    np.testing.assert_array_equal(
        refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0}, [0., 1.], unavailable, norm="native", grid_name="test"), expected,
    )
    with pytest.raises(pytest.fail.Exception, match='metadata does not match inputs'):
        refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0}, [0., 2.], unavailable, norm="native", grid_name="test")
    monkeypatch.setattr(refs, 'COEFFICIENT_HASH', 'different-input')
    with pytest.raises(pytest.fail.Exception, match='metadata does not match inputs'):
        refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0}, [0., 1.], unavailable, norm="native", grid_name="test")

    # Integrity checks are independent of input matching: keep metadata unchanged
    # while altering the numerical results in a valid NPZ archive.
    path = next(refs.PY_CACHE.glob('*.npz'))
    manifest_path = refs.PY_CACHE / 'python_manifest.json'
    original_manifest = manifest_path.read_bytes()
    with np.load(path, allow_pickle=False) as saved:
        contents = {key: saved[key].copy() for key in saved.files}
    contents['result'] = contents['result'] + 1
    np.savez(path, **contents)
    with pytest.raises(pytest.fail.Exception, match='modified cached Python result'):
        refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0},
                           [0., 1.], unavailable, norm='native', grid_name='test')
    assert manifest_path.read_bytes() == original_manifest
    # An existing file without a manifest entry must also be rejected.
    manifest_path.write_text(json.dumps({'schema': 1, 'outputs': {}}))
    with pytest.raises(pytest.fail.Exception, match='Missing or modified cached Python result'):
        refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0},
                           [0., 1.], unavailable, norm='native', grid_name='test')
    for invalid in ('not json', '{"schema": 99, "outputs": {}}'):
        manifest_path.write_text(invalid)
        with pytest.raises(pytest.fail.Exception, match='invalid Python reference manifest'):
            refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0},
                               [0., 1.], unavailable, norm='native', grid_name='test')
    manifest_path.unlink()
    with pytest.raises(pytest.fail.Exception, match='Missing or invalid Python reference manifest'):
        refs.python_result('screen', {'model': 'ZBL', 'z1': 5, 'z2': 14, 'rnorm': 1.0},
                           [0., 1.], unavailable, norm='native', grid_name='test')


def make_fortran_cache(directory):
    """Create a Fortran cache for testing purposes.

    This function creates a Fortran cache in the specified directory for testing purposes. 
    It generates a sample output file and a manifest file containing the coefficient hash and output digest. 
    The function returns the path to the generated output file.

    Args:
        directory (Path): The directory where the Fortran cache will be created.

    Returns:
        Path: The path to the generated output file.
    """

    # Create the specified directory and generate a sample output file and manifest file for the Fortran cache
    directory.mkdir(parents=True)
    output = directory / 'sample.txt'
    output.write_text('0.0 1.0\n')
    (directory / 'screening_manifest.json').write_text(json.dumps({
        'coefficient_sha256': refs.COEFFICIENT_HASH,
        'outputs': {output.name: generation.digest(output)},
    }))
    return output


def test_missing_fortran_source_preserves_reference(tmp_path, monkeypatch):
    """Verify that the Fortran cache preserves reference data when the source is missing.

    This test checks that the Fortran cache correctly preserves reference data when the source files are missing.
    It verifies that the cached output file remains unchanged and that the loader raises an exception when an attempt
    is made to generate the Fortran cache with a missing source file. The test also checks that the loader raises
    an exception when an attempt is made to generate the Fortran cache with a modified cached output file.

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        monkeypatch (pytest.MonkeyPatch): A pytest fixture for modifying environment variables and attributes.

    Returns:
        None

    Raises:
        pytest.fail.Exception: If the Fortran cache does not correctly preserve reference data,
                               or if the loader does not raise an exception when the source is missing or the cached output is modified.
    """

    # Set the TESTS attribute to the temporary path and the OPEN_SRIM_FORTRAN_SOURCE environment variable to a missing source directory
    monkeypatch.setattr(generation, 'TESTS', tmp_path)
    monkeypatch.setenv('OPEN_SRIM_FORTRAN_SOURCE', str(tmp_path / 'missing'))
    directory = tmp_path / 'data/f90/screening'
    output = make_fortran_cache(directory)
    original = output.stat().st_mtime_ns

    # Verify that the Fortran cache preserves reference data when the source is 
    # missing and raises an exception when attempting to generate the cache with a missing source file
    assert generation.generate_fortran('screening', ['missing.f90']) == directory
    assert output.stat().st_mtime_ns == original
    output.write_text('0.0 999.0\n')
    with pytest.raises(pytest.fail.Exception, match='modified cached'):
        generation.generate_fortran('screening', ['missing.f90'])


def test_missing_fortran_and_missing_reference_fail(tmp_path, monkeypatch):
    """Verify that the Fortran cache raises an exception when both the source and reference data are missing.

    This test checks that the Fortran cache correctly raises an exception when both the source files and reference data are missing.
    It verifies that the loader raises a pytest.fail.Exception with the expected error message when an attempt is made to generate 
    the Fortran cache with a missing source file and no reference data available.

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        monkeypatch (pytest.MonkeyPatch): A pytest fixture for modifying environment variables and attributes.

    Returns:
        None

    Raises:
        pytest.fail.Exception: If the Fortran cache does not raise an exception when both the source and reference data are missing.
    """

    # Set the TESTS attribute to the temporary path and the OPEN_SRIM_FORTRAN_SOURCE environment variable to a missing source directory
    monkeypatch.setattr(generation, 'TESTS', tmp_path)
    monkeypatch.setenv('OPEN_SRIM_FORTRAN_SOURCE', str(tmp_path / 'missing'))

    # Verify that the Fortran cache raises an exception when both the source files and reference data are missing
    with pytest.raises(pytest.fail.Exception, match='no saved Fortran'):
        generation.generate_fortran('screening', ['missing.f90'])


@pytest.mark.parametrize('group', ['cm', 'extended'])
def test_fortran_rejects_old_diagnostic_or_mehler_cache(tmp_path, group):
    """ Verify that the Fortran cache rejects old diagnostic or Mehler cache files and raises an exception.

    This test checks that the Fortran cache correctly rejects old diagnostic or Mehler cache files and raises an exception
    when an attempt is made to generate the Fortran cache with an outdated manifest file. It verifies that the loader raises
    a pytest.fail.Exception with the expected error message when an outdated manifest file is provided.

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        group (str): The name of the group for which the Fortran cache is being generated (e.g., 'cm' or 'extended').

    Returns:
        None

    Raises:
        pytest.fail.Exception: If the Fortran cache does not raise an exception when an outdated manifest file is provided.
    """

    # Create a temporary directory for the Fortran cache and generate a sample output file and manifest file
    output = tmp_path / 'sample.txt'
    output.write_text('0.0 1.0\n')
    manifest = tmp_path / f'{group}_manifest.json'

    # Verify that the Fortran cache rejects old diagnostic or Mehler cache files and raises an exception when an outdated manifest file is provided
    metadata = {'coefficient_sha256': refs.COEFFICIENT_HASH,
                'outputs': {output.name: generation.digest(output)}}
    manifest.write_text(json.dumps(metadata))
    with pytest.raises(pytest.fail.Exception, match='predate production chi/Legendre'):
        generation.cached_fortran(tmp_path, group, 'test')

    # Update the manifest file with the correct calculation API and verify that the Fortran cache correctly generates the cache for the specified group
    metadata['calculation_api'] = 'production-chi-legendre-v1'
    manifest.write_text(json.dumps(metadata))
    assert generation.cached_fortran(tmp_path, group, 'test') == tmp_path


def test_fortran_build_failure_does_not_fall_back_or_delete(tmp_path, monkeypatch):
    """Verify that the Fortran build failure does not fall back or delete the cache.

    This test checks that the Fortran build failure does not fall back to an alternative implementation or delete the existing cache.
    It verifies that the loader raises a pytest.fail.Exception with the expected error message when an attempt is made to generate the 
    Fortran cache and the build fails. The test also checks that the existing cache remains unchanged after the build failure.

    Args:
        tmp_path (Path): A temporary directory provided by pytest for creating test files.
        monkeypatch (pytest.MonkeyPatch): A pytest fixture for modifying environment variables and attributes.

    Returns:
        None

    Raises:
        pytest.fail.Exception: If the Fortran build failure falls back to an alternative implementation or deletes the existing cache.
    """

    # Set the TESTS attribute to the temporary path, create a temporary directory for the Fortran source files, and simulate a compiler failure
    monkeypatch.setattr(generation, 'TESTS', tmp_path)
    root = tmp_path / 'implementation'
    (root / 'code').mkdir(parents=True)
    (root / 'code/srcs.mk').write_text('${SRCDIR}/solver.f90\n')
    (root / 'code/solver.f90').write_text('! source\n')
    (root / 'driver.f90').write_text('! driver\n')
    monkeypatch.setenv('OPEN_SRIM_FORTRAN_SOURCE', str(root))
    monkeypatch.delenv('OPEN_SRIM_FORTRAN_MODE', raising=False)
    monkeypatch.setattr(generation.shutil, 'which', lambda command: '/fake/gfortran')

    def broken(*args, **kwargs):
        """Simulate a broken subprocess.run that raises a CalledProcessError to simulate a compiler failure.

        Args:
            *args: Positional arguments passed to subprocess.run.
            **kwargs: Keyword arguments passed to subprocess.run.

        Raises:
            subprocess.CalledProcessError: Always raised to simulate a compiler failure.
        """
        raise subprocess.CalledProcessError(1, args[0])

    # Verify that the Fortran build failure does not fall back to an alternative implementation or delete the cache, 
    # and raises an exception when attempting to generate the Fortran cache
    monkeypatch.setattr(generation.subprocess, 'run', broken)
    output = make_fortran_cache(tmp_path / 'data/f90/screening')
    old = output.read_bytes()
    with pytest.raises(pytest.fail.Exception, match='not falling back'):
        generation.generate_fortran('screening', ['driver.f90'])
    assert output.read_bytes() == old

