"""Infrastructure regressions: unavailable code must not erase/use stale data."""
import json
import subprocess

import numpy as np
import pytest

import reference_data as refs
import conftest as generation


def test_python_loader_rejects_alternative_sources(tmp_path, monkeypatch):
    """Reject another complete Python source tree before executing its code."""
    source = tmp_path / 'other-implementation'
    source.mkdir()
    for name in ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter'):
        (source / f'{name}.py').write_text('raise RuntimeError("wrong source executed")\n')
    monkeypatch.setattr(refs, 'PY_SOURCE', source)
    # Bypass the loader cache without changing modules used by other tests.
    with pytest.raises(pytest.fail.Exception, match='refusing alternative source directory'):
        refs.python_modules.__wrapped__()
    monkeypatch.setattr(refs, 'PY_SOURCE', tmp_path / 'missing')
    assert refs.python_modules.__wrapped__() is None


def test_python_cache_roundtrip_and_input_validation(tmp_path, monkeypatch):
    source = tmp_path / 'source'
    source.mkdir()
    for name in ('table1d', 'apsis', 'zbl', 'nlhlin', 'cm_scatter'):
        (source / f'{name}.py').write_text('# test source\n')
    monkeypatch.setattr(refs, 'PY_SOURCE', source)
    monkeypatch.setattr(refs, 'PY_CACHE', tmp_path / 'saved')
    monkeypatch.setattr(refs, 'python_modules', lambda: {})
    expected = refs.python_result('screen', {'model': 'ZBL'}, [0., 1.], lambda: [1., 0.5])
    monkeypatch.setattr(refs, 'python_modules', lambda: None)
    def unavailable():
        raise AssertionError('Cached mode must not call implementation')
    np.testing.assert_array_equal(
        refs.python_result('screen', {'model': 'ZBL'}, [0., 1.], unavailable), expected,
    )
    with pytest.raises(pytest.fail.Exception, match='no matching saved result'):
        refs.python_result('screen', {'model': 'ZBL'}, [0., 2.], unavailable)
    monkeypatch.setattr(refs, 'COEFFICIENT_HASH', 'different-input')
    with pytest.raises(pytest.fail.Exception, match='no matching saved result'):
        refs.python_result('screen', {'model': 'ZBL'}, [0., 1.], unavailable)


def make_fortran_cache(directory):
    directory.mkdir(parents=True)
    output = directory / 'sample.txt'
    output.write_text('0.0 1.0\n')
    (directory / 'screening_manifest.json').write_text(json.dumps({
        'coefficient_sha256': refs.COEFFICIENT_HASH,
        'outputs': {output.name: generation.digest(output)},
    }))
    return output


def test_missing_fortran_source_preserves_reference(tmp_path, monkeypatch):
    monkeypatch.setattr(generation, 'TESTS', tmp_path)
    monkeypatch.setenv('OPEN_SRIM_FORTRAN_SOURCE', str(tmp_path / 'missing'))
    directory = tmp_path / 'data/f90/screening'
    output = make_fortran_cache(directory)
    original = output.stat().st_mtime_ns
    assert generation.generate_fortran('screening', ['missing.f90']) == directory
    assert output.stat().st_mtime_ns == original
    output.write_text('0.0 999.0\n')
    with pytest.raises(pytest.fail.Exception, match='modified cached'):
        generation.generate_fortran('screening', ['missing.f90'])


def test_missing_fortran_and_missing_reference_fail(tmp_path, monkeypatch):
    monkeypatch.setattr(generation, 'TESTS', tmp_path)
    monkeypatch.setenv('OPEN_SRIM_FORTRAN_SOURCE', str(tmp_path / 'missing'))
    with pytest.raises(pytest.fail.Exception, match='no saved Fortran'):
        generation.generate_fortran('screening', ['missing.f90'])


def test_fortran_build_failure_does_not_fall_back_or_delete(tmp_path, monkeypatch):
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
        raise subprocess.CalledProcessError(1, args[0])
    monkeypatch.setattr(generation.subprocess, 'run', broken)
    output = make_fortran_cache(tmp_path / 'data/f90/screening')
    old = output.read_bytes()
    with pytest.raises(pytest.fail.Exception, match='not falling back'):
        generation.generate_fortran('screening', ['driver.f90'])
    assert output.read_bytes() == old
