"""Provide session-scoped native Fortran references for comparison tests.

Build complete source trees in temporary working directories. When code or the
compiler is unavailable, validate preserved results instead; a failed live build
or run remains an error. Numerical test modules validate the requested grids.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import pytest

from reference_data import TESTS, ROOT, COEFFICIENTS, COEFFICIENT_HASH


def digest(path):
    """Return the SHA256 of file contents for reference integrity/provenance.
    
    Timestamps are not used: identical bytes have the same digest after copying.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cached_fortran(directory, group, reason):
    """Validate an existing reference group and return its directory.
    
    With a manifest, require the current coefficient hash and every recorded
    output digest. Stored source hashes document provenance but are not compared
    with an expected code revision here. Legacy text files without a manifest
    remain supported; individual tests must validate their filenames and grids.
    Missing or modified references fail instead of skipping comparisons.
    """
    manifest = directory / f'{group}_manifest.json'
    if manifest.is_file():
        metadata = json.loads(manifest.read_text())
        if metadata['coefficient_sha256'] != COEFFICIENT_HASH:
            pytest.fail(f'Fortran references use different coefficients: {manifest}')
        for name, expected in metadata['outputs'].items():
            path = directory / name
            if not path.is_file() or digest(path) != expected:
                pytest.fail(f'Missing or modified cached Fortran result: {path}')
    elif not any(directory.glob('*.txt')):
        pytest.fail(f'{reason}; no saved Fortran data in {directory}')
    else:
        print(f'Legacy Fortran files without provenance: {directory}; grids are validated by tests')
    print(f'Fortran {group}: cached results ({reason})')
    return directory


def generate_fortran(group, drivers, common=(), runs=None):
    """Build native drivers and publish references, or validate saved results.
    
    Args:
        group: Output group; cm uses data/f90 directly, other groups a subfolder.
        drivers: Driver paths relative to the IMSIL root, or absolute paths.
        common: Shared Fortran source paths relative to the IMSIL root.
        runs: Argument tuples passed to each executable; defaults to one run.
    
    Sources are read from code/srcs.mk, excluding the application main program.
    Missing sources/compiler/drivers permit cached mode; compilation or execution
    errors do not. Required runtime tables must also exist for a live run.
    Builds run in temporary storage. Only successful runs publish output files
    and a manifest; replacement is atomic per file, not for the whole group.
    """
    imsil = Path(os.environ.get('OPEN_SRIM_FORTRAN_SOURCE', ROOT / 'imsil-omp'))
    directory = TESTS / 'data/f90' / (group if group != 'cm' else '')
    src_list = imsil / 'code/srcs.mk'
    # Explicit mode also permits exercising fallback without removing source.
    if os.environ.get('OPEN_SRIM_FORTRAN_MODE') == 'cached':
        return cached_fortran(directory, group, 'explicit cached mode')
    if not src_list.is_file() or shutil.which('gfortran') is None:
        return cached_fortran(directory, group, 'source tree or gfortran unavailable')
    # srcs.mk is the source inventory, not merely a marker that a directory exists.
    # Exclude main.f90 because each reference driver supplies its own program.
    sources = []
    for line in src_list.read_text().splitlines():
        entry = line.strip().rstrip('\\').strip()
        if entry.startswith('${SRCDIR}/') and not entry.endswith('/main.f90'):
            sources.append(imsil / entry.replace('${SRCDIR}', 'code'))
    # The screen/extended drivers live under tests; classic drivers live in IMSIL.
    driver_paths = [imsil / d if not Path(d).is_absolute() else Path(d) for d in drivers]
    common_paths = [imsil / d for d in common]
    if not sources or not all(p.is_file() for p in sources + driver_paths + common_paths):
        return cached_fortran(directory, group, 'Fortran implementation or driver incomplete')
    directory.mkdir(parents=True, exist_ok=True)
    # A failed build/run never destroys previously saved outputs.
    with tempfile.TemporaryDirectory(prefix=f'opensrim-{group}-') as tmp:
        stage = Path(tmp)
        (stage / 'mod').mkdir()
        tables = imsil / 'tables'
        required = ('gauss_legendre.dat', 'rhoe.dat', 'rhozbl.dat', 'scattab', 'screentab') if group in ('cm', 'extended') else ()
        for name in required:
            src = tables / name
            if src.is_dir():
                (stage / name).symlink_to(src.resolve(), target_is_directory=True)
            elif src.is_file():
                shutil.copy2(src, stage / name)
            else:
                pytest.fail(f'Required Fortran input missing: {src}')
        (stage / 'screencoefs').mkdir()
        # Only a format header differs from Numba; numerical rows are untouched.
        (stage / 'screencoefs/NLHLIN_python_coeffs.dat').write_text('# nexp=3\n' + COEFFICIENTS.read_text())
        for driver in driver_paths:
            executable = stage / driver.stem
            build_log = directory / f'build_{driver.stem}.log'
            run_log = directory / f'run_{driver.stem}.log'
            cmd = ['gfortran', '-g', '-O0', '-fcheck=all', '-fbacktrace', '-Wall',
                   '-J', str(stage / 'mod'), '-cpp', '-o', str(executable)]
            try:
                with build_log.open('wb') as log:
                    subprocess.run(cmd + [str(p) for p in sources + common_paths + [driver]],
                                   check=True, stdout=log, stderr=subprocess.STDOUT, timeout=600)
                with run_log.open('wb') as log:
                    for args in (runs or [()]):
                        subprocess.run([str(executable), *map(str, args)], cwd=stage, check=True, stdout=log,
                                       stderr=subprocess.STDOUT, timeout=600)
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                pytest.fail(f'Fortran {group} build/run failed, not falling back to stale data: '
                            f'{exc}; logs: {build_log}, {run_log}')
        outputs = list(stage.glob('*.txt'))
        if not outputs:
            pytest.fail(f'Fortran {group} generated no result files')
        # Publish only after all requested runs succeed. Each replacement is atomic;
        # an interrupted multi-file publication is detected by manifest hashes.
        for result in outputs:
            target = directory / result.name
            tmp_target = target.with_suffix('.tmp')
            shutil.copy2(result, tmp_target)
            tmp_target.replace(target)
        metadata = {'schema': 1, 'coefficient_sha256': COEFFICIENT_HASH,
                    'runs': [list(args) for args in (runs or [()])],
                    'sources': {str(p): digest(p) for p in sources + common_paths + driver_paths},
                    'outputs': {p.name: digest(p) for p in outputs}}
        manifest = directory / f'{group}_manifest.json'
        temporary = manifest.with_suffix('.tmp')
        temporary.write_text(json.dumps(metadata, indent=2))
        temporary.replace(manifest)
    print(f'Fortran {group}: generated {len(outputs)} files with shared NLHlin_3eV coefficients')
    return directory


@pytest.fixture(scope='session')
def build_and_run_f90():
    """Prepare classic B/Si ZBL and NLHlin references once per pytest session.
    
    Both drivers link the shared comparison-output module and use native solvers.
    """
    base = 'code/unittests/sn/cmscatter/'
    return generate_fortran('cm', [base+'test_cm_scatter_f90.f90', base+'test_cm_scatter_nlhlin_f90.f90'],
                            [base+'cm_scatter_comparison_data.f90'])


@pytest.fixture(scope='session')
def screening_fortran():
    """Prepare native phi/derivative grids for the screening tests once per session.
    
    The driver covers two models, three configurations and 92 pairs each.
    """
    return generate_fortran('screening', [TESTS / 'fortran/screen_reference.f90'])


@pytest.fixture(scope='session')
def extended_fortran():
    """Prepare pair, boundary and quadrature-convergence references.
    
    Each model/pair/order runs in a fresh process so native allocation and saved
    module state cannot leak between configurations. Only B/Si uses every Mehler
    order; the other element pairs use order 512.
    """
    from scattering_support import EXTENDED_PAIRS, MEHLER_ORDERS
    runs = [(model, z1, z2, order)
            for model in ('ZBL', 'NLHlin')
            for z1, z2 in EXTENDED_PAIRS
            for order in (MEHLER_ORDERS if (z1, z2) == (5, 14) else (512,))]
    return generate_fortran('extended', [TESTS / 'fortran/scatter_extended.f90'], runs=runs)
