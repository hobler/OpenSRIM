""" Checks if Fortran source is available and builds reference outputs for comparison 
    if not available reference data is used

    - digest()
        calculates the SHA256 of a file for provenance and integrity checks
    - cached_fortran()
        checks if reference data is available and valid, returns the directory of the reference data
    - generate_fortran()
        builds and runs Fortran drivers to generate reference data, or validates existing reference data
    - build_and_run_f90()
        pytest fixture that prepares ZBL and NLHlin references once per session
    - screening_fortran()
        pytest fixture that prepares screening references for screening tests once per session
    - extended_fortran()
        pytest fixture that prepares different pair and boundary references for extended tests once per session
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
    
    Args:
        path: The file path to calculate the digest for.

    Returns:
        The SHA256 digest as a hexadecimal string.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cached_fortran(directory, group, reason):
    """Check if reference data is available and valid, return the directory of the reference data.

    Args:
        directory: The directory where the reference data is stored.
        group: The output group name (e.g., 'cm', 'screening', 'extended').
        reason: The reason for using cached data (e.g., 'explicit cached mode', 'source tree or gfortran unavailable').

    Returns:
        The directory containing the cached Fortran reference data.
    """

    # sets the path to the manifest file for the given group
    manifest = directory / f'{group}_manifest.json'
    if manifest.is_file():
        metadata = json.loads(manifest.read_text())
        if group in ('cm', 'extended') and metadata.get('calculation_api') != 'production-chi-legendre-v1':
            pytest.fail(f'Fortran references predate production chi/Legendre comparison: {manifest}; regenerate with sources')

        # validates the coefficient hash in the manifest against the expected hash
        if metadata['coefficient_sha256'] != COEFFICIENT_HASH:
            pytest.fail(f'Fortran references use different coefficients: {manifest}')

        # goes through each output file listed in the manifest and checks if it exists and matches the expected hash
        for name, expected in metadata['outputs'].items():
            path = directory / name
            if not path.is_file() or digest(path) != expected:
                pytest.fail(f'Missing or modified cached Fortran result: {path}')
    elif group in ('cm', 'extended'):
        pytest.fail(f'Fortran {group} requires a production chi/Legendre manifest: {manifest}')
    elif not any(directory.glob('*.txt')):
        pytest.fail(f'{reason}; no saved Fortran data in {directory}')
    else:
        print(f'Legacy Fortran files without provenance: {directory}; grids are validated by tests')
    print(f'Fortran {group}: cached results ({reason})')
    return directory


def generate_fortran(group, drivers, common=(), runs=None):
    """Generate reference data by building and running Fortran drivers, or validate existing reference data.

    Args:
        group: The output group name (e.g., 'cm', 'screening', 'extended').
        drivers: A list of Fortran driver source files to compile and run.
        common: A list of common Fortran source files to include in the build (e.g., shared modules)
        runs: A list of argument tuples to pass to the drivers when running (e.g., [(arg1, arg2), ...]). If None, the drivers are run without arguments.

    Returns:
        The directory containing the cached Fortran reference data.
    """

    # get the path to the Fortran source code, either from an environment variable or a default location
    imsil = Path(os.environ.get('OPEN_SRIM_FORTRAN_SOURCE', ROOT / 'imsil-omp'))

    # The 'cm' group uses the base directory, while other groups get subdirectories.
    directory = TESTS / 'data/f90' / (group if group != 'cm' else '')

    # Define the path to the Makefile source list containing the Fortran files.
    src_list = imsil / 'code/srcs.mk'

    # Check whether cached mode was explicitly requested.
    # In this mode, existing reference data is used instead of compiling Fortran.
    if os.environ.get('OPEN_SRIM_FORTRAN_MODE') == 'cached':
        return cached_fortran(directory, group, 'explicit cached mode')

    # Check whether the source list exists and whether the GNU Fortran compiler
    # (gfortran) is available on the system.
    # If either is missing, use the existing cached reference data instead.
    if not src_list.is_file() or shutil.which('gfortran') is None:
        return cached_fortran(directory, group, 'source tree or gfortran unavailable')

    # read the source list and collect the Fortran source files, excluding the main driver files
    sources = []
    for line in src_list.read_text().splitlines():
        entry = line.strip().rstrip('\\').strip()
        if entry.startswith('${SRCDIR}/') and not entry.endswith('/main.f90'):
            sources.append(imsil / entry.replace('${SRCDIR}', 'code'))

    # Create the complete paths to the Fortran driver files.
    # Absolute paths are used directly; relative paths are interpreted
    # relative to the Fortran source directory.
    driver_paths = [imsil / d if not Path(d).is_absolute() else Path(d) for d in drivers]
    # Create the complete paths to additional shared Fortran source files.
    common_paths = [imsil / d for d in common]

    # Check whether source files were found and whether all required files exist.
    # If the implementation or any driver is incomplete, use cached reference data.
    if not sources or not all(p.is_file() for p in sources + driver_paths + common_paths):
        return cached_fortran(directory, group, 'Fortran implementation or driver incomplete')

    # Create the output directory if it does not already exist.
    # exist_ok=True prevents an error if the directory already exists.
    directory.mkdir(parents=True, exist_ok=True)

    # Create a temporary directory for compilation and execution.
    # The directory is automatically deleted when the with-block ends.
    # Previously saved reference data remains untouched during this process.
    with tempfile.TemporaryDirectory(prefix=f'opensrim-{group}-') as tmp:
        stage = Path(tmp)
        (stage / 'mod').mkdir()
        tables = imsil / 'tables'

        # Define the input files required for the 'cm' and 'extended' groups.
        # Other groups do not require these additional table files.
        required = ('gauss_legendre.dat',) if group in ('cm', 'extended') else ()
        for name in required:
            src = tables / name

            # Create a symbolic link to the original directory.
            # This avoids copying potentially large directories.
            if src.is_dir():
                (stage / name).symlink_to(src.resolve(), target_is_directory=True)

            # Copy the file into the temporary working directory.
            # copy2 also preserves metadata such as modification timestamps.
            elif src.is_file():
                shutil.copy2(src, stage / name)

            # If neither a file nor a directory exists, the input is missing.
            else:
                pytest.fail(f'Required Fortran input missing: {src}')

        # Create the directory for the screening coefficients.
        (stage / 'screencoefs').mkdir()

        # Add a format header to the coefficient file.
        # The numerical coefficient values themselves remain unchanged.
        (stage / 'screencoefs/NLHLIN_python_coeffs.dat').write_text('# nexp=3\n' + COEFFICIENTS.read_text())

        # Compile and execute each requested Fortran driver.
        for driver in driver_paths:
            executable = stage / driver.stem

            # Define the path where the compiler output will be saved.
            build_log = directory / f'build_{driver.stem}.log'

            # Define the path where the program execution output will be saved.
            run_log = directory / f'run_{driver.stem}.log'

            # Define the path where the program execution output will be saved.
            run_log = directory / f'run_{driver.stem}.log'

            # Create the compiler command with the required options.
            # -g: Include debugging information.
            # -O0: Disable compiler optimizations.
            # -fcheck=all: Enable runtime checks.
            # -fbacktrace: Show a traceback if a runtime error occurs.
            # -Wall: Enable compiler warnings.
            # -J: Specify the directory for generated module files.
            # -cpp: Enable the C preprocessor.
            # -o: Specify the output executable.
            cmd = ['gfortran', '-g', '-O0', '-fcheck=all', '-fbacktrace', '-Wall',
                   '-J', str(stage / 'mod'), '-cpp', '-o', str(executable)]

            # Try to compile the Fortran sources and execute the resulting program.
            try:
                with build_log.open('wb') as log:
                    subprocess.run(cmd + [str(p) for p in sources + common_paths + [driver]],
                                   check=True, stdout=log, stderr=subprocess.STDOUT, timeout=600)
                with run_log.open('wb') as log:
                    for args in (runs or [()]):
                        subprocess.run([str(executable), *map(str, args)], cwd=stage, check=True, stdout=log,
                                       stderr=subprocess.STDOUT, timeout=600)

            # Abort the test and report the error.
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                pytest.fail(f'Fortran {group} build/run failed, not falling back to stale data: '
                            f'{exc}; logs: {build_log}, {run_log}')

        # Collect all TXT result files generated in the temporary directory.
        outputs = list(stage.glob('*.txt'))
        
        # Check whether the Fortran drivers generated any result files.
        if not outputs:
            pytest.fail(f'Fortran {group} generated no result files')

        # Publish the generated results only after all requested runs succeeded.
        # Each file is first copied to a temporary target and then renamed.
        # This prevents partially written individual result files.
        for result in outputs:
            target = directory / result.name
            tmp_target = target.with_suffix('.tmp')
            shutil.copy2(result, tmp_target)
            tmp_target.replace(target)

        # Create a dictionary containing metadata about the generated reference data.
        # This includes the schema version, coefficient hash, run arguments, source file digests, and output file digests.
        metadata = {'schema': 1, 'calculation_api': 'production-chi-legendre-v1',
                    'coefficient_sha256': COEFFICIENT_HASH,
                    'runs': [list(args) for args in (runs or [()])],
                    'sources': {str(p): digest(p) for p in sources + common_paths + driver_paths},
                    'outputs': {p.name: digest(p) for p in outputs}}

        # Define the path to the manifest JSON file.
        # This file stores the metadata and hashes of the generated results.
        manifest = directory / f'{group}_manifest.json'

        # Create a temporary filename for writing the manifest.
        temporary = manifest.with_suffix('.tmp')

        # Convert the metadata dictionary into formatted JSON and save it.
        # indent=2 makes the JSON file easier to read.
        temporary.write_text(json.dumps(metadata, indent=2))

        # Replace the old manifest with the new one.
        # This ensures that the manifest is only published after it is complete.
        temporary.replace(manifest)
    print(f'Fortran {group}: generated {len(outputs)} files with shared NLHlin_3eV coefficients')
    return directory


@pytest.fixture(scope='session')
def build_and_run_f90():
    """Prepare ZBL and NLHlin references once per session.

    Args:
        None
    Returns:
        The directory containing the cached Fortran reference data.

    The driver covers two models, three configurations and 92 pairs each.

    The reference data is cached in tests/data/f90/cm, and validated by a manifest.


    """
    base = TESTS / 'fortran'

    # executes the Fortran drivers to generate reference data for the 'cm' group, using the specified source files and comparison data.
    return generate_fortran('cm', [base / 'generate_zbl_scattering_references.f90', base / 'generate_nlhlin_scattering_references.f90'],
                            [base / 'scattering_reference_writers.f90'])


@pytest.fixture(scope='session')
def screening_fortran():
    """Prepare screening references for screening tests once per session.

    Args:
        None
    Returns:
        The directory containing the cached Fortran reference data.
    """

    # executes the Fortran driver to generate reference data for the 'screening' group, using the specified source file.
    return generate_fortran('screening', [TESTS / 'fortran/generate_screening_references.f90'])


@pytest.fixture(scope='session')
def extended_fortran():
    """Prepare different pair and boundary references for extended tests once per session.

    Args:
        None
    Returns:
        The directory containing the cached Fortran reference data.
    """
    from scattering_support import EXTENDED_PAIRS

    # Define the runs for the extended tests, which include different models and pairs of atomic numbers (z1, z2) with four Legendre quadrature points.
    runs = [(model, z1, z2, 4)
            for model in ('ZBL', 'NLHlin')
            for z1, z2 in EXTENDED_PAIRS]

    # executes the Fortran driver to generate reference data for the 'extended' group, using the specified source file and the defined runs.
    return generate_fortran('extended', [TESTS / 'fortran/generate_extended_scattering_references.f90'], runs=runs)
