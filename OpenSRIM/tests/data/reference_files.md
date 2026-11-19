# Test reference files

Python and Fortran use the same readable naming scheme. All data paths in this
document are relative to this directory (`OpenSRIM/tests/data/`).
The linked generators are located one level above, under `tests/`.

| Directory | Contents |
| --- | --- |
| `py/reference_v1/` | Python results and `python_manifest.json` |
| `f90/` | Standard Chi/scattering references and `cm_manifest.json` |
| `f90/screening/` | Screening references and `screening_manifest.json` |
| `f90/extended/` | Extended scattering references and `extended_manifest.json` |

The Fortran group `cm` is stored directly in `f90/`, not in a `cm/` subdirectory.

```text
<kind>_<model>_z<z1>_z<z2>_norm_<norm>[_e<energy>][_r0_<radius>][_n<order>]_grid_<grid>.<extension>
```

| Component | Meaning |
| --- | --- |
| `kind` | `screen`, `scatter`, `chi` or `theta` (Fortran files containing only the angle) |
| `model` | `zbl` or `nlhlin` |
| `z1`, `z2` | Atomic numbers, in input order |
| `norm` | `native` for model-specific normalization; `c1`, `c2`, `c3` for the screening configuration in `test_extended_element_pairs.py` |
| `energy` | Reduced energy, e.g. `1e-04` |
| `radius` | Actual reference radius `r0`, not a list index |
| `order` | Legendre quadrature order, e.g. `4` or `32` |
| `grid` | Explicit label of the test grid used |
| `extension` | `npz` for Python, `txt` for Fortran |

Examples:

```text
screen_nlhlin_z5_z5_norm_c1_grid_screening.npz
screen_nlhlin_z5_z5_norm_c1_grid_screening.txt
scatter_nlhlin_z5_z14_norm_native_e1e-04_n4_grid_error.npz
scatter_nlhlin_z5_z14_norm_native_e1e-04_n4_grid_error.txt
chi_zbl_z5_z14_norm_native_r0_1e-02_grid_near_one.npz
chi_zbl_z5_z14_norm_native_r0_1e-02_grid_near_one.txt
```

The screening configuration identifies the entire comparison, including unit
conversions. For ZBL, the Python call internally uses native normalization;
`norm_c1` and similar labels identify the comparison configuration, not a
different internal normalization value.

## Test grids

| Grid name | Inputs |
| --- | --- |
| `screening` | 101 radii per element pair and screening configuration |
| `regular` | 999 evenly spaced Chi coordinates |
| `near_zero`, `near_one` | Chi grids near the integration boundaries |
| `head_on` | Zero impact parameter |
| `cutoff` | Impact parameters at or beyond the NLHlin range |
| `error` | 101 impact parameters with an energy-dependent upper bound |
| `ratio` | 201 impact parameters for Tau-versus-Theta comparisons |
| `sin` | 201 impact parameters for sine-squared comparisons |
| `pairs` | Seven impact parameters for additional element pairs |
| `boundary` | Seven impact parameters near zero or inside the potential range |

Each extended Fortran file contains only **one energy**, just like its Python
counterpart. The previous combined files are no longer needed. The TXT columns
remain `e, p, theta, tau`. Files with matching names across backends identify
the same test case; their arrays or columns may differ (e.g. Python:
`phi, chi`, Fortran: `coordinate, chi`).

## Filenames and generators

In the following table, `<basis>` stands for
`<model>_z<z1>_z<z2>_norm_<norm>`.

| Filename pattern | Location | Generator |
| --- | --- | --- |
| `screen_<basis>_grid_screening.npz` | `py/reference_v1/` | [reference_data.py](../reference_data.py): `PythonScreen.__call__()` → `python_result()` |
| `scatter_<basis>_e<energy>_n<order>_grid_<grid>.npz` | `py/reference_v1/` | [reference_data.py](../reference_data.py): `python_scatter()` → `python_result()` |
| `chi_<basis>_r0_<radius>_grid_<grid>.npz` | `py/reference_v1/` | [reference_data.py](../reference_data.py): `python_chi()` → `python_result()` |
| `screen_<basis>_grid_screening.txt` | `f90/screening/` | [generate_screening_references.f90](../fortran/generate_screening_references.f90): program `generate_screening_references` |
| `chi_<basis>_r0_<radius>_grid_<grid>.txt` | `f90/` | [scattering_reference_writers.f90](../fortran/scattering_reference_writers.f90): `write_chi_data()` |
| `scatter_<basis>_e<energy>_n<order>_grid_<grid>.txt` | `f90/` | [scattering_reference_writers.f90](../fortran/scattering_reference_writers.f90): `write_comparison_data()`; grids `cutoff`, `ratio`, `error` |
| `theta_<basis>_e<energy>_n4_grid_<grid>.txt` | `f90/` | [scattering_reference_writers.f90](../fortran/scattering_reference_writers.f90): `write_comparison_data()`; grids `head_on`, `sin`, `error` |
| `scatter_<basis>_e<energy>_n<order>_grid_<grid>.txt` | `f90/extended/` | [generate_extended_scattering_references.f90](../fortran/generate_extended_scattering_references.f90): program `generate_extended_scattering_references`; grids `pairs`, `boundary` |
| `python_manifest.json` | `py/reference_v1/` | [reference_data.py](../reference_data.py): `python_result()` |
| `<group>_manifest.json` | Fortran group directory | [conftest.py](../conftest.py): `generate_fortran()` |
| `build_<driver>.log`, `run_<driver>.log` | Fortran group directory | [conftest.py](../conftest.py): `generate_fortran()`; diagnostic logs, not reference data |

`reference_filename()` in [reference_data.py](../reference_data.py) constructs
Python filenames and the Fortran filenames expected by the tests. The Fortran
drivers construct their own output names using the same scheme; they do not
call this Python function.
`generate_fortran()` compiles and runs the drivers in a temporary directory
and publishes their TXT files only after successful execution.

Python also stores both scattering values under `scatter_…` for the `head_on`,
`sin` and `error` grids. Fortran angle files are named `theta_…` because they
contain only the angle. Consistent naming means a shared scheme, not an
identical filename for every comparison.

## Contents and input validation

| File type | Numerical contents |
| --- | --- |
| Python `screen_…npz` | Array containing Phi and its radial derivative |
| Python `scatter_…npz` | Array containing Theta and Tau |
| Python `chi_…npz` | Array containing Phi and Chi |
| Fortran `screen_…txt` | Columns: input radius, Phi, derivative with respect to the input radius |
| Fortran `chi_…txt` | Columns: coordinate, Chi; for `near_one`, the coordinate is `1-u` |
| Fortran `theta_…txt` | Columns: impact parameter, Theta |
| Fortran `scatter_…txt` in `f90/` | Columns: impact parameter, Theta, Tau |
| Fortran `scatter_…txt` in `f90/extended/` | Columns: energy, impact parameter, Theta, Tau |

Filenames contain no hash. Python still stores the full parameters, complete
input grid, normalization/grid labels and coefficient hash in the metadata.
Source hashes record provenance. Metadata must match exactly when loading.
Even during live calculations, a name must not silently be reused for different
parameters or grids: a collision causes a failure before overwriting the file.
Results may be recomputed with changed coefficients when original code is available.

The Python naming function formats numbers with up to twelve decimal places in
the scientific mantissa, omitting trailing zeros. If two inputs produce the
same name, exact metadata validation prevents them from being confused.
The Fortran drivers use corresponding predefined numeric labels for their fixed
energy and radius lists. New grids or normalization configurations require
new labels, and both filename implementations must be kept in sync.

## Numerical tolerances

`assert_arrays_equal(reference, result, rtol=..., atol=...)` checks that values
are finite and, for each element:

```text
abs(result - reference) <= atol + rtol * abs(reference)
```

`rtol` is relative and dimensionless; `atol` has the units of the quantity being
compared and is particularly important near zero. Smaller values tighten the check.
`TOLERANCES` in `scattering_support.py` stores `(rtol, atol)` by quantity and
backend combination. Screening tests use `SCREEN_TOLERANCES`.
`TOLERANCES["default"] = (1e-3, 1e-6)` holds the defaults; all current callers
explicitly supply both tolerances.

## Integrity checks and manifests

Under `outputs`, [`py/reference_v1/python_manifest.json`](py/reference_v1/python_manifest.json)
stores the SHA-256 of each complete NPZ reference file. `python_result()` only
adds freshly computed results and replaces the manifest through a temporary
file. In cached mode, the manifest must exist and be valid; files without an
entry or with a mismatched hash are rejected before loading. Exact input
metadata validation follows.

A missing manifest is rebuilt only through live calculations, not from
unverified old files. To regenerate it completely, all Python comparison
cases must run with the original sources available.

Fortran's `cached_fortran()` checks the coefficient hash and the SHA-256 values
of all result files listed in the group manifest. For `cm` and `extended`,
`calculation_api` must also equal `production-chi-legendre-v1`.
The tests subsequently check shapes and input grids using their specified
tolerances. The code still permits legacy screening data without a manifest
if TXT files exist; the current screening references have a manifest.

Python metadata validation checks the complete grid exactly. Fortran does not
store an identical per-file input metadata structure; its test grid checks
complement manifest validation. For both backends, a file hash confirms that
the stored bytes are unchanged, not that the results are physically correct
or up to date with modified sources.

## Source code changes

After changes to the calculation methods, Python and Fortran references must
be recomputed in live mode. A stored file hash confirms that a reference is
unchanged, not that it reflects the current source code. Source hashes remain
provenance information for both backends and are not checked against current
sources in cached mode.
