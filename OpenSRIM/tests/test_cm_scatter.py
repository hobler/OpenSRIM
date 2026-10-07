"""Black-box scattering comparisons; shared setup lives in scattering_support."""
import numpy as np
import pytest
from scattering_support import (
    get_python_screen_fun,
    build_numba_potential,
    all_models,
    chi_u_values,
    chi_near_zero_u_values,
    chi_near_one_values,
    calculate_scatter,
    calculate_scatter_python,
    calculate_scatter_numba,
    cm_scatter_nb,
    assert_arrays_equal,
    load_fortran_data,
    record_max_error,
    compare_three_way,
    finish_comparison,
    compare_chi_reference,
    E_VALUES_ERROR,
    P_VALUES_ERROR,
    E_VALUES_TAU_THETA,
    P_VALUES_TAU_THETA,
    E_VALUES_SIN,
    P_VALUES_SIN,
    ATOL_TAU_FORTRAN,
    TOLERANCES,
    EXTENDED_PAIRS,
    MEHLER_ORDERS,
    pair_inputs,
    extended_grid,
    extended_outputs,
)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_chi_python_numba_fortran(model):
    """Compare the screening-derived chi diagnostic on the original plot grid.

    All three backends supply native screening values. This does not test
    three independent internal chi implementations of the scattering solvers."""
    # The shared helper reports all three pairwise differences.
    u = chi_u_values()
    compare_chi_reference(model, "", u, u)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_chi_near_zero_python_numba_fortran(model):
    """Compare screening-derived chi near u=0 across all three backends.

    The logarithmic grid covers the transition to the Taylor approximation
    used to avoid cancellation of nearly equal screening values."""
    # The grid helper excludes the singular endpoint u=1.
    u = chi_near_zero_u_values()
    compare_chi_reference(model, "_near_zero", u, u)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_chi_near_one_python_numba_fortran(model):
    """Compare the screening-derived chi diagnostic near u=1.

    Small values of 1-u correspond to large evaluation radii and therefore
    also probe the screening tail."""
    # Reference coordinates are 1-u, while the calculation receives u.
    u1, u = chi_near_one_values()
    compare_chi_reference(model, "_near_one", u1, u)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_head_on_scattering_angle_python_numba_fortran(model):
    """Check head-on scattering at p=0 for both potential models.

    Besides pairwise agreement, each backend must return the physical
    backscattering angle pi. This test does not check tau."""
    # Pairwise agreement alone would not detect a common error in all backends.
    screen = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)
    summaries, failures = {}, []
    for e in E_VALUES_SIN:
        (theta_py, _), (theta_nb, _) = calculate_scatter(e, [0.0], screen, pot_model, pot_coefs)
        data = load_fortran_data(f"head_on_{model.lower()}_e{e:.0e}.txt", [0.0])
        values = dict(Python=theta_py, Numba=theta_nb, Fortran=data[:, 1])
        compare_three_way(summaries, failures, values, [0.0], f"e={e}", "boundary", tolerances=TOLERANCES)
        for label, theta in values.items():
            try:
                assert_arrays_equal([np.pi], theta, rtol=0, atol=1e-14)
            except AssertionError as exc:
                failures.append(f"{label} head-on angle: {exc}")
    finish_comparison(f"{model} head-on theta", summaries, failures)


@pytest.mark.usefixtures("build_and_run_f90")
def test_scattering_at_and_beyond_cutoff_python_numba_fortran():
    """Check vanishing NLHlin scattering at and beyond the cutoff.

    All backends must return theta=tau=0. Also check that the direct Numba
    apsis call returns p without performing a Newton iteration."""
    # nextafter selects the immediately adjacent representable value above rmax.
    screen = get_python_screen_fun("NLHlin")
    pot_model, pot_coefs = build_numba_potential("NLHlin")
    assert_arrays_equal(screen.rmax, pot_coefs.rmax, rtol=0, atol=1e-14)
    p_values = np.array([screen.rmax, np.nextafter(screen.rmax, np.inf), 30.0])
    summaries, failures = {}, []
    for e in E_VALUES_SIN:
        for p in p_values:
            r0, count = cm_scatter_nb.get_apsis(e, p, pot_model, pot_coefs)
            assert r0 == p and count == 0
        python, numba = calculate_scatter(e, p_values, screen, pot_model, pot_coefs)
        data = load_fortran_data(f"cutoff_nlhlin_e{e:.0e}.txt", p_values, columns=3)
        for index, quantity in enumerate(("theta", "tau")):
            values = dict(Python=python[index], Numba=numba[index], Fortran=data[:, index+1])
            compare_three_way(summaries, failures, values, p_values, f"e={e}", "boundary", tolerances=TOLERANCES,
                              prefix=f"{quantity}: ")
            for label, result in values.items():
                try:
                    assert_arrays_equal(np.zeros_like(p_values), result, rtol=0, atol=0)
                except AssertionError as exc:
                    failures.append(f"{label} {quantity} at cutoff: {exc}")
    finish_comparison("NLHlin cutoff", summaries, failures)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_scatter_integrals_python_numba_fortran(model):
    """Compare theta and tau on three energy-dependent impact-parameter grids.

    Python and Numba use four Legendre points. Fortran supplies the
    four-point Legendre angle and tau from 512 Mehler points."""
    # Angle and time integral have separate, explicit error budgets.
    screen = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)
    summaries, failures = {}, []
    for e, p_values in zip(E_VALUES_ERROR, P_VALUES_ERROR):
        python, numba = calculate_scatter(e, p_values, screen, pot_model, pot_coefs)
        data = load_fortran_data(f"scatter_native_{model.lower()}_e{e:.0e}.txt", p_values, columns=3)
        for index, quantity in enumerate(("theta", "tau")):
            compare_three_way(summaries, failures,
                              dict(Python=python[index], Numba=numba[index], Fortran=data[:, index+1]),
                              p_values, f"e={e}", quantity, tolerances=TOLERANCES, prefix=f"{quantity}: ")
    finish_comparison(f"{model} scattering", summaries, failures)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_tau_error_data_python_numba_fortran(model):
    """Compare native time integrals and report a finer Python calculation.

    Three-way comparisons are assertions. The additional comparison with
    Python at 32 points only diagnoses quadrature errors."""
    # The separate accuracy test supplies the actual convergence assertions.
    screen = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)
    summaries, failures = {}, []
    for e, p_values in zip(E_VALUES_ERROR, P_VALUES_ERROR):
        _, tau_ref = calculate_scatter_python(e, p_values, screen, n_absc=32)
        (_, tau_py), (_, tau_nb) = calculate_scatter(e, p_values, screen, pot_model, pot_coefs)
        data = load_fortran_data(f"scatter_native_{model.lower()}_e{e:.0e}.txt", p_values, columns=3)
        assert np.all(np.isfinite(tau_ref))
        values = dict(Python=tau_py, Numba=tau_nb, Fortran=data[:, 2])
        compare_three_way(summaries, failures, values, p_values, f"e={e}", "tau", tolerances=TOLERANCES)
        for label, result in values.items():
            record_max_error(summaries, f"{label} vs Python 32-point (diagnostic)",
                             result, tau_ref, p_values, context=f"e={e}",
                             labels=(label, "Python 32-point"))
    finish_comparison(f"{model} tau and convergence", summaries, failures)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_tau_over_theta_data_python_numba_fortran(model):
    """Compare tau as a function of theta on the original plotting grid.

    Both coordinates are checked separately across all three backends.
    Python and Numba use the four-point default configuration here."""
    # The plot uses theta as its axis; it does not divide tau by theta.
    screen = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)
    summaries, failures = {}, []
    p_values = P_VALUES_TAU_THETA
    for e in E_VALUES_TAU_THETA:
        python, numba = calculate_scatter(e, p_values, screen, pot_model, pot_coefs)
        data = load_fortran_data(f"scatter_ratio_{model.lower()}_e{e:.0e}.txt", p_values, columns=3)
        outputs = dict(Python=python, Numba=numba, Fortran=data[:, 1:].T)
        for index, quantity in enumerate(("theta", "tau")):
            compare_three_way(summaries, failures,
                              {label: result[index] for label, result in outputs.items()},
                              p_values, f"e={e}", quantity, tolerances=TOLERANCES, prefix=f"{quantity}: ")
    finish_comparison(f"{model} tau as a function of theta", summaries, failures)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.usefixtures("build_and_run_f90")
def test_sinhalftheta_over_p_data_python_numba_fortran(model):
    """Compare sin²(theta/2) as a function of impact parameter p.

    Each native scattering angle is converted to the quantity displayed
    by the original plotting function."""
    # This conversion defines the tested quantity; it does not correct results.
    screen = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)
    summaries, failures = {}, []
    for e in E_VALUES_SIN:
        (theta_py, _), (theta_nb, _) = calculate_scatter(e, P_VALUES_SIN, screen, pot_model, pot_coefs)
        data = load_fortran_data(f"theta_sin_{model.lower()}_e{e:.0e}.txt", P_VALUES_SIN)
        values = {label: np.sin(theta/2)**2 for label, theta in
                  (("Python", theta_py), ("Numba", theta_nb), ("Fortran", data[:, 1]))}
        compare_three_way(summaries, failures, values, P_VALUES_SIN, f"e={e}", "sin2", tolerances=TOLERANCES)
    finish_comparison(f"{model} sin²(theta/2)", summaries, failures)


@pytest.mark.parametrize("model", all_models())
@pytest.mark.parametrize("energy_index", range(len(E_VALUES_ERROR)))
@pytest.mark.usefixtures("build_and_run_f90")
def test_theta_error_data_python_numba_fortran(model, energy_index):
    """Compare four-point scattering angles at one parametrized energy.

    Each potential uses identical impact parameters for Python, Numba
    and the native Fortran LegendreChiCorr method."""
    # The energy selects the corresponding impact-parameter and reference grid.
    e = E_VALUES_ERROR[energy_index]
    p_values = P_VALUES_ERROR[energy_index]
    data = load_fortran_data(f"theta_error_{model.lower()}_e{e:.0e}_n4.txt", p_values)
    screen = get_python_screen_fun(model)
    pot_model, pot_coefs = build_numba_potential(model)
    (theta_py, _), (theta_nb, _) = calculate_scatter(e, p_values, screen, pot_model, pot_coefs)
    summaries, failures = {}, []
    compare_three_way(summaries, failures, dict(Python=theta_py, Numba=theta_nb, Fortran=data[:, 1]),
                      p_values, f"e={e}", "theta", tolerances=TOLERANCES)
    finish_comparison(f"{model} four-point theta", summaries, failures)


@pytest.mark.usefixtures("build_and_run_f90")
def test_scatter_tau_python_numba_fortran():
    """Compare ZBL tau with the Fortran table routines cmangl and timeint.

    Each Fortran routine is compared separately with Python and Numba.
    Table interpolation has its own absolute error budget."""
    # These legacy tables are evaluated independently of Mehler quadrature.
    screen = get_python_screen_fun("ZBL")
    pot_model, pot_coefs = build_numba_potential("ZBL")
    summaries, failures = {}, []
    for e, p_values in zip(E_VALUES_ERROR, P_VALUES_ERROR):
        data = load_fortran_data(f"tau_zbl_e{e:.0e}.txt", p_values)
        timeint = load_fortran_data(f"tau_timeint_zbl_e{e:.0e}.txt", p_values)
        assert_arrays_equal(data[:, 0], timeint[:, 0], rtol=0, atol=0)
        (_, tau_py), (_, tau_nb) = calculate_scatter(e, data[:, 0], screen, pot_model, pot_coefs)
        for api, tau_f90 in (("cmangl", data[:, 1]), ("timeint", timeint[:, 1])):
            compare_three_way(summaries, failures, dict(Python=tau_py, Numba=tau_nb, Fortran=tau_f90),
                              data[:, 0], f"e={e}", "tau", tolerances=TOLERANCES, prefix=f"{api}: ",
                              fortran_tolerance=(1e-3, ATOL_TAU_FORTRAN))
    finish_comparison("ZBL legacy tau", summaries, failures)


@pytest.mark.parametrize('model', all_models())
@pytest.mark.parametrize('pair', EXTENDED_PAIRS, ids=lambda pair: f'Z{pair[0]}-Z{pair[1]}')
def test_additional_element_pairs_python_numba_fortran(model, pair, extended_fortran):
    """Compare theta and tau for H–H, Au–Au, H–Au and B–Si across all backends.

    Three energies and seven impact parameters cover light, heavy and
    strongly asymmetric element pairs for both potential models."""
    # Each backend reads its own coefficients from the shared input dataset.
    screen, coefs = pair_inputs(model, *pair)
    energies, points = extended_grid('pairs', model, screen.rmax)
    fortran = extended_outputs(extended_fortran, 'pairs', model, pair, 512, energies, points)
    summaries, failures = {}, []
    for i, e in enumerate(energies):
        python, numba = calculate_scatter(e, points, screen, model, coefs)
        for j, quantity in enumerate(('theta', 'tau')):
            compare_three_way(summaries, failures,
                              dict(Python=python[j], Numba=numba[j], Fortran=fortran[i, j]),
                              points, f'e={e}, pair={pair}', quantity,
                              tolerances=TOLERANCES, prefix=f'{quantity}: ')
    finish_comparison(f'{model} pair={pair}', summaries, failures)


@pytest.mark.parametrize('model', all_models())
def test_small_impact_and_inside_cutoff_python_numba_fortran(model, extended_fortran):
    """Check small positive impact parameters and the approach to the NLHlin cutoff.

    Besides three-way comparisons, check the Numba apsis bounds and
    the physical angle range of every backend."""
    # Cutoff inputs must be bitwise identical; extended_outputs enforces this.
    screen, coefs = pair_inputs(model, 5, 14)
    energies, points = extended_grid('boundary', model, screen.rmax)
    fortran = extended_outputs(extended_fortran, 'boundary', model, (5, 14), 512, energies, points)
    summaries, failures = {}, []
    for i, e in enumerate(energies):
        for p in points:
            r0, iterations = cm_scatter_nb.get_apsis(e, p, model, coefs)
            assert np.isfinite(r0) and p <= r0 <= screen.rmax
            assert 0 <= iterations <= 100
        python, numba = calculate_scatter(e, points, screen, model, coefs)
        for j, quantity in enumerate(('theta', 'tau')):
            values = dict(Python=python[j], Numba=numba[j], Fortran=fortran[i, j])
            compare_three_way(summaries, failures, values, points, f'e={e}', quantity,
                              tolerances=TOLERANCES, prefix=f'{quantity}: ')
        for label, angles in (('Python', python[0]), ('Numba', numba[0]), ('Fortran', fortran[i, 0])):
            if not np.all((angles >= -1e-11) & (angles <= np.pi+1e-11)):
                failures.append(f'{label}, e={e}: angle outside [0, pi] or nonfinite')
    finish_comparison(f'{model} small p and approach to cutoff', summaries, failures)


# Accuracy target is intentionally separate from the loose legacy tau comparison.
# It requires approximately 0.1% accuracy, with a 1e-6 absolute floor near zero.
TAU_ACCURACY = {'tau': {'numba': (1e-3, 1e-6), 'fortran': (1e-3, 1e-6)}}


@pytest.mark.parametrize('model', all_models())
def test_tau_accuracy_and_convergence_python_numba_fortran(model, extended_fortran):
    """Check tau accuracy and quadrature convergence across all three backends.

    Python and Numba use 128 points for the accuracy comparison; Fortran
    uses its highest tested order. Adjacent high orders must agree, and
    the Numba default must remain equivalent to an explicit four-point call."""
    # A high order alone does not prove convergence; convergence is asserted separately.
    screen, coefs = pair_inputs(model, 5, 14)
    with pytest.raises(ValueError, match="Unsupported Legendre order"):
        cm_scatter_nb.scatter_integrals(1., 1., model, coefs, 7)
    energies, points = extended_grid('convergence', model, screen.rmax)
    native = {n: extended_outputs(extended_fortran, 'convergence', model, (5, 14), n, energies, points)
              for n in MEHLER_ORDERS}
    previous_order, final_order = MEHLER_ORDERS[-2:]
    summaries, failures = {}, []
    for i, e in enumerate(energies):
        _, numba_default = calculate_scatter(e, points, screen, model, coefs)
        numba = {n: calculate_scatter_numba(e, points, model, coefs, n_absc=n)
                 for n in (4, 16, 32, 64, 128, 256)}
        assert_arrays_equal(numba_default, numba[4], rtol=0, atol=0)
        python = {n: calculate_scatter_python(e, points, screen, n_absc=n)[1] for n in (4, 32, 64, 128)}
        # Assert that the final two orders agree, rather than calling a single
        # high-order result an exact reference.
        for label, low, high in (
                ('Python 64 vs 128', python[64], python[128]),
                ('Numba 64 vs 128', numba[64][1], numba[128][1]),
                ('Numba 128 vs 256', numba[128][1], numba[256][1]),
                (f'Fortran {previous_order} vs {final_order}',
                 native[previous_order][i, 1], native[final_order][i, 1])):
            record_max_error(summaries, label, low, high, points, context=f'e={e}',
                             labels=('lower order', 'higher order'))
            try:
                assert_arrays_equal(high, low, rtol=1e-3, atol=1e-6)
            except AssertionError as exc:
                failures.append(f'{label}, e={e}: convergence target not reached: {exc}')
        for n in (4, 32):
            record_max_error(summaries, f'Python {n} vs 128', python[n], python[128], points,
                             context=f'e={e}', labels=(f'Python {n}', 'Python 128'))
        for n in (4, 16, 32):
            record_max_error(summaries, f'Numba {n} vs 128', numba[n][1], numba[128][1], points,
                             context=f'e={e}', labels=(f'Numba {n}', 'Numba 128'))
        for n in MEHLER_ORDERS[:-1]:
            record_max_error(summaries, f'Fortran {n} vs {final_order}', native[n][i, 1], native[final_order][i, 1], points,
                             context=f'e={e}', labels=(f'Fortran {n}', f'Fortran {final_order}'))
        compare_three_way(summaries, failures,
                          dict(Python=python[128], Numba=numba[128][1], Fortran=native[final_order][i, 1]),
                          points, f'e={e}', 'tau', tolerances=TAU_ACCURACY,
                          prefix='tau accuracy: ')
    finish_comparison(f'{model} tau accuracy and quadrature convergence', summaries, failures,
                      note=f'Accuracy comparison: Python 128, Numba 128 (default remains 4), Fortran {final_order} points; convergence checks are mandatory.')
