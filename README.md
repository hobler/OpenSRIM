# CyTRIM
A Python/Cython implementation of TRIM.

## Directories:
    pytrim: pure Python code, kept for reference.
    cytrim: Cython code ("model")
    gui: (Most of) the GUI code ("view")
    doc: Documentation

## KORAL toml file structure:
    [params]
    method = 'ZBL'
    z_ion = 33
    m_ion = 74.992
    z_target = [14]
    m_target = [28.086]
    d_target = [0.04996]
    c_target = [1]
    s_e_corr = 1.0
    start_energy = 1
    stop_energy = 10e6
    nr_values = 100

    [settings]
    nr_iterations = 1
    integration_method = 'LSODA'
    rtol = 1e-6
    atol = 1e-6
