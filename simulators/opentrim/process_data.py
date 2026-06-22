import sys
from pathlib import Path
import numpy as np
from .stats import standardize_moments

def _get_element_names_from_input(input_params):
    beam = input_params["beam"]
    layers = input_params["layer"]

    elements = [(beam["symbol"], beam["name"], beam["Z"], beam["M"])]
    for layer in layers:
        for elem in layer["element"]:
            descriptor = (elem["symbol"], elem["name"], elem["Z"], elem["M"])
            if descriptor not in elements[1:]:
                elements.append(descriptor)

    elems = {}
    for _, name, z, _ in elements:
        elems[name] = int(z)
    nelem_target = len(elements) - 1
    return elems, nelem_target

def _build_distribution_headers(key, elems, follow_recoils, nelem_target):
    elem_names = list(elems.keys())

    if key in ["x", "y"]:
        mom_species = elem_names + elem_names[1:] if follow_recoils else elem_names[:]
        species_labels = mom_species[:]
        header_indexes = [0]
        for i, label in enumerate(species_labels):
            z = elems[label]
            if i == 0:
                species_labels[i] = f"{label} (ion)"
                header_indexes.append(z)
                continue
            if (i - 1) // nelem_target:
                species_labels[i] = f"{label} (vacancies)"
                header_indexes.append(-z)
            else:
                species_labels[i] = f"{label} (interstitials)"
                header_indexes.append(z)
        return header_indexes, [f"{key.upper()} Position"] + species_labels

    if key[0] in ["b", "t"]:
        kind_map = {"b": "Backscattered", "t": "Transmitted", "e": "energy", "a": "angle"}
        mom_species = elem_names[:] if follow_recoils else elem_names[:1]
        species_labels = mom_species[:]
        header_indexes = [0]
        for i, label in enumerate(species_labels):
            z = elems[label]
            if i == 0:
                species_labels[i] = f"{label} (ion)"
            header_indexes.append(z)
        return header_indexes, [f"{kind_map[key[0]]} {kind_map[key[1]]}"] + species_labels

    if key[1] in ["n", "e"]:
        kind_map = {"n": "NED", "e": "EED"}
        mom_species = elem_names[:]
        species_labels = mom_species[:]
        header_indexes = [0]
        for i, label in enumerate(species_labels):
            z = elems[label]
            if i == 0:
                species_labels[i] = f"{label} (ion)"
            header_indexes.append(z)
        return header_indexes, [f"{key[0].upper()} {kind_map[key[1]]}"] + species_labels

    return None, None

def _header_and_column_row(indexes, labels):
    column_row = ", ".join([str(el) for el in indexes])
    header = "\n".join([f"{k}:{v}" for k, v in zip(indexes, labels)])
    return header, column_row

def _write_histogram(path, val, x_vals, header_indexes, header_elems):
    header, column_row = _header_and_column_row(header_indexes, header_elems)
    data = np.vstack((x_vals, val["counts"][:, 1:-1])).T
    with open(path, "w") as f:
        np.savetxt(f, data, delimiter=", ", fmt="%d", header=header + "\n" + column_row)

def _write_histogram_binary_2d(path, val, species_labels):
    counts = np.asarray(val["counts"][:, 1:-1, 1:-1], dtype="<f8")
    n_species, _, _ = counts.shape
    nx = int(val["x_nbins"])
    ny = int(val["y_nbins"])
    x_values = np.linspace(val["x_limits"][0], val["x_limits"][1], nx, dtype="<f8")
    y_values = np.linspace(val["y_limits"][0], val["y_limits"][1], ny, dtype="<f8")
    labels = np.asarray(species_labels, dtype="S32")
    header = [
        np.array([0x00fa], dtype="<u2"),
        np.array([n_species, nx, ny], dtype="<u4"),
        x_values,
        y_values,
        labels,
    ]
    with open(path, "wb") as f:
        for item in header:
            item.tofile(f)
        counts.tofile(f)


def _write_moments(path, val, species_labels, species_indexes):
    metric_names = ["total", "mean", "std", "skewness", "kurtosis"]
    max_order_local = (val["power_sums"].shape[1] - 1) // 2
    n_metrics = min(max_order_local + 1, len(metric_names))

    header_indexes = [0]
    header_elems = ["Moment order"]
    data_cols = [np.arange(n_metrics, dtype=np.float64)]
    for ivar, label in enumerate(species_labels):
        std_moments, std_err_moments = standardize_moments(val, ivar)
        val_idx = species_indexes[ivar]

        header_indexes.append(val_idx)
        header_elems.append(f"{label} moments")
        data_cols.append(std_moments[:n_metrics])

        header_indexes.append(val_idx)
        header_elems.append(f"{label} moments errors")
        data_cols.append(std_err_moments[:n_metrics])

    header, column_row = _header_and_column_row(header_indexes, header_elems)
    data = np.column_stack(data_cols)
    with open(path, "w") as f:
        np.savetxt(f, data, delimiter=", ", fmt="%.8e", header=header + "\n" + column_row)

def _write_raw_moments(path, val, species_labels, species_indexes):
    raw_nmom = val["power_sums"].shape[1]
    header_indexes = [0]
    header_elems = ["Moment order"]
    data_cols = [np.arange(raw_nmom, dtype=np.float64)]
    for ivar, label in enumerate(species_labels):
        header_indexes.append(species_indexes[ivar])
        header_elems.append(f"{label} raw moments")
        data_cols.append(val["power_sums"][ivar, :])

    header, column_row = _header_and_column_row(header_indexes, header_elems)
    data = np.column_stack(data_cols)
    with open(path, "w") as f:
        np.savetxt(f, data, delimiter=", ", fmt="%.8e", header=header + "\n" + column_row)

def write_stats(input_params, stats):
    """Write histograms into a directoty of the current simulation
    
    Parametrers:
        input_params (dict): Input patameters used to start the simulation (from TOML)
        stats (np.recarray[STATS_DTYPE]): Statistics to save
    
    Returns:
        str: Absolute path to the output directory
    """
    stats = stats[0]
    workdir = input_params["simulation"]["workdir"]
    out_path = Path(__file__).parent / workdir
    out_path.mkdir(parents=True, exist_ok=True)
    
    follow_recoils = input_params["simulation"]["follow_recoils"]
    elems, nelem_target = _get_element_names_from_input(input_params)
    for key, val in zip(stats.dtype.names, stats):
        if not val["score"]:
            continue
        if key in ["xy", "xyn", "xye"]:
            source_key = {"xy": "x", "xyn": "xn", "xye": "xe"}[key]
            _, header_elems = _build_distribution_headers(
                source_key, elems, follow_recoils, nelem_target
            )
            if header_elems is not None:
                _write_histogram_binary_2d(out_path / f"{key}.hisb", val, header_elems[1:])
            continue

        x_vals = np.linspace(val["limits"][0], val["limits"][1], val["nbins"])
        header_indexes, header_elems = _build_distribution_headers(
            key, elems, follow_recoils, nelem_target
        )
        if header_indexes is None or header_elems is None:
            continue

        species_labels = header_elems[1:]
        species_indexes = header_indexes[1:]
        _write_histogram(out_path / f"{key}.his", val, x_vals, header_indexes, header_elems)
        _write_moments(out_path / f"{key}.mom", val, species_labels, species_indexes)
        # _write_raw_moments(out_path / "moments_raw" / f"{key}.mom", val, species_labels, species_indexes)
    return out_path


def _resolve_out_path(input_params):
    workdir = input_params["simulation"]["workdir"]
    return Path(__file__).parent / workdir


def save_progress(input_params, ions_done, ions_total, status="running"):
    out_path = _resolve_out_path(input_params)
    out_path.mkdir(parents=True, exist_ok=True)
    with open(out_path / "progress", "w") as f:
        f.write(f"{int(ions_done)} / {int(ions_total)}")
    with open(out_path / "status", "w") as f:
        f.write(str(status))


def is_stop_requested(input_params):
    """Return True if the UI has requested a cooperative stop.

    The UI signals an abort by creating a ``stop_requested`` file inside the
    working directory; the chunked simulation loop polls for it and exits
    gracefully after the current chunk, leaving valid partial results.
    """
    try:
        return (_resolve_out_path(input_params) / "stop_requested").exists()
    except Exception:
        return False

# if __name__ == "__main__":
#     with open(Path(__file__).parent / "pickle_dump.p", "rb") as f:
#         _, stats, input_params, _ = pickle.load(f)
#     write_stats(input_params, stats)
