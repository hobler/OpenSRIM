import sys
from pathlib import Path
import numpy as np
from .stats import standardize_moments


def _build_distribution_headers(field_name, params):
    """Build headers for statistics data files.
    
    Parameters:
        field_name: (str) Stats field name
        params: (PARAMS_DTYPE) Simulation parameters
    """
    elem_names = [params.elements[ielem]["name"] 
                  for ielem in range(params.nelem)]

    if field_name in ["i", "x", "y"]:
        if params.cascade.follow_recoils:
            elem_names += elem_names[1:]  # add target elements for vacancies

        species_labels = elem_names[:]
        header_indexes = [0]
        for ielem, elem_name in enumerate(species_labels):
            if ielem == 0:
                header_indexes.append(params.elements[ielem]["Z"])
                species_labels[ielem] = f"{elem_name} (ion)"
            elif ielem <= params.nelem_target:
                header_indexes.append(params.elements[ielem]["Z"])
                species_labels[ielem] = f"{elem_name} (interstitial)"
            else:
                jelem = ielem - params.nelem_target
                header_indexes.append(-params.elements[jelem]["Z"])
                species_labels[ielem] = f"{elem_name} (vacancy)"
        if field_name == "i":
            return header_indexes, ["Yield"] + species_labels
        else:
            return (header_indexes, [f"{field_name.upper()} Position"] 
                                     + species_labels)

    if field_name[0] in ["b", "t"]:
        follow_recoils = params.cascade.follow_recoils
        species_labels = elem_names[:] if follow_recoils else elem_names[:1]
        header_indexes = [0]
        for ielem, species_label in enumerate(species_labels):
            header_indexes.append(params.elements[ielem]["Z"])
            if ielem == 0:
                species_label += " (ion)"
            else:
                species_label += " (target atom)"
        if field_name[0] == "b":
            first_label = "Backscattered " 
        else:
            first_label = "Transmitted "
        if len(field_name) == 1:
            first_label += " yield"
        elif field_name[1] == "e":
            first_label += "energy"
        elif field_name[1] == "a":
            first_label += "angle"
        return header_indexes, [first_label] + species_labels

    if field_name[1] in ["n", "e"]:
        species_labels = elem_names[:]
        header_indexes = [0]
        for ielem, elem_name in enumerate(species_labels):
            header_indexes.append(params.elements[ielem]["Z"])
            if ielem == 0:
                species_labels[ielem] = f"{elem_name} (ion)"
            else:
                species_labels[ielem] = f"{elem_name} (target atom)"
        return (header_indexes, [f"{field_name[0].upper()} Position"] 
                                 + species_labels)

    return None, None


def _header_and_column_row(indexes, labels):
    column_row = ", ".join([str(el) for el in indexes])
    header = "\n".join([f"{k}: {v}" for k, v in zip(indexes, labels)])
    return header, column_row


def _write_histogram(path, stats_field, nion, header_indexes, header_labels):
    """Write a 1D histogram file.
    
    Units of densities are [1/A], [1/eV], or [1/deg] depending on the histogram 
    type.

    The number of ions processed is saved in the ion column of the last row of 
    the histogram file. I may be used to estimate the statistical error of the 
    histogram data.
    
    Parameters:
        path: (str) The path to the output file.
        stats_field: (dict) The statistics field containing the histogram data.
        nion: (int) The number of projectiles already processed.
        header_indexes: (list) The indexes for the header.
        header_labels: (list) The labels for the header.
    """
    x_vals = np.linspace(stats_field["limits"][0], stats_field["limits"][1], 
                         stats_field["nbins"] + 1)
    header, column_row = _header_and_column_row(header_indexes, header_labels)
    densities = stats_field["counts"][:, 1:] / (nion * stats_field["bin_width"])
    densities[0, -1] = nion
    densities[1:, -1] = 0.0
    data = np.vstack((x_vals, densities)).T
    np.savetxt(path, data, delimiter=", ", fmt="%13.6e", 
               header=header + "\n" + column_row)


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


def _write_moments(path, val, species_indexes, species_labels):
    metric_names = ["total", "mean", "std", "skewness", "kurtosis"]
    max_order_local = (val["power_sums"].shape[1] - 1) // 2
    n_metrics = min(max_order_local + 1, len(metric_names))

    header_indexes = [0]
    header_labels = ["Moment order"]
    data_cols = [np.arange(n_metrics, dtype=np.float64)]
    for ivar, label in enumerate(species_labels):
        std_moments, std_err_moments = standardize_moments(val, ivar)
        val_idx = species_indexes[ivar]

        header_indexes.append(val_idx)
        header_labels.append(f"{label} moments")
        data_cols.append(std_moments[:n_metrics])

        header_indexes.append(val_idx)
        header_labels.append(f"{label} moments errors")
        data_cols.append(std_err_moments[:n_metrics])

    header, column_row = _header_and_column_row(header_indexes, header_labels)
    data = np.column_stack(data_cols)
    fmt = "%2i" + " %13.6e" * (len(data_cols) - 1)
    np.savetxt(path, data, fmt=fmt, header=header + "\n" + column_row)


def _write_raw_moments(path, val, species_indexes, species_labels):
    raw_nmom = val["power_sums"].shape[1]
    header_indexes = [0]
    header_labels = ["Moment order"]
    data_cols = [np.arange(raw_nmom, dtype=np.float64)]
    for ivar, label in enumerate(species_labels):
        header_indexes.append(species_indexes[ivar])
        header_labels.append(f"{label} raw moments")
        data_cols.append(val["power_sums"][ivar, :])

    header, column_row = _header_and_column_row(header_indexes, header_labels)
    data = np.column_stack(data_cols)
    fmt = "%2i" + " %13.6e" * (len(data_cols) - 1)
    np.savetxt(path, data, fmt=fmt, header=header + "\n" + column_row)


def write_stats(params, stats, nion, workdir):
    """Write histograms into a directory of the current simulation
    
    Parameters:
        params:  (PARAMS_DTYPE) Simulation parameters
        stats: (np.recarray[STATS_DTYPE]) Statistics to save
        nion: (int) Number of projectiles already processed
        workdir: (str) Directory to save the output files
    
    Returns:
        str: Absolute path to the output directory
    """
    out_path = Path(__file__).parent / workdir
    out_path.mkdir(parents=True, exist_ok=True)
    
    for stats_name in stats.dtype.names:
        stats_field = stats[stats_name]

        if not stats_field["score"]:
            continue

        if stats_name not in ["xy", "xyn", "xye"]:  # 0d and 1d
            header_indexes, header_labels = (
                _build_distribution_headers(stats_name, params)
            )
            if header_indexes is None or header_labels is None:
                continue
            species_labels = header_labels[1:]
            species_indexes = header_indexes[1:]
            _write_moments(out_path / f"{stats_name}.mom", stats_field, 
                           species_indexes, species_labels)
            # _write_raw_moments(out_path / "moments_raw" / f"{stats_name}.mom", 
            #                    stats_field, species_indexes, species_labels)
            if stats_name not in ["i", "b", "t"]:  # 1d
                _write_histogram(out_path / f"{stats_name}.his", stats_field, 
                                 nion, header_indexes, header_labels)
        else:  # 2d
            map = {"xy": "x", "xyn": "xn", "xye": "xe"}
            source_stats_name = map[stats_name]
            _, header_labels = (
                _build_distribution_headers(source_stats_name, params)
            )
            if header_labels is not None:
                _write_histogram_binary_2d(out_path / f"{stats_name}.hisb", 
                                           stats_field, header_labels[1:])
    
    return out_path


def save_progress(workdir, ions_done, ions_total):
    out_path = Path(__file__).parent / workdir
    out_path.mkdir(parents=True, exist_ok=True)
    with open(out_path / "progress", "w") as f:
        f.write(f"{int(ions_done)} / {int(ions_total)}")

# if __name__ == "__main__":
#     with open(Path(__file__).parent / "pickle_dump.p", "rb") as f:
#         _, stats, input_params, _ = pickle.load(f)
#     write_stats(input_params, stats)
