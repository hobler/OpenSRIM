from pathlib import Path
import sys
if __package__ is None:
    # Running as a script: add parent dir to sys.path
    project_root = str(Path(__file__).parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    __package__ = str(Path(__file__).parent.name)
import numpy as np
from .stats import max_order
from .read_params import read_params
from .process_data import _get_element_names_from_input


def _build_measure_configs(input_params, include_unscored=False):
    short_names = {
        "depth_distribution": "x",
        "lateral_distribution": "y",
        "backscattered_atoms": "b",
        "transmitted_atoms": "t",
        "ion_recoils": "",
        "nuclear_energy_deposition": "n",
        "electronic_energy_deposition": "e",
        "energy": "e",
        "angle": "a",
    }

    _, nelem_target = _get_element_names_from_input(input_params)
    nelem_ion = 1
    nelem = nelem_ion + nelem_target
    follow_recoils = input_params["simulation"]["follow_recoils"]
    nvar = {
        "x": nelem + nelem_target if follow_recoils else nelem,
        "y": nelem + nelem_target if follow_recoils else nelem,
        "xn": nelem,
        "yn": nelem,
        "xe": nelem,
        "ye": nelem,
        "be": nelem if follow_recoils else nelem_ion,
        "ba": nelem if follow_recoils else nelem_ion,
        "te": nelem if follow_recoils else nelem_ion,
        "ta": nelem if follow_recoils else nelem_ion,
    }

    configs = {}
    for name, group in input_params["output"].items():
        if name == "trajectories":
            continue
        for subname, cfg in group.items():
            key = f"{short_names[name]}{short_names[subname]}"
            score = bool(cfg["score"])
            if not include_unscored and not score:
                continue
            configs[key] = {
                "score": score,
                "nbins": int(cfg["nbins"]),
                "limits": tuple(cfg["limits"]),
                "nvar": int(nvar[key]),
            }
    return configs


def _build_stats_dtype(configs):
    nmom = max_order + 1
    fields = []
    for key, cfg in configs.items():
        measure_dtype = np.dtype(
            [
                ("hist", np.int32, (cfg["nbins"], cfg["nvar"] + 1)),
                ("mom", np.float64, (nmom, 1 + 2 * cfg["nvar"])),
            ],
            align=True,
        )
        fields.append((key, measure_dtype))
    return np.dtype(fields, align=True)


def _measure_label(measure_key):
    labels = {
        "x": "Depth distribution",
        "y": "Lateral distribution",
        "xn": "Depth nuclear energy deposition",
        "yn": "Lateral nuclear energy deposition",
        "xe": "Depth electronic energy deposition",
        "ye": "Lateral electronic energy deposition",
        "be": "Backscattered energy distribution",
        "ba": "Backscattered angular distribution",
        "te": "Transmitted energy distribution",
        "ta": "Transmitted angular distribution",
    }
    return labels.get(measure_key, measure_key)


def _read_column_labels(path):
    labels = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if ":" in text and "," not in text:
                labels.append(text)
    return labels


def _read_histogram_binary_2d(path, expected_nvar, expected_nx, expected_ny):
    with open(path, "rb") as f:
        n_species = np.fromfile(f, dtype="<u4", count=1)
        if n_species.size != 1:
            raise ValueError(f"Invalid binary histogram header in {path}")
        version = np.fromfile(f, dtype="<u2", count=1)
        shape = np.fromfile(f, dtype="<u4", count=2)
        x_limits = np.fromfile(f, dtype="<f8", count=2)
        y_limits = np.fromfile(f, dtype="<f8", count=2)
        bin_widths = np.fromfile(f, dtype="<f8", count=2)
        species_ids = np.fromfile(f, dtype="<i4", count=int(n_species[0]))
        counts = np.fromfile(f, dtype="<f8")

    if int(version[0]) != 0x00fa:
        raise ValueError(f"Unsupported binary histogram version in {path}: {int(version[0])}")
    if int(n_species[0]) != expected_nvar:
        raise ValueError(
            f"Invalid species count in {path}: expected {expected_nvar}, got {int(n_species[0])}"
        )
    if tuple(shape) != (expected_nx, expected_ny):
        raise ValueError(
            f"Invalid binary histogram shape in {path}: expected {(expected_nx, expected_ny)}, got {tuple(shape)}"
        )
    expected_size = expected_nvar * (expected_nx + 2) * (expected_ny + 2)
    if counts.size != expected_size:
        raise ValueError(
            f"Invalid binary histogram payload size in {path}: expected {expected_size}, got {counts.size}"
        )
    return counts.reshape(expected_nvar, expected_nx + 2, expected_ny + 2), x_limits, y_limits, bin_widths, species_ids


def _hist_from_binary_pair(base_path, key, configs):
    pairs = {
        "x": ("xy", "x", "y"),
        "y": ("xy", "x", "y"),
        "xn": ("xyn", "xn", "yn"),
        "yn": ("xyn", "xn", "yn"),
        "xe": ("xye", "xe", "ye"),
        "ye": ("xye", "xe", "ye"),
    }
    if key not in pairs:
        return None

    binary_key, x_key, y_key = pairs[key]
    path = base_path / f"{binary_key}.hisb"
    if not path.exists() or x_key not in configs or y_key not in configs:
        return None

    x_cfg = configs[x_key]
    y_cfg = configs[y_key]
    counts, _, _, _, _ = _read_histogram_binary_2d(
        path, x_cfg["nvar"], x_cfg["nbins"], y_cfg["nbins"]
    )
    if key == x_key:
        values = np.linspace(x_cfg["limits"][0], x_cfg["limits"][1], x_cfg["nbins"])
        marginal = counts[:, 1:-1, :].sum(axis=2)
    else:
        values = np.linspace(y_cfg["limits"][0], y_cfg["limits"][1], y_cfg["nbins"])
        marginal = counts[:, :, 1:-1].sum(axis=1)
    return np.vstack((values, marginal)).T


def _build_helper_labels(configs, base_path):
    helper = []
    for key in configs:
        hist_path = base_path / f"{key}.his"
        mom_path = base_path / f"{key}.mom"
        hist_labels = _read_column_labels(hist_path) if hist_path.exists() else []
        mom_labels = _read_column_labels(mom_path) if mom_path.exists() else []
        helper.append(
            (
                _measure_label(key),
                [("hist", hist_labels), ("mom", mom_labels)],
            )
        )
    return helper


def read_stats(input_params, include_unscored=False):
    """Read histograms and moments for configured measures."""
    path = Path(input_params["simulation"]["workdir"])
    if path.is_absolute():
        base_path = path
    else:
        base_path = Path(__file__).parent / path

    if not base_path.exists():
        raise FileNotFoundError(f"Statistics path does not exist: {base_path}")
    if not base_path.is_dir():
        raise NotADirectoryError(f"Statistics path is not a directory: {base_path}")

    configs = _build_measure_configs(
        input_params, include_unscored=include_unscored
    )
    stats = np.zeros(1, dtype=_build_stats_dtype(configs))
    nmom = max_order + 1

    for key, cfg in configs.items():
        if not cfg["score"]:
            continue

        hist_path = base_path / f"{key}.his"
        mom_path = base_path / f"{key}.mom"
        if not mom_path.exists():
            raise FileNotFoundError(
                f"Missing moments file for active metric '{key}': {mom_path}"
            )

        hist_data = _hist_from_binary_pair(base_path, key, configs)
        if hist_data is None:
            if not hist_path.exists():
                raise FileNotFoundError(
                    f"Missing histogram file for active metric '{key}': {hist_path}"
                )
            hist_data = np.loadtxt(hist_path, delimiter=",")
            if hist_data.ndim == 1:
                hist_data = hist_data.reshape(1, -1)
        expected_hist_shape = (cfg["nbins"], cfg["nvar"] + 1)
        if hist_data.shape != expected_hist_shape:
            raise ValueError(
                f"Invalid histogram data shape for '{key}' in {hist_path}: "
                f"expected {expected_hist_shape}, got {hist_data.shape}"
            )
        stats[0][key]["hist"] = hist_data

        mom_data = np.loadtxt(mom_path, delimiter=",")
        if mom_data.ndim == 1:
            mom_data = mom_data.reshape(1, -1)
        expected_mom_shape = (nmom, 1 + 2 * cfg["nvar"])
        if mom_data.shape != expected_mom_shape:
            raise ValueError(
                f"Invalid moments data shape for '{key}' in {mom_path}: "
                f"expected {expected_mom_shape}, got {mom_data.shape}"
            )
        stats[0][key]["mom"] = mom_data

    return stats, _build_helper_labels(configs, base_path)

if __name__ == "__main__":
    input_params = read_params("input.toml")
    stats, helper_labels = read_stats(input_params)
    import pprint
    pprint.pprint(stats.dtype.descr)
    pprint.pprint(helper_labels)
