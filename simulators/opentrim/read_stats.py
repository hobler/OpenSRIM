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
        if name in ["trajectories", "distribution_2d"]:
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


def _build_2d_measure_configs(input_params):
    names = {
        "ion_recoils": "xy",
        "nuclear_energy_deposition": "xyn",
        "electronic_energy_deposition": "xye",
    }
    _, nelem_target = _get_element_names_from_input(input_params)
    nelem = 1 + nelem_target
    follow_recoils = input_params["simulation"]["follow_recoils"]
    nvar = {
        "xy": nelem + nelem_target if follow_recoils else nelem,
        "xyn": nelem,
        "xye": nelem,
    }
    configs = {}
    distribution_2d = input_params["output"]["distribution_2d"]
    for name, key in names.items():
        cfg = distribution_2d[name]
        configs[key] = {
            "score": bool(cfg["score"]),
            "x_nbins": int(cfg["nbins"][0]),
            "y_nbins": int(cfg["nbins"][1]),
            "x_limits": tuple(cfg["limits"][0]),
            "y_limits": tuple(cfg["limits"][1]),
            "nvar": int(nvar[key]),
        }
    return configs


def _read_histogram_binary_2d(path, expected_nvar=None, expected_nx=None, expected_ny=None):
    with open(path, "rb") as f:
        version = np.fromfile(f, dtype="<u2", count=1)
        if version.size != 1:
            raise ValueError(f"Invalid binary histogram header in {path}")
        shape = np.fromfile(f, dtype="<u4", count=3)
        if shape.size != 3:
            raise ValueError(f"Invalid binary histogram shape header in {path}")
        n_species, nx, ny = [int(value) for value in shape]
        species_labels = np.fromfile(f, dtype="S32", count=n_species)
        tables = np.fromfile(f, dtype="<f8")

    if int(version[0]) != 0x00fa:
        raise ValueError(f"Unsupported binary histogram version in {path}: {int(version[0])}")
    if expected_nvar is not None and n_species != expected_nvar:
        raise ValueError(
            f"Invalid species count in {path}: expected {expected_nvar}, got {n_species}"
        )
    if expected_nx is not None and expected_ny is not None and (nx, ny) != (expected_nx, expected_ny):
        raise ValueError(
            f"Invalid binary histogram shape in {path}: expected {(expected_nx, expected_ny)}, got {(nx, ny)}"
        )
    if species_labels.size != n_species:
        raise ValueError(f"Invalid binary histogram species metadata in {path}")
    expected_size = n_species * (nx + 1) * (ny + 1)
    if tables.size != expected_size:
        raise ValueError(
            f"Invalid binary histogram payload size in {path}: expected {expected_size}, got {tables.size}"
        )
    tables = tables.reshape(n_species, nx + 1, ny + 1)
    labels = [label.rstrip(b"\x00").decode("utf-8") for label in species_labels]
    return tables, labels


def _build_stats_2d_dtype(configs_2d):
    fields = []
    for key, cfg in configs_2d.items():
        if not cfg["score"]:
            continue
        measure_dtype = np.dtype(
            [
                ("hist", np.float64, (cfg["nvar"], cfg["x_nbins"] + 1, cfg["y_nbins"] + 1)),
            ],
            align=True,
        )
        fields.append((key, measure_dtype))
    return np.dtype(fields, align=True)


def _measure_label_2d(measure_key):
    labels = {
        "xy": "Depth-lateral ion/recoil distribution",
        "xyn": "Depth-lateral nuclear energy deposition",
        "xye": "Depth-lateral electronic energy deposition",
    }
    return labels.get(measure_key, measure_key)


def _read_stats_2d(base_path, configs_2d):
    stats_2d = np.zeros(1, dtype=_build_stats_2d_dtype(configs_2d))
    labels_2d = {
        "x_values": {"pos": [1, 0], "type": "column"},
        "y_values": {"pos": [0, 1], "type": "row"},
        "metrics": {},
    }

    for key, cfg in configs_2d.items():
        if not cfg["score"]:
            continue
        path = base_path / f"{key}.hisb"
        if not path.exists():
            raise FileNotFoundError(
                f"Missing binary histogram file for active 2D metric '{key}': {path}"
            )
        tables, species = _read_histogram_binary_2d(
            path, cfg["nvar"], cfg["x_nbins"], cfg["y_nbins"]
        )
        stats_2d[0][key]["hist"] = tables
        labels_2d["metrics"][_measure_label_2d(key)] = {
            "hist": {
                "shape": list(tables.shape),
                "species": species,
            }
        }

    return stats_2d, labels_2d


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
    configs_2d = _build_2d_measure_configs(input_params)
    stats_1d = np.zeros(1, dtype=_build_stats_dtype(configs))
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
        stats_1d[0][key]["hist"] = hist_data

        mom_data = np.loadtxt(mom_path, delimiter=",")
        if mom_data.ndim == 1:
            mom_data = mom_data.reshape(1, -1)
        expected_mom_shape = (nmom, 1 + 2 * cfg["nvar"])
        if mom_data.shape != expected_mom_shape:
            raise ValueError(
                f"Invalid moments data shape for '{key}' in {mom_path}: "
                f"expected {expected_mom_shape}, got {mom_data.shape}"
            )
        stats_1d[0][key]["mom"] = mom_data

    labels_1d = _build_helper_labels(configs, base_path)
    stats_2d, labels_2d = _read_stats_2d(base_path, configs_2d)
    return stats_1d, labels_1d, stats_2d, labels_2d

if __name__ == "__main__":
    input_params = read_params("defaults.toml")
    stats_1d, labels_1d, stats_2d, labels_2d = read_stats(input_params, include_unscored=True)
    import pprint
    pprint.pprint(stats_1d.dtype.descr)
    pprint.pprint(labels_1d)
    pprint.pprint(stats_2d.dtype.descr)
    pprint.pprint(labels_2d)
