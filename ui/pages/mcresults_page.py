from __future__ import annotations

import os
from pathlib import Path
from typing import List, Dict, Any, Optional

import numpy as np
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QWidget, QVBoxLayout

try:
    from ui.pages.simulation.simulation_page import MCResultsWidget
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.pages.simulation.simulation_page import MCResultsWidget  # type: ignore


def _read_histogram_file(path: Path) -> tuple:
    """Read a .his file and return (column_labels, x_values, data_columns).

    Returns
    -------
    labels : list[str]
        Human-readable column labels from the header comments.
    x : np.ndarray
        First column (bin positions).
    cols : list[np.ndarray]
        Remaining columns (one per species).
    """
    labels: list[str] = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if ":" in text and "," not in text:
                labels.append(text.split(":", 1)[1].strip())

    data = np.loadtxt(path, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    x = data[:, 0]
    cols = [data[:, i] for i in range(1, data.shape[1])]
    return labels, x, cols


def _read_histogram_binary_2d(path: Path) -> tuple:
    """Read a 2D binary histogram (.hisb) file written by OpenTRIM.

    Returns
    -------
    counts : np.ndarray of shape (n_species, nx, ny)
    x_values : np.ndarray of shape (nx,)  – bin centers along the x-axis
    y_values : np.ndarray of shape (ny,)  – bin centers along the y-axis
    species : list[str]                   – per-species labels (length n_species)
    """
    with open(path, "rb") as f:
        version = np.fromfile(f, dtype="<u2", count=1)
        if version.size != 1 or int(version[0]) != 0x00fa:
            raise ValueError(f"Unsupported 2D histogram version in {path}")
        shape = np.fromfile(f, dtype="<u4", count=3)
        if shape.size != 3:
            raise ValueError(f"Invalid 2D histogram shape in {path}")
        n_species, nx, ny = (int(v) for v in shape)
        x_values = np.fromfile(f, dtype="<f8", count=nx)
        y_values = np.fromfile(f, dtype="<f8", count=ny)
        species_raw = np.fromfile(f, dtype="S32", count=n_species)
        counts = np.fromfile(f, dtype="<f8", count=n_species * nx * ny)

    if (x_values.size, y_values.size, species_raw.size) != (nx, ny, n_species):
        raise ValueError(f"Truncated 2D histogram metadata in {path}")
    if counts.size != n_species * nx * ny:
        raise ValueError(f"Truncated 2D histogram payload in {path}")

    counts = counts.reshape(n_species, nx, ny)
    species = [b.rstrip(b"\x00").decode("utf-8", errors="replace") for b in species_raw]
    return counts, x_values, y_values, species


def _read_moments_file(path: Path) -> dict:
    """Read a .mom file and return moments with uncertainties per species."""
    labels: list[str] = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if ":" in text and "," not in text:
                labels.append(text.split(":", 1)[1].strip())

    data = np.loadtxt(path, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    moment_names = ["total", "mean", "std", "skewness", "kurtosis"]
    results = {}
    # Columns are: order, (moments, errors) pairs per species.
    # labels[0] = "Moment order", then pairs
    col = 1
    label_idx = 1
    while col < data.shape[1] and label_idx < len(labels):
        species_name = labels[label_idx].replace(" moments", "")
        vals = {}
        for i, mname in enumerate(moment_names):
            if i < data.shape[0]:
                vals[mname] = float(data[i, col])
                if (col + 1) < data.shape[1]:
                    vals[f"{mname}_err"] = float(data[i, col + 1])
        results[species_name] = vals
        col += 2  # skip errors column
        label_idx += 2
    return results


# Mapping from file prefix to human-readable name and axis labels
_MEASURE_META = {
    "x":   ("Depth Distribution: Ion/Recoil",           "Depth (Å)",            "Counts"),
    "xn":  ("Depth: Nuclear Energy Deposition",         "Depth (Å)",            "Energy (eV)"),
    "xe":  ("Depth: Electronic Energy Deposition",      "Depth (Å)",            "Energy (eV)"),
    "y":   ("Lateral Distribution: Ion/Recoil",         "Lateral Position (Å)", "Counts"),
    "yn":  ("Lateral: Nuclear Energy Deposition",       "Lateral Position (Å)", "Energy (eV)"),
    "ye":  ("Lateral: Electronic Energy Deposition",    "Lateral Position (Å)", "Energy (eV)"),
    "be":  ("Backscattered: Energy",                    "Energy (keV)",         "Counts"),
    "ba":  ("Backscattered: Angle",                     "Angle (deg)",          "Counts"),
    "te":  ("Transmitted: Energy",                      "Energy (keV)",         "Counts"),
    "ta":  ("Transmitted: Angle",                       "Angle (deg)",          "Counts"),
}

# 2D distribution metadata: (display_prefix, x_label, y_label, colorbar_label)
_MEASURE_2D_META = {
    "xy":  ("2D Distribution: Ion/Recoil",            "Depth (Å)", "Lateral (Å)", "Counts"),
    "xyn": ("2D Distribution: Nuclear Energy Dep.",   "Depth (Å)", "Lateral (Å)", "Energy (eV)"),
    "xye": ("2D Distribution: Electronic Energy Dep.", "Depth (Å)", "Lateral (Å)", "Energy (eV)"),
}

_COLORS = [
    "#3274A1", "#E1812C", "#3A923A", "#C03D3E", "#9372B2",
    "#8E6C8A", "#D97706", "#2E86AB", "#4C956C", "#B85C38",
]


def _rebin_step(x: np.ndarray, cols, factor: int):
    """Combine adjacent bins by *factor* (counts are summed, x is averaged).

    Trailing bins that don't fill a full group are dropped. Returns
    ``(x_new, [col_new, ...])`` even if no change was made (``factor <= 1``).
    """
    if factor is None or factor <= 1:
        return x, list(cols)
    n = (len(x) // factor) * factor
    if n == 0:
        return x, list(cols)
    x_arr = np.asarray(x[:n], dtype=float)
    x_new = x_arr.reshape(-1, factor).mean(axis=1)
    cols_new = [
        np.asarray(c[:n], dtype=float).reshape(-1, factor).sum(axis=1)
        for c in cols
    ]
    return x_new, cols_new


def _parse_species(label: str) -> tuple:
    """Extract (element, kind) from species labels like 'Si (vacancies)'.

    Returns (element_name, kind) where kind is one of:
    'ion', 'vacancies', 'interstitials', or '' if no parenthetical suffix.
    """
    text = label.strip()
    if "(" in text and text.endswith(")"):
        elem, _, rest = text.partition("(")
        return elem.strip(), rest[:-1].strip().lower()
    return text, ""


def _build_atom_table(all_moments: Dict[str, dict], total_ions: float) -> Dict[str, Any]:
    """Build a 2D table: rows are parameter names, columns are atom species.

    The grouping follows the user's request:
      Backscattered:  reflected (ion column), sputtered (target atom columns)
      Transmitted:    transmitted (ion column), forward-sputtered (target atoms)
      Inside:         stopped (ion column), vacancies, interstitials (target atoms)

    Returns a dict with keys: "atoms", "rows" where
       atoms : ordered list of atom labels (ion first, then targets)
       rows  : list of (row_name, {atom_label: value_str}) tuples
    """
    ion_label: Optional[str] = None
    target_atoms: List[str] = []

    def _register(label: str, *, is_ion: bool) -> None:
        nonlocal ion_label
        if is_ion and ion_label is None:
            ion_label = label
            return
        if not is_ion and label not in target_atoms and label != ion_label:
            target_atoms.append(label)

    for key in ("x", "be", "te"):
        for species in all_moments.get(key, {}).keys():
            elem, kind = _parse_species(species)
            _register(elem, is_ion=(kind == "ion"))

    atoms: List[str] = []
    if ion_label:
        atoms.append(ion_label)
    atoms.extend(target_atoms)
    if not atoms:
        return {"atoms": [], "rows": []}

    def _zeros() -> Dict[str, float]:
        return {a: 0.0 for a in atoms}

    refl = _zeros()
    sput = _zeros()
    trans = _zeros()
    fsput = _zeros()
    stopped = _zeros()
    vac = _zeros()
    inter = _zeros()

    for species, vals in all_moments.get("be", {}).items():
        elem, kind = _parse_species(species)
        total = float(vals.get("total", 0.0) or 0.0)
        if elem not in atoms:
            continue
        if kind == "ion":
            refl[elem] += total
        else:
            sput[elem] += total

    for species, vals in all_moments.get("te", {}).items():
        elem, kind = _parse_species(species)
        total = float(vals.get("total", 0.0) or 0.0)
        if elem not in atoms:
            continue
        if kind == "ion":
            trans[elem] += total
        else:
            fsput[elem] += total

    for species, vals in all_moments.get("x", {}).items():
        elem, kind = _parse_species(species)
        total = float(vals.get("total", 0.0) or 0.0)
        if elem not in atoms:
            continue
        if kind == "ion":
            stopped[elem] += total
        elif kind == "vacancies":
            vac[elem] += total
        elif kind == "interstitials":
            inter[elem] += total

    def _row(name: str, values: Dict[str, float], *, as_pct: bool = False) -> tuple:
        cells = {}
        for a in atoms:
            v = values.get(a, 0.0)
            if as_pct and total_ions > 0:
                cells[a] = f"{(100.0 * v / total_ions):.3f} %"
            elif v == 0.0:
                cells[a] = "—"
            else:
                cells[a] = f"{v:.0f}"
        return (name, cells)

    rows: List[tuple] = [
        ("— Backscattered —", {a: "" for a in atoms}),
        _row("Reflected", refl),
        _row("Sputtered", sput),
        ("— Transmitted —", {a: "" for a in atoms}),
        _row("Transmitted", trans),
        _row("Forward sputtered", fsput),
        ("— Inside (1 − backscatt. − transm.) —", {a: "" for a in atoms}),
        _row("Stopped in target", stopped),
        _row("Vacancies", vac),
        _row("Interstitials", inter),
    ]

    if total_ions > 0:
        rows.append(("— Yields (% of total ions) —", {a: "" for a in atoms}))
        rows.append(_row("Reflected yield", refl, as_pct=True))
        rows.append(_row("Sputter yield", sput, as_pct=True))
        rows.append(_row("Transmitted yield", trans, as_pct=True))
        rows.append(_row("Forward sputter yield", fsput, as_pct=True))

    return {"atoms": atoms, "rows": rows}


def _build_plots_from_directory(results_dir: str):
    """Scan a results directory for .his/.mom files and build plot dicts.

    Returns
    -------
    plots : dict  – keyed by plot_id, compatible with MCResultsWidget.set_results()
    numerical_values : dict with keys 'rows' (single-column param/value rows) and
                       'atom_table' (per-atom 2D structure)
    """
    base = Path(results_dir)
    plots: Dict[str, dict] = {}
    all_moments: Dict[str, dict] = {}

    # Discover all histogram files
    his_files = sorted(base.glob("*.his"))
    for his_path in his_files:
        key = his_path.stem  # e.g. "x", "xn", "be"
        mom_path = base / f"{key}.mom"

        meta = _MEASURE_META.get(key)
        if meta is None:
            continue
        plot_name, xlabel, ylabel = meta

        try:
            labels, x, cols = _read_histogram_file(his_path)
        except Exception:
            continue

        if mom_path.exists():
            try:
                all_moments[key] = _read_moments_file(mom_path)
            except Exception:
                pass

        # Build plot function (closure)
        def _make_plot_func(_x, _cols, _labels, _xlabel, _ylabel, _title):
            def plot_func(ax, bin_factor: int = 1):
                x_p, cols_p = _rebin_step(_x, _cols, int(bin_factor or 1))
                for i, col_data in enumerate(cols_p):
                    label = _labels[i + 1] if (i + 1) < len(_labels) else f"Series {i}"
                    color = _COLORS[i % len(_COLORS)]
                    ax.step(x_p, col_data, where="mid", color=color,
                            linewidth=1.2, label=label)
                ax.set_xlabel(_xlabel)
                ax.set_ylabel(_ylabel)
                ax.set_title(_title)
                if len(cols_p) > 1:
                    ax.legend(fontsize=8)
            return plot_func

        # Build series labels (one per column)
        series_labels = [
            labels[i + 1] if (i + 1) < len(labels) else f"Series {i}"
            for i in range(len(cols))
        ]

        def _make_stats_func(_key, _moments_dict):
            def stats_func():
                moms = _moments_dict.get(_key, {})
                result = []
                for species, vals in moms.items():
                    total = float(vals.get("total", 0.0) or 0.0)
                    total_err = float(vals.get("total_err", 0.0) or 0.0)
                    mean = float(vals.get("mean", 0.0) or 0.0)
                    mean_err = float(vals.get("mean_err", 0.0) or 0.0)
                    std = float(vals.get("std", 0.0) or 0.0)
                    std_err = float(vals.get("std_err", 0.0) or 0.0)
                    result.append((f"{species} total", f"{total:.0f} ± {total_err:.0f}"))
                    result.append((f"{species} mean", f"{mean:.2f} ± {mean_err:.2f}"))
                    result.append((f"{species} std", f"{std:.2f} ± {std_err:.2f}"))
                    if "skewness" in vals:
                        result.append((f"{species} skewness", f"{float(vals.get('skewness', 0.0)):.4f}"))
                    if "kurtosis" in vals:
                        result.append((f"{species} kurtosis", f"{float(vals.get('kurtosis', 0.0)):.4f}"))
                if not result:
                    result = [("No data", "")]
                return result
            return stats_func

        plots[key] = {
            "name": plot_name,
            "plot_func": _make_plot_func(x, cols, labels, xlabel, ylabel, plot_name),
            "projection": None,
            "stats_func": _make_stats_func(key, all_moments),
            "x_values": x,
            "columns": cols,
            "series_labels": series_labels,
            "x_label": xlabel,
            "y_label": ylabel,
            "is_depth_profile": key in ("x", "xn", "xe"),
            "colors": [_COLORS[i % len(_COLORS)] for i in range(len(cols))],
        }

    # Discover all 2D binary histograms (.hisb)
    for hisb_path in sorted(base.glob("*.hisb")):
        key = hisb_path.stem  # "xy", "xyn", "xye"
        meta = _MEASURE_2D_META.get(key)
        if meta is None:
            continue
        plot_prefix, xlabel, ylabel, cbar_label = meta

        try:
            counts2d, x_edges, y_edges, species = _read_histogram_binary_2d(hisb_path)
        except Exception:
            continue

        for ispec, spec_label in enumerate(species):
            data = counts2d[ispec]
            if not np.any(data):
                continue

            plot_id = f"{key}__{ispec}"
            plot_name = f"{plot_prefix} ({spec_label})"

            def _make_2d_plot_func(_x, _y, _data, _xlabel, _ylabel, _title, _cbar):
                def plot_func(ax, bin_factor: int = 1, conv=None):
                    # x_edges/y_edges from the binary file are bin centers
                    # (linspace over limits with `nbins` points). pcolormesh
                    # treats them as cell-edge coordinates with shading="auto".
                    data = _data
                    if conv:
                        sigma = float(conv.get("sigma", 0.0) or 0.0)
                        sa = conv.get("scan_area") or [0.0, 0.0]
                        a = float(sa[0]) if len(sa) > 0 else 0.0
                        b = float(sa[1]) if len(sa) > 1 else 0.0
                        if sigma > 0 and _x.size > 1 and _y.size > 1:
                            from scipy.ndimage import gaussian_filter1d, convolve1d
                            from scipy.special import erf
                            dx = (_x[-1] - _x[0]) / (_x.size - 1)
                            dy = (_y[-1] - _y[0]) / (_y.size - 1)
                            # x (depth) direction: plain Gaussian smear.
                            data = gaussian_filter1d(
                                _data, sigma=sigma / dx, axis=0, mode="constant"
                            )
                            # y (lateral) direction: erf-rectangle if a scan
                            # area is set (prof's formula), otherwise Gauss.
                            if abs(b - a) < 1e-12:
                                data = gaussian_filter1d(
                                    data, sigma=sigma / dy, axis=1, mode="constant"
                                )
                            else:
                                sqrt2s = np.sqrt(2.0) * sigma
                                reach = 5.0 * sigma + max(abs(a), abs(b))
                                max_half = max(3, (_y.size - 1) // 2)
                                half = min(max(3, int(np.ceil(reach / dy))), max_half)
                                d = np.arange(-half, half + 1, dtype=float) * dy
                                kernel = 0.5 * (
                                    erf((d - a) / sqrt2s) - erf((d - b) / sqrt2s)
                                ) * dy
                                data = convolve1d(data, kernel, axis=1, mode="constant")
                    mesh = ax.pcolormesh(
                        _x, _y, data.T,
                        shading="auto", cmap="viridis",
                    )
                    ax.set_xlabel(_xlabel)
                    ax.set_ylabel(_ylabel)
                    ax.set_title(_title)
                    cbar = ax.figure.colorbar(mesh, ax=ax)
                    cbar.set_label(_cbar)
                return plot_func

            def _make_2d_stats_func(_data, _spec_label):
                def stats_func():
                    total = float(_data.sum())
                    peak = float(_data.max())
                    if total <= 0.0:
                        return [(f"{_spec_label} total", "0")]
                    return [
                        (f"{_spec_label} total", f"{total:.0f}"),
                        (f"{_spec_label} peak bin", f"{peak:.3g}"),
                    ]
                return stats_func

            plots[plot_id] = {
                "name": plot_name,
                "plot_func": _make_2d_plot_func(
                    x_edges, y_edges, data, xlabel, ylabel, plot_name, cbar_label,
                ),
                "projection": None,
                "stats_func": _make_2d_stats_func(data, spec_label),
                # Marker: indicates a 2D plot. Used by callers (e.g. the single
                # plot tab) to decide whether the dataset is overlayable.
                "is_2d": True,
                # Raw 2-D axes for downstream convolution / σ-default heuristics.
                "x_values_2d": x_edges,
                "y_values_2d": y_edges,
                "x_label": xlabel,
                "y_label": ylabel,
                "colors": [],
            }

    # Build numerical values from moments.
    summary_rows: List[Dict[str, str]] = []

    # Read progress for total ions simulated
    progress_path = base / "progress"
    total_ions = 0.0
    if progress_path.exists():
        try:
            text = progress_path.read_text(encoding="utf-8").strip()
            parts = text.split("/")
            if len(parts) == 2:
                summary_rows.append({"name": "Total Ions Simulated", "value": parts[1].strip()})
                summary_rows.append({"name": "Ions Completed", "value": parts[0].strip()})
                total_ions = float(parts[1].strip()) if parts[1].strip() else 0.0
        except Exception:
            pass

    # Build the per-atom 2D table (Backscattered/Transmitted/Inside × atom columns).
    atom_table = _build_atom_table(all_moments, total_ions)

    # Add per-species moments for the depth distribution as flat summary rows
    # (mean range, straggling, etc.).
    x_moms = all_moments.get("x", {})
    for species, vals in x_moms.items():
        mean = float(vals.get("mean", 0.0) or 0.0)
        mean_err = float(vals.get("mean_err", 0.0) or 0.0)
        std = float(vals.get("std", 0.0) or 0.0)
        std_err = float(vals.get("std_err", 0.0) or 0.0)
        summary_rows.append({"name": f"{species} – Mean Range",
                             "value": f"{mean:.2f} ± {mean_err:.2f} Å"})
        summary_rows.append({"name": f"{species} – Straggling",
                             "value": f"{std:.2f} ± {std_err:.2f} Å"})
        if "skewness" in vals:
            summary_rows.append({"name": f"{species} – Skewness",
                                 "value": f"{float(vals['skewness']):.4f}"})
        if "kurtosis" in vals:
            summary_rows.append({"name": f"{species} – Kurtosis",
                                 "value": f"{float(vals['kurtosis']):.4f}"})

    numerical = {"rows": summary_rows, "atom_table": atom_table}
    return plots, numerical


class MCResultsPage(QWidget):
    advanced_requested = pyqtSignal(str)
    plot_open_in_single = pyqtSignal(str, object)

    def __init__(self, parent=None):
        super().__init__(parent)
        # Match the darker / thicker group-box borders used on the MC Setup
        # and KORAL pages.
        self.setStyleSheet(
            "QGroupBox { border: 2px solid palette(shadow); border-radius: 4px;"
            " margin-top: 6px; padding-top: 6px; }"
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        self._results_widget = MCResultsWidget()
        self._results_widget.advanced_requested.connect(self.advanced_requested)
        self._results_widget.plot_open_in_single.connect(self.plot_open_in_single)
        self._results_widget.load_directory_requested.connect(
            self.load_results_from_directory
        )
        layout.addWidget(self._results_widget)

    def get_results_widget(self):
        return self._results_widget

    def load_results_from_directory(self, results_dir: str, *, silent: bool = False) -> None:
        """Read OpenTRIM output files from *results_dir* and display them.

        Parameters
        ----------
        silent : bool
            If *True*, suppress error/info dialogs (used for live updates
            during a running simulation).
        """
        try:
            plots, numerical = _build_plots_from_directory(results_dir)
        except Exception as exc:
            if not silent:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Load Results",
                                    f"Failed to load results from:\n{results_dir}\n\n{exc}")
            return

        if not plots:
            if not silent:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.information(self, "Load Results",
                                        "No result files found in the output directory.")
            return

        self._results_widget.set_results(plots, numerical)

        # If only a single histogram is present, skip straight to the single
        # plot tab — there's nothing to compare or arrange. Only on real
        # loads (not the silent live updates during a running simulation),
        # otherwise the view would keep jumping mid-run.
        if not silent and len(plots) == 1:
            only_id, only_info = next(iter(plots.items()))
            self.plot_open_in_single.emit(only_id, only_info)
