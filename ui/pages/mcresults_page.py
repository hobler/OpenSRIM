from __future__ import annotations

import os
from pathlib import Path
from typing import List, Dict, Any

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


def _read_moments_file(path: Path) -> dict:
    """Read a .mom file and return a dict of species → {total, mean, std, skewness, kurtosis}."""
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
    # Columns are: order, (moments, errors) pairs per species
    # labels[0] = "Moment order", then pairs
    col = 1
    label_idx = 1
    while col < data.shape[1] and label_idx < len(labels):
        species_name = labels[label_idx].replace(" moments", "")
        vals = {}
        for i, mname in enumerate(moment_names):
            if i < data.shape[0]:
                vals[mname] = float(data[i, col])
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

_COLORS = [
    "#3274A1", "#E1812C", "#3A923A", "#C03D3E", "#9372B2",
    "#8E6C8A", "#D97706", "#2E86AB", "#4C956C", "#B85C38",
]


def _build_plots_from_directory(results_dir: str):
    """Scan a results directory for .his/.mom files and build plot dicts.

    Returns
    -------
    plots : dict  – keyed by plot_id, compatible with MCResultsWidget.set_results()
    numerical_values : list[dict]
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
            def plot_func(ax):
                for i, col_data in enumerate(_cols):
                    label = _labels[i + 1] if (i + 1) < len(_labels) else f"Series {i}"
                    color = _COLORS[i % len(_COLORS)]
                    ax.step(_x, col_data, where="mid", color=color,
                            linewidth=1.2, label=label)
                ax.set_xlabel(_xlabel)
                ax.set_ylabel(_ylabel)
                ax.set_title(_title)
                if len(_cols) > 1:
                    ax.legend(fontsize=8)
            return plot_func

        def _make_stats_func(_key, _moments_dict):
            def stats_func():
                moms = _moments_dict.get(_key, {})
                result = []
                for species, vals in moms.items():
                    result.append((f"{species} total", f"{vals.get('total', 0):.0f}"))
                    result.append((f"{species} mean", f"{vals.get('mean', 0):.2f}"))
                    result.append((f"{species} std", f"{vals.get('std', 0):.2f}"))
                if not result:
                    result = [("No data", "")]
                return result[:4]
            return stats_func

        plots[key] = {
            "name": plot_name,
            "plot_func": _make_plot_func(x, cols, labels, xlabel, ylabel, plot_name),
            "projection": None,
            "stats_func": _make_stats_func(key, all_moments),
        }

    # Build numerical values from moments
    numerical: List[Dict[str, str]] = []

    # Read progress for total ions simulated
    progress_path = base / "progress"
    if progress_path.exists():
        try:
            text = progress_path.read_text(encoding="utf-8").strip()
            parts = text.split("/")
            if len(parts) == 2:
                numerical.append({"name": "Total Ions Simulated", "value": parts[1].strip()})
                numerical.append({"name": "Ions Completed", "value": parts[0].strip()})
        except Exception:
            pass

    # Add moment statistics for the primary distribution (x = depth)
    x_moms = all_moments.get("x", {})
    for species, vals in x_moms.items():
        numerical.append({"name": f"{species} – Total", "value": f"{vals.get('total', 0):.0f}"})
        numerical.append({"name": f"{species} – Mean Range", "value": f"{vals.get('mean', 0):.2f} Å"})
        numerical.append({"name": f"{species} – Straggling", "value": f"{vals.get('std', 0):.2f} Å"})
        if "skewness" in vals:
            numerical.append({"name": f"{species} – Skewness", "value": f"{vals['skewness']:.4f}"})
        if "kurtosis" in vals:
            numerical.append({"name": f"{species} – Kurtosis", "value": f"{vals['kurtosis']:.4f}"})

    return plots, numerical


class MCResultsPage(QWidget):
    advanced_requested = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        self._results_widget = MCResultsWidget()
        self._results_widget.advanced_requested.connect(self.advanced_requested)
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
