import csv
import json
import os
import sys
import numpy as np
from pathlib import Path
from typing import List, Dict, Any, Callable, Optional

from PyQt6.QtCore import Qt, pyqtSignal, QTimer, QMimeData, QPoint, QSize
from PyQt6.QtGui import QAction, QDrag, QIcon, QPixmap, QPainter, QColor, QPen, QCursor
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QSplitter, QScrollArea, QFrame, QGridLayout, QDialog,
    QListWidget, QListWidgetItem, QTableWidget, QTableWidgetItem,
    QHeaderView, QSizePolicy, QPushButton, QGroupBox, QAbstractItemView,
    QToolButton, QStyle, QFileDialog, QStackedWidget, QComboBox,
    QLineEdit, QMessageBox, QCheckBox, QDoubleSpinBox, QSpinBox,
    QColorDialog, QFormLayout, QToolTip,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.image as mpimg
from state import get_last_used_directory, remember_last_used_path

try:
    from ui.widgets.toggle_switch import ToggleSwitch
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.widgets.toggle_switch import ToggleSwitch  # type: ignore

try:
    from ui.widgets.hints_popup import HintSystem
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.widgets.hints_popup import HintSystem  # type: ignore

try:
    from ui.widgets.advanced_settings_button import AdvancedSettingsButton
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.widgets.advanced_settings_button import AdvancedSettingsButton  # type: ignore


# =====================================================================
#  HARDCODED DEFINITIONS — replace these with real simulation results
# =====================================================================

def _generate_hardcoded_data():
    """Generate OpenTRIM-like placeholder outputs and numerical values.

    Returns two dicts:
        AVAILABLE_PLOTS  – keyed by plot_id, each value contains:
            "name"      : display name
            "plot_func" : callable(ax) drawing the plot
            "stats_func": callable() -> list of 4 (name, value) pairs
        AVAILABLE_NUMERICAL_VALUES – global result key/value rows
    """
    rng = np.random.default_rng(42)

    n_total = 5000
    n_traj = 60

    # --- base trajectory cloud used by 3 trajectory outputs ---
    trajectories_xy: List[tuple[np.ndarray, np.ndarray]] = []
    collision_counts: List[int] = []
    segment_lengths: List[float] = []

    for _ in range(n_traj):
        n_pts = int(rng.integers(18, 70))
        steps = rng.normal(loc=0.0, scale=14.0, size=(n_pts, 2))
        xy = np.cumsum(steps, axis=0)
        x = xy[:, 0]
        y = xy[:, 1]
        trajectories_xy.append((x, y))

        collisions = max(0, n_pts - 1)
        collision_counts.append(collisions)
        if n_pts > 1:
            seg = np.hypot(np.diff(x), np.diff(y))
            segment_lengths.extend(seg.tolist())

    start_x = np.array([x[0] for x, _ in trajectories_xy if len(x)], dtype=float)
    start_y = np.array([y[0] for _, y in trajectories_xy if len(y)], dtype=float)
    end_x = np.array([x[-1] for x, _ in trajectories_xy if len(x)], dtype=float)
    end_y = np.array([y[-1] for _, y in trajectories_xy if len(y)], dtype=float)

    collision_x = np.concatenate(
        [x[1:-1] for x, _ in trajectories_xy if len(x) > 2],
        dtype=float,
    ) if any(len(x) > 2 for x, _ in trajectories_xy) else np.array([], dtype=float)
    collision_y = np.concatenate(
        [y[1:-1] for _, y in trajectories_xy if len(y) > 2],
        dtype=float,
    ) if any(len(y) > 2 for _, y in trajectories_xy) else np.array([], dtype=float)

    # --- range distributions ---
    ion_depths = rng.normal(loc=1500.0, scale=460.0, size=n_total)
    ion_depths = ion_depths[(ion_depths >= 0.0) & (ion_depths <= 4000.0)]

    recoil_depths = rng.exponential(scale=620.0, size=int(1.4 * n_total))
    recoil_depths = recoil_depths[recoil_depths <= 4000.0]

    photon_depths = rng.gamma(shape=3.5, scale=360.0, size=3200)
    photon_depths = photon_depths[(photon_depths >= 0.0) & (photon_depths <= 4000.0)]

    ionization_depth_axis = np.linspace(0.0, 4000.0, 140)
    ionization_range_profile = (
        23.0 * np.exp(-((ionization_depth_axis - 1350.0) ** 2) / (2.0 * 430.0**2))
        + 7.0 * np.exp(-((ionization_depth_axis - 2100.0) ** 2) / (2.0 * 700.0**2))
        + rng.normal(0.0, 0.6, size=ionization_depth_axis.size)
    )
    ionization_range_profile = np.clip(ionization_range_profile, 0.0, None)

    # --- lateral range distributions ---
    lateral_ion = rng.normal(loc=0.0, scale=220.0, size=n_total)
    lateral_recoil = rng.normal(loc=0.0, scale=310.0, size=int(1.2 * n_total))
    lateral_photons = rng.normal(loc=0.0, scale=360.0, size=3200)

    lateral_axis = np.linspace(-900.0, 900.0, 180)
    lateral_ionization_profile = (
        92.0 * np.exp(-(lateral_axis**2) / (2.0 * 250.0**2))
        + 16.0 * np.exp(-((lateral_axis - 260.0) ** 2) / (2.0 * 180.0**2))
        + rng.normal(0.0, 1.2, size=lateral_axis.size)
    )
    lateral_ionization_profile = np.clip(lateral_ionization_profile, 0.0, None)

    # --- boundary outputs ---
    backscatter_count = max(80, int(0.032 * n_total))
    backscatter_energy = rng.gamma(shape=2.2, scale=1700.0, size=backscatter_count)
    backscatter_energy = np.clip(backscatter_energy, 0.0, 50000.0)
    backscatter_angle = 90.0 + 90.0 * rng.beta(2.6, 1.8, size=backscatter_count)

    transmitted_count = max(80, int(0.018 * n_total))
    transmitted_energy = rng.normal(loc=32000.0, scale=5200.0, size=transmitted_count)
    transmitted_energy = np.clip(transmitted_energy, 0.0, 50000.0)
    transmitted_angle = 75.0 * rng.beta(1.8, 5.2, size=transmitted_count)

    def _set_equal_limits(ax, x_values: np.ndarray, y_values: np.ndarray) -> None:
        """Set a square data window so x/y units are not visually distorted."""
        if x_values.size == 0 or y_values.size == 0:
            ax.set_aspect("equal", adjustable="box")
            return

        x_min, x_max = float(np.min(x_values)), float(np.max(x_values))
        y_min, y_max = float(np.min(y_values)), float(np.max(y_values))

        span = max(x_max - x_min, y_max - y_min, 1.0)
        half = 0.5 * span * 1.05
        x_mid = 0.5 * (x_min + x_max)
        y_mid = 0.5 * (y_min + y_max)

        ax.set_xlim(x_mid - half, x_mid + half)
        ax.set_ylim(y_mid - half, y_mid + half)
        ax.set_aspect("equal", adjustable="box")

    def _safe_mean(values: np.ndarray) -> float:
        return float(np.mean(values)) if values.size else 0.0

    def _safe_std(values: np.ndarray) -> float:
        return float(np.std(values)) if values.size else 0.0

    def _safe_max(values: np.ndarray) -> float:
        return float(np.max(values)) if values.size else 0.0

    def _safe_percentile(values: np.ndarray, q: float) -> float:
        return float(np.percentile(values, q)) if values.size else 0.0

    def _safe_peak_center(values: np.ndarray, bins: int, limits: tuple[float, float]) -> float:
        if values.size == 0:
            return 0.0
        hist, edges = np.histogram(values, bins=bins, range=limits)
        idx = int(np.argmax(hist))
        return float(0.5 * (edges[idx] + edges[idx + 1]))

    # ---- plot functions ------------------------------------------------
    def plot_traj_start(ax):
        ax.scatter(start_x, start_y, s=18, alpha=0.75, color="#3274A1")
        _set_equal_limits(ax, start_x, start_y)
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        ax.set_title("Trajectories: Start")

    def plot_traj_end(ax):
        ax.scatter(end_x, end_y, s=18, alpha=0.75, color="#E1812C")
        _set_equal_limits(ax, end_x, end_y)
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        ax.set_title("Trajectories: End")

    def plot_traj_collisions(ax):
        ax.scatter(collision_x, collision_y, s=9, alpha=0.35, color="#C03D3E")
        _set_equal_limits(ax, collision_x, collision_y)
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        ax.set_title("Trajectories: Collisions")

    def plot_range_ion_recoil(ax):
        ax.hist(
            ion_depths,
            bins=60,
            range=(0.0, 4000.0),
            alpha=0.70,
            color="#3274A1",
            label="Ion",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.hist(
            recoil_depths,
            bins=60,
            range=(0.0, 4000.0),
            alpha=0.55,
            color="#E1812C",
            label="Recoil",
            edgecolor="black",
            linewidth=0.3,
        )
        ax.set_xlabel("Depth (A)")
        ax.set_ylabel("Counts")
        ax.set_title("Range Distribution: Ion/Recoil")
        ax.legend()

    def plot_range_photons(ax):
        ax.hist(
            photon_depths,
            bins=60,
            range=(0.0, 4000.0),
            alpha=0.82,
            color="#4C956C",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.set_xlabel("Depth (A)")
        ax.set_ylabel("Photon Counts")
        ax.set_title("Range Distribution: Photons")

    def plot_range_ionization(ax):
        ax.plot(
            ionization_depth_axis,
            ionization_range_profile,
            color="#8E6C8A",
            linewidth=1.6,
            label="Ionization",
        )
        ax.fill_between(
            ionization_depth_axis,
            ionization_range_profile,
            color="#8E6C8A",
            alpha=0.25,
        )
        ax.set_xlabel("Depth (A)")
        ax.set_ylabel("Ionization (arb. units)")
        ax.set_title("Range Distribution: Ionization")
        ax.legend()

    def plot_lateral_ion_recoil(ax):
        ax.hist(
            lateral_ion,
            bins=60,
            range=(-1000.0, 1000.0),
            alpha=0.70,
            color="#3274A1",
            label="Ion",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.hist(
            lateral_recoil,
            bins=60,
            range=(-1000.0, 1000.0),
            alpha=0.55,
            color="#E1812C",
            label="Recoil",
            edgecolor="black",
            linewidth=0.3,
        )
        ax.set_xlabel("Lateral Position (A)")
        ax.set_ylabel("Counts")
        ax.set_title("Lateral Range Distribution: Ion/Recoil")
        ax.legend()

    def plot_lateral_photons(ax):
        ax.hist(
            lateral_photons,
            bins=60,
            range=(-1200.0, 1200.0),
            alpha=0.82,
            color="#4C956C",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.set_xlabel("Lateral Position (A)")
        ax.set_ylabel("Photon Counts")
        ax.set_title("Lateral Range Distribution: Photons")

    def plot_lateral_ionization(ax):
        ax.plot(
            lateral_axis,
            lateral_ionization_profile,
            color="#8E6C8A",
            linewidth=1.6,
            label="Ionization",
        )
        ax.fill_between(
            lateral_axis,
            lateral_ionization_profile,
            color="#8E6C8A",
            alpha=0.25,
        )
        ax.set_xlabel("Lateral Position (A)")
        ax.set_ylabel("Ionization (arb. units)")
        ax.set_title("Lateral Range Distribution: Ionization")
        ax.legend()

    def plot_backscatter_energy(ax):
        ax.hist(
            backscatter_energy,
            bins=50,
            range=(0.0, 50000.0),
            alpha=0.85,
            color="#B85C38",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Counts")
        ax.set_title("Backscattered: Energy")

    def plot_backscatter_angle(ax):
        ax.hist(
            backscatter_angle,
            bins=45,
            range=(90.0, 180.0),
            alpha=0.85,
            color="#D97706",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.set_xlabel("Angle (deg)")
        ax.set_ylabel("Counts")
        ax.set_title("Backscattered: Angle")

    def plot_transmitted_energy(ax):
        ax.hist(
            transmitted_energy,
            bins=50,
            range=(0.0, 50000.0),
            alpha=0.85,
            color="#2E86AB",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Counts")
        ax.set_title("Transmitted: Energy")

    def plot_transmitted_angle(ax):
        ax.hist(
            transmitted_angle,
            bins=45,
            range=(0.0, 90.0),
            alpha=0.85,
            color="#3B82A0",
            edgecolor="black",
            linewidth=0.4,
        )
        ax.set_xlabel("Angle (deg)")
        ax.set_ylabel("Counts")
        ax.set_title("Transmitted: Angle")

    # ---- per-plot 4-value stats (shown in zoom popup) ---------------
    def stats_traj_start() -> List[tuple[str, str]]:
        start_r = np.hypot(start_x, start_y)
        return [
            ("Trajectories", f"{len(trajectories_xy)}"),
            ("Mean X(start)", f"{_safe_mean(start_x):.1f} A"),
            ("Mean Y(start)", f"{_safe_mean(start_y):.1f} A"),
            ("Max start radius", f"{_safe_max(start_r):.1f} A"),
        ]

    def stats_traj_end() -> List[tuple[str, str]]:
        end_r = np.hypot(end_x, end_y)
        return [
            ("Trajectories", f"{len(trajectories_xy)}"),
            ("Mean X(end)", f"{_safe_mean(end_x):.1f} A"),
            ("Mean Y(end)", f"{_safe_mean(end_y):.1f} A"),
            ("Mean end radius", f"{_safe_mean(end_r):.1f} A"),
        ]

    def stats_traj_collisions() -> List[tuple[str, str]]:
        collisions_arr = np.array(collision_counts, dtype=float)
        return [
            ("Trajectories", f"{len(trajectories_xy)}"),
            ("Collision points", f"{len(collision_x)}"),
            ("Mean collisions / ion", f"{_safe_mean(collisions_arr):.1f}"),
            ("Mean segment length", f"{_safe_mean(np.array(segment_lengths, dtype=float)):.1f} A"),
        ]

    def stats_range_ion_recoil() -> List[tuple[str, str]]:
        return [
            ("Ion samples", f"{len(ion_depths)}"),
            ("Recoil samples", f"{len(recoil_depths)}"),
            ("Mean ion depth", f"{_safe_mean(ion_depths):.1f} A"),
            ("Mean recoil depth", f"{_safe_mean(recoil_depths):.1f} A"),
        ]

    def stats_range_photons() -> List[tuple[str, str]]:
        peak = _safe_peak_center(photon_depths, bins=60, limits=(0.0, 4000.0))
        return [
            ("Photon samples", f"{len(photon_depths)}"),
            ("Mean depth", f"{_safe_mean(photon_depths):.1f} A"),
            ("Peak depth", f"{peak:.1f} A"),
            ("P90 depth", f"{_safe_percentile(photon_depths, 90.0):.1f} A"),
        ]

    def stats_range_ionization() -> List[tuple[str, str]]:
        peak_idx = int(np.argmax(ionization_range_profile))
        return [
            ("Depth bins", f"{len(ionization_depth_axis)}"),
            ("Peak depth", f"{ionization_depth_axis[peak_idx]:.1f} A"),
            ("Peak ionization", f"{_safe_max(ionization_range_profile):.2f}"),
            (
                "Integrated ionization",
                f"{float(np.trapz(ionization_range_profile, ionization_depth_axis)):.1f}",
            ),
        ]

    def stats_lateral_ion_recoil() -> List[tuple[str, str]]:
        return [
            ("Ion samples", f"{len(lateral_ion)}"),
            ("Recoil samples", f"{len(lateral_recoil)}"),
            ("Std ion lateral", f"{_safe_std(lateral_ion):.1f} A"),
            ("Std recoil lateral", f"{_safe_std(lateral_recoil):.1f} A"),
        ]

    def stats_lateral_photons() -> List[tuple[str, str]]:
        abs_lat = np.abs(lateral_photons)
        return [
            ("Photon samples", f"{len(lateral_photons)}"),
            ("Mean lateral", f"{_safe_mean(lateral_photons):.1f} A"),
            ("Std lateral", f"{_safe_std(lateral_photons):.1f} A"),
            ("P95 |lateral|", f"{_safe_percentile(abs_lat, 95.0):.1f} A"),
        ]

    def stats_lateral_ionization() -> List[tuple[str, str]]:
        peak_idx = int(np.argmax(lateral_ionization_profile))
        return [
            ("Lateral bins", f"{len(lateral_axis)}"),
            ("Peak lateral", f"{lateral_axis[peak_idx]:.1f} A"),
            ("Peak ionization", f"{_safe_max(lateral_ionization_profile):.2f}"),
            (
                "Integrated ionization",
                f"{float(np.trapz(lateral_ionization_profile, lateral_axis)):.1f}",
            ),
        ]

    def stats_backscatter_energy() -> List[tuple[str, str]]:
        peak = _safe_peak_center(backscatter_energy, bins=50, limits=(0.0, 50000.0))
        return [
            ("Backscattered ions", f"{len(backscatter_energy)}"),
            ("Mean energy", f"{_safe_mean(backscatter_energy):.1f} eV"),
            ("Peak energy", f"{peak:.1f} eV"),
            ("P90 energy", f"{_safe_percentile(backscatter_energy, 90.0):.1f} eV"),
        ]

    def stats_backscatter_angle() -> List[tuple[str, str]]:
        peak = _safe_peak_center(backscatter_angle, bins=45, limits=(90.0, 180.0))
        return [
            ("Backscattered ions", f"{len(backscatter_angle)}"),
            ("Mean angle", f"{_safe_mean(backscatter_angle):.1f} deg"),
            ("Std angle", f"{_safe_std(backscatter_angle):.1f} deg"),
            ("Peak angle", f"{peak:.1f} deg"),
        ]

    def stats_transmitted_energy() -> List[tuple[str, str]]:
        peak = _safe_peak_center(transmitted_energy, bins=50, limits=(0.0, 50000.0))
        return [
            ("Transmitted ions", f"{len(transmitted_energy)}"),
            ("Mean energy", f"{_safe_mean(transmitted_energy):.1f} eV"),
            ("Peak energy", f"{peak:.1f} eV"),
            ("P90 energy", f"{_safe_percentile(transmitted_energy, 90.0):.1f} eV"),
        ]

    def stats_transmitted_angle() -> List[tuple[str, str]]:
        peak = _safe_peak_center(transmitted_angle, bins=45, limits=(0.0, 90.0))
        return [
            ("Transmitted ions", f"{len(transmitted_angle)}"),
            ("Mean angle", f"{_safe_mean(transmitted_angle):.1f} deg"),
            ("Std angle", f"{_safe_std(transmitted_angle):.1f} deg"),
            ("P90 angle", f"{_safe_percentile(transmitted_angle, 90.0):.1f} deg"),
        ]

    plots = {
        "traj_start": {
            "name": "Trajectories: Start",
            "plot_func": plot_traj_start,
            "projection": None,
            "stats_func": stats_traj_start,
        },
        "traj_end": {
            "name": "Trajectories: End",
            "plot_func": plot_traj_end,
            "projection": None,
            "stats_func": stats_traj_end,
        },
        "traj_collisions": {
            "name": "Trajectories: Collisions",
            "plot_func": plot_traj_collisions,
            "projection": None,
            "stats_func": stats_traj_collisions,
        },
        "range_ion_recoil": {
            "name": "Range Distributions: Ion/Recoil",
            "plot_func": plot_range_ion_recoil,
            "projection": None,
            "stats_func": stats_range_ion_recoil,
        },
        "range_photons": {
            "name": "Range Distributions: Photons",
            "plot_func": plot_range_photons,
            "projection": None,
            "stats_func": stats_range_photons,
        },
        "range_ionization": {
            "name": "Range Distributions: Ionization",
            "plot_func": plot_range_ionization,
            "projection": None,
            "stats_func": stats_range_ionization,
        },
        "lateral_ion_recoil": {
            "name": "Lateral Range Distribution: Ion/Recoil",
            "plot_func": plot_lateral_ion_recoil,
            "projection": None,
            "stats_func": stats_lateral_ion_recoil,
        },
        "lateral_photons": {
            "name": "Lateral Range Distribution: Photons",
            "plot_func": plot_lateral_photons,
            "projection": None,
            "stats_func": stats_lateral_photons,
        },
        "lateral_ionization": {
            "name": "Lateral Range Distribution: Ionization",
            "plot_func": plot_lateral_ionization,
            "projection": None,
            "stats_func": stats_lateral_ionization,
        },
        "backscatter_energy": {
            "name": "Backscattered: Energy",
            "plot_func": plot_backscatter_energy,
            "projection": None,
            "stats_func": stats_backscatter_energy,
        },
        "backscatter_angle": {
            "name": "Backscattered: Angle",
            "plot_func": plot_backscatter_angle,
            "projection": None,
            "stats_func": stats_backscatter_angle,
        },
        "transmitted_energy": {
            "name": "Transmitted: Energy",
            "plot_func": plot_transmitted_energy,
            "projection": None,
            "stats_func": stats_transmitted_energy,
        },
        "transmitted_angle": {
            "name": "Transmitted: Angle",
            "plot_func": plot_transmitted_angle,
            "projection": None,
            "stats_func": stats_transmitted_angle,
        },
    }

    ion_mean = _safe_mean(ion_depths)
    ion_std = _safe_std(ion_depths)
    ion_centered = ion_depths - ion_mean
    ion_var = float(np.mean(ion_centered**2)) if ion_depths.size else 0.0
    if ion_var > 1e-12:
        ion_skew = float(np.mean(ion_centered**3) / (ion_var**1.5))
        ion_kurt = float(np.mean(ion_centered**4) / (ion_var**2) - 3.0)
    else:
        ion_skew = 0.0
        ion_kurt = 0.0

    numerical = [
        {"name": "Total Ions Simulated", "value": f"{n_total}"},
        {"name": "Trajectory Count", "value": f"{n_traj}"},
        {"name": "Mean Collisions / Ion", "value": f"{_safe_mean(np.array(collision_counts, dtype=float)):.1f}"},
        {"name": "Ions Stopped in Target", "value": f"{len(ion_depths)}"},
        {"name": "Backscattered Ions", "value": f"{len(backscatter_energy)}"},
        {"name": "Transmitted Ions", "value": f"{len(transmitted_energy)}"},
        {"name": "Mean Ion Range", "value": f"{ion_mean:.1f} A"},
        {"name": "Range Straggling", "value": f"{ion_std:.1f} A"},
        {"name": "Skewness", "value": f"{ion_skew:.3f}"},
        {"name": "Kurtosis", "value": f"{ion_kurt:.3f}"},
    ]

    return plots, numerical


# Generate once at module load
AVAILABLE_PLOTS, AVAILABLE_NUMERICAL_VALUES = _generate_hardcoded_data()


# =====================================================================
#  Zoomed-plot dialog  (shown on double-click)
# =====================================================================

class ZoomedPlotDialog(QDialog):
    """Full-size matplotlib dialog with NavigationToolbar and stats table."""

    def __init__(self, plot_id: str, plot_info: dict, parent=None, font_size: float = 10.0):
        super().__init__(parent)
        self._font_size = float(font_size)
        self.setWindowTitle(plot_info["name"])
        self.resize(1100, 680)

        root = QHBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(8)

        plot_widget = QWidget(self)
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(0)

        projection = plot_info.get("projection")
        self.figure = Figure(figsize=(8, 5), dpi=100)
        if projection:
            self.ax = self.figure.add_subplot(111, projection=projection)
        else:
            self.ax = self.figure.add_subplot(111)

        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        root.addWidget(plot_widget, 1)

        stats_group = QGroupBox("Numerical Values")
        stats_group.setMinimumWidth(320)
        stats_layout = QVBoxLayout(stats_group)
        stats_layout.setContentsMargins(6, 6, 6, 6)

        self._stats_table = QTableWidget(0, 2)
        self._stats_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self._stats_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self._stats_table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.Stretch
        )
        self._stats_table.verticalHeader().setVisible(False)
        self._stats_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._stats_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self._stats_table.setAlternatingRowColors(True)
        stats_layout.addWidget(self._stats_table)
        root.addWidget(stats_group, 0)

        plot_info["plot_func"](self.ax)
        _apply_axes_font_size(self.ax, self._font_size)
        self._populate_stats_table(self._resolve_stats_rows(plot_info))

        # 3D: remap rotation to right mouse button
        if projection:
            try:
                self.ax._rotate_btn = np.array([3])
                self.ax._zoom_btn = np.array([2])
            except AttributeError:
                pass

        self.figure.tight_layout()
        self.canvas.draw()

    def _resolve_stats_rows(self, plot_info: dict) -> List[tuple[str, str]]:
        rows: List[tuple[str, str]] = []
        stats_func = plot_info.get("stats_func")
        if callable(stats_func):
            try:
                rows = list(stats_func())
            except Exception:
                rows = []

        if not rows:
            rows = [
                (str(entry.get("name", "Value")), str(entry.get("value", "-")))
                for entry in AVAILABLE_NUMERICAL_VALUES
            ]
        return rows

    def _populate_stats_table(self, rows: List[tuple[str, str]]) -> None:
        self._stats_table.setRowCount(len(rows))
        for row_idx, (name, value) in enumerate(rows):
            name_item = QTableWidgetItem(name)
            value_item = QTableWidgetItem(value)
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            value_item.setFlags(value_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._stats_table.setItem(row_idx, 0, name_item)
            self._stats_table.setItem(row_idx, 1, value_item)


# =====================================================================
#  Tile style constants
# =====================================================================

_STYLE_TILE_NOBORDER = "border: none;"

_STYLE_TILE = "border: 1px solid #bbb; border-radius: 4px;"
_STYLE_DROP_BEFORE = (
    "border: 1px solid #bbb; border-radius: 4px; border-left: 4px solid #3274A1;"
)
_STYLE_DROP_AFTER = (
    "border: 1px solid #bbb; border-radius: 4px; border-right: 4px solid #3274A1;"
)


def _apply_axes_font_size(ax, font_size: float) -> None:
    """Apply a consistent base font size to titles, labels, ticks and legends."""
    base = max(6.0, float(font_size))
    # Keep the title at the same size as the axis labels (a noticeably larger
    # title was reported as distracting on publication plots).
    ax.title.set_fontsize(base)
    ax.xaxis.label.set_fontsize(base)
    ax.yaxis.label.set_fontsize(base)
    ax.tick_params(axis="both", labelsize=max(6.0, base - 1.0))

    if hasattr(ax, "zaxis"):
        try:
            ax.zaxis.label.set_fontsize(base)
            ax.tick_params(axis="z", labelsize=max(6.0, base - 1.0))
        except Exception:
            pass

    legend = ax.get_legend()
    if legend is not None:
        for text in legend.get_texts():
            text.set_fontsize(max(6.0, base - 1.0))
        title = legend.get_title()
        if title is not None:
            title.set_fontsize(base)


def _make_drag_cursor_pixmap() -> QPixmap:
    """Return a small high-contrast pixmap used as active drag cursor."""
    size = 24
    pix = QPixmap(size, size)
    pix.fill(Qt.GlobalColor.transparent)

    p = QPainter(pix)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor("#2c7be5"))
    p.drawEllipse(2, 2, size - 4, size - 4)

    p.setPen(QPen(QColor("white"), 2))
    p.drawLine(size // 2, 6, size // 2, size - 6)
    p.drawLine(6, size // 2, size - 6, size // 2)
    p.end()
    return pix


def _apply_drag_cursors(drag: QDrag, cursor_pixmap: QPixmap) -> None:
    """Set the same cursor pixmap for all drag actions."""
    drag.setDragCursor(cursor_pixmap, Qt.DropAction.MoveAction)
    drag.setDragCursor(cursor_pixmap, Qt.DropAction.CopyAction)
    drag.setDragCursor(cursor_pixmap, Qt.DropAction.LinkAction)
    drag.setDragCursor(cursor_pixmap, Qt.DropAction.IgnoreAction)
    drag.setDragCursor(cursor_pixmap, Qt.DropAction.TargetMoveAction)


class ReorderListWidget(QListWidget):
    """QListWidget with reliable reorder commit + drag cursor feedback."""

    reordered = pyqtSignal()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.MouseButton.LeftButton:
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        super().mouseReleaseEvent(event)

    def startDrag(self, supportedActions):
        indexes = self.selectedIndexes()
        if not indexes:
            return

        mime_data = self.model().mimeData(indexes)
        if mime_data is None:
            return

        drag = QDrag(self)
        drag.setMimeData(mime_data)

        rect = self.visualRect(indexes[0])
        if rect.isValid():
            row_pixmap = self.viewport().grab(rect)
            if not row_pixmap.isNull():
                drag.setPixmap(row_pixmap)
                drag.setHotSpot(QPoint(10, max(8, row_pixmap.height() // 2)))

        drag_cursor = _make_drag_cursor_pixmap()
        _apply_drag_cursors(drag, drag_cursor)
        QApplication.setOverrideCursor(QCursor(drag_cursor, 4, 4))
        try:
            drag.exec(supportedActions, Qt.DropAction.MoveAction)
        finally:
            QApplication.restoreOverrideCursor()
            self.setCursor(Qt.CursorShape.OpenHandCursor)
            QTimer.singleShot(0, self.reordered.emit)

    def dragEnterEvent(self, event):
        self.setCursor(Qt.CursorShape.ClosedHandCursor)
        super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        self.setCursor(Qt.CursorShape.ClosedHandCursor)
        super().dragMoveEvent(event)

    def dragLeaveEvent(self, event):
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        super().dragLeaveEvent(event)

    def dropEvent(self, event):
        super().dropEvent(event)
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        QTimer.singleShot(0, self.reordered.emit)


# =====================================================================
#  Single plot tile  (canvas + toolbar, full-canvas drag, drop target)
# =====================================================================

class PlotTile(QFrame):
    """A framed widget holding one matplotlib plot with its toolbar.

    Features:
    - Left-click + drag on the canvas to reorder (when no mpl tool active).
    - Accepts drops from other tiles; shows left/right drop indicator.
    - Double-click on the canvas opens a zoomed dialog.
    - 3D plots use right mouse button for rotation.
    """

    zoom_requested = pyqtSignal(str)  # emits plot_id
    tile_reorder_requested = pyqtSignal(str, str, bool)  # source, target, before

    def __init__(self, plot_id: str, plot_info: dict, parent=None, font_size: float = 10.0,
                 bin_factor: int = 1):
        super().__init__(parent)
        self.plot_id = plot_id
        self.plot_info = plot_info
        self._is_3d = plot_info.get("projection") == "3d"
        self._has_border = True
        self._font_size = float(font_size)
        self._bin_factor = max(1, int(bin_factor))
        self._beam_from_top = True
        self.setAcceptDrops(True)

        self.setFrameStyle(QFrame.Shape.Box | QFrame.Shadow.Plain)
        self.setLineWidth(1)
        self.setStyleSheet(_STYLE_TILE)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(0)

        # --- matplotlib canvas ---
        projection = plot_info.get("projection")
        self.figure = Figure(tight_layout=True)
        if projection:
            self.ax = self.figure.add_subplot(111, projection=projection)
        else:
            self.ax = self.figure.add_subplot(111)

        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self.canvas.setCursor(Qt.CursorShape.OpenHandCursor)

        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setStyleSheet("QToolBar { border: none; spacing: 2px; }")

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        # Draw the plot
        self._call_plot_func(self.ax)
        _apply_axes_font_size(self.ax, self._font_size)
        self.figure.tight_layout()
        self.canvas.draw()

        # 3D: remap rotation to right mouse button, zoom to middle
        if self._is_3d:
            try:
                self.ax._rotate_btn = np.array([3])
                self.ax._zoom_btn = np.array([2])
            except AttributeError:
                pass

        # Double-click on canvas opens zoomed dialog
        self.canvas.mpl_connect("button_press_event", self._on_mpl_press)

        # Drag initiation via canvas event filter
        self._drag_start_pos: Optional[QPoint] = None
        self.canvas.installEventFilter(self)

    # ---- drag initiation via canvas ---------------------------------

    def eventFilter(self, obj, event):
        """Intercept left-click + drag on the canvas to start tile reorder.

        Drag is only initiated when no matplotlib toolbar tool (zoom/pan)
        is active.  The initial press is still forwarded to matplotlib so
        that double-click detection keeps working.
        """
        if obj is not self.canvas:
            return super().eventFilter(obj, event)

        etype = event.type()

        if (
            etype == event.Type.MouseButtonPress
            and event.button() == Qt.MouseButton.LeftButton
            and not self.toolbar.mode
        ):
            self._drag_start_pos = event.position().toPoint()
            self.canvas.setCursor(Qt.CursorShape.ClosedHandCursor)
            return False  # let matplotlib see the press (for dblclick)

        if etype == event.Type.MouseMove and self._drag_start_pos is not None:
            if event.buttons() & Qt.MouseButton.LeftButton:
                self.canvas.setCursor(Qt.CursorShape.ClosedHandCursor)
                dist = (
                    event.position().toPoint() - self._drag_start_pos
                ).manhattanLength()
                if dist >= QApplication.startDragDistance():
                    self._drag_start_pos = None
                    self._start_drag()
                    return True
            else:
                self._drag_start_pos = None

        if etype == event.Type.MouseButtonRelease:
            self._drag_start_pos = None
            self.canvas.setCursor(Qt.CursorShape.OpenHandCursor)

        return False

    def _start_drag(self):
        drag = QDrag(self)
        mime = QMimeData()
        mime.setText(self.plot_id)
        drag.setMimeData(mime)
        pixmap = self.grab().scaled(
            180, 130,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        drag.setPixmap(pixmap)
        drag.setHotSpot(QPoint(pixmap.width() // 2, pixmap.height() // 2))
        drag_cursor = _make_drag_cursor_pixmap()
        _apply_drag_cursors(drag, drag_cursor)
        QApplication.setOverrideCursor(QCursor(drag_cursor, 4, 4))
        drag.exec(Qt.DropAction.MoveAction)
        QApplication.restoreOverrideCursor()
        self.canvas.setCursor(Qt.CursorShape.OpenHandCursor)

    # ---- drop target -----------------------------------------------

    def dragEnterEvent(self, event):
        if (
            event.mimeData().hasText()
            and event.mimeData().text() in AVAILABLE_PLOTS
            and event.mimeData().text() != self.plot_id
        ):
            event.acceptProposedAction()
            self._update_drop_style(event.position().toPoint())
            self.setCursor(Qt.CursorShape.DragMoveCursor)

    def dragMoveEvent(self, event):
        self._update_drop_style(event.position().toPoint())
        event.acceptProposedAction()

    def dragLeaveEvent(self, event):
        self.setStyleSheet(
            _STYLE_TILE if self._has_border else _STYLE_TILE_NOBORDER
        )
        self.unsetCursor()

    def dropEvent(self, event):
        source_id = event.mimeData().text()
        before = event.position().toPoint().x() < self.width() / 2
        self.setStyleSheet(
            _STYLE_TILE if self._has_border else _STYLE_TILE_NOBORDER
        )
        self.unsetCursor()
        self.tile_reorder_requested.emit(source_id, self.plot_id, before)
        event.acceptProposedAction()

    def _update_drop_style(self, pos):
        if pos.x() < self.width() / 2:
            self.setStyleSheet(_STYLE_DROP_BEFORE)
        else:
            self.setStyleSheet(_STYLE_DROP_AFTER)

    # ---- display settings -------------------------------------------

    def set_border_visible(self, visible: bool):
        """Show or hide the tile border frame."""
        self._has_border = visible
        if visible:
            self.setFrameStyle(QFrame.Shape.Box | QFrame.Shadow.Plain)
            self.setStyleSheet(_STYLE_TILE)
        else:
            self.setFrameStyle(QFrame.Shape.NoFrame)
            self.setStyleSheet(_STYLE_TILE_NOBORDER)

    # ---- other -----------------------------------------------------

    def _on_mpl_press(self, event):
        if event.dblclick:
            self.zoom_requested.emit(self.plot_id)

    def redraw(self):
        """Redraw the plot (e.g. after data update)."""
        # The figure can hold extra artists (colorbars, twin axes) added by a
        # previous plot_func call — clear at the figure level so they vanish too.
        self.figure.clf()
        self.ax = self.figure.add_subplot(
            111, projection=self.plot_info.get("projection")
        )
        self._call_plot_func(self.ax)
        _apply_axes_font_size(self.ax, self._font_size)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def set_font_size(self, font_size: float):
        self._font_size = float(font_size)
        self.redraw()

    def set_bin_factor(self, factor: int) -> None:
        new = max(1, int(factor))
        if new == self._bin_factor:
            return
        self._bin_factor = new
        self.redraw()

    def set_beam_from_top(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if self._beam_from_top == enabled:
            return
        self._beam_from_top = enabled
        self.redraw()

    def _call_plot_func(self, ax) -> None:
        """Invoke the plot_func, passing bin_factor when the function accepts it."""
        fn = self.plot_info["plot_func"]
        try:
            fn(ax, bin_factor=self._bin_factor, beam_from_top=self._beam_from_top)
        except TypeError:
            try:
                fn(ax, bin_factor=self._bin_factor)
            except TypeError:
                fn(ax)


# =====================================================================
#  Plot area  (right side — nested splitters, drag-and-drop)
# =====================================================================

class PlotArea(QWidget):
    """Right-side area with resizable plot tiles arranged via nested
    QSplitters (vertical outer, horizontal per row).

    Tiles can be **resized** by dragging the splitter handles and
    **reordered** via drag-and-drop on the plot canvases.
    """

    order_changed = pyqtSignal(list)  # emits new plot_id order after D&D
    plot_double_clicked = pyqtSignal(str, object)  # plot_id, plot_info dict

    def __init__(self, parent=None):
        super().__init__(parent)

        self._plot_order: List[str] = []
        self._tiles: Dict[str, PlotTile] = {}

        # Display settings
        self._col_override: int = 0   # 0 = auto
        self._toolbar_visible: bool = True
        self._borders_visible: bool = True
        self._font_size: float = 10.0
        self._bin_factor: int = 1     # 1 = no combining; up to 6
        self._per_plot_bin_factor: Dict[str, int] = {}
        self._beam_from_top: bool = True

        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setSpacing(0)

        self._main_splitter = QSplitter(Qt.Orientation.Vertical)
        self._layout.addWidget(self._main_splitter)

        self._row_splitters: List[QSplitter] = []

    # --- public API ---------------------------------------------------

    def set_visible_plots(self, plot_ids: List[str]):
        """Rebuild the splitter grid to show exactly the requested plots."""
        self._plot_order = list(plot_ids)
        self._rebuild()

    def reorder(self, new_order: List[str]):
        """Reorder tiles to match *new_order* (list of plot_ids)."""
        self._plot_order = [pid for pid in new_order if pid in self._plot_order]
        self._rebuild()

    def set_columns(self, n: int):
        """Set column count.  0 = auto, 1/2/3 = fixed."""
        self._col_override = n
        self._rebuild()

    def set_toolbar_visible(self, visible: bool):
        """Show or hide matplotlib navigation toolbars on every tile."""
        self._toolbar_visible = visible
        for tile in self._tiles.values():
            tile.toolbar.setVisible(visible)

    def set_borders_visible(self, visible: bool):
        """Show or hide plot tile borders."""
        self._borders_visible = visible
        for tile in self._tiles.values():
            tile.set_border_visible(visible)

    def set_font_size(self, size: float):
        """Set the base plot font size for all tiles and zoom dialogs."""
        self._font_size = float(size)
        for tile in self._tiles.values():
            tile.set_font_size(self._font_size)

    def set_bin_factor(self, factor: int) -> None:
        """Combine *factor* adjacent histogram bins on every tile (1–6)."""
        new = max(1, min(6, int(factor)))
        self._bin_factor = new
        for tile in self._tiles.values():
            tile.set_bin_factor(new)

    def set_plot_bin_factor(self, plot_id: str, factor: int) -> None:
        """Set per-plot bin factor for one tile (1–6)."""
        new = max(1, min(6, int(factor)))
        self._per_plot_bin_factor[str(plot_id)] = new
        tile = self._tiles.get(str(plot_id))
        if tile is not None:
            tile.set_bin_factor(new)

    def set_beam_from_top(self, enabled: bool) -> None:
        self._beam_from_top = bool(enabled)
        for tile in self._tiles.values():
            tile.set_beam_from_top(self._beam_from_top)

    # --- internal rebuild ---------------------------------------------

    def _rebuild(self):
        # Destroy old tiles
        for tile in self._tiles.values():
            tile.setParent(None)
            tile.deleteLater()
        self._tiles.clear()

        # Remove old row splitters
        while self._main_splitter.count():
            w = self._main_splitter.widget(0)
            w.setParent(None)
        for sp in self._row_splitters:
            sp.deleteLater()
        self._row_splitters.clear()

        if not self._plot_order:
            return

        cols = self._get_columns(len(self._plot_order))
        for start in range(0, len(self._plot_order), cols):
            chunk = self._plot_order[start : start + cols]
            row_sp = QSplitter(Qt.Orientation.Horizontal)
            for pid in chunk:
                info = AVAILABLE_PLOTS.get(pid)
                if info is None:
                    continue
                tile = PlotTile(pid, info, font_size=self._font_size,
                                bin_factor=self._per_plot_bin_factor.get(pid, self._bin_factor))
                tile.zoom_requested.connect(self._open_zoomed)
                tile.tile_reorder_requested.connect(self._handle_tile_reorder)
                tile.setMinimumSize(250, 200)
                # Apply current display settings
                tile.toolbar.setVisible(self._toolbar_visible)
                tile.set_border_visible(self._borders_visible)
                tile.set_beam_from_top(self._beam_from_top)
                self._tiles[pid] = tile
                row_sp.addWidget(tile)
            self._row_splitters.append(row_sp)
            self._main_splitter.addWidget(row_sp)

    # --- drag-and-drop reorder ----------------------------------------

    def _handle_tile_reorder(self, source_id: str, target_id: str, before: bool):
        if source_id not in self._plot_order or target_id not in self._plot_order:
            return
        if source_id == target_id:
            return

        order = list(self._plot_order)
        order.remove(source_id)
        idx = order.index(target_id)
        if not before:
            idx += 1
        order.insert(idx, source_id)

        self._plot_order = order
        self._rebuild()
        self.order_changed.emit(list(self._plot_order))

    # --- helpers ------------------------------------------------------

    def _get_columns(self, n: int) -> int:
        if self._col_override > 0:
            return self._col_override
        return self._auto_columns(n)

    @staticmethod
    def _auto_columns(n: int) -> int:
        if n <= 1:
            return 1
        if n <= 4:
            return 2
        return 3

    def _open_zoomed(self, plot_id: str):
        info = AVAILABLE_PLOTS.get(plot_id)
        if info is None:
            return
        self.plot_double_clicked.emit(plot_id, info)


def _to_float(value: Any) -> Optional[float]:
    if isinstance(value, (int, float, np.number)):
        v = float(value)
        return v if np.isfinite(v) else None
    if isinstance(value, str):
        txt = value.strip()
        if not txt:
            return None
        txt = txt.replace(",", ".")
        try:
            v = float(txt)
        except ValueError:
            return None
        return v if np.isfinite(v) else None
    return None


def _to_numeric_vector(values: Any) -> Optional[np.ndarray]:
    if isinstance(values, np.ndarray):
        try:
            arr = np.asarray(values, dtype=float).reshape(-1)
        except (TypeError, ValueError):
            return None
    elif isinstance(values, (list, tuple)):
        if not values:
            return None
        try:
            arr = np.asarray(values, dtype=float).reshape(-1)
        except (TypeError, ValueError):
            return None
    else:
        return None

    if arr.size < 2:
        return None
    if not np.all(np.isfinite(arr)):
        return None
    return arr


def _default_preview_name(dataset_name: str, data_path: str) -> str:
    base = os.path.basename(data_path) if data_path else "data"
    return f"{dataset_name} - {base}"


def _extract_json_datasets(node: Any, prefix: str = "") -> List[Dict[str, Any]]:
    datasets: List[Dict[str, Any]] = []

    if isinstance(node, dict):
        x_direct = _to_numeric_vector(node.get("x"))
        y_direct = _to_numeric_vector(node.get("y"))
        if x_direct is not None and y_direct is not None and x_direct.size == y_direct.size:
            label = str(node.get("name") or node.get("label") or prefix or "dataset")
            datasets.append(
                {
                    "name": label,
                    "x": x_direct,
                    "y": y_direct,
                    "x_label": "x",
                    "y_label": "y",
                }
            )

        numeric_series: Dict[str, np.ndarray] = {}
        for key, value in node.items():
            if key in {"x", "y", "name", "label"}:
                continue
            arr = _to_numeric_vector(value)
            if arr is not None:
                numeric_series[str(key)] = arr

        preferred_x = next(
            (
                key
                for key in ("energy", "depth", "time", "position", "range", "z", "index")
                if key in numeric_series
            ),
            None,
        )

        if preferred_x is not None:
            x_ref = numeric_series[preferred_x]
            for key, arr in numeric_series.items():
                if key == preferred_x or arr.size != x_ref.size:
                    continue
                label = f"{prefix}.{key}" if prefix else key
                datasets.append(
                    {
                        "name": label,
                        "x": x_ref,
                        "y": arr,
                        "x_label": preferred_x,
                        "y_label": key,
                    }
                )
        else:
            for key, arr in numeric_series.items():
                label = f"{prefix}.{key}" if prefix else key
                datasets.append(
                    {
                        "name": label,
                        "x": np.arange(arr.size, dtype=float),
                        "y": arr,
                        "x_label": "Index",
                        "y_label": key,
                    }
                )

        processed_vector_keys = set(numeric_series.keys()) | {"x", "y"}
        for key, value in node.items():
            if key in processed_vector_keys:
                continue
            if isinstance(value, (dict, list)):
                child_prefix = f"{prefix}.{key}" if prefix else str(key)
                datasets.extend(_extract_json_datasets(value, child_prefix))

    elif isinstance(node, list):
        arr = _to_numeric_vector(node)
        if arr is not None:
            datasets.append(
                {
                    "name": prefix or "dataset",
                    "x": np.arange(arr.size, dtype=float),
                    "y": arr,
                    "x_label": "Index",
                    "y_label": prefix or "Value",
                }
            )
        else:
            for idx, value in enumerate(node):
                if isinstance(value, (dict, list)):
                    child_prefix = f"{prefix}[{idx}]" if prefix else f"dataset_{idx + 1}"
                    datasets.extend(_extract_json_datasets(value, child_prefix))

    return datasets


def _extract_table_rows(path: str) -> List[List[str]]:
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        lines = [line.rstrip("\n") for line in fh]

    if not lines:
        return []

    suffix = os.path.splitext(path)[1].lower()
    rows: List[List[str]] = []
    if suffix == ".txt":
        for line in lines:
            if not line.strip() or line.strip().startswith("#"):
                continue
            rows.append(line.split())
        return rows

    sample = "\n".join(lines[:20])
    delimiter = ","
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        delimiter = dialect.delimiter
    except csv.Error:
        if ";" in sample:
            delimiter = ";"
        elif "\t" in sample:
            delimiter = "\t"

    for line in lines:
        if not line.strip() or line.strip().startswith("#"):
            continue
        rows.append([part.strip() for part in line.split(delimiter)])
    return rows


def _extract_table_datasets(path: str) -> List[Dict[str, Any]]:
    rows = _extract_table_rows(path)
    if not rows:
        return []

    header_row = rows[0]
    has_header = any(_to_float(token) is None for token in header_row)
    data_rows = rows[1:] if has_header else rows
    if not data_rows:
        return []

    n_cols = max(len(row) for row in data_rows)
    if n_cols == 0:
        return []

    col_names: List[str] = []
    for idx in range(n_cols):
        if has_header and idx < len(header_row) and header_row[idx].strip():
            col_names.append(header_row[idx].strip())
        else:
            col_names.append(f"Column {idx + 1}")

    matrix = np.full((len(data_rows), n_cols), np.nan, dtype=float)
    for r_idx, row in enumerate(data_rows):
        for c_idx, token in enumerate(row[:n_cols]):
            v = _to_float(token)
            if v is not None:
                matrix[r_idx, c_idx] = v

    datasets: List[Dict[str, Any]] = []
    if n_cols == 1:
        y_vals = matrix[:, 0]
        mask = np.isfinite(y_vals)
        if int(np.count_nonzero(mask)) >= 2:
            y = y_vals[mask]
            datasets.append(
                {
                    "name": col_names[0],
                    "x": np.arange(y.size, dtype=float),
                    "y": y,
                    "x_label": "Index",
                    "y_label": col_names[0],
                }
            )
        return datasets

    x_raw = matrix[:, 0]
    x_label = col_names[0]
    for c_idx in range(1, n_cols):
        y_raw = matrix[:, c_idx]
        mask = np.isfinite(x_raw) & np.isfinite(y_raw)
        if int(np.count_nonzero(mask)) < 2:
            continue
        datasets.append(
            {
                "name": col_names[c_idx],
                "x": x_raw[mask],
                "y": y_raw[mask],
                "x_label": x_label,
                "y_label": col_names[c_idx],
            }
        )

    return datasets


def load_datasets_from_path(path: str) -> List[Dict[str, Any]]:
    suffix = os.path.splitext(path)[1].lower()
    raw: List[Dict[str, Any]] = []

    if suffix == ".json":
        with open(path, "r", encoding="utf-8") as fh:
            payload = json.load(fh)
        raw = _extract_json_datasets(payload)
    elif suffix in {".csv", ".txt"}:
        raw = _extract_table_datasets(path)
    else:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
            raw = _extract_json_datasets(payload)
        except Exception:
            raw = _extract_table_datasets(path)

    normalized: List[Dict[str, Any]] = []
    for idx, entry in enumerate(raw):
        x_arr = _to_numeric_vector(entry.get("x"))
        y_arr = _to_numeric_vector(entry.get("y"))
        if x_arr is None or y_arr is None or x_arr.size != y_arr.size:
            continue

        name = str(entry.get("name") or f"Dataset {idx + 1}")
        normalized.append(
            {
                "id": f"dataset_{idx + 1}",
                "name": name,
                "x": x_arr,
                "y": y_arr,
                "x_label": str(entry.get("x_label") or "X"),
                "y_label": str(entry.get("y_label") or "Y"),
                "path": path,
            }
        )

    return normalized


class SinglePlotArea(QWidget):
    """Single-plot renderer for multiple user-selected datasets."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._curves: List[Dict[str, Any]] = []
        self._toolbar_visible: bool = True
        self._borders_visible: bool = True
        self._font_size: float = 10.0

        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        self._plot_frame = QFrame(self)
        self._plot_frame.setFrameStyle(QFrame.Shape.Box | QFrame.Shadow.Plain)
        self._plot_frame.setLineWidth(1)
        self._plot_frame.setStyleSheet(_STYLE_TILE)

        frame_layout = QVBoxLayout(self._plot_frame)
        frame_layout.setContentsMargins(2, 2, 2, 2)
        frame_layout.setSpacing(0)

        self.figure = Figure(tight_layout=True)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setStyleSheet("QToolBar { border: none; spacing: 2px; }")

        frame_layout.addWidget(self.toolbar)
        frame_layout.addWidget(self.canvas)
        root.addWidget(self._plot_frame, 1)

        self._render_plot()

    def set_curves(self, curves: List[Dict[str, Any]]) -> None:
        self._curves = list(curves)
        self._render_plot()

    # Backward-compatible convenience for legacy single-dataset calls.
    def set_dataset(self, dataset: Optional[Dict[str, Any]], preview_name: str = "") -> None:
        if dataset is None:
            self.set_curves([])
            return
        preview = preview_name.strip() or _default_preview_name(
            str(dataset.get("name") or "dataset"),
            str(dataset.get("path") or ""),
        )
        self.set_curves(
            [
                {
                    "preview_name": preview,
                    "dataset_name": str(dataset.get("name") or "dataset"),
                    "path": str(dataset.get("path") or ""),
                    "dataset": dataset,
                }
            ]
        )

    def set_toolbar_visible(self, visible: bool) -> None:
        self._toolbar_visible = bool(visible)
        self.toolbar.setVisible(self._toolbar_visible)

    def set_borders_visible(self, visible: bool) -> None:
        self._borders_visible = bool(visible)
        if self._borders_visible:
            self._plot_frame.setFrameStyle(QFrame.Shape.Box | QFrame.Shadow.Plain)
            self._plot_frame.setStyleSheet(_STYLE_TILE)
        else:
            self._plot_frame.setFrameStyle(QFrame.Shape.NoFrame)
            self._plot_frame.setStyleSheet(_STYLE_TILE_NOBORDER)

    def set_font_size(self, size: float) -> None:
        self._font_size = float(size)
        self._render_plot()

    def _render_plot(self) -> None:
        self.ax.clear()
        if not self._curves:
            self.ax.text(
                0.5,
                0.5,
                "No curves selected",
                horizontalalignment="center",
                verticalalignment="center",
                transform=self.ax.transAxes,
            )
            self.ax.set_axis_off()
        else:
            self.ax.set_axis_on()
            x_labels: List[str] = []
            y_labels: List[str] = []
            plotted = 0
            for curve in self._curves:
                dataset = curve.get("dataset") or {}
                x_vals = np.asarray(dataset.get("x", []), dtype=float)
                y_vals = np.asarray(dataset.get("y", []), dtype=float)
                if x_vals.size < 2 or y_vals.size < 2 or x_vals.size != y_vals.size:
                    continue

                preview = str(curve.get("preview_name") or dataset.get("name") or "curve")
                self.ax.plot(x_vals, y_vals, marker="o", linewidth=1.6, label=preview)
                x_labels.append(str(dataset.get("x_label") or "X"))
                y_labels.append(str(dataset.get("y_label") or "Y"))
                plotted += 1

            if plotted == 0:
                self.ax.text(
                    0.5,
                    0.5,
                    "No valid numeric curve data",
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=self.ax.transAxes,
                )
                self.ax.set_axis_off()
            else:
                x_label = x_labels[0] if x_labels and all(lbl == x_labels[0] for lbl in x_labels) else "X"
                y_label = y_labels[0] if y_labels and all(lbl == y_labels[0] for lbl in y_labels) else "Y"
                self.ax.set_xlabel(x_label)
                self.ax.set_ylabel(y_label)
                self.ax.set_title("Single Plot (Multiple Curves)")
                self.ax.legend()
            self.ax.grid(alpha=0.25)

        _apply_axes_font_size(self.ax, self._font_size)
        self.figure.tight_layout()
        self.canvas.draw_idle()


# =====================================================================
#  Single-Plot Page  (top-level tab, populated on tile double-click)
# =====================================================================

_COLORS = [
    "#3274A1", "#E1812C", "#3A923A", "#C03D3E", "#9372B2",
    "#8E6C8A", "#D97706", "#2E86AB", "#4C956C", "#B85C38",
]

_DRAWSTYLES = [
    ("Line",         "default"),
    ("Step (mid)",   "steps-mid"),
    ("Step (pre)",   "steps-pre"),
    ("Step (post)",  "steps-post"),
]

_LINESTYLES = [
    ("Solid",        "-"),
    ("Dashed",       "--"),
    ("Dotted",       ":"),
    ("Dash-dot",     "-."),
    ("None",         "None"),
]

_MARKERS = [
    ("None",         "None"),
    ("Circle",       "o"),
    ("Square",       "s"),
    ("Triangle ▲",   "^"),
    ("Triangle ▼",   "v"),
    ("Diamond",      "D"),
    ("Cross",        "x"),
    ("Plus",         "+"),
    ("Point",        "."),
    ("Star",         "*"),
]

_LEGEND_POSITIONS = [
    "best", "upper right", "upper left", "lower left", "lower right",
    "right", "center left", "center right", "lower center", "upper center", "center",
]


def _color_icon(color: str, size: int = 14) -> QIcon:
    """Build a simple solid-color square icon for the curves list."""
    pix = QPixmap(size, size)
    pix.fill(Qt.GlobalColor.transparent)
    p = QPainter(pix)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
    qc = QColor(color)
    if not qc.isValid():
        qc = QColor("#888888")
    p.setBrush(qc)
    p.setPen(QPen(QColor("#444"), 1))
    p.drawRect(0, 0, size - 1, size - 1)
    p.end()
    return QIcon(pix)


class CollapsibleBox(QGroupBox):
    """A bordered box whose body collapses/expands via a clickable header.

    The title and an optional hint button live in the header row. ``setEnabled``
    is overridden to grey out only the body, so the header arrow stays usable
    even when the section's controls are disabled.
    """

    def __init__(self, title: str, *, hint_btn: Optional[QWidget] = None,
                 expanded: bool = False, parent: Optional[QWidget] = None):
        super().__init__("", parent)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        self._toggle = QToolButton(self)
        self._toggle.setText(title)
        self._toggle.setCheckable(True)
        self._toggle.setChecked(expanded)
        self._toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._toggle.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self._toggle.setCursor(Qt.CursorShape.PointingHandCursor)
        self._toggle.setStyleSheet("QToolButton { border: none; font-weight: 600; }")
        self._toggle.toggled.connect(self._on_toggled)

        header = QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.setSpacing(6)
        header.addWidget(self._toggle)
        if hint_btn is not None:
            header.addWidget(hint_btn)
        header.addStretch(1)
        outer.addLayout(header)

        self._content = QWidget(self)
        self._content.setVisible(expanded)
        outer.addWidget(self._content)

    @property
    def content(self) -> QWidget:
        return self._content

    def set_expanded(self, expanded: bool) -> None:
        self._toggle.setChecked(expanded)

    def _on_toggled(self, checked: bool) -> None:
        self._toggle.setArrowType(
            Qt.ArrowType.DownArrow if checked else Qt.ArrowType.RightArrow
        )
        self._content.setVisible(checked)

    def setEnabled(self, enabled: bool) -> None:  # type: ignore[override]
        # Only the body greys out; the header toggle stays clickable.
        self._content.setEnabled(enabled)


class SinglePlotPage(QWidget):
    """Multi-curve overlay plot with full per-curve matplotlib styling.

    Workflow:
      - Double-click a tile in any MC Results tab → all visible series of
        that tile are *appended* as curves here (existing curves stay).
      - Each curve has independent visibility, color, line/draw style,
        line width, marker, alpha and (for depth profiles) Gaussian
        convolution.
      - Axes (title, labels, log scale, grid, legend) are global and
        editable.
    """

    def _init_hints(self) -> None:
        """Initialize hint system and lock it to this page."""
        try:
            hints_path = Path(__file__).resolve().parents[2] / "widgets" / "hints.json"
            self._hint_system = HintSystem(repo_path=hints_path, parent=self)
            self._hint_system.set_current_page("Single Plot")
        except Exception:
            self._hint_system = None

    def _hint_btn(self, hint_id: str, parent: Optional[QWidget] = None) -> QPushButton:
        """Create a '?' button for this page; disabled if hint unavailable."""
        if self._hint_system is None:
            btn = QPushButton("?", parent)
            btn.setEnabled(False)
            btn.setFixedSize(22, 22)
            btn.setToolTip("Hints not available")
            return btn
        try:
            return self._hint_system.make_hint_button(
                page_id="Single Plot", hint_id=hint_id, parent=parent
            )
        except Exception:
            btn = QPushButton("?", parent)
            btn.setEnabled(False)
            btn.setFixedSize(22, 22)
            btn.setToolTip(f"Hint '{hint_id}' not available")
            return btn

    def _groupbox_header(self, title: str, hint_id: Optional[str] = None,
                         parent: Optional[QWidget] = None) -> QWidget:
        """Header row: bold title label + optional '?' hint button."""
        w = QWidget(parent)
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(6)
        lbl = QLabel(title, w)
        lbl.setStyleSheet("font-weight: 600;")
        h.addWidget(lbl)
        if hint_id:
            h.addWidget(self._hint_btn(hint_id, parent=w))
        h.addStretch(1)
        return w

    def __init__(self, parent=None, font_size: float = 10.0):
        super().__init__(parent)
        # Match the darker / thicker group-box borders used on the MC Setup
        # and KORAL pages.
        self.setStyleSheet(
            "QGroupBox { border: 2px solid palette(shadow); border-radius: 4px;"
            " margin-top: 6px; padding-top: 6px; }"
        )
        self._hint_system: Optional[HintSystem] = None
        self._init_hints()
        self._curves: List[Dict[str, Any]] = []
        self._next_curve_id = 1
        self._font_size = float(font_size)
        self._color_cycle_idx = 0

        # Reference image overlaid under the curves (e.g. a JPEG/PNG from an
        # earlier simulation, used to compare edge structures).
        self._ref_image_path: Optional[str] = None
        self._ref_image_data: Optional[np.ndarray] = None
        self._ref_image_alpha: float = 0.5

        # Axes-level state
        self._axes_state: Dict[str, Any] = {
            "title":        "",
            "xlabel":       "",
            "xlabel_top":   "",
            "ylabel":       "",
            "ylabel_right": "",
            "log_x":         False,
            "log_y":         False,
            "log_x_top":     False,
            "log_y_right":   False,
            "grid":     True,
            "legend":   True,
            "legend_loc": "best",
        }
        self._axes_user_overrides: set = set()
        # View-preservation across re-renders (see _render_plot / _render_plot_keep_view).
        self._preserve_view_once: bool = False
        self._saved_view_limits = None

        root = QHBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(8)

        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        root.addWidget(splitter, 1)

        # ===== left: plot =====
        plot_box = QWidget()
        plot_layout = QVBoxLayout(plot_box)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(0)

        self.figure = Figure(figsize=(8, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self, coordinates=False)
        if hasattr(self.toolbar, "locLabel"):
            try:
                self.toolbar.locLabel.setVisible(False)
            except Exception:
                pass
        self._toolbar_spacer = QWidget()
        self._toolbar_spacer.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Preferred,
        )
        self.toolbar.addWidget(self._toolbar_spacer)
        self._cursor_status = QLabel("x: -, y: -")
        self._cursor_status.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._cursor_status.setStyleSheet("color: #000; font-size: 20px; padding: 0 6px;")
        self.toolbar.addWidget(self._cursor_status)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, 1)
        self.canvas.mpl_connect("motion_notify_event", self._on_canvas_motion)
        self.canvas.mpl_connect("figure_leave_event", self._on_canvas_leave)

        splitter.addWidget(plot_box)

        # ===== right: side panel (scrollable) =====
        side_scroll = QScrollArea()
        side_scroll.setWidgetResizable(True)
        side_scroll.setFrameShape(QFrame.Shape.NoFrame)
        side = QWidget()
        side.setMinimumWidth(360)
        side_layout = QVBoxLayout(side)
        side_layout.setContentsMargins(6, 6, 6, 6)
        side_layout.setSpacing(8)
        side_scroll.setWidget(side)

        # ---- Curves group ----
        curves_group = QGroupBox("")
        cg_layout = QVBoxLayout(curves_group)
        cg_layout.setContentsMargins(6, 6, 6, 6)
        cg_layout.setSpacing(4)
        cg_layout.addWidget(self._groupbox_header("Curves", hint_id="curves", parent=curves_group))

        self._curve_list = QListWidget()
        self._curve_list.setAlternatingRowColors(True)
        self._curve_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._curve_list.itemChanged.connect(self._on_curve_item_changed)
        self._curve_list.itemSelectionChanged.connect(self._on_selection_changed)
        cg_layout.addWidget(self._curve_list)

        list_btn_row = QHBoxLayout()
        self._btn_remove = QPushButton("Remove")
        self._btn_remove.setToolTip("Remove the selected curve")
        self._btn_remove.clicked.connect(self._remove_selected_curve)
        self._btn_clear = QPushButton("Clear All")
        self._btn_clear.clicked.connect(self._clear_all_curves)
        list_btn_row.addWidget(self._btn_remove)
        list_btn_row.addWidget(self._btn_clear)
        cg_layout.addLayout(list_btn_row)

        config_btn_row = QHBoxLayout()
        self._btn_save_cfg = QPushButton("Save Config…")
        self._btn_save_cfg.setToolTip(
            "Save this single-plot configuration (curves, styles, axes) with "
            "relative paths to the data sources"
        )
        self._btn_save_cfg.clicked.connect(self._save_plot_config)
        self._btn_load_cfg = QPushButton("Load Config…")
        self._btn_load_cfg.setToolTip(
            "Load a saved single-plot configuration; missing data sources are "
            "reported and skipped"
        )
        self._btn_load_cfg.clicked.connect(self._load_plot_config)
        config_btn_row.addWidget(self._btn_save_cfg)
        config_btn_row.addWidget(self._btn_load_cfg)
        cg_layout.addLayout(config_btn_row)

        side_layout.addWidget(curves_group)

        # ---- Load previous simulation output ----
        load_group = CollapsibleBox("Load Output Directory",
                                    hint_btn=self._hint_btn("load_directory"))
        load_layout = QVBoxLayout(load_group.content)
        load_layout.setContentsMargins(0, 0, 0, 0)
        load_layout.setSpacing(4)

        self._btn_load_dir = QPushButton("Load Output Directory…")
        self._btn_load_dir.setToolTip(
            "Pick a previous simulation output directory; all .his series are "
            "listed in the Dataset dropdown below"
        )
        self._btn_load_dir.clicked.connect(self._open_load_directory_dialog)
        load_layout.addWidget(self._btn_load_dir)

        ds_row = QHBoxLayout()
        ds_row.setSpacing(4)
        ds_row.addWidget(QLabel("Dataset:"))
        self._load_dataset_cmb = QComboBox()
        self._load_dataset_cmb.setEnabled(False)
        self._load_dataset_cmb.currentIndexChanged.connect(self._on_loaded_dataset_changed)
        ds_row.addWidget(self._load_dataset_cmb, 1)
        load_layout.addLayout(ds_row)

        lbl_row = QHBoxLayout()
        lbl_row.setSpacing(4)
        lbl_row.addWidget(QLabel("Label:"))
        self._load_label_edit = QLineEdit()
        self._load_label_edit.setEnabled(False)
        self._load_label_edit.setPlaceholderText("Display label for the new curve")
        lbl_row.addWidget(self._load_label_edit, 1)
        load_layout.addLayout(lbl_row)

        self._btn_add_loaded = QPushButton("Add as Curve")
        self._btn_add_loaded.setEnabled(False)
        self._btn_add_loaded.clicked.connect(self._add_loaded_as_curve)
        load_layout.addWidget(self._btn_add_loaded)

        self._load_status = QLabel("No file loaded.")
        self._load_status.setStyleSheet("color: #888; font-size: 10px;")
        self._load_status.setWordWrap(True)
        load_layout.addWidget(self._load_status)

        side_layout.addWidget(load_group)

        # State for loaded datasets
        self._loaded_datasets: List[Dict[str, Any]] = []
        self._loaded_data_path: Optional[str] = None

        # ---- Selected curve style ----
        self._style_group = CollapsibleBox("Curve Properties",
                                           hint_btn=self._hint_btn("selected_curve"))
        style_form = QFormLayout(self._style_group.content)
        style_form.setContentsMargins(0, 0, 0, 0)
        style_form.setHorizontalSpacing(6)
        style_form.setVerticalSpacing(4)
        style_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self._label_edit = QLineEdit()
        self._label_edit.textEdited.connect(self._on_label_edited)
        style_form.addRow("Label:", self._label_edit)

        color_row = QHBoxLayout()
        color_row.setSpacing(4)
        self._color_btn = QPushButton("")
        self._color_btn.setFixedSize(36, 22)
        self._color_btn.clicked.connect(self._pick_color)
        self._color_text = QLabel("—")
        self._color_text.setStyleSheet("color: #666; font-family: monospace;")
        color_row.addWidget(self._color_btn)
        color_row.addWidget(self._color_text)
        color_row.addStretch(1)
        color_holder = QWidget()
        color_holder.setLayout(color_row)
        style_form.addRow("Color:", color_holder)

        self._drawstyle_cmb = QComboBox()
        for name, val in _DRAWSTYLES:
            self._drawstyle_cmb.addItem(name, val)
        self._drawstyle_cmb.currentIndexChanged.connect(self._on_drawstyle_changed)
        style_form.addRow("Draw:", self._drawstyle_cmb)

        self._linestyle_cmb = QComboBox()
        for name, val in _LINESTYLES:
            self._linestyle_cmb.addItem(name, val)
        self._linestyle_cmb.currentIndexChanged.connect(self._on_linestyle_changed)
        style_form.addRow("Line:", self._linestyle_cmb)

        self._linewidth_spin = QDoubleSpinBox()
        self._linewidth_spin.setRange(0.1, 20.0)
        self._linewidth_spin.setDecimals(2)
        self._linewidth_spin.setSingleStep(0.2)
        self._linewidth_spin.setValue(1.4)
        self._linewidth_spin.valueChanged.connect(self._on_linewidth_changed)
        style_form.addRow("Width:", self._linewidth_spin)

        self._marker_cmb = QComboBox()
        for name, val in _MARKERS:
            self._marker_cmb.addItem(name, val)
        self._marker_cmb.currentIndexChanged.connect(self._on_marker_changed)
        style_form.addRow("Marker:", self._marker_cmb)

        self._marker_size_spin = QDoubleSpinBox()
        self._marker_size_spin.setRange(0.0, 40.0)
        self._marker_size_spin.setDecimals(1)
        self._marker_size_spin.setSingleStep(0.5)
        self._marker_size_spin.setValue(4.0)
        self._marker_size_spin.valueChanged.connect(self._on_markersize_changed)
        style_form.addRow("M. size:", self._marker_size_spin)

        self._alpha_spin = QDoubleSpinBox()
        self._alpha_spin.setRange(0.0, 1.0)
        self._alpha_spin.setDecimals(2)
        self._alpha_spin.setSingleStep(0.05)
        self._alpha_spin.setValue(1.0)
        self._alpha_spin.valueChanged.connect(self._on_alpha_changed)
        style_form.addRow("Alpha:", self._alpha_spin)

        self._zorder_spin = QSpinBox()
        self._zorder_spin.setRange(0, 999)
        self._zorder_spin.setValue(2)
        self._zorder_spin.valueChanged.connect(self._on_zorder_changed)
        style_form.addRow("Z-order:", self._zorder_spin)

        self._axis_y_cmb = QComboBox()
        self._axis_y_cmb.addItem("Left",  "left")
        self._axis_y_cmb.addItem("Right", "right")
        self._axis_y_cmb.currentIndexChanged.connect(self._on_axis_y_changed)
        style_form.addRow("Y axis:", self._axis_y_cmb)

        self._axis_x_cmb = QComboBox()
        self._axis_x_cmb.addItem("Bottom", "bottom")
        self._axis_x_cmb.addItem("Top",    "top")
        self._axis_x_cmb.currentIndexChanged.connect(self._on_axis_x_changed)
        style_form.addRow("X axis:", self._axis_x_cmb)

        self._bin_combine_spin = QSpinBox()
        self._bin_combine_spin.setRange(1, 6)
        self._bin_combine_spin.setValue(1)
        self._bin_combine_spin.setToolTip(
            "Sum N adjacent histogram bins for this curve (1 = no combining)."
        )
        self._bin_combine_spin.valueChanged.connect(self._on_bin_combine_changed)
        style_form.addRow("Combine bins:", self._bin_combine_spin)

        self._twod_render_label = QLabel("2D render:")
        self._twod_mode_cmb = QComboBox()
        self._twod_mode_cmb.addItem("Color Mesh", "mesh")
        self._twod_mode_cmb.addItem("Contours", "contour")
        self._twod_mode_cmb.currentIndexChanged.connect(self._on_twod_render_mode_changed)
        style_form.addRow(self._twod_render_label, self._twod_mode_cmb)

        self._twod_contour_label = QLabel("Contours:")
        contour_row = QHBoxLayout()
        contour_row.setContentsMargins(0, 0, 0, 0)
        contour_row.setSpacing(4)
        contour_row.addWidget(QLabel("Levels"))
        self._twod_levels_spin = QSpinBox()
        self._twod_levels_spin.setRange(2, 40)
        self._twod_levels_spin.setValue(8)
        self._twod_levels_spin.valueChanged.connect(self._on_twod_levels_changed)
        contour_row.addWidget(self._twod_levels_spin)
        self._twod_label_contours_check = QCheckBox("Label contours")
        self._twod_label_contours_check.toggled.connect(self._on_twod_contour_labels_toggled)
        contour_row.addWidget(self._twod_label_contours_check)
        contour_row.addStretch(1)
        self._twod_contour_holder = QWidget()
        self._twod_contour_holder.setLayout(contour_row)
        style_form.addRow(self._twod_contour_label, self._twod_contour_holder)

        side_layout.addWidget(self._style_group)
        self._style_group.setEnabled(False)
        self._twod_render_label.setVisible(False)
        self._twod_mode_cmb.setVisible(False)
        self._twod_contour_label.setVisible(False)
        self._twod_contour_holder.setVisible(False)

        # ---- Convolution / scan area (per-curve) ----
        self._conv_group = CollapsibleBox("Gaussian Convolution",
                                          hint_btn=self._hint_btn("convolution"))
        conv_layout = QVBoxLayout(self._conv_group.content)
        conv_layout.setContentsMargins(0, 0, 0, 0)
        conv_layout.setSpacing(4)

        # Depth-profile convolution controls
        self._conv_depth_widget = QWidget()
        conv_depth_layout = QVBoxLayout(self._conv_depth_widget)
        conv_depth_layout.setContentsMargins(0, 0, 0, 0)
        conv_depth_layout.setSpacing(4)

        self._conv_check = QCheckBox("Apply convolution")
        self._conv_check.toggled.connect(self._on_conv_toggled)

        self._conv_info = QLabel(
            "Lateral resolution corresponds physically to the former beam sigma. "
            "Use depth resolution for projected/depth broadening."
        )
        self._conv_info.setWordWrap(True)
        self._conv_info.setStyleSheet("color: #666; font-size: 10px;")
        conv_depth_layout.addWidget(self._conv_info)

        self._sigma_row_widget = QWidget()
        sigma_row = QHBoxLayout(self._sigma_row_widget)
        sigma_row.setContentsMargins(0, 0, 0, 0)
        sigma_row.addWidget(QLabel("Lateral resolution (Å):"))
        self._sigma_spin = QDoubleSpinBox()
        self._sigma_spin.setRange(0.0, 1_000_000.0)
        self._sigma_spin.setDecimals(2)
        self._sigma_spin.setValue(0.0)
        self._sigma_spin.setSingleStep(5.0)
        self._sigma_spin.valueChanged.connect(self._on_sigma_changed)
        sigma_row.addWidget(self._sigma_spin, 1)
        conv_depth_layout.addWidget(self._sigma_row_widget)

        self._depth_sigma_row_widget = QWidget()
        depth_sigma_row = QHBoxLayout(self._depth_sigma_row_widget)
        depth_sigma_row.setContentsMargins(0, 0, 0, 0)
        depth_sigma_row.addWidget(QLabel("Depth resolution (Å):"))
        self._depth_sigma_spin = QDoubleSpinBox()
        self._depth_sigma_spin.setRange(0.0, 1_000_000.0)
        self._depth_sigma_spin.setDecimals(2)
        self._depth_sigma_spin.setValue(0.0)
        self._depth_sigma_spin.setSingleStep(5.0)
        self._depth_sigma_spin.valueChanged.connect(self._on_depth_sigma_changed)
        depth_sigma_row.addWidget(self._depth_sigma_spin, 1)
        conv_depth_layout.addWidget(self._depth_sigma_row_widget)
        conv_layout.addWidget(self._conv_depth_widget)

        # Scan-area controls (2-D / lateral curves only)
        self._scan_area_widget = QWidget()
        scan_layout = QVBoxLayout(self._scan_area_widget)
        scan_layout.setContentsMargins(0, 0, 0, 0)
        scan_layout.setSpacing(4)

        scan_row = QHBoxLayout()
        scan_row.addWidget(QLabel("Lateral mask / scan area:"))
        self._scan_area_min = QDoubleSpinBox()
        self._scan_area_max = QDoubleSpinBox()
        for _sp in (self._scan_area_min, self._scan_area_max):
            _sp.setRange(-1_000_000.0, 1_000_000.0)
            _sp.setDecimals(2)
            _sp.setValue(0.0)
            _sp.setSingleStep(5.0)
            _sp.valueChanged.connect(self._on_scan_area_changed)
        scan_row.addWidget(self._scan_area_min, 1)
        scan_row.addWidget(QLabel("–"))
        scan_row.addWidget(self._scan_area_max, 1)
        scan_layout.addLayout(scan_row)

        self._sigma_2d_x_widget = QWidget()
        sigma_2d_x_row = QHBoxLayout(self._sigma_2d_x_widget)
        sigma_2d_x_row.setContentsMargins(0, 0, 0, 0)
        sigma_2d_x_row.addWidget(QLabel("2D depth resolution X (Å):"))
        self._sigma_2d_x_spin = QDoubleSpinBox()
        self._sigma_2d_x_spin.setRange(0.0, 1_000_000.0)
        self._sigma_2d_x_spin.setDecimals(2)
        self._sigma_2d_x_spin.setValue(0.0)
        self._sigma_2d_x_spin.setSingleStep(5.0)
        self._sigma_2d_x_spin.valueChanged.connect(self._on_sigma_2d_changed)
        sigma_2d_x_row.addWidget(self._sigma_2d_x_spin, 1)
        scan_layout.addWidget(self._sigma_2d_x_widget)

        self._sigma_2d_y_widget = QWidget()
        sigma_2d_y_row = QHBoxLayout(self._sigma_2d_y_widget)
        sigma_2d_y_row.setContentsMargins(0, 0, 0, 0)
        sigma_2d_y_row.addWidget(QLabel("2D lateral resolution Y (Å):"))
        self._sigma_2d_y_spin = QDoubleSpinBox()
        self._sigma_2d_y_spin.setRange(0.0, 1_000_000.0)
        self._sigma_2d_y_spin.setDecimals(2)
        self._sigma_2d_y_spin.setValue(0.0)
        self._sigma_2d_y_spin.setSingleStep(5.0)
        self._sigma_2d_y_spin.valueChanged.connect(self._on_sigma_2d_changed)
        sigma_2d_y_row.addWidget(self._sigma_2d_y_spin, 1)
        scan_layout.addWidget(self._sigma_2d_y_widget)
        conv_layout.addWidget(self._scan_area_widget)
        self._sigma_2d_x_widget.setVisible(False)
        self._sigma_2d_y_widget.setVisible(False)

        self._btn_add_convolution_curve = QPushButton("Add convolution as dataset")
        self._btn_add_convolution_curve.clicked.connect(self._add_convolution_as_dataset)
        conv_layout.addWidget(self._btn_add_convolution_curve)

        # Multi-σ comparison: convolve the same source data with several σ
        # values at once and overlay them for direct comparison (per prof's
        # idea). For 2-D curves the comparison is rendered as contour lines.
        multi_sigma_row = QHBoxLayout()
        multi_sigma_row.setContentsMargins(0, 0, 0, 0)
        multi_sigma_row.setSpacing(4)
        multi_sigma_row.addWidget(QLabel("Compare σ (Å):"))
        self._multi_sigma_edit = QLineEdit()
        self._multi_sigma_edit.setPlaceholderText("e.g. 10, 25, 50")
        self._multi_sigma_edit.setToolTip(
            "Enter several σ values separated by commas. Each value convolves "
            "the selected curve's original data and is added as a separate "
            "curve (2-D curves are added as contour lines)."
        )
        multi_sigma_row.addWidget(self._multi_sigma_edit, 1)
        conv_layout.addLayout(multi_sigma_row)

        self._btn_add_multi_sigma = QPushButton("Add σ comparison")
        self._btn_add_multi_sigma.setToolTip(
            "Convolve the selected curve with each σ value above and overlay "
            "the results for comparison."
        )
        self._btn_add_multi_sigma.clicked.connect(self._add_multi_sigma_comparison)
        conv_layout.addWidget(self._btn_add_multi_sigma)

        # "Apply convolution" toggle sits at the bottom of the panel, after the
        # σ and scan-area inputs, so the user sets parameters first and enables
        # the effect last.
        conv_layout.addWidget(self._conv_check)

        # Keep convolution below the style/axes panels; curves and loading stay highest-priority.

        # ---- Axes & legend ----
        axes_group = CollapsibleBox("Axes & Legend", hint_btn=self._hint_btn("axes"))
        axes_form = QFormLayout(axes_group.content)
        axes_form.setContentsMargins(0, 0, 0, 0)
        axes_form.setHorizontalSpacing(6)
        axes_form.setVerticalSpacing(4)
        axes_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        self._title_edit = QLineEdit()
        self._title_edit.textEdited.connect(lambda t: self._on_axes_edit("title", t))
        axes_form.addRow("Title:", self._title_edit)

        self._xlabel_edit = QLineEdit()
        self._xlabel_edit.textEdited.connect(lambda t: self._on_axes_edit("xlabel", t))
        self._xlabel_label = QLabel("X label:")
        axes_form.addRow(self._xlabel_label, self._xlabel_edit)

        self._xlabel_top_edit = QLineEdit()
        self._xlabel_top_edit.textEdited.connect(lambda t: self._on_axes_edit("xlabel_top", t))
        self._xlabel_top_label = QLabel("X label (top):")
        axes_form.addRow(self._xlabel_top_label, self._xlabel_top_edit)

        self._ylabel_edit = QLineEdit()
        self._ylabel_edit.textEdited.connect(lambda t: self._on_axes_edit("ylabel", t))
        self._ylabel_label = QLabel("Y label:")
        axes_form.addRow(self._ylabel_label, self._ylabel_edit)

        self._ylabel_right_edit = QLineEdit()
        self._ylabel_right_edit.textEdited.connect(lambda t: self._on_axes_edit("ylabel_right", t))
        self._ylabel_right_label = QLabel("Y label (right):")
        axes_form.addRow(self._ylabel_right_label, self._ylabel_right_edit)

        scale_row = QHBoxLayout()
        self._log_x_check = QCheckBox("Log X")
        self._log_y_check = QCheckBox("Log Y")
        self._log_x_check.toggled.connect(lambda v: self._on_axes_flag("log_x", v))
        self._log_y_check.toggled.connect(lambda v: self._on_axes_flag("log_y", v))
        scale_row.addWidget(self._log_x_check)
        scale_row.addWidget(self._log_y_check)
        scale_row.addStretch(1)
        scale_holder = QWidget()
        scale_holder.setLayout(scale_row)
        axes_form.addRow("Scale:", scale_holder)

        scale_row2 = QHBoxLayout()
        self._log_x_top_check = QCheckBox("Log X (top)")
        self._log_y_right_check = QCheckBox("Log Y (right)")
        self._log_x_top_check.toggled.connect(lambda v: self._on_axes_flag("log_x_top", v))
        self._log_y_right_check.toggled.connect(lambda v: self._on_axes_flag("log_y_right", v))
        scale_row2.addWidget(self._log_x_top_check)
        scale_row2.addWidget(self._log_y_right_check)
        scale_row2.addStretch(1)
        self._scale2_holder = QWidget()
        self._scale2_holder.setLayout(scale_row2)
        self._scale2_label = QLabel("Scale (sec.):")
        axes_form.addRow(self._scale2_label, self._scale2_holder)

        misc_row = QHBoxLayout()
        self._grid_check = QCheckBox("Grid")
        self._grid_check.setChecked(True)
        self._grid_check.toggled.connect(lambda v: self._on_axes_flag("grid", v))
        self._legend_check = QCheckBox("Legend")
        self._legend_check.setChecked(True)
        self._legend_check.toggled.connect(lambda v: self._on_axes_flag("legend", v))
        misc_row.addWidget(self._grid_check)
        misc_row.addWidget(self._legend_check)
        misc_row.addStretch(1)
        misc_holder = QWidget()
        misc_holder.setLayout(misc_row)
        axes_form.addRow("Display:", misc_holder)

        self._legend_loc_cmb = QComboBox()
        for loc in _LEGEND_POSITIONS:
            self._legend_loc_cmb.addItem(loc)
        self._legend_loc_cmb.currentTextChanged.connect(
            lambda t: self._on_axes_edit("legend_loc", t)
        )
        axes_form.addRow("Legend at:", self._legend_loc_cmb)

        side_layout.addWidget(axes_group)
        side_layout.addWidget(self._conv_group)

        # ---- Reference image overlay ----
        ref_group = CollapsibleBox("Reference Image", hint_btn=self._hint_btn("reference"))
        self._ref_group = ref_group
        ref_layout = QVBoxLayout(ref_group.content)
        ref_layout.setContentsMargins(0, 0, 0, 0)
        ref_layout.setSpacing(4)

        self._ref_load_btn = QPushButton("Load Image…")
        self._ref_load_btn.setToolTip(
            "Overlay a JPEG/PNG (e.g. from an earlier simulation) under the curves."
        )
        self._ref_load_btn.clicked.connect(self._pick_reference_image)
        ref_layout.addWidget(self._ref_load_btn)

        self._ref_clear_btn = QPushButton("Clear Image")
        self._ref_clear_btn.clicked.connect(self._clear_reference_image)
        self._ref_clear_btn.setEnabled(False)
        ref_layout.addWidget(self._ref_clear_btn)

        alpha_row = QHBoxLayout()
        alpha_row.addWidget(QLabel("Opacity:"))
        self._ref_alpha_spin = QDoubleSpinBox()
        self._ref_alpha_spin.setRange(0.0, 1.0)
        self._ref_alpha_spin.setDecimals(2)
        self._ref_alpha_spin.setSingleStep(0.05)
        self._ref_alpha_spin.setValue(self._ref_image_alpha)
        self._ref_alpha_spin.valueChanged.connect(self._on_ref_alpha_changed)
        alpha_row.addWidget(self._ref_alpha_spin, 1)
        ref_layout.addLayout(alpha_row)

        self._ref_status = QLabel("No image loaded.")
        self._ref_status.setStyleSheet("color: #888; font-size: 10px;")
        self._ref_status.setWordWrap(True)
        ref_layout.addWidget(self._ref_status)

        side_layout.addWidget(ref_group)

        # ---- Stats for selected curve ----
        self._stats_group = CollapsibleBox("Stats (Selected)")
        stats_layout = QVBoxLayout(self._stats_group.content)
        stats_layout.setContentsMargins(0, 0, 0, 0)
        self._stats_table = QTableWidget(0, 2)
        self._stats_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self._stats_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self._stats_table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.Stretch
        )
        self._stats_table.verticalHeader().setVisible(False)
        self._stats_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._stats_table.setAlternatingRowColors(True)
        stats_layout.addWidget(self._stats_table)
        side_layout.addWidget(self._stats_group)
        self._stats_group.setVisible(False)

        side_layout.addStretch(1)

        splitter.addWidget(side_scroll)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([760, 400])

        self._update_secondary_axes_visibility()
        self._render_plot()

    # =====================================================================
    # Public API
    # =====================================================================

    def set_font_size(self, size: float) -> None:
        self._font_size = float(size)
        self._render_plot()

    def set_beam_from_top(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if getattr(self, "_beam_from_top", True) == enabled:
            return
        self._beam_from_top = enabled
        self._render_plot()

    # Backwards-compat shim — callers that used set_plot now add instead.
    def set_plot(self, plot_id: str, plot_info: Dict[str, Any]) -> None:
        self._curves = []
        self._refresh_curve_list()
        self.add_plot(plot_id, plot_info)

    def add_plot(self, plot_id: str, plot_info: Dict[str, Any]) -> int:
        """Append every visible series of *plot_info* as a new curve."""
        if not isinstance(plot_info, dict):
            return 0

        requested_kind = "2d" if plot_info.get("is_2d") else "1d"
        if self._curves:
            existing_kind = "2d" if any(c.get("is_2d") for c in self._curves) else "1d"
            if existing_kind != requested_kind:
                QMessageBox.information(
                    self,
                    "Single Plot",
                    "1D and 2D histograms cannot be mixed in Single Plot. Clear the current curves first.",
                )
                return 0

        if plot_info.get("is_2d") and callable(plot_info.get("plot_func")):
            source_name = str(plot_info.get("name", plot_id))
            x_label = str(plot_info.get("x_label", "") or "")
            y_label = str(plot_info.get("y_label", "") or "")
            curve = {
                "id":          self._next_curve_id,
                "source_id":   plot_id,
                "source_name": source_name,
                "series_idx":  0,
                "label":       source_name,
                "visible":     True,
                "is_2d":       True,
                "plot_func":   plot_info["plot_func"],
                "stats_func":  plot_info.get("stats_func"),
                "x_values_2d": np.asarray(
                    plot_info["x_values_2d"] if plot_info.get("x_values_2d") is not None else [],
                    dtype=float,
                ),
                "y_values_2d": np.asarray(
                    plot_info["y_values_2d"] if plot_info.get("y_values_2d") is not None else [],
                    dtype=float,
                ),
                "z_values_2d": np.asarray(
                    plot_info["z_values_2d"] if plot_info.get("z_values_2d") is not None else [],
                    dtype=float,
                ),
                "x_label":     x_label,
                "y_label":     y_label,
                "colorbar_label": str(plot_info.get("colorbar_label", "") or ""),
                "color":       self._next_cycle_color(),
                "drawstyle":   "default",
                "linestyle":   "-",
                "linewidth":   1.0,
                "marker":      "None",
                "markersize":  4.0,
                "alpha":       1.0,
                "zorder":      2 + len(self._curves),
                "conv_enabled": False,
                "conv_sigma":   0.0,
                "conv_sigma_lateral": 0.0,
                "conv_sigma_depth": 0.0,
                "scan_area":    [0.0, 0.0],
                "bin_factor":   1,
                "is_depth_profile": False,
                "axis_x":       "bottom",
                "axis_y":       "left",
                "render_mode":  "mesh" if not self._curves else "contour",
                "contour_levels": 8,
                "label_contours": False,
            }
            self._next_curve_id += 1
            self._curves.append(curve)
            # Take title/labels from the 2D plot unless the user already
            # customized them.
            if "title" not in self._axes_user_overrides:
                self._axes_state["title"] = source_name
                self._title_edit.blockSignals(True)
                self._title_edit.setText(source_name)
                self._title_edit.blockSignals(False)
            if "xlabel" not in self._axes_user_overrides:
                self._axes_state["xlabel"] = x_label
                self._xlabel_edit.blockSignals(True)
                self._xlabel_edit.setText(x_label)
                self._xlabel_edit.blockSignals(False)
            if "ylabel" not in self._axes_user_overrides:
                self._axes_state["ylabel"] = y_label
                self._ylabel_edit.blockSignals(True)
                self._ylabel_edit.setText(y_label)
                self._ylabel_edit.blockSignals(False)
            self._refresh_curve_list(select_last=True)
            self._render_plot()
            return 1

        x_vals = plot_info.get("x_values")
        cols = plot_info.get("columns")
        if x_vals is None or cols is None:
            return 0
        x_arr = np.asarray(x_vals, dtype=float)

        labels = list(plot_info.get("series_labels", []) or [])
        src_colors = list(plot_info.get("colors", []) or [])
        source_name = str(plot_info.get("name", plot_id))
        is_depth = bool(plot_info.get("is_depth_profile", False))
        x_label = str(plot_info.get("x_label", "") or "")
        y_label = str(plot_info.get("y_label", "") or "")
        plot_stats_func = plot_info.get("stats_func")

        # Adopt the first plot's axes labels and title as defaults — only if
        # the user has not already customized them.
        if not self._curves:
            if "title" not in self._axes_user_overrides:
                self._axes_state["title"] = source_name
                self._title_edit.blockSignals(True)
                self._title_edit.setText(source_name)
                self._title_edit.blockSignals(False)
            if "xlabel" not in self._axes_user_overrides:
                self._axes_state["xlabel"] = x_label
                self._xlabel_edit.blockSignals(True)
                self._xlabel_edit.setText(x_label)
                self._xlabel_edit.blockSignals(False)
            if "ylabel" not in self._axes_user_overrides:
                self._axes_state["ylabel"] = y_label
                self._ylabel_edit.blockSignals(True)
                self._ylabel_edit.setText(y_label)
                self._ylabel_edit.blockSignals(False)

        added = 0
        for i, col_data in enumerate(cols):
            y_arr = np.asarray(col_data, dtype=float)
            if y_arr.size != x_arr.size or y_arr.size < 2:
                continue
            series_label = labels[i] if i < len(labels) else f"Series {i}"
            display_label = f"{series_label} — {source_name}"
            color = src_colors[i] if i < len(src_colors) else self._next_cycle_color()
            curve = {
                "id":          self._next_curve_id,
                "source_id":   plot_id,
                "source_name": source_name,
                "series_idx":  i,
                "x":           x_arr,
                "y_original":  y_arr,
                "label":       display_label,
                "visible":     True,
                "color":       color,
                "drawstyle":   "steps-mid",
                "linestyle":   "-",
                "linewidth":   1.4,
                "marker":      "None",
                "markersize":  4.0,
                "alpha":       1.0,
                "zorder":      2 + len(self._curves),
                "x_label":     x_label,
                "y_label":     y_label,
                "is_depth_profile": is_depth,
                "conv_enabled": False,
                "conv_sigma":   0.0,
                "conv_sigma_lateral": 0.0,
                "conv_sigma_depth": 0.0,
                "scan_area":    [0.0, 0.0],
                "bin_factor":   1,
                "stats_func":   plot_stats_func,
                "axis_x":       "bottom",
                "axis_y":       "left",
            }
            self._next_curve_id += 1
            self._curves.append(curve)
            added += 1

        if added > 0:
            self._refresh_curve_list(select_last=True)
            self._render_plot()
        return added

    def add_dataset_as_curve(self, dataset: Dict[str, Any], label: str = "") -> bool:
        """Append a single externally-loaded dataset as a curve."""
        requested_kind = "2d" if bool(dataset.get("is_2d")) else "1d"
        if self._curves:
            existing_kind = "2d" if any(c.get("is_2d") for c in self._curves) else "1d"
            if existing_kind != requested_kind:
                QMessageBox.information(
                    self,
                    "Single Plot",
                    "1D and 2D histograms cannot be mixed in Single Plot. Clear the current curves first.",
                )
                return False

        if bool(dataset.get("is_2d")):
            plot_info = dataset.get("plot_info")
            if not isinstance(plot_info, dict):
                return False
            added = self.add_plot(str(dataset.get("id") or "loaded_2d"), dict(plot_info))
            if added > 0 and label.strip():
                self._curves[-1]["label"] = label.strip()
                self._refresh_curve_list(select_last=True)
                self._render_plot()
            return added > 0

        x = _to_numeric_vector(dataset.get("x"))
        y = _to_numeric_vector(dataset.get("y"))
        if x is None or y is None or x.size != y.size:
            return False

        x_label = str(dataset.get("x_label") or "")
        y_label = str(dataset.get("y_label") or "")
        ds_name = str(dataset.get("name") or "dataset")
        tag = (
            dataset.get("_dir_name")
            or os.path.basename(str(dataset.get("path") or "loaded data"))
        )
        display_label = (label or "").strip() or f"{ds_name} ({tag})"

        # Set initial axes labels if this is the very first curve.
        if not self._curves:
            if "title" not in self._axes_user_overrides:
                self._axes_state["title"] = ds_name
                self._title_edit.blockSignals(True)
                self._title_edit.setText(ds_name)
                self._title_edit.blockSignals(False)
            if "xlabel" not in self._axes_user_overrides and x_label:
                self._axes_state["xlabel"] = x_label
                self._xlabel_edit.blockSignals(True)
                self._xlabel_edit.setText(x_label)
                self._xlabel_edit.blockSignals(False)
            if "ylabel" not in self._axes_user_overrides and y_label:
                self._axes_state["ylabel"] = y_label
                self._ylabel_edit.blockSignals(True)
                self._ylabel_edit.setText(y_label)
                self._ylabel_edit.blockSignals(False)

        curve = {
            "id":          self._next_curve_id,
            "source_id":   str(dataset.get("id") or "loaded"),
            "source_name": tag,
            "source_path": str(dataset.get("path") or ""),
            "series_idx":  0,
            "x":           x,
            "y_original":  y,
            "label":       display_label,
            "visible":     True,
            "color":       self._next_cycle_color(),
            "drawstyle":   "default",
            "linestyle":   "-",
            "linewidth":   1.4,
            "marker":      "None",
            "markersize":  4.0,
            "alpha":       1.0,
            "zorder":      2 + len(self._curves),
            "x_label":     x_label,
            "y_label":     y_label,
            # If the dataset came from a known depth profile (e.g. via
            # directory load) use that; otherwise default to True so the
            # user can still try convolution on it.
            "is_depth_profile": bool(dataset.get("is_depth_profile", True)),
            "conv_enabled": False,
            "conv_sigma":   0.0,
            "conv_sigma_lateral": 0.0,
            "conv_sigma_depth": 0.0,
            "scan_area":    [0.0, 0.0],
            "stats_func":   None,
            "axis_x":       "bottom",
            "axis_y":       "left",
        }
        self._next_curve_id += 1
        self._curves.append(curve)
        self._refresh_curve_list(select_last=True)
        self._render_plot()
        return True

    # =====================================================================
    # Load-from-directory helpers
    # =====================================================================

    def _open_load_directory_dialog(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Simulation Output Directory",
            get_last_used_directory(),
            QFileDialog.Option.ShowDirsOnly,
        )
        if not directory:
            return
        remember_last_used_path(directory)
        # Lazy import to avoid a circular dependency with mcresults_page.
        try:
            from ui.pages.mcresults_page import _build_plots_from_directory
        except ImportError:
            try:
                from OpenSRIM.ui.pages.mcresults_page import _build_plots_from_directory  # type: ignore
            except ImportError as exc:
                QMessageBox.warning(self, "Load Directory",
                                    f"Could not load directory reader:\n{exc}")
                return
        try:
            plots, _numerical = _build_plots_from_directory(directory)
        except Exception as exc:
            QMessageBox.warning(self, "Load Directory",
                                f"Failed to read directory:\n{exc}")
            return
        if not plots:
            QMessageBox.information(
                self, "Load Directory",
                "No .his/.mom result files were found in the selected directory.",
            )
            return

        # Flatten: each (plot, series) pair becomes one dataset entry so
        # users pick a specific curve from the combo box.
        dir_name = os.path.basename(os.path.normpath(directory)) or directory
        datasets: List[Dict[str, Any]] = []
        for plot_id, info in plots.items():
            if info.get("is_2d"):
                datasets.append({
                    "id": plot_id,
                    "name": str(info.get("name", plot_id)),
                    "path": directory,
                    "is_2d": True,
                    "plot_info": dict(info),
                    "_dir_name": dir_name,
                })
                continue
            x_vals = info.get("x_values")
            cols = info.get("columns")
            if x_vals is None or cols is None:
                continue
            x_arr = np.asarray(x_vals, dtype=float)
            labels = list(info.get("series_labels", []) or [])
            x_label = str(info.get("x_label", "") or "")
            y_label = str(info.get("y_label", "") or "")
            plot_name = str(info.get("name", plot_id))
            is_depth = bool(info.get("is_depth_profile", False))
            for i, col in enumerate(cols):
                y_arr = np.asarray(col, dtype=float)
                if y_arr.size != x_arr.size or y_arr.size < 2:
                    continue
                series_label = labels[i] if i < len(labels) else f"Series {i}"
                datasets.append({
                    "id":      f"{plot_id}__{i}",
                    "name":    f"{plot_name} — {series_label}",
                    "x":       x_arr,
                    "y":       y_arr,
                    "x_label": x_label,
                    "y_label": y_label,
                    "path":    directory,
                    "is_depth_profile": is_depth,
                    "_dir_name": dir_name,
                })

        if not datasets:
            QMessageBox.information(
                self, "Load Directory",
                "No plottable series were found in the result files.",
            )
            return

        self._loaded_data_path = directory
        self._loaded_datasets = datasets

        self._load_dataset_cmb.blockSignals(True)
        self._load_dataset_cmb.clear()
        for ds in datasets:
            self._load_dataset_cmb.addItem(str(ds["name"]))
        self._load_dataset_cmb.blockSignals(False)

        self._load_dataset_cmb.setEnabled(True)
        self._load_label_edit.setEnabled(True)
        self._btn_add_loaded.setEnabled(True)
        self._load_dataset_cmb.setCurrentIndex(0)
        self._load_status.setText(
            f"Loaded {len(datasets)} series from directory “{dir_name}”"
        )

    def _on_loaded_dataset_changed(self, idx: int) -> None:
        if 0 <= idx < len(self._loaded_datasets):
            ds = self._loaded_datasets[idx]
            tag = (
                ds.get("_dir_name")
                or os.path.basename(self._loaded_data_path or "")
                or "loaded"
            )
            base = str(ds.get("name") or "dataset")
            default_label = f"{base} ({tag})"
            self._load_label_edit.blockSignals(True)
            self._load_label_edit.setText(default_label)
            self._load_label_edit.blockSignals(False)

    def _add_loaded_as_curve(self) -> None:
        idx = self._load_dataset_cmb.currentIndex()
        if not (0 <= idx < len(self._loaded_datasets)):
            return
        ds = self._loaded_datasets[idx]
        label = self._load_label_edit.text().strip()
        if not self.add_dataset_as_curve(ds, label=label):
            QMessageBox.warning(
                self, "Add Curve",
                "Dataset does not contain valid plottable data.",
            )

    # =====================================================================
    # Curve list helpers
    # =====================================================================

    def _refresh_curve_list(self, *, select_last: bool = False) -> None:
        self._curve_list.blockSignals(True)
        self._curve_list.clear()
        for curve in self._curves:
            item = QListWidgetItem(str(curve["label"]))
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsEnabled
            )
            item.setCheckState(
                Qt.CheckState.Checked if curve["visible"] else Qt.CheckState.Unchecked
            )
            item.setIcon(_color_icon(curve["color"]))
            item.setData(Qt.ItemDataRole.UserRole, int(curve["id"]))
            self._curve_list.addItem(item)
        self._curve_list.blockSignals(False)

        self._update_secondary_axes_visibility()

        if select_last and self._curves:
            self._curve_list.setCurrentRow(len(self._curves) - 1)
        else:
            self._on_selection_changed()

    def _selected_curve(self) -> Optional[Dict[str, Any]]:
        item = self._curve_list.currentItem()
        if item is None:
            return None
        cid = int(item.data(Qt.ItemDataRole.UserRole) or -1)
        for c in self._curves:
            if int(c["id"]) == cid:
                return c
        return None

    def _on_curve_item_changed(self, item: QListWidgetItem) -> None:
        cid = int(item.data(Qt.ItemDataRole.UserRole) or -1)
        for c in self._curves:
            if int(c["id"]) != cid:
                continue
            c["visible"] = (item.checkState() == Qt.CheckState.Checked)
            new_label = item.text().strip()
            if new_label and new_label != c["label"]:
                c["label"] = new_label
            break
        self._render_plot()
        sel = self._selected_curve()
        if sel is not None and cid == int(sel["id"]):
            self._label_edit.blockSignals(True)
            self._label_edit.setText(sel["label"])
            self._label_edit.blockSignals(False)

    def _on_selection_changed(self) -> None:
        curve = self._selected_curve()
        has = curve is not None
        self._style_group.setEnabled(has)
        self._stats_group.setVisible(has)
        # Convolution panel stays interactive even without a selection — the
        # prof reported it looking "broken" when greyed out. The individual
        # handlers (_on_conv_toggled, _on_sigma_changed, ...) silently no-op
        # or auto-select the first curve when the user interacts.
        self._conv_group.setEnabled(True)
        if not has:
            self._set_twod_controls_visible(False)
            self._stats_table.setRowCount(0)
            return

        self._label_edit.blockSignals(True)
        self._label_edit.setText(str(curve["label"]))
        self._label_edit.blockSignals(False)

        self._update_color_button(curve["color"])

        self._set_combo_data(self._linestyle_cmb, curve.get("linestyle", "-"))
        self._linewidth_spin.blockSignals(True)
        self._linewidth_spin.setValue(float(curve.get("linewidth", 1.4)))
        self._linewidth_spin.blockSignals(False)
        self._alpha_spin.blockSignals(True)
        self._alpha_spin.setValue(float(curve.get("alpha", 1.0)))
        self._alpha_spin.blockSignals(False)
        self._zorder_spin.blockSignals(True)
        self._zorder_spin.setValue(int(curve.get("zorder", 2)))
        self._zorder_spin.blockSignals(False)

        if curve.get("is_2d"):
            self._set_twod_controls_visible(True)
            self._drawstyle_cmb.setEnabled(False)
            self._marker_cmb.setEnabled(False)
            self._marker_size_spin.setEnabled(False)
            self._axis_y_cmb.setEnabled(False)
            self._axis_x_cmb.setEnabled(False)
            self._bin_combine_spin.setEnabled(False)
            self._conv_depth_widget.setVisible(True)
            self._scan_area_widget.setVisible(True)
            self._conv_check.blockSignals(True)
            self._conv_check.setChecked(bool(curve.get("conv_enabled", False)))
            self._conv_check.blockSignals(False)
            self._sigma_spin.blockSignals(True)
            self._sigma_spin.setValue(float(curve.get("conv_sigma_lateral", 0.0)))
            self._sigma_spin.blockSignals(False)
            self._depth_sigma_spin.blockSignals(True)
            self._depth_sigma_spin.setValue(float(curve.get("conv_sigma_depth", 0.0)))
            self._depth_sigma_spin.blockSignals(False)
            self._sigma_2d_x_spin.blockSignals(True)
            self._sigma_2d_x_spin.setValue(float(curve.get("conv_sigma_depth", 0.0)))
            self._sigma_2d_x_spin.blockSignals(False)
            self._sigma_2d_y_spin.blockSignals(True)
            self._sigma_2d_y_spin.setValue(float(curve.get("conv_sigma_lateral", 0.0)))
            self._sigma_2d_y_spin.blockSignals(False)
            sa = curve.get("scan_area") or [0.0, 0.0]
            self._scan_area_min.blockSignals(True)
            self._scan_area_max.blockSignals(True)
            self._scan_area_min.setValue(float(sa[0]) if len(sa) > 0 else 0.0)
            self._scan_area_max.setValue(float(sa[1]) if len(sa) > 1 else 0.0)
            self._scan_area_min.blockSignals(False)
            self._scan_area_max.blockSignals(False)
            self._set_combo_data(self._twod_mode_cmb, curve.get("render_mode", "mesh"))
            self._twod_levels_spin.blockSignals(True)
            self._twod_levels_spin.setValue(int(curve.get("contour_levels", 8) or 8))
            self._twod_levels_spin.blockSignals(False)
            self._twod_label_contours_check.blockSignals(True)
            self._twod_label_contours_check.setChecked(bool(curve.get("label_contours", False)))
            self._twod_label_contours_check.blockSignals(False)
            self._populate_stats_for(curve)
            return

        self._set_twod_controls_visible(False)
        self._drawstyle_cmb.setEnabled(True)
        self._marker_cmb.setEnabled(True)
        self._marker_size_spin.setEnabled(True)
        self._axis_y_cmb.setEnabled(True)
        self._axis_x_cmb.setEnabled(True)
        self._bin_combine_spin.setEnabled(True)
        self._set_combo_data(self._drawstyle_cmb, curve["drawstyle"])
        self._set_combo_data(self._marker_cmb, curve["marker"])

        self._marker_size_spin.blockSignals(True)
        self._marker_size_spin.setValue(float(curve["markersize"]))
        self._marker_size_spin.blockSignals(False)

        self._set_combo_data(self._axis_y_cmb, curve.get("axis_y", "left"))
        self._set_combo_data(self._axis_x_cmb, curve.get("axis_x", "bottom"))

        self._bin_combine_spin.blockSignals(True)
        self._bin_combine_spin.setValue(int(curve.get("bin_factor", 1)))
        self._bin_combine_spin.blockSignals(False)

        # Convolution panel: σ + checkbox always available; the scan-area row
        # only makes sense for lateral / 2-D curves (per prof's note).
        is_depth = bool(curve.get("is_depth_profile", False))
        self._conv_group.setEnabled(True)
        self._conv_depth_widget.setVisible(True)
        self._scan_area_widget.setVisible(not is_depth)
        self._conv_check.blockSignals(True)
        self._conv_check.setChecked(bool(curve.get("conv_enabled", False)))
        self._conv_check.blockSignals(False)
        self._sigma_spin.blockSignals(True)
        self._sigma_spin.setValue(float(curve.get("conv_sigma_lateral", 0.0)))
        self._sigma_spin.blockSignals(False)
        self._depth_sigma_spin.blockSignals(True)
        self._depth_sigma_spin.setValue(float(curve.get("conv_sigma_depth", 0.0)))
        self._depth_sigma_spin.blockSignals(False)
        if not is_depth:
            sa = curve.get("scan_area") or [0.0, 0.0]
            self._scan_area_min.blockSignals(True)
            self._scan_area_max.blockSignals(True)
            self._scan_area_min.setValue(float(sa[0]) if len(sa) > 0 else 0.0)
            self._scan_area_max.setValue(float(sa[1]) if len(sa) > 1 else 0.0)
            self._scan_area_min.blockSignals(False)
            self._scan_area_max.blockSignals(False)

        self._populate_stats_for(curve)

    @staticmethod
    def _set_combo_data(combo: QComboBox, value: Any) -> None:
        combo.blockSignals(True)
        idx = combo.findData(value)
        if idx < 0:
            idx = 0
        combo.setCurrentIndex(idx)
        combo.blockSignals(False)

    def _update_color_button(self, color: str) -> None:
        qc = QColor(color)
        if not qc.isValid():
            qc = QColor("#888888")
        self._color_btn.setStyleSheet(
            f"background-color: {qc.name()}; border: 1px solid #555;"
        )
        self._color_text.setText(qc.name())

    def _next_cycle_color(self) -> str:
        c = _COLORS[self._color_cycle_idx % len(_COLORS)]
        self._color_cycle_idx += 1
        return c

    def _set_twod_controls_visible(self, visible: bool) -> None:
        self._twod_render_label.setVisible(visible)
        self._twod_mode_cmb.setVisible(visible)
        self._twod_contour_label.setVisible(visible)
        self._twod_contour_holder.setVisible(visible)
        self._sigma_2d_x_widget.setVisible(visible)
        self._sigma_2d_y_widget.setVisible(visible)
        self._sigma_row_widget.setVisible(not visible)
        self._depth_sigma_row_widget.setVisible(not visible)

    @staticmethod
    def _curve_sigma_key(curve: Dict[str, Any]) -> str:
        return "conv_sigma_depth" if bool(curve.get("is_depth_profile", False)) else "conv_sigma_lateral"

    def _effective_curve_sigma(self, curve: Dict[str, Any]) -> float:
        key = self._curve_sigma_key(curve)
        return float(curve.get(key, curve.get("conv_sigma", 0.0)) or 0.0)

    def _curve_scan_area(self, curve: Dict[str, Any]) -> Optional[List[float]]:
        if bool(curve.get("is_depth_profile", False)):
            return None
        return curve.get("scan_area")

    # -----------------------------------------------------------------
    # Reference-image availability
    #
    # A reference image overlay only makes physical sense where the axes
    # carry a spatial meaning (2-D maps, depth/lateral profiles). For pure
    # distribution histograms (energy / angle) the axis scaling is arbitrary
    # and an underlaid image is misleading, so the feature is hidden there
    # (per prof's note).
    # -----------------------------------------------------------------
    @staticmethod
    def _curve_is_distribution_histogram(curve: Dict[str, Any]) -> bool:
        if curve.get("is_2d"):
            return False
        label = str(curve.get("x_label", "")).strip().lower()
        return any(
            tok in label
            for tok in ("energy", "angle", "(ev", "(kev", "(mev", "(deg", "(rad")
        )

    def _curve_supports_reference(self, curve: Dict[str, Any]) -> bool:
        return not self._curve_is_distribution_histogram(curve)

    def _reference_supported(self) -> bool:
        visible = [c for c in self._curves if c.get("visible", True)]
        if not visible:
            # Nothing plotted yet — keep the panel available (neutral state).
            return True
        return any(self._curve_supports_reference(c) for c in visible)

    def _update_reference_availability(self) -> None:
        if not hasattr(self, "_ref_group"):
            return
        self._ref_group.setVisible(self._reference_supported())

    def _apply_curve_convolution_1d(
        self, curve: Dict[str, Any], x: np.ndarray, y: np.ndarray
    ) -> np.ndarray:
        if not curve.get("conv_enabled"):
            return y
        return self._gauss_convolve(
            x,
            y,
            self._effective_curve_sigma(curve),
            scan_area=self._curve_scan_area(curve),
        )

    @staticmethod
    def _gauss_convolve_2d(
        x: np.ndarray,
        y: np.ndarray,
        data: np.ndarray,
        sigma_x: float,
        sigma_y: float,
        scan_area: Optional[List[float]] = None,
    ) -> np.ndarray:
        if data.ndim != 2 or x.size < 2 or y.size < 2:
            return data
        if sigma_x <= 0.0 and sigma_y <= 0.0:
            return data
        try:
            from scipy.ndimage import gaussian_filter1d, convolve1d
            from scipy.special import erf
        except ImportError:
            return data

        out = np.asarray(data, dtype=float)
        dx = float(np.mean(np.diff(x))) if x.size > 1 else 0.0
        dy = float(np.mean(np.diff(y))) if y.size > 1 else 0.0
        if dx <= 0.0 or dy <= 0.0:
            return out

        if sigma_x > 0.0:
            out = gaussian_filter1d(out, sigma=sigma_x / dx, axis=0, mode="constant")

        if sigma_y > 0.0:
            a = float(scan_area[0]) if scan_area and len(scan_area) > 0 else 0.0
            b = float(scan_area[1]) if scan_area and len(scan_area) > 1 else 0.0
            if abs(b - a) < 1e-12:
                out = gaussian_filter1d(out, sigma=sigma_y / dy, axis=1, mode="constant")
            else:
                sqrt2s = np.sqrt(2.0) * sigma_y
                reach = 5.0 * sigma_y + max(abs(a), abs(b))
                max_half = max(3, (y.size - 1) // 2)
                half = min(max(3, int(np.ceil(reach / dy))), max_half)
                d = np.arange(-half, half + 1, dtype=float) * dy
                kernel = 0.5 * (erf((d - a) / sqrt2s) - erf((d - b) / sqrt2s)) * dy
                out = convolve1d(out, kernel, axis=1, mode="constant")
        return out

    def _current_2d_data(self, curve: Dict[str, Any]) -> np.ndarray:
        data = np.asarray(curve.get("z_values_2d", []), dtype=float)
        if not curve.get("conv_enabled"):
            return data
        return self._gauss_convolve_2d(
            np.asarray(curve.get("x_values_2d", []), dtype=float),
            np.asarray(curve.get("y_values_2d", []), dtype=float),
            data,
            float(curve.get("conv_sigma_depth", 0.0) or 0.0),
            float(curve.get("conv_sigma_lateral", 0.0) or 0.0),
            scan_area=curve.get("scan_area"),
        )

    # =====================================================================
    # Curve action slots
    # =====================================================================

    def _remove_selected_curve(self) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        cid = int(curve["id"])
        self._curves = [c for c in self._curves if int(c["id"]) != cid]
        self._refresh_curve_list()
        self._render_plot()

    def _clear_all_curves(self) -> None:
        if not self._curves:
            return
        self._curves = []
        self._refresh_curve_list()
        self._render_plot()

    # =====================================================================
    # Single-plot configuration save / load
    # =====================================================================

    # Style/state keys persisted per curve (data arrays and callables excluded).
    _CONFIG_STYLE_KEYS = (
        "label", "visible", "color", "drawstyle", "linestyle", "linewidth",
        "marker", "markersize", "alpha", "zorder", "axis_x", "axis_y",
        "bin_factor", "conv_enabled", "conv_sigma", "conv_sigma_depth",
        "conv_sigma_lateral", "scan_area", "is_depth_profile", "render_mode",
        "contour_levels", "label_contours", "x_label", "y_label",
        "colorbar_label",
    )

    def _save_plot_config(self) -> None:
        if not self._curves:
            QMessageBox.information(self, "Save Plot Config",
                                    "There are no curves to save.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Single Plot Configuration", get_last_used_directory(),
            "Single Plot Config (*.spc.json);;JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".spc.json"
        config_dir = os.path.dirname(os.path.abspath(path))

        curves_out: List[Dict[str, Any]] = []
        skipped: List[str] = []
        for curve in self._curves:
            entry: Dict[str, Any] = {k: curve[k] for k in self._CONFIG_STYLE_KEYS if k in curve}
            entry["is_2d"] = bool(curve.get("is_2d"))
            entry["source_id"] = str(curve.get("source_id", ""))
            src = str(curve.get("source_path") or "")
            if src and os.path.isdir(src):
                # Reference the data source by a path relative to the config file.
                try:
                    entry["source_path_rel"] = os.path.relpath(src, config_dir)
                except ValueError:
                    entry["source_path_rel"] = src
            elif curve.get("is_2d"):
                # 2-D curves without an external directory (e.g. convolved
                # overlays) cannot be re-read from disk and are skipped.
                skipped.append(str(curve.get("label", "?")))
                continue
            else:
                # No external source (derived/convolved 1-D curve): embed data.
                entry["data"] = {
                    "x": np.asarray(curve.get("x", []), dtype=float).tolist(),
                    "y": np.asarray(curve.get("y_original", []), dtype=float).tolist(),
                }
            curves_out.append(entry)

        payload = {
            "format": "OpenSRIM.singleplot",
            "version": 1,
            "axes": dict(self._axes_state),
            "font_size": float(self._font_size),
            "beam_from_top": bool(getattr(self, "_beam_from_top", True)),
            "curves": curves_out,
        }
        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
        except OSError as exc:
            QMessageBox.warning(self, "Save Plot Config",
                                f"Unable to save configuration:\n{exc}")
            return
        remember_last_used_path(path)
        msg = f"Saved {len(curves_out)} curve(s) to {os.path.basename(path)}."
        if skipped:
            msg += "\n\nSkipped (no re-loadable data source):\n- " + "\n- ".join(skipped)
        QMessageBox.information(self, "Save Plot Config", msg)

    def _load_plot_config(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Single Plot Configuration", get_last_used_directory(),
            "Single Plot Config (*.spc.json);;JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            QMessageBox.warning(self, "Load Plot Config",
                                f"Unable to read configuration:\n{exc}")
            return
        if not isinstance(payload, dict) or payload.get("format") != "OpenSRIM.singleplot":
            QMessageBox.warning(self, "Load Plot Config",
                                "This file is not a single-plot configuration.")
            return

        config_dir = os.path.dirname(os.path.abspath(path))
        try:
            from ui.pages.mcresults_page import _build_plots_from_directory
        except ImportError:
            from OpenSRIM.ui.pages.mcresults_page import _build_plots_from_directory  # type: ignore

        self._curves = []
        self._refresh_curve_list()

        dir_cache: Dict[str, Any] = {}
        missing: List[str] = []

        for entry in payload.get("curves", []):
            if not isinstance(entry, dict):
                continue
            is_2d = bool(entry.get("is_2d"))
            source_id = str(entry.get("source_id", ""))
            label = str(entry.get("label", "")) or source_id
            dataset = None

            rel = entry.get("source_path_rel")
            if rel:
                abs_dir = os.path.normpath(os.path.join(config_dir, rel))
                if abs_dir not in dir_cache:
                    if os.path.isdir(abs_dir):
                        try:
                            dir_cache[abs_dir] = _build_plots_from_directory(abs_dir)[0]
                        except Exception:
                            dir_cache[abs_dir] = None
                    else:
                        dir_cache[abs_dir] = None
                plots = dir_cache[abs_dir]
                if not plots:
                    missing.append(f"{label}  ←  {rel}")
                    continue
                dataset = self._dataset_from_plots(plots, source_id, is_2d, abs_dir)
                if dataset is None:
                    missing.append(f"{label}  (series not found in {rel})")
                    continue
            elif "data" in entry and not is_2d:
                d = entry["data"]
                dataset = {
                    "id": source_id or "loaded",
                    "name": label,
                    "x": np.asarray(d.get("x", []), dtype=float),
                    "y": np.asarray(d.get("y", []), dtype=float),
                    "x_label": str(entry.get("x_label", "")),
                    "y_label": str(entry.get("y_label", "")),
                    "is_depth_profile": bool(entry.get("is_depth_profile", True)),
                }
            else:
                missing.append(label)
                continue

            if not self.add_dataset_as_curve(dataset, label=label):
                missing.append(label)
                continue
            self._apply_curve_style(self._curves[-1], entry)

        # Restore axes / global view state.
        axes = payload.get("axes")
        if isinstance(axes, dict):
            self._axes_state.update(axes)
            self._sync_axes_inputs()
        try:
            self.set_font_size(float(payload.get("font_size", self._font_size)))
        except (TypeError, ValueError):
            pass
        if "beam_from_top" in payload:
            self.set_beam_from_top(bool(payload.get("beam_from_top")))

        self._refresh_curve_list()
        self._render_plot()

        if missing:
            QMessageBox.warning(
                self, "Load Plot Config",
                "The following curves could not be loaded because their data "
                "source is missing:\n\n- " + "\n- ".join(missing),
            )

    def _dataset_from_plots(self, plots: Dict[str, Any], source_id: str,
                            is_2d: bool, path: str) -> Optional[Dict[str, Any]]:
        """Reconstruct a dataset entry for *source_id* from freshly-read plots."""
        if is_2d:
            info = plots.get(source_id)
            if not isinstance(info, dict):
                return None
            return {
                "id": source_id,
                "name": str(info.get("name", source_id)),
                "path": path,
                "is_2d": True,
                "plot_info": dict(info),
            }
        plot_id, sep, idx_s = source_id.rpartition("__")
        if not sep:
            return None
        try:
            i = int(idx_s)
        except ValueError:
            return None
        info = plots.get(plot_id)
        if not isinstance(info, dict):
            return None
        cols = info.get("columns")
        x_vals = info.get("x_values")
        if cols is None or x_vals is None or i < 0 or i >= len(cols):
            return None
        labels = list(info.get("series_labels", []) or [])
        series_label = labels[i] if i < len(labels) else f"Series {i}"
        return {
            "id": source_id,
            "name": f"{info.get('name', plot_id)} — {series_label}",
            "x": np.asarray(x_vals, dtype=float),
            "y": np.asarray(cols[i], dtype=float),
            "x_label": str(info.get("x_label", "")),
            "y_label": str(info.get("y_label", "")),
            "path": path,
            "is_depth_profile": bool(info.get("is_depth_profile", False)),
        }

    def _apply_curve_style(self, curve: Dict[str, Any], entry: Dict[str, Any]) -> None:
        """Apply persisted style/state values from *entry* onto *curve*."""
        for key in self._CONFIG_STYLE_KEYS:
            if key in entry:
                curve[key] = entry[key]

    def _sync_axes_inputs(self) -> None:
        """Push the current ``_axes_state`` into the axes input widgets."""
        text_widgets = {
            "title": getattr(self, "_title_edit", None),
            "xlabel": getattr(self, "_xlabel_edit", None),
            "xlabel_top": getattr(self, "_xlabel_top_edit", None),
            "ylabel": getattr(self, "_ylabel_edit", None),
            "ylabel_right": getattr(self, "_ylabel_right_edit", None),
        }
        for key, widget in text_widgets.items():
            if widget is not None:
                widget.blockSignals(True)
                widget.setText(str(self._axes_state.get(key, "")))
                widget.blockSignals(False)
        flag_widgets = {
            "log_x": getattr(self, "_log_x_check", None),
            "log_y": getattr(self, "_log_y_check", None),
            "log_x_top": getattr(self, "_log_x_top_check", None),
            "log_y_right": getattr(self, "_log_y_right_check", None),
            "grid": getattr(self, "_grid_check", None),
            "legend": getattr(self, "_legend_check", None),
        }
        for key, widget in flag_widgets.items():
            if widget is not None:
                widget.blockSignals(True)
                widget.setChecked(bool(self._axes_state.get(key, False)))
                widget.blockSignals(False)
        loc_widget = getattr(self, "_legend_loc_cmb", None)
        if loc_widget is not None:
            idx = loc_widget.findText(str(self._axes_state.get("legend_loc", "best")))
            if idx >= 0:
                loc_widget.blockSignals(True)
                loc_widget.setCurrentIndex(idx)
                loc_widget.blockSignals(False)

    def _add_convolution_as_dataset(self) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            QMessageBox.information(
                self,
                "Add Convolution",
                "Load or select a curve first.",
            )
            return

        if curve.get("is_2d"):
            sigma_x = float(curve.get("conv_sigma_depth", 0.0) or 0.0)
            sigma_y = float(curve.get("conv_sigma_lateral", 0.0) or 0.0)
            if sigma_x <= 0.0 and sigma_y <= 0.0:
                QMessageBox.information(
                    self,
                    "Add Convolution",
                    "Set a non-zero 2D depth or lateral resolution first.",
                )
                return
            new_curve = dict(curve)
            new_curve["id"] = self._next_curve_id
            new_curve["label"] = f"{curve['label']} (convolved)"
            new_curve["z_values_2d"] = self._gauss_convolve_2d(
                np.asarray(curve.get("x_values_2d", []), dtype=float),
                np.asarray(curve.get("y_values_2d", []), dtype=float),
                np.asarray(curve.get("z_values_2d", []), dtype=float),
                sigma_x,
                sigma_y,
                scan_area=curve.get("scan_area"),
            )
            new_curve["conv_enabled"] = False
            new_curve["conv_sigma"] = 0.0
            new_curve["conv_sigma_depth"] = 0.0
            new_curve["conv_sigma_lateral"] = 0.0
            # The convolved copy inherits the source's render mode, so a
            # colormesh source yields a colormesh result (switchable per curve
            # via the "2D render" dropdown). Toggle the source's visibility to
            # compare the two meshes, or set one to contour to overlay them.
            new_curve["render_mode"] = curve.get("render_mode") or "mesh"
            new_curve["color"] = self._next_cycle_color()
            new_curve["zorder"] = int(curve.get("zorder", 2)) + 1
            self._next_curve_id += 1
            self._curves.append(new_curve)
            self._refresh_curve_list(select_last=True)
            self._render_plot_keep_view()
            return

        sigma = self._effective_curve_sigma(curve)
        if sigma <= 0.0:
            QMessageBox.information(
                self,
                "Add Convolution",
                "Set a non-zero depth or lateral resolution first.",
            )
            return

        x = np.asarray(curve.get("x", []), dtype=float)
        y = np.asarray(curve.get("y_original", []), dtype=float)
        if x.size == 0 or y.size != x.size:
            return
        bf = int(curve.get("bin_factor", 1) or 1)
        if bf > 1:
            x, y = self._rebin_curve(x, y, bf)
        y_conv = self._gauss_convolve(
            x,
            y,
            sigma,
            scan_area=self._curve_scan_area(curve),
        )
        new_curve = dict(curve)
        new_curve["id"] = self._next_curve_id
        new_curve["label"] = f"{curve['label']} (convolved)"
        new_curve["x"] = x.copy()
        new_curve["y_original"] = y_conv.copy()
        new_curve["drawstyle"] = "default"
        new_curve["linestyle"] = "--"
        new_curve["conv_enabled"] = False
        new_curve["conv_sigma"] = 0.0
        new_curve["conv_sigma_depth"] = 0.0
        new_curve["conv_sigma_lateral"] = 0.0
        new_curve["bin_factor"] = 1
        new_curve["color"] = self._next_cycle_color()
        new_curve["zorder"] = int(curve.get("zorder", 2)) + 1
        self._next_curve_id += 1
        self._curves.append(new_curve)
        self._refresh_curve_list(select_last=True)
        self._render_plot_keep_view()

    @staticmethod
    def _fmt_sigma(value: float) -> str:
        return f"{float(value):g}"

    @staticmethod
    def _parse_sigma_list(text: str) -> List[float]:
        """Parse a comma/space separated list of positive σ values, keeping
        order and dropping duplicates / non-positive / unparsable tokens."""
        cleaned = str(text).replace(";", " ").replace(",", " ")
        values: List[float] = []
        seen: set = set()
        for tok in cleaned.split():
            try:
                v = float(tok)
            except ValueError:
                continue
            if v <= 0.0:
                continue
            key = round(v, 9)
            if key in seen:
                continue
            seen.add(key)
            values.append(v)
        return values

    def _add_multi_sigma_comparison(self) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            QMessageBox.information(
                self,
                "σ Comparison",
                "Load or select a curve first.",
            )
            return
        sigmas = self._parse_sigma_list(self._multi_sigma_edit.text())
        if not sigmas:
            QMessageBox.information(
                self,
                "σ Comparison",
                "Enter one or more positive σ values, e.g. 10, 25, 50.",
            )
            return
        if curve.get("is_2d"):
            added = self._add_multi_sigma_2d(curve, sigmas)
        else:
            added = self._add_multi_sigma_1d(curve, sigmas)
        if added:
            self._refresh_curve_list(select_last=True)
            self._render_plot_keep_view()

    def _add_multi_sigma_1d(self, curve: Dict[str, Any], sigmas: List[float]) -> int:
        x = np.asarray(curve.get("x", []), dtype=float)
        y = np.asarray(curve.get("y_original", []), dtype=float)
        if x.size == 0 or y.size != x.size:
            return 0
        bf = int(curve.get("bin_factor", 1) or 1)
        if bf > 1:
            x, y = self._rebin_curve(x, y, bf)
        scan_area = self._curve_scan_area(curve)
        base_label = str(curve["label"])
        base_zorder = int(curve.get("zorder", 2))
        added = 0
        for sigma in sigmas:
            y_conv = self._gauss_convolve(x, y, sigma, scan_area=scan_area)
            new_curve = dict(curve)
            new_curve["id"] = self._next_curve_id
            new_curve["label"] = f"{base_label} (σ={self._fmt_sigma(sigma)} Å)"
            new_curve["x"] = x.copy()
            new_curve["y_original"] = y_conv.copy()
            new_curve["drawstyle"] = "default"
            new_curve["linestyle"] = "--"
            new_curve["conv_enabled"] = False
            new_curve["conv_sigma"] = 0.0
            new_curve["conv_sigma_depth"] = 0.0
            new_curve["conv_sigma_lateral"] = 0.0
            new_curve["bin_factor"] = 1
            new_curve["color"] = self._next_cycle_color()
            new_curve["zorder"] = base_zorder + 1 + added
            self._next_curve_id += 1
            self._curves.append(new_curve)
            added += 1
        return added

    def _add_multi_sigma_2d(self, curve: Dict[str, Any], sigmas: List[float]) -> int:
        x = np.asarray(curve.get("x_values_2d", []), dtype=float)
        y = np.asarray(curve.get("y_values_2d", []), dtype=float)
        z = np.asarray(curve.get("z_values_2d", []), dtype=float)
        if x.size < 2 or y.size < 2 or z.size == 0:
            return 0
        scan_area = curve.get("scan_area")
        base_label = str(curve["label"])
        base_zorder = int(curve.get("zorder", 2))
        # The source keeps its render mode (typically a colormesh, so the
        # colour scale stays); each σ variant is overlaid as contour lines.
        added = 0
        for sigma in sigmas:
            z_conv = self._gauss_convolve_2d(x, y, z, sigma, sigma, scan_area=scan_area)
            new_curve = dict(curve)
            new_curve["id"] = self._next_curve_id
            new_curve["label"] = f"{base_label} (σ={self._fmt_sigma(sigma)} Å)"
            new_curve["z_values_2d"] = z_conv
            new_curve["conv_enabled"] = False
            new_curve["conv_sigma"] = 0.0
            new_curve["conv_sigma_depth"] = 0.0
            new_curve["conv_sigma_lateral"] = 0.0
            new_curve["render_mode"] = "contour"
            new_curve["color"] = self._next_cycle_color()
            new_curve["zorder"] = base_zorder + 1 + added
            self._next_curve_id += 1
            self._curves.append(new_curve)
            added += 1
        return added

    def _on_label_edited(self, text: str) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        text = text.strip() or curve["label"]
        curve["label"] = text
        # Sync list row
        item = self._curve_list.currentItem()
        if item is not None:
            self._curve_list.blockSignals(True)
            item.setText(text)
            self._curve_list.blockSignals(False)
        self._render_plot()

    def _pick_color(self) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        current = QColor(curve["color"])
        color = QColorDialog.getColor(current, self, "Select Curve Color")
        if not color.isValid():
            return
        curve["color"] = color.name()
        self._update_color_button(curve["color"])
        item = self._curve_list.currentItem()
        if item is not None:
            self._curve_list.blockSignals(True)
            item.setIcon(_color_icon(curve["color"]))
            self._curve_list.blockSignals(False)
        self._render_plot()

    def _on_drawstyle_changed(self, _idx: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["drawstyle"] = self._drawstyle_cmb.currentData()
        self._render_plot()

    def _on_linestyle_changed(self, _idx: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["linestyle"] = self._linestyle_cmb.currentData()
        self._render_plot()

    def _on_linewidth_changed(self, v: float) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["linewidth"] = float(v)
        self._render_plot()

    def _on_marker_changed(self, _idx: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["marker"] = self._marker_cmb.currentData()
        self._render_plot()

    def _on_markersize_changed(self, v: float) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["markersize"] = float(v)
        self._render_plot()

    def _on_alpha_changed(self, v: float) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["alpha"] = float(v)
        self._render_plot()

    def _on_zorder_changed(self, v: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["zorder"] = int(v)
        self._render_plot()

    def _on_axis_y_changed(self, _idx: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["axis_y"] = self._axis_y_cmb.currentData() or "left"
        self._update_secondary_axes_visibility()
        self._render_plot()

    def _on_axis_x_changed(self, _idx: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        curve["axis_x"] = self._axis_x_cmb.currentData() or "bottom"
        self._update_secondary_axes_visibility()
        self._render_plot()

    def _ensure_curve_selected(self) -> Optional[Dict[str, Any]]:
        """Return the selected curve, auto-selecting the first 1-D curve if
        the user has not picked one yet. Returns ``None`` only when *no*
        curve has been loaded at all."""
        curve = self._selected_curve()
        if curve is not None:
            return curve
        if self._curves:
            self._curve_list.setCurrentRow(0)
            return self._selected_curve()
        return None

    def _on_conv_toggled(self, checked: bool) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            # No curves loaded yet — undo the toggle so the UI doesn't lie
            # about the state, and nudge the user via a transient tooltip.
            self._conv_check.blockSignals(True)
            self._conv_check.setChecked(False)
            self._conv_check.blockSignals(False)
            QToolTip.showText(
                self._conv_check.mapToGlobal(self._conv_check.rect().bottomLeft()),
                "Load a curve first: double-click a plot tile on the MC Results tab.",
                self._conv_check,
            )
            return
        curve["conv_enabled"] = bool(checked)
        # When the user first turns convolution on with σ still at zero, pick
        # a sensible non-zero default so the effect is immediately visible.
        # Heuristic: ~3 % of the curve's x-range.
        needs_default = False
        if curve.get("is_2d"):
            needs_default = (
                float(curve.get("conv_sigma_depth", 0.0) or 0.0) <= 0.0
                and float(curve.get("conv_sigma_lateral", 0.0) or 0.0) <= 0.0
            )
        else:
            needs_default = self._effective_curve_sigma(curve) <= 0.0
        if checked and needs_default:
            if curve.get("is_2d"):
                x = np.asarray(curve.get("x_values_2d", []), dtype=float)
                y = np.asarray(curve.get("y_values_2d", []), dtype=float)
            else:
                x = np.asarray(curve.get("x", []), dtype=float)
                y = np.array([], dtype=float)
            if curve.get("is_2d"):
                if x.size >= 2 and float(curve.get("conv_sigma_depth", 0.0) or 0.0) <= 0.0:
                    sigma_x = max(float(x[-1] - x[0]) * 0.03, 0.0)
                    curve["conv_sigma_depth"] = sigma_x
                    self._depth_sigma_spin.blockSignals(True)
                    self._depth_sigma_spin.setValue(sigma_x)
                    self._depth_sigma_spin.blockSignals(False)
                    self._sigma_2d_x_spin.blockSignals(True)
                    self._sigma_2d_x_spin.setValue(sigma_x)
                    self._sigma_2d_x_spin.blockSignals(False)
                if y.size >= 2 and float(curve.get("conv_sigma_lateral", 0.0) or 0.0) <= 0.0:
                    sigma_y = max(float(y[-1] - y[0]) * 0.03, 0.0)
                    curve["conv_sigma_lateral"] = sigma_y
                    curve["conv_sigma"] = sigma_y
                    self._sigma_spin.blockSignals(True)
                    self._sigma_spin.setValue(sigma_y)
                    self._sigma_spin.blockSignals(False)
                    self._sigma_2d_y_spin.blockSignals(True)
                    self._sigma_2d_y_spin.setValue(sigma_y)
                    self._sigma_2d_y_spin.blockSignals(False)
            elif x.size >= 2 and self._effective_curve_sigma(curve) <= 0.0:
                sigma_default = max(float(x[-1] - x[0]) * 0.03, 0.0)
                key = self._curve_sigma_key(curve)
                curve[key] = sigma_default
                curve["conv_sigma"] = sigma_default
                target_spin = self._depth_sigma_spin if key == "conv_sigma_depth" else self._sigma_spin
                target_spin.blockSignals(True)
                target_spin.setValue(sigma_default)
                target_spin.blockSignals(False)
        # Toggling convolution must not snap the view back to autoscale.
        self._render_plot_keep_view()

    def _on_sigma_changed(self, v: float) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            return
        curve["conv_sigma_lateral"] = float(v)
        curve["conv_sigma"] = float(v)
        if curve.get("conv_enabled"):
            self._render_plot_keep_view()

    def _on_depth_sigma_changed(self, v: float) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            return
        curve["conv_sigma_depth"] = float(v)
        if curve.get("is_2d"):
            self._sigma_2d_x_spin.blockSignals(True)
            self._sigma_2d_x_spin.setValue(float(v))
            self._sigma_2d_x_spin.blockSignals(False)
        if curve.get("conv_enabled"):
            self._render_plot_keep_view()

    def _on_sigma_2d_changed(self, _v: float) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            return
        curve["conv_sigma_depth"] = float(self._sigma_2d_x_spin.value())
        curve["conv_sigma_lateral"] = float(self._sigma_2d_y_spin.value())
        curve["conv_sigma"] = float(self._sigma_2d_y_spin.value())
        self._depth_sigma_spin.blockSignals(True)
        self._depth_sigma_spin.setValue(curve["conv_sigma_depth"])
        self._depth_sigma_spin.blockSignals(False)
        self._sigma_spin.blockSignals(True)
        self._sigma_spin.setValue(curve["conv_sigma_lateral"])
        self._sigma_spin.blockSignals(False)
        if curve.get("conv_enabled"):
            self._render_plot_keep_view()

    def _on_scan_area_changed(self, _v: float) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            return
        curve["scan_area"] = [
            float(self._scan_area_min.value()),
            float(self._scan_area_max.value()),
        ]
        if curve.get("conv_enabled"):
            self._render_plot_keep_view()

    def _on_twod_render_mode_changed(self, _idx: int) -> None:
        curve = self._selected_curve()
        if curve is None or not curve.get("is_2d"):
            return
        curve["render_mode"] = self._twod_mode_cmb.currentData() or "mesh"
        self._render_plot()

    def _on_twod_levels_changed(self, value: int) -> None:
        curve = self._selected_curve()
        if curve is None or not curve.get("is_2d"):
            return
        curve["contour_levels"] = int(value)
        self._render_plot()

    def _on_twod_contour_labels_toggled(self, checked: bool) -> None:
        curve = self._selected_curve()
        if curve is None or not curve.get("is_2d"):
            return
        curve["label_contours"] = bool(checked)
        self._render_plot()

    def _on_bin_combine_changed(self, value: int) -> None:
        curve = self._selected_curve()
        if curve is None:
            return
        new = max(1, min(6, int(value)))
        if int(curve.get("bin_factor", 1)) == new:
            return
        curve["bin_factor"] = new
        self._render_plot()
        self._populate_stats_for(curve)

    # =====================================================================
    # Reference-image slots
    # =====================================================================

    def _pick_reference_image(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Reference Image", get_last_used_directory(),
            "Images (*.png *.jpg *.jpeg);;All files (*)",
        )
        if not path:
            return
        remember_last_used_path(path)
        try:
            data = mpimg.imread(path)
        except Exception as exc:  # pragma: no cover - depends on file
            QMessageBox.warning(
                self, "Reference Image", f"Could not load image:\n{exc}"
            )
            return
        self._ref_image_path = path
        self._ref_image_data = data
        self._ref_clear_btn.setEnabled(True)
        self._ref_status.setText(os.path.basename(path))
        self._render_plot()

    def _clear_reference_image(self) -> None:
        self._ref_image_path = None
        self._ref_image_data = None
        self._ref_clear_btn.setEnabled(False)
        self._ref_status.setText("No image loaded.")
        self._render_plot()

    def _on_ref_alpha_changed(self, v: float) -> None:
        self._ref_image_alpha = float(v)
        if self._ref_image_data is not None:
            self._render_plot()

    def _on_canvas_leave(self, _event) -> None:
        self._cursor_status.setText("x: -, y: -")

    def _on_canvas_motion(self, event) -> None:
        if event is None or event.inaxes is None or event.xdata is None or event.ydata is None:
            self._cursor_status.setText("x: -, y: -")
            return

        x_val = float(event.xdata)
        y_val = float(event.ydata)
        twod_curve = self._selected_curve()
        if twod_curve is None or not twod_curve.get("is_2d") or not twod_curve.get("visible", True):
            twod_curve = next(
                (c for c in self._curves if c.get("is_2d") and c.get("visible", True)),
                None,
            )
        if twod_curve is not None:
            x_grid = np.asarray(twod_curve.get("x_values_2d", []), dtype=float)
            y_grid = np.asarray(twod_curve.get("y_values_2d", []), dtype=float)
            data = self._current_2d_data(twod_curve)
            if x_grid.size and y_grid.size and data.ndim == 2 and data.shape == (x_grid.size, y_grid.size):
                if getattr(self, "_beam_from_top", True):
                    ix = int(np.argmin(np.abs(x_grid - y_val)))
                    iy = int(np.argmin(np.abs(y_grid - x_val)))
                else:
                    ix = int(np.argmin(np.abs(x_grid - x_val)))
                    iy = int(np.argmin(np.abs(y_grid - y_val)))
                z_val = float(data[ix, iy])
                self._cursor_status.setText(
                    f"x: {x_val:.3f}, y: {y_val:.3f}, f(x,y): {z_val:.6g}"
                )
                return

        self._cursor_status.setText(f"x: {x_val:.3f}, y: {y_val:.3f}")

    # =====================================================================
    # Axes-level slots
    # =====================================================================

    def _on_axes_edit(self, key: str, text: str) -> None:
        self._axes_state[key] = text
        self._axes_user_overrides.add(key)
        self._render_plot()

    def _on_axes_flag(self, key: str, val: bool) -> None:
        self._axes_state[key] = bool(val)
        self._axes_user_overrides.add(key)
        self._render_plot()

    def _has_curve_on(self, axis: str, value: str) -> bool:
        return any(
            c.get(axis, "left" if axis == "axis_y" else "bottom") == value
            for c in self._curves
        )

    def _update_secondary_axes_visibility(self) -> None:
        """Show secondary-axis label/scale controls only when a curve uses them."""
        has_right_y = self._has_curve_on("axis_y", "right")
        has_top_x = self._has_curve_on("axis_x", "top")

        for w in (self._ylabel_right_label, self._ylabel_right_edit):
            w.setVisible(has_right_y)
        for w in (self._xlabel_top_label, self._xlabel_top_edit):
            w.setVisible(has_top_x)
        for w in (self._scale2_label, self._scale2_holder):
            w.setVisible(has_right_y or has_top_x)
        # Within the scale row, hide individual checkboxes that aren't usable.
        self._log_y_right_check.setVisible(has_right_y)
        self._log_x_top_check.setVisible(has_top_x)

    # =====================================================================
    # Statistics
    # =====================================================================

    def _populate_stats_for(self, curve: Dict[str, Any]) -> None:
        if curve.get("is_2d"):
            data = self._current_2d_data(curve)
            x = np.asarray(curve.get("x_values_2d", []), dtype=float)
            y = np.asarray(curve.get("y_values_2d", []), dtype=float)
            rows: List[tuple] = [
                ("Source", str(curve.get("source_name", ""))),
                ("Series", str(curve.get("label", ""))),
                ("Grid", f"{x.size} × {y.size}"),
                ("Sum", f"{float(np.sum(data)):.3f}"),
                ("Peak value", f"{float(np.max(data)):.3f}" if data.size else "0.000"),
            ]
            if data.size:
                peak_idx = np.unravel_index(int(np.argmax(data)), data.shape)
                rows.append(("Peak X", f"{float(x[peak_idx[0]]):.3f}" if x.size else "0.000"))
                rows.append(("Peak Y", f"{float(y[peak_idx[1]]):.3f}" if y.size else "0.000"))
            self._stats_table.setRowCount(len(rows))
            for r, (name, value) in enumerate(rows):
                ni = QTableWidgetItem(str(name))
                ni.setFlags(ni.flags() & ~Qt.ItemFlag.ItemIsEditable)
                vi = QTableWidgetItem(str(value))
                vi.setFlags(vi.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._stats_table.setItem(r, 0, ni)
                self._stats_table.setItem(r, 1, vi)
            return

        x = np.asarray(curve["x"], dtype=float)
        y = np.asarray(curve["y_original"], dtype=float)
        bf = int(curve.get("bin_factor", 1) or 1)
        if bf > 1:
            x, y = self._rebin_curve(x, y, bf)
        y = self._apply_curve_convolution_1d(curve, x, y)

        rows: List[tuple] = []
        if y.size:
            total = float(np.sum(y))
            rows.append(("Source", str(curve["source_name"])))
            rows.append(("Series", str(curve["label"])))
            rows.append(("N samples", f"{y.size}"))
            rows.append(("Sum", f"{total:.3f}"))
            if total > 0:
                mean = float(np.sum(x * y) / total)
                var = float(np.sum(((x - mean) ** 2) * y) / total)
                std = float(np.sqrt(max(0.0, var)))
                rows.append(("Mean (weighted)", f"{mean:.3f}"))
                rows.append(("Std  (weighted)", f"{std:.3f}"))
            rows.append(("Peak value", f"{float(np.max(y)):.3f}"))
            peak_idx = int(np.argmax(y))
            rows.append(("Peak X",     f"{float(x[peak_idx]):.3f}"))
            rows.append(("X range", f"{float(np.min(x)):.3f} … {float(np.max(x)):.3f}"))

        self._stats_table.setRowCount(len(rows))
        for r, (name, value) in enumerate(rows):
            ni = QTableWidgetItem(str(name))
            ni.setFlags(ni.flags() & ~Qt.ItemFlag.ItemIsEditable)
            vi = QTableWidgetItem(str(value))
            vi.setFlags(vi.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._stats_table.setItem(r, 0, ni)
            self._stats_table.setItem(r, 1, vi)

    # =====================================================================
    # Render
    # =====================================================================

    @staticmethod
    def _rebin_curve(x: np.ndarray, y: np.ndarray, factor: int) -> tuple:
        """Combine adjacent bins by *factor* (sum y, average x). Trailing
        partial groups are dropped. ``factor <= 1`` returns the inputs."""
        f = int(factor)
        if f <= 1 or x.size < 2 or y.size != x.size:
            return x, y
        n = (x.size // f) * f
        if n == 0:
            return x, y
        x_new = x[:n].reshape(-1, f).mean(axis=1)
        y_new = y[:n].reshape(-1, f).sum(axis=1)
        return x_new, y_new

    @staticmethod
    def _gauss_convolve(x: np.ndarray, y: np.ndarray, sigma: float,
                        scan_area: Optional[List[float]] = None) -> np.ndarray:
        """Convolve y(x) with a Gaussian, optionally smeared over a scan area.

        For the depth-profile case (or any call without a scan area) this is
        the classical Gaussian smearing: kernel ∝ exp(−x²/2σ²), normalized so
        the kernel sum is 1 (i.e. total counts are preserved).

        For lateral / 2-D profiles with a scan area [a, b], it uses the
        bin-integral form from the prof's slides:

            G(d) = ½ · [erf((d − a)/(√2 σ)) − erf((d − b)/(√2 σ))]
            h(x_i) = Σⱼ f(x_j) · G(x_i − x_j) · Δx

        i.e. the Gaussian convolved with the rectangular scan-area indicator,
        evaluated on the bin grid.
        """
        if sigma <= 0 or x.size < 2 or y.size != x.size:
            return y
        dx = float(np.mean(np.diff(x))) if x.size > 1 else 1.0
        if dx <= 0:
            return y

        a = float(scan_area[0]) if scan_area and len(scan_area) > 0 else 0.0
        b = float(scan_area[1]) if scan_area and len(scan_area) > 1 else 0.0
        sqrt2sigma = float(np.sqrt(2.0) * sigma)
        # Kernel extends ±5σ beyond the scan area on either side. Cap the
        # kernel size below the data length — np.convolve(mode="same") would
        # otherwise return an array in kernel length, not in data length.
        reach = 5.0 * sigma + max(abs(a), abs(b))
        max_half = max(3, (x.size - 1) // 2)
        half = min(max(3, int(np.ceil(reach / dx))), max_half)
        d = np.arange(-half, half + 1, dtype=float) * dx

        if abs(b - a) < 1e-12:
            # No (or degenerate) scan area: plain Gaussian smearing, normalized.
            kernel = np.exp(-0.5 * (d / sigma) ** 2)
            s = float(kernel.sum())
            if s <= 0:
                return y
            kernel /= s
        else:
            try:
                from scipy.special import erf
            except ImportError:
                # NumPy ≥ 2.0 ships an erf in numpy.polynomial; fall back to
                # math.erf vectorized via numpy.vectorize as a last resort.
                import math
                _vec_erf = np.vectorize(math.erf, otypes=[float])
                erf = _vec_erf  # type: ignore[assignment]
            kernel = 0.5 * (erf((d - a) / sqrt2sigma) - erf((d - b) / sqrt2sigma)) * dx

        return np.convolve(y, kernel, mode="same")

    def _render_plot(self) -> None:
        # Keep the Reference Image panel hidden where it makes no sense
        # (energy / angle histograms) and visible for spatial / 2-D curves.
        self._update_reference_availability()

        # Optionally preserve the current view (zoom/pan section) across this
        # re-render. Set by operations that should not snap back to autoscale,
        # e.g. applying a convolution after the user changed the plot section.
        preserve_view = getattr(self, "_preserve_view_once", False)
        self._preserve_view_once = False
        self._saved_view_limits = None
        if preserve_view:
            try:
                cur_xlim = self.ax.get_xlim()
                cur_ylim = self.ax.get_ylim()
                # Ignore the matplotlib default (0,1)x(0,1) of an empty axes.
                if cur_xlim != (0.0, 1.0) or cur_ylim != (0.0, 1.0):
                    self._saved_view_limits = (cur_xlim, cur_ylim)
            except Exception:
                self._saved_view_limits = None

        # Recreate axes from scratch every render — this is the cleanest way
        # to handle twin axes appearing / disappearing as curves change.
        self.figure.clf()
        self.ax = self.figure.add_subplot(111)

        ax_BL = self.ax       # bottom-x + left-y  (primary)
        ax_BR = None          # bottom-x + right-y (twinx of primary)
        ax_TL = None          # top-x    + left-y  (twiny of primary)
        ax_TR = None          # top-x    + right-y (twinx of top-x)

        if not self._curves:
            self.ax.text(0.5, 0.5,
                         "No curves yet.\nDouble-click a plot tile in MC Results,\n"
                         "or use Load File / Load Directory.",
                         ha="center", va="center", transform=self.ax.transAxes,
                         color="#888")
            self.ax.set_axis_off()
            self.canvas.draw_idle()
            return
        self.ax.set_axis_on()

        visible_2d = [c for c in self._curves if c.get("is_2d") and c.get("visible", True)]
        if visible_2d:
            try:
                # Each 2-D curve renders by its own mode: a "mesh" curve draws
                # a filled colormesh (with colour scale), a "contour" curve
                # draws overlaid contour lines. This lets several 2-D datasets
                # coexist — e.g. one colormesh background plus convolved
                # contour overlays — instead of collapsing everything to lines.
                colorbar = None
                colorbar_added = False
                legend_lines: list = []
                if self._ref_image_data is not None:
                    base = visible_2d[0]
                    x0 = np.asarray(base.get("x_values_2d", []), dtype=float)
                    y0 = np.asarray(base.get("y_values_2d", []), dtype=float)
                    if x0.size and y0.size:
                        self.ax.imshow(
                            self._ref_image_data,
                            extent=(float(np.min(x0)), float(np.max(x0)), float(np.min(y0)), float(np.max(y0))),
                            aspect="auto",
                            zorder=0,
                            alpha=float(self._ref_image_alpha),
                        )

                for idx, curve in enumerate(visible_2d):
                    x = np.asarray(curve.get("x_values_2d", []), dtype=float)
                    y = np.asarray(curve.get("y_values_2d", []), dtype=float)
                    data = self._current_2d_data(curve)
                    if x.size == 0 or y.size == 0 or data.size == 0:
                        continue
                    if (curve.get("render_mode") or "mesh") == "contour":
                        levels_count = max(2, int(curve.get("contour_levels", 8) or 8))
                        vmax = float(np.max(data))
                        if vmax <= 0.0:
                            continue
                        levels = np.linspace(vmax / levels_count, vmax, levels_count)
                        contour_x = y if getattr(self, "_beam_from_top", True) else x
                        contour_y = x if getattr(self, "_beam_from_top", True) else y
                        contour_z = data if getattr(self, "_beam_from_top", True) else data.T
                        contour = self.ax.contour(
                            contour_x,
                            contour_y,
                            contour_z,
                            levels=levels,
                            colors=[curve.get("color", "#444444")],
                            linewidths=float(curve.get("linewidth", 1.0)),
                            linestyles=curve.get("linestyle", "-") or "-",
                            alpha=float(curve.get("alpha", 1.0)),
                            zorder=int(curve.get("zorder", 2)) + idx,
                        )
                        if curve.get("label_contours"):
                            try:
                                self.ax.clabel(contour, inline=True, fontsize=max(6, int(self._font_size - 1)))
                            except Exception:
                                pass
                        proxy_line, = self.ax.plot(
                            [], [],
                            color=curve.get("color", "#444444"),
                            linestyle=curve.get("linestyle", "-") or "-",
                            linewidth=float(curve.get("linewidth", 1.0)),
                            label=str(curve.get("label", "")),
                        )
                        legend_lines.append(proxy_line)
                    else:
                        mesh_x = y if getattr(self, "_beam_from_top", True) else x
                        mesh_y = x if getattr(self, "_beam_from_top", True) else y
                        mesh_z = data if getattr(self, "_beam_from_top", True) else data.T
                        mesh = self.ax.pcolormesh(
                            mesh_x,
                            mesh_y,
                            mesh_z,
                            shading="auto",
                            cmap="viridis",
                            alpha=float(curve.get("alpha", 1.0)),
                        )
                        # Only the first mesh owns the colour scale; further
                        # meshes still draw but don't stack extra colorbars.
                        if not colorbar_added:
                            colorbar = self.figure.colorbar(mesh, ax=self.ax)
                            if curve.get("colorbar_label"):
                                colorbar.set_label(str(curve.get("colorbar_label")))
                            colorbar_added = True
                if getattr(self, "_beam_from_top", True):
                    self.ax.invert_yaxis()
                self.ax.set_aspect("equal", adjustable="box")
                if self._axes_state.get("grid"):
                    self.ax.grid(True, alpha=0.25)
                else:
                    self.ax.grid(False)
                if legend_lines and self._axes_state.get("legend"):
                    self.ax.legend(
                        legend_lines,
                        [ln.get_label() for ln in legend_lines],
                        loc=str(self._axes_state.get("legend_loc", "best")),
                        fontsize=8,
                    )
            except Exception as exc:
                self.ax.text(
                    0.5, 0.5,
                    f"Failed to render 2D plot:\n{exc}",
                    ha="center", va="center", transform=self.ax.transAxes,
                    color="#a33",
                )
                self.ax.set_axis_off()
                self.canvas.draw_idle()
                return
            if self._axes_state.get("title"):
                self.ax.set_title(str(self._axes_state["title"]))
            if visible_2d and getattr(self, "_beam_from_top", True):
                if "xlabel" not in self._axes_user_overrides:
                    self.ax.set_xlabel(str(visible_2d[0].get("y_label", "")))
                elif self._axes_state.get("xlabel"):
                    self.ax.set_xlabel(str(self._axes_state["xlabel"]))
                if "ylabel" not in self._axes_user_overrides:
                    self.ax.set_ylabel(str(visible_2d[0].get("x_label", "")))
                elif self._axes_state.get("ylabel"):
                    self.ax.set_ylabel(str(self._axes_state["ylabel"]))
            else:
                if self._axes_state.get("xlabel"):
                    self.ax.set_xlabel(str(self._axes_state["xlabel"]))
                if self._axes_state.get("ylabel"):
                    self.ax.set_ylabel(str(self._axes_state["ylabel"]))
            _apply_axes_font_size(self.ax, self._font_size)
            try:
                self.figure.tight_layout()
            except Exception:
                pass
            self._restore_saved_view(self.ax)
            self.canvas.draw_idle()
            sel = self._selected_curve()
            if sel is not None:
                self._populate_stats_for(sel)
            return

        any_plotted = False
        # Collect (axis, line) for a combined legend across all axes
        legend_lines: list = []

        for curve in self._curves:
            if not curve.get("visible", True):
                continue
            if curve.get("is_2d"):
                continue
            x = np.asarray(curve["x"], dtype=float)
            y = np.asarray(curve["y_original"], dtype=float)
            bf = int(curve.get("bin_factor", 1) or 1)
            if bf > 1:
                x, y = self._rebin_curve(x, y, bf)
            y = self._apply_curve_convolution_1d(curve, x, y)

            axis_x = curve.get("axis_x", "bottom")
            axis_y = curve.get("axis_y", "left")
            if axis_x == "bottom" and axis_y == "left":
                target = ax_BL
            elif axis_x == "bottom" and axis_y == "right":
                if ax_BR is None:
                    ax_BR = ax_BL.twinx()
                target = ax_BR
            elif axis_x == "top" and axis_y == "left":
                if ax_TL is None:
                    ax_TL = ax_BL.twiny()
                target = ax_TL
            else:  # top + right
                if ax_TR is None:
                    # Build a twin that shares y with right-y and x with top-x.
                    if ax_BR is None:
                        ax_BR = ax_BL.twinx()
                    ax_TR = ax_BR.twiny()
                target = ax_TR

            ls = curve.get("linestyle", "-") or "None"
            marker = curve.get("marker", "None") or "None"
            kwargs = dict(
                color=curve.get("color"),
                drawstyle=curve.get("drawstyle", "default") or "default",
                linestyle=(ls if ls != "None" else "None"),
                linewidth=float(curve.get("linewidth", 1.4)),
                marker=(marker if marker != "None" else None),
                markersize=float(curve.get("markersize", 4.0)),
                alpha=float(curve.get("alpha", 1.0)),
                zorder=int(curve.get("zorder", 2)),
                label=str(curve.get("label", "")),
            )
            line, = target.plot(x, y, **kwargs)
            legend_lines.append(line)
            any_plotted = True

        # Reference image: drawn under the curves, spanning the data extent.
        # Only overlaid for spatial profiles — never under energy/angle
        # distribution histograms, where axis scaling is arbitrary.
        if self._ref_image_data is not None and any_plotted and self._reference_supported():
            try:
                xlim = ax_BL.get_xlim()
                ylim = ax_BL.get_ylim()
                ax_BL.imshow(
                    self._ref_image_data,
                    extent=(xlim[0], xlim[1], ylim[0], ylim[1]),
                    aspect="auto",
                    zorder=0,
                    alpha=float(self._ref_image_alpha),
                )
                ax_BL.set_xlim(xlim)
                ax_BL.set_ylim(ylim)
            except Exception:
                pass

        # Title goes on the primary axes
        ax_BL.set_title(str(self._axes_state.get("title", "")))

        # Axis labels
        ax_BL.set_xlabel(str(self._axes_state.get("xlabel", "")))
        ax_BL.set_ylabel(str(self._axes_state.get("ylabel", "")))
        if ax_BR is not None:
            ax_BR.set_ylabel(str(self._axes_state.get("ylabel_right", "")))
        if ax_TL is not None:
            ax_TL.set_xlabel(str(self._axes_state.get("xlabel_top", "")))

        # Log / linear scaling
        def _set_scale(ax, axis: str, log: bool) -> None:
            if ax is None:
                return
            try:
                if axis == "x":
                    ax.set_xscale("log" if log else "linear")
                else:
                    ax.set_yscale("log" if log else "linear")
            except Exception:
                pass

        log_x = bool(self._axes_state.get("log_x"))
        log_y = bool(self._axes_state.get("log_y"))
        log_x_top = bool(self._axes_state.get("log_x_top"))
        log_y_right = bool(self._axes_state.get("log_y_right"))
        _set_scale(ax_BL, "x", log_x)
        _set_scale(ax_BL, "y", log_y)
        _set_scale(ax_BR, "y", log_y_right)
        _set_scale(ax_TL, "x", log_x_top)
        if ax_TR is not None:
            _set_scale(ax_TR, "x", log_x_top)
            _set_scale(ax_TR, "y", log_y_right)

        # Grid (only on primary to avoid double grids)
        if self._axes_state.get("grid"):
            ax_BL.grid(True, alpha=0.3)
        else:
            ax_BL.grid(False)

        # Combined legend across all (twin) axes, placed on the primary axis.
        if any_plotted and self._axes_state.get("legend"):
            try:
                labels = [ln.get_label() for ln in legend_lines]
                ax_BL.legend(
                    legend_lines, labels,
                    loc=str(self._axes_state.get("legend_loc", "best")),
                    fontsize=8,
                )
            except Exception:
                pass

        # Font sizes on all axes that exist
        for ax in (ax_BL, ax_BR, ax_TL, ax_TR):
            if ax is not None:
                _apply_axes_font_size(ax, self._font_size)

        try:
            self.figure.tight_layout()
        except Exception:
            pass
        self._restore_saved_view(ax_BL)
        self.canvas.draw_idle()

        # Stats follow the currently selected curve (if any)
        sel = self._selected_curve()
        if sel is not None:
            self._populate_stats_for(sel)

    def _restore_saved_view(self, ax) -> None:
        """Re-apply the view limits captured before a preserve-view re-render."""
        limits = getattr(self, "_saved_view_limits", None)
        if not limits or ax is None:
            return
        xlim, ylim = limits
        try:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        except Exception:
            pass
        self._saved_view_limits = None

    def _render_plot_keep_view(self) -> None:
        """Re-render while preserving the current zoom/pan section."""
        self._preserve_view_once = True
        self._render_plot()


# =====================================================================
#  Sidebar  (left side — plot selector + numerical values)
# =====================================================================

class ResultsSidebar(QWidget):
    """Left column: plot checkboxes on top, numerical values table below.

    The plot list supports **drag-and-drop reordering** — the display
    order in the plot area follows the list order.
    """

    selection_changed = pyqtSignal(list)  # emits list of selected plot_ids
    mode_changed = pyqtSignal(bool)       # True = multiple plot mode
    single_plot_config_changed = pyqtSignal(list)
    plot_bin_factor_changed = pyqtSignal(str, int)
    advanced_requested = pyqtSignal(str)
    load_directory_requested = pyqtSignal(str)  # emits the chosen directory path
    unload_requested = pyqtSignal()             # clear the current results tab

    def __init__(self, parent=None):
        super().__init__(parent)
        self._loaded_data_path: Optional[str] = None
        self._single_datasets: List[Dict[str, Any]] = []
        self._single_curve_entries: List[Dict[str, Any]] = []
        self._single_curve_id_seq: int = 1
        self._plot_row_controls: Dict[str, Dict[str, Any]] = {}
        self._plot_bin_factors: Dict[str, int] = {}
        self._num_groups: Dict[str, Dict[str, Any]] = {}
        self.setMinimumWidth(320)
        self.setMaximumWidth(640)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(8)

        # Single/Multiple plot toggle removed — Single Plot is now its own
        # top-level tab in the main window. Keep an invisible ToggleSwitch
        # so the rest of this widget's code (which still references
        # ``self._mode_switch.isChecked()`` for various branches) keeps
        # behaving as if "Multiple Plots" is always selected.
        # mode_row = QHBoxLayout()
        # mode_row.setContentsMargins(0, 0, 0, 0)
        # mode_row.setSpacing(6)
        # mode_row.addWidget(QLabel("Single Plot"))
        self._mode_switch = ToggleSwitch()
        self._mode_switch.setChecked(True)  # True = multiple mode
        self._mode_switch.toggled.connect(self._on_mode_toggled)
        self._mode_switch.hide()
        # mode_row.addWidget(self._mode_switch)
        # mode_row.addWidget(QLabel("Multiple Plots"))
        # mode_row.addStretch(1)
        # outer.addLayout(mode_row)

        # ---- Plot selection group ----
        plot_group = QGroupBox("Select Plots")
        plot_layout = QVBoxLayout(plot_group)
        plot_layout.setContentsMargins(4, 4, 4, 4)
        plot_layout.setSpacing(4)

        self._mode_stack = QStackedWidget()
        plot_layout.addWidget(self._mode_stack)

        # -- Multiple mode page (existing behavior) --
        multi_page = QWidget()
        multi_layout = QVBoxLayout(multi_page)
        multi_layout.setContentsMargins(0, 0, 0, 0)
        multi_layout.setSpacing(2)

        self._plot_list = ReorderListWidget()
        self._plot_list.setAlternatingRowColors(True)
        self._plot_list.setCursor(Qt.CursorShape.OpenHandCursor)

        # Enable drag-and-drop reordering inside the list
        self._plot_list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._plot_list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._plot_list.setDragDropOverwriteMode(False)

        # Start with an empty plot list. Populated only when results are
        # actually loaded (via update_available_plots).
        self._plot_id_order: List[str] = []

        self._plot_list.reordered.connect(self._on_reordered)

        list_header = QWidget(multi_page)
        list_header_layout = QHBoxLayout(list_header)
        list_header_layout.setContentsMargins(8, 0, 8, 0)
        list_header_layout.setSpacing(6)
        lbl_plots = QLabel("Plots", list_header)
        lbl_plots.setStyleSheet("font-weight: 600;")
        lbl_bins = QLabel("Bins", list_header)
        lbl_bins.setStyleSheet("font-weight: 600;")
        lbl_bins.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl_bins.setFixedWidth(58)
        list_header_layout.addWidget(lbl_plots, 1)
        list_header_layout.addWidget(lbl_bins)
        multi_layout.addWidget(list_header)
        self._plot_list_header = list_header

        multi_layout.addWidget(self._plot_list)

        # Placeholder shown while the list is empty.
        self._plot_empty_label = QLabel(
            "No plots available.\nRun a simulation or load results to see plots here."
        )
        self._plot_empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._plot_empty_label.setWordWrap(True)
        self._plot_empty_label.setStyleSheet(
            "color: #888; font-style: italic; padding: 8px;"
        )
        multi_layout.addWidget(self._plot_empty_label)

        # Reflect empty state in the visibility of list vs placeholder.
        self._plot_list_header.setVisible(False)
        self._plot_list.setVisible(False)

        btn_row = QHBoxLayout()
        self._btn_all = QPushButton("Select All")
        self._btn_none = QPushButton("Deselect All")
        self._btn_all.clicked.connect(self._select_all)
        self._btn_none.clicked.connect(self._deselect_all)
        btn_row.addWidget(self._btn_all)
        btn_row.addWidget(self._btn_none)
        multi_layout.addLayout(btn_row)

        load_row = QHBoxLayout()
        self._btn_load_data = QPushButton("Load Directory…")
        self._btn_load_data.setToolTip("Pick a previous simulation output directory")
        self._btn_load_data.clicked.connect(self._open_load_data_dialog)
        load_row.addWidget(self._btn_load_data, 1)

        self._btn_unload_data = QPushButton("Unload")
        self._btn_unload_data.setToolTip("Clear the results currently shown in this tab")
        self._btn_unload_data.clicked.connect(self.unload_requested)
        load_row.addWidget(self._btn_unload_data)
        multi_layout.addLayout(load_row)

        self._btn_display_settings = AdvancedSettingsButton("Open Advanced Options > Display Settings")
        self._btn_display_settings.clicked.connect(
            lambda: self.advanced_requested.emit("display_settings")
        )
        settings_row = QHBoxLayout()
        settings_row.addStretch(1)
        settings_row.addWidget(self._btn_display_settings)
        multi_layout.addLayout(settings_row)

        self._mode_stack.addWidget(multi_page)

        # -- Single mode page --
        single_page = QWidget()
        single_layout = QVBoxLayout(single_page)
        single_layout.setContentsMargins(0, 0, 0, 0)
        single_layout.setSpacing(6)

        self._btn_single_load_data = QPushButton("Load Data...")
        self._btn_single_load_data.setToolTip("Open a previously simulated data file")
        self._btn_single_load_data.clicked.connect(self._open_load_data_dialog)
        single_layout.addWidget(self._btn_single_load_data)

        ds_row = QHBoxLayout()
        ds_row.addWidget(QLabel("Dataset:"))
        self._single_dataset_combo = QComboBox()
        self._single_dataset_combo.setEnabled(False)
        self._single_dataset_combo.currentIndexChanged.connect(self._on_single_dataset_changed)
        ds_row.addWidget(self._single_dataset_combo, 1)
        single_layout.addLayout(ds_row)

        preview_row = QHBoxLayout()
        preview_row.addWidget(QLabel("Preview Name:"))
        self._single_preview_edit = QLineEdit()
        self._single_preview_edit.setEnabled(False)
        self._single_preview_edit.setPlaceholderText("name - filename")
        self._single_preview_edit.textChanged.connect(self._on_single_preview_changed)
        preview_row.addWidget(self._single_preview_edit, 1)
        single_layout.addLayout(preview_row)

        add_row = QHBoxLayout()
        self._btn_single_add_curve = QPushButton("Add Curve")
        self._btn_single_add_curve.setEnabled(False)
        self._btn_single_add_curve.clicked.connect(self._add_current_single_curve)
        add_row.addWidget(self._btn_single_add_curve)
        single_layout.addLayout(add_row)

        self._single_curve_list = QListWidget()
        self._single_curve_list.setAlternatingRowColors(True)
        self._single_curve_list.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
            | QAbstractItemView.EditTrigger.SelectedClicked
        )
        self._single_curve_list.itemChanged.connect(self._on_single_curve_item_changed)
        single_layout.addWidget(self._single_curve_list)

        curve_btn_row = QHBoxLayout()
        self._btn_single_remove_curve = QPushButton("Remove Selected")
        self._btn_single_remove_curve.clicked.connect(self._remove_selected_single_curve)
        self._btn_single_clear_curves = QPushButton("Clear Curves")
        self._btn_single_clear_curves.clicked.connect(self._clear_single_curves)
        curve_btn_row.addWidget(self._btn_single_remove_curve)
        curve_btn_row.addWidget(self._btn_single_clear_curves)
        single_layout.addLayout(curve_btn_row)

        self._single_status = QLabel("Load calculated data, select a dataset and add one or more curves.")
        self._single_status.setWordWrap(True)
        self._single_status.setStyleSheet("color: #666;")
        single_layout.addWidget(self._single_status)
        single_layout.addStretch(1)

        self._mode_stack.addWidget(single_page)

        outer.addWidget(plot_group)

        # ---- Numerical values group ----
        num_group = QGroupBox("Numerical Results")
        num_layout = QVBoxLayout(num_group)
        num_layout.setContentsMargins(4, 4, 4, 4)

        self._num_table = QTableWidget(0, 2)
        self._num_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        # Interactive resize: user can drag the column separators.
        header = self._num_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(40)
        self._num_table.verticalHeader().setVisible(False)
        self._num_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._num_table.setAlternatingRowColors(True)
        self._num_table.setStyleSheet(
            "QTableWidget { gridline-color: #ccc; }"
        )
        self._num_table.cellClicked.connect(self._on_num_table_cell_clicked)
        self._num_table.cellDoubleClicked.connect(self._on_num_table_double_clicked)

        # Start empty — no default sample data.
        self._num_table.setRowCount(0)

        # Placeholder shown until results are loaded.
        self._num_empty_label = QLabel(
            "No numerical results yet.\nRun a simulation or load results to populate this table."
        )
        self._num_empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._num_empty_label.setWordWrap(True)
        self._num_empty_label.setStyleSheet(
            "color: #888; font-style: italic; padding: 8px;"
        )

        num_layout.addWidget(self._num_table)
        num_layout.addWidget(self._num_empty_label)
        # Initially empty → hide the table, show the placeholder.
        self._num_table.setVisible(False)
        outer.addWidget(num_group)

        self._on_mode_toggled(self._mode_switch.isChecked())

    # --- helpers ------------------------------------------------------

    def selected_plot_ids(self) -> List[str]:
        """Return checked plot IDs in *list order* (determines display order)."""
        ids = []
        for i in range(self._plot_list.count()):
            item = self._plot_list.item(i)
            if item is None:
                continue
            pid = str(item.data(Qt.ItemDataRole.UserRole) or "")
            ctrl = self._plot_row_controls.get(pid, {})
            check = ctrl.get("check")
            if check is not None and check.isChecked():
                ids.append(pid)
        return ids

    def all_plot_ids_in_order(self) -> List[str]:
        """Return all plot IDs (checked or not) in the current list order."""
        out: List[str] = []
        for i in range(self._plot_list.count()):
            item = self._plot_list.item(i)
            if item is None:
                continue
            out.append(str(item.data(Qt.ItemDataRole.UserRole) or ""))
        return out

    def _on_reordered(self):
        """Fired after a drag-and-drop reorder inside the list."""
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _select_all(self):
        for ctrl in self._plot_row_controls.values():
            check = ctrl.get("check")
            if check is not None:
                check.blockSignals(True)
                check.setChecked(True)
                check.blockSignals(False)
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _deselect_all(self):
        for ctrl in self._plot_row_controls.values():
            check = ctrl.get("check")
            if check is not None:
                check.blockSignals(True)
                check.setChecked(False)
                check.blockSignals(False)
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _on_mode_toggled(self, multiple_mode: bool):
        self._mode_stack.setCurrentIndex(0 if multiple_mode else 1)
        self.mode_changed.emit(bool(multiple_mode))
        if multiple_mode:
            self.selection_changed.emit(self.selected_plot_ids())
        else:
            self._emit_single_plot_config()
            self._update_single_status()

    def current_single_dataset(self) -> Optional[Dict[str, Any]]:
        idx = self._single_dataset_combo.currentIndex()
        if 0 <= idx < len(self._single_datasets):
            return self._single_datasets[idx]
        return None

    def _on_single_dataset_changed(self, _index: int):
        dataset = self.current_single_dataset()
        self._btn_single_add_curve.setEnabled(dataset is not None)
        if dataset is not None:
            default_name = _default_preview_name(
                str(dataset.get("name") or "dataset"),
                str(dataset.get("path") or self._loaded_data_path or ""),
            )
            self._single_preview_edit.blockSignals(True)
            self._single_preview_edit.setText(default_name)
            self._single_preview_edit.blockSignals(False)
        self._update_single_status()

    def _on_single_preview_changed(self, _text: str):
        self._update_single_status()

    def _next_single_curve_id(self) -> int:
        cid = self._single_curve_id_seq
        self._single_curve_id_seq += 1
        return cid

    def _refresh_single_curve_list(self) -> None:
        self._single_curve_list.blockSignals(True)
        self._single_curve_list.clear()
        for entry in self._single_curve_entries:
            preview = str(entry.get("preview_name") or "curve")
            item = QListWidgetItem(preview)
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsEditable
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsEnabled
            )
            item.setData(Qt.ItemDataRole.UserRole, int(entry.get("id", -1)))
            item.setToolTip(str(entry.get("path") or ""))
            self._single_curve_list.addItem(item)
        self._single_curve_list.blockSignals(False)

    def _add_current_single_curve(self) -> None:
        dataset = self.current_single_dataset()
        if dataset is None:
            return

        preview = self._single_preview_edit.text().strip()
        if not preview:
            preview = _default_preview_name(
                str(dataset.get("name") or "dataset"),
                str(dataset.get("path") or self._loaded_data_path or ""),
            )

        self._single_curve_entries.append(
            {
                "id": self._next_single_curve_id(),
                "preview_name": preview,
                "dataset_name": str(dataset.get("name") or "dataset"),
                "path": str(dataset.get("path") or ""),
                "dataset": dataset,
            }
        )
        self._refresh_single_curve_list()
        self._emit_single_plot_config()
        self._update_single_status()

    def _on_single_curve_item_changed(self, item: QListWidgetItem) -> None:
        curve_id = int(item.data(Qt.ItemDataRole.UserRole) or -1)
        new_preview = item.text().strip()

        for entry in self._single_curve_entries:
            if int(entry.get("id", -2)) != curve_id:
                continue
            if not new_preview:
                new_preview = str(entry.get("preview_name") or "curve")
                self._single_curve_list.blockSignals(True)
                item.setText(new_preview)
                self._single_curve_list.blockSignals(False)
            entry["preview_name"] = new_preview
            break

        self._emit_single_plot_config()
        self._update_single_status()

    def _remove_selected_single_curve(self) -> None:
        selected = self._single_curve_list.selectedItems()
        if not selected:
            return

        remove_ids = {int(item.data(Qt.ItemDataRole.UserRole) or -1) for item in selected}
        self._single_curve_entries = [
            entry for entry in self._single_curve_entries if int(entry.get("id", -2)) not in remove_ids
        ]
        self._refresh_single_curve_list()
        self._emit_single_plot_config()
        self._update_single_status()

    def _clear_single_curves(self) -> None:
        self._single_curve_entries = []
        self._refresh_single_curve_list()
        self._emit_single_plot_config()
        self._update_single_status()

    def _emit_single_plot_config(self):
        payload = [dict(entry) for entry in self._single_curve_entries]
        self.single_plot_config_changed.emit(payload)

    def _update_single_status(self) -> None:
        file_part = "No file loaded"
        if self._loaded_data_path:
            file_name = os.path.basename(self._loaded_data_path)
            file_part = f"Loaded {len(self._single_datasets)} dataset(s) from {file_name}"
        curve_part = f"Curves in diagram: {len(self._single_curve_entries)}"
        self._single_status.setText(f"{file_part}. {curve_part}.")

    def _open_load_data_dialog(self):
        """Pick a previous simulation output directory and load all plots from it."""
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Simulation Output Directory",
            get_last_used_directory(),
            QFileDialog.Option.ShowDirsOnly,
        )
        if not directory:
            return
        remember_last_used_path(directory)
        # Bubble the chosen path up to MCResultsPage, which knows how to
        # actually read the .his/.mom files and populate this view.
        self.load_directory_requested.emit(directory)

    # --- public API for syncing with external reorder -----------------

    def set_plot_order(self, ordered_ids: List[str]):
        """Reorder the list items to match *ordered_ids*.

        Called when the plot area is reordered via tile drag-and-drop so
        that the sidebar stays in sync.  Check states are preserved.
        """
        self._plot_list.blockSignals(True)

        # Snapshot current items
        data_map: Dict[str, dict] = {}
        for i in range(self._plot_list.count()):
            item = self._plot_list.item(i)
            if item is None:
                continue
            pid = str(item.data(Qt.ItemDataRole.UserRole) or "")
            ctrl = self._plot_row_controls.get(pid, {})
            check = ctrl.get("check")
            data_map[pid] = {
                "name": str(check.text() if check is not None else item.text() or ""),
                "checked": bool(check.isChecked()) if check is not None else True,
                "bin_factor": int(self._plot_bin_factors.get(pid, 1)),
            }

        self._plot_list.clear()
        self._plot_row_controls = {}

        added: set = set()
        for pid in ordered_ids:
            if pid in data_map:
                self._add_list_item(pid, data_map[pid])
                added.add(pid)
        # Append remaining items that were not in ordered_ids
        for pid, d in data_map.items():
            if pid not in added:
                self._add_list_item(pid, d)

        self._plot_list.blockSignals(False)

    def _add_list_item(self, pid: str, data: dict):
        item = QListWidgetItem("")
        item.setFlags(
            item.flags()
            | Qt.ItemFlag.ItemIsDragEnabled
            | Qt.ItemFlag.ItemIsDropEnabled
            | Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
        )
        item.setData(Qt.ItemDataRole.UserRole, pid)
        self._plot_list.addItem(item)

        row = QWidget()
        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(4, 1, 4, 1)
        row_layout.setSpacing(6)

        check = QCheckBox(str(data.get("name", pid)), row)
        check.setChecked(bool(data.get("checked", True)))
        check.toggled.connect(lambda _v, this_pid=pid: self._on_plot_check_changed(this_pid))
        row_layout.addWidget(check, 1)

        bins = QComboBox(row)
        bins.setMinimumWidth(58)
        bins.setMaximumWidth(58)
        for factor, bins_value in ((1, 120), (2, 60), (3, 40), (4, 30), (5, 24), (6, 20)):
            bins.addItem(f"{bins_value}", factor)
        selected_factor = max(1, min(6, int(data.get("bin_factor", self._plot_bin_factors.get(pid, 1)))))
        idx = bins.findData(selected_factor)
        bins.setCurrentIndex(idx if idx >= 0 else 0)
        bins.currentIndexChanged.connect(
            lambda _idx, this_pid=pid, cmb=bins: self._on_plot_bins_changed(this_pid, cmb)
        )
        row_layout.addWidget(bins)

        self._plot_row_controls[pid] = {"check": check, "bins": bins}
        self._plot_list.setItemWidget(item, row)

    def _on_plot_check_changed(self, _pid: str) -> None:
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _on_plot_bins_changed(self, pid: str, combo: QComboBox) -> None:
        factor = int(combo.currentData() or 1)
        self._plot_bin_factors[pid] = factor
        self.plot_bin_factor_changed.emit(pid, factor)

    # --- public API for updating values later -------------------------

    def update_numerical_values(self, values: Any):
        """Replace the numerical values table contents.

        Accepts either:
        - A flat ``list[{"name": str, "value": str}]`` (legacy single-column)
        - A dict ``{"rows": [...], "atom_table": {...}}`` (atom-column layout)
        """
        rows: List[Dict[str, str]] = []
        if isinstance(values, dict):
            rows = list(values.get("rows", []) or [])
        elif isinstance(values, list):
            rows = list(values)

        # Toggle empty-state placeholder vs table visibility.
        has_data = bool(rows)
        if hasattr(self, "_num_empty_label"):
            self._num_empty_label.setVisible(not has_data)
        self._num_table.setVisible(has_data)
        if not has_data:
            self._num_table.setRowCount(0)
            self._num_groups = {}
            return

        # Parameter | Value layout with optional accordion rows.
        self._num_table.setColumnCount(2)
        self._num_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        header = self._num_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(40)
        self._num_table.setColumnWidth(0, 190)
        self._num_table.setRowCount(len(rows))
        self._num_groups = {}

        for row, entry in enumerate(rows):
            row_type = str(entry.get("row_type", ""))
            group_id = str(entry.get("group", ""))
            name_item = QTableWidgetItem(str(entry.get("name", "")))
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            if row_type in {"section", "group"}:
                f = name_item.font()
                f.setBold(True)
                name_item.setFont(f)
                if row_type == "group":
                    prefix = "▸ " if bool(entry.get("collapsed", True)) else "▾ "
                    name_item.setText(prefix + str(entry.get("name", "")))
                    name_item.setToolTip("Click to expand/collapse")

            value_item = QTableWidgetItem(str(entry.get("value", "")))
            value_item.setFlags(value_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            tooltip = str(entry.get("tooltip", "") or "").strip()
            if tooltip:
                value_item.setToolTip(tooltip)
            if row_type == "child":
                name_item.setText("   " + name_item.text())

            self._num_table.setItem(row, 0, name_item)
            self._num_table.setItem(row, 1, value_item)

            if row_type == "group" and group_id:
                self._num_groups[group_id] = {
                    "header": row,
                    "rows": [],
                    "collapsed": bool(entry.get("collapsed", True)),
                }
            elif row_type == "child" and group_id in self._num_groups:
                self._num_groups[group_id]["rows"].append(row)

        for info in self._num_groups.values():
            collapsed = bool(info.get("collapsed", True))
            for r in info.get("rows", []):
                self._num_table.setRowHidden(int(r), collapsed)

    def _on_num_table_cell_clicked(self, row: int, column: int) -> None:
        if column != 0:
            return
        for info in self._num_groups.values():
            if int(info.get("header", -1)) != int(row):
                continue

            collapsed = not bool(info.get("collapsed", True))
            info["collapsed"] = collapsed
            for child_row in info.get("rows", []):
                self._num_table.setRowHidden(int(child_row), collapsed)

            header_item = self._num_table.item(row, 0)
            if header_item is not None:
                text = header_item.text()
                if text.startswith("▸ ") or text.startswith("▾ "):
                    header_item.setText(("▸ " if collapsed else "▾ ") + text[2:])
            return

    def _on_num_table_double_clicked(self, _row: int, _column: int) -> None:
        self._open_num_table_popup()

    def _open_num_table_popup(self) -> None:
        if self._num_table.rowCount() == 0:
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Numerical Results")
        dlg.resize(900, 620)
        layout = QVBoxLayout(dlg)

        table = QTableWidget(self._num_table.rowCount(), self._num_table.columnCount(), dlg)
        labels: List[str] = []
        for c in range(self._num_table.columnCount()):
            hdr = self._num_table.horizontalHeaderItem(c)
            labels.append(hdr.text() if hdr is not None else "")
        table.setHorizontalHeaderLabels(labels)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setAlternatingRowColors(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        table.horizontalHeader().setStretchLastSection(True)
        table.setColumnWidth(0, 240)

        popup_groups: Dict[str, Dict[str, Any]] = {
            gid: {
                "header": int(info.get("header", -1)),
                "rows": [int(r) for r in info.get("rows", [])],
                "collapsed": bool(info.get("collapsed", True)),
            }
            for gid, info in self._num_groups.items()
        }

        for r in range(self._num_table.rowCount()):
            if self._num_table.isRowHidden(r):
                table.setRowHidden(r, True)
            for c in range(self._num_table.columnCount()):
                src = self._num_table.item(r, c)
                if src is None:
                    continue
                clone = QTableWidgetItem(src)
                clone.setFlags(clone.flags() & ~Qt.ItemFlag.ItemIsEditable)
                table.setItem(r, c, clone)

        def _toggle_popup_group(row: int, column: int) -> None:
            if column != 0:
                return
            for info in popup_groups.values():
                if int(info.get("header", -1)) != int(row):
                    continue

                collapsed = not bool(info.get("collapsed", True))
                info["collapsed"] = collapsed
                for child_row in info.get("rows", []):
                    table.setRowHidden(int(child_row), collapsed)

                header_item = table.item(row, 0)
                if header_item is not None:
                    text = header_item.text()
                    if text.startswith("▸ ") or text.startswith("▾ "):
                        header_item.setText(("▸ " if collapsed else "▾ ") + text[2:])
                return

        table.cellClicked.connect(_toggle_popup_group)

        layout.addWidget(table)
        dlg.exec()

    def update_available_plots(self, plots: Dict[str, dict]):
        """Replace the plot selection list."""
        self._plot_list.blockSignals(True)
        self._plot_list.clear()
        self._plot_row_controls = {}
        self._plot_id_order = list(plots.keys())
        for pid, info in plots.items():
            self._add_list_item(
                str(pid),
                {
                    "name": str(info.get("name", pid)),
                    "checked": True,
                    "bin_factor": int(self._plot_bin_factors.get(str(pid), 1)),
                },
            )
        self._plot_list.blockSignals(False)

        # Toggle the empty-state placeholder.
        has_plots = bool(plots)
        if hasattr(self, "_plot_empty_label"):
            self._plot_empty_label.setVisible(not has_plots)
        if hasattr(self, "_plot_list_header"):
            self._plot_list_header.setVisible(has_plots)
        self._plot_list.setVisible(has_plots)

        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())


# =====================================================================
#  MC Results Widget  (public API — imported by mcresults_page.py)
# =====================================================================

class MCResultsWidget(QWidget):
    """Embeddable results widget: sidebar (plot select + values) | plot area.

    The sidebar plot list and the tile grid are kept in sync:
    - Checking / reordering items in the sidebar updates the plot grid.
    - Drag-and-drop reordering of tiles updates the sidebar order.
    """

    advanced_requested = pyqtSignal(str)
    # Emitted on tile double-click — the parent should show the plot in a
    # dedicated Single Plot tab.
    plot_open_in_single = pyqtSignal(str, object)
    # Bubbled up when the user picks a directory from the sidebar's
    # "Load Output Directory" button.
    load_directory_requested = pyqtSignal(str)

    def __init__(self, parent=None, show_toolbar=True):
        super().__init__(parent)

        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        self._splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout.addWidget(self._splitter)

        # Left: sidebar
        self.sidebar = ResultsSidebar()
        self._splitter.addWidget(self.sidebar)

        # Right: stacked plot views (multiple-plot area / single-plot area)
        self.plot_area = PlotArea()
        self.single_plot_area = SinglePlotArea()
        self._plot_stack = QStackedWidget()
        self._plot_stack.addWidget(self.plot_area)
        self._plot_stack.addWidget(self.single_plot_area)
        self._splitter.addWidget(self._plot_stack)

        # Sidebar occupies ~25 % initially
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setSizes([420, 860])

        # Wire selection changes (sidebar → plot area)
        self.sidebar.selection_changed.connect(self._on_sidebar_selection_changed)
        # Wire tile reorder (plot area → sidebar)
        self.plot_area.order_changed.connect(self.sidebar.set_plot_order)
        # Tile double-click → forward to parent (main window) for Single Plot tab
        self.plot_area.plot_double_clicked.connect(self.plot_open_in_single)
        # Wire advanced-options navigation (sidebar → parent)
        self.sidebar.advanced_requested.connect(self.advanced_requested)
        # Wire mode + single-plot selection
        self.sidebar.mode_changed.connect(self._on_mode_changed)
        self.sidebar.single_plot_config_changed.connect(self._on_single_plot_config_changed)
        self.sidebar.plot_bin_factor_changed.connect(self.plot_area.set_plot_bin_factor)
        # Sidebar "Load Output Directory" → bubble to parent page
        self.sidebar.load_directory_requested.connect(self.load_directory_requested)
        # Sidebar "Unload" → clear this tab's results
        self.sidebar.unload_requested.connect(self.unload_results)

        # Default mode: multiple plots
        self._on_mode_changed(True)

    def _on_sidebar_selection_changed(self, plot_ids: List[str]) -> None:
        self.plot_area.set_visible_plots(plot_ids)
        if len(plot_ids) == 1:
            plot_id = str(plot_ids[0])
            info = AVAILABLE_PLOTS.get(plot_id)
            if info is not None:
                self.plot_open_in_single.emit(plot_id, info)

    def _on_mode_changed(self, multiple_mode: bool) -> None:
        self._plot_stack.setCurrentIndex(0 if multiple_mode else 1)

    def _on_single_plot_config_changed(self, curves: List[Dict[str, Any]]) -> None:
        self.single_plot_area.set_curves(curves)

    def set_plot_toolbar_visible(self, visible: bool) -> None:
        self.plot_area.set_toolbar_visible(visible)
        self.single_plot_area.set_toolbar_visible(visible)

    def set_plot_columns(self, count: int) -> None:
        self.plot_area.set_columns(int(count))

    def set_plot_borders_visible(self, visible: bool) -> None:
        self.plot_area.set_borders_visible(visible)
        self.single_plot_area.set_borders_visible(visible)

    def set_plot_font_size(self, size: float) -> None:
        self.plot_area.set_font_size(float(size))
        self.single_plot_area.set_font_size(float(size))

    def set_plot_bin_combine(self, factor: int) -> None:
        self.plot_area.set_bin_factor(int(factor))

    def set_beam_from_top(self, enabled: bool) -> None:
        self.plot_area.set_beam_from_top(bool(enabled))
        if hasattr(self.single_plot_area, "set_beam_from_top"):
            self.single_plot_area.set_beam_from_top(bool(enabled))

    # --- public API for future real-data integration ------------------

    def set_results(
        self,
        plots: Dict[str, dict],
        numerical_values: Any,
    ):
        """Replace both plots and numerical values with new data.

        Parameters
        ----------
        plots : dict
            Same structure as AVAILABLE_PLOTS.
        numerical_values : list[dict] | dict
            Either a flat list of ``{"name", "value"}`` rows, or a dict with
            ``{"rows": [...], "atom_table": {...}}`` for the per-atom layout.
        """
        global AVAILABLE_PLOTS, AVAILABLE_NUMERICAL_VALUES
        AVAILABLE_PLOTS = plots
        AVAILABLE_NUMERICAL_VALUES = numerical_values
        self.sidebar.update_available_plots(plots)
        self.sidebar.update_numerical_values(numerical_values)
        self.plot_area.set_visible_plots(list(plots.keys()))

    def unload_results(self) -> None:
        """Clear all plots and numerical values shown in this results tab."""
        self.set_results({}, {"rows": [], "atom_table": None})


# =====================================================================
#  Standalone main (for quick testing)
# =====================================================================

def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    pal = app.palette()
    pal.setColor(pal.ColorRole.Window, Qt.GlobalColor.white)
    pal.setColor(pal.ColorRole.Base, Qt.GlobalColor.white)
    pal.setColor(pal.ColorRole.Text, Qt.GlobalColor.black)
    pal.setColor(pal.ColorRole.WindowText, Qt.GlobalColor.black)
    app.setPalette(pal)

    w = QMainWindow()
    w.setWindowTitle("OpenSRIM – Results")
    w.resize(1366, 768)
    w.setCentralWidget(MCResultsWidget())
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
