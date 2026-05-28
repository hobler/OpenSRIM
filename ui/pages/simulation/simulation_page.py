import csv
import json
import os
import sys
import numpy as np
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
    QColorDialog, QFormLayout,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.image as mpimg

try:
    from ui.widgets.toggle_switch import ToggleSwitch
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.widgets.toggle_switch import ToggleSwitch  # type: ignore


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
    ax.title.set_fontsize(base + 1.0)
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

    def _call_plot_func(self, ax) -> None:
        """Invoke the plot_func, passing bin_factor when the function accepts it."""
        fn = self.plot_info["plot_func"]
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
                                bin_factor=self._bin_factor)
                tile.zoom_requested.connect(self._open_zoomed)
                tile.tile_reorder_requested.connect(self._handle_tile_reorder)
                tile.setMinimumSize(250, 200)
                # Apply current display settings
                tile.toolbar.setVisible(self._toolbar_visible)
                tile.set_border_visible(self._borders_visible)
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

    def __init__(self, parent=None, font_size: float = 10.0):
        super().__init__(parent)
        # Match the darker / thicker group-box borders used on the MC Setup
        # and KORAL pages.
        self.setStyleSheet(
            "QGroupBox { border: 2px solid palette(shadow); border-radius: 4px;"
            " margin-top: 6px; padding-top: 6px; }"
        )
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

        self._hint_label = QLabel(
            "Double-click a plot tile on the MC Results tab to load its curves here. "
            "Then pick a curve in the right-hand list and enable Gaussian Convolution "
            "to smooth it (σ in plot units, scan area for lateral curves)."
        )
        self._hint_label.setStyleSheet(
            "color: #555; padding: 4px 6px; background: #f5f5f5; "
            "border-bottom: 1px solid #ddd;"
        )
        self._hint_label.setWordWrap(True)
        plot_layout.addWidget(self._hint_label)

        self.figure = Figure(figsize=(8, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, 1)

        splitter.addWidget(plot_box)

        # ===== right: side panel (scrollable) =====
        side_scroll = QScrollArea()
        side_scroll.setWidgetResizable(True)
        side_scroll.setFrameShape(QFrame.Shape.NoFrame)
        side = QWidget()
        side.setMinimumWidth(300)
        side_layout = QVBoxLayout(side)
        side_layout.setContentsMargins(6, 6, 6, 6)
        side_layout.setSpacing(8)
        side_scroll.setWidget(side)

        # ---- Curves group ----
        curves_group = QGroupBox("Curves")
        cg_layout = QVBoxLayout(curves_group)
        cg_layout.setContentsMargins(6, 6, 6, 6)
        cg_layout.setSpacing(4)

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

        side_layout.addWidget(curves_group)

        # ---- Load previous simulation output ----
        load_group = QGroupBox("Load Output Directory")
        load_layout = QVBoxLayout(load_group)
        load_layout.setContentsMargins(6, 6, 6, 6)
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
        self._style_group = QGroupBox("Selected Curve")
        style_form = QFormLayout(self._style_group)
        style_form.setContentsMargins(6, 6, 6, 6)
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

        side_layout.addWidget(self._style_group)
        self._style_group.setEnabled(False)

        # ---- Convolution / scan area (per-curve) ----
        self._conv_group = QGroupBox("Gaussian Convolution")
        conv_layout = QVBoxLayout(self._conv_group)
        conv_layout.setContentsMargins(6, 6, 6, 6)
        conv_layout.setSpacing(4)

        # Depth-profile convolution controls
        self._conv_depth_widget = QWidget()
        conv_depth_layout = QVBoxLayout(self._conv_depth_widget)
        conv_depth_layout.setContentsMargins(0, 0, 0, 0)
        conv_depth_layout.setSpacing(4)

        self._conv_check = QCheckBox("Apply convolution")
        self._conv_check.toggled.connect(self._on_conv_toggled)
        conv_depth_layout.addWidget(self._conv_check)

        sigma_row = QHBoxLayout()
        sigma_row.addWidget(QLabel("Beam σ (Å):"))
        self._sigma_spin = QDoubleSpinBox()
        self._sigma_spin.setRange(0.0, 1_000_000.0)
        self._sigma_spin.setDecimals(2)
        self._sigma_spin.setValue(0.0)
        self._sigma_spin.setSingleStep(5.0)
        self._sigma_spin.valueChanged.connect(self._on_sigma_changed)
        sigma_row.addWidget(self._sigma_spin, 1)
        conv_depth_layout.addLayout(sigma_row)

        ch = QLabel("Folds the selected depth-profile curve with a Gaussian of width σ.")
        ch.setWordWrap(True)
        ch.setStyleSheet("color: #888; font-size: 10px;")
        conv_depth_layout.addWidget(ch)
        conv_layout.addWidget(self._conv_depth_widget)

        # Scan-area controls (2-D / lateral curves only)
        self._scan_area_widget = QWidget()
        scan_layout = QVBoxLayout(self._scan_area_widget)
        scan_layout.setContentsMargins(0, 0, 0, 0)
        scan_layout.setSpacing(4)

        scan_row = QHBoxLayout()
        scan_row.addWidget(QLabel("Scan area:"))
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

        sh = QLabel("Lateral scan area applied to the selected 2-D curve.")
        sh.setWordWrap(True)
        sh.setStyleSheet("color: #888; font-size: 10px;")
        scan_layout.addWidget(sh)
        conv_layout.addWidget(self._scan_area_widget)

        # Place the convolution panel high in the side bar (right after
        # "Curves") so it doesn't get lost below the larger style form.
        side_layout.insertWidget(1, self._conv_group)
        # Always interactive — even without a selected curve. The handlers
        # auto-pick the first 1-D curve as soon as the user touches anything,
        # so the panel never sits "broken-looking" in a greyed-out state.

        # ---- Axes & legend ----
        axes_group = QGroupBox("Axes & Legend")
        axes_form = QFormLayout(axes_group)
        axes_form.setContentsMargins(6, 6, 6, 6)
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

        # ---- Reference image overlay ----
        ref_group = QGroupBox("Reference Image")
        ref_layout = QVBoxLayout(ref_group)
        ref_layout.setContentsMargins(6, 6, 6, 6)
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
        self._stats_group = QGroupBox("Stats (Selected)")
        stats_layout = QVBoxLayout(self._stats_group)
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
        self._stats_table.setAlternatingRowColors(True)
        stats_layout.addWidget(self._stats_table)
        side_layout.addWidget(self._stats_group, 1)
        self._stats_group.setVisible(False)

        side_layout.addStretch(0)

        splitter.addWidget(side_scroll)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([800, 320])

        self._update_secondary_axes_visibility()
        self._render_plot()

    # =====================================================================
    # Public API
    # =====================================================================

    def set_font_size(self, size: float) -> None:
        self._font_size = float(size)
        self._render_plot()

    # Backwards-compat shim — callers that used set_plot now add instead.
    def set_plot(self, plot_id: str, plot_info: Dict[str, Any]) -> None:
        self.add_plot(plot_id, plot_info)

    def add_plot(self, plot_id: str, plot_info: Dict[str, Any]) -> int:
        """Append every visible series of *plot_info* as a new curve."""
        if not isinstance(plot_info, dict):
            return 0

        # 2D plots can't overlay with line curves — replace the curve list
        # with the single heatmap entry. The actual drawing is delegated to
        # the plot_func supplied by the results builder.
        if plot_info.get("is_2d") and callable(plot_info.get("plot_func")):
            source_name = str(plot_info.get("name", plot_id))
            x_label = str(plot_info.get("x_label", "") or "")
            y_label = str(plot_info.get("y_label", "") or "")
            self._curves = []
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
                "x_label":     x_label,
                "y_label":     y_label,
                # Placeholders so the curve-list/style panel code paths don't
                # blow up when they read these fields.
                "color":       "#444444",
                "drawstyle":   "default",
                "linestyle":   "-",
                "linewidth":   1.0,
                "marker":      "None",
                "markersize":  4.0,
                "alpha":       1.0,
                "zorder":      1,
                "conv_enabled": False,
                "conv_sigma":   0.0,
                "scan_area":    [0.0, 0.0],
                "bin_factor":   1,
                "is_depth_profile": False,
                "axis_x":       "bottom",
                "axis_y":       "left",
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
            "",
            QFileDialog.Option.ShowDirsOnly,
        )
        if not directory:
            return
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
                "Dataset does not contain valid numeric x/y data.",
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
            self._stats_table.setRowCount(0)
            return

        # 2D curves: line-style controls don't apply, but the Gauss filter
        # (and the erf scan-area along the lateral axis) does — populate the
        # convolution fields from the curve.
        if curve.get("is_2d"):
            self._style_group.setEnabled(False)
            self._conv_depth_widget.setVisible(True)
            self._scan_area_widget.setVisible(True)
            self._conv_check.blockSignals(True)
            self._conv_check.setChecked(bool(curve.get("conv_enabled", False)))
            self._conv_check.blockSignals(False)
            self._sigma_spin.blockSignals(True)
            self._sigma_spin.setValue(float(curve.get("conv_sigma", 0.0)))
            self._sigma_spin.blockSignals(False)
            sa = curve.get("scan_area") or [0.0, 0.0]
            self._scan_area_min.blockSignals(True)
            self._scan_area_max.blockSignals(True)
            self._scan_area_min.setValue(float(sa[0]) if len(sa) > 0 else 0.0)
            self._scan_area_max.setValue(float(sa[1]) if len(sa) > 1 else 0.0)
            self._scan_area_min.blockSignals(False)
            self._scan_area_max.blockSignals(False)
            self._bin_combine_spin.blockSignals(True)
            self._bin_combine_spin.setValue(1)
            self._bin_combine_spin.blockSignals(False)
            self._label_edit.blockSignals(True)
            self._label_edit.setText(str(curve["label"]))
            self._label_edit.blockSignals(False)
            self._populate_stats_for(curve)
            return

        self._label_edit.blockSignals(True)
        self._label_edit.setText(str(curve["label"]))
        self._label_edit.blockSignals(False)

        self._update_color_button(curve["color"])

        self._set_combo_data(self._drawstyle_cmb, curve["drawstyle"])
        self._set_combo_data(self._linestyle_cmb, curve["linestyle"])
        self._set_combo_data(self._marker_cmb, curve["marker"])

        self._linewidth_spin.blockSignals(True)
        self._linewidth_spin.setValue(float(curve["linewidth"]))
        self._linewidth_spin.blockSignals(False)

        self._marker_size_spin.blockSignals(True)
        self._marker_size_spin.setValue(float(curve["markersize"]))
        self._marker_size_spin.blockSignals(False)

        self._alpha_spin.blockSignals(True)
        self._alpha_spin.setValue(float(curve["alpha"]))
        self._alpha_spin.blockSignals(False)

        self._zorder_spin.blockSignals(True)
        self._zorder_spin.setValue(int(curve["zorder"]))
        self._zorder_spin.blockSignals(False)

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
        self._sigma_spin.setValue(float(curve.get("conv_sigma", 0.0)))
        self._sigma_spin.blockSignals(False)
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
        for idx, c in enumerate(self._curves):
            if not c.get("is_2d"):
                self._curve_list.setCurrentRow(idx)
                return self._selected_curve()
        return None

    def _on_conv_toggled(self, checked: bool) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            # No curves loaded yet — undo the toggle so the UI doesn't lie
            # about the state, and nudge the user via the status hint.
            self._conv_check.blockSignals(True)
            self._conv_check.setChecked(False)
            self._conv_check.blockSignals(False)
            self._hint_label.setText(
                "Load a curve first: double-click a plot tile on the MC Results tab."
            )
            return
        curve["conv_enabled"] = bool(checked)
        # When the user first turns convolution on with σ still at zero, pick
        # a sensible non-zero default so the effect is immediately visible.
        # Heuristic: ~3 % of the curve's x-range.
        if checked and float(curve.get("conv_sigma", 0.0)) <= 0.0:
            if curve.get("is_2d"):
                x = np.asarray(curve.get("x_values_2d", []), dtype=float)
            else:
                x = np.asarray(curve.get("x", []), dtype=float)
            if x.size >= 2:
                span = float(x[-1] - x[0])
                sigma_default = max(span * 0.03, 0.0)
                curve["conv_sigma"] = sigma_default
                self._sigma_spin.blockSignals(True)
                self._sigma_spin.setValue(sigma_default)
                self._sigma_spin.blockSignals(False)
        self._render_plot()

    def _on_sigma_changed(self, v: float) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            return
        curve["conv_sigma"] = float(v)
        if curve.get("conv_enabled"):
            self._render_plot()

    def _on_scan_area_changed(self, _v: float) -> None:
        curve = self._ensure_curve_selected()
        if curve is None:
            return
        curve["scan_area"] = [
            float(self._scan_area_min.value()),
            float(self._scan_area_max.value()),
        ]
        if curve.get("conv_enabled"):
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
            self, "Select Reference Image", "",
            "Images (*.png *.jpg *.jpeg);;All files (*)",
        )
        if not path:
            return
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
        # 2D plots: defer to the stats_func provided by the results builder.
        if curve.get("is_2d"):
            rows: List[tuple] = [("Source", str(curve.get("source_name", "")))]
            stats_func = curve.get("stats_func")
            if callable(stats_func):
                try:
                    rows.extend((str(n), str(v)) for n, v in stats_func())
                except Exception:
                    pass
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
        if curve.get("conv_enabled"):
            y = self._gauss_convolve(x, y, float(curve.get("conv_sigma", 0.0)),
                                     scan_area=curve.get("scan_area"))

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

        # 2D plot: a single heatmap drives the whole figure. Delegate to the
        # plot_func provided by the results builder and skip the curve loop.
        twod_curve = next(
            (c for c in self._curves if c.get("is_2d") and c.get("visible", True)),
            None,
        )
        if twod_curve is not None:
            conv_kwargs = {}
            if twod_curve.get("conv_enabled"):
                conv_kwargs["conv"] = {
                    "sigma": float(twod_curve.get("conv_sigma", 0.0) or 0.0),
                    "scan_area": twod_curve.get("scan_area"),
                }
            try:
                try:
                    twod_curve["plot_func"](self.ax, **conv_kwargs)
                except TypeError:
                    # Older plot_func without conv kwarg.
                    twod_curve["plot_func"](self.ax)
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
            # Override the title/labels with the page's editable axes state.
            if self._axes_state.get("title"):
                self.ax.set_title(str(self._axes_state["title"]))
            if self._axes_state.get("xlabel"):
                self.ax.set_xlabel(str(self._axes_state["xlabel"]))
            if self._axes_state.get("ylabel"):
                self.ax.set_ylabel(str(self._axes_state["ylabel"]))
            _apply_axes_font_size(self.ax, self._font_size)
            try:
                self.figure.tight_layout()
            except Exception:
                pass
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
            if curve.get("conv_enabled"):
                y = self._gauss_convolve(x, y, float(curve.get("conv_sigma", 0.0)),
                                     scan_area=curve.get("scan_area"))

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
        if self._ref_image_data is not None and any_plotted:
            try:
                xlim = ax_BL.get_xlim()
                ylim = ax_BL.get_ylim()
                ax_BL.imshow(
                    self._ref_image_data,
                    extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
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
        self.canvas.draw_idle()

        # Stats follow the currently selected curve (if any)
        sel = self._selected_curve()
        if sel is not None:
            self._populate_stats_for(sel)


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
    advanced_requested = pyqtSignal(str)
    load_directory_requested = pyqtSignal(str)  # emits the chosen directory path

    def __init__(self, parent=None):
        super().__init__(parent)
        self._loaded_data_path: Optional[str] = None
        self._single_datasets: List[Dict[str, Any]] = []
        self._single_curve_entries: List[Dict[str, Any]] = []
        self._single_curve_id_seq: int = 1
        self.setMinimumWidth(260)
        self.setMaximumWidth(380)

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

        # Start with an empty plot list. Populated only when results are
        # actually loaded (via update_available_plots).
        self._plot_id_order: List[str] = []

        self._plot_list.itemChanged.connect(self._on_item_changed)
        self._plot_list.reordered.connect(self._on_reordered)

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

        hint = QLabel("\u21c5 Drag items to reorder plots")
        hint.setStyleSheet("color: #888; font-size: 10px;")
        hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        multi_layout.addWidget(hint)
        self._plot_list_hint = hint

        # Reflect empty state in the visibility of list/hint vs placeholder.
        self._plot_list.setVisible(False)
        self._plot_list_hint.setVisible(False)

        btn_row = QHBoxLayout()
        self._btn_all = QPushButton("Select All")
        self._btn_none = QPushButton("Deselect All")
        self._btn_all.clicked.connect(self._select_all)
        self._btn_none.clicked.connect(self._deselect_all)
        btn_row.addWidget(self._btn_all)
        btn_row.addWidget(self._btn_none)
        multi_layout.addLayout(btn_row)

        self._btn_load_data = QPushButton("Load Output Directory…")
        self._btn_load_data.setToolTip("Pick a previous simulation output directory")
        self._btn_load_data.clicked.connect(self._open_load_data_dialog)
        multi_layout.addWidget(self._btn_load_data)

        self._btn_display_settings = QToolButton()
        self._btn_display_settings.setToolTip("Open Advanced Options > Display Settings")
        self._btn_display_settings.setAutoRaise(True)
        self._btn_display_settings.setIconSize(QSize(18, 18))
        icon = QIcon.fromTheme("preferences-system")
        if icon.isNull():
            icon = self.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogDetailedView)
        self._btn_display_settings.setIcon(icon)
        self._btn_display_settings.setCursor(Qt.CursorShape.PointingHandCursor)
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
            if item.checkState() == Qt.CheckState.Checked:
                ids.append(item.data(Qt.ItemDataRole.UserRole))
        return ids

    def all_plot_ids_in_order(self) -> List[str]:
        """Return all plot IDs (checked or not) in the current list order."""
        return [
            self._plot_list.item(i).data(Qt.ItemDataRole.UserRole)
            for i in range(self._plot_list.count())
        ]

    def _on_item_changed(self, _item):
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _on_reordered(self):
        """Fired after a drag-and-drop reorder inside the list."""
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _select_all(self):
        self._plot_list.blockSignals(True)
        for i in range(self._plot_list.count()):
            self._plot_list.item(i).setCheckState(Qt.CheckState.Checked)
        self._plot_list.blockSignals(False)
        if self._mode_switch.isChecked():
            self.selection_changed.emit(self.selected_plot_ids())

    def _deselect_all(self):
        self._plot_list.blockSignals(True)
        for i in range(self._plot_list.count()):
            self._plot_list.item(i).setCheckState(Qt.CheckState.Unchecked)
        self._plot_list.blockSignals(False)
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
            "",
            QFileDialog.Option.ShowDirsOnly,
        )
        if not directory:
            return
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
            pid = item.data(Qt.ItemDataRole.UserRole)
            data_map[pid] = {
                "name": item.text(),
                "checked": item.checkState(),
            }

        self._plot_list.clear()

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
        item = QListWidgetItem(data["name"])
        item.setFlags(
            item.flags()
            | Qt.ItemFlag.ItemIsUserCheckable
            | Qt.ItemFlag.ItemIsDragEnabled
            | Qt.ItemFlag.ItemIsDropEnabled
        )
        item.setCheckState(data["checked"])
        item.setData(Qt.ItemDataRole.UserRole, pid)
        self._plot_list.addItem(item)

    # --- public API for updating values later -------------------------

    def update_numerical_values(self, values: Any):
        """Replace the numerical values table contents.

        Accepts either:
        - A flat ``list[{"name": str, "value": str}]`` (legacy single-column)
        - A dict ``{"rows": [...], "atom_table": {...}}`` (atom-column layout)
        """
        rows: List[Dict[str, str]] = []
        atom_table: Optional[Dict[str, Any]] = None
        if isinstance(values, dict):
            rows = list(values.get("rows", []) or [])
            atom_table = values.get("atom_table")
        elif isinstance(values, list):
            rows = list(values)

        atoms: List[str] = []
        atom_rows: List[tuple] = []
        if isinstance(atom_table, dict):
            atoms = list(atom_table.get("atoms", []) or [])
            atom_rows = list(atom_table.get("rows", []) or [])

        # Toggle empty-state placeholder vs table visibility.
        has_data = bool(atoms) or bool(rows)
        if hasattr(self, "_num_empty_label"):
            self._num_empty_label.setVisible(not has_data)
        self._num_table.setVisible(has_data)
        if not has_data:
            self._num_table.setRowCount(0)
            return

        if atoms:
            n_cols = 1 + len(atoms)
            self._num_table.setColumnCount(n_cols)
            headers = ["Parameter"] + atoms
            self._num_table.setHorizontalHeaderLabels(headers)
            header = self._num_table.horizontalHeader()
            # Interactive: user can drag column dividers freely.
            for c in range(n_cols):
                header.setSectionResizeMode(c, QHeaderView.ResizeMode.Interactive)
            header.setStretchLastSection(True)
            header.setMinimumSectionSize(40)
            # Sensible initial widths — narrow Parameter column, compact atoms.
            self._num_table.setColumnWidth(0, 130)
            for c in range(1, n_cols):
                self._num_table.setColumnWidth(c, 70)

            total_rows = len(atom_rows) + len(rows)
            self._num_table.setRowCount(total_rows)

            for r, (name, cells) in enumerate(atom_rows):
                name_item = QTableWidgetItem(str(name))
                name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._num_table.setItem(r, 0, name_item)
                # Visually distinguish section headers (rows whose value cells
                # are empty strings).
                is_header = all(not str(cells.get(a, "")).strip() for a in atoms)
                if is_header:
                    font = name_item.font()
                    font.setBold(True)
                    name_item.setFont(font)
                for c, atom in enumerate(atoms, start=1):
                    cell = QTableWidgetItem(str(cells.get(atom, "")))
                    cell.setFlags(cell.flags() & ~Qt.ItemFlag.ItemIsEditable)
                    cell.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                    self._num_table.setItem(r, c, cell)

            for r, entry in enumerate(rows, start=len(atom_rows)):
                name_item = QTableWidgetItem(str(entry.get("name", "")))
                name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self._num_table.setItem(r, 0, name_item)
                value_item = QTableWidgetItem(str(entry.get("value", "")))
                value_item.setFlags(value_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                # Span the value across the remaining columns.
                for c in range(1, n_cols):
                    self._num_table.setItem(r, c, QTableWidgetItem(""))
                self._num_table.setItem(r, 1, value_item)
                if n_cols > 2:
                    self._num_table.setSpan(r, 1, 1, n_cols - 1)
            return

        # Legacy / no atom table: keep simple Parameter | Value layout.
        self._num_table.setColumnCount(2)
        self._num_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        header = self._num_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(40)
        self._num_table.setColumnWidth(0, 130)
        self._num_table.setRowCount(len(rows))
        for row, entry in enumerate(rows):
            name_item = QTableWidgetItem(str(entry.get("name", "")))
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            value_item = QTableWidgetItem(str(entry.get("value", "")))
            value_item.setFlags(value_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._num_table.setItem(row, 0, name_item)
            self._num_table.setItem(row, 1, value_item)

    def update_available_plots(self, plots: Dict[str, dict]):
        """Replace the plot selection list."""
        self._plot_list.blockSignals(True)
        self._plot_list.clear()
        self._plot_id_order = list(plots.keys())
        for pid, info in plots.items():
            item = QListWidgetItem(info["name"])
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsDragEnabled
                | Qt.ItemFlag.ItemIsDropEnabled
            )
            item.setCheckState(Qt.CheckState.Checked)
            item.setData(Qt.ItemDataRole.UserRole, pid)
            self._plot_list.addItem(item)
        self._plot_list.blockSignals(False)

        # Toggle the empty-state placeholder.
        has_plots = bool(plots)
        if hasattr(self, "_plot_empty_label"):
            self._plot_empty_label.setVisible(not has_plots)
        self._plot_list.setVisible(has_plots)
        if hasattr(self, "_plot_list_hint"):
            self._plot_list_hint.setVisible(has_plots)

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
        self._splitter.setSizes([280, 900])

        # Wire selection changes (sidebar → plot area)
        self.sidebar.selection_changed.connect(self.plot_area.set_visible_plots)
        # Wire tile reorder (plot area → sidebar)
        self.plot_area.order_changed.connect(self.sidebar.set_plot_order)
        # Tile double-click → forward to parent (main window) for Single Plot tab
        self.plot_area.plot_double_clicked.connect(self.plot_open_in_single)
        # Wire advanced-options navigation (sidebar → parent)
        self.sidebar.advanced_requested.connect(self.advanced_requested)
        # Wire mode + single-plot selection
        self.sidebar.mode_changed.connect(self._on_mode_changed)
        self.sidebar.single_plot_config_changed.connect(self._on_single_plot_config_changed)
        # Sidebar "Load Output Directory" → bubble to parent page
        self.sidebar.load_directory_requested.connect(self.load_directory_requested)

        # Default mode: multiple plots
        self._on_mode_changed(True)

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
