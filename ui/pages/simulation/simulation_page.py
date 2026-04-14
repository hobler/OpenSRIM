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
    QLineEdit, QMessageBox,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

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

        self._stats_table = QTableWidget(4, 2)
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
                for entry in AVAILABLE_NUMERICAL_VALUES[:4]
            ]

        rows = rows[:4]
        while len(rows) < 4:
            rows.append((f"Value {len(rows) + 1}", "-"))
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

    def __init__(self, plot_id: str, plot_info: dict, parent=None, font_size: float = 10.0):
        super().__init__(parent)
        self.plot_id = plot_id
        self.plot_info = plot_info
        self._is_3d = plot_info.get("projection") == "3d"
        self._has_border = True
        self._font_size = float(font_size)
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
        plot_info["plot_func"](self.ax)
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
        self.ax.clear()
        self.plot_info["plot_func"](self.ax)
        _apply_axes_font_size(self.ax, self._font_size)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def set_font_size(self, font_size: float):
        self._font_size = float(font_size)
        self.redraw()


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

    def __init__(self, parent=None):
        super().__init__(parent)

        self._plot_order: List[str] = []
        self._tiles: Dict[str, PlotTile] = {}

        # Display settings
        self._col_override: int = 0   # 0 = auto
        self._toolbar_visible: bool = True
        self._borders_visible: bool = True
        self._font_size: float = 10.0

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
                tile = PlotTile(pid, info, font_size=self._font_size)
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
        dlg = ZoomedPlotDialog(plot_id, info, self, font_size=self._font_size)
        dlg.exec()


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

        mode_row = QHBoxLayout()
        mode_row.setContentsMargins(0, 0, 0, 0)
        mode_row.setSpacing(6)
        mode_row.addWidget(QLabel("Single Plot"))
        self._mode_switch = ToggleSwitch()
        self._mode_switch.setChecked(True)  # True = multiple mode
        self._mode_switch.toggled.connect(self._on_mode_toggled)
        mode_row.addWidget(self._mode_switch)
        mode_row.addWidget(QLabel("Multiple Plots"))
        mode_row.addStretch(1)
        outer.addLayout(mode_row)

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

        self._plot_id_order: List[str] = []
        for pid, info in AVAILABLE_PLOTS.items():
            item = QListWidgetItem(info["name"])
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsDragEnabled
                | Qt.ItemFlag.ItemIsDropEnabled
            )
            item.setCheckState(Qt.CheckState.Unchecked)
            item.setData(Qt.ItemDataRole.UserRole, pid)
            self._plot_list.addItem(item)
            self._plot_id_order.append(pid)

        self._plot_list.itemChanged.connect(self._on_item_changed)
        self._plot_list.reordered.connect(self._on_reordered)

        multi_layout.addWidget(self._plot_list)

        hint = QLabel("\u21c5 Drag items to reorder plots")
        hint.setStyleSheet("color: #888; font-size: 10px;")
        hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        multi_layout.addWidget(hint)

        btn_row = QHBoxLayout()
        self._btn_all = QPushButton("Select All")
        self._btn_none = QPushButton("Deselect All")
        self._btn_all.clicked.connect(self._select_all)
        self._btn_none.clicked.connect(self._deselect_all)
        btn_row.addWidget(self._btn_all)
        btn_row.addWidget(self._btn_none)
        multi_layout.addLayout(btn_row)

        self._btn_load_data = QPushButton("Load Data...")
        self._btn_load_data.setToolTip("Open a previously simulated data file")
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

        self._num_table = QTableWidget(len(AVAILABLE_NUMERICAL_VALUES), 2)
        self._num_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self._num_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self._num_table.verticalHeader().setVisible(False)
        self._num_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._num_table.setAlternatingRowColors(True)
        self._num_table.setStyleSheet(
            "QTableWidget { gridline-color: #ccc; }"
        )

        for row, entry in enumerate(AVAILABLE_NUMERICAL_VALUES):
            name_item = QTableWidgetItem(entry["name"])
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            val_item = QTableWidgetItem(entry["value"])
            val_item.setFlags(val_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._num_table.setItem(row, 0, name_item)
            self._num_table.setItem(row, 1, val_item)

        num_layout.addWidget(self._num_table)
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
        """Open file picker for already simulated data and parse datasets."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Simulated Data",
            "",
            "JSON Files (*.json);;CSV Files (*.csv);;Text Files (*.txt);;All Files (*)",
        )
        if not path:
            return

        try:
            datasets = load_datasets_from_path(path)
        except Exception as exc:
            QMessageBox.warning(self, "Load Data", f"Unable to load data file:\n{exc}")
            return

        if not datasets:
            QMessageBox.information(
                self,
                "Load Data",
                "No plottable numeric datasets were found in this file.",
            )
            return

        self._loaded_data_path = path
        self._single_datasets = datasets
        self._single_dataset_combo.blockSignals(True)
        self._single_dataset_combo.clear()
        for ds in datasets:
            self._single_dataset_combo.addItem(str(ds.get("name") or "dataset"))
        self._single_dataset_combo.blockSignals(False)

        self._single_dataset_combo.setEnabled(True)
        self._single_preview_edit.setEnabled(True)
        self._btn_single_add_curve.setEnabled(True)
        self._single_dataset_combo.setCurrentIndex(0)
        self._update_single_status()

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

    def update_numerical_values(self, values: List[Dict[str, str]]):
        """Replace the numerical values table contents."""
        self._num_table.setRowCount(len(values))
        for row, entry in enumerate(values):
            self._num_table.setItem(row, 0, QTableWidgetItem(entry["name"]))
            self._num_table.setItem(row, 1, QTableWidgetItem(entry["value"]))

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
        # Wire advanced-options navigation (sidebar → parent)
        self.sidebar.advanced_requested.connect(self.advanced_requested)
        # Wire mode + single-plot selection
        self.sidebar.mode_changed.connect(self._on_mode_changed)
        self.sidebar.single_plot_config_changed.connect(self._on_single_plot_config_changed)

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

    # --- public API for future real-data integration ------------------

    def set_results(
        self,
        plots: Dict[str, dict],
        numerical_values: List[Dict[str, str]],
    ):
        """Replace both plots and numerical values with new data.

        Parameters
        ----------
        plots : dict
            Same structure as AVAILABLE_PLOTS.
        numerical_values : list[dict]
            Same structure as AVAILABLE_NUMERICAL_VALUES.
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
