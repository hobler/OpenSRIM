from __future__ import annotations

import json
import subprocess
import threading
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

from PyQt6.QtCore import Qt, QTimer, QRect, QRectF, pyqtSignal
from PyQt6.QtGui import QFontMetrics, QPainter, QTextDocument, QTextOption
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QLabel,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QComboBox,
    QFrame, QDialog, QListWidget, QDialogButtonBox, QMessageBox,
    QSplitter, QStyle, QStyleOptionHeader, QScrollArea, QToolButton,
    QFileDialog,
    QSizePolicy,
)

from state import AppState


class _AdaptiveDecimalSpinBox(QDoubleSpinBox):
    """QDoubleSpinBox that shows only significant decimal places (strips trailing zeros)."""

    def textFromValue(self, value: float) -> str:
        text = f"{value:.10f}".rstrip("0").rstrip(".")
        prefix = self.prefix()
        suffix = self.suffix()
        return f"{prefix}{text}{suffix}"
from ui.widgets.periodic_table_picker import PeriodicTableButton, PeriodicTableDialog

try:
    from simulators.opentrim.read_params import read_params as _read_opentrim_params
except Exception:
    _read_opentrim_params = None  # type: ignore
from ui.dialogs.compound_dictionary_dialog import CompoundDictionaryDialog


def _load_density_table_ang() -> tuple[dict[str, float], dict[str, float]]:
    """Load atomic densities and isotope masses from material_densities.csv.

    Returns
    -------
    densities : dict
        Upper-case symbol → atomic number density in atoms/Å³.
    masses : dict
        Upper-case symbol → dominant isotope mass (mass1 column, amu).
    """
    csv_path = Path(__file__).parents[2] / "data" / "material_densities" / "material_densities.csv"
    densities: dict[str, float] = {}
    masses: dict[str, float] = {}
    try:
        with open(csv_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = [p.strip() for p in line.split(",")]
                if len(parts) < 4:
                    continue
                symbol = parts[0].upper()
                try:
                    masses[symbol] = float(parts[1])
                except ValueError:
                    pass
                try:
                    dens_csv = float(parts[3])
                except ValueError:
                    continue
                densities[symbol] = dens_csv * 0.01  # atoms/Å³
    except OSError:
        pass
    return densities, masses


_DENSITY_TABLE_ANG, _MASS_TABLE = _load_density_table_ang()


def _load_element_energy_defaults() -> dict[int, dict[str, float]]:
    """Load SRIM-compatible default energy parameters per element.

    Returns a dict mapping atomic number Z to
    {"disp": Ed, "latt": Eb, "surf": Es}  (all in eV).
    """
    csv_path = Path(__file__).parents[2] / "data" / "element_energy_defaults" / "element_energy_defaults.csv"
    table: dict[int, dict[str, float]] = {}
    try:
        with open(csv_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = [p.strip() for p in line.split(",")]
                if len(parts) < 5:
                    continue
                try:
                    z = int(parts[0])
                    disp = float(parts[2])
                    latt = float(parts[3])
                    surf = float(parts[4])
                except ValueError:
                    continue
                table[z] = {"disp": disp, "latt": latt, "surf": surf}
    except OSError:
        pass
    return table


_ELEMENT_ENERGY_DEFAULTS: dict[int, dict[str, float]] = _load_element_energy_defaults()

try:
    from ui.logging import log as emit_log
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.logging import log as emit_log  # type: ignore


# HintSystem import (support both possible module locations)
try:
    from OpenSRIM.ui.hints.hints_popup import HintSystem  # type: ignore
except ModuleNotFoundError:
    from OpenSRIM.ui.widgets.hints_popup import HintSystem  # type: ignore


class MCSetupPage(QWidget):
    advanced_requested = pyqtSignal(str)
    save_requested = pyqtSignal()
    load_requested = pyqtSignal()
    simulation_finished = pyqtSignal(str)  # emits results directory path
    results_update = pyqtSignal(str)      # emits results dir for live refresh

    def __init__(self, state: AppState, on_log: Optional[Callable[[str], None]] = None, parent=None):
        super().__init__(parent)
        self.state = state
        self._on_log = on_log or emit_log

        self.layer_elements = []
        self._updating_elements_table = False
        self.latest_log_button = None
        self._logs_dialog = None
        self._logs_list_widget = None
        self.mc_progress = None
        self.run_button = None
        self._progress_timer = None
        self._last_update_ions = 0  # tracks last ion count that triggered a live update
        self.no_of_ions_spin = None
        self.update_after_ions_spin = None
        # Column indices for atoms-per-layer energy columns (accounts for delete column at index 0).
        # Damage column removed; Disp.E./Latt/Surf are now cols 7/8/9.
        self._atoms_disp_col = 7
        self._atoms_latt_col = 8
        self._atoms_surf_col = 9
        self._working_directory: Optional[str] = None
        self._density_user_override: set[int] = set()
        self._updating_layers_table = False
        self._histogram_settings: dict = {}

        # Hints
        self._hint_system: Optional[HintSystem] = None
        self._init_hints()

        # Darker, thicker borders on all container group boxes.
        self.setStyleSheet(
            "QGroupBox { border: 2px solid palette(shadow); border-radius: 4px;"
            " margin-top: 6px; padding-top: 6px; }"
        )

        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        ion_box = self.build_ion_data()

        # Resizable top tables (horizontal)
        top_split = QSplitter(Qt.Orientation.Horizontal)
        top_split.setChildrenCollapsible(False)
        top_split.addWidget(self.build_target_layers())
        top_split.addWidget(self.build_input_elements())
        # Make layer container narrower so atoms-per-layer table gets more space.
        top_split.setStretchFactor(0, 5)
        top_split.setStretchFactor(1, 5)

        # Resizable model/simulator (horizontal)
        model_sim_split = QSplitter(Qt.Orientation.Horizontal)
        model_sim_split.setChildrenCollapsible(False)
        model_sim_split.addWidget(self.build_model_selection())

        # Resizable model/simulator vs output options (vertical)
        bottom_split = QSplitter(Qt.Orientation.Vertical)
        bottom_split.setChildrenCollapsible(False)
        bottom_split.addWidget(model_sim_split)
        bottom_split.addWidget(self.build_trajectories_output())

        # Resizable tables vs bottom section (vertical)
        center = QSplitter(Qt.Orientation.Vertical)
        center.setChildrenCollapsible(False)
        center.addWidget(top_split)
        center.addWidget(bottom_split)

        # Resizable ion selection vs the rest (vertical)
        main_split = QSplitter(Qt.Orientation.Vertical)
        main_split.setChildrenCollapsible(False)
        main_split.addWidget(ion_box)
        main_split.addWidget(center)
        layout.addWidget(main_split, 1)

        layout.addWidget(self._build_mc_setup_footer())

        self._refresh_element_table()

    def _init_hints(self) -> None:
        """Initialize hint system and lock it to this page."""
        try:
            hints_path = Path(__file__).resolve().parents[1] / "widgets" / "hints.json"
            self._hint_system = HintSystem(repo_path=hints_path, parent=self)
            self._hint_system.set_current_page("MC Setup")
        except Exception:
            self._hint_system = None

    def _hint_btn(self, hint_id: str, parent: Optional[QWidget] = None) -> QPushButton:
        """Create a '?' button for this page; disabled if hint system unavailable / hint missing."""
        if self._hint_system is None:
            btn = QPushButton("?", parent)
            btn.setEnabled(False)
            btn.setFixedSize(22, 22)
            btn.setToolTip("Hints not available")
            return btn
        try:
            return self._hint_system.make_hint_button(page_id="MC Setup", hint_id=hint_id, parent=parent)
        except Exception:
            btn = QPushButton("?", parent)
            btn.setEnabled(False)
            btn.setFixedSize(22, 22)
            btn.setToolTip(f"Hint '{hint_id}' not available")
            return btn

    _LOG_PREVIEW_MAX = 120  # max characters shown in footer button

    # -------- external bridge ----------
    def update_latest_log(self, entry: str):
        if self.latest_log_button:
            preview = entry if len(entry) <= self._LOG_PREVIEW_MAX else entry[:self._LOG_PREVIEW_MAX] + "…"
            self.latest_log_button.setText(preview)
            self.latest_log_button.setToolTip(entry)
        if self._logs_list_widget:
            if self._logs_list_widget.count() == 1 and self._logs_list_widget.item(0).text() == "No logs available.":
                self._logs_list_widget.clear()
            self._logs_list_widget.addItem(entry)
            while self._logs_list_widget.count() > 1000:
                self._logs_list_widget.takeItem(0)
            self._logs_list_widget.scrollToBottom()

    def add_log_entry(self, message: str):
        if self._on_log:
            self._on_log(message)

    # -------- configuration (used by MainWindow) ----------
    def collect_simulation_config(self) -> dict:
        try:
            _z = int(self.ion_z.text() or 0)
        except (TypeError, ValueError):
            _z = 0
        ion_data = {
            "symbol": self._ion_symbol_value,
            "name": self.ion_name.text(),
            "number": _z,
            "mass": float(self.ion_mass.text() or 0),
            "energy": self.ion_energy.value(),
            "angle": self.ion_angle.value(),
        }
        # NEW (optional): persist these footer fields
        ions_meta = {
            "no_of_ions": int(self.no_of_ions_spin.value()) if self.no_of_ions_spin else 0,
            "update_after_ions": int(self.update_after_ions_spin.value()) if self.update_after_ions_spin else 0,
        }

        layers = []
        layer_rows = max(self.layers_table.rowCount() - 1, 0) if hasattr(self, "layers_table") else 0
        global_width_unit = self.width_unit_combo.currentText() if hasattr(self, "width_unit_combo") else ""
        global_density_unit = self.density_unit_combo.currentText() if hasattr(self, "density_unit_combo") else "g/cm³"
        for row in range(layer_rows):
            gas_widget = self.layers_table.cellWidget(row, 6)
            gas_checkbox = gas_widget.findChild(QCheckBox) if gas_widget else None
            entries = self.layer_elements[row] if row < len(self.layer_elements) else []
            layer = {
                "name": (self.layers_table.item(row, 1).text() if self.layers_table.item(row, 1) else ""),
                "width": (self.layers_table.item(row, 2).text() if self.layers_table.item(row, 2) else ""),
                "unit": global_width_unit,
                "density": (self.layers_table.item(row, 4).text() if self.layers_table.item(row, 4) else ""),
                "density_unit": global_density_unit,
                "compound_corr": (self.layers_table.item(row, 5).text() if self.layers_table.item(row, 5) else ""),
                "gas": gas_checkbox.isChecked() if gas_checkbox else False,
                "elements": [
                    {
                        "Z": entry["element"]["number"],
                        "symbol": entry["element"]["symbol"],
                        "name": entry["element"]["name"],
                        "mass": entry["element"].get("atomic_mass", 0.0),
                        "ratio": entry["ratio"],
                        "damage": entry["damage"],
                        "disp": entry["disp"],
                        "latt": entry["latt"],
                        "surf": entry["surf"],
                    }
                    for entry in entries
                ],
            }
            layers.append(layer)
        model_meta = {
            "cascade": self.cascade_combo.currentText() if hasattr(self, "cascade_combo") else "",
            "nuclear_stopping": self.nuclear_stopping_combo.currentText() if hasattr(self, "nuclear_stopping_combo") else "",
            "electronic_stopping": self.electronic_stopping_combo.currentText() if hasattr(self, "electronic_stopping_combo") else "",
            "simulator": self.simulator_combo.currentText() if hasattr(self, "simulator_combo") else "",
        }
        output_meta = {
            "traj_start": bool(self.chk_traj_start.isChecked()) if hasattr(self, "chk_traj_start") else False,
            "traj_end": bool(self.chk_traj_end.isChecked()) if hasattr(self, "chk_traj_end") else False,
            "traj_collisions": bool(self.chk_traj_coll.isChecked()) if hasattr(self, "chk_traj_coll") else False,
            "traj_preview": bool(self.chk_traj_preview.isChecked()) if hasattr(self, "chk_traj_preview") else False,
            "range_ion_recoil": bool(self.chk_range_ion_recoil.isChecked()) if hasattr(self, "chk_range_ion_recoil") else False,
            "range_phonons": bool(self.chk_range_phonons.isChecked()) if hasattr(self, "chk_range_phonons") else False,
            "range_ionization": bool(self.chk_range_ionization.isChecked()) if hasattr(self, "chk_range_ionization") else False,
            "lateral_ion_recoil": bool(self.chk_lateral_ion_recoil.isChecked()) if hasattr(self, "chk_lateral_ion_recoil") else False,
            "lateral_phonons": bool(self.chk_lateral_phonons.isChecked()) if hasattr(self, "chk_lateral_phonons") else False,
            "lateral_ionization": bool(self.chk_lateral_ionization.isChecked()) if hasattr(self, "chk_lateral_ionization") else False,
            "dist2d_ion_recoil": bool(self.chk_2d_ion_recoil.isChecked()) if hasattr(self, "chk_2d_ion_recoil") else False,
            "dist2d_phonons": bool(self.chk_2d_phonons.isChecked()) if hasattr(self, "chk_2d_phonons") else False,
            "dist2d_ionization": bool(self.chk_2d_ionization.isChecked()) if hasattr(self, "chk_2d_ionization") else False,
            "dist2d_all": bool(self.chk_2d_all.isChecked()) if hasattr(self, "chk_2d_all") else False,
            "backscattered_energy": bool(self.chk_backscattered_energy.isChecked()) if hasattr(self, "chk_backscattered_energy") else False,
            "backscattered_angle": bool(self.chk_backscattered_angle.isChecked()) if hasattr(self, "chk_backscattered_angle") else False,
            "transmitted_energy": bool(self.chk_transmitted_energy.isChecked()) if hasattr(self, "chk_transmitted_energy") else False,
            "transmitted_angle": bool(self.chk_transmitted_angle.isChecked()) if hasattr(self, "chk_transmitted_angle") else False,
        }

        return {
            "ion": ion_data,
            "ions": ions_meta,
            "layers": layers,
            "selection": model_meta,
            "output": output_meta,
        }

    def apply_simulation_config(self, payload: dict):
        ion = payload.get("ion", {})
        if ion:
            self._ion_symbol_value = ion.get("symbol", "")
            self.ion_name.setText(ion.get("name", ""))
            try:
                self.ion_z.setText(str(int(ion.get("number", self.ion_z.text() or 0))))
            except (TypeError, ValueError):
                pass
            try:
                mass_val = float(ion.get("mass", self.ion_mass.text() or 0))
                self.ion_mass.setText(f"{mass_val:.10f}".rstrip("0").rstrip("."))
            except (TypeError, ValueError):
                pass
            for spin, key in ((self.ion_energy, "energy"), (self.ion_angle, "angle")):
                try:
                    spin.setValue(float(ion.get(key, spin.value())))
                except (TypeError, ValueError):
                    pass

        ions_meta = payload.get("ions") or {}
        if self.no_of_ions_spin is not None:
            try:
                self.no_of_ions_spin.setValue(int(ions_meta.get("no_of_ions", self.no_of_ions_spin.value())))
            except (TypeError, ValueError):
                pass
        if self.update_after_ions_spin is not None:
            try:
                self.update_after_ions_spin.setValue(int(ions_meta.get("update_after_ions", self.update_after_ions_spin.value())))
            except (TypeError, ValueError):
                pass

        layers = payload.get("layers") or []
        self.layers_table.setRowCount(0)
        self.layer_elements = []
        self._density_user_override.clear()

        # Set global unit combos from the first layer's data (if available).
        if layers:
            first = layers[0]
            if hasattr(self, "width_unit_combo"):
                unit = first.get("unit", "Ång")
                idx = self.width_unit_combo.findText(unit)
                if idx >= 0:
                    self.width_unit_combo.setCurrentIndex(idx)
            if hasattr(self, "density_unit_combo"):
                dunit = first.get("density_unit", "g/cm³")
                didx = self.density_unit_combo.findText(dunit)
                if didx >= 0:
                    self.density_unit_combo.setCurrentIndex(didx)

        for idx, layer_data in enumerate(layers):
            self.layers_table.insertRow(idx)
            self.seed_layer_row(idx)
            self.layer_elements.append([])
            self._apply_layer_data(idx, layer_data)

        if not layers:
            self.layers_table.insertRow(0)
            self.seed_layer_row(0)
            self.layer_elements = [[]]

        # Keep a final action row inside the table.
        self.layers_table.insertRow(self.layers_table.rowCount())
        self._ensure_layers_action_row()

        if self.layers_table.rowCount():
            self.layers_table.selectRow(0)
        self._refresh_element_table()

        selection = payload.get("selection") or {}
        if hasattr(self, "cascade_combo"):
            cascade = selection.get("cascade")
            if isinstance(cascade, str) and cascade:
                idx = self.cascade_combo.findText(cascade)
                if idx >= 0:
                    self.cascade_combo.setCurrentIndex(idx)
        if hasattr(self, "nuclear_stopping_combo"):
            nuclear = selection.get("nuclear_stopping")
            if isinstance(nuclear, str) and nuclear:
                idx = self.nuclear_stopping_combo.findText(nuclear)
                if idx >= 0:
                    self.nuclear_stopping_combo.setCurrentIndex(idx)
        if hasattr(self, "electronic_stopping_combo"):
            electronic = selection.get("electronic_stopping")
            if isinstance(electronic, str) and electronic:
                idx = self.electronic_stopping_combo.findText(electronic)
                if idx >= 0:
                    self.electronic_stopping_combo.setCurrentIndex(idx)

        if hasattr(self, "simulator_combo"):
            simulator = selection.get("simulator")
            if isinstance(simulator, str) and simulator:
                idx = self.simulator_combo.findText(simulator)
                if idx >= 0:
                    self.simulator_combo.setCurrentIndex(idx)

        output = payload.get("output") or {}
        for attr, key in (
            ("chk_traj_start", "traj_start"),
            ("chk_traj_end", "traj_end"),
            ("chk_traj_coll", "traj_collisions"),
            ("chk_traj_preview", "traj_preview"),
            ("chk_range_ion_recoil", "range_ion_recoil"),
            ("chk_range_phonons", "range_phonons"),
            ("chk_range_ionization", "range_ionization"),
            ("chk_lateral_ion_recoil", "lateral_ion_recoil"),
            ("chk_lateral_phonons", "lateral_phonons"),
            ("chk_lateral_ionization", "lateral_ionization"),
            ("chk_2d_ion_recoil", "dist2d_ion_recoil"),
            ("chk_2d_phonons", "dist2d_phonons"),
            ("chk_2d_ionization", "dist2d_ionization"),
            ("chk_2d_all", "dist2d_all"),
            ("chk_backscattered_energy", "backscattered_energy"),
            ("chk_backscattered_angle", "backscattered_angle"),
            ("chk_transmitted_energy", "transmitted_energy"),
            ("chk_transmitted_angle", "transmitted_angle"),
        ):
            if hasattr(self, attr) and key in output:
                getattr(self, attr).setChecked(bool(output.get(key)))


    def _apply_layer_data(self, row: int, data: dict):
        for col, key in enumerate(["name", "width", None, "density", "compound_corr"], start=1):
            if key is None:
                continue
            self.layers_table.setItem(row, col, QTableWidgetItem(str(data.get(key, ""))))

        # Apply global width unit from first layer's data (done once in apply_simulation_config).
        # No per-row unit widget needed anymore.

        gas_widget = self.layers_table.cellWidget(row, 6)
        gas_checkbox = gas_widget.findChild(QCheckBox) if gas_widget else None
        if gas_checkbox:
            gas_checkbox.setChecked(bool(data.get("gas", False)))

        self.layer_elements[row] = []
        for entry in data.get("elements", []):
            number = entry.get("Z") or entry.get("number")
            element = self.state.elements_by_number.get(int(number)) if number else None
            if not element:
                continue
            overrides = {k: entry.get(k) for k in ("damage", "disp", "latt", "surf")}
            self._add_element_to_layer(row, element, entry.get("ratio", 0.0), overrides=overrides, refresh=False)

        # Enforce gas/solid-state constraint after loading elements.
        self._enforce_gas_rule_for_layer(row, show_message=False)

    # --------- UI builders ----------
    def _groupbox_header(self, title: str, hint_id: Optional[str] = None, parent: Optional[QWidget] = None) -> QWidget:
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

    def build_ion_data(self) -> QGroupBox:
        box = QGroupBox("")
        v = QVBoxLayout(box)
        v.setSpacing(10)

        v.addWidget(self._groupbox_header("Ion Selection", hint_id="ion", parent=box))

        # Hidden symbol storage (not displayed, used for config).
        self._ion_symbol_value = ""

        grid = QGridLayout()
        grid.setHorizontalSpacing(6)
        grid.setVerticalSpacing(6)

        self.pick_btn = PeriodicTableButton(
            "Select",
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=True,
        )
        self.pick_btn.element_selected.connect(self.on_element_selected)

        # Row 0: all labels in one horizontal line
        grid.addWidget(QLabel("Element"), 0, 0, Qt.AlignmentFlag.AlignLeft)
        grid.addWidget(QLabel("Name"), 0, 1, Qt.AlignmentFlag.AlignLeft)
        grid.addWidget(QLabel("Atomic Number"), 0, 2, Qt.AlignmentFlag.AlignLeft)
        grid.addWidget(QLabel("Mass (amu)"), 0, 3, Qt.AlignmentFlag.AlignLeft)
        grid.addWidget(QLabel("Energy"), 0, 4, Qt.AlignmentFlag.AlignLeft)
        grid.addWidget(QLabel("Tilt"), 0, 5, Qt.AlignmentFlag.AlignLeft)

        # Row 1: all input fields in one horizontal line
        pick_row = QWidget(box)
        pick_row_l = QHBoxLayout(pick_row)
        pick_row_l.setContentsMargins(0, 0, 0, 0)
        pick_row_l.setSpacing(6)
        pick_row_l.addWidget(self.pick_btn)
        pick_row_l.addStretch(1)
        grid.addWidget(pick_row, 1, 0, Qt.AlignmentFlag.AlignLeft)

        self.ion_name = QLineEdit()
        self.ion_name.setReadOnly(True)
        self.ion_name.setPlaceholderText("Element Name")
        self.ion_name.setMaximumWidth(180)
        self.ion_name.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_name, 1, 1, Qt.AlignmentFlag.AlignLeft)

        self.ion_z = QLineEdit()
        self.ion_z.setReadOnly(True)
        self.ion_z.setPlaceholderText("Z")
        self.ion_z.setMaximumWidth(120)
        self.ion_z.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_z, 1, 2, Qt.AlignmentFlag.AlignLeft)

        self.ion_mass = QLineEdit()
        self.ion_mass.setPlaceholderText("Mass")
        self.ion_mass.setMaximumWidth(120)
        self.ion_mass.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_mass, 1, 3, Qt.AlignmentFlag.AlignLeft)

        self.ion_energy = QDoubleSpinBox()
        self.ion_energy.setRange(0.001, 1_000_000)
        self.ion_energy.setValue(10.0)
        self.ion_energy.setMaximumWidth(160)
        self.ion_energy.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_energy, 1, 4, Qt.AlignmentFlag.AlignLeft)

        # Tilt (angle of incidence) — visible next to Energy
        self.ion_angle = QDoubleSpinBox()
        self.ion_angle.setRange(0.0, 90.0)
        self.ion_angle.setDecimals(1)
        self.ion_angle.setValue(0.0)
        self.ion_angle.setMaximumWidth(120)
        self.ion_angle.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_angle, 1, 5, Qt.AlignmentFlag.AlignLeft)

        # Ensure all input fields have the same height.
        _target_h = max(
            self.ion_name.sizeHint().height(),
            self.ion_z.sizeHint().height(),
            self.ion_mass.sizeHint().height(),
            self.ion_energy.sizeHint().height(),
            self.ion_angle.sizeHint().height(),
        )
        for _w in (self.ion_name, self.ion_z, self.ion_mass, self.ion_energy, self.ion_angle):
            _w.setFixedHeight(_target_h)

        ion_settings_btn = QToolButton(box)
        ion_settings_btn.setText("⚙")
        ion_settings_btn.setToolTip("Open ion advanced options")
        ion_settings_btn.clicked.connect(lambda: self.advanced_requested.emit("ion_selection_mc"))
        grid.addWidget(ion_settings_btn, 1, 6, Qt.AlignmentFlag.AlignLeft)

        # Force extra horizontal space to the far right (keeps fields packed on the left)
        for c in range(0, 7):
            grid.setColumnStretch(c, 0)
        grid.setColumnStretch(7, 1)
        v.addLayout(grid)
        return box

    def set_ion_angle(self, angle: float) -> None:
        try:
            self.ion_angle.setValue(float(angle))
        except (TypeError, ValueError):
            return

    def get_ion_angle(self) -> float:
        return float(self.ion_angle.value())

    def on_element_selected(self, element: dict):
        self._ion_symbol_value = element.get("symbol", "")
        self.ion_name.setText(element.get("name", ""))
        try:
            self.ion_z.setText(str(int(element.get("number", self.ion_z.text() or 0))))
        except (TypeError, ValueError):
            pass
        # Prefer isotope mass (mass1) from SRIM density table; fall back to
        # periodic-table atomic_mass.
        sym_upper = self._ion_symbol_value.upper()
        mass_val = _MASS_TABLE.get(sym_upper)
        if mass_val is None:
            try:
                mass_val = float(element.get("atomic_mass", 0))
            except (TypeError, ValueError):
                mass_val = None
        if mass_val is not None and mass_val > 0:
            self.ion_mass.setText(f"{mass_val:.10f}".rstrip("0").rstrip("."))

    def build_target_layers(self) -> QGroupBox:
        box = QGroupBox("")
        v = QVBoxLayout(box)

        # Header row (like KORAL): title + hint on the left, '+' button right-aligned.
        header = QWidget(box)
        header_l = QHBoxLayout(header)
        header_l.setContentsMargins(0, 0, 0, 0)
        header_l.setSpacing(6)
        title_lbl = QLabel("Target layer selection", header)
        title_lbl.setStyleSheet("font-weight: 600;")
        header_l.addWidget(title_lbl)
        header_l.addWidget(self._hint_btn("target_layers", parent=header))
        header_l.addStretch(1)
        v.addWidget(header)

        # Table: first column is a per-row delete button.
        # Extra last row is reserved for the "Add layer" action.
        self.layers_table = QTableWidget(2, 7)
        self.layers_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.layers_table.setHorizontalHeaderLabels(["", "Layer", "Width", "", "Density", "Compound Corr", "Gas"])
        layers_hdr = _WrapHeaderView(self.layers_table)
        self.layers_table.setHorizontalHeader(layers_hdr)
        hdr = self.layers_table.horizontalHeader()
        hdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        try:
            hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
            hdr.setSectionResizeMode(4, QHeaderView.ResizeMode.Fixed)  # Density wider
            hdr.setSectionResizeMode(6, QHeaderView.ResizeMode.Fixed)  # Gas narrow
        except Exception:
            pass
        self.layers_table.setColumnWidth(0, 28)
        self.layers_table.setColumnWidth(4, 110)
        self.layers_table.setColumnWidth(6, 38)
        # Hide the old per-row "Units" column (col 3) — unit is now in the Width header.
        self.layers_table.setColumnHidden(3, True)
        self.layers_table.verticalHeader().setVisible(False)
        self.layers_table.setAlternatingRowColors(True)
        self.layers_table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.layers_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        # Embed unit combo boxes directly inside the header cells.
        # Size to widest entry instead of stretching to full column width.
        self.width_unit_combo = QComboBox()
        self.width_unit_combo.addItems(self.state.unit_options)
        self.width_unit_combo.setCurrentText("Ång")
        self.width_unit_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        layers_hdr.set_header_widget(2, self.width_unit_combo)  # "Width" column

        self.density_unit_combo = QComboBox()
        self.density_unit_combo.addItems(["g/cm³", "kg/m³", "atoms/cm³"])
        self.density_unit_combo.setCurrentText("g/cm³")
        self.density_unit_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        layers_hdr.set_header_widget(4, self.density_unit_combo)  # "Density" column

        self.seed_layer_row(0)
        self._ensure_layers_action_row()
        self.layers_table.itemChanged.connect(self._handle_layer_item_changed)
        v.addWidget(self.layers_table, 1)

        self.layer_elements.append([])
        self.layers_table.selectRow(0)
        self.layers_table.itemSelectionChanged.connect(self._handle_layer_selection_changed)

        # Button is created in the action row.
        return box

    def _ensure_layers_action_row(self) -> None:
        """Ensure the last row of layers_table is an in-table action row."""
        if not hasattr(self, "layers_table"):
            return
        if self.layers_table.columnCount() != 7:
            return

        if self.layers_table.rowCount() == 0:
            self.layers_table.setRowCount(1)

        action_row = self.layers_table.rowCount() - 1

        # Clear any prior widgets in this row.
        for c in range(self.layers_table.columnCount()):
            self.layers_table.setCellWidget(action_row, c, None)
            self.layers_table.setItem(action_row, c, None)

        # Ensure span is applied reliably.
        try:
            self.layers_table.clearSpans()
        except Exception:
            pass

        # Span across all columns.
        try:
            self.layers_table.setSpan(action_row, 0, 1, self.layers_table.columnCount())
        except Exception:
            pass

        cell = QWidget(self.layers_table)
        cell_l = QHBoxLayout(cell)
        cell_l.setContentsMargins(4, 0, 0, 0)
        cell_l.setSpacing(6)
        add_layer = QPushButton("+", cell)
        add_layer.setToolTip("Add layer")
        add_layer.setFixedSize(20, 20)
        add_layer.setStyleSheet("padding: 0px;")
        add_layer.clicked.connect(self.add_layer_row)
        cell_l.addWidget(add_layer)
        cell_l.addStretch(1)
        self.layers_table.setCellWidget(action_row, 0, cell)

    def seed_layer_row(self, r: int):
        # Guard: don't seed the action row.
        if hasattr(self, "layers_table") and r == self.layers_table.rowCount() - 1:
            return
        # Per-row delete button (column 0)
        btn = QPushButton("-")
        btn.setFixedSize(20, 20)
        btn.setStyleSheet("padding: 0px;")
        btn.setToolTip("Delete layer")
        btn.clicked.connect(lambda _=False, b=btn: self._delete_layer_row_for_button(b))
        cell = QWidget(self.layers_table)
        cell_l = QHBoxLayout(cell)
        cell_l.setContentsMargins(0, 0, 0, 0)
        cell_l.setSpacing(0)
        cell_l.addStretch(1)
        cell_l.addWidget(btn)
        cell_l.addStretch(1)
        self.layers_table.setCellWidget(r, 0, cell)

        self.layers_table.setItem(r, 1, QTableWidgetItem(f"Layer {r + 1}"))
        self.layers_table.setItem(r, 2, QTableWidgetItem("10000" if r == 0 else ""))
        # Column 3 (units) is hidden; global unit combo is used instead.
        self._updating_layers_table = True
        self.layers_table.setItem(r, 4, QTableWidgetItem(""))  # auto-filled from elements
        self._updating_layers_table = False
        self.layers_table.setItem(r, 5, QTableWidgetItem("1"))
        gas_chk = QCheckBox()
        gas_chk.toggled.connect(self._handle_layer_gas_toggled)
        gas_widget = QWidget(); lay = QHBoxLayout(gas_widget); lay.setContentsMargins(0, 0, 0, 0)
        lay.addStretch(1); lay.addWidget(gas_chk); lay.addStretch(1)
        self.layers_table.setCellWidget(r, 6, gas_widget)

    def _delete_layer_row_for_button(self, button: QPushButton) -> None:
        if not hasattr(self, "layers_table"):
            return
        for r in range(max(self.layers_table.rowCount() - 1, 0)):
            cell = self.layers_table.cellWidget(r, 0)
            btn = cell.findChild(QPushButton) if cell else None
            if btn is button:
                self._delete_layer_row(r)
                return

    def _delete_layer_row(self, row: int) -> None:
        if not hasattr(self, "layers_table"):
            return
        data_rows = max(self.layers_table.rowCount() - 1, 0)
        if row < 0 or row >= data_rows:
            return
        self.layers_table.removeRow(row)
        if 0 <= row < len(self.layer_elements):
            self.layer_elements.pop(row)
        # Update density override indices after deletion.
        self._density_user_override.discard(row)
        self._density_user_override = {i - 1 if i > row else i for i in self._density_user_override}
        if max(self.layers_table.rowCount() - 1, 0) == 0:
            # Keep one empty layer row + action row.
            if self.layers_table.rowCount() == 0:
                self.layers_table.setRowCount(2)
            else:
                # Ensure we have at least 2 rows (one data + action).
                while self.layers_table.rowCount() < 2:
                    self.layers_table.insertRow(self.layers_table.rowCount())
            self.seed_layer_row(0)
            self.layer_elements = [[]]
            self._density_user_override.clear()

        self._ensure_layers_action_row()
        # Reselect a valid row and refresh the atoms table.
        data_rows = max(self.layers_table.rowCount() - 1, 0)
        if data_rows:
            self.layers_table.selectRow(min(max(row - 1, 0), data_rows - 1))
        self._refresh_element_table()

    def add_layer_row(self):
        action_row = max(self.layers_table.rowCount() - 1, 0)
        self.layers_table.insertRow(action_row)
        self.seed_layer_row(action_row)
        self.layer_elements.insert(action_row, [])
        self._ensure_layers_action_row()
        self.layers_table.selectRow(action_row)

    def delete_selected_layers(self):
        data_rows = max(self.layers_table.rowCount() - 1, 0)
        rows = sorted({idx.row() for idx in self.layers_table.selectedIndexes() if idx.row() < data_rows}, reverse=True)
        for r in rows:
            self.layers_table.removeRow(r)
            if 0 <= r < len(self.layer_elements):
                self.layer_elements.pop(r)
        self._density_user_override.clear()
        if max(self.layers_table.rowCount() - 1, 0) == 0:
            # Keep one empty layer row + action row.
            if self.layers_table.rowCount() == 0:
                self.layers_table.setRowCount(2)
            else:
                while self.layers_table.rowCount() < 2:
                    self.layers_table.insertRow(self.layers_table.rowCount())
            self.seed_layer_row(0)
            self.layer_elements = [[]]
        self._ensure_layers_action_row()
        data_rows = max(self.layers_table.rowCount() - 1, 0)
        if data_rows:
            self.layers_table.selectRow(min(data_rows - 1, 0))
        self._refresh_element_table()

    def _handle_layer_gas_toggled(self, checked: bool) -> None:
        """Prevent gas layers from containing solid-state energy parameters."""
        sender = self.sender()
        if not isinstance(sender, QCheckBox):
            return
        row = self._find_layer_row_for_gas_checkbox(sender)
        if row < 0:
            return
        self._enforce_gas_rule_for_layer(row, show_message=True)
        # If the currently selected layer changed, refresh displayed atoms table.
        if row == self._current_layer_index():
            self._refresh_element_table()

    def _find_layer_row_for_gas_checkbox(self, checkbox: QCheckBox) -> int:
        if not hasattr(self, "layers_table"):
            return -1
        for r in range(max(self.layers_table.rowCount() - 1, 0)):
            gas_widget = self.layers_table.cellWidget(r, 6)
            gas_cb = gas_widget.findChild(QCheckBox) if gas_widget else None
            if gas_cb is checkbox:
                return r
        return -1

    def _is_layer_gas(self, layer_idx: int) -> bool:
        if not hasattr(self, "layers_table"):
            return False
        data_rows = max(self.layers_table.rowCount() - 1, 0)
        if layer_idx < 0 or layer_idx >= data_rows:
            return False
        gas_widget = self.layers_table.cellWidget(layer_idx, 6)
        gas_cb = gas_widget.findChild(QCheckBox) if gas_widget else None
        return bool(gas_cb.isChecked()) if gas_cb else False

    def _enforce_gas_rule_for_layer(self, layer_idx: int, *, show_message: bool) -> None:
        """If layer is gas: store solid params and set energies to 0. If not gas: restore."""
        entries = self._get_layer_entries(layer_idx)
        is_gas = self._is_layer_gas(layer_idx)

        if is_gas:
            had_solid = False
            for entry in entries:
                # Keep a backup so unchecking gas restores prior values.
                if "_solid_energy_backup" not in entry:
                    entry["_solid_energy_backup"] = {
                        "damage": entry.get("damage"),
                        "disp": entry.get("disp"),
                        "latt": entry.get("latt"),
                        "surf": entry.get("surf"),
                    }
                for key in ("damage", "disp", "latt", "surf"):
                    val = entry.get(key)
                    try:
                        had_solid = had_solid or (float(val) != 0.0)
                    except (TypeError, ValueError):
                        had_solid = True
                    entry[key] = "0"
            if show_message and had_solid and entries:
                QMessageBox.information(
                    self,
                    "Gas layer",
                    "Gas layers cannot contain solid-state energy parameters. "
                    "Damage/Disp/Latt/Surf were reset to 0.",
                )
            return

        # Not gas: restore previous values if available, else defaults.
        for entry in entries:
            backup = entry.pop("_solid_energy_backup", None)
            if isinstance(backup, dict) and all(k in backup for k in ("damage", "disp", "latt", "surf")):
                for key in ("damage", "disp", "latt", "surf"):
                    if backup.get(key) is not None:
                        entry[key] = str(backup[key])
            else:
                defaults = self._get_default_energy_params(entry["element"])
                for key in ("damage", "disp", "latt", "surf"):
                    entry[key] = defaults[key]

    def _handle_layer_selection_changed(self):
        self._refresh_element_table()

    def build_input_elements(self) -> QGroupBox:
        box = QGroupBox("")
        v = QVBoxLayout(box)
        # Header row (like KORAL): title + hint left; actions right.
        header = QWidget(box)
        header_l = QHBoxLayout(header)
        header_l.setContentsMargins(0, 0, 0, 0)
        header_l.setSpacing(6)
        title_lbl = QLabel("Atoms per layer", header)
        title_lbl.setStyleSheet("font-weight: 600;")
        header_l.addWidget(title_lbl)
        header_l.addWidget(self._hint_btn("elements", parent=header))
        header_l.addStretch(1)

        dict_btn = QPushButton("Compound Dictionary")
        dict_btn.clicked.connect(self._open_compound_dictionary)
        header_l.addWidget(dict_btn)

        settings_btn = QToolButton()
        settings_btn.setText("⚙")
        settings_btn.setToolTip("Open advanced options")
        settings_btn.clicked.connect(lambda: self.advanced_requested.emit("atoms_per_layer"))
        header_l.addWidget(settings_btn)

        v.addWidget(header)

        # Table: first column is a per-row delete button.
        self.elem_table = QTableWidget(0, 10)
        self.elem_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.elem_table.setHorizontalHeaderLabels([
            "", "Symbol", "Name", "Atomic No.", "Weight\n(amu)",
            "Atom Stoich", "Atom Stoich\n(%)", "Disp", "Latt", "Surf"
        ])
        elem_hdr = _WrapHeaderView(self.elem_table)
        elem_hdr.set_group_header("Damage Energy (eV)", 7, 9)
        self.elem_table.setHorizontalHeader(elem_hdr)
        hdr = self.elem_table.horizontalHeader()
        hdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        try:
            hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)   # delete btn
            hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.Fixed)   # Symbol
            hdr.setSectionResizeMode(7, QHeaderView.ResizeMode.Fixed)   # Disp
            hdr.setSectionResizeMode(8, QHeaderView.ResizeMode.Fixed)   # Latt
            hdr.setSectionResizeMode(9, QHeaderView.ResizeMode.Fixed)   # Surf
        except Exception:
            pass
        self.elem_table.setColumnWidth(0, 28)
        self.elem_table.setColumnWidth(1, 62)
        self.elem_table.setColumnWidth(7, 52)
        self.elem_table.setColumnWidth(8, 52)
        self.elem_table.setColumnWidth(9, 52)
        self.elem_table.verticalHeader().setVisible(False)
        self.elem_table.setAlternatingRowColors(True)

        self.elem_table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.elem_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

        self._set_atoms_energy_columns_visible(show_disp=True, show_latt=True, show_surf=True)

        v.addWidget(self.elem_table, 1)

        # Create the in-table action row (last row).
        self.elem_pick_btn = PeriodicTableButton(
            "+",
            parent=box,
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=False,
        )
        self.elem_pick_btn.setToolTip("Add element")
        self.elem_pick_btn.setFixedSize(20, 20)
        self.elem_pick_btn.setStyleSheet("padding: 0px;")
        self.elem_pick_btn.element_selected.connect(self.on_target_element_selected)

        self.elem_table.itemChanged.connect(self._handle_element_item_changed)
        self.elem_table.cellDoubleClicked.connect(self._handle_element_cell_double_clicked)

        return box

    def _set_atoms_energy_columns_visible(self, *, show_disp: bool, show_latt: bool, show_surf: bool) -> None:
        if not hasattr(self, "elem_table"):
            return
        self.elem_table.setColumnHidden(self._atoms_disp_col, not show_disp)
        self.elem_table.setColumnHidden(self._atoms_latt_col, not show_latt)
        self.elem_table.setColumnHidden(self._atoms_surf_col, not show_surf)

    def on_target_element_selected(self, element: dict):
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        self._add_element_to_layer(layer_idx, element, 1.0)

    def delete_selected_elements(self):
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        data_rows = max(self.elem_table.rowCount() - 1, 0)
        rows = sorted({idx.row() for idx in self.elem_table.selectedIndexes() if idx.row() < data_rows}, reverse=True)
        entries = self._get_layer_entries(layer_idx)
        for r in rows:
            if 0 <= r < len(entries):
                entries.pop(r)
        self._refresh_element_table()

    def build_model_selection(self) -> QGroupBox:
        box = QGroupBox("")
        v = QVBoxLayout(box)
        v.addWidget(self._groupbox_header("Model selection", hint_id="model", parent=box))

        def _combo_with_gear(items: list[str], signal_id: str) -> tuple[QComboBox, QWidget]:
            """Return (combo, container) where container holds combo + gear button flush together."""
            combo = QComboBox()
            combo.addItems(items)
            combo.setMinimumWidth(120)
            btn = QToolButton(box)
            btn.setText("\u2699")
            btn.setToolTip(f"Open {signal_id.replace('_', ' ')} advanced options")
            btn.clicked.connect(lambda: self.advanced_requested.emit(signal_id))
            container = QWidget(box)
            h = QHBoxLayout(container)
            h.setContentsMargins(0, 0, 0, 0)
            h.setSpacing(2)
            h.addWidget(combo)
            h.addWidget(btn)
            h.addStretch(1)
            return combo, container

        grid = QGridLayout()
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(4)

        # Row 0: labels
        grid.addWidget(QLabel("Simulator:"), 0, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignBottom)
        grid.addWidget(QLabel("Cascade option:"), 0, 2, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignBottom)
        grid.addWidget(QLabel("Nuclear stopping:"), 0, 4, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignBottom)
        grid.addWidget(QLabel("Electronic stopping:"), 0, 6, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignBottom)

        # Row 1: combo + gear in a container so button is always flush to combo
        self.simulator_combo = QComboBox()
        self.simulator_combo.addItems(["OpenTRIM", "PyTRIM"])
        self.simulator_combo.setMinimumWidth(120)
        grid.addWidget(self.simulator_combo, 1, 0, Qt.AlignmentFlag.AlignLeft)

        # Vertical separator between simulator and the rest
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.VLine)
        sep.setFrameShadow(QFrame.Shadow.Sunken)
        grid.addWidget(sep, 0, 1, 2, 1)

        self.cascade_combo, _cascade_w = _combo_with_gear(["Ions only", "Full cascade"], "cascade_options")
        grid.addWidget(_cascade_w, 1, 2, Qt.AlignmentFlag.AlignLeft)

        self.nuclear_stopping_combo, _nuclear_w = _combo_with_gear(["ZBL", "NLHlin"], "nuclear_stopping")
        grid.addWidget(_nuclear_w, 1, 4, Qt.AlignmentFlag.AlignLeft)

        self.electronic_stopping_combo, _electronic_w = _combo_with_gear(["SRIM"], "electronic_stopping")
        grid.addWidget(_electronic_w, 1, 6, Qt.AlignmentFlag.AlignLeft)

        grid.setColumnStretch(7, 1)
        v.addLayout(grid)
        v.addStretch(1)
        return box

    def build_trajectories_output(self) -> QGroupBox:
        box = QGroupBox("")
        v = QVBoxLayout(box)
        v.setContentsMargins(8, 8, 8, 8)

        # Header: title + hint.
        header = QWidget(box)
        header_l = QHBoxLayout(header)
        header_l.setContentsMargins(0, 0, 0, 0)
        header_l.setSpacing(6)

        title_lbl = QLabel("Output Options", header)
        title_lbl.setStyleSheet("font-weight: 600;")
        header_l.addWidget(title_lbl)
        header_l.addWidget(self._hint_btn("output", parent=header))
        output_settings_btn = QToolButton(header)
        output_settings_btn.setText("⚙")
        output_settings_btn.setToolTip("Open histogram limits in Advanced Options")
        output_settings_btn.clicked.connect(lambda: self.advanced_requested.emit("histogram_settings"))
        header_l.addWidget(output_settings_btn)
        header_l.addStretch(1)

        v.addWidget(header)

        # Left block (QGridLayout keeps checkboxes vertically aligned within the block)
        self.chk_traj_start = QCheckBox("Start")
        self.chk_traj_end = QCheckBox("End")
        self.chk_traj_coll = QCheckBox("Collisions")
        self.chk_traj_coll.toggled.connect(lambda checked: (
            self.chk_traj_start.setChecked(True),
            self.chk_traj_end.setChecked(True),
        ) if checked else None)
        self.chk_traj_preview = QCheckBox("Preview")

        self.chk_range_ion_recoil = QCheckBox("Ion/Recoil")
        self.chk_range_phonons = QCheckBox("Phonons")
        self.chk_range_ionization = QCheckBox("Ionization")

        self.chk_lateral_ion_recoil = QCheckBox("Ion/Recoil")
        self.chk_lateral_phonons = QCheckBox("Phonons")
        self.chk_lateral_ionization = QCheckBox("Ionization")

        self.chk_2d_ion_recoil = QCheckBox("Ion/Recoil")
        self.chk_2d_phonons = QCheckBox("Phonons")
        self.chk_2d_ionization = QCheckBox("Ionization")

        def _make_all_checkbox(children: list[QCheckBox]) -> QCheckBox:
            """Returns an 'All' checkbox that toggles all children, and updates itself when children change."""
            chk_all = QCheckBox("All")
            def _on_all_clicked():
                # Skip PartiallyChecked when user clicks — jump straight to Checked.
                if chk_all.checkState() == Qt.CheckState.PartiallyChecked:
                    chk_all.setCheckState(Qt.CheckState.Checked)
                checked = chk_all.checkState() == Qt.CheckState.Checked
                for c in children:
                    c.blockSignals(True)
                    c.setChecked(checked)
                    c.blockSignals(False)
            def _on_child_changed():
                all_checked = all(c.isChecked() for c in children)
                any_checked = any(c.isChecked() for c in children)
                chk_all.blockSignals(True)
                chk_all.setCheckState(
                    Qt.CheckState.Checked if all_checked else
                    Qt.CheckState.PartiallyChecked if any_checked else
                    Qt.CheckState.Unchecked
                )
                chk_all.blockSignals(False)
            chk_all.setTristate(True)
            chk_all.clicked.connect(_on_all_clicked)
            for c in children:
                c.toggled.connect(lambda _: _on_child_changed())
            return chk_all

        self.chk_range_all = _make_all_checkbox([self.chk_range_ion_recoil, self.chk_range_phonons, self.chk_range_ionization])
        self.chk_lateral_all = _make_all_checkbox([self.chk_lateral_ion_recoil, self.chk_lateral_phonons, self.chk_lateral_ionization])
        self.chk_2d_all = _make_all_checkbox([self.chk_2d_ion_recoil, self.chk_2d_phonons, self.chk_2d_ionization])

        # Right block widgets
        self.chk_backscattered_energy = QCheckBox("Energy")
        self.chk_backscattered_angle = QCheckBox("Angle")
        self.chk_transmitted_energy = QCheckBox("Energy")
        self.chk_transmitted_angle = QCheckBox("Angle")

        self.chk_backscattered_all = _make_all_checkbox([self.chk_backscattered_energy, self.chk_backscattered_angle])
        self.chk_transmitted_all = _make_all_checkbox([self.chk_transmitted_energy, self.chk_transmitted_angle])

        # Single grid: cols 0-4 = left block, col 5 = gap, cols 6-9 = right block.
        # Row indices are shared so Backscattered aligns with Range (row 1),
        # Transmitted aligns with Lateral Range (row 2).
        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(8)
        grid.setColumnMinimumWidth(5, 24)   # gap between left and right block

        # Row 0: Trajectories (left only)
        grid.addWidget(QLabel("Trajectories:"),       0, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(self.chk_traj_start,           0, 1)
        grid.addWidget(self.chk_traj_end,             0, 2)
        grid.addWidget(self.chk_traj_coll,            0, 3)
        grid.addWidget(self.chk_traj_preview,         0, 4)

        # Row 1: Range (left) + Backscattered (right)
        grid.addWidget(QLabel("Range:"),              1, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(self.chk_range_ion_recoil,     1, 1)
        grid.addWidget(self.chk_range_phonons,        1, 2)
        grid.addWidget(self.chk_range_ionization,     1, 3)
        grid.addWidget(self.chk_range_all,            1, 4)
        grid.addWidget(QLabel("Backscattered:"),      1, 6, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(self.chk_backscattered_energy, 1, 7)
        grid.addWidget(self.chk_backscattered_angle,  1, 8)
        grid.addWidget(self.chk_backscattered_all,    1, 9)

        # Row 2: Lateral Range (left) + Transmitted (right)
        grid.addWidget(QLabel("Lateral Range:"),      2, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(self.chk_lateral_ion_recoil,   2, 1)
        grid.addWidget(self.chk_lateral_phonons,      2, 2)
        grid.addWidget(self.chk_lateral_ionization,   2, 3)
        grid.addWidget(self.chk_lateral_all,          2, 4)
        grid.addWidget(QLabel("Transmitted:"),        2, 6, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(self.chk_transmitted_energy,   2, 7)
        grid.addWidget(self.chk_transmitted_angle,    2, 8)
        grid.addWidget(self.chk_transmitted_all,      2, 9)

        # Row 3: 2D Distribution (left)
        grid.addWidget(QLabel("2D-Distribution:"),   3, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(self.chk_2d_ion_recoil,        3, 1)
        grid.addWidget(self.chk_2d_phonons,           3, 2)
        grid.addWidget(self.chk_2d_ionization,        3, 3)
        grid.addWidget(self.chk_2d_all,               3, 4)

        grid.setRowStretch(4, 1)

        # --- Right side: Simulation controls ---
        _default_nions = 10000
        _default_update = 5000
        if _read_opentrim_params is not None:
            try:
                _p = _read_opentrim_params()
                _sim = _p.get("simulation", {})
                _default_nions = int(_sim.get("nions", 10000))
                _default_update = int(_sim.get("nions_update", 100))
            except Exception:
                pass

        sim_right = QHBoxLayout()

        sim_grid = QGridLayout()
        sim_grid.setHorizontalSpacing(8)
        sim_grid.setVerticalSpacing(8)

        sim_grid.addWidget(QLabel("Working Directory:"), 0, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self._wd_btn = QPushButton("Select", box)
        self._wd_btn.setToolTip("Select the working directory used for outputs")
        self._wd_btn.clicked.connect(self._choose_working_directory)
        sim_grid.addWidget(self._wd_btn, 0, 1, Qt.AlignmentFlag.AlignLeft)

        sim_grid.addWidget(QLabel("No. of Ions:"), 1, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.no_of_ions_spin = QSpinBox()
        self.no_of_ions_spin.setRange(1, 1_000_000_000)
        self.no_of_ions_spin.setSingleStep(1000)
        self.no_of_ions_spin.setValue(_default_nions)
        sim_grid.addWidget(self.no_of_ions_spin, 1, 1, Qt.AlignmentFlag.AlignLeft)

        sim_grid.addWidget(QLabel("Update after Ions:"), 2, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.update_after_ions_spin = QSpinBox()
        self.update_after_ions_spin.setRange(1, 1_000_000_000)
        self.update_after_ions_spin.setSingleStep(100)
        self.update_after_ions_spin.setValue(_default_update)
        sim_grid.addWidget(self.update_after_ions_spin, 2, 1, Qt.AlignmentFlag.AlignLeft)
        sim_grid.setRowStretch(3, 1)

        sim_right.addLayout(sim_grid)

        # Assemble: output options left, simulation controls right.
        content = QHBoxLayout()
        content.setSpacing(80)
        content.addLayout(grid)
        _right_sep = QFrame()
        _right_sep.setFrameShape(QFrame.Shape.VLine)
        _right_sep.setFrameShadow(QFrame.Shadow.Sunken)
        content.addWidget(_right_sep)
        content.addLayout(sim_right)
        content.addStretch(1)

        v.addLayout(content, 1)

        return box

    # -------- footer / logs / progress ----------
    def _build_mc_setup_footer(self) -> QFrame:
        from PyQt6.QtWidgets import QProgressBar  # local import

        footer = QFrame()
        layout = QHBoxLayout(footer)
        layout.setContentsMargins(16, 10, 16, 10)
        layout.setSpacing(12)

        log_container = QWidget(footer)
        log_container_l = QVBoxLayout(log_container)
        log_container_l.setContentsMargins(0, 0, 0, 0)
        log_container_l.setSpacing(2)

        log_btn = QPushButton("No updates yet")
        log_btn.setToolTip("Click to open update notifications")
        log_btn.clicked.connect(self._show_logs_dialog)
        self.latest_log_button = log_btn
        log_container_l.addWidget(log_btn)

        self.mc_progress = QProgressBar()
        self.mc_progress.setRange(0, 100)
        self.mc_progress.setValue(0)
        self.mc_progress.setFormat("Ready")

        load_btn = QPushButton("Load")
        load_btn.setToolTip("Load configuration")
        load_btn.clicked.connect(self.load_requested.emit)

        save_btn = QPushButton("Save")
        save_btn.setToolTip("Save configuration")
        save_btn.clicked.connect(self.save_requested.emit)

        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._handle_run_clicked)

        layout.addWidget(log_container, 3)
        layout.addWidget(self.mc_progress, 1)
        layout.addWidget(load_btn)
        layout.addWidget(save_btn)
        layout.addWidget(self.run_button)
        return footer

    def _choose_working_directory(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Select working directory",
            self._working_directory or str(Path.home()),
        )
        if not path:
            return
        self._working_directory = str(path)
        # Show last path component on button; full path on hover.
        tail = Path(path).name or path
        if hasattr(self, "_wd_btn") and self._wd_btn is not None:
            self._wd_btn.setText(tail)
            self._wd_btn.setToolTip(path)
        self.add_log_entry(f"Working directory set to: {self._working_directory}")

    def _show_logs_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Update Notifications")
        dialog.resize(800, 400)
        layout = QVBoxLayout(dialog)
        list_widget = QListWidget()
        list_widget.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)

        if self.state.log_entries:
            list_widget.addItems(list(self.state.log_entries))
            list_widget.scrollToBottom()
        else:
            list_widget.addItem("No logs available.")

        layout.addWidget(list_widget)

        btn_row = QHBoxLayout()

        copy_btn = QPushButton("Copy Selected")
        copy_btn.setToolTip("Copy selected log entries to clipboard")
        copy_btn.clicked.connect(lambda: self._copy_selected_logs(list_widget))
        btn_row.addWidget(copy_btn)

        clear_btn = QPushButton("Clear Logs")
        clear_btn.clicked.connect(lambda: self._clear_logs(list_widget))
        btn_row.addWidget(clear_btn)

        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)

        self._logs_dialog = dialog
        self._logs_list_widget = list_widget

        def _cleanup(_result: int) -> None:
            self._logs_dialog = None
            self._logs_list_widget = None

        dialog.finished.connect(_cleanup)
        dialog.exec()

    @staticmethod
    def _copy_selected_logs(list_widget: QListWidget) -> None:
        items = list_widget.selectedItems()
        if not items:
            # Nothing selected → copy all
            texts = [list_widget.item(i).text() for i in range(list_widget.count())]
        else:
            texts = [item.text() for item in items]
        from PyQt6.QtWidgets import QApplication
        clipboard = QApplication.clipboard()
        if clipboard is not None:
            clipboard.setText("\n".join(texts))

    def _clear_logs(self, list_widget: QListWidget):
        self.state.clear_logs()
        list_widget.clear()
        list_widget.addItem("No logs available.")
        if self.latest_log_button:
            self.latest_log_button.setText("No updates yet")
            self.latest_log_button.setToolTip("")

    # -------- unit conversion helpers ----------
    _WIDTH_TO_ANGSTROM = {
        "Ång": 1.0,
        "nm": 10.0,
        "µm": 1e4,
        "mm": 1e7,
        "cm": 1e8,
        "m": 1e10,
        "km": 1e13,
    }

    _DENSITY_TO_ATOMS_PER_ANG3 = {
        "atoms/cm³": 1e-24,   # 1 atoms/cm³ = 1e-24 atoms/Å³
    }

    # Avogadro constant
    _N_A = 6.02214076e23

    def _width_in_angstrom(self, value: float, unit: str) -> float:
        return value * self._WIDTH_TO_ANGSTROM.get(unit, 1.0)

    @staticmethod
    def _avg_atomic_mass(elements: list[dict]) -> float:
        """Stoichiometry-weighted average atomic mass for a layer's elements."""
        total_stoich = 0.0
        weighted_mass = 0.0
        for elem in elements:
            try:
                stoich = float(elem.get("ratio", 1.0))
            except (TypeError, ValueError):
                stoich = 1.0
            try:
                mass = float(elem.get("mass", 0.0))
            except (TypeError, ValueError):
                mass = 0.0
            weighted_mass += stoich * mass
            total_stoich += stoich
        return weighted_mass / total_stoich if total_stoich > 0 else 1.0

    @staticmethod
    def _density_from_table(elements: list[dict]) -> float | None:
        """Compute stoichiometry-weighted atomic density (atoms/Å³) from the
        material_densities.csv table.  Returns None if any element is missing."""
        if not elements:
            return None
        total_stoich = 0.0
        weighted = 0.0
        for elem in elements:
            sym = str(elem.get("symbol", "")).upper()
            dens = _DENSITY_TABLE_ANG.get(sym)
            if dens is None:
                return None  # element not in table
            try:
                stoich = float(elem.get("ratio", 1.0))
            except (TypeError, ValueError):
                stoich = 1.0
            weighted += stoich * dens
            total_stoich += stoich
        return weighted / total_stoich if total_stoich > 0 else None

    def _density_in_atoms_per_ang3(self, value: float, unit: str,
                                    elements: list[dict] | None = None) -> float:
        """Convert density to atoms/Å³ which OpenTRIM expects.

        First tries using the material_densities table (most reliable).
        Falls back to converting the user-entered value."""
        # Prefer table density
        if elements:
            table_d = self._density_from_table(elements)
            if table_d is not None:
                return table_d

        # Fallback: convert manually entered value
        if unit == "atoms/cm³":
            return value * 1e-24
        if unit == "kg/m³":
            value = value * 1e-3
            unit = "g/cm³"
        if unit == "g/cm³":
            m_avg = self._avg_atomic_mass(elements) if elements else 1.0
            if m_avg <= 0:
                m_avg = 1.0
            return value * self._N_A / m_avg * 1e-24
        return value

    def set_histogram_settings(self, settings: dict) -> None:
        if isinstance(settings, dict):
            self._histogram_settings = dict(settings)

    # -------- TOML generation ----------
    def _build_input_toml(self, results_dir: str) -> str:
        """Build the contents of an OpenTRIM input.toml from the current UI state."""
        cfg = self.collect_simulation_config()
        ion = cfg["ion"]
        ions = cfg["ions"]
        layers = cfg["layers"]
        sel = cfg["selection"]
        out = cfg["output"]

        width_unit = (self.width_unit_combo.currentText()
                      if hasattr(self, "width_unit_combo") else "Ång")
        density_unit = (self.density_unit_combo.currentText()
                        if hasattr(self, "density_unit_combo") else "g/cm³")

        follow_recoils = sel.get("cascade", "") == "Full cascade"

        lines = [
            "# Auto-generated by OpenSRIM MC Setup",
            "",
            "[simulation]",
            f"follow_recoils = {'true' if follow_recoils else 'false'}",
            f"nions = {ions['no_of_ions']}",
            f"nions_update = {ions['update_after_ions']}",
            f"rng_seed = 12345",
            f'workdir = "{results_dir}"',
            "",
            "[beam]",
            f'symbol = "{ion["symbol"]}"',
            f'name = "{ion["name"]}"',
            f"Z = {ion['number']}",
            f"M = {ion['mass']}",
            f"energy = {ion['energy']}",
            f"tilt = {ion['angle']}",
            "",
        ]

        total_width = 0.0
        for layer in layers:
            try:
                w_raw = float(layer["width"])
            except (TypeError, ValueError):
                w_raw = 0.0
            w_ang = self._width_in_angstrom(w_raw, width_unit)
            total_width += w_ang

            try:
                d_raw = float(layer["density"])
            except (TypeError, ValueError):
                d_raw = 0.0
            d_val = self._density_in_atoms_per_ang3(d_raw, density_unit,
                                                       elements=layer.get("elements", []))

            try:
                cc = float(layer.get("compound_corr") or layer.get("compound_correction", 1.0))
            except (TypeError, ValueError):
                cc = 1.0

            lines += [
                "[[layer]]",
                f'name = "{layer["name"]}"',
                f"width = {w_ang}",
                f"density = {d_val}",
                f"compound_correction = {cc}",
                f"gas = {'true' if layer.get('gas') else 'false'}",
                "",
            ]
            for elem in layer.get("elements", []):
                try:
                    disp_e = float(elem.get("disp", 25.0))
                except (TypeError, ValueError):
                    disp_e = 25.0
                lines += [
                    "    [[layer.element]]",
                    f'    symbol = "{elem["symbol"]}"',
                    f'    name = "{elem["name"]}"',
                    f"    Z = {elem['Z']}",
                    f"    M = {elem.get('mass', 0.0)}",
                    f"    stoichiometry = {elem.get('ratio', 1.0)}",
                    f"    displacement_energy = {disp_e}",
                    "",
                ]

        # Models
        potential = sel.get("nuclear_stopping", "ZBL")
        estop = sel.get("electronic_stopping", "SRIM")
        lines += [
            "[models]",
            f'potential = "{potential}"',
            f'electronic_stopping = "{estop}"',
            "",
            "[models.scattering_integrals]",
            'algorithm = "magic"',
            "n_absc = 4",
            "",
        ]

        # Output section
        energy_kev = ion["energy"]
        half_width = total_width / 2.0 if total_width > 0 else 2000.0
        hs = self._histogram_settings if isinstance(self._histogram_settings, dict) else {}

        def _get_hs(key: str, default: float) -> float:
            if not hs.get("enabled"):
                return float(default)
            try:
                return float(hs.get(key, default))
            except (TypeError, ValueError):
                return float(default)

        depth_min = _get_hs("depth_min", 0.0)
        depth_max = _get_hs("depth_max", total_width)
        lateral_min = _get_hs("lateral_min", -half_width)
        lateral_max = _get_hs("lateral_max", half_width)
        energy_min = _get_hs("energy_min", 0.0)
        energy_max = _get_hs("energy_max", energy_kev)
        angle_min = _get_hs("angle_min", -90.0)
        angle_max = _get_hs("angle_max", 90.0)

        if depth_max <= depth_min:
            depth_min, depth_max = 0.0, total_width
        if lateral_max <= lateral_min:
            lateral_min, lateral_max = -half_width, half_width
        if energy_max <= energy_min:
            energy_min, energy_max = 0.0, energy_kev
        if angle_max <= angle_min:
            angle_min, angle_max = -90.0, 90.0

        try:
            nbins = int(hs.get("nbins", 120)) if hs.get("enabled") else 120
        except (TypeError, ValueError):
            nbins = 120
        if nbins < 10:
            nbins = 10

        lines += [
            "[output]",
            "",
            "[output.trajectories]",
            f"start = {'true' if out.get('traj_start') else 'false'}",
            f"collisions = {'true' if out.get('traj_collisions') else 'false'}",
            f"stopped = {'true' if out.get('traj_end') else 'false'}",
            "backscattered = false",
            "transmitted = false",
            "",
            "[output.depth_distribution.ion_recoils]",
            f"score = {'true' if out.get('range_ion_recoil') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{depth_min}, {depth_max}]",
            "",
            "[output.depth_distribution.nuclear_energy_deposition]",
            f"score = {'true' if out.get('range_phonons') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{depth_min}, {depth_max}]",
            "",
            "[output.depth_distribution.electronic_energy_deposition]",
            f"score = {'true' if out.get('range_ionization') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{depth_min}, {depth_max}]",
            "",
            "[output.lateral_distribution.ion_recoils]",
            f"score = {'true' if out.get('lateral_ion_recoil') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{lateral_min}, {lateral_max}]",
            "",
            "[output.lateral_distribution.nuclear_energy_deposition]",
            f"score = {'true' if out.get('lateral_phonons') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{lateral_min}, {lateral_max}]",
            "",
            "[output.lateral_distribution.electronic_energy_deposition]",
            f"score = {'true' if out.get('lateral_ionization') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{lateral_min}, {lateral_max}]",
            "",
            "[output.backscattered_atoms.energy]",
            f"score = {'true' if out.get('backscattered_energy') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{energy_min}, {energy_max}]",
            "",
            "[output.backscattered_atoms.angle]",
            f"score = {'true' if out.get('backscattered_angle') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{angle_min}, {angle_max}]",
            "",
            "[output.transmitted_atoms.energy]",
            f"score = {'true' if out.get('transmitted_energy') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{energy_min}, {energy_max}]",
            "",
            "[output.transmitted_atoms.angle]",
            f"score = {'true' if out.get('transmitted_angle') else 'false'}",
            f"nbins = {nbins}",
            f"limits = [{angle_min}, {angle_max}]",
            "",
        ]

        return "\n".join(lines)

    # -------- simulation execution ----------
    def _validate_before_run(self) -> str | None:
        """Return an error message if the config is incomplete, else None."""
        errors: list[str] = []

        # Ion
        if not getattr(self, "_ion_symbol_value", ""):
            errors.append("No ion selected.")
        if hasattr(self, "ion_mass") and float(self.ion_mass.text() or 0) <= 0:
            errors.append("Ion mass must be > 0.")
        if hasattr(self, "ion_energy") and self.ion_energy.value() <= 0:
            errors.append("Ion energy must be > 0.")

        # Layers
        layer_count = max(self.layers_table.rowCount() - 1, 0) if hasattr(self, "layers_table") else 0
        if layer_count == 0:
            errors.append("No target layers defined.")
        else:
            for i in range(layer_count):
                entries = self.layer_elements[i] if i < len(self.layer_elements) else []
                if not entries:
                    errors.append(f"Layer {i + 1} has no elements.")
                # Width
                w_item = self.layers_table.item(i, 2) if hasattr(self, "layers_table") else None
                try:
                    w = float(w_item.text()) if w_item else 0.0
                except (TypeError, ValueError):
                    w = 0.0
                if w <= 0:
                    errors.append(f"Layer {i + 1} has no valid width.")
                # Density
                d_item = self.layers_table.item(i, 4) if hasattr(self, "layers_table") else None
                try:
                    d = float(d_item.text()) if d_item else 0.0
                except (TypeError, ValueError):
                    d = 0.0
                if d <= 0:
                    errors.append(f"Layer {i + 1} has no valid density.")

        # Output – at least one scoring option
        cfg = self.collect_simulation_config()
        out = cfg.get("output", {})
        has_output = any(v for k, v in out.items() if isinstance(v, bool))
        if not has_output:
            errors.append("No output options selected (check at least one distribution).")

        # Working directory
        if not self._working_directory:
            errors.append("No working directory set. Click 'Set working directory' first.")

        return "\n".join(errors) if errors else None

    def _handle_run_clicked(self):
        if not self.mc_progress or not self.run_button:
            return

        # Validate
        err = self._validate_before_run()
        if err:
            QMessageBox.warning(self, "Cannot Start Simulation", err)
            return

        base_dir = Path(self._working_directory)  # type: ignore[arg-type]  # validated above

        results_dir = base_dir / datetime.now().strftime("run_%Y%m%d_%H%M%S")
        results_dir.mkdir(parents=True, exist_ok=True)

        toml_path = results_dir / "input.toml"
        toml_content = self._build_input_toml(str(results_dir))
        toml_path.write_text(toml_content, encoding="utf-8")

        self._current_results_dir = str(results_dir)
        self._current_nions = int(self.no_of_ions_spin.value()) if self.no_of_ions_spin else 10000

        self.run_button.setEnabled(False)
        self.mc_progress.setValue(0)
        self.mc_progress.setFormat("0% – Starting…")
        self.add_log_entry(f"Simulation started. Config: {toml_path}")

        self._sim_process = None
        self._sim_error = None

        def _run_in_thread():
            try:
                import sys
                import os as _os
                env = _os.environ.copy()
                env["MPLBACKEND"] = "Agg"  # prevent matplotlib GUI windows
                proc = subprocess.Popen(
                    [sys.executable, "-m", "simulators.run_opentrim", str(toml_path)],
                    cwd=str(Path(__file__).resolve().parents[2]),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    env=env,
                )
                self._sim_process = proc
                proc.wait()
                if proc.returncode != 0:
                    stderr = proc.stderr.read().decode(errors="replace") if proc.stderr else ""
                    self._sim_error = stderr or f"Process exited with code {proc.returncode}"
            except Exception as exc:
                self._sim_error = str(exc)

        self._sim_thread = threading.Thread(target=_run_in_thread, daemon=True)
        self._sim_thread.start()
        self._last_update_ions = 0

        if not self._progress_timer:
            self._progress_timer = QTimer(self)
            self._progress_timer.timeout.connect(self._poll_progress)
        self._progress_timer.start(500)

    def _poll_progress(self):
        if not self.mc_progress:
            return

        # Check if thread is still alive
        thread_alive = hasattr(self, "_sim_thread") and self._sim_thread.is_alive()

        # Read progress file
        progress_file = Path(self._current_results_dir) / "progress"
        done = 0
        if progress_file.exists():
            try:
                text = progress_file.read_text(encoding="utf-8").strip()
                parts = text.split("/")
                if len(parts) == 2:
                    done = int(parts[0].strip())
                    total = int(parts[1].strip())
                    if total > 0:
                        pct = min(int(done * 100 / total), 100)
                        self.mc_progress.setValue(pct)
                        self.mc_progress.setFormat(f"{pct}% – {done}/{total} ions")
            except (ValueError, OSError):
                pass

        # Emit live update when enough new ions have been processed
        update_interval = self.update_after_ions_spin.value() if self.update_after_ions_spin else 5000
        if done > 0 and done - self._last_update_ions >= update_interval:
            self._last_update_ions = done
            self.results_update.emit(self._current_results_dir)

        if not thread_alive:
            if self._progress_timer:
                self._progress_timer.stop()

            # Re-enable button BEFORE any modal dialog
            if self.run_button:
                self.run_button.setEnabled(True)

            if self._sim_error:
                self.mc_progress.setFormat("Error")
                self.mc_progress.setValue(0)
                self.add_log_entry(f"Simulation error: {self._sim_error}")
                # Defer dialog so timer callback returns first
                err_msg = self._sim_error or ""
                QTimer.singleShot(0, lambda: QMessageBox.warning(
                    self, "Simulation Error",
                    f"OpenTRIM failed:\n{err_msg[:500]}"))
            else:
                self.mc_progress.setValue(100)
                self.mc_progress.setFormat("Complete")
                # List output files for transparency
                res_path = Path(self._current_results_dir)
                his_files = sorted(res_path.glob("*.his"))
                mom_files = sorted(res_path.glob("*.mom"))
                self.add_log_entry(
                    f"Simulation completed. {len(his_files)} histogram(s), "
                    f"{len(mom_files)} moment file(s) in: {self._current_results_dir}"
                )
                self.simulation_finished.emit(self._current_results_dir)

    # -------- element/layer logic ----------
    def _open_compound_dictionary(self):
        dialog = CompoundDictionaryDialog(self)
        dialog.compound_selected.connect(self._add_compound_to_layer)
        dialog.exec()

    def _add_compound_to_layer(self, compound: dict):
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        components = [
            part for part in compound.get("composition", [])
            if isinstance(part, dict) and part.get("Z") and part.get("fraction") is not None
        ]
        if not components:
            return
        for part in components:
            element = self.state.elements_by_number.get(int(part["Z"]))
            if not element:
                continue
            ratio = float(part["fraction"])
            self._add_element_to_layer(layer_idx, element, ratio)

    def _add_element_to_layer(self, layer_idx, element, ratio, overrides=None, refresh=True):
        entries = self._get_layer_entries(layer_idx)
        energy_defaults = self._get_default_energy_params(element)
        if overrides:
            for key, value in overrides.items():
                if value is not None:
                    energy_defaults[key] = str(value)
        try:
            ratio_value = float(ratio)
        except (TypeError, ValueError):
            ratio_value = 0.0
        if ratio_value < 0:
            ratio_value = 0.0

        # If the element is already present in this layer, increment its stoichiometry instead of adding a new row.
        try:
            element_number = int(element.get("number"))
        except Exception:
            element_number = None

        if element_number is not None:
            for entry in entries:
                try:
                    existing_number = int(entry.get("element", {}).get("number"))
                except Exception:
                    continue
                if existing_number != element_number:
                    continue
                try:
                    current = float(entry.get("ratio", 0.0) or 0.0)
                except (TypeError, ValueError):
                    current = 0.0
                entry["ratio"] = max(current + ratio_value, 0.0)
                # Enforce gas/solid-state constraint immediately for gas layers.
                self._enforce_gas_rule_for_layer(layer_idx, show_message=False)
                if refresh:
                    self._refresh_element_table()
                return
        entries.append({
            "element": element,
            "ratio": ratio_value,
            "damage": energy_defaults["damage"],
            "disp": energy_defaults["disp"],
            "latt": energy_defaults["latt"],
            "surf": energy_defaults["surf"],
        })

        # Enforce gas/solid-state constraint immediately for gas layers.
        self._enforce_gas_rule_for_layer(layer_idx, show_message=False)
        if refresh:
            self._refresh_element_table()

    def _get_default_energy_params(self, element: dict) -> dict:
        params = {}
        # Try element-specific defaults from the SRIM-compatible CSV table
        try:
            z = int(element.get("number", 0))
        except (TypeError, ValueError):
            z = 0
        table_entry = _ELEMENT_ENERGY_DEFAULTS.get(z, {})
        for key, fallback in self.state.energy_defaults.items():
            if key in table_entry:
                params[key] = str(table_entry[key])
            else:
                candidate = element.get(f"{key}_eV", element.get(key, fallback))
                params[key] = str(candidate)
        return params

    def _refresh_element_table(self):
        if not hasattr(self, "elem_table"):
            return
        layer_idx = self._current_layer_index()
        entries = self._get_layer_entries(layer_idx) if layer_idx >= 0 else []
        self._updating_elements_table = True

        for entry in entries:
            ratio_src = entry.get("ratio", entry.get("stoich", 0.0) or 0.0)
            try:
                entry["ratio"] = float(ratio_src)
            except (TypeError, ValueError):
                entry["ratio"] = 0.0
            defaults = self._get_default_energy_params(entry["element"])
            for key in ("damage", "disp", "latt", "surf"):
                entry.setdefault(key, defaults[key])

        total_ratio = sum(entry["ratio"] for entry in entries)
        # Last row is reserved for the in-table "Add element" action.
        self.elem_table.setRowCount(len(entries) + 1)

        for row, entry in enumerate(entries):
            element = entry["element"]
            mass_raw = element.get("atomic_mass")
            try:
                mass_text = f"{float(mass_raw):.3f}"
            except (TypeError, ValueError):
                mass_text = str(mass_raw)

            def ro_item(text: str):
                it = QTableWidgetItem(text)
                it.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
                return it

            # Per-row delete button (column 0)
            btn = QPushButton("-")
            btn.setFixedSize(20, 20)
            btn.setStyleSheet("padding: 0px;")
            btn.setToolTip("Delete element")
            btn.clicked.connect(lambda _=False, r=row: self._delete_element_row(r))
            cell = QWidget(self.elem_table)
            cell_l = QHBoxLayout(cell)
            cell_l.setContentsMargins(0, 0, 0, 0)
            cell_l.setSpacing(0)
            cell_l.addStretch(1)
            cell_l.addWidget(btn)
            cell_l.addStretch(1)
            self.elem_table.setCellWidget(row, 0, cell)

            self.elem_table.setItem(row, 1, ro_item(element["symbol"]))
            self.elem_table.setItem(row, 2, ro_item(element["name"]))
            self.elem_table.setItem(row, 3, ro_item(str(element["number"])))
            self.elem_table.setItem(row, 4, ro_item(mass_text))

            ratio_item = QTableWidgetItem(f"{entry['ratio']:.4f}")
            ratio_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
            self.elem_table.setItem(row, 5, ratio_item)

            percent = (entry["ratio"] / total_ratio * 100.0) if total_ratio else 0.0
            self.elem_table.setItem(row, 6, ro_item(f"{percent:.2f}"))

            for offset, key in enumerate(("disp", "latt", "surf"), start=7):
                energy_item = QTableWidgetItem(str(entry[key]))
                energy_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
                self.elem_table.setItem(row, offset, energy_item)

        # Action row: span across all columns and place the add-element button.
        action_row = self.elem_table.rowCount() - 1
        for c in range(self.elem_table.columnCount()):
            self.elem_table.setCellWidget(action_row, c, None)
            self.elem_table.setItem(action_row, c, None)

        # Ensure span is applied reliably.
        try:
            self.elem_table.clearSpans()
        except Exception:
            pass
        try:
            self.elem_table.setSpan(action_row, 0, 1, self.elem_table.columnCount())
        except Exception:
            pass

        action_cell = QWidget(self.elem_table)
        action_l = QHBoxLayout(action_cell)
        action_l.setContentsMargins(4, 0, 0, 0)
        action_l.setSpacing(6)
        if hasattr(self, "elem_pick_btn") and self.elem_pick_btn is not None:
            action_l.addWidget(self.elem_pick_btn)
        action_l.addStretch(1)
        self.elem_table.setCellWidget(action_row, 0, action_cell)

        self._updating_elements_table = False
        self._update_layer_density_from_elements()

    def _update_layer_density_from_elements(self) -> None:
        """Auto-fill the layer density cell from the SRIM density table.

        Computes a stoichiometry-weighted density and converts it to the
        currently selected UI unit (g/cm³, kg/m³, or atoms/cm³).
        Skips if the user has manually edited the density for this layer.
        """
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        # Skip if user manually overrode the density for this layer.
        if layer_idx in self._density_user_override:
            return
        entries = self._get_layer_entries(layer_idx)
        if not entries:
            return
        # Build elements list for _density_from_table
        elements = []
        for e in entries:
            elem = e.get("element", {})
            elements.append({
                "symbol": elem.get("symbol", ""),
                "ratio": e.get("ratio", 1.0),
                "mass": elem.get("atomic_mass", 0.0),
            })
        dens_ang = self._density_from_table(elements)
        if dens_ang is None or dens_ang <= 0:
            return
        # Convert atoms/Å³ to the current UI unit
        unit = self.density_unit_combo.currentText() if hasattr(self, "density_unit_combo") else "g/cm³"
        if unit == "atoms/cm³":
            display_val = dens_ang * 1e24
        else:
            m_avg = self._avg_atomic_mass(elements)
            if m_avg <= 0:
                return
            rho_g_cm3 = dens_ang * m_avg / 0.6022
            if unit == "kg/m³":
                display_val = rho_g_cm3 * 1e3
            else:  # g/cm³
                display_val = rho_g_cm3
        # Update the density cell in the layers table
        item = self.layers_table.item(layer_idx, 4)
        if item is None:
            item = QTableWidgetItem()
            self.layers_table.setItem(layer_idx, 4, item)
        self._updating_layers_table = True
        item.setText(f"{display_val:.4f}")
        self._updating_layers_table = False

    def _handle_layer_item_changed(self, item):
        """Track manual edits to the density column (col 4)."""
        if self._updating_layers_table:
            return
        if item.column() == 4:
            self._density_user_override.add(item.row())

    def _handle_element_item_changed(self, item):
        if self._updating_elements_table:
            return
        # Ignore edits in the action row.
        if hasattr(self, "elem_table") and item.row() == self.elem_table.rowCount() - 1:
            return
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        entries = self._get_layer_entries(layer_idx)
        row = item.row()
        if row < 0 or row >= len(entries):
            return
        entry = entries[row]
        col = item.column()
        if col == 5:
            try:
                entry["ratio"] = max(float(item.text()), 0.0)
            except ValueError:
                entry["ratio"] = 0.0
            self._refresh_element_table()
        elif col in (7, 8, 9):
            key = ("disp", "latt", "surf")[col - 7]
            try:
                entry[key] = str(max(float(item.text()), 0.0))
            except ValueError:
                pass

    def _handle_element_cell_double_clicked(self, row, column):
        # Ignore the action row.
        if hasattr(self, "elem_table") and row == self.elem_table.rowCount() - 1:
            return
        if column != 1:
            return
        layer_idx = self._current_layer_index()
        entries = self._get_layer_entries(layer_idx)
        if row < 0 or row >= len(entries):
            return
        dialog = PeriodicTableDialog(self, compact=True, show_hover_info=True, bordered=True)
        dialog.element_selected.connect(lambda element, r=row, idx=layer_idx: self._replace_layer_element(idx, r, element))
        dialog.exec()

    def _replace_layer_element(self, layer_idx, row, element):
        entries = self._get_layer_entries(layer_idx)
        if row < 0 or row >= len(entries):
            return
        entries[row]["element"] = element
        entries[row].update(self._get_default_energy_params(element))
        self._refresh_element_table()

    def _delete_element_row(self, row: int) -> None:
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        entries = self._get_layer_entries(layer_idx)
        if row < 0 or row >= len(entries):
            return
        entries.pop(row)
        self._refresh_element_table()

    def _current_layer_index(self):
        if not hasattr(self, "layers_table"):
            return -1
        row = self.layers_table.currentRow()
        data_rows = max(self.layers_table.rowCount() - 1, 0)
        if data_rows <= 0:
            return -1
        if row < 0:
            row = 0
            self.layers_table.selectRow(0)
        if row >= data_rows:
            row = data_rows - 1
            self.layers_table.selectRow(row)
        return row

    def _get_layer_entries(self, layer_idx):
        while len(self.layer_elements) <= layer_idx:
            self.layer_elements.append([])
        return self.layer_elements[layer_idx] if layer_idx >= 0 else []


class _WrapHeaderView(QHeaderView):
    """QHeaderView that wraps header text instead of eliding it.

    Optionally supports embedding QComboBox widgets below the text of specific
    header sections via ``set_header_widget(logical_index, widget)``.

    Supports group headers that span across multiple columns via
    ``set_group_header(text, first_col, last_col)``.
    """

    def __init__(self, parent=None):
        super().__init__(Qt.Orientation.Horizontal, parent)
        self.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setTextElideMode(Qt.TextElideMode.ElideNone)

        self._sync_scheduled = False
        self._header_widgets: dict[int, QWidget] = {}  # logical_index -> widget
        self._group_headers: list[tuple[str, int, int]] = []  # (text, first_col, last_col)

        # Recompute after geometry changes so Stretch-mode columns are accounted for.
        self.geometriesChanged.connect(self._schedule_sync)
        self.sectionResized.connect(lambda *_: self._schedule_sync())

    # --- public API for embedding widgets in headers ---
    def set_header_widget(self, logical_index: int, widget: QWidget) -> None:
        """Place *widget* inside the header section at *logical_index*."""
        widget.setParent(self.viewport())
        self._header_widgets[logical_index] = widget
        self._schedule_sync()

    def set_group_header(self, text: str, first_col: int, last_col: int) -> None:
        """Register a group header *text* spanning from *first_col* to *last_col*."""
        self._group_headers.append((text, first_col, last_col))
        self._schedule_sync()

    def _bold_font(self):
        f = self.font()
        f.setBold(True)
        return f

    def _schedule_sync(self) -> None:
        if self._sync_scheduled:
            return
        self._sync_scheduled = True
        QTimer.singleShot(0, self._sync_height)

    def _sync_height(self) -> None:
        self._sync_scheduled = False
        h = self.sizeHint().height()
        if h > 0:
            self.setMinimumHeight(h)
            self.setMaximumHeight(h)
        self._reposition_header_widgets()

    def _reposition_header_widgets(self) -> None:
        """Move overlay widgets to their correct position inside the header."""
        for col, widget in self._header_widgets.items():
            if self.isSectionHidden(col):
                widget.setVisible(False)
                continue
            widget.setVisible(True)
            x = self.sectionPosition(col) - self.offset()
            w = self.sectionSize(col)
            total_h = self.height()

            widget_h = widget.sizeHint().height()
            padding = 2
            widget_y = total_h - widget_h - padding

            # Use content-width (sizeHint) instead of full column width, centred.
            content_w = widget.sizeHint().width()
            if content_w < w - 2 * padding:
                widget_x = x + (w - content_w) // 2
                widget.setGeometry(widget_x, widget_y, content_w, widget_h)
            else:
                widget.setGeometry(x + padding, widget_y, w - 2 * padding, widget_h)

    def _header_text(self, logical_index: int) -> str:
        model = self.model()
        if model is None:
            return ""
        value = model.headerData(logical_index, self.orientation(), Qt.ItemDataRole.DisplayRole)
        return "" if value is None else str(value)

    def _text_height_for_section(self, logical_index: int) -> int:
        """Return the text height for a given section."""
        text = self._header_text(logical_index)
        if not text:
            return 0
        bold = self._bold_font()
        doc = QTextDocument()
        doc.setDefaultFont(bold)
        doc.setDocumentMargin(0)
        option = QTextOption()
        option.setAlignment(Qt.AlignmentFlag.AlignCenter)
        option.setWrapMode(QTextOption.WrapMode.WrapAtWordBoundaryOrAnywhere)
        doc.setDefaultTextOption(option)
        doc.setPlainText(text)
        width = max(40, self.sectionSize(logical_index))
        doc.setTextWidth(float(width))
        return int(doc.size().height())

    def _group_header_height(self) -> int:
        """Height needed for the group-header row (0 if no groups)."""
        if not self._group_headers:
            return 0
        fm = QFontMetrics(self._bold_font())
        return fm.height() + 4

    def _group_for_section(self, logical_index: int):
        """Return (text, first, last) if this section belongs to a group, else None."""
        for text, first, last in self._group_headers:
            if first <= logical_index <= last:
                return text, first, last
        return None

    def sizeHint(self):
        base = super().sizeHint()
        model = self.model()
        if model is None:
            return base

        group_h = self._group_header_height()
        height = 0
        for i in range(model.columnCount()):
            if self.isSectionHidden(i):
                continue
            th = self._text_height_for_section(i)
            extra = 0
            if i in self._header_widgets:
                extra = self._header_widgets[i].sizeHint().height() + 4
            # Grouped columns need group_h + their own text height;
            # non-grouped columns need just their text height.
            if self._group_for_section(i) is not None:
                height = max(height, th + extra + group_h)
            else:
                height = max(height, th + extra)

        if height:
            base.setHeight(max(base.height(), height + 8))
        return base

    def paintSection(self, painter, rect, logicalIndex):
        if not rect.isValid():
            return
        painter.save()

        # If this section belongs to a group, reserve the top strip for the group header
        # by shrinking the section drawing rect downward.
        in_group = self._group_for_section(logicalIndex) is not None
        group_h = self._group_header_height() if in_group else 0
        section_rect = rect.adjusted(0, group_h, 0, 0) if group_h > 0 else rect

        opt = QStyleOptionHeader()
        self.initStyleOption(opt)
        opt.rect = section_rect
        opt.section = logicalIndex
        opt.text = ""

        # draw background/section (only in the lower portion when grouped)
        self.style().drawControl(QStyle.ControlElement.CE_HeaderSection, opt, painter, self)

        text_rect = self.style().subElementRect(QStyle.SubElement.SE_HeaderLabel, opt, self)
        text = self._header_text(logicalIndex)

        # If this section has an embedded widget, only use the top portion for text
        if logicalIndex in self._header_widgets:
            widget_h = self._header_widgets[logicalIndex].sizeHint().height() + 4
            text_rect.setHeight(max(text_rect.height() - widget_h, 10))

        doc = QTextDocument()
        doc.setDefaultFont(self._bold_font())
        doc.setDocumentMargin(0)
        option = QTextOption()
        option.setAlignment(Qt.AlignmentFlag.AlignCenter)
        option.setWrapMode(QTextOption.WrapMode.WrapAtWordBoundaryOrAnywhere)
        doc.setDefaultTextOption(option)
        doc.setPlainText(text)
        doc.setTextWidth(float(text_rect.width()))

        painter.translate(text_rect.topLeft())
        clip = QRectF(0.0, 0.0, float(text_rect.width()), float(text_rect.height()))
        doc.drawContents(painter, clip)
        painter.restore()

    def paintEvent(self, event):
        super().paintEvent(event)
        if not self._group_headers:
            return
        painter = QPainter(self.viewport())
        painter.setFont(self._bold_font())
        group_h = self._group_header_height()
        for text, first_col, last_col in self._group_headers:
            # Find bounding rect across visible columns in the group.
            x_start = None
            total_w = 0
            for c in range(first_col, last_col + 1):
                if self.isSectionHidden(c):
                    continue
                pos = self.sectionPosition(c) - self.offset()
                if x_start is None:
                    x_start = pos
                total_w = pos + self.sectionSize(c) - x_start
            if x_start is None or total_w <= 0:
                continue

            # Draw a real merged header cell (background + border) for the group.
            cell_rect = QRect(int(x_start), 0, int(total_w), int(group_h))
            opt = QStyleOptionHeader()
            self.initStyleOption(opt)
            opt.rect = cell_rect
            opt.section = first_col
            opt.text = ""
            self.style().drawControl(QStyle.ControlElement.CE_HeaderSection, opt, painter, self)

            # Draw the group text centered inside the merged cell.
            painter.drawText(QRectF(cell_rect), Qt.AlignmentFlag.AlignCenter, text)
        painter.end()
