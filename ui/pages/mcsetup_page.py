from __future__ import annotations

from datetime import datetime
from typing import Callable, Optional

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QLabel,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QPushButton,
    QTableWidget, QTableWidgetItem, QHeaderView, QComboBox,
    QFrame, QDialog, QListWidget, QDialogButtonBox, QMessageBox,
)

from app.state import AppState
from app.ui.widgets.periodic_table_picker import PeriodicTableButton, PeriodicTableDialog
from app.ui.dialogs.compound_dictionary_dialog import CompoundDictionaryDialog


class MCSetupPage(QWidget):
    def __init__(self, state: AppState, on_log: Optional[Callable[[str], None]] = None, parent=None):
        super().__init__(parent)
        self.state = state
        self._on_log = on_log

        self.layer_elements = []
        self._updating_elements_table = False
        self.latest_log_button = None
        self.mc_progress = None
        self.run_button = None
        self._progress_timer = None
        self.no_of_ions_spin = None
        self.update_after_ions_spin = None

        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        ion_box = self.build_ion_data()
        ion_box.setTitle("Ion Selection")
        layout.addWidget(ion_box)

        middle = QHBoxLayout()
        middle.setSpacing(10)
        middle.addWidget(self.build_target_layers(), stretch=1)
        middle.addWidget(self.build_input_elements(), stretch=1)
        layout.addLayout(middle)

        layout.addWidget(self.build_model_selection())
        layout.addWidget(self.build_trajectories_output())

        layout.addStretch(1)
        layout.addWidget(self._build_mc_setup_footer())

        self._refresh_element_table()

    # -------- external bridge ----------
    def update_latest_log(self, entry: str):
        if self.latest_log_button:
            self.latest_log_button.setText(entry)
            self.latest_log_button.setToolTip(entry)

    def add_log_entry(self, message: str):
        if self._on_log:
            self._on_log(message)

    # -------- configuration (used by MainWindow) ----------
    def collect_simulation_config(self) -> dict:
        ion_data = {
            "symbol": self.ion_symbol.text(),
            "name": self.ion_name.text(),
            "number": self.ion_z.value(),
            "mass": self.ion_mass.value(),
            "energy": self.ion_energy.value(),
            "angle": self.ion_angle.value(),
        }
        # NEW (optional): persist these footer fields
        ions_meta = {
            "no_of_ions": int(self.no_of_ions_spin.value()) if self.no_of_ions_spin else 0,
            "update_after_ions": int(self.update_after_ions_spin.value()) if self.update_after_ions_spin else 0,
        }

        layers = []
        for row in range(self.layers_table.rowCount()):
            unit_widget = self.layers_table.cellWidget(row, 2)
            unit = unit_widget.currentText() if isinstance(unit_widget, QComboBox) else ""
            gas_widget = self.layers_table.cellWidget(row, 5)
            gas_checkbox = gas_widget.findChild(QCheckBox) if gas_widget else None
            entries = self.layer_elements[row] if row < len(self.layer_elements) else []
            layer = {
                "name": (self.layers_table.item(row, 0).text() if self.layers_table.item(row, 0) else ""),
                "width": (self.layers_table.item(row, 1).text() if self.layers_table.item(row, 1) else ""),
                "unit": unit,
                "density": (self.layers_table.item(row, 3).text() if self.layers_table.item(row, 3) else ""),
                "compound_corr": (self.layers_table.item(row, 4).text() if self.layers_table.item(row, 4) else ""),
                "gas": gas_checkbox.isChecked() if gas_checkbox else False,
                "elements": [
                    {
                        "Z": entry["element"]["number"],
                        "symbol": entry["element"]["symbol"],
                        "name": entry["element"]["name"],
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
        return {"ion": ion_data, "ions": ions_meta, "layers": layers}

    def apply_simulation_config(self, payload: dict):
        ion = payload.get("ion", {})
        if ion:
            self.ion_symbol.setText(ion.get("symbol", ""))
            self.ion_name.setText(ion.get("name", ""))
            try:
                self.ion_z.setValue(int(ion.get("number", self.ion_z.value())))
            except (TypeError, ValueError):
                pass
            for spin, key in ((self.ion_mass, "mass"), (self.ion_energy, "energy"), (self.ion_angle, "angle")):
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
        for idx, layer_data in enumerate(layers):
            self.layers_table.insertRow(idx)
            self.seed_layer_row(idx)
            self.layer_elements.append([])
            self._apply_layer_data(idx, layer_data)

        if not layers:
            self.layers_table.insertRow(0)
            self.seed_layer_row(0)
            self.layer_elements = [[]]

        if self.layers_table.rowCount():
            self.layers_table.selectRow(0)
        self._refresh_element_table()

    def _apply_layer_data(self, row: int, data: dict):
        for col, key in enumerate(["name", "width", None, "density", "compound_corr"]):
            if key is None:
                continue
            self.layers_table.setItem(row, col, QTableWidgetItem(str(data.get(key, ""))))

        unit_widget = self.layers_table.cellWidget(row, 2)
        if isinstance(unit_widget, QComboBox):
            idx = unit_widget.findText(data.get("unit", "Ång"))
            if idx >= 0:
                unit_widget.setCurrentIndex(idx)

        gas_widget = self.layers_table.cellWidget(row, 5)
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

    # --------- UI builders ----------
    def build_ion_data(self) -> QGroupBox:
        box = QGroupBox("Ion Data")
        grid = QGridLayout(box)
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(10)

        self.pick_btn = PeriodicTableButton(
            "Click to Select Element",
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=True,
        )
        self.pick_btn.element_selected.connect(self.on_element_selected)

        self.ion_symbol = QLineEdit(); self.ion_symbol.setReadOnly(True); self.ion_symbol.setPlaceholderText("Symbol")
        self.ion_name = QLineEdit(); self.ion_name.setReadOnly(True); self.ion_name.setPlaceholderText("Element Name")
        self.ion_z = QSpinBox(); self.ion_z.setRange(1, 120); self.ion_z.setReadOnly(True)
        self.ion_mass = QDoubleSpinBox(); self.ion_mass.setRange(0.01, 1000.0); self.ion_mass.setDecimals(3); self.ion_mass.setReadOnly(True)
        self.ion_energy = QDoubleSpinBox(); self.ion_energy.setRange(0.001, 1_000_000); self.ion_energy.setSuffix(" keV"); self.ion_energy.setValue(10.0)
        self.ion_angle = QDoubleSpinBox(); self.ion_angle.setRange(0.0, 90.0); self.ion_angle.setDecimals(1); self.ion_angle.setSuffix(" °"); self.ion_angle.setValue(0.0)

        grid.addWidget(self.pick_btn, 0, 0, 1, 2)
        grid.addWidget(QLabel("Symbol"), 0, 2); grid.addWidget(self.ion_symbol, 0, 3)
        grid.addWidget(QLabel("Name of Element"), 0, 4); grid.addWidget(self.ion_name, 0, 5)
        grid.addWidget(QLabel("Atomic Number"), 0, 6); grid.addWidget(self.ion_z, 0, 7)
        grid.addWidget(QLabel("Mass (amu)"), 0, 8); grid.addWidget(self.ion_mass, 0, 9)
        grid.addWidget(QLabel("Energy (keV)"), 0, 10); grid.addWidget(self.ion_energy, 0, 11)
        grid.addWidget(QLabel("Angle of Incidence"), 0, 12); grid.addWidget(self.ion_angle, 0, 13)

        grid.setColumnStretch(5, 1)
        grid.setColumnStretch(13, 1)
        return box

    def on_element_selected(self, element: dict):
        self.ion_symbol.setText(element.get("symbol", ""))
        self.ion_name.setText(element.get("name", ""))
        try:
            self.ion_z.setValue(int(element.get("number", self.ion_z.value())))
        except (TypeError, ValueError):
            pass
        try:
            self.ion_mass.setValue(float(element.get("atomic_mass", self.ion_mass.value())))
        except (TypeError, ValueError):
            pass

    def build_target_layers(self) -> QGroupBox:
        box = QGroupBox("Target layer selection")
        v = QVBoxLayout(box)

        controls = QHBoxLayout()
        add_layer = QPushButton("Add New Layer")
        del_layer = QPushButton("Delete Selected Layer(s)")
        controls.addWidget(add_layer); controls.addWidget(del_layer); controls.addStretch(1)
        v.addLayout(controls)

        self.layers_table = QTableWidget(1, 6)
        self.layers_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.layers_table.setHorizontalHeaderLabels(["Layer", "Width", "Units", "Density (g/cm³)", "Compound Corr", "Gas"])
        self.layers_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.layers_table.verticalHeader().setVisible(False)
        self.layers_table.setAlternatingRowColors(True)

        self.seed_layer_row(0)
        v.addWidget(self.layers_table)

        self.layer_elements.append([])
        self.layers_table.selectRow(0)
        self.layers_table.itemSelectionChanged.connect(self._handle_layer_selection_changed)

        add_layer.clicked.connect(self.add_layer_row)
        del_layer.clicked.connect(self.delete_selected_layers)
        return box

    def seed_layer_row(self, r: int):
        self.layers_table.setItem(r, 0, QTableWidgetItem(f"Layer {r + 1}"))
        self.layers_table.setItem(r, 1, QTableWidgetItem("10000" if r == 0 else ""))
        unit_combo = QComboBox(); unit_combo.addItems(self.state.unit_options); unit_combo.setCurrentText("Ång")
        self.layers_table.setCellWidget(r, 2, unit_combo)
        self.layers_table.setItem(r, 3, QTableWidgetItem("1.0" if r == 0 else ""))
        self.layers_table.setItem(r, 4, QTableWidgetItem("0"))
        gas_chk = QCheckBox()
        gas_widget = QWidget(); lay = QHBoxLayout(gas_widget); lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(gas_chk); lay.addStretch(1)
        self.layers_table.setCellWidget(r, 5, gas_widget)

    def add_layer_row(self):
        r = self.layers_table.rowCount()
        self.layers_table.insertRow(r)
        self.seed_layer_row(r)
        self.layer_elements.append([])
        self.layers_table.selectRow(r)

    def delete_selected_layers(self):
        rows = sorted({idx.row() for idx in self.layers_table.selectedIndexes()}, reverse=True)
        for r in rows:
            self.layers_table.removeRow(r)
            if 0 <= r < len(self.layer_elements):
                self.layer_elements.pop(r)
        if self.layers_table.rowCount() == 0:
            self.layers_table.insertRow(0)
            self.seed_layer_row(0)
            self.layer_elements = [[]]
        self.layers_table.selectRow(min(self.layers_table.rowCount() - 1, 0))
        self._refresh_element_table()

    def _handle_layer_selection_changed(self):
        self._refresh_element_table()

    def build_input_elements(self) -> QGroupBox:
        box = QGroupBox("Atoms per layer")
        v = QVBoxLayout(box)

        controls = QHBoxLayout()
        self.elem_pick_btn = PeriodicTableButton(
            "Click to Select Element",
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=False,
        )
        self.elem_pick_btn.element_selected.connect(self.on_target_element_selected)

        del_elem = QPushButton("Delete Selected Element(s)")
        dict_btn = QPushButton("Compound Dictionary")
        dict_btn.clicked.connect(self._open_compound_dictionary)

        controls.addWidget(self.elem_pick_btn)
        controls.addWidget(del_elem)
        controls.addStretch(1)
        controls.addWidget(dict_btn)
        v.addLayout(controls)

        self.elem_table = QTableWidget(0, 10)
        self.elem_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.elem_table.setHorizontalHeaderLabels([
            "Symbol", "Name", "Atomic No.", "Weight (amu)",
            "Atom Stoich", "Atom Stoich %", "Damage (eV)", "Disp (eV)", "Latt (eV)", "Surf (eV)"
        ])
        self.elem_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.elem_table.verticalHeader().setVisible(False)
        self.elem_table.setAlternatingRowColors(True)

        v.addWidget(self.elem_table)

        self.elem_table.itemChanged.connect(self._handle_element_item_changed)
        self.elem_table.cellDoubleClicked.connect(self._handle_element_cell_double_clicked)
        del_elem.clicked.connect(self.delete_selected_elements)

        return box

    def on_target_element_selected(self, element: dict):
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        self._add_element_to_layer(layer_idx, element, 1.0)

    def delete_selected_elements(self):
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        rows = sorted({idx.row() for idx in self.elem_table.selectedIndexes()}, reverse=True)
        entries = self._get_layer_entries(layer_idx)
        for r in rows:
            if 0 <= r < len(entries):
                entries.pop(r)
        self._refresh_element_table()

    def build_model_selection(self) -> QGroupBox:
        box = QGroupBox("Model selection")
        h = QHBoxLayout(box)
        self.model_combo = QComboBox()
        self.model_combo.addItems(["Sample Model 1", "Sample Model 2", "Sample Model 3"])
        h.addWidget(self.model_combo, 1)
        h.addStretch(1)
        btn_advanced = QPushButton("Advanced")
        btn_advanced.clicked.connect(lambda: QMessageBox.information(self, "Advanced", "Advanced tab was split out; implement as needed."))
        h.addWidget(btn_advanced)
        return box

    def build_trajectories_output(self) -> QGroupBox:
        box = QGroupBox("Output Options")
        v = QVBoxLayout(box)

        v.addWidget(QLabel("Trajectories:"))
        row = QHBoxLayout(); row.addSpacing(40)
        self.chk_traj_start = QCheckBox("Start"); self.chk_traj_end = QCheckBox("End"); self.chk_traj_coll = QCheckBox("Collisions")
        row.addWidget(self.chk_traj_start); row.addWidget(self.chk_traj_end); row.addWidget(self.chk_traj_coll); row.addStretch(1)
        v.addLayout(row)

        v.addWidget(QLabel("Range Distributions:"))
        row_range = QHBoxLayout(); row_range.addSpacing(40)
        self.chk_range_ion_recoil = QCheckBox("Ion/Recoil"); self.chk_range_phonons = QCheckBox("Phonons"); self.chk_range_ionization = QCheckBox("Ionization")
        row_range.addWidget(self.chk_range_ion_recoil); row_range.addWidget(self.chk_range_phonons); row_range.addWidget(self.chk_range_ionization); row_range.addStretch(1)
        v.addLayout(row_range)

        v.addWidget(QLabel("Lateral Range Distributions:"))
        row_lat = QHBoxLayout(); row_lat.addSpacing(40)
        self.chk_lateral_ion_recoil = QCheckBox("Ion/Recoil"); self.chk_lateral_phonons = QCheckBox("Phonons"); self.chk_lateral_ionization = QCheckBox("Ionization")
        row_lat.addWidget(self.chk_lateral_ion_recoil); row_lat.addWidget(self.chk_lateral_phonons); row_lat.addWidget(self.chk_lateral_ionization); row_lat.addStretch(1)
        v.addLayout(row_lat)

        v.addWidget(QLabel("Backscattered:"))
        row_back = QHBoxLayout(); row_back.addSpacing(40)
        self.chk_backscattered_energy = QCheckBox("Energy"); self.chk_backscattered_angle = QCheckBox("Angle")
        row_back.addWidget(self.chk_backscattered_energy); row_back.addWidget(self.chk_backscattered_angle); row_back.addStretch(1)
        v.addLayout(row_back)

        v.addWidget(QLabel("Transmitted:"))
        row_trans = QHBoxLayout(); row_trans.addSpacing(40)
        self.chk_transmitted_energy = QCheckBox("Energy"); self.chk_transmitted_angle = QCheckBox("Angle")
        row_trans.addWidget(self.chk_transmitted_energy); row_trans.addWidget(self.chk_transmitted_angle); row_trans.addStretch(1)
        v.addLayout(row_trans)

        v.addStretch(1)
        return box

    # -------- footer / logs / progress ----------
    def _build_mc_setup_footer(self) -> QFrame:
        from PyQt6.QtWidgets import QProgressBar  # local import

        footer = QFrame()
        layout = QHBoxLayout(footer)
        layout.setContentsMargins(16, 10, 16, 10)
        layout.setSpacing(12)

        # NEW: integer inputs (left side)
        ions_row = QHBoxLayout()
        ions_row.setSpacing(8)

        ions_row.addWidget(QLabel("No. of Ions"))
        self.no_of_ions_spin = QSpinBox()
        self.no_of_ions_spin.setRange(1, 1_000_000_000)
        self.no_of_ions_spin.setValue(10000)
        ions_row.addWidget(self.no_of_ions_spin)

        ions_row.addSpacing(8)
        ions_row.addWidget(QLabel("Update after Ions"))
        self.update_after_ions_spin = QSpinBox()
        self.update_after_ions_spin.setRange(1, 1_000_000_000)
        self.update_after_ions_spin.setValue(100)
        ions_row.addWidget(self.update_after_ions_spin)

        ions_row.addStretch(1)
        layout.addLayout(ions_row, 2)

        log_btn = QPushButton("No updates yet")
        log_btn.clicked.connect(self._show_logs_dialog)
        self.latest_log_button = log_btn

        self.mc_progress = QProgressBar()
        self.mc_progress.setRange(0, 100)
        self.mc_progress.setValue(0)
        self.mc_progress.setFormat("Ready")

        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._handle_run_clicked)

        layout.addWidget(log_btn, 2)
        layout.addWidget(self.mc_progress, 2)
        layout.addWidget(self.run_button)
        return footer

    def _show_logs_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Update Notifications")
        dialog.resize(800, 400)
        layout = QVBoxLayout(dialog)
        list_widget = QListWidget()

        if self.state.log_entries:
            list_widget.addItems(list(reversed(self.state.log_entries)))
        else:
            list_widget.addItem("No logs available.")

        layout.addWidget(list_widget)

        clear_btn = QPushButton("Clear Logs")
        clear_btn.clicked.connect(lambda: self._clear_logs(list_widget))
        layout.addWidget(clear_btn)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)

        dialog.exec()

    def _clear_logs(self, list_widget: QListWidget):
        self.state.clear_logs()
        list_widget.clear()
        list_widget.addItem("No logs available.")
        if self.latest_log_button:
            self.latest_log_button.setText("No updates yet")
            self.latest_log_button.setToolTip("")

    def _handle_run_clicked(self):
        if not self.mc_progress or not self.run_button:
            return
        self.run_button.setEnabled(False)
        self.mc_progress.setValue(0)
        self.mc_progress.setFormat("%p% Complete")
        if not self._progress_timer:
            self._progress_timer = QTimer(self)
            self._progress_timer.timeout.connect(self._advance_progress)
        self._progress_timer.start(120)
        self.add_log_entry("Simulation run started.")

    def _advance_progress(self):
        if not self.mc_progress:
            return
        value = self.mc_progress.value() + 5
        self.mc_progress.setValue(min(value, 100))
        if self.mc_progress.value() >= 100:
            if self._progress_timer:
                self._progress_timer.stop()
            self.mc_progress.setFormat("Complete")
            if self.run_button:
                self.run_button.setEnabled(True)
            self.add_log_entry("Simulation run completed.")

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
        entries.append({
            "element": element,
            "ratio": ratio_value,
            "damage": energy_defaults["damage"],
            "disp": energy_defaults["disp"],
            "latt": energy_defaults["latt"],
            "surf": energy_defaults["surf"],
        })
        if refresh:
            self._refresh_element_table()

    def _get_default_energy_params(self, element: dict) -> dict:
        params = {}
        for key, fallback in self.state.energy_defaults.items():
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
        self.elem_table.setRowCount(len(entries))

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

            self.elem_table.setItem(row, 0, ro_item(element["symbol"]))
            self.elem_table.setItem(row, 1, ro_item(element["name"]))
            self.elem_table.setItem(row, 2, ro_item(str(element["number"])))
            self.elem_table.setItem(row, 3, ro_item(mass_text))

            ratio_item = QTableWidgetItem(f"{entry['ratio']:.4f}")
            ratio_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
            self.elem_table.setItem(row, 4, ratio_item)

            percent = (entry["ratio"] / total_ratio * 100.0) if total_ratio else 0.0
            self.elem_table.setItem(row, 5, ro_item(f"{percent:.2f}"))

            for offset, key in enumerate(("damage", "disp", "latt", "surf"), start=6):
                self.elem_table.setItem(row, offset, ro_item(str(entry[key])))

        self._updating_elements_table = False

    def _handle_element_item_changed(self, item):
        if self._updating_elements_table:
            return
        layer_idx = self._current_layer_index()
        if layer_idx < 0:
            return
        entries = self._get_layer_entries(layer_idx)
        row = item.row()
        if row < 0 or row >= len(entries):
            return
        entry = entries[row]
        if item.column() == 4:
            try:
                entry["ratio"] = max(float(item.text()), 0.0)
            except ValueError:
                entry["ratio"] = 0.0
            self._refresh_element_table()

    def _handle_element_cell_double_clicked(self, row, column):
        if column != 0:
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

    def _current_layer_index(self):
        if not hasattr(self, "layers_table"):
            return -1
        row = self.layers_table.currentRow()
        if row < 0 and self.layers_table.rowCount() > 0:
            row = 0
            self.layers_table.selectRow(0)
        return row

    def _get_layer_entries(self, layer_idx):
        while len(self.layer_elements) <= layer_idx:
            self.layer_elements.append([])
        return self.layer_elements[layer_idx] if layer_idx >= 0 else []
