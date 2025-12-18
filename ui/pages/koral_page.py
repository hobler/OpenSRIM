from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QLabel,
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTableWidget, QTableWidgetItem,
    QHeaderView, QStackedWidget, QCheckBox, QPushButton, QMessageBox, QDialog,
    QSplitter, QScrollArea,
    QSizePolicy,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from app.state import AppState
from app.ui.widgets.toggle_switch import ToggleSwitch
from app.ui.widgets.periodic_table_picker import PeriodicTableButton, PeriodicTableDialog
from app.ui.dialogs.compound_dictionary_dialog import CompoundDictionaryDialog


class KoralPage(QWidget):
    def __init__(self, state: AppState, parent=None):
        super().__init__(parent)
        self.state = state

        # NOTE: target layers removed -> keep a flat list of element entries
        self.element_entries: list[dict] = []

        self.layer_elements: list[list[dict]] = []  # ...existing code (now unused)...
        self._updating_elements_table = False

        # Plot fullscreen toggle bookkeeping
        self._plot_fullscreen_dialog: Optional[QDialog] = None
        self._plot_widget: Optional[QWidget] = None
        self._plot_grid: Optional(QGridLayout) = None
        self._plot_placeholder: Optional[QWidget] = None

        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        ion_box = self.build_ion_data()
        ion_box.setTitle("Ion Selection")
        layout.addWidget(ion_box)

        middle = QSplitter(Qt.Orientation.Vertical)  # Use QSplitter for adjustable height
        middle.addWidget(self.build_input_elements())  # Atoms per layer
        middle.addWidget(self._build_koral_bottom_section())  # Ausgabe Optionen / Range / Straggling
        layout.addWidget(middle)

        self._refresh_element_table()

    def _build_koral_bottom_section(self) -> QSplitter:
        splitter = QSplitter(Qt.Orientation.Horizontal)  # Use QSplitter for horizontal resizing
        splitter.addWidget(self._build_model_selection())  # Add "Model Selection" as its own column
        splitter.addWidget(self._build_koral_left_options())  # Add "Ausgabe Optionen"
        splitter.addWidget(self._build_koral_plot_list_section())  # Add "Range / Straggling"
        return splitter

    def _build_model_selection(self) -> QGroupBox:
        box = QGroupBox("Model Selection")
        layout = QVBoxLayout(box)
        self.model_button = QPushButton("Select Models")
        self.model_button.clicked.connect(self._open_model_selection_dialog)
        self.selected_models_label = QLabel("Selected Models:\nNone")
        self.selected_models_label.setWordWrap(True)

        layout.addWidget(self.model_button)
        layout.addWidget(self.selected_models_label)
        return box

    def _open_model_selection_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Select Models")
        dialog.setModal(True)
        dialog.resize(300, 300)  # Adjust height to display all models

        layout = QVBoxLayout(dialog)
        scroll_area = QScrollArea(dialog)
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        # Example model options
        self.model_checkboxes = []
        for model in ["Model A", "Model B", "Model C", "Model D"]:
            checkbox = QCheckBox(model)
            self.model_checkboxes.append(checkbox)
            scroll_layout.addWidget(checkbox)

        scroll_content.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area)

        button_box = QHBoxLayout()
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(lambda: self._update_selected_models(dialog))
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(dialog.reject)
        button_box.addWidget(ok_button)
        button_box.addWidget(cancel_button)
        layout.addLayout(button_box)

        dialog.setLayout(layout)
        dialog.exec()

    def _update_selected_models(self, dialog: QDialog):
        selected_models = [cb.text() for cb in self.model_checkboxes if cb.isChecked()]
        if selected_models:
            self.selected_models_label.setText("Selected Models:\n" + "\n".join(f"• {model}" for model in selected_models))
        else:
            self.selected_models_label.setText("Selected Models:\nNone")
        dialog.accept()

    # --- identical to MCSetup: Ion Selection ---
    def build_ion_data(self) -> QGroupBox:
        box = QGroupBox("Ion Data")
        layout = QVBoxLayout(box)
        layout.setSpacing(10)

        grid = QGridLayout()
        grid.setHorizontalSpacing(6)   # was 10
        grid.setVerticalSpacing(6)     # was 10
        # Optional: slightly tighter margins inside the box
        # grid.setContentsMargins(0, 0, 0, 0)

        # Add the "Click to Select Element" button as the most left element
        self.pick_btn = PeriodicTableButton(
            "Click to Select Element",
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=True,
        )
        self.pick_btn.element_selected.connect(self.on_element_selected)
        grid.addWidget(self.pick_btn, 1, 0, Qt.AlignmentFlag.AlignLeft)

        # Symbol
        grid.addWidget(QLabel("Symbol"), 0, 1, Qt.AlignmentFlag.AlignLeft)
        self.ion_symbol = QLineEdit()
        self.ion_symbol.setReadOnly(True)
        self.ion_symbol.setPlaceholderText("Symbol")
        self.ion_symbol.setMaximumWidth(90)
        self.ion_symbol.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_symbol, 1, 1, Qt.AlignmentFlag.AlignLeft)

        # Name of Element
        grid.addWidget(QLabel("Name of Element"), 0, 2, Qt.AlignmentFlag.AlignLeft)
        self.ion_name = QLineEdit()
        self.ion_name.setReadOnly(True)
        self.ion_name.setPlaceholderText("Element Name")
        self.ion_name.setMaximumWidth(180)
        self.ion_name.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_name, 1, 2, Qt.AlignmentFlag.AlignLeft)

        # Atomic Number
        grid.addWidget(QLabel("Atomic Number"), 0, 3, Qt.AlignmentFlag.AlignLeft)
        self.ion_z = QSpinBox()
        self.ion_z.setRange(1, 120)
        self.ion_z.setReadOnly(True)
        self.ion_z.setFixedWidth(120)  # ensures it is actually wider (maxWidth alone doesn't grow it)
        self.ion_z.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_z, 1, 3, Qt.AlignmentFlag.AlignLeft)

        # Mass (amu)
        grid.addWidget(QLabel("Mass (amu)"), 0, 4, Qt.AlignmentFlag.AlignLeft)
        self.ion_mass = QDoubleSpinBox()
        self.ion_mass.setRange(0.01, 1000.0)
        self.ion_mass.setDecimals(3)
        self.ion_mass.setReadOnly(False)
        self.ion_mass.setMaximumWidth(120)
        self.ion_mass.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_mass, 1, 4, Qt.AlignmentFlag.AlignLeft)

        # Energy min/max (keV)
        grid.addWidget(QLabel("Energy min (keV)"), 0, 5, Qt.AlignmentFlag.AlignLeft)
        self.energy_min = QDoubleSpinBox()
        self.energy_min.setRange(0.0, 1e6)
        self.energy_min.setDecimals(2)
        self.energy_min.setMaximumWidth(130)
        self.energy_min.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.energy_min, 1, 5, Qt.AlignmentFlag.AlignLeft)

        grid.addWidget(QLabel("Energy max (keV)"), 0, 6, Qt.AlignmentFlag.AlignLeft)
        self.energy_max = QDoubleSpinBox()
        self.energy_max.setRange(0.0, 1e6)
        self.energy_max.setDecimals(2)
        self.energy_max.setMaximumWidth(130)
        self.energy_max.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.energy_max, 1, 6, Qt.AlignmentFlag.AlignLeft)

        # Force extra horizontal space to the far right (keeps fields packed on the left)
        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 0)
        grid.setColumnStretch(2, 0)
        grid.setColumnStretch(3, 0)
        grid.setColumnStretch(4, 0)
        grid.setColumnStretch(5, 0)
        grid.setColumnStretch(6, 0)
        grid.setColumnStretch(7, 1)  # trailing stretch column

        layout.addLayout(grid)
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

    def build_input_elements(self) -> QGroupBox:
        box = QGroupBox("Atoms per layer")
        v = QVBoxLayout(box)

        controls = QHBoxLayout()
        self.elem_pick_btn = PeriodicTableButton(
            "+",
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=False,
        )
        self.elem_pick_btn.element_selected.connect(self.on_target_element_selected)

        del_elem = QPushButton("-")
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

        # Hide the last three columns
        self.elem_table.setColumnHidden(7, True)  # "Disp (eV)"
        self.elem_table.setColumnHidden(8, True)  # "Latt (eV)"
        self.elem_table.setColumnHidden(9, True)  # "Surf (eV)"

        v.addWidget(self.elem_table)

        self.elem_table.itemChanged.connect(self._handle_element_item_changed)
        self.elem_table.cellDoubleClicked.connect(self._handle_element_cell_double_clicked)
        del_elem.clicked.connect(self.delete_selected_elements)

        return box

    def on_target_element_selected(self, element: dict):
        self._add_element_to_table(element, 1.0)

    def delete_selected_elements(self):
        if not hasattr(self, "elem_table"):
            return
        rows = sorted({idx.row() for idx in self.elem_table.selectedIndexes()}, reverse=True)
        for r in rows:
            if 0 <= r < len(self.element_entries):
                self.element_entries.pop(r)
        self._refresh_element_table()

    def _open_compound_dictionary(self):
        dialog = CompoundDictionaryDialog(self, editable=True)
        dialog.compound_selected.connect(self._add_compound_to_table)
        dialog.exec()

    def _add_compound_to_table(self, compound: dict):
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
            self._add_element_to_table(element, float(part["fraction"]), refresh=False)
        self._refresh_element_table()

    def _add_element_to_table(self, element: dict, ratio: float, overrides: Optional[dict] = None, refresh: bool = True):
        energy_defaults = self._get_default_energy_params(element)
        if overrides:
            for key, value in overrides.items():
                if value is not None:
                    energy_defaults[key] = str(value)

        try:
            ratio_value = float(ratio)
        except (TypeError, ValueError):
            ratio_value = 0.0

        self.element_entries.append({
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

        entries = self.element_entries
        self._updating_elements_table = True

        # normalize + fill defaults
        for entry in entries:
            ratio_src = entry.get("ratio", entry.get("stoich", 0.0) or 0.0)
            try:
                entry["ratio"] = float(ratio_src)
            except (TypeError, ValueError):
                entry["ratio"] = 0.0
            defaults = self._get_default_energy_params(entry["element"])
            for key in ("damage", "disp", "latt", "surf"):
                entry.setdefault(key, defaults[key])

        total_ratio = sum(e["ratio"] for e in entries)
        self.elem_table.setRowCount(len(entries))

        def ro_item(text: str):
            it = QTableWidgetItem(text)
            it.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            return it

        for row, entry in enumerate(entries):
            element = entry["element"]
            mass_raw = element.get("atomic_mass")
            try:
                mass_text = f"{float(mass_raw):.3f}"
            except (TypeError, ValueError):
                mass_text = str(mass_raw)

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
        row = item.row()
        if not (0 <= row < len(self.element_entries)):
            return

        if item.column() == 4:
            try:
                self.element_entries[row]["ratio"] = max(float(item.text()), 0.0)
            except ValueError:
                self.element_entries[row]["ratio"] = 0.0
            self._refresh_element_table()

    def _handle_element_cell_double_clicked(self, row, column):
        if column != 0:
            return
        if not (0 <= row < len(self.element_entries)):
            return
        dialog = PeriodicTableDialog(self, compact=True, show_hover_info=True, bordered=True)
        dialog.element_selected.connect(lambda element, r=row: self._replace_element_row(r, element))
        dialog.exec()

    def _replace_element_row(self, row: int, element: dict):
        if not (0 <= row < len(self.element_entries)):
            return
        self.element_entries[row]["element"] = element
        self.element_entries[row].update(self._get_default_energy_params(element))
        self._refresh_element_table()

    # --- existing KORAL-specific options/plot section ---
    def _build_koral_left_options(self) -> QGroupBox:
        box = QGroupBox("Ausgabe Optionen")
        v = QVBoxLayout(box)
        v.setSpacing(6)

        self.all_none_btn = QPushButton("All")
        self.all_none_btn.clicked.connect(self._toggle_all_options)

        row_prange = QHBoxLayout()
        self.chk_prange = QCheckBox("Projectile Range")
        self.cmb_prange = QComboBox()
        self.cmb_prange.clear()
        self.cmb_prange.addItems(self.state.unit_options)
        row_prange.addWidget(self.chk_prange)
        row_prange.addStretch(1)
        row_prange.addWidget(self.cmb_prange)
        v.addLayout(row_prange)

        self.chk_long_strag = QCheckBox("Long. Straggling")
        v.addWidget(self.chk_long_strag)
        self.chk_lat_strag = QCheckBox("Lat. Straggling")
        v.addWidget(self.chk_lat_strag)

        row_nucl = QHBoxLayout()
        self.chk_nucl_strag = QCheckBox("Nuclear Stopping")
        self.cmb_nucl_unit = QComboBox()
        self.cmb_nucl_unit.addItems(self.state.unit_options)
        row_nucl.addWidget(self.chk_nucl_strag)
        row_nucl.addStretch(1)
        row_nucl.addWidget(self.cmb_nucl_unit)
        v.addLayout(row_nucl)

        row_elect = QHBoxLayout()
        self.chk_elec_hop = QCheckBox("Electron Stopping")
        self.cmb_elect_unit = QComboBox()
        self.cmb_elect_unit.addItems(self.state.unit_options)
        row_elect.addWidget(self.chk_elec_hop)
        row_elect.addStretch(1)
        row_elect.addWidget(self.cmb_elect_unit)
        v.addLayout(row_elect)

        row_switch = QHBoxLayout()
        row_switch.addWidget(QLabel("Plot"))
        self.sw_koral_mode = ToggleSwitch()
        self.sw_koral_mode.setChecked(False)  # False = Plot, True = List
        row_switch.addWidget(self.sw_koral_mode)
        row_switch.addWidget(QLabel("List"))
        row_switch.addStretch(1)
        v.addLayout(row_switch)

        # Add the "All/None" button at the bottom
        v.addWidget(self.all_none_btn)

        self.sw_koral_mode.toggled.connect(self._update_koral_plot_view)

        v.addStretch(1)
        return box

    def _toggle_all_options(self):
        # Toggle all checkboxes between checked and unchecked
        all_checked = all(
            checkbox.isChecked()
            for checkbox in [self.chk_prange, self.chk_long_strag, self.chk_lat_strag, self.chk_nucl_strag, self.chk_elec_hop]
        )
        new_state = not all_checked

        # Update all checkboxes
        for checkbox in [self.chk_prange, self.chk_long_strag, self.chk_lat_strag, self.chk_nucl_strag, self.chk_elec_hop]:
            checkbox.setChecked(new_state)

        # Update button text
        self.all_none_btn.setText("None" if new_state else "All")

    def _build_koral_plot_list_section(self) -> QGroupBox:
        box = QGroupBox("Range / Straggling vs Energy")
        v = QVBoxLayout(box)

        self.koral_plot_list_stack = QStackedWidget()
        v.addWidget(self.koral_plot_list_stack)

        # Plot page
        plot_page = QWidget()
        grid = QGridLayout(plot_page)
        grid.setContentsMargins(4, 4, 4, 4)

        # Create the Matplotlib figure and canvas (+ toolbar)
        self.figure = Figure(figsize=(5, 3))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        ax = self.figure.add_subplot(111)
        ax.set_title("Range and Straggling vs Energy")
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Range / Straggling (µm)")
        energy = [10, 20, 30]
        ax.plot(energy, [0.5, 1.0, 1.5], label="Range", marker="o")
        ax.plot(energy, [0.1, 0.2, 0.3], label="Straggling", marker="s")
        ax.legend()

        # Wrap toolbar+canvas into a single movable widget
        self._plot_grid = grid
        self._plot_widget = QWidget()
        plot_v = QVBoxLayout(self._plot_widget)
        plot_v.setContentsMargins(0, 0, 0, 0)
        plot_v.setSpacing(2)
        plot_v.addWidget(self.toolbar)
        plot_v.addWidget(self.canvas)

        # Placeholder used while plot is re-parented to fullscreen dialog
        self._plot_placeholder = QWidget()
        self._plot_placeholder.setMinimumSize(1, 1)

        # Connect double-click event for fullscreen toggle
        self.canvas.mpl_connect("button_press_event", self._handle_plot_double_click)

        # Add the plot widget to the grid layout
        grid.addWidget(self._plot_widget, 0, 0, 1, 3)  # Span all columns

        grid.setColumnStretch(1, 1)
        self.koral_plot_list_stack.addWidget(plot_page)

        # List page
        list_page = QWidget()
        lv = QVBoxLayout(list_page)
        self.koral_result_table = QTableWidget(3, 3)
        self.koral_result_table.setHorizontalHeaderLabels(["Energy (keV)", "Range (µm)", "Straggling (µm)"])
        self.koral_result_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.koral_result_table.verticalHeader().setVisible(False)
        self.koral_result_table.setAlternatingRowColors(True)

        for row, (e, r, s) in enumerate([("10", "0.5", "0.1"), ("20", "1.0", "0.2"), ("30", "1.5", "0.3")]):
            self.koral_result_table.setItem(row, 0, QTableWidgetItem(e))
            self.koral_result_table.setItem(row, 1, QTableWidgetItem(r))
            self.koral_result_table.setItem(row, 2, QTableWidgetItem(s))

        lv.addWidget(self.koral_result_table)
        self.koral_plot_list_stack.addWidget(list_page)

        self.koral_plot_list_stack.setCurrentIndex(0)
        return box

    def _handle_plot_double_click(self, event):
        if getattr(event, "dblclick", False):
            self._toggle_plot_fullscreen()

    def _toggle_plot_fullscreen(self):
        if not self._plot_widget or not self._plot_grid:
            return

        # If currently fullscreen -> restore
        if self._plot_fullscreen_dialog and self._plot_fullscreen_dialog.isVisible():
            self._plot_fullscreen_dialog.close()
            return

        # Move plot widget into a fullscreen dialog
        dlg = QDialog(self)
        dlg.setWindowTitle("Plot")
        dlg.setModal(False)

        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        # Keep grid cell occupied while plot is detached
        self._plot_grid.removeWidget(self._plot_widget)
        if self._plot_placeholder:
            self._plot_grid.addWidget(self._plot_placeholder, 1, 1)

        self._plot_widget.setParent(dlg)
        lay.addWidget(self._plot_widget)

        def _restore():
            if not self._plot_widget or not self._plot_grid:
                return
            # Remove placeholder and put plot widget back into its original cell
            if self._plot_placeholder:
                self._plot_grid.removeWidget(self._plot_placeholder)
                self._plot_placeholder.setParent(None)
            self._plot_widget.setParent(self._plot_grid.parentWidget())
            self._plot_grid.addWidget(self._plot_widget, 1, 1)

        dlg.finished.connect(_restore)
        self._plot_fullscreen_dialog = dlg
        dlg.showFullScreen()

    def _update_koral_plot_view(self, checked: bool):
        self.koral_plot_list_stack.setCurrentIndex(1 if checked else 0)
