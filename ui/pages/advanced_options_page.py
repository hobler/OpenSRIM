from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal, QPropertyAnimation, QEasingCurve
from PyQt6.QtGui import QCursor, QIcon
from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QCheckBox,
    QComboBox,
    QLabel,
    QLineEdit,
    QScrollArea,
    QDoubleSpinBox,
    QSpinBox,
    QFrame,
    QPushButton,
    QStyle,
)


class AccordionItem(QFrame):
    """Accordion item with animated expand/collapse (based on test.py sample)."""

    toggled = pyqtSignal(bool)

    def __init__(self, title: str, content: QWidget, expanded: bool = True, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._expanded = expanded
        self._anim: Optional[QPropertyAnimation] = None

        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setObjectName("accordion-item")
        # Use palette roles (no hard-coded colors) but keep the sample's "modern" layout.
        self.setStyleSheet(
            """
            QFrame#accordion-item {
                border: 1px solid palette(dark);
                border-radius: 10px;
                background: palette(window);
            }
            QPushButton#accordion-header {
                text-align: left;
                padding: 10px 14px;
                border: none;
                font-weight: 600;
                background: palette(alternate-base);
                border-top-left-radius: 10px;
                border-top-right-radius: 10px;
                border-bottom: 1px solid palette(dark);
            }
            QPushButton#accordion-header:hover { background: palette(light); }
            QFrame#accordion-body {
                background: palette(base);
                border-bottom-left-radius: 10px;
                border-bottom-right-radius: 10px;
            }
            """
        )

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        self.headerBtn = QPushButton(title)
        self.headerBtn.setObjectName("accordion-header")
        self.headerBtn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.headerBtn.setIcon(self._icon_for_state(expanded))
        self.headerBtn.clicked.connect(self.toggle)
        lay.addWidget(self.headerBtn)

        self.body = QFrame()
        self.body.setObjectName("accordion-body")
        self.bodyLay = QVBoxLayout(self.body)
        self.bodyLay.setContentsMargins(12, 12, 12, 12)
        self.bodyLay.setSpacing(8)
        self.bodyLay.addWidget(content)
        lay.addWidget(self.body)

        self._apply_state(expanded, animate=False)

    def _icon_for_state(self, expanded: bool) -> QIcon:
        icon = QIcon()
        icon.addPixmap(
            self.style().standardPixmap(
                QStyle.StandardPixmap.SP_ArrowDown if expanded else QStyle.StandardPixmap.SP_ArrowRight
            )
        )
        return icon

    def toggle(self) -> None:
        self._apply_state(not self._expanded, animate=True)

    def set_expanded(self, expanded: bool, *, animate: bool = True) -> None:
        if self._expanded == expanded:
            return
        self._apply_state(expanded, animate=animate)

    def _apply_state(self, expanded: bool, animate: bool) -> None:
        self._expanded = expanded
        self.headerBtn.setIcon(self._icon_for_state(expanded))

        if animate:
            start = int(self.body.maximumHeight())

            if expanded:
                self.body.setMaximumHeight(10**6)
                target = int(self.body.sizeHint().height())
                self.body.setMaximumHeight(start)
                self.body.setVisible(True)
            else:
                target = 0
                self.body.setVisible(True)

            anim = QPropertyAnimation(self.body, b"maximumHeight", self)
            anim.setDuration(160)
            anim.setStartValue(max(start, 0))
            anim.setEndValue(target)
            anim.setEasingCurve(QEasingCurve.Type.InOutCubic)
            anim.finished.connect(lambda: self._on_anim_finished(expanded))
            anim.start()
            self._anim = anim
        else:
            if expanded:
                self.body.setMaximumHeight(16777215)
                self.body.setVisible(True)
            else:
                self.body.setMaximumHeight(0)
                self.body.setVisible(False)

        self.toggled.emit(expanded)

    def _on_anim_finished(self, expanded: bool) -> None:
        if expanded:
            self.body.setMaximumHeight(16777215)
            self.body.setVisible(True)
        else:
            self.body.setMaximumHeight(0)
            self.body.setVisible(False)


class AdvancedOptionsPage(QWidget):
    atoms_columns_visibility_changed = pyqtSignal(bool, bool, bool)
    mc_ion_angle_changed = pyqtSignal(float)
    toolbar_visibility_changed = pyqtSignal(bool)
    columns_changed = pyqtSignal(int)        # 0=auto, 1, 2, 3
    borders_visibility_changed = pyqtSignal(bool)
    plot_font_size_changed = pyqtSignal(float)
    bin_combine_changed = pyqtSignal(int)    # 1..6 adjacent bins summed per tile
    koral_solver_changed = pyqtSignal(dict)
    histogram_settings_changed = pyqtSignal(dict)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)

        root = QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(10)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)

        content = QWidget(scroll)
        content_l = QVBoxLayout(content)
        content_l.setContentsMargins(0, 0, 0, 0)
        content_l.setSpacing(10)

        # --- Ion selection content ---
        ion = QWidget(content)
        ion_l = QVBoxLayout(ion)
        ion_l.setContentsMargins(0, 0, 0, 0)
        ion_l.setSpacing(8)
        ion_l.addWidget(QLabel("Angle of Incidence"))

        ion_l.addWidget(QLabel("MC Setup"))
        self.spin_mc_ion_angle = QDoubleSpinBox()
        self.spin_mc_ion_angle.setRange(0.0, 90.0)
        self.spin_mc_ion_angle.setDecimals(1)
        self.spin_mc_ion_angle.setSuffix(" °")
        self.spin_mc_ion_angle.setValue(0.0)
        ion_l.addWidget(self.spin_mc_ion_angle)
        ion_l.addStretch(1)

        # --- Atoms per layer content ---
        atoms = QWidget(content)
        atoms_l = QVBoxLayout(atoms)
        atoms_l.setContentsMargins(0, 0, 0, 0)
        atoms_l.setSpacing(8)
        atoms_l.addWidget(QLabel("Atoms per layer: Table columns"))
        self.chk_disp = QCheckBox("Show Disp. E. (eV)")
        self.chk_latt = QCheckBox("Show Latt (eV)")
        self.chk_surf = QCheckBox("Show Surf (eV)")
        self.chk_disp.setChecked(True)
        self.chk_latt.setChecked(True)
        self.chk_surf.setChecked(True)
        atoms_l.addWidget(self.chk_disp)
        atoms_l.addWidget(self.chk_latt)
        atoms_l.addWidget(self.chk_surf)
        atoms_l.addStretch(1)

        # --- Model selection content ---
        model = QWidget(content)
        model_l = QVBoxLayout(model)
        model_l.setContentsMargins(0, 0, 0, 0)
        model_l.setSpacing(8)
        nbins_row = QHBoxLayout()
        nbins_row.addWidget(QLabel("Number of bins:"))
        self.spin_nbins = QSpinBox()
        self.spin_nbins.setRange(10, 10000)
        self.spin_nbins.setSingleStep(10)
        self.spin_nbins.setValue(120)
        self.spin_nbins.setToolTip("Number of bins for depth and lateral distributions")
        nbins_row.addWidget(self.spin_nbins)
        nbins_row.addStretch(1)
        model_l.addLayout(nbins_row)
        model_l.addStretch(1)

        # --- Cascade Options content (placeholder) ---
        cascade_opts = QWidget(content)
        cascade_opts_l = QVBoxLayout(cascade_opts)
        cascade_opts_l.setContentsMargins(0, 0, 0, 0)
        cascade_opts_l.setSpacing(8)
        cascade_opts_l.addWidget(QLabel("Cascade options (placeholder)"))
        cascade_opts_l.addStretch(1)

        # --- Nuclear Stopping content (placeholder) ---
        nuclear_opts = QWidget(content)
        nuclear_opts_l = QVBoxLayout(nuclear_opts)
        nuclear_opts_l.setContentsMargins(0, 0, 0, 0)
        nuclear_opts_l.setSpacing(8)
        nuclear_opts_l.addWidget(QLabel("Nuclear stopping options (placeholder)"))
        nuclear_opts_l.addStretch(1)

        # --- Electronic Stopping content (placeholder) ---
        electronic_opts = QWidget(content)
        electronic_opts_l = QVBoxLayout(electronic_opts)
        electronic_opts_l.setContentsMargins(0, 0, 0, 0)
        electronic_opts_l.setSpacing(8)
        electronic_opts_l.addWidget(QLabel("Electronic stopping options (placeholder)"))
        electronic_opts_l.addStretch(1)

        # --- KORAL Solver content ---
        koral_solver = QWidget(content)
        koral_solver_l = QVBoxLayout(koral_solver)
        koral_solver_l.setContentsMargins(0, 0, 0, 0)
        koral_solver_l.setSpacing(8)

        iter_row = QHBoxLayout()
        iter_row.addWidget(QLabel("Iterations:"))
        self.spin_koral_nr_iter = QSpinBox()
        self.spin_koral_nr_iter.setRange(1, 20)
        self.spin_koral_nr_iter.setValue(1)
        self.spin_koral_nr_iter.setToolTip("Number of KORAL iterations")
        iter_row.addWidget(self.spin_koral_nr_iter)
        iter_row.addStretch(1)
        koral_solver_l.addLayout(iter_row)

        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Integration Method:"))
        self.cmb_koral_method = QComboBox()
        self.cmb_koral_method.addItems(["LSODA", "RK45", "RK23", "DOP853", "Radau", "BDF"])
        self.cmb_koral_method.setCurrentIndex(0)
        self.cmb_koral_method.setToolTip("ODE solver method (SciPy solve_ivp)")
        method_row.addWidget(self.cmb_koral_method)
        method_row.addStretch(1)
        koral_solver_l.addLayout(method_row)

        rtol_row = QHBoxLayout()
        rtol_row.addWidget(QLabel("Rel. Tolerance (rtol):"))
        self.le_koral_rtol = QLineEdit("1e-6")
        self.le_koral_rtol.setToolTip("Relative integration tolerance (e.g. 1e-6)")
        rtol_row.addWidget(self.le_koral_rtol)
        koral_solver_l.addLayout(rtol_row)

        atol_row = QHBoxLayout()
        atol_row.addWidget(QLabel("Abs. Tolerance (atol):"))
        self.le_koral_atol = QLineEdit("1e-6")
        self.le_koral_atol.setToolTip("Absolute integration tolerance (e.g. 1e-6)")
        atol_row.addWidget(self.le_koral_atol)
        koral_solver_l.addLayout(atol_row)

        koral_solver_l.addStretch(1)

        # --- Output histogram settings ---
        hist = QWidget(content)
        hist_l = QVBoxLayout(hist)
        hist_l.setContentsMargins(0, 0, 0, 0)
        hist_l.setSpacing(8)

        self.chk_hist_override = QCheckBox("Override histogram limits")
        self.chk_hist_override.setChecked(False)
        hist_l.addWidget(self.chk_hist_override)

        self.spin_hist_nbins = QSpinBox()
        self.spin_hist_nbins.setRange(10, 10000)
        self.spin_hist_nbins.setSingleStep(10)
        self.spin_hist_nbins.setValue(120)
        row_nbins = QHBoxLayout()
        row_nbins.addWidget(QLabel("Bins:"))
        row_nbins.addWidget(self.spin_hist_nbins)
        row_nbins.addStretch(1)
        hist_l.addLayout(row_nbins)

        def _range_row(label: str, default_min: float, default_max: float):
            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            s_min = QDoubleSpinBox()
            s_min.setRange(-1e9, 1e9)
            s_min.setDecimals(4)
            s_min.setValue(default_min)
            s_max = QDoubleSpinBox()
            s_max.setRange(-1e9, 1e9)
            s_max.setDecimals(4)
            s_max.setValue(default_max)
            row.addWidget(QLabel("Min"))
            row.addWidget(s_min)
            row.addWidget(QLabel("Max"))
            row.addWidget(s_max)
            row.addStretch(1)
            hist_l.addLayout(row)
            return s_min, s_max

        self.spin_depth_min, self.spin_depth_max = _range_row("Depth (Å):", 0.0, 4000.0)
        self.spin_lateral_min, self.spin_lateral_max = _range_row("Lateral (Å):", -2000.0, 2000.0)
        self.spin_energy_min, self.spin_energy_max = _range_row("Energy (keV):", 0.0, 1000.0)
        self.spin_angle_min, self.spin_angle_max = _range_row("Angle (deg):", -90.0, 90.0)
        hist_l.addStretch(1)

        # --- Display Settings content (MC Results plots) ---
        display = QWidget(content)
        display_l = QVBoxLayout(display)
        display_l.setContentsMargins(0, 0, 0, 0)
        display_l.setSpacing(8)

        self.chk_toolbars = QCheckBox("Show Plot Toolbars")
        self.chk_toolbars.setChecked(True)
        self.chk_toolbars.toggled.connect(self.toolbar_visibility_changed)
        display_l.addWidget(self.chk_toolbars)

        self.chk_borders = QCheckBox("Show Plot Borders")
        self.chk_borders.setChecked(True)
        self.chk_borders.toggled.connect(self.borders_visibility_changed)
        display_l.addWidget(self.chk_borders)

        font_row = QHBoxLayout()
        font_row.addWidget(QLabel("Plot Font Size:"))
        self.spin_plot_font_size = QDoubleSpinBox()
        self.spin_plot_font_size.setRange(6.0, 30.0)
        self.spin_plot_font_size.setDecimals(1)
        self.spin_plot_font_size.setSingleStep(0.5)
        self.spin_plot_font_size.setSuffix(" pt")
        self.spin_plot_font_size.setValue(10.0)
        self.spin_plot_font_size.valueChanged.connect(self.plot_font_size_changed)
        font_row.addWidget(self.spin_plot_font_size)
        display_l.addLayout(font_row)

        col_row = QHBoxLayout()
        col_row.addWidget(QLabel("Grid Columns:"))
        self.col_combo = QComboBox()
        self.col_combo.addItems(["Auto", "1", "2", "3"])
        self.col_combo.currentIndexChanged.connect(self.columns_changed)
        col_row.addWidget(self.col_combo)
        display_l.addLayout(col_row)

        bin_row = QHBoxLayout()
        bin_row.addWidget(QLabel("Combine Bins:"))
        self.spin_bin_combine = QSpinBox()
        self.spin_bin_combine.setRange(1, 6)
        self.spin_bin_combine.setValue(1)
        self.spin_bin_combine.setToolTip(
            "Sum N adjacent histogram bins per tile (1 = no combining, up to 6)."
        )
        self.spin_bin_combine.valueChanged.connect(self.bin_combine_changed)
        bin_row.addWidget(self.spin_bin_combine)
        bin_row.addStretch(1)
        display_l.addLayout(bin_row)
        display_l.addStretch(1)

        self._acc_ion = AccordionItem("Ion selection", ion, expanded=False, parent=content)
        self._acc_atoms = AccordionItem("Atoms per layer", atoms, expanded=True, parent=content)
        self._acc_model = AccordionItem("Model selection", model, expanded=False, parent=content)
        self._acc_cascade = AccordionItem("Cascade Options", cascade_opts, expanded=False, parent=content)
        self._acc_nuclear = AccordionItem("Nuclear Stopping", nuclear_opts, expanded=False, parent=content)
        self._acc_electronic = AccordionItem("Electronic Stopping", electronic_opts, expanded=False, parent=content)
        self._acc_koral_solver = AccordionItem("KORAL Solver", koral_solver, expanded=False, parent=content)
        self._acc_hist = AccordionItem("Output Histograms", hist, expanded=False, parent=content)
        self._acc_display = AccordionItem("Display Settings", display, expanded=False, parent=content)

        self._accordion_by_id = {
            "ion_selection_mc": self._acc_ion,
            "atoms_per_layer": self._acc_atoms,
            "model_selection": self._acc_model,
            "cascade_options": self._acc_cascade,
            "nuclear_stopping": self._acc_nuclear,
            "electronic_stopping": self._acc_electronic,
            "koral_solver": self._acc_koral_solver,
            "histogram_settings": self._acc_hist,
            "display_settings": self._acc_display,
        }

        self._all_accordions = (
            self._acc_ion, self._acc_atoms, self._acc_model,
            self._acc_cascade, self._acc_nuclear, self._acc_electronic,
            self._acc_koral_solver, self._acc_hist, self._acc_display,
        )
        for item in self._all_accordions:
            item.toggled.connect(lambda on, src=item: self._handle_item_toggled(src, on))

        content_l.addWidget(self._acc_ion)
        content_l.addWidget(self._acc_atoms)
        content_l.addWidget(self._acc_model)
        content_l.addWidget(self._acc_cascade)
        content_l.addWidget(self._acc_nuclear)
        content_l.addWidget(self._acc_electronic)
        content_l.addWidget(self._acc_koral_solver)
        content_l.addWidget(self._acc_hist)
        content_l.addWidget(self._acc_display)
        content_l.addStretch(1)
        content.setLayout(content_l)
        scroll.setWidget(content)
        root.addWidget(scroll)

        self.chk_disp.toggled.connect(self._emit_atoms_visibility)
        self.chk_latt.toggled.connect(self._emit_atoms_visibility)
        self.chk_surf.toggled.connect(self._emit_atoms_visibility)
        self.spin_mc_ion_angle.valueChanged.connect(self._emit_mc_ion_angle)
        self.spin_koral_nr_iter.valueChanged.connect(self._emit_koral_solver)
        self.cmb_koral_method.currentIndexChanged.connect(self._emit_koral_solver)
        self.le_koral_rtol.textChanged.connect(self._emit_koral_solver)
        self.le_koral_atol.textChanged.connect(self._emit_koral_solver)
        self.chk_hist_override.toggled.connect(self._emit_histogram_settings)
        self.spin_hist_nbins.valueChanged.connect(self._emit_histogram_settings)
        for spin in (
            self.spin_depth_min,
            self.spin_depth_max,
            self.spin_lateral_min,
            self.spin_lateral_max,
            self.spin_energy_min,
            self.spin_energy_max,
            self.spin_angle_min,
            self.spin_angle_max,
        ):
            spin.valueChanged.connect(self._emit_histogram_settings)

    def _handle_item_toggled(self, source: AccordionItem, expanded: bool) -> None:
        if not expanded:
            return
        for item in self._all_accordions:
            if item is not source:
                item.set_expanded(False, animate=True)

    def _emit_atoms_visibility(self) -> None:
        self.atoms_columns_visibility_changed.emit(
            bool(self.chk_disp.isChecked()),
            bool(self.chk_latt.isChecked()),
            bool(self.chk_surf.isChecked()),
        )

    def _emit_mc_ion_angle(self, value: float) -> None:
        self.mc_ion_angle_changed.emit(float(value))

    def _emit_koral_solver(self, *_args) -> None:
        self.koral_solver_changed.emit(self._collect_koral_solver())

    def _collect_histogram_settings(self) -> dict:
        return {
            "enabled": bool(self.chk_hist_override.isChecked()),
            "nbins": int(self.spin_hist_nbins.value()),
            "depth_min": float(self.spin_depth_min.value()),
            "depth_max": float(self.spin_depth_max.value()),
            "lateral_min": float(self.spin_lateral_min.value()),
            "lateral_max": float(self.spin_lateral_max.value()),
            "energy_min": float(self.spin_energy_min.value()),
            "energy_max": float(self.spin_energy_max.value()),
            "angle_min": float(self.spin_angle_min.value()),
            "angle_max": float(self.spin_angle_max.value()),
        }

    def _emit_histogram_settings(self, *_args) -> None:
        self.histogram_settings_changed.emit(self._collect_histogram_settings())

    # Backwards-compat helper: treat set_ion_angle as the MC angle.
    def set_ion_angle(self, value: float) -> None:
        self.set_mc_ion_angle(value)

    def set_mc_ion_angle(self, value: float) -> None:
        self.spin_mc_ion_angle.blockSignals(True)
        try:
            self.spin_mc_ion_angle.setValue(float(value))
        finally:
            self.spin_mc_ion_angle.blockSignals(False)

    def _collect_koral_solver(self) -> dict:
        return {
            "nr_iterations": int(self.spin_koral_nr_iter.value()),
            "integration_method": str(self.cmb_koral_method.currentText()),
            "rtol": str(self.le_koral_rtol.text()),
            "atol": str(self.le_koral_atol.text()),
        }

    def collect_config(self) -> dict:
        return {
            "ion": {
                "mc_angle": float(self.spin_mc_ion_angle.value()),
            },
            "atoms_per_layer": {
                "show_disp": bool(self.chk_disp.isChecked()),
                "show_latt": bool(self.chk_latt.isChecked()),
                "show_surf": bool(self.chk_surf.isChecked()),
            },
            "model_selection": {
                "nbins": int(self.spin_nbins.value()),
            },
            "koral_solver": self._collect_koral_solver(),
            "histogram_settings": self._collect_histogram_settings(),
            "display_settings": {
                "show_toolbars": bool(self.chk_toolbars.isChecked()),
                "show_borders": bool(self.chk_borders.isChecked()),
                "plot_font_size": float(self.spin_plot_font_size.value()),
                "grid_columns": int(self.col_combo.currentIndex()),
            },
        }

    def apply_config(self, payload: dict) -> None:
        if not isinstance(payload, dict):
            return
        ion = payload.get("ion") or {}
        if isinstance(ion, dict):
            if "mc_angle" in ion:
                try:
                    self.set_mc_ion_angle(float(ion.get("mc_angle", 0.0)))
                except (TypeError, ValueError):
                    pass

            # Backwards-compat: older configs used ion.angle (apply to MC).
            if ("mc_angle" not in ion) and ("angle" in ion):
                try:
                    v = float(ion.get("angle", 0.0))
                except (TypeError, ValueError):
                    v = None
                if v is not None:
                    self.set_mc_ion_angle(v)

        atoms = payload.get("atoms_per_layer") or {}
        if isinstance(atoms, dict):
            # Keep signals enabled so UI reacts (column visibility etc.)
            self.chk_disp.setChecked(bool(atoms.get("show_disp", self.chk_disp.isChecked())))
            self.chk_latt.setChecked(bool(atoms.get("show_latt", self.chk_latt.isChecked())))
            self.chk_surf.setChecked(bool(atoms.get("show_surf", self.chk_surf.isChecked())))

        model_sel = payload.get("model_selection") or {}
        if isinstance(model_sel, dict) and "nbins" in model_sel:
            try:
                self.spin_nbins.setValue(int(model_sel["nbins"]))
            except (TypeError, ValueError):
                pass

        ks = payload.get("koral_solver") or {}
        if isinstance(ks, dict):
            if "nr_iterations" in ks:
                try:
                    self.spin_koral_nr_iter.setValue(int(ks["nr_iterations"]))
                except (TypeError, ValueError):
                    pass
            if "integration_method" in ks:
                idx = self.cmb_koral_method.findText(str(ks["integration_method"]))
                if idx >= 0:
                    self.cmb_koral_method.setCurrentIndex(idx)
            if "rtol" in ks:
                self.le_koral_rtol.setText(str(ks["rtol"]))
            if "atol" in ks:
                self.le_koral_atol.setText(str(ks["atol"]))

        hist = payload.get("histogram_settings") or {}
        if isinstance(hist, dict):
            self.chk_hist_override.setChecked(bool(hist.get("enabled", self.chk_hist_override.isChecked())))
            if "nbins" in hist:
                try:
                    self.spin_hist_nbins.setValue(int(hist["nbins"]))
                except (TypeError, ValueError):
                    pass
            for key, spin in (
                ("depth_min", self.spin_depth_min),
                ("depth_max", self.spin_depth_max),
                ("lateral_min", self.spin_lateral_min),
                ("lateral_max", self.spin_lateral_max),
                ("energy_min", self.spin_energy_min),
                ("energy_max", self.spin_energy_max),
                ("angle_min", self.spin_angle_min),
                ("angle_max", self.spin_angle_max),
            ):
                if key in hist:
                    try:
                        spin.setValue(float(hist[key]))
                    except (TypeError, ValueError):
                        pass

        disp = payload.get("display_settings") or {}
        if isinstance(disp, dict):
            self.chk_toolbars.setChecked(bool(disp.get("show_toolbars", self.chk_toolbars.isChecked())))
            self.chk_borders.setChecked(bool(disp.get("show_borders", self.chk_borders.isChecked())))
            if "plot_font_size" in disp:
                try:
                    self.spin_plot_font_size.setValue(float(disp.get("plot_font_size", self.spin_plot_font_size.value())))
                except (TypeError, ValueError):
                    pass
            idx = int(disp.get("grid_columns", self.col_combo.currentIndex()))
            if 0 <= idx < self.col_combo.count():
                self.col_combo.setCurrentIndex(idx)

    def open_section(self, section_id: str) -> None:
        item = self._accordion_by_id.get(section_id)
        if not item:
            return
        item.set_expanded(True, animate=True)
