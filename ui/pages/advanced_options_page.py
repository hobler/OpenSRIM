from __future__ import annotations

import json
import os
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
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QAbstractItemView,
    QSizePolicy,
)

try:
    from ui.widgets.periodic_table_picker import PeriodicTableDialog
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.widgets.periodic_table_picker import PeriodicTableDialog  # type: ignore


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
    mc_follow_recoils_changed = pyqtSignal(bool)
    mc_rng_seed_changed = pyqtSignal(int)
    mc_electronic_stopping_changed = pyqtSignal(str)
    mc_electronic_straggling_changed = pyqtSignal(str)
    mc_scattering_algorithm_changed = pyqtSignal(str)
    mc_n_absc_changed = pyqtSignal(int)
    mc_lindhard_correction_changed = pyqtSignal(dict)
    mc_pmax_min_changed = pyqtSignal(float)
    mc_pmax_max_changed = pyqtSignal(float)
    mc_psi_min_changed = pyqtSignal(float)
    mc_de_min_changed = pyqtSignal(float)
    mc_psi_min_surface_changed = pyqtSignal(float)
    mc_de_min_surface_changed = pyqtSignal(float)
    mc_replacement_collisions_changed = pyqtSignal(bool)
    mc_nthreads_changed = pyqtSignal(int)
    mc_target_roughness_changed = pyqtSignal(float)
    toolbar_visibility_changed = pyqtSignal(bool)
    columns_changed = pyqtSignal(int)        # 0=auto, 1, 2, 3
    borders_visibility_changed = pyqtSignal(bool)
    plot_font_size_changed = pyqtSignal(float)
    beam_from_top_changed = pyqtSignal(bool)
    koral_solver_changed = pyqtSignal(dict)
    histogram_settings_changed = pyqtSignal(dict)

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)

        # Optional callback (set via set_project_elements_provider) returning
        # the project's current ion + target element symbols, so pickers like
        # the Lindhard-correction "Add Row" can restrict to those elements.
        self._project_elements_provider = None

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

        # --- Target Layers content ---
        target_layers = QWidget(content)
        target_layers_l = QVBoxLayout(target_layers)
        target_layers_l.setContentsMargins(0, 0, 0, 0)
        target_layers_l.setSpacing(8)

        roughness_row = QHBoxLayout()
        roughness_row.addWidget(QLabel("Top layer roughness (Å):"))
        self.spin_target_roughness = QDoubleSpinBox()
        self.spin_target_roughness.setRange(0.0, 1.0e9)
        self.spin_target_roughness.setDecimals(4)
        self.spin_target_roughness.setValue(0.0)
        roughness_row.addWidget(self.spin_target_roughness)
        roughness_row.addStretch(1)
        target_layers_l.addLayout(roughness_row)
        target_layers_l.addStretch(1)

        # --- Cascade Options content ---
        cascade_opts = QWidget(content)
        cascade_opts_l = QVBoxLayout(cascade_opts)
        cascade_opts_l.setContentsMargins(0, 0, 0, 0)
        cascade_opts_l.setSpacing(8)

        self.chk_follow_recoils = QCheckBox("Follow recoils")
        self.chk_follow_recoils.setChecked(True)
        cascade_opts_l.addWidget(self.chk_follow_recoils)

        self.chk_replacement_collisions = QCheckBox("Replacement collisions")
        self.chk_replacement_collisions.setChecked(False)
        self.chk_replacement_collisions.setToolTip(
            "Allow a slow recoil of the same species as the struck lattice atom "
            "to replace it in place, instead of creating a separate interstitial."
        )
        cascade_opts_l.addWidget(self.chk_replacement_collisions)

        def _cascade_param_row(label: str, minimum: float, maximum: float,
                                value: float, decimals: int = 3) -> QDoubleSpinBox:
            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            spin = QDoubleSpinBox()
            spin.setRange(minimum, maximum)
            spin.setDecimals(decimals)
            spin.setValue(value)
            row.addWidget(spin)
            row.addStretch(1)
            cascade_opts_l.addLayout(row)
            return spin

        self.spin_pmax_max = _cascade_param_row("Maximum impact parameter (Å):", 0.0, 1.0e6, 4.0)
        self.spin_pmax_min = _cascade_param_row("Minimum impact parameter (Å):", 0.0, 1.0e6, 0.0)
        self.spin_de_min = _cascade_param_row("Minimum energy transfer (eV):", 0.0, 1.0e9, 15.0)
        self.spin_psi_min = _cascade_param_row("Minimum scattering angle (°):", 0.0, 180.0, 5.0)
        self.spin_psi_min_surface = _cascade_param_row("Minimum scattering angle above surface (°):", 0.0, 180.0, 5.0)
        self.spin_de_min_surface = _cascade_param_row("Minimum energy transfer above surface (eV):", 0.0, 1.0e9, 15.0)

        seed_row = QHBoxLayout()
        seed_row.addWidget(QLabel("RNG Seed:"))
        self.spin_rng_seed = QSpinBox()
        self.spin_rng_seed.setRange(0, 2_147_483_647)
        self.spin_rng_seed.setValue(12345)
        seed_row.addWidget(self.spin_rng_seed)
        seed_row.addStretch(1)
        cascade_opts_l.addLayout(seed_row)

        # Default to ~80% of logical cores: at 100% the simulation can make
        # the whole UI (including the Stop button) sluggish to unresponsive.
        _cpu_count = os.cpu_count() or 1
        _default_threads = max(1, round(0.8 * _cpu_count))
        threads_row = QHBoxLayout()
        threads_row.addWidget(QLabel(f"Simulation threads (of {_cpu_count}):"))
        self.spin_nthreads = QSpinBox()
        self.spin_nthreads.setRange(1, _cpu_count)
        self.spin_nthreads.setValue(_default_threads)
        self.spin_nthreads.setToolTip(
            "Number of CPU threads used by the Monte Carlo simulation. "
            "Defaults to ~80% of available cores so the UI (including Stop) "
            "stays responsive while a simulation is running."
        )
        threads_row.addWidget(self.spin_nthreads)
        threads_row.addStretch(1)
        cascade_opts_l.addLayout(threads_row)
        cascade_opts_l.addStretch(1)

        # --- Nuclear Stopping content ---
        nuclear_opts = QWidget(content)
        nuclear_opts_l = QVBoxLayout(nuclear_opts)
        nuclear_opts_l.setContentsMargins(0, 0, 0, 0)
        nuclear_opts_l.setSpacing(8)

        alg_row = QHBoxLayout()
        alg_row.addWidget(QLabel("Scattering Algorithm:"))
        self.cmb_scattering_algorithm = QComboBox()
        self.cmb_scattering_algorithm.addItems(["magic", "Legendre"])
        self.cmb_scattering_algorithm.setCurrentText("Legendre")
        alg_row.addWidget(self.cmb_scattering_algorithm)
        alg_row.addStretch(1)
        nuclear_opts_l.addLayout(alg_row)

        absc_row = QHBoxLayout()
        absc_row.addWidget(QLabel("n_absc:"))
        self.spin_n_absc = QSpinBox()
        self.spin_n_absc.setRange(1, 1000)
        self.spin_n_absc.setValue(4)
        absc_row.addWidget(self.spin_n_absc)
        absc_row.addStretch(1)
        nuclear_opts_l.addLayout(absc_row)
        nuclear_opts_l.addStretch(1)

        # --- Electronic Stopping content ---
        electronic_opts = QWidget(content)
        electronic_opts_l = QVBoxLayout(electronic_opts)
        electronic_opts_l.setContentsMargins(0, 0, 0, 0)
        electronic_opts_l.setSpacing(8)

        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("Electronic Stopping Model:"))
        self.cmb_electronic_stopping = QComboBox()
        self.cmb_electronic_stopping.addItems(["SRIM", "Lindhard"])
        self.cmb_electronic_stopping.setCurrentText("SRIM")
        model_row.addWidget(self.cmb_electronic_stopping)
        model_row.addStretch(1)
        electronic_opts_l.addLayout(model_row)

        self.chk_electronic_straggling = QCheckBox("Electronic straggling")
        self.chk_electronic_straggling.setChecked(False)
        self.chk_electronic_straggling.setToolTip(
            "When enabled, electronic energy-loss straggling (Bohr model) is applied.\n"
            "Maps to models.electronic_straggling in the simulation TOML "
            '("Off" when disabled, "Bohr" when enabled).'
        )
        electronic_opts_l.addWidget(self.chk_electronic_straggling)

        electronic_opts_l.addWidget(QLabel("Lindhard correction factors:"))
        self.tbl_lindhard_correction = QTableWidget(0, 3)
        self.tbl_lindhard_correction.setHorizontalHeaderLabels(["Projectile", "Target", "Factor"])
        self.tbl_lindhard_correction.setFrameShape(QFrame.Shape.Box)
        self.tbl_lindhard_correction.setFrameShadow(QFrame.Shadow.Plain)
        self.tbl_lindhard_correction.setLineWidth(1)
        self.tbl_lindhard_correction.verticalHeader().setVisible(False)
        self.tbl_lindhard_correction.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tbl_lindhard_correction.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.tbl_lindhard_correction.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        self.tbl_lindhard_correction.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Interactive)
        self.tbl_lindhard_correction.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Interactive)
        self.tbl_lindhard_correction.horizontalHeader().setStretchLastSection(False)
        self.tbl_lindhard_correction.setColumnWidth(0, 120)
        self.tbl_lindhard_correction.setColumnWidth(1, 120)
        self.tbl_lindhard_correction.setColumnWidth(2, 240)
        self.tbl_lindhard_correction.setMinimumWidth(500)
        self.tbl_lindhard_correction.setMaximumWidth(500)
        self.tbl_lindhard_correction.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)
        self.tbl_lindhard_correction.setAlternatingRowColors(True)
        self.tbl_lindhard_correction.setShowGrid(True)
        self.tbl_lindhard_correction.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.tbl_lindhard_correction.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.tbl_lindhard_correction.setStyleSheet(
            "QTableWidget {"
            " border: 1px solid palette(mid);"
            " border-radius: 6px;"
            " background: palette(base);"
            " alternate-background-color: palette(alternate-base);"
            " selection-background-color: palette(highlight);"
            " gridline-color: palette(mid);"
            " }"
            " QHeaderView::section {"
            " background: palette(alternate-base);"
            " border: 1px solid palette(mid);"
            " padding: 6px 8px;"
            " font-weight: 600;"
            " }"
        )
        electronic_opts_l.addWidget(self.tbl_lindhard_correction)

        lindhard_btn_row = QHBoxLayout()
        self.btn_lindhard_add = QPushButton("Add Row")
        self.btn_lindhard_remove = QPushButton("Remove Row")
        lindhard_btn_row.addWidget(self.btn_lindhard_add)
        lindhard_btn_row.addWidget(self.btn_lindhard_remove)
        lindhard_btn_row.addStretch(1)
        electronic_opts_l.addLayout(lindhard_btn_row)
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

        self.chk_hist_override = QCheckBox("Use custom histogram limits for simulation output")
        self.chk_hist_override.setChecked(False)
        self.chk_hist_override.setToolTip(
            "When enabled, the custom bins and min/max ranges below are written into the simulation TOML. "
            "When disabled, histogram limits are derived from the current layer width and ion energy."
        )
        hist_l.addWidget(self.chk_hist_override)

        self.lbl_hist_info = QLabel(
            "These settings affect generated simulation histograms, not the display of existing MC Results plots."
        )
        self.lbl_hist_info.setWordWrap(True)
        self.lbl_hist_info.setStyleSheet("color: palette(mid);")
        hist_l.addWidget(self.lbl_hist_info)

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

        self.chk_beam_from_top = QCheckBox("Show beam from top in 2D plots")
        self.chk_beam_from_top.setChecked(True)
        self.chk_beam_from_top.toggled.connect(self.beam_from_top_changed)
        display_l.addWidget(self.chk_beam_from_top)

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

        display_l.addStretch(1)

        self._acc_ion = AccordionItem("Ion selection", ion, expanded=False, parent=content)
        self._acc_target_layers = AccordionItem("Target Layers", target_layers, expanded=False, parent=content)
        self._acc_atoms = AccordionItem("Atoms per layer", atoms, expanded=True, parent=content)
        self._acc_model = AccordionItem("Model selection", model, expanded=False, parent=content)
        self._acc_cascade = AccordionItem("Cascade Options", cascade_opts, expanded=False, parent=content)
        self._acc_nuclear = AccordionItem("Nuclear Stopping", nuclear_opts, expanded=False, parent=content)
        self._acc_electronic = AccordionItem("Electronic Stopping", electronic_opts, expanded=False, parent=content)
        self._acc_koral_solver = AccordionItem("KORAL Solver", koral_solver, expanded=False, parent=content)
        self._acc_hist = AccordionItem("Simulation Output Histograms", hist, expanded=False, parent=content)
        self._acc_display = AccordionItem("Display Settings", display, expanded=False, parent=content)

        self._accordion_by_id = {
            "ion_selection_mc": self._acc_ion,
            "target_layers": self._acc_target_layers,
            "atoms_per_layer": self._acc_atoms,
            "model_selection": self._acc_model,
            "mc_setup_advanced": self._acc_cascade,
            "cascade_options": self._acc_cascade,
            "nuclear_stopping": self._acc_nuclear,
            "electronic_stopping": self._acc_electronic,
            "koral_solver": self._acc_koral_solver,
            "histogram_settings": self._acc_hist,
            "display_settings": self._acc_display,
        }

        self._all_accordions = (
            self._acc_ion, self._acc_target_layers, self._acc_atoms, self._acc_model,
            self._acc_cascade, self._acc_nuclear, self._acc_electronic,
            self._acc_koral_solver, self._acc_hist, self._acc_display,
        )
        for item in self._all_accordions:
            item.toggled.connect(lambda on, src=item: self._handle_item_toggled(src, on))

        content_l.addWidget(self._acc_ion)
        content_l.addWidget(self._acc_target_layers)
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
        self.chk_follow_recoils.toggled.connect(self._emit_mc_setup_advanced)
        self.chk_replacement_collisions.toggled.connect(self._emit_mc_setup_advanced)
        self.spin_pmax_max.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_pmax_min.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_de_min.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_psi_min.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_psi_min_surface.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_de_min_surface.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_nthreads.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_rng_seed.valueChanged.connect(self._emit_mc_setup_advanced)
        self.spin_target_roughness.valueChanged.connect(self._emit_mc_setup_advanced)
        self.cmb_electronic_stopping.currentIndexChanged.connect(self._emit_mc_setup_advanced)
        self.chk_electronic_straggling.toggled.connect(self._emit_mc_setup_advanced)
        self.cmb_scattering_algorithm.currentIndexChanged.connect(self._emit_mc_setup_advanced)
        self.spin_n_absc.valueChanged.connect(self._emit_mc_setup_advanced)
        self.btn_lindhard_add.clicked.connect(self._on_add_lindhard_row)
        self.btn_lindhard_remove.clicked.connect(self._on_remove_lindhard_row)
        self.chk_hist_override.toggled.connect(self._emit_histogram_settings)
        self.chk_hist_override.toggled.connect(self._update_histogram_override_ui)
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

        self._elements_by_symbol = self._load_elements_lookup()
        self._element_symbols = sorted(self._elements_by_symbol.keys())
        self._updating_lindhard = False
        self._suppress_incoming_mc_setup_apply = False

        # Default examples to avoid empty free-text parsing by users.
        self._apply_lindhard_correction({"B->Si": 1.5, "Si->Si": 1.0})
        self._update_lindhard_table_height()
        self._update_histogram_override_ui(bool(self.chk_hist_override.isChecked()))

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

    def _emit_mc_setup_advanced(self, *_args) -> None:
        self._suppress_incoming_mc_setup_apply = True
        try:
            self.mc_follow_recoils_changed.emit(bool(self.chk_follow_recoils.isChecked()))
            self.mc_replacement_collisions_changed.emit(bool(self.chk_replacement_collisions.isChecked()))
            self.mc_pmax_min_changed.emit(float(self.spin_pmax_min.value()))
            self.mc_pmax_max_changed.emit(float(self.spin_pmax_max.value()))
            self.mc_psi_min_changed.emit(float(self.spin_psi_min.value()))
            self.mc_de_min_changed.emit(float(self.spin_de_min.value()))
            self.mc_psi_min_surface_changed.emit(float(self.spin_psi_min_surface.value()))
            self.mc_de_min_surface_changed.emit(float(self.spin_de_min_surface.value()))
            self.mc_nthreads_changed.emit(int(self.spin_nthreads.value()))
            self.mc_rng_seed_changed.emit(int(self.spin_rng_seed.value()))
            self.mc_target_roughness_changed.emit(float(self.spin_target_roughness.value()))
            self.mc_electronic_stopping_changed.emit(str(self.cmb_electronic_stopping.currentText()))
            self.mc_electronic_straggling_changed.emit(
                "Bohr" if self.chk_electronic_straggling.isChecked() else "Off"
            )
            self.mc_scattering_algorithm_changed.emit(str(self.cmb_scattering_algorithm.currentText()))
            self.mc_n_absc_changed.emit(int(self.spin_n_absc.value()))
            self.mc_lindhard_correction_changed.emit(self._collect_lindhard_correction())
        finally:
            self._suppress_incoming_mc_setup_apply = False

    def _on_lindhard_editor_changed(self, *_args) -> None:
        if self._updating_lindhard:
            return
        self._emit_mc_setup_advanced()

    def _on_add_lindhard_row(self) -> None:
        self._add_lindhard_row("", "", 1.0, emit=False)
        row = self.tbl_lindhard_correction.rowCount() - 1
        if row >= 0:
            self.tbl_lindhard_correction.selectRow(row)

    def _on_remove_lindhard_row(self) -> None:
        row = self.tbl_lindhard_correction.currentRow()
        if row < 0:
            row = self.tbl_lindhard_correction.rowCount() - 1
        if row < 0:
            return
        self.tbl_lindhard_correction.removeRow(row)
        self._update_lindhard_table_height()
        self._emit_mc_setup_advanced()

    def set_project_elements_provider(self, provider) -> None:
        """Register a no-arg callable returning the project's current ion +
        target element symbols (see MCSetupPage.get_project_element_symbols)."""
        self._project_elements_provider = provider

    def _on_lindhard_pick_clicked(self, button: QPushButton) -> None:
        allowed_symbols = None
        if self._project_elements_provider is not None:
            try:
                allowed_symbols = self._project_elements_provider() or None
            except Exception:
                allowed_symbols = None
        dialog = PeriodicTableDialog(self, compact=True, show_hover_info=True, bordered=True,
                                      allowed_symbols=allowed_symbols)

        def _set_selected_element(element: dict) -> None:
            symbol = str(element.get("symbol", "")).strip()
            name = str(element.get("name", "")).strip()
            if not symbol:
                return
            button.setProperty("element_symbol", symbol)
            button.setProperty("element_name", name or symbol)
            button.setText(name or symbol)
            self._emit_mc_setup_advanced()

        dialog.element_selected.connect(_set_selected_element)
        dialog.exec()

    def _load_elements_lookup(self) -> dict[str, str]:
        json_path = os.path.join(os.path.dirname(__file__), "..", "widgets", "PeriodicTableJSON.json")
        json_path = os.path.normpath(json_path)
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except OSError:
            return {}
        except json.JSONDecodeError:
            return {}

        elements = payload.get("elements", []) if isinstance(payload, dict) else []
        lookup: dict[str, str] = {}
        for elem in elements:
            if not isinstance(elem, dict):
                continue
            symbol = str(elem.get("symbol", "")).strip()
            name = str(elem.get("name", "")).strip()
            if symbol:
                lookup[symbol] = name or symbol
        return lookup

    def _create_lindhard_pick_button(self, symbol: str) -> QPushButton:
        btn = QPushButton(self.tbl_lindhard_correction)
        sym = symbol.strip()
        name = self._elements_by_symbol.get(sym, sym)
        btn.setProperty("element_symbol", sym)
        btn.setProperty("element_name", name)
        btn.setText(name if name else "Pick element…")
        btn.clicked.connect(lambda _checked=False, b=btn: self._on_lindhard_pick_clicked(b))
        return btn

    def _create_lindhard_factor_spin(self, factor: float) -> QDoubleSpinBox:
        spin = QDoubleSpinBox(self.tbl_lindhard_correction)
        spin.setRange(0.0, 1000.0)
        spin.setDecimals(6)
        spin.setSingleStep(0.1)
        spin.setValue(float(factor))
        spin.valueChanged.connect(self._on_lindhard_editor_changed)
        return spin

    def _update_lindhard_table_height(self) -> None:
        header_h = self.tbl_lindhard_correction.horizontalHeader().height()
        frame = self.tbl_lindhard_correction.frameWidth() * 2
        rows_h = 0
        for row in range(self.tbl_lindhard_correction.rowCount()):
            rows_h += self.tbl_lindhard_correction.rowHeight(row)
        if self.tbl_lindhard_correction.rowCount() == 0:
            rows_h = self.tbl_lindhard_correction.verticalHeader().defaultSectionSize()
        total = header_h + rows_h + frame + 2
        self.tbl_lindhard_correction.setFixedHeight(total)

    def _add_lindhard_row(self, projectile: str, target: str, factor: float, *, emit: bool = True) -> None:
        row = self.tbl_lindhard_correction.rowCount()
        self.tbl_lindhard_correction.insertRow(row)
        self.tbl_lindhard_correction.setCellWidget(row, 0, self._create_lindhard_pick_button(str(projectile)))
        self.tbl_lindhard_correction.setCellWidget(row, 1, self._create_lindhard_pick_button(str(target)))
        self.tbl_lindhard_correction.setCellWidget(row, 2, self._create_lindhard_factor_spin(float(factor)))
        self._update_lindhard_table_height()
        if emit:
            self._emit_mc_setup_advanced()

    def _collect_lindhard_correction(self) -> dict[str, float]:
        parsed: dict[str, float] = {}
        for row in range(self.tbl_lindhard_correction.rowCount()):
            p_widget = self.tbl_lindhard_correction.cellWidget(row, 0)
            t_widget = self.tbl_lindhard_correction.cellWidget(row, 1)
            f_widget = self.tbl_lindhard_correction.cellWidget(row, 2)

            projectile = str(p_widget.property("element_symbol") or "").strip() if isinstance(p_widget, QPushButton) else ""
            target = str(t_widget.property("element_symbol") or "").strip() if isinstance(t_widget, QPushButton) else ""
            if not projectile or not target:
                continue
            factor = float(f_widget.value()) if isinstance(f_widget, QDoubleSpinBox) else 1.0
            parsed[f"{projectile}->{target}"] = factor
        return parsed

    def _apply_lindhard_correction(self, correction: dict[str, float]) -> None:
        self._updating_lindhard = True
        self.tbl_lindhard_correction.setRowCount(0)
        for key, value in correction.items():
            if "->" not in str(key):
                continue
            projectile, target = str(key).split("->", 1)
            try:
                factor = float(value)
            except (TypeError, ValueError):
                continue
            self._add_lindhard_row(projectile.strip(), target.strip(), factor, emit=False)
        self._update_lindhard_table_height()
        self._updating_lindhard = False

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

    def _update_histogram_override_ui(self, enabled: bool) -> None:
        for widget in (
            self.spin_hist_nbins,
            self.spin_depth_min,
            self.spin_depth_max,
            self.spin_lateral_min,
            self.spin_lateral_max,
            self.spin_energy_min,
            self.spin_energy_max,
            self.spin_angle_min,
            self.spin_angle_max,
        ):
            widget.setEnabled(bool(enabled))

    def apply_histogram_defaults(self, defaults: dict) -> None:
        if not isinstance(defaults, dict):
            return
        if self.chk_hist_override.isChecked():
            return
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
            if key not in defaults:
                continue
            spin.blockSignals(True)
            try:
                spin.setValue(float(defaults[key]))
            except (TypeError, ValueError):
                pass
            finally:
                spin.blockSignals(False)
        if "nbins" in defaults:
            self.spin_hist_nbins.blockSignals(True)
            try:
                self.spin_hist_nbins.setValue(int(defaults["nbins"]))
            except (TypeError, ValueError):
                pass
            finally:
                self.spin_hist_nbins.blockSignals(False)

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

    def _collect_mc_setup_advanced(self) -> dict:
        return {
            "follow_recoils": bool(self.chk_follow_recoils.isChecked()),
            "replacement_collisions": bool(self.chk_replacement_collisions.isChecked()),
            "pmax_min": float(self.spin_pmax_min.value()),
            "pmax_max": float(self.spin_pmax_max.value()),
            "psi_min": float(self.spin_psi_min.value()),
            "de_min": float(self.spin_de_min.value()),
            "psi_min_surface": float(self.spin_psi_min_surface.value()),
            "de_min_surface": float(self.spin_de_min_surface.value()),
            "nthreads": int(self.spin_nthreads.value()),
            "target_roughness": float(self.spin_target_roughness.value()),
            "rng_seed": int(self.spin_rng_seed.value()),
            "electronic_stopping": str(self.cmb_electronic_stopping.currentText()),
            "electronic_straggling": "Bohr" if self.chk_electronic_straggling.isChecked() else "Off",
            "scattering_algorithm": str(self.cmb_scattering_algorithm.currentText()),
            "n_absc": int(self.spin_n_absc.value()),
            "lindhard_correction": self._collect_lindhard_correction(),
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
            "mc_setup_advanced": self._collect_mc_setup_advanced(),
            "koral_solver": self._collect_koral_solver(),
            "histogram_settings": self._collect_histogram_settings(),
            "display_settings": {
                "show_toolbars": bool(self.chk_toolbars.isChecked()),
                "show_borders": bool(self.chk_borders.isChecked()),
                "beam_from_top": bool(self.chk_beam_from_top.isChecked()),
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

        mc_setup_adv = payload.get("mc_setup_advanced") or {}
        if isinstance(mc_setup_adv, dict):
            self.apply_mc_setup_advanced_config(mc_setup_adv)

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
            self.chk_beam_from_top.setChecked(bool(disp.get("beam_from_top", self.chk_beam_from_top.isChecked())))
            if "plot_font_size" in disp:
                try:
                    self.spin_plot_font_size.setValue(float(disp.get("plot_font_size", self.spin_plot_font_size.value())))
                except (TypeError, ValueError):
                    pass
            idx = int(disp.get("grid_columns", self.col_combo.currentIndex()))
            if 0 <= idx < self.col_combo.count():
                self.col_combo.setCurrentIndex(idx)

    def apply_mc_setup_advanced_config(self, payload: dict) -> None:
        if self._suppress_incoming_mc_setup_apply:
            return
        if not isinstance(payload, dict):
            return
        if "follow_recoils" in payload:
            self.chk_follow_recoils.blockSignals(True)
            self.chk_follow_recoils.setChecked(bool(payload.get("follow_recoils", self.chk_follow_recoils.isChecked())))
            self.chk_follow_recoils.blockSignals(False)
        if "replacement_collisions" in payload:
            self.chk_replacement_collisions.blockSignals(True)
            self.chk_replacement_collisions.setChecked(
                bool(payload.get("replacement_collisions", self.chk_replacement_collisions.isChecked())))
            self.chk_replacement_collisions.blockSignals(False)
        for key, widget in (
            ("pmax_min", self.spin_pmax_min),
            ("pmax_max", self.spin_pmax_max),
            ("psi_min", self.spin_psi_min),
            ("de_min", self.spin_de_min),
            ("psi_min_surface", self.spin_psi_min_surface),
            ("de_min_surface", self.spin_de_min_surface),
        ):
            if key in payload:
                try:
                    widget.blockSignals(True)
                    widget.setValue(float(payload[key]))
                except (TypeError, ValueError):
                    pass
                finally:
                    widget.blockSignals(False)
        if "nthreads" in payload:
            try:
                self.spin_nthreads.blockSignals(True)
                self.spin_nthreads.setValue(int(payload["nthreads"]))
            except (TypeError, ValueError):
                pass
            finally:
                self.spin_nthreads.blockSignals(False)
        if "rng_seed" in payload:
            try:
                self.spin_rng_seed.blockSignals(True)
                self.spin_rng_seed.setValue(int(payload["rng_seed"]))
            except (TypeError, ValueError):
                pass
            finally:
                self.spin_rng_seed.blockSignals(False)
        if "target_roughness" in payload:
            try:
                self.spin_target_roughness.blockSignals(True)
                self.spin_target_roughness.setValue(float(payload["target_roughness"]))
            except (TypeError, ValueError):
                pass
            finally:
                self.spin_target_roughness.blockSignals(False)
        if "electronic_stopping" in payload:
            idx = self.cmb_electronic_stopping.findText(str(payload.get("electronic_stopping", self.cmb_electronic_stopping.currentText())))
            if idx >= 0:
                self.cmb_electronic_stopping.blockSignals(True)
                self.cmb_electronic_stopping.setCurrentIndex(idx)
                self.cmb_electronic_stopping.blockSignals(False)
        if "electronic_straggling" in payload:
            raw = payload.get("electronic_straggling")
            if isinstance(raw, bool):
                enabled = raw
            else:
                enabled = str(raw).strip().lower() not in ("", "off", "none", "false", "0")
            self.chk_electronic_straggling.blockSignals(True)
            self.chk_electronic_straggling.setChecked(enabled)
            self.chk_electronic_straggling.blockSignals(False)
        if "scattering_algorithm" in payload:
            idx = self.cmb_scattering_algorithm.findText(str(payload.get("scattering_algorithm", self.cmb_scattering_algorithm.currentText())))
            if idx >= 0:
                self.cmb_scattering_algorithm.blockSignals(True)
                self.cmb_scattering_algorithm.setCurrentIndex(idx)
                self.cmb_scattering_algorithm.blockSignals(False)
        if "n_absc" in payload:
            try:
                self.spin_n_absc.blockSignals(True)
                self.spin_n_absc.setValue(int(payload["n_absc"]))
            except (TypeError, ValueError):
                pass
            finally:
                self.spin_n_absc.blockSignals(False)
        if "lindhard_correction" in payload and isinstance(payload.get("lindhard_correction"), dict):
            self._apply_lindhard_correction(payload.get("lindhard_correction") or {})

    def open_section(self, section_id: str) -> None:
        item = self._accordion_by_id.get(section_id)
        if not item:
            return
        item.set_expanded(True, animate=True)

    # Which page each advanced section belongs to. Used to show only the
    # relevant sections depending on where the user opened Advanced Settings.
    _SECTION_CONTEXT = {
        "ion_selection_mc": "mc_setup",
        "target_layers": "mc_setup",
        "atoms_per_layer": "mc_setup",
        "model_selection": "mc_setup",
        "mc_setup_advanced": "mc_setup",
        "cascade_options": "mc_setup",
        "nuclear_stopping": "mc_setup",
        "electronic_stopping": "mc_setup",
        "histogram_settings": "mc_setup",
        "koral_solver": "koral",
        "display_settings": "mc_results",
    }

    def context_for_section(self, section_id: str) -> Optional[str]:
        """Return the page context ('mc_setup', 'mc_results', 'koral') of a section."""
        return self._SECTION_CONTEXT.get(section_id)

    def _accordion_context(self, item: "AccordionItem") -> Optional[str]:
        for sid, acc in self._accordion_by_id.items():
            if acc is item:
                return self._SECTION_CONTEXT.get(sid)
        return None

    def set_visible_context(self, context: Optional[str]) -> None:
        """Show only the accordion sections relevant to *context*.

        ``None`` shows every section. Sections with no known context are always
        shown so nothing is accidentally hidden.
        """
        for item in self._all_accordions:
            if context is None:
                item.setVisible(True)
                continue
            ctx = self._accordion_context(item)
            item.setVisible(ctx is None or ctx == context)
