from __future__ import annotations

from datetime import datetime
from typing import Callable, Optional
from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal, QObject, QThread
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QLabel,
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTableWidget, QTableWidgetItem,
    QHeaderView, QStackedWidget, QCheckBox, QPushButton, QMessageBox, QDialog,
    QListWidget, QDialogButtonBox,
    QRadioButton, QButtonGroup,
    QFileDialog,
    QSplitter, QScrollArea,
    QFrame,
    QSizePolicy,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
try:
    from matplotlib.backends.backend_qtagg import NavigationToolbar2QT  # type: ignore
except Exception:  # pragma: no cover
    from matplotlib.backends.backend_qt import NavigationToolbar2QT  # type: ignore
from matplotlib.figure import Figure

from state import AppState
from ui.widgets.toggle_switch import ToggleSwitch
from ui.widgets.periodic_table_picker import PeriodicTableButton, PeriodicTableDialog
from ui.dialogs.compound_dictionary_dialog import CompoundDictionaryDialog

try:
    from ui.logging import log as emit_log
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.logging import log as emit_log  # type: ignore

try:
    from simulators.koral_interface import (
        discover_models,
        filter_outputs,
        load_model,
        load_ui_parameters,
    )
except Exception:  # pragma: no cover
    discover_models = None  # type: ignore
    load_model = None  # type: ignore
    filter_outputs = None  # type: ignore
    load_ui_parameters = None  # type: ignore


class _KoralWorker(QObject):
    finished = pyqtSignal(list)
    error = pyqtSignal(str)

    def __init__(self, *, model_ids: list[str], request: dict, parent: Optional[QObject] = None):
        super().__init__(parent)
        self._model_ids = list(model_ids)
        self._request = dict(request)

    def run(self) -> None:
        if load_model is None:
            self.error.emit("KORAL interface not available")
            return

        requested_outputs = self._request.get("output", {}).get("requested")
        if not isinstance(requested_outputs, list) or not requested_outputs:
            self.error.emit("No output options selected")
            return

        results: list[dict] = []
        try:
            for mid in self._model_ids:
                emit_log(f"KORAL: running model '{mid}'")
                model = load_model(str(mid))
                res = model.run(self._request)
                if filter_outputs is not None:
                    res = filter_outputs(res, [str(x) for x in requested_outputs])
                results.append(res)
        except Exception as exc:
            emit_log(f"KORAL: model run failed: {exc}")
            self.error.emit(str(exc))
            return

        self.finished.emit(results)


# HintSystem import (support both workspace and packaged layouts)
try:
    from ui.widgets.hints_popup import HintSystem  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.widgets.hints_popup import HintSystem  # type: ignore


class KoralPage(QWidget):
    advanced_requested = pyqtSignal(str)

    def __init__(self, state: AppState, on_log: Optional[Callable[[str], None]] = None, parent=None):
        super().__init__(parent)
        self.state = state
        self._on_log = on_log or emit_log

        self._ion_angle: float = 0.0
        self._selected_models: list[str] = []
        self._available_models: list[str] = []
        self._primary_models: list[str] = []

        self._discovered_models: list[object] = []
        self._model_supported_outputs: dict[str, set[str]] = {}
        self._model_display_label: dict[str, str] = {}

        self.latest_log_button = None
        self._logs_dialog = None
        self._logs_list_widget = None
        self.koral_progress = None
        self.run_button = None

        self._output_option_widgets: dict[str, list[QWidget]] = {}

        # Hints (KORAL-only)
        self._hint_system: Optional[HintSystem] = None
        self._init_hints()

        # NOTE: target layers removed -> keep a flat list of element entries
        self.element_entries: list[dict] = []

        self.layer_elements: list[list[dict]] = []  # ...existing code (now unused)...
        self._updating_elements_table = False

        # Plot fullscreen toggle bookkeeping
        self._plot_fullscreen_dialog: Optional[QDialog] = None
        self._plot_widget: Optional[QWidget] = None
        self._plot_grid: Optional[QGridLayout] = None
        self._plot_placeholder: Optional[QWidget] = None

        self._table_fullscreen_dialog: Optional[QDialog] = None
        self._table_layout: Optional[QVBoxLayout] = None
        self._table_placeholder: Optional[QWidget] = None

        self._calc_thread: Optional[QThread] = None
        self._calc_worker: Optional[_KoralWorker] = None

        self._last_request: Optional[dict] = None
        self._last_results: Optional[list] = None
        self._ui_param_base_specs: dict[str, dict] = {}

        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        ion_box = self.build_ion_data()
        target_box = self.build_input_elements()

        # Top row: Ion Selection | Target Data (horizontally resizable)
        top = QSplitter(Qt.Orientation.Horizontal)
        top.addWidget(ion_box)
        top.addWidget(target_box)
        try:
            top.setStretchFactor(0, 0)
            top.setStretchFactor(1, 1)
        except Exception:
            pass

        middle = QSplitter(Qt.Orientation.Vertical)  # adjustable height
        middle.addWidget(top)
        middle.addWidget(self._build_koral_bottom_section())
        layout.addWidget(middle)

        layout.addWidget(self._build_koral_footer())

        self._refresh_element_table()

        # Models: allow running without opening the selection dialog first.
        self._refresh_available_models()
        if not self._selected_models and self._primary_models:
            # Enforce single-model selection.
            self._selected_models = [str(self._primary_models[0])]
        self._update_selected_models_label()
        self._capture_base_ui_param_specs()
        self._apply_selected_models_ui_param_specs()
        self._apply_selected_model_supported_outputs()

    def _ui_param_widgets(self) -> dict[str, object]:
        widgets: dict[str, object] = {}
        if hasattr(self, "ion_mass"):
            widgets["ion_mass_amu"] = self.ion_mass
        if hasattr(self, "energy_min"):
            widgets["energy_min_keV"] = self.energy_min
        if hasattr(self, "energy_max"):
            widgets["energy_max_keV"] = self.energy_max
        if hasattr(self, "spin_compound_corr"):
            widgets["compound_correction"] = self.spin_compound_corr
        return widgets

    def _capture_base_ui_param_specs(self) -> None:
        """Remember initial widget ranges so model constraints can be applied/reverted."""

        base: dict[str, dict] = {}
        for pid, w in self._ui_param_widgets().items():
            if isinstance(w, QSpinBox):
                base[pid] = {
                    "min": int(w.minimum()),
                    "max": int(w.maximum()),
                    "step": int(w.singleStep()),
                    "default": int(w.value()),
                }
            elif isinstance(w, QDoubleSpinBox):
                base[pid] = {
                    "min": float(w.minimum()),
                    "max": float(w.maximum()),
                    "decimals": int(w.decimals()),
                    "step": float(w.singleStep()),
                    "default": float(w.value()),
                }
        self._ui_param_base_specs = base

    def _apply_selected_models_ui_param_specs(self) -> None:
        if not self._ui_param_base_specs:
            self._capture_base_ui_param_specs()

        widgets = self._ui_param_widgets()
        if not widgets:
            return

        # Start from baseline and intersect constraints of all selected models.
        merged: dict[str, dict] = {k: dict(v) for k, v in self._ui_param_base_specs.items()}

        selected = [str(m) for m in (self._selected_models or []) if str(m).strip()]
        if selected and load_ui_parameters is not None:
            per_model: list[dict[str, dict]] = []
            for mid in selected:
                try:
                    spec = load_ui_parameters(mid)
                except Exception:
                    spec = {}
                if isinstance(spec, dict) and spec:
                    per_model.append(spec)

            for pid, base in merged.items():
                if not isinstance(base, dict):
                    continue
                try:
                    cur_min = float(base.get("min"))
                    cur_max = float(base.get("max"))
                except (TypeError, ValueError):
                    continue

                cur_decimals = base.get("decimals")
                try:
                    cur_decimals_i = int(cur_decimals) if cur_decimals is not None else None
                except (TypeError, ValueError):
                    cur_decimals_i = None

                cur_step = base.get("step")
                try:
                    cur_step_f = float(cur_step) if cur_step is not None else None
                except (TypeError, ValueError):
                    cur_step_f = None

                for model_spec in per_model:
                    p = model_spec.get(pid)
                    if not isinstance(p, dict):
                        continue
                    try:
                        cur_min = max(cur_min, float(p.get("min", cur_min)))
                        cur_max = min(cur_max, float(p.get("max", cur_max)))
                    except (TypeError, ValueError):
                        continue
                    if "decimals" in p:
                        try:
                            d = int(p.get("decimals"))
                            cur_decimals_i = d if cur_decimals_i is None else max(cur_decimals_i, d)
                        except (TypeError, ValueError):
                            pass
                    if "step" in p:
                        try:
                            s = float(p.get("step"))
                            if s > 0:
                                cur_step_f = s if cur_step_f is None else min(cur_step_f, s)
                        except (TypeError, ValueError):
                            pass

                # If the intersection is empty, fall back to the baseline.
                if cur_min <= cur_max:
                    base["min"] = cur_min
                    base["max"] = cur_max
                    if cur_decimals_i is not None:
                        base["decimals"] = cur_decimals_i
                    if cur_step_f is not None:
                        base["step"] = cur_step_f

        for pid, w in widgets.items():
            spec = merged.get(pid)
            if not isinstance(spec, dict):
                continue
            if isinstance(w, QSpinBox):
                try:
                    vmin = int(spec.get("min"))
                    vmax = int(spec.get("max"))
                except (TypeError, ValueError):
                    continue
                if vmin > vmax:
                    continue
                old = int(w.value())
                w.setRange(vmin, vmax)
                step = spec.get("step")
                try:
                    step_i = int(step) if step is not None else None
                except (TypeError, ValueError):
                    step_i = None
                if step_i is not None and step_i > 0:
                    w.setSingleStep(step_i)
                if vmin <= old <= vmax:
                    w.setValue(old)
            elif isinstance(w, QDoubleSpinBox):
                try:
                    vmin = float(spec.get("min"))
                    vmax = float(spec.get("max"))
                except (TypeError, ValueError):
                    continue
                if vmin > vmax:
                    continue
                old = float(w.value())
                w.setRange(vmin, vmax)
                decimals = spec.get("decimals")
                try:
                    decimals_i = int(decimals) if decimals is not None else None
                except (TypeError, ValueError):
                    decimals_i = None
                if decimals_i is not None and decimals_i >= 0:
                    w.setDecimals(decimals_i)
                step = spec.get("step")
                try:
                    step_f = float(step) if step is not None else None
                except (TypeError, ValueError):
                    step_f = None
                if step_f is not None and step_f > 0:
                    w.setSingleStep(step_f)
                if vmin <= old <= vmax:
                    w.setValue(old)

    def _update_selected_models_label(self) -> None:
        if not hasattr(self, "selected_models_label"):
            return
        mid = self._selected_models[0] if self._selected_models else ""
        if mid:
            label = self._model_display_label.get(mid) or str(mid)
            self.selected_models_label.setText("Selected Model:\n" + label)
        else:
            self.selected_models_label.setText("Selected Model:\nNone")

    # --- logging bridge (same pattern as MC Setup) ---
    def update_latest_log(self, entry: str) -> None:
        if self.latest_log_button:
            self.latest_log_button.setText(entry)
            self.latest_log_button.setToolTip(entry)
        if self._logs_list_widget:
            if self._logs_list_widget.count() == 1 and self._logs_list_widget.item(0).text() == "No logs available.":
                self._logs_list_widget.clear()
            self._logs_list_widget.addItem(entry)
            while self._logs_list_widget.count() > 1000:
                self._logs_list_widget.takeItem(0)
            self._logs_list_widget.scrollToBottom()

    def add_log_entry(self, message: str) -> None:
        if self._on_log:
            self._on_log(message)

    def _show_logs_dialog(self) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("Update Notifications")
        dialog.resize(800, 400)
        layout = QVBoxLayout(dialog)
        list_widget = QListWidget()

        if self.state.log_entries:
            list_widget.addItems(list(self.state.log_entries))
            list_widget.scrollToBottom()
        else:
            list_widget.addItem("No logs available.")

        layout.addWidget(list_widget)

        clear_btn = QPushButton("Clear Logs")
        clear_btn.clicked.connect(lambda: self._clear_logs(list_widget))
        layout.addWidget(clear_btn)

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

    def _clear_logs(self, list_widget: QListWidget) -> None:
        self.state.clear_logs()
        list_widget.clear()
        list_widget.addItem("No logs available.")
        if self.latest_log_button:
            self.latest_log_button.setText("No updates yet")
            self.latest_log_button.setToolTip("")

    def _build_koral_footer(self) -> QWidget:
        from PyQt6.QtWidgets import QProgressBar  # local import

        footer = QFrame(self)
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

        self.koral_progress = QProgressBar()
        self.koral_progress.setRange(0, 100)
        self.koral_progress.setValue(0)
        self.koral_progress.setFormat("Ready")

        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._handle_run_clicked)

        layout.addWidget(log_container, 2)
        layout.addWidget(self.koral_progress, 2)
        layout.addWidget(self.run_button)

        return footer

    def _handle_run_clicked(self) -> None:
        self._start_calculation_async()

    def _start_calculation_async(self) -> None:
        model_ids = list(self._selected_models)
        if not model_ids:
            QMessageBox.warning(self, "KORAL", "No calculation model selected.\nSelect at least one model in 'Model Selection'.")
            return

        request = self._collect_calculation_request()
        self._last_request = request

        # Basic input validation
        ion_symbol = str(request.get("ion", {}).get("symbol") or "").strip()
        if not ion_symbol:
            QMessageBox.warning(self, "KORAL", "No ion selected.\nPlease select an ion in 'Ion Selection'.")
            self.add_log_entry("KORAL: run blocked (no ion selected)")
            return

        elements = request.get("target", {}).get("elements")
        if not isinstance(elements, list) or not elements:
            QMessageBox.warning(self, "KORAL", "No target defined.\nPlease add at least one target element in 'Target Data'.")
            self.add_log_entry("KORAL: run blocked (no target defined)")
            return
        try:
            ratio_sum = sum(float(e.get("ratio", 0.0) or 0.0) for e in elements if isinstance(e, dict))
        except Exception:
            ratio_sum = 0.0
        if ratio_sum <= 0.0:
            QMessageBox.warning(self, "KORAL", "Target stoichiometry is zero.\nPlease set a positive ratio for at least one target element.")
            self.add_log_entry("KORAL: run blocked (target stoichiometry is zero)")
            return

        requested = request.get("output", {}).get("requested")
        if not isinstance(requested, list) or not requested:
            QMessageBox.warning(self, "KORAL", "No output option selected.\nSelect at least one output option (Projectile Range, Straggling, Stopping, ...).")
            return

        if self.run_button:
            self.run_button.setEnabled(False)
        if self.koral_progress:
            self.koral_progress.setRange(0, 0)  # busy
            self.koral_progress.setFormat("Running")

        self.add_log_entry("Starting KORAL calculation…")

        thread = QThread(self)
        worker = _KoralWorker(model_ids=model_ids, request=request)
        worker.moveToThread(thread)

        thread.started.connect(worker.run)
        worker.finished.connect(self._on_calc_finished)
        worker.error.connect(self._on_calc_error)
        worker.finished.connect(thread.quit)
        worker.error.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        worker.error.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)

        self._calc_thread = thread
        self._calc_worker = worker
        thread.start()

    def _on_calc_error(self, message: str) -> None:
        if self.koral_progress:
            self.koral_progress.setRange(0, 100)
            self.koral_progress.setValue(0)
            self.koral_progress.setFormat("Ready")
        if self.run_button:
            self.run_button.setEnabled(True)
        QMessageBox.warning(self, "KORAL", str(message))
        self.add_log_entry(f"KORAL calculation failed: {message}")

    def _on_calc_finished(self, results: list) -> None:
        if self.koral_progress:
            self.koral_progress.setRange(0, 100)
            self.koral_progress.setValue(100)
            self.koral_progress.setFormat("Complete")
        if self.run_button:
            self.run_button.setEnabled(True)

        self._last_results = results
        try:
            self._render_calculation_results(results)
        except Exception as exc:
            QMessageBox.warning(self, "KORAL", f"Unable to render results:\n{exc}")
            self.add_log_entry(f"Render failed: {exc}")
            return

        self.add_log_entry("KORAL calculation finished.")

    def _collect_calculation_request(self) -> dict:
        """Collect all KORAL parameters to send to calculation models."""

        def _srim_like_energy_grid(min_keV: float, max_keV: float) -> list[float]:
            """Generate an energy grid similar to SRIM stopping/range tables.

            Uses SRIM-like preferred steps per decade:
            1.0, 1.1, ..., 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75,
            4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0 (times 10^n).
            """

            try:
                a = float(min_keV)
                b = float(max_keV)
            except (TypeError, ValueError):
                return []
            if a < 0:
                a = 0.0
            if b < 0:
                b = 0.0
            if b < a:
                a, b = b, a
            if a == b:
                return [a]
            if a <= 0.0 or b <= 0.0:
                return []

            import math

            mantissas = [
                1.0,
                1.1,
                1.2,
                1.3,
                1.4,
                1.5,
                1.6,
                1.7,
                1.8,
                2.0,
                2.25,
                2.5,
                2.75,
                3.0,
                3.25,
                3.5,
                3.75,
                4.0,
                4.5,
                5.0,
                5.5,
                6.0,
                6.5,
                7.0,
                8.0,
                9.0,
            ]

            start_decade = int(math.floor(math.log10(a)))
            end_decade = int(math.ceil(math.log10(b)))

            out: list[float] = []
            for p in range(start_decade, end_decade + 1):
                base = 10.0**p
                for m in mantissas:
                    e = m * base
                    if a <= e <= b:
                        out.append(float(f"{e:.6g}"))

            out = sorted(set(out))
            if not out:
                out = [a, b]
            return out

        # Ion
        ion = {
            "symbol": self.ion_symbol.text() if hasattr(self, "ion_symbol") else "",
            "name": self.ion_name.text() if hasattr(self, "ion_name") else "",
            "Z": int(self.ion_z.value()) if hasattr(self, "ion_z") else 0,
            "mass_amu": float(self.ion_mass.value()) if hasattr(self, "ion_mass") else 0.0,
            "angle_deg": float(self.get_ion_angle()),
        }

        # Energy grid: SRIM-like table energies
        e0 = float(self.energy_min.value()) if hasattr(self, "energy_min") else 0.0
        e1 = float(self.energy_max.value()) if hasattr(self, "energy_max") else 0.0
        if e1 < e0:
            e0, e1 = e1, e0
        energies = _srim_like_energy_grid(e0, e1)
        if not energies and e0 > 0.0:
            energies = [e0]

        energy = {
            "min_keV": e0,
            "max_keV": e1,
            "energies_keV": energies,
        }

        # Also provide SI energies (J) for models that want strict SI input.
        eV_J = 1.602176634e-19
        energies_J = [float(e) * 1e3 * eV_J for e in energies]
        energy["energies_J"] = energies_J

        # Target composition
        elements = []
        for entry in getattr(self, "element_entries", []):
            element = entry.get("element") or {}
            if not element:
                continue

            ratio = float(entry.get("ratio", 0.0) or 0.0)
            mass_override = entry.get("mass_override")
            try:
                mass_override_amu = float(mass_override) if mass_override is not None else None
            except (TypeError, ValueError):
                mass_override_amu = None

            try:
                mass_amu = float(mass_override_amu) if mass_override_amu is not None else float(element.get("atomic_mass", 0.0) or 0.0)
            except (TypeError, ValueError):
                mass_amu = 0.0

            elements.append(
                {
                    "Z": int(element.get("number") or 0),
                    "symbol": element.get("symbol"),
                    "name": element.get("name"),
                    "mass_amu": mass_amu,
                    "ratio": ratio,
                    "damage_eV": entry.get("damage"),
                    "disp_eV": entry.get("disp"),
                    "latt_eV": entry.get("latt"),
                    "surf_eV": entry.get("surf"),
                    "mass_override_amu": mass_override_amu,
                }
            )

        # Provide SRIM-like defaults for atomic density so physics models can produce non-zero linear stopping.
        # If/when a real density input exists in the UI, this should use that.
        number_density_atoms_cm3 = 6.1238e22
        # Estimate mass density from atomic density and average atomic mass.
        total_ratio = sum(float(e.get("ratio", 0.0) or 0.0) for e in elements) if elements else 0.0
        avg_mass_amu = 0.0
        if total_ratio > 0.0:
            for e in elements:
                try:
                    avg_mass_amu += (float(e.get("ratio", 0.0) or 0.0) / total_ratio) * float(e.get("mass_amu", 0.0) or 0.0)
                except Exception:
                    continue
        amu_g = 1.66053906660e-24
        density_g_cm3 = float(number_density_atoms_cm3) * float(avg_mass_amu) * amu_g

        number_density_atoms_m3 = float(number_density_atoms_cm3) * 1.0e6
        density_kg_m3 = float(density_g_cm3) * 1000.0

        # Output selection
        requested: list[str] = []
        if hasattr(self, "chk_prange") and self.chk_prange.isChecked():
            requested.append("prange")
        if hasattr(self, "chk_long_strag") and self.chk_long_strag.isChecked():
            requested.append("long_strag")
        if hasattr(self, "chk_lat_strag") and self.chk_lat_strag.isChecked():
            requested.append("lat_strag")
        if hasattr(self, "chk_nucl_strag") and self.chk_nucl_strag.isChecked():
            requested.append("nucl_stop")
        if hasattr(self, "chk_elec_hop") and self.chk_elec_hop.isChecked():
            requested.append("elec_stop")

        units = {
            "prange": self.cmb_prange.currentText() if hasattr(self, "cmb_prange") else "",
            "long_strag": self.cmb_long_strag.currentText() if hasattr(self, "cmb_long_strag") else "",
            "lat_strag": self.cmb_lat_strag.currentText() if hasattr(self, "cmb_lat_strag") else "",
            # Separate stopping units for nuclear vs electronic.
            "nucl_stop": self.cmb_nucl_stop_unit.currentText() if hasattr(self, "cmb_nucl_stop_unit") else "",
            "elec_stop": self.cmb_elec_stop_unit.currentText() if hasattr(self, "cmb_elec_stop_unit") else "",
        }

        compound_correction = float(self.spin_compound_corr.value()) if hasattr(self, "spin_compound_corr") else 1.0

        output = {
            "requested": requested,
            "units": units,
            "compound_correction": compound_correction,
        }

        return {
            "ion": ion,
            "energy": energy,
            "target": {
                "elements": elements,
                "number_density_atoms_cm3": number_density_atoms_cm3,
                "number_density_atoms_m3": number_density_atoms_m3,
                "density_g_cm3": density_g_cm3,
                "density_kg_m3": density_kg_m3,
                "avg_mass_amu": avg_mass_amu,
                "compound_correction": compound_correction,
            },
            "output": output,
        }

    def _convert_stopping(self, value_J_per_m: float, unit: str, *, number_density_atoms_cm3: float, density_g_cm3: float) -> float:
        """Convert stopping from SI base (J/m) to the selected unit."""
        u = str(unit)
        v = float(value_J_per_m)

        eV_J = 1.602176634e-19
        # J/m -> eV/Å
        ev_per_A = v * (1.0e-10 / eV_J)

        # LSS reduced stopping unit (SRIM-style):
        # Convert linear stopping (eV/Å) -> stopping cross section (eV/(1e15 atoms/cm^2))
        # then divide by the nuclear-prefactor so that the result matches the universal
        # reduced stopping scale (dimensionless).
        if u == "L.S.S. reduced units":
            nd = float(number_density_atoms_cm3) if number_density_atoms_cm3 > 0 else 0.0
            if nd <= 0:
                return 0.0

            # linear -> cross section
            s_cs = ev_per_A / (nd * 1.0e-23)

            Z1 = 0
            M1 = 0.0
            Z2_eff = 0.0
            M2_eff = 0.0
            if isinstance(self._last_request, dict):
                ion = self._last_request.get("ion", {})
                if isinstance(ion, dict):
                    try:
                        Z1 = int(ion.get("Z", 0) or 0)
                    except (TypeError, ValueError):
                        Z1 = 0
                    try:
                        M1 = float(ion.get("mass_amu", ion.get("mass", 0.0)) or 0.0)
                    except (TypeError, ValueError):
                        M1 = 0.0

                tgt = self._last_request.get("target", {})
                if isinstance(tgt, dict):
                    els = tgt.get("elements")
                    if isinstance(els, list) and els:
                        total = 0.0
                        for e in els:
                            if not isinstance(e, dict):
                                continue
                            try:
                                total += float(e.get("ratio", 0.0) or 0.0)
                            except (TypeError, ValueError):
                                continue
                        if total > 0:
                            for e in els:
                                if not isinstance(e, dict):
                                    continue
                                try:
                                    fi = float(e.get("ratio", 0.0) or 0.0) / total
                                except (TypeError, ValueError):
                                    continue
                                if fi <= 0:
                                    continue
                                try:
                                    Z2_eff += fi * float(int(e.get("Z", 0) or 0))
                                except (TypeError, ValueError):
                                    pass
                                try:
                                    M2_eff += fi * float(e.get("mass_amu", 0.0) or 0.0)
                                except (TypeError, ValueError):
                                    pass

            if Z1 <= 0 or M1 <= 0 or Z2_eff <= 0 or M2_eff <= 0:
                return 0.0

            denom = (M1 + M2_eff) * (Z1 ** 0.23 + Z2_eff ** 0.23)
            if denom <= 0:
                return 0.0

            pref = 8.462 * (Z1 * Z2_eff * M1) / denom
            return (s_cs / pref) if pref > 0 else 0.0

        # Pure unit conversions (do not depend on target density)
        if u == "eV/Å":
            return ev_per_A
        if u == "keV/µm":
            # 1 eV/Å = 10 keV/µm
            return ev_per_A * 10.0
        if u == "MeV/mm":
            # 1 eV/Å = 10 MeV/mm
            return ev_per_A * 10.0

        # Conversions requiring atomic density / mass density
        nd = float(number_density_atoms_cm3) if number_density_atoms_cm3 > 0 else 0.0
        rho = float(density_g_cm3) if density_g_cm3 > 0 else 0.0

        if u == "eV/(10¹⁵ atoms/cm²)":
            # dE/dx[eV/Å] = S[eV/(1e15 at/cm^2)] * n[at/cm^3] * 1e-23
            return ev_per_A / (nd * 1.0e-23) if nd > 0 else 0.0

        # dE/dx (eV/Å) -> (eV/cm)
        ev_per_cm = ev_per_A * 1.0e8
        # (eV/cm) / (g/cm^3) = eV/(g/cm^2)
        ev_per_gcm2 = (ev_per_cm / rho) if rho > 0 else 0.0

        if u == "keV/(µg/cm²)":
            # eV/(g/cm2) -> keV/(ug/cm2): /1e9
            return ev_per_gcm2 / 1.0e9
        if u == "MeV/(mg/cm²)":
            # eV/(g/cm2) -> MeV/(mg/cm2): /1e9
            return ev_per_gcm2 / 1.0e9
        if u == "keV/(mg/cm²)":
            # eV/(g/cm2) -> keV/(mg/cm2): /1e6
            return ev_per_gcm2 / 1.0e6

        return v

    def _length_from_m(self, value_m: float, unit: str) -> float:
        u = str(unit)
        # AppState.unit_options default: ["Ång", "nm", "µm", "mm", "cm", "m", "km"]
        factors = {
            "Ång": 1.0e10,
            "nm": 1.0e9,
            "µm": 1.0e6,
            "mm": 1.0e3,
            "cm": 1.0e2,
            "m": 1.0,
            "km": 1.0e-3,
        }
        return float(value_m) * float(factors.get(u, 1.0))

    def _format_energy_label(self, e_keV: float, *, use_comma: bool = False) -> str:
        try:
            e = float(e_keV)
        except (TypeError, ValueError):
            e = 0.0
        if e >= 1000.0:
            s = f"{e/1000.0:.2f} MeV"
        else:
            s = f"{e:.2f} keV"
        return s.replace(".", ",") if use_comma else s

    def _format_sci(self, v: float, *, use_comma: bool = False) -> str:
        try:
            x = float(v)
        except (TypeError, ValueError):
            x = 0.0
        s = f"{x:.3E}"
        return s.replace(".", ",") if use_comma else s

    def _format_length(self, v_m: float, unit: str, *, use_comma: bool = False) -> str:
        val = self._length_from_m(v_m, unit)
        u = str(unit)
        if u == "Ång":
            s = f"{int(round(val))} A"
        else:
            s = f"{val:.3g} {u}"
        return s.replace(".", ",") if use_comma else s

    def _render_calculation_results(self, results: list) -> None:
        """Render returned model results into the plot and list table."""
        # Normalize
        norm: list[dict] = [r for r in results if isinstance(r, dict)]
        if not norm:
            raise ValueError("No results returned")

        # Plot
        if hasattr(self, "figure"):
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.set_title("KORAL results")
            ax.set_xlabel("Energy (keV)")

            # Use the same units as the table/export (UI converts from SI base units)
            range_unit = "Ång"
            if isinstance(self._last_request, dict):
                range_unit = str(self._last_request.get("output", {}).get("units", {}).get("prange") or range_unit)
            range_unit_long = "Ång"
            range_unit_lat = "Ång"
            if isinstance(self._last_request, dict):
                units = self._last_request.get("output", {}).get("units", {})
                if isinstance(units, dict):
                    range_unit_long = str(units.get("long_strag") or range_unit_long)
                    range_unit_lat = str(units.get("lat_strag") or range_unit_lat)
            stop_unit_elec = "eV/Å"
            stop_unit_nucl = "eV/Å"
            if isinstance(self._last_request, dict):
                units = self._last_request.get("output", {}).get("units", {})
                if isinstance(units, dict):
                    stop_unit_elec = str(units.get("elec_stop") or stop_unit_elec)
                    stop_unit_nucl = str(units.get("nucl_stop") or stop_unit_nucl)

            number_density_atoms_cm3 = 0.0
            density_g_cm3 = 0.0
            if isinstance(self._last_request, dict):
                try:
                    number_density_atoms_cm3 = float(self._last_request.get("target", {}).get("number_density_atoms_cm3", 0.0) or 0.0)
                except (TypeError, ValueError):
                    number_density_atoms_cm3 = 0.0
                try:
                    density_g_cm3 = float(self._last_request.get("target", {}).get("density_g_cm3", 0.0) or 0.0)
                except (TypeError, ValueError):
                    density_g_cm3 = 0.0

            plotted_any = False
            for res in norm:
                mid = str(res.get("model_id", "model"))
                energies = res.get("energies_keV")
                outputs = res.get("outputs")
                if not isinstance(energies, list) or not isinstance(outputs, dict):
                    continue
                for key, vals in outputs.items():
                    if not isinstance(vals, list):
                        continue
                    label = f"{mid}: {key}"

                    # Convert for plotting
                    plot_vals = vals
                    if key in ("elec_stop", "nucl_stop"):
                        try:
                            unit_for_stop = stop_unit_elec if key == "elec_stop" else stop_unit_nucl
                            plot_vals = [
                                self._convert_stopping(
                                    float(v),
                                    unit_for_stop,
                                    number_density_atoms_cm3=number_density_atoms_cm3,
                                    density_g_cm3=density_g_cm3,
                                )
                                for v in vals
                            ]
                        except Exception:
                            plot_vals = vals
                    elif key in ("prange", "range_csda", "long_strag", "lat_strag"):
                        try:
                            unit_for_key = range_unit
                            if key == "long_strag":
                                unit_for_key = range_unit_long
                            elif key == "lat_strag":
                                unit_for_key = range_unit_lat
                            plot_vals = [self._length_from_m(float(v), unit_for_key) for v in vals]
                        except Exception:
                            plot_vals = vals

                    ax.plot(energies, plot_vals, marker="o", label=label)
                    plotted_any = True

            if plotted_any:
                ax.legend()
            self.canvas.draw_idle()

        # Table (SRIM-like ordering)
        requested = []
        if isinstance(self._last_request, dict):
            req_out = self._last_request.get("output", {}).get("requested")
            if isinstance(req_out, list):
                requested = [str(x) for x in req_out]

        range_unit = "Ång"
        range_unit_long = "Ång"
        range_unit_lat = "Ång"
        if isinstance(self._last_request, dict):
            units = self._last_request.get("output", {}).get("units", {})
            if isinstance(units, dict):
                range_unit = str(units.get("prange") or range_unit)
                range_unit_long = str(units.get("long_strag") or range_unit_long)
                range_unit_lat = str(units.get("lat_strag") or range_unit_lat)
        stop_unit_elec = "eV/Å"
        stop_unit_nucl = "eV/Å"
        if isinstance(self._last_request, dict):
            units = self._last_request.get("output", {}).get("units", {})
            if isinstance(units, dict):
                stop_unit_elec = str(units.get("elec_stop") or stop_unit_elec)
                stop_unit_nucl = str(units.get("nucl_stop") or stop_unit_nucl)

        number_density_atoms_cm3 = 0.0
        density_g_cm3 = 0.0
        if isinstance(self._last_request, dict):
            try:
                number_density_atoms_cm3 = float(self._last_request.get("target", {}).get("number_density_atoms_cm3", 0.0) or 0.0)
            except (TypeError, ValueError):
                number_density_atoms_cm3 = 0.0
            try:
                density_g_cm3 = float(self._last_request.get("target", {}).get("density_g_cm3", 0.0) or 0.0)
            except (TypeError, ValueError):
                density_g_cm3 = 0.0

        order = ["elec_stop", "nucl_stop", "prange", "long_strag", "lat_strag"]
        keys_present: list[str] = []
        for k in order:
            if requested and k not in requested:
                continue
            # include if any result has it
            def _has_key(r: dict) -> bool:
                outs = r.get("outputs")
                return isinstance(outs, dict) and k in outs

            if any(_has_key(r) for r in norm):
                keys_present.append(k)

        multi = len(norm) > 1
        headers: list[str] = []
        if multi:
            headers.append("Model")
        headers.append("Ion\nEnergy")
        for k in keys_present:
            if k == "elec_stop":
                headers.append(f"dE/dx\nElec.\n({stop_unit_elec})")
            elif k == "nucl_stop":
                headers.append(f"dE/dx\nNuclear\n({stop_unit_nucl})")
            elif k == "prange":
                headers.append(f"Projected\nRange\n({range_unit})")
            elif k == "long_strag":
                headers.append(f"Long.\nStraggling\n({range_unit_long})")
            elif k == "lat_strag":
                headers.append(f"Lat.\nStraggling\n({range_unit_lat})")
            else:
                headers.append(k)

        rows: list[list[str]] = []
        for res in norm:
            mid = str(res.get("model_id", "model"))
            energies = res.get("energies_keV")
            outputs = res.get("outputs")
            if not isinstance(energies, list) or not isinstance(outputs, dict):
                continue
            n = len(energies)
            for i in range(n):
                base: list[str] = []
                if multi:
                    base.append(mid)
                base.append(self._format_energy_label(energies[i]))

                for k in keys_present:
                    vals = outputs.get(k)
                    if not (isinstance(vals, list) and i < len(vals)):
                        base.append("")
                        continue
                    try:
                        v = float(vals[i])
                    except (TypeError, ValueError):
                        v = 0.0

                    if k in {"elec_stop", "nucl_stop"}:
                        unit_for_stop = stop_unit_elec if k == "elec_stop" else stop_unit_nucl
                        conv = self._convert_stopping(
                            v,
                            unit_for_stop,
                            number_density_atoms_cm3=number_density_atoms_cm3,
                            density_g_cm3=density_g_cm3,
                        )
                        base.append(self._format_sci(conv))
                    elif k == "prange":
                        base.append(self._format_length(v, range_unit))
                    elif k == "long_strag":
                        base.append(self._format_length(v, range_unit_long))
                    elif k == "lat_strag":
                        base.append(self._format_length(v, range_unit_lat))
                    else:
                        base.append(str(vals[i]))
                rows.append(base)

        self.koral_result_table.clear()
        self.koral_result_table.setColumnCount(len(headers))
        self.koral_result_table.setRowCount(len(rows))
        self.koral_result_table.setHorizontalHeaderLabels(headers)
        for r, row in enumerate(rows):
            for c, txt in enumerate(row):
                self.koral_result_table.setItem(r, c, QTableWidgetItem(txt))
        hdr = self.koral_result_table.horizontalHeader()
        if hdr is not None:
            # Ensure multi-line headers are fully visible.
            try:
                max_lines = max(1, max(h.count("\n") + 1 for h in headers))
                hdr.setMinimumHeight(hdr.fontMetrics().height() * max_lines + 12)
            except Exception:
                pass
            hdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

    def _export_results_txt(self) -> None:
        if not self._last_results or not isinstance(self._last_results, list):
            QMessageBox.information(self, "KORAL", "No results to export.\nRun a calculation first.")
            return

        path, _ = QFileDialog.getSaveFileName(self, "Export SRIM-style table", "koral_results.txt", "Text Files (*.txt)")
        if not path:
            return

        text = self._build_srim_export_text(self._last_results, export_path=path)
        try:
            Path(path).write_text(text, encoding="utf-8")
        except Exception as exc:
            QMessageBox.warning(self, "KORAL", f"Failed to write file:\n{exc}")
            return

        self.add_log_entry(f"Exported results to: {path}")

    def _build_srim_export_text(self, results: list, *, export_path: str | None = None) -> str:
        norm: list[dict] = [r for r in results if isinstance(r, dict)]
        if not norm:
            return ""

        def _fmt_E4(x: float) -> str:
            try:
                s = f"{float(x):.4E}"
            except (TypeError, ValueError):
                s = f"{0.0:.4E}"
            return s.replace(".", ",")

        def _fmt_2(x: float) -> str:
            try:
                s = f"{float(x):.2f}"
            except (TypeError, ValueError):
                s = f"{0.0:.2f}"
            return s.replace(".", ",")

        def _srim_month_name(m: int) -> str:
            # Match the example's German month spellings.
            return {
                1: "Jänner",
                2: "Februar",
                3: "März",
                4: "April",
                5: "Mai",
                6: "Juni",
                7: "Juli",
                8: "August",
                9: "September",
                10: "Oktober",
                11: "November",
                12: "Dezember",
            }.get(int(m), "")

        # Units
        range_unit = "Ång"
        range_unit_long = "Ång"
        range_unit_lat = "Ång"
        # TXT export requirement: always use SRIM default stopping unit
        stop_unit_elec = "eV/Å"
        stop_unit_nucl = "eV/Å"
        compound_corr = 1.0
        requested: list[str] = []
        number_density_atoms_cm3 = 0.0
        density_g_cm3 = 0.0
        ion_name = ""
        ion_Z = 0
        ion_mass_amu = 0.0
        elements_for_header: list[dict] = []
        if isinstance(self._last_request, dict):
            out = self._last_request.get("output", {})
            requested_raw = out.get("requested")
            if isinstance(requested_raw, list):
                requested = [str(x) for x in requested_raw]
            units = out.get("units", {})
            if isinstance(units, dict):
                range_unit = str(units.get("prange") or range_unit)
                range_unit_long = str(units.get("long_strag") or range_unit_long)
                range_unit_lat = str(units.get("lat_strag") or range_unit_lat)
                # ignore stopping units for TXT export
            try:
                compound_corr = float(out.get("compound_correction", 1.0))
            except (TypeError, ValueError):
                compound_corr = 1.0

            try:
                number_density_atoms_cm3 = float(self._last_request.get("target", {}).get("number_density_atoms_cm3", 0.0) or 0.0)
            except (TypeError, ValueError):
                number_density_atoms_cm3 = 0.0
            try:
                density_g_cm3 = float(self._last_request.get("target", {}).get("density_g_cm3", 0.0) or 0.0)
            except (TypeError, ValueError):
                density_g_cm3 = 0.0

            ion = self._last_request.get("ion", {})
            if isinstance(ion, dict):
                ion_name = str(ion.get("name") or ion.get("symbol") or "").strip()
                try:
                    ion_Z = int(ion.get("Z", 0) or 0)
                except (TypeError, ValueError):
                    ion_Z = 0
                try:
                    ion_mass_amu = float(ion.get("mass_amu", 0.0) or 0.0)
                except (TypeError, ValueError):
                    ion_mass_amu = 0.0

            tgt = self._last_request.get("target", {})
            if isinstance(tgt, dict):
                els = tgt.get("elements")
                if isinstance(els, list):
                    elements_for_header = [e for e in els if isinstance(e, dict)]

        bragg_corr_pct = (compound_corr - 1.0) * 100.0
        # Column order exactly like SRIM table snippet
        order = ["elec_stop", "nucl_stop", "prange", "long_strag", "lat_strag"]
        keys_present: list[str] = []
        for k in order:
            if requested and k not in requested:
                continue
            def _has_key(r: dict) -> bool:
                outs = r.get("outputs")
                return isinstance(outs, dict) and k in outs

            if any(_has_key(r) for r in norm):
                keys_present.append(k)

        def line_for_point(e_keV: float, outs: dict, i: int) -> str:
            energy = self._format_energy_label(e_keV, use_comma=True).rjust(12)
            vals: list[str] = []
            for k in keys_present:
                series = outs.get(k)
                if not isinstance(series, list) or i >= len(series):
                    vals.append("".rjust(10))
                    continue
                try:
                    v = float(series[i])
                except (TypeError, ValueError):
                    v = 0.0
                if k in {"elec_stop", "nucl_stop"}:
                    # Always export stopping in eV/Å (SRIM default)
                    unit_for_stop = "eV/Å"
                    conv = self._convert_stopping(
                        v,
                        unit_for_stop,
                        number_density_atoms_cm3=number_density_atoms_cm3,
                        density_g_cm3=density_g_cm3,
                    )
                    vals.append(self._format_sci(conv, use_comma=True).rjust(10))
                else:
                    unit_for_k = range_unit
                    if k == "long_strag":
                        unit_for_k = range_unit_long
                    elif k == "lat_strag":
                        unit_for_k = range_unit_lat
                    vals.append(self._format_length(v, unit_for_k, use_comma=True).rjust(10))
            return f"  {energy}  " + " ".join(vals)

        now_dt = datetime.now()
        calc_date = f"{_srim_month_name(now_dt.month)} {now_dt.day:02d}, {now_dt.year}"
        out_lines: list[str] = []

        # SRIM-like header block
        out_lines.append(" ==================================================================")
        out_lines.append("              OpenSrim version ---> OpenSrim-2025.00")
        out_lines.append(f"              Calc. date   ---> {calc_date} ")
        out_lines.append(" ==================================================================")
        out_lines.append("")

        disk_name = "koral_results.txt"
        if export_path:
            try:
                disk_name = Path(str(export_path)).name
            except Exception:
                disk_name = str(export_path)
        out_lines.append(f" Disk File Name = {disk_name}")
        out_lines.append("")

        if ion_name or ion_Z or ion_mass_amu:
            out_lines.append(f" Ion = {ion_name} [{ion_Z}] , Mass = {ion_mass_amu:g} amu")
            out_lines.append("")

        out_lines.append(
            f" Target Density =  {_fmt_E4(density_g_cm3)} g/cm3 = {_fmt_E4(number_density_atoms_cm3)} atoms/cm3"
        )
        out_lines.append(" ======= Target  Composition ========")
        out_lines.append("    Atom   Atom   Atomic    Mass     ")
        out_lines.append("    Name   Numb   Percent   Percent  ")
        out_lines.append("    ----   ----   -------   -------  ")

        total_ratio = 0.0
        for e in elements_for_header:
            try:
                total_ratio += float(e.get("ratio", 0.0) or 0.0)
            except (TypeError, ValueError):
                continue

        atomic_fracs: list[tuple[str, int, float, float]] = []
        if total_ratio > 0.0:
            for e in elements_for_header:
                sym = str(e.get("symbol") or "").strip()
                try:
                    z = int(e.get("Z", 0) or 0)
                except (TypeError, ValueError):
                    z = 0
                try:
                    r = float(e.get("ratio", 0.0) or 0.0)
                except (TypeError, ValueError):
                    r = 0.0
                try:
                    mass = float(e.get("mass_amu", 0.0) or 0.0)
                except (TypeError, ValueError):
                    mass = 0.0
                if r <= 0:
                    continue
                atomic_fracs.append((sym, z, r / total_ratio, mass))

        mass_total = 0.0
        for _sym, _z, af, mass in atomic_fracs:
            mass_total += af * mass

        for sym, z, af, mass in atomic_fracs:
            ap = af * 100.0
            mp = (af * mass / mass_total * 100.0) if mass_total > 0 else 0.0
            out_lines.append(f"{sym:>6} {z:>6} {_fmt_2(ap):>9} {_fmt_2(mp):>9}")

        out_lines.append(" ====================================")
        out_lines.append(f" Bragg Correction = {bragg_corr_pct:.2f}%".replace(".", ","))
        out_lines.append(" Stopping Units =   eV / Angstrom ")
        out_lines.append(" See bottom of Table for other Stopping units ")
        out_lines.append("")
        out_lines.append("        Ion        dE/dx      dE/dx     Projected  Longitudinal   Lateral")
        out_lines.append("       Energy      Elec.      Nuclear     Range     Straggling   Straggling")
        out_lines.append("  --------------  ---------- ---------- ----------  ----------  ----------")

        multi = len(norm) > 1
        for res in norm:
            mid = str(res.get("model_id", "model"))
            energies = res.get("energies_keV")
            outputs = res.get("outputs")
            if not isinstance(energies, list) or not isinstance(outputs, dict):
                continue
            if multi:
                out_lines.append("")
                out_lines.append(f" Model = {mid}")
            for i, e in enumerate(energies):
                try:
                    ee = float(e)
                except (TypeError, ValueError):
                    ee = 0.0
                out_lines.append(line_for_point(ee, outputs, i))

        out_lines.append("-----------------------------------------------------------")
        out_lines.append(" Multiply Stopping by        for Stopping Units")
        out_lines.append(" -------------------        ------------------")
        # Multipliers relative to 1 eV/Å
        conv_table: list[tuple[float, str]] = []
        conv_table.append((1.0, "eV / Angstrom"))
        conv_table.append((10.0, "keV / micron"))
        conv_table.append((10.0, "MeV / mm"))

        # Density-dependent conversions (only meaningful if density is known)
        if density_g_cm3 > 0:
            # See _convert_stopping: keV/(ug/cm2) and MeV/(mg/cm2) are /1e9 after eV/(g/cm2)
            factor = (1.0e8 / density_g_cm3) / 1.0e9
            conv_table.append((factor, "keV / (ug/cm2)"))
            conv_table.append((factor, "MeV / (mg/cm2)"))
            conv_table.append(((1.0e8 / density_g_cm3) / 1.0e6, "keV / (mg/cm2)"))
        if number_density_atoms_cm3 > 0:
            cs_factor = 1.0 / (number_density_atoms_cm3 * 1.0e-23)
            conv_table.append((cs_factor, "eV / (1E15 atoms/cm2)"))

            # LSS reduced units (SRIM-style): divide the stopping cross section by
            # the ZBL nuclear prefactor so the result is on the universal Sn(ε) scale.
            Z2_eff = 0.0
            M2_eff = 0.0
            for _sym, _z, af, mass in atomic_fracs:
                Z2_eff += af * float(_z)
                M2_eff += af * float(mass)
            if ion_Z > 0 and ion_mass_amu > 0 and Z2_eff > 0 and M2_eff > 0:
                denom = (ion_mass_amu + M2_eff) * (ion_Z ** 0.23 + Z2_eff ** 0.23)
                if denom > 0:
                    pref = 8.462 * (ion_Z * Z2_eff * ion_mass_amu) / denom
                    if pref > 0:
                        conv_table.append((cs_factor / pref, "L.S.S. reduced units"))

        for mult, name in conv_table:
            out_lines.append(f" {mult: .4E}                {name}".replace(".", ","))

        out_lines.append("")
        out_lines.append("  Program name: OpenSrim")
        return "\n".join(out_lines) + "\n"


    def _init_hints(self) -> None:
        """Initialize hint system and lock it to this page."""
        try:
            hints_path = Path(__file__).resolve().parents[1] / "widgets" / "hints.json"  # app/ui/widgets/hints.json
            hs = HintSystem(repo_path=hints_path, parent=self)
            hs.set_current_page("KORAL")
            self._hint_system = hs
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
            return self._hint_system.make_hint_button(page_id="KORAL", hint_id=hint_id, parent=parent)
        except Exception:
            btn = QPushButton("?", parent)
            btn.setEnabled(False)
            btn.setFixedSize(22, 22)
            btn.setToolTip(f"Hint '{hint_id}' not available")
            return btn

    def _build_koral_bottom_section(self) -> QSplitter:
        splitter = QSplitter(Qt.Orientation.Horizontal)  # Use QSplitter for horizontal resizing
        splitter.addWidget(self._build_model_selection())  # Add "Model Selection" as its own column
        splitter.addWidget(self._build_koral_left_options())  # Add "Ausgabe Optionen"
        splitter.addWidget(self._build_koral_plot_list_section())  # Add "Range / Straggling"
        return splitter

    def _build_model_selection(self) -> QGroupBox:
        box = QGroupBox("")
        layout = QVBoxLayout(box)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.setContentsMargins(9, 0, 9, 9)

        # Header (title + hint button)
        layout.addWidget(self._groupbox_header("Model Selection", hint_id="model", parent=box))

        self.model_button = QPushButton("Select Model")
        self.model_button.clicked.connect(self._open_model_selection_dialog)
        self.selected_models_label = QLabel("Selected Model:\nNone")
        self.selected_models_label.setWordWrap(True)

        layout.addWidget(self.model_button)
        layout.addWidget(self.selected_models_label)
        layout.addStretch(1)

        return box

    def _refresh_available_models(self) -> None:
        self._available_models = []
        self._primary_models = []
        self._discovered_models = []
        self._model_supported_outputs = {}
        self._model_display_label = {}
        if discover_models is None:
            return
        try:
            discovered = discover_models()
        except Exception:
            discovered = []
        self._discovered_models = list(discovered)

        avail: list[str] = []
        primary: list[str] = []
        for m in discovered:
            mid = getattr(m, "id", None)
            if not isinstance(mid, str) or not mid.strip():
                continue
            avail.append(mid)
            if bool(getattr(m, "is_primary", False)):
                primary.append(mid)

            pkg = getattr(m, "package_display_name", None) or getattr(m, "package_id", None) or ""
            mdl = getattr(m, "model_display_name", None) or getattr(m, "model_id", None) or mid
            if isinstance(pkg, str) and pkg.strip():
                self._model_display_label[mid] = f"{pkg.strip()}: {str(mdl).strip()}"
            else:
                self._model_display_label[mid] = str(mdl).strip()

            outs = getattr(m, "supported_outputs", None)
            if isinstance(outs, (list, tuple)):
                self._model_supported_outputs[mid] = {str(x) for x in outs if str(x).strip()}

        self._available_models = avail
        self._primary_models = primary

    def _open_model_selection_dialog(self):
        self._refresh_available_models()

        dialog = QDialog(self)
        dialog.setWindowTitle("Select Model")
        dialog.setModal(True)
        dialog.resize(380, 380)

        layout = QVBoxLayout(dialog)
        scroll_area = QScrollArea(dialog)
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        radio_group = QButtonGroup(dialog)
        radio_group.setExclusive(True)
        self._model_radio_buttons: list[QRadioButton] = []

        if not self._available_models:
            scroll_layout.addWidget(QLabel("No KORAL models found."))
        else:
            # Group models by package display name.
            packages: dict[str, list[object]] = {}
            for m in (self._discovered_models or []):
                pkg = getattr(m, "package_display_name", None) or getattr(m, "package_id", None) or "KORAL"
                packages.setdefault(str(pkg), []).append(m)

            # Determine which model should be selected when opening.
            current = self._selected_models[0] if self._selected_models else ""
            if current and current in self._available_models:
                preselect_id = current
            elif self._primary_models:
                preselect_id = str(self._primary_models[0])
            else:
                preselect_id = str(self._available_models[0])

            for pkg_name in sorted(packages.keys()):
                pkg_label = QLabel(str(pkg_name))
                pkg_label.setStyleSheet("font-weight: 600;")
                scroll_layout.addWidget(pkg_label)

                for m in packages[pkg_name]:
                    mid = getattr(m, "id", None)
                    if not isinstance(mid, str) or not mid.strip():
                        continue
                    model_display = getattr(m, "model_display_name", None) or getattr(m, "model_id", None) or mid
                    rb = QRadioButton(str(model_display))
                    rb.setProperty("koral_model_id", mid)
                    rb.setChecked(mid == preselect_id)
                    radio_group.addButton(rb)
                    self._model_radio_buttons.append(rb)
                    scroll_layout.addWidget(rb)

                scroll_layout.addSpacing(6)

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
        selected: list[str] = []
        for rb in getattr(self, "_model_radio_buttons", []) or []:
            if rb.isChecked():
                mid = rb.property("koral_model_id")
                if isinstance(mid, str) and mid.strip():
                    selected = [mid]
                break
        self._selected_models = selected
        self._update_selected_models_label()
        self._apply_selected_models_ui_param_specs()
        self._apply_selected_model_supported_outputs()
        dialog.accept()

    def _apply_selected_model_supported_outputs(self) -> None:
        """Hide output options that the selected model does not support."""

        if not self._output_option_widgets:
            return

        mid = self._selected_models[0] if self._selected_models else ""
        supported = self._model_supported_outputs.get(str(mid)) if mid else None
        if not supported:
            # If unknown, do not hide anything.
            for widgets in self._output_option_widgets.values():
                for w in widgets:
                    w.setVisible(True)
            return

        for out_id, widgets in self._output_option_widgets.items():
            visible = out_id in supported
            for w in widgets:
                w.setVisible(bool(visible))
            # If hidden, also uncheck (avoid "invisible selected" outputs).
            if not visible:
                for w in widgets:
                    if isinstance(w, QCheckBox):
                        w.setChecked(False)

    # -------- configuration persistence ----------
    def collect_config(self) -> dict:
        ion = {
            "symbol": self.ion_symbol.text() if hasattr(self, "ion_symbol") else "",
            "name": self.ion_name.text() if hasattr(self, "ion_name") else "",
            "number": int(self.ion_z.value()) if hasattr(self, "ion_z") else 0,
            "mass": float(self.ion_mass.value()) if hasattr(self, "ion_mass") else 0.0,
            "energy_min": float(self.energy_min.value()) if hasattr(self, "energy_min") else 0.0,
            "energy_max": float(self.energy_max.value()) if hasattr(self, "energy_max") else 0.0,
            "angle": float(self.get_ion_angle()),
        }

        # Selected model (single-choice)
        models = list(self._selected_models)

        elements = []
        for entry in getattr(self, "element_entries", []):
            element = entry.get("element") or {}
            z = element.get("number")
            if z is None:
                continue
            elements.append(
                {
                    "Z": int(z),
                    "ratio": entry.get("ratio", 0.0),
                    "damage": entry.get("damage"),
                    "disp": entry.get("disp"),
                    "latt": entry.get("latt"),
                    "surf": entry.get("surf"),
                    "mass_override": entry.get("mass_override"),
                }
            )

        output = {}
        for name in (
            "chk_prange",
            "chk_long_strag",
            "chk_lat_strag",
            "chk_nucl_strag",
            "chk_elec_hop",
        ):
            if hasattr(self, name):
                output[name] = bool(getattr(self, name).isChecked())
        for name in ("cmb_prange", "cmb_long_strag", "cmb_lat_strag", "cmb_elect_unit"):
            if hasattr(self, name):
                output[name] = getattr(self, name).currentText()
        if hasattr(self, "cmb_elec_stop_unit"):
            output["cmb_elec_stop_unit"] = self.cmb_elec_stop_unit.currentText()
        if hasattr(self, "cmb_nucl_stop_unit"):
            output["cmb_nucl_stop_unit"] = self.cmb_nucl_stop_unit.currentText()
        if hasattr(self, "spin_compound_corr"):
            output["compound_corr"] = float(self.spin_compound_corr.value())
        if hasattr(self, "sw_koral_mode"):
            output["sw_koral_mode"] = bool(self.sw_koral_mode.isChecked())
        if hasattr(self, "all_none_chk"):
            output["all_none_chk"] = bool(self.all_none_chk.isChecked())

        return {
            "ion": ion,
            "models": models,
            "elements": elements,
            "output": output,
        }

    def apply_config(self, payload: dict) -> None:
        if not isinstance(payload, dict):
            return

        ion = payload.get("ion") or {}
        if isinstance(ion, dict):
            if hasattr(self, "ion_symbol"):
                self.ion_symbol.setText(str(ion.get("symbol", "")))
            if hasattr(self, "ion_name"):
                self.ion_name.setText(str(ion.get("name", "")))
            if hasattr(self, "ion_z"):
                try:
                    self.ion_z.setValue(int(ion.get("number", self.ion_z.value())))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "ion_mass"):
                try:
                    self.ion_mass.setValue(float(ion.get("mass", self.ion_mass.value())))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "energy_min"):
                try:
                    self.energy_min.setValue(float(ion.get("energy_min", self.energy_min.value())))
                except (TypeError, ValueError):
                    pass
            if hasattr(self, "energy_max"):
                try:
                    self.energy_max.setValue(float(ion.get("energy_max", self.energy_max.value())))
                except (TypeError, ValueError):
                    pass
            if "angle" in ion:
                try:
                    self.set_ion_angle(float(ion.get("angle", 0.0)))
                except (TypeError, ValueError):
                    pass

        models = payload.get("models") or []
        if isinstance(models, list):
            # Enforce single-model selection and migrate old IDs like "ZBL" -> "koral:ZBL".
            self._refresh_available_models()
            cand = [str(m) for m in models if isinstance(m, (str, int, float)) and str(m).strip()]
            chosen = cand[0] if cand else ""
            if chosen and ":" not in chosen:
                namespaced = f"koral:{chosen}"
                if namespaced in self._available_models:
                    chosen = namespaced

            if chosen and self._available_models and chosen not in self._available_models:
                chosen = ""

            if not chosen and self._primary_models:
                chosen = str(self._primary_models[0])

            self._selected_models = [chosen] if chosen else []
            self._update_selected_models_label()
            self._apply_selected_models_ui_param_specs()
            self._apply_selected_model_supported_outputs()

        # elements table
        if hasattr(self, "element_entries"):
            self.element_entries = []
            for e in payload.get("elements") or []:
                if not isinstance(e, dict):
                    continue
                z = e.get("Z")
                if z is None:
                    continue
                element = self.state.elements_by_number.get(int(z))
                if not element:
                    continue
                overrides = {k: e.get(k) for k in ("damage", "disp", "latt", "surf")}
                added = self._add_element_to_table(element, e.get("ratio", 0.0), overrides=overrides, refresh=False)
                if isinstance(added, dict):
                    added["mass_override"] = e.get("mass_override")
            self._refresh_element_table()

        output = payload.get("output") or {}
        if isinstance(output, dict):
            for name in (
                "chk_prange",
                "chk_long_strag",
                "chk_lat_strag",
                "chk_nucl_strag",
                "chk_elec_hop",
            ):
                if hasattr(self, name) and name in output:
                    getattr(self, name).setChecked(bool(output.get(name)))

            for name in ("cmb_prange", "cmb_long_strag", "cmb_lat_strag", "cmb_elect_unit"):
                if hasattr(self, name) and name in output:
                    val = output.get(name)
                    if isinstance(val, str):
                        cb = getattr(self, name)
                        idx = cb.findText(val)
                        if idx >= 0:
                            cb.setCurrentIndex(idx)

            # Backward compatibility:
            # - older configs used a single stopping-unit dropdown (cmb_stop_power_unit)
            #   -> apply to both new dropdowns.
            legacy_stop = output.get("cmb_stop_power_unit")
            if isinstance(legacy_stop, str) and legacy_stop:
                for name in ("cmb_elec_stop_unit", "cmb_nucl_stop_unit"):
                    if hasattr(self, name):
                        cb = getattr(self, name)
                        idx = cb.findText(legacy_stop)
                        if idx >= 0:
                            cb.setCurrentIndex(idx)

            for name in ("cmb_elec_stop_unit", "cmb_nucl_stop_unit"):
                if hasattr(self, name) and name in output:
                    val = output.get(name)
                    if isinstance(val, str):
                        cb = getattr(self, name)
                        idx = cb.findText(val)
                        if idx >= 0:
                            cb.setCurrentIndex(idx)

            if hasattr(self, "spin_compound_corr") and "compound_corr" in output:
                try:
                    self.spin_compound_corr.setValue(float(output.get("compound_corr", self.spin_compound_corr.value())))
                except (TypeError, ValueError):
                    pass

            if hasattr(self, "sw_koral_mode") and "sw_koral_mode" in output:
                self.sw_koral_mode.setChecked(bool(output.get("sw_koral_mode")))
            if hasattr(self, "all_none_chk") and "all_none_chk" in output:
                self.all_none_chk.setChecked(bool(output.get("all_none_chk")))

    def _groupbox_header(self, title: str, hint_id: Optional[str] = None, parent: Optional[QWidget] = None) -> QWidget:
        """Header row: title label + optional hint button directly to the right."""
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

    # --- identical to MCSetup: Ion Selection ---
    def build_ion_data(self) -> QGroupBox:
        box = QGroupBox("")
        layout = QVBoxLayout(box)
        # Keep header close to the top edge (like "Target Data").
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        # Header (title + hint button)
        layout.addWidget(self._groupbox_header("Ion Selection", hint_id="ion", parent=box))

        grid = QGridLayout()
        grid.setHorizontalSpacing(6)   # was 10
        grid.setVerticalSpacing(6)     # was 10
        # Optional: slightly tighter margins inside the box
        # grid.setContentsMargins(0, 0, 0, 0)

        # --- Periodic table button (ion) ---
        self.pick_btn = PeriodicTableButton(
            "Click to Select Element",
            compact=True,
            show_hover_info=True,
            bordered=True,
            update_button_text=True,
        )
        self.pick_btn.element_selected.connect(self.on_element_selected)

        pick_row = QWidget()
        pick_row_l = QHBoxLayout(pick_row)
        pick_row_l.setContentsMargins(0, 0, 0, 0)
        pick_row_l.setSpacing(6)
        pick_row_l.addWidget(self.pick_btn)
        pick_row_l.addStretch(1)

        # Keep the selector in the input row, but top-align it so it doesn't look vertically centered.
        grid.addWidget(pick_row, 1, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        # Keep these widgets for data/config, but do not show them (PeriodicTableButton already displays them).
        self.ion_symbol = QLineEdit()
        self.ion_symbol.setReadOnly(True)
        self.ion_symbol.setVisible(False)

        self.ion_name = QLineEdit()
        self.ion_name.setReadOnly(True)
        self.ion_name.setVisible(False)

        self.ion_z = QSpinBox()
        self.ion_z.setRange(1, 120)
        self.ion_z.setReadOnly(True)
        self.ion_z.setVisible(False)

        # Mass (amu)
        grid.addWidget(QLabel("Mass (amu)"), 0, 1, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.ion_mass = QDoubleSpinBox()
        self.ion_mass.setRange(0.01, 1000.0)
        self.ion_mass.setDecimals(3)
        self.ion_mass.setReadOnly(False)
        self.ion_mass.setMaximumWidth(120)
        self.ion_mass.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_mass, 1, 1, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        # Energy min/max (keV)
        grid.addWidget(QLabel("Energy min (keV)"), 0, 2, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.energy_min = QDoubleSpinBox()
        self.energy_min.setRange(0.0, 1e6)
        self.energy_min.setDecimals(2)
        self.energy_min.setValue(10.0)
        self.energy_min.setMaximumWidth(130)
        self.energy_min.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.energy_min, 1, 2, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        grid.addWidget(QLabel("Energy max (keV)"), 0, 3, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.energy_max = QDoubleSpinBox()
        self.energy_max.setRange(0.0, 1e6)
        self.energy_max.setDecimals(2)
        self.energy_max.setValue(10000.0)
        self.energy_max.setMaximumWidth(130)
        self.energy_max.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.energy_max, 1, 3, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        # Keep any extra space below the inputs.
        grid.setRowStretch(2, 1)

        # Force extra horizontal space to the far right (keeps fields packed on the left)
        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 0)
        grid.setColumnStretch(2, 0)
        grid.setColumnStretch(3, 0)
        grid.setColumnStretch(4, 1)  # trailing stretch column

        layout.addLayout(grid)
        return box

    def set_ion_angle(self, angle: float) -> None:
        try:
            self._ion_angle = float(angle)
        except (TypeError, ValueError):
            return

    def get_ion_angle(self) -> float:
        return float(self._ion_angle)

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
        box = QGroupBox("")
        v = QVBoxLayout(box)
        # Keep header close to the top edge (like Ion Selection).
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(6)

        # Header row: title + hint + (Add Element | Compound Dictionary)
        header = QWidget(box)
        header_l = QHBoxLayout(header)
        header_l.setContentsMargins(0, 0, 0, 0)
        header_l.setSpacing(6)
        title_lbl = QLabel("Target Data", header)
        title_lbl.setStyleSheet("font-weight: 600;")
        header_l.addWidget(title_lbl)
        header_l.addWidget(self._hint_btn("elements", parent=header))
        header_l.addStretch(1)

        dict_btn = QPushButton("Compound Dictionary")
        dict_btn.clicked.connect(self._open_compound_dictionary)
        header_l.addWidget(dict_btn)

        v.addWidget(header)

        # Table
        self.elem_table = QTableWidget(0, 11)
        self.elem_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.elem_table.setHorizontalHeaderLabels([
            "", "Symbol", "Name", "Atomic No.", "Weight (amu)",
            "Atom Stoich", "Atom Stoich %", "Damage (eV)", "Disp (eV)", "Latt (eV)", "Surf (eV)"
        ])
        hdr = self.elem_table.horizontalHeader()
        if hdr is not None:
            hdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            # keep delete column compact
            try:
                hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
            except Exception:
                pass
            self.elem_table.setColumnWidth(0, 28)
            f = hdr.font()
            f.setBold(True)
            hdr.setFont(f)
        vh = self.elem_table.verticalHeader()
        if vh is not None:
            vh.setVisible(False)
        self.elem_table.setAlternatingRowColors(True)
        self.elem_table.setColumnHidden(8, True)
        self.elem_table.setColumnHidden(9, True)
        self.elem_table.setColumnHidden(10, True)

        v.addWidget(self.elem_table)

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

        # Ensure the action row is visible even when empty.
        self._refresh_element_table()

        return box

    def _delete_element_row(self, row: int) -> None:
        if not (0 <= row < len(self.element_entries)):
            return
        self.element_entries.pop(row)
        self._refresh_element_table()

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
        if ratio_value < 0:
            ratio_value = 0.0

        # If the element is already present, increment its stoichiometry instead of adding a new row.
        try:
            element_number = int(element.get("number"))
        except Exception:
            element_number = None

        if element_number is not None:
            for entry in self.element_entries:
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
                if refresh:
                    self._refresh_element_table()
                return entry

        entry = {
            "element": element,
            "ratio": ratio_value,
            "damage": energy_defaults["damage"],
            "disp": energy_defaults["disp"],
            "latt": energy_defaults["latt"],
            "surf": energy_defaults["surf"],
        }
        self.element_entries.append(entry)
        if refresh:
            self._refresh_element_table()
        return entry

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
        # Last row is reserved for the in-table "Add element" action.
        self.elem_table.setRowCount(len(entries) + 1)

        def ro_item(text: str):
            it = QTableWidgetItem(text)
            it.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled)
            return it

        for row, entry in enumerate(entries):
            element = entry["element"]
            # Delete button (first column)
            btn = QPushButton("-")
            btn.setFixedSize(20, 20)
            btn.setStyleSheet("padding: 0px;")
            btn.setToolTip("Delete element")
            btn.clicked.connect(lambda _=False, r=row: self._delete_element_row(r))
            # Center the button inside the cell (independent of table styles/alignment).
            cell = QWidget(self.elem_table)
            cell_l = QHBoxLayout(cell)
            cell_l.setContentsMargins(0, 0, 0, 0)
            cell_l.setSpacing(0)
            cell_l.addStretch(1)
            cell_l.addWidget(btn)
            cell_l.addStretch(1)
            self.elem_table.setCellWidget(row, 0, cell)

            # Symbol, Name, Atomic No.
            self.elem_table.setItem(row, 1, ro_item(element["symbol"]))
            self.elem_table.setItem(row, 2, ro_item(element["name"]))
            self.elem_table.setItem(row, 3, ro_item(str(element["number"])))

            # Mass (amu) -- now editable
            mass_val = entry.get("mass_override")
            if mass_val is not None:
                mass_text = f"{float(mass_val):.3f}"
            else:
                mass_raw = element.get("atomic_mass")
                try:
                    mass_text = f"{float(mass_raw):.3f}"
                except (TypeError, ValueError):
                    mass_text = str(mass_raw)
            mass_item = QTableWidgetItem(mass_text)
            mass_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
            self.elem_table.setItem(row, 4, mass_item)

            # Ratio (editable)
            ratio_item = QTableWidgetItem(f"{entry['ratio']:.4f}")
            ratio_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
            self.elem_table.setItem(row, 5, ratio_item)

            percent = (entry["ratio"] / total_ratio * 100.0) if total_ratio else 0.0
            self.elem_table.setItem(row, 6, ro_item(f"{percent:.2f}"))

            for offset, key in enumerate(("damage", "disp", "latt", "surf"), start=7):
                self.elem_table.setItem(row, offset, ro_item(str(entry[key])))

        # Action row: span across all columns and place the add-element button.
        action_row = self.elem_table.rowCount() - 1
        for c in range(self.elem_table.columnCount()):
            self.elem_table.setCellWidget(action_row, c, None)
            self.elem_table.setItem(action_row, c, None)
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
        action_l.setContentsMargins(0, 0, 0, 0)
        action_l.setSpacing(6)
        if hasattr(self, "elem_pick_btn") and self.elem_pick_btn is not None:
            action_l.addWidget(self.elem_pick_btn)
        action_l.addStretch(1)
        self.elem_table.setCellWidget(action_row, 0, action_cell)

        self._updating_elements_table = False

    def _handle_element_item_changed(self, item):
        if self._updating_elements_table:
            return
        row = item.row()
        if not (0 <= row < len(self.element_entries)):
            return

        # Ratio column
        if item.column() == 5:
            try:
                self.element_entries[row]["ratio"] = max(float(item.text()), 0.0)
            except ValueError:
                self.element_entries[row]["ratio"] = 0.0
            self._refresh_element_table()
        # Mass column (editable)
        elif item.column() == 4:
            try:
                val = float(item.text())
                self.element_entries[row]["mass_override"] = val
            except ValueError:
                self.element_entries[row]["mass_override"] = None
            self._refresh_element_table()

    def _handle_element_cell_double_clicked(self, row, column):
        # Replace element via periodic table: double-click Symbol column
        if column != 1:
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
        box = QGroupBox("")
        v = QVBoxLayout(box)
        v.setSpacing(6)

        # Map of output-id -> widgets to show/hide together.
        # Keys match the engine output ids used in requests/results.
        self._output_option_widgets = {}

        # Header (title + hint button)
        v.addWidget(self._groupbox_header("Output Options", hint_id="output", parent=box))

        row_prange = QHBoxLayout()
        self.chk_prange = QCheckBox("Projectile Range")
        self.cmb_prange = QComboBox()
        self.cmb_prange.clear()
        self.cmb_prange.addItems(self.state.unit_options)
        row_prange.addWidget(self.chk_prange)
        row_prange.addStretch(1)
        row_prange.addWidget(self.cmb_prange)
        v.addLayout(row_prange)
        self._output_option_widgets["prange"] = [self.chk_prange, self.cmb_prange]

        row_long = QHBoxLayout()
        self.chk_long_strag = QCheckBox("Long. Straggling")
        self.cmb_long_strag = QComboBox()
        self.cmb_long_strag.clear()
        self.cmb_long_strag.addItems(self.state.unit_options)
        row_long.addWidget(self.chk_long_strag)
        row_long.addStretch(1)
        row_long.addWidget(self.cmb_long_strag)
        v.addLayout(row_long)
        self._output_option_widgets["long_strag"] = [self.chk_long_strag, self.cmb_long_strag]

        row_lat = QHBoxLayout()
        self.chk_lat_strag = QCheckBox("Lat. Straggling")
        self.cmb_lat_strag = QComboBox()
        self.cmb_lat_strag.clear()
        self.cmb_lat_strag.addItems(self.state.unit_options)
        row_lat.addWidget(self.chk_lat_strag)
        row_lat.addStretch(1)
        row_lat.addWidget(self.cmb_lat_strag)
        v.addLayout(row_lat)
        self._output_option_widgets["lat_strag"] = [self.chk_lat_strag, self.cmb_lat_strag]

        row_nucl = QHBoxLayout()
        self.chk_nucl_strag = QCheckBox("Nuclear Stopping")
        row_nucl.addWidget(self.chk_nucl_strag)
        row_nucl.addStretch(1)
        self.cmb_nucl_stop_unit = QComboBox()
        self.cmb_nucl_stop_unit.addItems(
            [
                "L.S.S. reduced units",
                "eV/Å",
                "keV/µm",
                "MeV/mm",
                "keV/(µg/cm²)",
                "MeV/(mg/cm²)",
                "keV/(mg/cm²)",
                "eV/(10¹⁵ atoms/cm²)",
            ]
        )
        row_nucl.addWidget(self.cmb_nucl_stop_unit)
        v.addLayout(row_nucl)
        self._output_option_widgets["nucl_stop"] = [self.chk_nucl_strag, self.cmb_nucl_stop_unit]

        row_elect = QHBoxLayout()
        self.chk_elec_hop = QCheckBox("Electron Stopping")
        row_elect.addWidget(self.chk_elec_hop)
        row_elect.addStretch(1)
        self.cmb_elec_stop_unit = QComboBox()
        self.cmb_elec_stop_unit.addItems(
            [
                "L.S.S. reduced units",
                "eV/Å",
                "keV/µm",
                "MeV/mm",
                "keV/(µg/cm²)",
                "MeV/(mg/cm²)",
                "keV/(mg/cm²)",
                "eV/(10¹⁵ atoms/cm²)",
            ]
        )
        row_elect.addWidget(self.cmb_elec_stop_unit)
        v.addLayout(row_elect)
        self._output_option_widgets["elec_stop"] = [self.chk_elec_hop, self.cmb_elec_stop_unit]

        # Compound correction input
        row_corr = QHBoxLayout()
        row_corr.addWidget(QLabel("Compound correction"))
        self.spin_compound_corr = QDoubleSpinBox()
        self.spin_compound_corr.setRange(0.0, 10.0)
        self.spin_compound_corr.setDecimals(4)
        self.spin_compound_corr.setSingleStep(0.01)
        self.spin_compound_corr.setValue(1.0)
        row_corr.addStretch(1)
        row_corr.addWidget(self.spin_compound_corr)
        v.addLayout(row_corr)

        row_switch = QHBoxLayout()
        row_switch.addWidget(QLabel("Plot"))
        self.sw_koral_mode = ToggleSwitch()
        self.sw_koral_mode.setChecked(False)  # False = Plot, True = List
        row_switch.addWidget(self.sw_koral_mode)
        row_switch.addWidget(QLabel("List"))
        row_switch.addStretch(1)
        v.addLayout(row_switch)

        self.sw_koral_mode.toggled.connect(self._update_koral_plot_view)

        # --- "All" checkbox ganz unten, leicht links eingerückt ---
        v.addStretch(1)
        all_row = QHBoxLayout()
        all_row.addSpacing(10)
        self.all_none_chk = QCheckBox("All")
        self.all_none_chk.setTristate(False)
        self.all_none_chk.setFixedWidth(50)
        self.all_none_chk.setChecked(False)
        self.all_none_chk.stateChanged.connect(self._toggle_all_options)
        all_row.addWidget(self.all_none_chk)
        all_row.addStretch(1)
        v.addLayout(all_row)

        return box

    def _toggle_all_options(self, state):
        # Toggle all checkboxes between checked and unchecked based on all_none_chk
        new_state = self.all_none_chk.isChecked()
        for checkbox in [self.chk_prange, self.chk_long_strag, self.chk_lat_strag, self.chk_nucl_strag, self.chk_elec_hop]:
            if checkbox.isVisible():
                checkbox.setChecked(new_state)

    def _build_koral_plot_list_section(self) -> QGroupBox:
        box = QGroupBox("")
        v = QVBoxLayout(box)

        # Intentionally no header here to maximize vertical space for the plot.

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

        export_row = QHBoxLayout()
        export_row.addStretch(1)
        btn_export = QPushButton("Export TXT")
        btn_export.clicked.connect(self._export_results_txt)
        export_row.addWidget(btn_export)
        lv.addLayout(export_row)

        self.koral_result_table = QTableWidget(3, 3)
        self.koral_result_table.setHorizontalHeaderLabels(["Energy\n(keV)", "Range\n(µm)", "Straggling\n(µm)"])
        hh = self.koral_result_table.horizontalHeader()
        if hh is not None:
            try:
                hh.setMinimumHeight(hh.fontMetrics().height() * 2 + 12)
            except Exception:
                pass
            hh.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        vh = self.koral_result_table.verticalHeader()
        if vh is not None:
            vh.setVisible(False)
        self.koral_result_table.setAlternatingRowColors(True)

        for row, (e, r, s) in enumerate([("10", "0.5", "0.1"), ("20", "1.0", "0.2"), ("30", "1.5", "0.3")]):
            self.koral_result_table.setItem(row, 0, QTableWidgetItem(e))
            self.koral_result_table.setItem(row, 1, QTableWidgetItem(r))
            self.koral_result_table.setItem(row, 2, QTableWidgetItem(s))

        lv.addWidget(self.koral_result_table)
        self._table_layout = lv

        # Placeholder used while table is re-parented to fullscreen dialog
        self._table_placeholder = QWidget()
        self._table_placeholder.setMinimumSize(1, 1)

        # Double-click table to fullscreen
        self.koral_result_table.cellDoubleClicked.connect(lambda *_: self._toggle_table_fullscreen())
        self.koral_plot_list_stack.addWidget(list_page)

        self.koral_plot_list_stack.setCurrentIndex(0)
        return box

    def _toggle_table_fullscreen(self) -> None:
        if not self.koral_result_table or not self._table_layout:
            return
        if self._table_fullscreen_dialog and self._table_fullscreen_dialog.isVisible():
            self._table_fullscreen_dialog.close()
            return

        # Keep list layout occupied while table is detached
        try:
            self._table_layout.removeWidget(self.koral_result_table)
            if self._table_placeholder:
                self._table_layout.addWidget(self._table_placeholder)
        except Exception:
            pass

        dlg = QDialog(self)
        dlg.setWindowTitle("Results")
        dlg.setModal(False)
        lay = QVBoxLayout(dlg)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        self.koral_result_table.setParent(dlg)
        lay.addWidget(self.koral_result_table)

        def _restore() -> None:
            if not self._table_layout:
                return
            try:
                if self._table_placeholder:
                    self._table_layout.removeWidget(self._table_placeholder)
                    self._table_placeholder.setParent(None)
                self.koral_result_table.setParent(self._table_layout.parentWidget())
                self._table_layout.addWidget(self.koral_result_table)
            except Exception:
                pass

        dlg.finished.connect(_restore)
        self._table_fullscreen_dialog = dlg
        dlg.showFullScreen()

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
