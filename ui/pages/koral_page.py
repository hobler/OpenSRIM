from __future__ import annotations

from datetime import datetime
from typing import Callable, Optional
from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal, QObject, QThread, QEvent, QTimer
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


def _load_density_table() -> dict[str, float]:
    """Load atomic densities from material_densities.csv.

    Returns a dict mapping upper-case element symbol to atomic number density
    in atoms/cm³.  The CSV column ``dens`` is stored as atoms/Å³ × 100
    (SRIM convention), so we convert: dens_atoms_cm3 = dens_csv × 1e22.
    """
    csv_path = Path(__file__).parents[2] / "data" / "material_densities" / "material_densities.csv"
    table: dict[str, float] = {}
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
                    dens_csv = float(parts[3])
                except ValueError:
                    continue
                # dens column is in units of 0.01 atoms/Å³  →  multiply by 1e22 to get atoms/cm³
                table[symbol] = dens_csv * 1.0e22
    except OSError:
        pass
    return table


_DENSITY_TABLE: dict[str, float] = _load_density_table()


class _ProportionalColumnFilter(QObject):
    """Event filter that rescales table columns proportionally on resize."""

    def __init__(self, table: QTableWidget, fixed_cols: dict[int, int],
                 stretch_weights: dict[int, int]) -> None:
        super().__init__(table)
        self._table = table
        self._fixed = fixed_cols          # col_idx -> fixed px width
        self._weights = stretch_weights   # col_idx -> relative weight
        self._total_weight = sum(stretch_weights.values())
        table.installEventFilter(self)

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:  # type: ignore[override]
        if obj is self._table and event.type() == QEvent.Type.Resize:
            self._resize_columns()
        return False

    def _resize_columns(self) -> None:
        vp = self._table.viewport()
        if vp is None:
            return
        available = vp.width() - sum(self._fixed.values())
        if self._total_weight <= 0 or available <= 0:
            return
        for col, weight in self._weights.items():
            if not self._table.isColumnHidden(col):
                self._table.setColumnWidth(
                    col, max(1, int(available * weight / self._total_weight))
                )


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

        # Match MC Setup container appearance.
        self.setStyleSheet(
            "QGroupBox { border: 2px solid palette(shadow); border-radius: 4px;"
            " margin-top: 6px; padding-top: 6px; }"
        )

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
        self._koral_solver_settings: dict = {}
        self._ui_param_base_specs: dict[str, dict] = {}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 4)
        layout.setSpacing(0)

        ion_box = self.build_ion_data()
        target_box = self.build_input_elements()
        model_box = self._build_model_selection()

        # Left column: Ion Selection / Target Data / Model Selection
        left_col = QSplitter(Qt.Orientation.Vertical)
        left_col.addWidget(ion_box)
        left_col.addWidget(target_box)
        left_col.addWidget(model_box)
        try:
            left_col.setStretchFactor(0, 0)
            left_col.setStretchFactor(1, 1)
            left_col.setStretchFactor(2, 0)
        except Exception:
            pass

        # Right column: Output Options / Plot+List
        right_col = QSplitter(Qt.Orientation.Vertical)
        right_col.addWidget(self._build_koral_left_options())
        right_col.addWidget(self._build_koral_plot_list_section())
        try:
            right_col.setStretchFactor(0, 0)
            right_col.setStretchFactor(1, 1)
        except Exception:
            pass

        # Main two-column splitter – equal initial widths
        main_splitter = QSplitter(Qt.Orientation.Horizontal)
        main_splitter.addWidget(left_col)
        main_splitter.addWidget(right_col)
        try:
            main_splitter.setStretchFactor(0, 1)
            main_splitter.setStretchFactor(1, 1)
        except Exception:
            pass
        layout.addWidget(main_splitter, 1)
        # Equalize after the layout has resolved real pixel sizes
        def _equalize_splitter() -> None:
            total = main_splitter.width()
            if total > 0:
                half = total // 2
                main_splitter.setSizes([half, total - half])
        QTimer.singleShot(0, _equalize_splitter)

        def _prioritize_koral_plot_height() -> None:
            total = right_col.height()
            if total > 0:
                top = min(max(170, total // 4), 240)
                right_col.setSizes([top, max(1, total - top)])
        QTimer.singleShot(0, _prioritize_koral_plot_height)

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

        # Clear Logs and Close share one row (same layout as MC Setup).
        btn_row = QHBoxLayout()
        clear_btn = QPushButton("Clear Logs")
        clear_btn.clicked.connect(lambda: self._clear_logs(list_widget))
        btn_row.addWidget(clear_btn)
        btn_row.addStretch(1)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button_box.rejected.connect(dialog.reject)
        btn_row.addWidget(button_box)
        layout.addLayout(btn_row)

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
        footer = QFrame(self)
        footer.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
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

        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._handle_run_clicked)

        load_btn = QPushButton("Load")
        load_btn.setToolTip("Load KORAL configuration and results")
        load_btn.clicked.connect(self._handle_load_koral)

        save_btn = QPushButton("Save")
        save_btn.setToolTip("Save KORAL configuration")
        save_btn.clicked.connect(self._handle_save_koral)

        layout.addWidget(log_container, 2)
        layout.addWidget(load_btn)
        layout.addWidget(save_btn)
        layout.addWidget(self.run_button)

        return footer

    def _handle_run_clicked(self) -> None:
        self._start_calculation_async()

    # -------- Save / Load KORAL (.toml) ----------

    @staticmethod
    def _to_toml_str(data: dict) -> str:
        """Minimal TOML serializer for the KORAL save format."""
        import math

        def _val(v):
            # Unwrap numpy scalars to plain Python types
            if hasattr(v, "item"):
                v = v.item()
            if v is None:
                return '""'
            if isinstance(v, bool):
                return "true" if v else "false"
            if isinstance(v, int):
                return str(v)
            if isinstance(v, float):
                if math.isnan(v) or math.isinf(v):
                    return str(v)
                return repr(v)
            if isinstance(v, str):
                return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
            if isinstance(v, list):
                return "[" + ", ".join(_val(x) for x in v) + "]"
            if isinstance(v, dict):
                return "{" + ", ".join(f"{k} = {_val(vv)}" for k, vv in v.items() if vv is not None) + "}"
            return '"' + str(v) + '"'

        def _section(prefix: str, d: dict, lines: list):
            inline, nested = {}, {}
            for k, v in d.items():
                if isinstance(v, dict):
                    nested[k] = v
                else:
                    inline[k] = v
            for k, v in inline.items():
                lines.append(f"{k} = {_val(v)}")
            for k, v in nested.items():
                key = f"{prefix}.{k}" if prefix else k
                lines.append(f"\n[{key}]")
                _section(key, v, lines)

        lines: list[str] = []
        _section("", data, lines)
        return "\n".join(lines)

    def _handle_save_koral(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save KORAL", self.state.get_dialog_start_directory(), "KORAL Files (*.koral);;All Files (*)"
        )
        if not path:
            return
        if not path.endswith(".koral"):
            path += ".koral"

        config = self.collect_config()

        # Serialize results: list of dicts with energies_keV + outputs
        results_list = []
        if self._last_results:
            for r in self._last_results:
                if not isinstance(r, dict):
                    continue
                entry: dict = {"model_id": str(r.get("model_id", ""))}
                energies = r.get("energies_keV")
                if isinstance(energies, list):
                    entry["energies_keV"] = energies
                outputs = r.get("outputs")
                if isinstance(outputs, dict):
                    for key, vals in outputs.items():
                        entry[f"out_{key}"] = vals if isinstance(vals, list) else []
                results_list.append(entry)

        request_section: dict = {}
        if isinstance(self._last_request, dict):
            out = self._last_request.get("output") or {}
            request_section["requested"] = out.get("requested") or []
            request_section["units"] = out.get("units") or {}

        payload = {
            "config": config,
            "request": request_section,
        }
        toml_str = self._to_toml_str(payload)

        # Append results as TOML array of tables (manual, since writer is not available)
        if results_list:
            toml_str += "\n"
            for entry in results_list:
                toml_str += "\n[[results]]\n"
                toml_str += self._to_toml_str(entry).lstrip("\n")

        try:
            Path(path).write_text(toml_str, encoding="utf-8")
            self.state.remember_dialog_path(path)
            self.add_log_entry(f"Saved to {Path(path).name}")
        except OSError as exc:
            QMessageBox.warning(self, "KORAL", f"Unable to save file:\n{exc}")

    def _handle_load_koral(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load KORAL", self.state.get_dialog_start_directory(), "KORAL Files (*.koral);;All Files (*)"
        )
        if not path:
            return
        try:
            try:
                import tomllib
            except ImportError:
                import tomli as tomllib  # type: ignore
            with open(path, "rb") as fh:
                payload = tomllib.load(fh)
        except Exception as exc:
            QMessageBox.warning(self, "KORAL", f"Unable to read file:\n{exc}")
            return

        config = payload.get("config")
        if isinstance(config, dict):
            self.apply_config(config)
        self.state.remember_dialog_path(path)

        # Restore _last_request units/requested for correct plot axis labels
        request_section = payload.get("request") or {}
        if request_section:
            self._last_request = {
                "output": {
                    "requested": request_section.get("requested") or [],
                    "units": request_section.get("units") or {},
                }
            }

        # Restore results and re-render
        raw_results = payload.get("results")
        if isinstance(raw_results, list) and raw_results:
            restored = []
            for entry in raw_results:
                if not isinstance(entry, dict):
                    continue
                model_id = entry.get("model_id", "model")
                energies_keV = entry.get("energies_keV") or []
                outputs: dict = {}
                for k, v in entry.items():
                    if k.startswith("out_"):
                        outputs[k[4:]] = v
                restored.append({
                    "model_id": model_id,
                    "energies_keV": energies_keV,
                    "outputs": outputs,
                })
            self._last_results = restored
            try:
                self._render_calculation_results(restored)
            except Exception as exc:
                QMessageBox.warning(self, "KORAL", f"Unable to render results:\n{exc}")

        self.add_log_entry(f"Loaded from {Path(path).name}")

    def _start_calculation_async(self) -> None:
        model_ids = list(self._selected_models)
        if not model_ids:
            QMessageBox.warning(self, "KORAL", "No calculation model selected.\nSelect at least one model in 'Model Selection'.")
            return

        request = self._collect_calculation_request()
        self._last_request = request

        nd = float(request.get("target", {}).get("number_density_atoms_cm3", 0.0) or 0.0)
        if nd <= 0:
            QMessageBox.warning(
                self, "KORAL",
                "Target atomic density is zero.\n\n"
                "Add at least one target element — the density will be auto-filled from the database.\n"
                "You can also enter it manually in the 'Target Density' field."
            )
            return

        # Print request in __main__.py style for debugging
        ion = request.get("ion", {})
        tgt = request.get("target", {})
        elements = tgt.get("elements", [])
        z_targets = [e.get("Z") for e in elements]
        m_targets = [e.get("mass_amu") for e in elements]
        d_target = tgt.get("number_density_atoms_cm3", 0.0)
        compound_corr = tgt.get("compound_correction", 1.0)
        s_e_f = [compound_corr] * len(elements)
        print("--- KORAL run settings ---")
        print(f"  method='ZBL',")
        print(f"  z_ion={ion.get('Z')},")
        print(f"  m_ion={ion.get('mass_amu')},")
        print(f"  z_target={z_targets},")
        print(f"  m_target={m_targets},")
        print(f"  d_target={d_target / 1e24:.5f},  # atoms/Å³  ({d_target:.4e} atoms/cm³)")
        print(f"  s_e_f={s_e_f},")
        energy = request.get("energy", {})
        energy_min = energy.get("min_keV", "?")
        energy_max = energy.get("max_keV", "?")
        print(f"  start_energy={energy_min},  # keV")
        print(f"  stop_energy={energy_max},  # keV")
        print(f"  nr_values=100)")
        print("--------------------------")

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
        if self.run_button:
            self.run_button.setEnabled(True)
        QMessageBox.warning(self, "KORAL", str(message))
        self.add_log_entry(f"KORAL calculation failed: {message}")

    def _on_calc_finished(self, results: list) -> None:
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

        # Energy grid: always compute from 1 eV (0.001 keV),
        # but keep the UI minimum as display/export lower bound.
        e0 = float(self.energy_min.value()) if hasattr(self, "energy_min") else 0.0
        e1 = float(self.energy_max.value()) if hasattr(self, "energy_max") else 0.0
        if e1 < e0:
            e0, e1 = e1, e0
        calc_min_keV = 1.0e-3
        calc_max_keV = max(e1, calc_min_keV)
        energies = _srim_like_energy_grid(calc_min_keV, calc_max_keV)
        if not energies:
            energies = [calc_min_keV, calc_max_keV] if calc_max_keV > calc_min_keV else [calc_min_keV]
        if energies and energies[0] > calc_min_keV:
            energies.insert(0, calc_min_keV)

        energy = {
            "min_keV": e0,
            "max_keV": e1,
            "display_min_keV": e0,
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

        # Read atomic number density from the UI spinner (weighted from DB, user-overridable).
        number_density_atoms_cm3 = (
            self._get_target_density()
            if hasattr(self, "spin_target_density")
            else 0.0
        )
        # Estimate mass density from atomic density and stoichiometry-weighted atomic mass.
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
                "gas": bool(self.chk_gas.isChecked()) if hasattr(self, "chk_gas") else False,
            },
            "output": output,
            "solver": self._koral_solver_settings,
        }

    # The target density value is stored canonically as an atomic number
    # density in atoms/cm³; the UI may display it in a different unit.
    _AMU_G = 1.66053906660e-24

    def _target_avg_mass_amu(self) -> float:
        """Stoichiometry-weighted average atomic mass of the target (amu)."""
        entries = getattr(self, "element_entries", [])
        total_ratio = sum(float(e.get("ratio", 0.0) or 0.0) for e in entries)
        if total_ratio <= 0.0:
            return 0.0
        avg = 0.0
        for e in entries:
            ratio = float(e.get("ratio", 0.0) or 0.0)
            override = e.get("mass_override")
            mass = 0.0
            if override is not None:
                try:
                    mass = float(override)
                except (ValueError, TypeError):
                    mass = 0.0
            if mass <= 0.0:
                mass = float((e.get("element") or {}).get("atomic_mass", 0.0) or 0.0)
            avg += (ratio / total_ratio) * mass
        return avg

    def _density_to_canonical(self, value: float, unit: str) -> float:
        """Convert a displayed density value (in `unit`) to atoms/cm³."""
        if unit == "atoms/cm³":
            return value
        if unit == "atoms/m³":
            return value * 1.0e-6
        denom = self._target_avg_mass_amu() * self._AMU_G
        if denom <= 0.0:
            return 0.0
        if unit == "g/cm³":
            return value / denom
        if unit == "kg/m³":
            return value / 1000.0 / denom
        return value

    def _density_from_canonical(self, value: float, unit: str) -> float:
        """Convert an atoms/cm³ value to the displayed `unit`."""
        if unit == "atoms/cm³":
            return value
        if unit == "atoms/m³":
            return value * 1.0e6
        g_cm3 = value * self._target_avg_mass_amu() * self._AMU_G
        if unit == "g/cm³":
            return g_cm3
        if unit == "kg/m³":
            return g_cm3 * 1000.0
        return value

    def _on_density_unit_changed(self) -> None:
        """Re-display the target density value in the newly selected unit."""
        if not hasattr(self, "cmb_density_unit"):
            return
        new_unit = self.cmb_density_unit.currentText()
        old_unit = getattr(self, "_density_unit", new_unit)
        if new_unit == old_unit:
            return
        try:
            val = float(self.spin_target_density.text())
        except (ValueError, TypeError):
            val = 0.0
        canonical = self._density_to_canonical(val, old_unit)
        self._density_unit = new_unit
        self.spin_target_density.setText(
            f"{self._density_from_canonical(canonical, new_unit):.1e}"
        )

    def _get_target_density(self) -> float:
        """Return the target density in atoms/cm³, regardless of display unit."""
        if not hasattr(self, "spin_target_density"):
            return 0.0
        try:
            val = float(self.spin_target_density.text())
        except (ValueError, TypeError):
            return 0.0
        return self._density_to_canonical(val, getattr(self, "_density_unit", "atoms/cm³"))

    def _set_target_density(self, value: float) -> None:
        """Set the target density (given in atoms/cm³) in the current display unit."""
        if not hasattr(self, "spin_target_density"):
            return
        disp = self._density_from_canonical(value, getattr(self, "_density_unit", "atoms/cm³"))
        self.spin_target_density.setText(f"{disp:.1e}")

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
            s = f"{e/1000.0:.1f} MeV"
        else:
            s = f"{e:.1f} keV"
        return s.replace(".", ",") if use_comma else s

    def _format_sci(self, v: float, *, use_comma: bool = False) -> str:
        try:
            x = float(v)
        except (TypeError, ValueError):
            x = 0.0
        s = f"{x:.1E}"
        return s.replace(".", ",") if use_comma else s

    def _format_length(self, v_m: float, unit: str, *, use_comma: bool = False) -> str:
        val = self._length_from_m(v_m, unit)
        u = str(unit)
        if u == "Ång":
            s = f"{int(round(val))} A"
        else:
            s = f"{val:.1f} {u}"
        return s.replace(".", ",") if use_comma else s

    def _filter_results_for_display(self, results: list[dict]) -> list[dict]:
        """Filter energies below the UI minimum for table/plot/export display."""
        display_min_keV = 0.0
        if isinstance(self._last_request, dict):
            try:
                display_min_keV = float(self._last_request.get("energy", {}).get("display_min_keV", 0.0) or 0.0)
            except (TypeError, ValueError):
                display_min_keV = 0.0

        if display_min_keV <= 0.0:
            return [r for r in results if isinstance(r, dict)]

        filtered: list[dict] = []
        for res in results:
            if not isinstance(res, dict):
                continue
            energies = res.get("energies_keV")
            outputs = res.get("outputs")
            if not isinstance(energies, list) or not isinstance(outputs, dict):
                continue

            keep_idx: list[int] = []
            for i, e in enumerate(energies):
                try:
                    ee = float(e)
                except (TypeError, ValueError):
                    continue
                if ee >= display_min_keV:
                    keep_idx.append(i)
            if not keep_idx:
                continue

            out_new: dict[str, list[float]] = {}
            for key, vals in outputs.items():
                if not isinstance(vals, list):
                    continue
                out_new[str(key)] = [vals[i] for i in keep_idx if i < len(vals)]

            filtered.append(
                {
                    "model_id": res.get("model_id", "model"),
                    "energies_keV": [energies[i] for i in keep_idx],
                    "outputs": out_new,
                }
            )
        return filtered

    def _render_calculation_results(self, results: list) -> None:
        """Render returned model results into the plot and list table."""
        # Normalize
        norm: list[dict] = self._filter_results_for_display([r for r in results if isinstance(r, dict)])
        if not norm:
            raise ValueError("No results returned")

        # Plot
        if hasattr(self, "figure"):
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.set_title("KORAL results")
            ax.set_xlabel("Energy (eV)")
            ax.set_xscale("log")
            ax.set_yscale("log")

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
            plotted_stop_units: list[str] = []
            plotted_range_units: list[str] = []
            for res in norm:
                mid = str(res.get("model_id", "model"))
                energies = [e * 1e3 for e in (res.get("energies_keV") or [])]  # keV -> eV
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
                            if unit_for_stop not in plotted_stop_units:
                                plotted_stop_units.append(unit_for_stop)
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
                            if unit_for_key not in plotted_range_units:
                                plotted_range_units.append(unit_for_key)
                        except Exception:
                            plot_vals = vals

                    ax.plot(energies, plot_vals, linewidth=1.0, label=label)
                    plotted_any = True

            # Build ylabel from plotted quantities
            ylabel_parts: list[str] = []
            if plotted_range_units:
                units_str = " / ".join(plotted_range_units)
                ylabel_parts.append(f"Range / Straggling ({units_str})")
            if plotted_stop_units:
                units_str = " / ".join(plotted_stop_units)
                ylabel_parts.append(f"Stopping power ({units_str})")
            ax.set_ylabel(" | ".join(ylabel_parts) if ylabel_parts else "Value")

            if plotted_any:
                ax.legend()
            self.figure.tight_layout()
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

        default_path = str(Path(self.state.get_dialog_start_directory()) / "koral_results.txt")
        path, _ = QFileDialog.getSaveFileName(self, "Export SRIM-style table", default_path, "Text Files (*.txt)")
        if not path:
            return

        text = self._build_srim_export_text(self._last_results, export_path=path)
        try:
            Path(path).write_text(text, encoding="utf-8")
        except Exception as exc:
            QMessageBox.warning(self, "KORAL", f"Failed to write file:\n{exc}")
            return

        self.state.remember_dialog_path(path)
        self.add_log_entry(f"Exported results to: {path}")

    def _export_results_csv(self) -> None:
        if not self._last_results or not isinstance(self._last_results, list):
            QMessageBox.information(self, "KORAL", "No results to export.\nRun a calculation first.")
            return

        default_path = str(Path(self.state.get_dialog_start_directory()) / "koral_results.csv")
        path, _ = QFileDialog.getSaveFileName(self, "Export KORAL CSV", default_path, "CSV Files (*.csv)")
        if not path:
            return

        text = self._build_csv_export_text(self._last_results)
        try:
            Path(path).write_text(text, encoding="utf-8")
        except Exception as exc:
            QMessageBox.warning(self, "KORAL", f"Failed to write file:\n{exc}")
            return

        self.state.remember_dialog_path(path)
        self.add_log_entry(f"Exported CSV results to: {path}")

    def _build_csv_export_text(self, results: list) -> str:
        norm: list[dict] = self._filter_results_for_display([r for r in results if isinstance(r, dict)])
        if not norm:
            return ""

        requested = []
        units: dict = {}
        if isinstance(self._last_request, dict):
            req_out = self._last_request.get("output", {}).get("requested")
            if isinstance(req_out, list):
                requested = [str(x) for x in req_out]
            raw_units = self._last_request.get("output", {}).get("units")
            if isinstance(raw_units, dict):
                units = raw_units

        range_unit = str(units.get("prange") or "Ång")
        range_unit_long = str(units.get("long_strag") or "Ång")
        range_unit_lat = str(units.get("lat_strag") or "Ång")
        stop_unit_elec = str(units.get("elec_stop") or "eV/Å")
        stop_unit_nucl = str(units.get("nucl_stop") or "eV/Å")

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

        keys_order = ["elec_stop", "nucl_stop", "prange", "long_strag", "lat_strag", "nucl_strag"]
        keys_present: list[str] = []
        for k in keys_order:
            if requested and k not in requested:
                continue
            if any(isinstance(r.get("outputs"), dict) and k in r.get("outputs", {}) for r in norm):
                keys_present.append(k)

        name_map = {
            "elec_stop": "dE_dx_elec",
            "nucl_stop": "dE_dx_nuclear",
            "prange": "projected_range",
            "long_strag": "longitudinal_straggling",
            "lat_strag": "lateral_straggling",
            "nucl_strag": "nuclear_straggling_qn",
        }
        unit_map = {
            "elec_stop": stop_unit_elec,
            "nucl_stop": stop_unit_nucl,
            "prange": range_unit,
            "long_strag": range_unit_long,
            "lat_strag": range_unit_lat,
            "nucl_strag": "eV²/Å",
        }

        lines: list[str] = []
        lines.append("# OpenSRIM KORAL export")
        if isinstance(self._last_request, dict):
            tgt = self._last_request.get("target", {})
            if isinstance(tgt, dict):
                elems = tgt.get("elements")
                if isinstance(elems, list):
                    lines.append("# damage_energies_eV: symbol,damage,disp,latt,surf")
                    for e in elems:
                        if not isinstance(e, dict):
                            continue
                        sym = str(e.get("symbol") or "")
                        lines.append(
                            f"# {sym};{e.get('damage_eV', '')};{e.get('disp_eV', '')};{e.get('latt_eV', '')};{e.get('surf_eV', '')}"
                        )

        lines.append("# units: energy_keV=keV; " + "; ".join(f"{name_map[k]}={unit_map[k]}" for k in keys_present))

        header = ["model_id", "energy_keV"] + [name_map[k] for k in keys_present]
        lines.append(";".join(header))

        for res in norm:
            model_id = str(res.get("model_id", "model"))
            energies = res.get("energies_keV") or []
            outputs = res.get("outputs") or {}
            if not isinstance(energies, list) or not isinstance(outputs, dict):
                continue

            for i, e_keV in enumerate(energies):
                row = [model_id, f"{float(e_keV):.6g}"]
                for k in keys_present:
                    vals = outputs.get(k)
                    if not isinstance(vals, list) or i >= len(vals):
                        row.append("")
                        continue
                    try:
                        v = float(vals[i])
                    except (TypeError, ValueError):
                        row.append("")
                        continue

                    if k in {"elec_stop", "nucl_stop"}:
                        unit_for_stop = stop_unit_elec if k == "elec_stop" else stop_unit_nucl
                        conv = self._convert_stopping(
                            v,
                            unit_for_stop,
                            number_density_atoms_cm3=number_density_atoms_cm3,
                            density_g_cm3=density_g_cm3,
                        )
                        row.append(f"{conv:.6g}")
                    elif k == "prange":
                        row.append(f"{self._length_from_m(v, range_unit):.6g}")
                    elif k == "long_strag":
                        row.append(f"{self._length_from_m(v, range_unit_long):.6g}")
                    elif k == "lat_strag":
                        row.append(f"{self._length_from_m(v, range_unit_lat):.6g}")
                    else:
                        row.append(f"{v:.6g}")
                lines.append(";".join(row))

        return "\n".join(lines) + "\n"

    def _build_srim_export_text(self, results: list, *, export_path: str | None = None) -> str:
        norm: list[dict] = self._filter_results_for_display([r for r in results if isinstance(r, dict)])
        if not norm:
            return ""

        def _fmt_E4(x: float) -> str:
            try:
                s = f"{float(x):.4E}"
            except (TypeError, ValueError):
                s = f"{0.0:.4E}"
            return s

        def _fmt_2(x: float) -> str:
            try:
                s = f"{float(x):.2f}"
            except (TypeError, ValueError):
                s = f"{0.0:.2f}"
            return s

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
        order = ["elec_stop", "nucl_stop", "prange", "long_strag", "lat_strag", "nucl_strag"]
        keys_present: list[str] = []
        for k in order:
            if requested and k not in requested:
                continue
            def _has_key(r: dict) -> bool:
                outs = r.get("outputs")
                return isinstance(outs, dict) and k in outs

            if any(_has_key(r) for r in norm):
                keys_present.append(k)

        def _format_length_no_unit(v_m: float, unit: str) -> str:
            """Format length value without the unit suffix."""
            val = self._length_from_m(v_m, unit)
            u = str(unit)
            if u == "Ång":
                s = f"{int(round(val))}"
            else:
                s = f"{val:.3g}"
            return s

        def _format_energy_no_unit(e_keV: float) -> str:
            """Format energy in keV without the unit suffix."""
            try:
                e = float(e_keV)
            except (TypeError, ValueError):
                e = 0.0
            return f"{e:.2f}"

        unit_by_key = {
            "elec_stop": stop_unit_elec,
            "nucl_stop": stop_unit_nucl,
            "prange": range_unit,
            "long_strag": range_unit_long,
            "lat_strag": range_unit_lat,
            "nucl_strag": "eV²/Å",
        }
        title_by_key = {
            "elec_stop": ("dE/dx", "Elec."),
            "nucl_stop": ("dE/dx", "Nuclear"),
            "prange": ("Projected", "Range"),
            "long_strag": ("Longitudinal", "Straggling"),
            "lat_strag": ("Lateral", "Straggling"),
            "nucl_strag": ("Nuclear", "Straggling (Qn)"),
        }

        columns: list[tuple[str, int]] = [("energy", 12)]
        for k in keys_present:
            columns.append((k, 12))

        def line_for_point(e_keV: float, outs: dict, i: int) -> str:
            energy = _format_energy_no_unit(e_keV).rjust(12)
            vals: list[str] = []
            for k in keys_present:
                series = outs.get(k)
                if not isinstance(series, list) or i >= len(series):
                    vals.append("".rjust(12))
                    continue
                try:
                    v = float(series[i])
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
                    vals.append(self._format_sci(conv).rjust(12))
                elif k == "nucl_strag":
                    vals.append(self._format_sci(v).rjust(12))
                else:
                    unit_for_k = range_unit
                    if k == "long_strag":
                        unit_for_k = range_unit_long
                    elif k == "lat_strag":
                        unit_for_k = range_unit_lat
                    vals.append(_format_length_no_unit(v, unit_for_k).rjust(12))
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
        out_lines.append(f" Bragg Correction = {bragg_corr_pct:.2f}%")
        out_lines.append(" Unit Definitions (used in table below):")
        out_lines.append("   energy: keV")
        for k in keys_present:
            t1, t2 = title_by_key.get(k, (k, ""))
            out_lines.append(f"   {t1} {t2}: {unit_by_key.get(k, '')}")
        out_lines.append("")

        header_l1 = ["Ion"]
        header_l2 = ["Energy"]
        header_l3 = ["[keV]"]
        for k in keys_present:
            t1, t2 = title_by_key.get(k, (k, ""))
            header_l1.append(t1)
            header_l2.append(t2)
            header_l3.append(f"[{unit_by_key.get(k, '')}]")

        def _fmt_header(parts: list[str]) -> str:
            cells = []
            for idx, (_, width) in enumerate(columns):
                txt = parts[idx] if idx < len(parts) else ""
                cells.append(f"{txt:^{width}}")
            return "  " + "  ".join(cells)

        out_lines.append(_fmt_header(header_l1))
        out_lines.append(_fmt_header(header_l2))
        out_lines.append(_fmt_header(header_l3))
        out_lines.append("  " + "  ".join("-" * width for _, width in columns))

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

        solver_settings_btn = QPushButton("⚙ Solver Settings")
        solver_settings_btn.setToolTip("Open KORAL solver settings in Advanced Options")
        solver_settings_btn.clicked.connect(lambda: self.advanced_requested.emit("koral_solver"))

        layout.addWidget(self.model_button)
        layout.addWidget(self.selected_models_label)
        layout.addWidget(solver_settings_btn)
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

    def set_koral_solver_settings(self, settings: dict) -> None:
        if isinstance(settings, dict):
            self._koral_solver_settings = settings

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
        if hasattr(self, "spin_target_density"):
            output["target_density"] = self._get_target_density()
        if hasattr(self, "cmb_density_unit"):
            output["target_density_unit"] = self.cmb_density_unit.currentText()
        if hasattr(self, "chk_gas"):
            output["gas"] = bool(self.chk_gas.isChecked())
        if hasattr(self, "sw_koral_mode"):
            output["sw_koral_mode"] = bool(self.sw_koral_mode.isChecked())

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
                def _to_float_or_none(v):
                    if v is None or v == "":
                        return None
                    try:
                        return float(v)
                    except (TypeError, ValueError):
                        return None

                overrides = {k: _to_float_or_none(e.get(k)) for k in ("damage", "disp", "latt", "surf")}
                added = self._add_element_to_table(element, e.get("ratio", 0.0), overrides=overrides, refresh=False)
                if isinstance(added, dict):
                    added["mass_override"] = _to_float_or_none(e.get("mass_override"))
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

            if hasattr(self, "cmb_density_unit") and "target_density_unit" in output:
                idx = self.cmb_density_unit.findText(str(output["target_density_unit"]))
                if idx >= 0:
                    self.cmb_density_unit.setCurrentIndex(idx)
                    self._density_unit = self.cmb_density_unit.currentText()

            if hasattr(self, "spin_target_density") and "target_density" in output:
                try:
                    self._set_target_density(float(output["target_density"]))
                    self._density_user_override = True
                except (TypeError, ValueError):
                    pass

            if hasattr(self, "chk_gas") and "gas" in output:
                self.chk_gas.setChecked(bool(output.get("gas")))

            if hasattr(self, "sw_koral_mode") and "sw_koral_mode" in output:
                self.sw_koral_mode.setChecked(bool(output.get("sw_koral_mode")))

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
        # Keep consistent inner spacing from the frame on all sides.
        layout.setContentsMargins(8, 6, 8, 8)
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
        self.ion_mass.setDecimals(1)
        self.ion_mass.setReadOnly(False)
        self.ion_mass.setMaximumWidth(120)
        self.ion_mass.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.ion_mass, 1, 1, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        # Energy min/max (keV)
        grid.addWidget(QLabel("Energy min (keV)"), 0, 2, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.energy_min = QDoubleSpinBox()
        self.energy_min.setRange(0.0, 1e6)
        self.energy_min.setDecimals(1)
        self.energy_min.setValue(10.0)
        self.energy_min.setMaximumWidth(130)
        self.energy_min.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        grid.addWidget(self.energy_min, 1, 2, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        grid.addWidget(QLabel("Energy max (keV)"), 0, 3, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.energy_max = QDoubleSpinBox()
        self.energy_max.setRange(0.0, 1e6)
        self.energy_max.setDecimals(1)
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
        # Keep consistent inner spacing from the frame on all sides.
        v.setContentsMargins(8, 6, 8, 8)
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
            "", "Symbol", "Name", "Atomic\nNo.", "Weight\n(amu)",
            "Atom\nStoich", "Atom Stoich\n(%)", "Damage\n(eV)", "Disp\n(eV)", "Latt\n(eV)", "Surf\n(eV)"
        ])
        hdr = self.elem_table.horizontalHeader()
        if hdr is not None:
            # Keep columns readable and prevent header overflow beyond column bounds.
            hdr.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
            hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Fixed)
            # cols 1-6 remain Interactive so setColumnWidth() from the resize filter takes effect
            hdr.setSectionResizeMode(7, QHeaderView.ResizeMode.Fixed)
            hdr.setSectionResizeMode(8, QHeaderView.ResizeMode.Fixed)
            hdr.setSectionResizeMode(9, QHeaderView.ResizeMode.Fixed)
            hdr.setSectionResizeMode(10, QHeaderView.ResizeMode.Fixed)
            self.elem_table.setColumnWidth(0, 28)
            self.elem_table.setColumnWidth(1, 72)
            self.elem_table.setColumnWidth(2, 100)
            self.elem_table.setColumnWidth(3, 76)
            self.elem_table.setColumnWidth(4, 86)
            self.elem_table.setColumnWidth(5, 88)
            self.elem_table.setColumnWidth(6, 92)
            self.elem_table.setColumnWidth(7, 78)
            self.elem_table.setColumnWidth(8, 72)
            self.elem_table.setColumnWidth(9, 72)
            self.elem_table.setColumnWidth(10, 72)
            # Install proportional-resize filter: col 0 stays fixed, cols 1-6 scale
            _ProportionalColumnFilter(
                self.elem_table,
                fixed_cols={0: 28},
                stretch_weights={1: 72, 2: 100, 3: 76, 4: 86, 5: 88, 6: 92},
            )
            try:
                max_lines = max(1, max(h.count("\n") + 1 for h in [
                    "", "Symbol", "Name", "Atomic\nNo.", "Weight\n(amu)",
                    "Atom\nStoich", "Atom Stoich\n(%)", "Damage\n(eV)",
                    "Disp\n(eV)", "Latt\n(eV)", "Surf\n(eV)"
                ]))
                hdr.setMinimumHeight(hdr.fontMetrics().height() * max_lines + 12)
            except Exception:
                pass
            f = hdr.font()
            f.setBold(True)
            hdr.setFont(f)
        vh = self.elem_table.verticalHeader()
        if vh is not None:
            vh.setVisible(False)
        self.elem_table.setAlternatingRowColors(True)
        self.elem_table.setColumnHidden(7, True)
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

        # Compound correction + gas on one row (50% each); target density below.
        params_box = QWidget(box)
        params_v = QVBoxLayout(params_box)
        params_v.setContentsMargins(0, 0, 0, 0)
        params_v.setSpacing(4)

        inputs_row = QHBoxLayout()
        inputs_row.setContentsMargins(0, 0, 0, 0)
        inputs_row.setSpacing(12)

        # --- compound correction (left half) ---
        self.spin_compound_corr = QDoubleSpinBox()
        self.spin_compound_corr.setRange(0.0, 10.0)
        self.spin_compound_corr.setDecimals(1)
        self.spin_compound_corr.setSingleStep(0.01)
        self.spin_compound_corr.setValue(1.0)
        self.spin_compound_corr.setToolTip(
            "Compound stopping-power correction factor applied to the target."
        )
        cc_holder = QWidget()
        cc_l = QHBoxLayout(cc_holder)
        cc_l.setContentsMargins(0, 0, 0, 0)
        cc_l.setSpacing(4)
        cc_l.addWidget(QLabel("Compound correction"))
        cc_l.addWidget(self.spin_compound_corr, 1)

        # --- gas (right half) ---
        self.chk_gas = QCheckBox("Gas")
        self.chk_gas.setToolTip("Treat the target material as a gas.")
        gas_holder = QWidget()
        gas_l = QHBoxLayout(gas_holder)
        gas_l.setContentsMargins(0, 0, 0, 0)
        gas_l.setSpacing(4)
        gas_l.addWidget(self.chk_gas)
        gas_l.addStretch(1)

        inputs_row.addWidget(cc_holder, 1)
        inputs_row.addWidget(gas_holder, 1)
        params_v.addLayout(inputs_row)

        # --- target density (own row, with selectable unit) ---
        self.spin_target_density = QLineEdit()
        self.spin_target_density.setText("0.0")
        self.spin_target_density.textEdited.connect(lambda: setattr(self, "_density_user_override", True))
        self.spin_target_density.setToolTip(
            "Weighted density of the target compound.\n"
            "Auto-filled from the material database when elements are added.\n"
            "You can override this value manually.\n"
            "Supports scientific notation, e.g. 5.0e22."
        )
        self.cmb_density_unit = QComboBox()
        self.cmb_density_unit.addItems(["atoms/cm³", "g/cm³", "atoms/m³", "kg/m³"])
        self._density_unit = "atoms/cm³"
        self.cmb_density_unit.setToolTip("Unit for the target density value.")
        self.cmb_density_unit.currentIndexChanged.connect(self._on_density_unit_changed)

        td_holder = QWidget()
        td_l = QHBoxLayout(td_holder)
        td_l.setContentsMargins(0, 0, 0, 0)
        td_l.setSpacing(4)
        td_l.addWidget(QLabel("Target Density"))
        td_l.addWidget(self.spin_target_density, 1)
        td_l.addWidget(self.cmb_density_unit)
        params_v.addWidget(td_holder)

        v.addWidget(params_box)

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
                mass_text = f"{float(mass_val):.1f}"
            else:
                mass_raw = element.get("atomic_mass")
                try:
                    mass_text = f"{float(mass_raw):.1f}"
                except (TypeError, ValueError):
                    mass_text = str(mass_raw)
            mass_item = QTableWidgetItem(mass_text)
            mass_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
            self.elem_table.setItem(row, 4, mass_item)

            # Ratio (editable)
            ratio_item = QTableWidgetItem(f"{entry['ratio']:.1f}")
            ratio_item.setFlags(Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable | Qt.ItemFlag.ItemIsEnabled)
            self.elem_table.setItem(row, 5, ratio_item)

            percent = (entry["ratio"] / total_ratio * 100.0) if total_ratio else 0.0
            self.elem_table.setItem(row, 6, ro_item(f"{percent:.1f}"))

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
        action_l.setContentsMargins(4, 0, 0, 0)
        action_l.setSpacing(6)
        if hasattr(self, "elem_pick_btn") and self.elem_pick_btn is not None:
            action_l.addWidget(self.elem_pick_btn)
        action_l.addStretch(1)
        self.elem_table.setCellWidget(action_row, 0, action_cell)

        self._updating_elements_table = False

        # Update the weighted target density from the database.
        # Only auto-fill when the user has not manually edited the spinner
        # (we detect this by checking _density_user_override).
        if not getattr(self, "_density_user_override", False):
            self._update_density_from_elements()

    def _update_density_from_elements(self) -> None:
        """Compute stoichiometry-weighted atomic density and update the spinner."""
        if not hasattr(self, "spin_target_density"):
            return
        entries = self.element_entries
        total_ratio = sum(float(e.get("ratio", 0.0) or 0.0) for e in entries)
        if total_ratio <= 0.0 or not entries:
            return
        weighted_density = 0.0
        for e in entries:
            symbol = (e.get("element") or {}).get("symbol", "")
            w = float(e.get("ratio", 0.0) or 0.0) / total_ratio
            dens = _DENSITY_TABLE.get(symbol.upper(), 0.0)
            weighted_density += w * dens
        if weighted_density > 0.0:
            self._set_target_density(weighted_density)

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

        # Grid: left column = projectile range + straggling options,
        #        right column = nuclear / electron stopping (own column).
        #        col 6 = Plot/List switch (row 0 only, right-aligned).
        _CMB_W = 120  # uniform combo width
        opt_grid = QGridLayout()
        opt_grid.setHorizontalSpacing(8)
        opt_grid.setVerticalSpacing(4)
        opt_grid.setColumnStretch(0, 0)        # left checkbox
        opt_grid.setColumnStretch(1, 0)        # left combo
        opt_grid.setColumnMinimumWidth(2, 24)  # gap between the two columns
        opt_grid.setColumnStretch(2, 0)
        opt_grid.setColumnStretch(3, 0)        # right checkbox
        opt_grid.setColumnStretch(4, 0)        # right combo
        opt_grid.setColumnStretch(5, 1)        # takes all extra space
        opt_grid.setColumnStretch(6, 0)        # switch (row 0)

        # --- left column: projectile range + straggling ---
        self.chk_prange = QCheckBox("Projectile Range")
        self.cmb_prange = QComboBox()
        self.cmb_prange.setFixedWidth(_CMB_W)
        self.cmb_prange.addItems(self.state.unit_options)
        opt_grid.addWidget(self.chk_prange, 0, 0)
        opt_grid.addWidget(self.cmb_prange, 0, 1)
        self._output_option_widgets["prange"] = [self.chk_prange, self.cmb_prange]

        # Plot / List switch – below all output checkboxes, right-aligned
        self.sw_koral_mode = ToggleSwitch()
        self.sw_koral_mode.setChecked(False)  # False = Plot, True = List
        sw_container = QWidget()
        sw_layout = QHBoxLayout(sw_container)
        sw_layout.setContentsMargins(0, 0, 8, 0)
        sw_layout.setSpacing(4)
        sw_layout.addWidget(QLabel("Plot"))
        sw_layout.addWidget(self.sw_koral_mode)
        sw_layout.addWidget(QLabel("List"))
        opt_grid.addWidget(sw_container, 3, 0, 1, 2, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.sw_koral_mode.toggled.connect(self._update_koral_plot_view)

        # "All" toggle, directly under nuclear/electronic stopping
        self.all_none_chk = QCheckBox("All")
        self.all_none_chk.setTristate(True)
        self.all_none_chk.setCheckState(Qt.CheckState.Unchecked)
        self.all_none_chk.clicked.connect(self._toggle_all_options)
        opt_grid.addWidget(
            self.all_none_chk, 2, 3, 1, 2,
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
        )

        self.chk_long_strag = QCheckBox("Long. Straggling (σ_x)")
        self.cmb_long_strag = QComboBox()
        self.cmb_long_strag.setFixedWidth(_CMB_W)
        self.cmb_long_strag.addItems(self.state.unit_options)
        opt_grid.addWidget(self.chk_long_strag, 1, 0)
        opt_grid.addWidget(self.cmb_long_strag, 1, 1)
        self._output_option_widgets["long_strag"] = [self.chk_long_strag, self.cmb_long_strag]

        self.chk_lat_strag = QCheckBox("Lat. Straggling (σ_z)")
        self.cmb_lat_strag = QComboBox()
        self.cmb_lat_strag.setFixedWidth(_CMB_W)
        self.cmb_lat_strag.addItems(self.state.unit_options)
        opt_grid.addWidget(self.chk_lat_strag, 2, 0)
        opt_grid.addWidget(self.cmb_lat_strag, 2, 1)
        self._output_option_widgets["lat_strag"] = [self.chk_lat_strag, self.cmb_lat_strag]

        # --- right column: nuclear / electron stopping ---
        _STOP_UNITS = [
            "eV/Å",
            "keV/µm",
            "MeV/mm",
            "keV/(µg/cm²)",
            "MeV/(mg/cm²)",
            "keV/(mg/cm²)",
            "eV/(10¹⁵ atoms/cm²)",
            "L.S.S. reduced units",
        ]

        self.chk_nucl_strag = QCheckBox("Nuclear Stopping")
        self.cmb_nucl_stop_unit = QComboBox()
        self.cmb_nucl_stop_unit.setFixedWidth(_CMB_W)
        self.cmb_nucl_stop_unit.addItems(_STOP_UNITS)
        opt_grid.addWidget(self.chk_nucl_strag, 0, 3)
        opt_grid.addWidget(self.cmb_nucl_stop_unit, 0, 4)
        self._output_option_widgets["nucl_stop"] = [self.chk_nucl_strag, self.cmb_nucl_stop_unit]

        self.chk_elec_hop = QCheckBox("Electron Stopping")
        self.cmb_elec_stop_unit = QComboBox()
        self.cmb_elec_stop_unit.setFixedWidth(_CMB_W)
        self.cmb_elec_stop_unit.addItems(_STOP_UNITS)
        opt_grid.addWidget(self.chk_elec_hop, 1, 3)
        opt_grid.addWidget(self.cmb_elec_stop_unit, 1, 4)
        self._output_option_widgets["elec_stop"] = [self.chk_elec_hop, self.cmb_elec_stop_unit]

        for checkbox in self._koral_output_checkboxes():
            checkbox.toggled.connect(self._sync_all_options_checkbox)

        self._sync_all_options_checkbox()

        v.addLayout(opt_grid)

        return box

    def _koral_output_checkboxes(self) -> list[QCheckBox]:
        return [
            self.chk_prange,
            self.chk_long_strag,
            self.chk_lat_strag,
            self.chk_nucl_strag,
            self.chk_elec_hop,
        ]

    def _sync_all_options_checkbox(self) -> None:
        checkboxes = [c for c in self._koral_output_checkboxes() if c.isEnabled() and c.isVisible()]
        if not checkboxes:
            checkboxes = self._koral_output_checkboxes()
        all_checked = all(c.isChecked() for c in checkboxes)
        any_checked = any(c.isChecked() for c in checkboxes)
        self.all_none_chk.blockSignals(True)
        self.all_none_chk.setCheckState(
            Qt.CheckState.Checked if all_checked else
            Qt.CheckState.PartiallyChecked if any_checked else
            Qt.CheckState.Unchecked
        )
        self.all_none_chk.blockSignals(False)

    def _toggle_all_options(self) -> None:
        if self.all_none_chk.checkState() == Qt.CheckState.PartiallyChecked:
            self.all_none_chk.setCheckState(Qt.CheckState.Checked)
        new_state = self.all_none_chk.checkState() == Qt.CheckState.Checked
        for checkbox in self._koral_output_checkboxes():
            if checkbox.isEnabled() and checkbox.isVisible():
                checkbox.blockSignals(True)
                checkbox.setChecked(new_state)
                checkbox.blockSignals(False)
        self._sync_all_options_checkbox()

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
        ax.set_title("KORAL results")
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Value")
        self.figure.tight_layout()

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
        btn_export_csv = QPushButton("Export CSV")
        btn_export_csv.clicked.connect(self._export_results_csv)
        export_row.addWidget(btn_export_csv)
        btn_export = QPushButton("Export TXT")
        btn_export.clicked.connect(self._export_results_txt)
        export_row.addWidget(btn_export)
        lv.addLayout(export_row)

        self.koral_result_table = QTableWidget(0, 3)
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
