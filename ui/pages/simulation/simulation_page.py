import sys
import json
import base64
import time
from pathlib import Path
from typing import Dict, List

from PyQt6.QtCore import Qt, QByteArray, QPoint, pyqtSignal, QEvent
from PyQt6.QtGui import QAction, QMouseEvent
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QDockWidget, QVBoxLayout, QLabel,
    QSplitter, QTableWidget, QTableWidgetItem,
    QFormLayout, QLineEdit, QSpinBox, QDoubleSpinBox, QToolBar, QPushButton,
    QMessageBox, QToolButton, QHBoxLayout, QMenu, QInputDialog, QHeaderView,
    QFrame, QGridLayout, QStyle, QDialog, QListWidget, QDialogButtonBox, QTabWidget
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# --- Directory for views ---
BASE_DIR = Path(__file__).resolve().parent
VIEWS_DIR = BASE_DIR / "views"
VIEWS_DIR.mkdir(exist_ok=True)


class LayoutStore:
    def __init__(self, directory: Path):
        self.directory = directory
        self.order_file = self.directory / ".views_order.json"
        self.views: Dict[str, dict] = {}
        self.order: List[str] = []
        self.reload()

    def reload(self):
        self.views.clear()
        for file in self.directory.glob("*.json"):
            try:
                self.views[file.stem] = json.loads(file.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
        stored_order = []
        if self.order_file.exists():
            try:
                stored_order = json.loads(self.order_file.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                stored_order = []
        existing = list(self.views.keys())
        ordered = [name for name in stored_order if name in self.views]
        leftovers = [name for name in existing if name not in ordered]
        ordered.extend(sorted(leftovers))
        self.order = ordered

    def names(self):
        return list(self.order)

    def get(self, name: str):
        return self.views.get(name)

    def save(self, name: str, payload: dict):
        (self.directory / f"{name}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        self.views[name] = payload
        if name not in self.order:
            self.order.append(name)
        self._persist_order()

    def remove(self, name: str):
        path = self.directory / f"{name}.json"
        if path.exists():
            path.unlink()
        self.views.pop(name, None)
        if name in self.order:
            self.order.remove(name)
        self._persist_order()

    def reorder(self, names: List[str]):
        filtered = [n for n in names if n in self.views]
        leftovers = [n for n in self.order if n not in filtered]
        self.order = filtered + leftovers
        self._persist_order()

    def apply_management(self, ordered_names: List[str], removed: List[str]):
        for name in removed:
            self.remove(name)
        self.reorder(ordered_names)

    def _persist_order(self):
        self.order_file.write_text(json.dumps(self.order, indent=2), encoding="utf-8")


# ----------------------------- Panels ----------------------------- #

class StandardPanel(QWidget):
    fullscreen_requested = pyqtSignal(object)

    def __init__(self):
        super().__init__()
        self.content_layout = QVBoxLayout(self)
        self.content_layout.setContentsMargins(0, 0, 0, 0)
        self.content_layout.setSpacing(8)
        self._fullscreen_toggle_enabled = True
        self._tracked_widgets = set()
        self._register_child_for_fullscreen(self)

    def set_fullscreen_toggle_enabled(self, enabled: bool):
        self._fullscreen_toggle_enabled = enabled

    def _register_child_for_fullscreen(self, widget: QWidget):
        if widget in self._tracked_widgets:
            return
        widget.installEventFilter(self)
        self._tracked_widgets.add(widget)
        for child in widget.findChildren(QWidget):
            self._register_child_for_fullscreen(child)

    def _unregister_child_for_fullscreen(self, widget: QWidget):
        if widget in self._tracked_widgets:
            widget.removeEventFilter(self)
            self._tracked_widgets.remove(widget)
        for child in widget.findChildren(QWidget):
            self._unregister_child_for_fullscreen(child)

    def childEvent(self, event: QEvent):
        if event.type() == QEvent.Type.ChildAdded:
            child = event.child()
            if isinstance(child, QWidget):
                self._register_child_for_fullscreen(child)
        elif event.type() == QEvent.Type.ChildRemoved:
            child = event.child()
            if isinstance(child, QWidget):
                self._unregister_child_for_fullscreen(child)
        super().childEvent(event)

    def eventFilter(self, obj, event):
        if (
            self._fullscreen_toggle_enabled
            and event.type() == QEvent.Type.MouseButtonDblClick
            and isinstance(event, QMouseEvent)
            and event.button() == Qt.MouseButton.LeftButton
        ):
            self.fullscreen_requested.emit(self)
            return True
        return super().eventFilter(obj, event)

    def mouseDoubleClickEvent(self, event: QMouseEvent):
        if self._fullscreen_toggle_enabled and event.button() == Qt.MouseButton.LeftButton:
            self.fullscreen_requested.emit(self)
        super().mouseDoubleClickEvent(event)


class UnconfiguredPanel(StandardPanel):
    def __init__(self):
        super().__init__()
        icon = self.style().standardIcon(QStyle.StandardPixmap.SP_MessageBoxInformation)
        pic = QLabel(); pic.setPixmap(icon.pixmap(48, 48)); pic.setAlignment(Qt.AlignmentFlag.AlignCenter)
        text = QLabel("This panel is not configured yet.\nUse ⚙ to select a view.")
        text.setAlignment(Qt.AlignmentFlag.AlignCenter); text.setWordWrap(True)
        self.content_layout.addStretch(); self.content_layout.addWidget(pic); self.content_layout.addWidget(text); self.content_layout.addStretch()


class IonPanel(StandardPanel):
    def __init__(self):
        super().__init__()
        form = QFormLayout()
        form.addRow("Ion Type", QLineEdit())
        e = QDoubleSpinBox(); e.setRange(0, 1e6); form.addRow("Energy [keV]", e)
        a = QDoubleSpinBox(); a.setRange(0, 90); form.addRow("Angle [°]", a)
        self.content_layout.addLayout(form); self.content_layout.addStretch()


class PlotControlPanel(StandardPanel):
    def __init__(self):
        super().__init__()
        self.content_layout.addWidget(QLabel("Plot Controls"))
        self.content_layout.addStretch()


class CalcPanel(StandardPanel):
    def __init__(self):
        super().__init__()
        form = QFormLayout()
        form.addRow("Backscattered", QSpinBox())
        form.addRow("Transmitted", QSpinBox())
        form.addRow("Vacancies/Ion", QSpinBox())
        self.content_layout.addLayout(form)
        self.content_layout.addStretch()


class TargetDataTable(QTableWidget):
    HEADERS = ["Layer", "Name", "Width (Å)", "Density", "Material", "Stop Corr."]
    def __init__(self):
        super().__init__(0, len(self.HEADERS))
        self.setHorizontalHeaderLabels(self.HEADERS)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.verticalHeader().setVisible(False)
        self.setShowGrid(True)
        self.setGridStyle(Qt.PenStyle.SolidLine)
        self.setFrameStyle(QFrame.Shape.Box | QFrame.Shadow.Plain)
        self.setLineWidth(1)
        self.setStyleSheet("""
            QTableWidget {
                gridline-color: black;
                background: white;
                border: 1px solid black;
            }
            QTableWidget::item { border: 1px solid black; }
        """)
        for row_data in [
            ["1", "Be", "100", "1.85", "Solid", "1"],
            ["2", "Au", "100", "19.32", "Solid", "1"],
            ["3", "Si", "100", "2.32", "Solid", "1"],
            ["4", "Vacuum", "0", "0", "-", "0"],
        ]:
            self._insert_row_with_values(self.rowCount(), row_data)

    def _insert_row_with_values(self, row: int, values=None):
        self.insertRow(row)
        for col in range(self.columnCount()):
            value = values[col] if values and col < len(values) else ""
            self.setItem(row, col, QTableWidgetItem(value))


class TargetTablePanel(StandardPanel):
    def __init__(self):
        super().__init__()
        self.content_layout.addWidget(TargetDataTable())


class ProjectionViewPanel(StandardPanel):
    LABELS = ["XY", "XZ", "YZ", "3D"]

    def __init__(self):
        super().__init__()
        self.tabs = QTabWidget()  # Store reference to tabs
        for label in self.LABELS:
            figure = Figure(figsize=(4, 3))
            if label == "3D":
                ax = figure.add_subplot(111, projection="3d")
                ax.plot([0, 1], [0, 1], [0, 1])
                ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
            else:
                ax = figure.add_subplot(111)
                ax.plot([0, 1], [0, 1])
                ax.set_xlabel(label[0]); ax.set_ylabel(label[1])
            canvas = FigureCanvas(figure)
            page = QWidget()
            page_layout = QVBoxLayout(page)
            page_layout.setContentsMargins(0, 0, 0, 0)
            page_layout.addWidget(canvas)
            self.tabs.addTab(page, label)
        self.content_layout.addWidget(self.tabs)

    def get_current_tab_index(self) -> int:
        return self.tabs.currentIndex()

    def set_current_tab_index(self, index: int):
        if 0 <= index < self.tabs.count():
            self.tabs.setCurrentIndex(index)

PANEL_TYPES = {
    "Unconfigured": UnconfiguredPanel,
    "Ion Parameters": IonPanel,
    "Plot Controls": PlotControlPanel,
    "Calculation Parameters": CalcPanel,
    "Target Table": TargetTablePanel,
    "Projection Views": ProjectionViewPanel,
}

# ------------------ Configurable Panel Container ------------------- #

class ConfigurablePanel(QFrame):
    fullscreen_toggle_requested = pyqtSignal(object)

    def __init__(self, panel_key: str = "Unconfigured"):
        super().__init__()
        self.setObjectName("ConfigurablePanel")
        self.setStyleSheet("#ConfigurablePanel { border: 2px solid black; border-radius: 6px; }")
        self.current_panel_key = panel_key
        self.panel_state = {}  # Store panel-specific state

        header = QWidget()
        header.setStyleSheet("background: #f5f5f5; border-bottom: 1px solid black;")
        h = QHBoxLayout(header)
        h.setContentsMargins(8, 4, 8, 4)
        self.title_label = QLabel(panel_key)
        self.title_label.setStyleSheet("font-weight: 600;")
        h.addWidget(self.title_label)
        h.addStretch()
        self.gear_button = QToolButton()
        self.gear_button.setText("⚙")
        self.gear_button.clicked.connect(self._show_menu)
        h.addWidget(self.gear_button)

        container = QWidget()
        self.body_layout = QVBoxLayout(container)
        self.body_layout.setContentsMargins(8, 8, 8, 8)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(header)
        outer.addWidget(container)

        self.panel_widget = None
        self.set_panel(panel_key)

    def _show_menu(self):
        menu = QMenu(self)
        for key in [k for k in PANEL_TYPES if k != "Unconfigured"]:
            act = menu.addAction(key)
            act.setCheckable(True)
            act.setChecked(key == self.current_panel_key)
        action = menu.exec(self.gear_button.mapToGlobal(self.gear_button.rect().bottomRight()))
        if action:
            self.set_panel(action.text())

    def _bind_fullscreen_signal(self, panel_widget):
        if hasattr(panel_widget, "fullscreen_requested"):
            panel_widget.fullscreen_requested.connect(self._handle_panel_fullscreen_request)

    def _unbind_fullscreen_signal(self, panel_widget):
        if hasattr(panel_widget, "fullscreen_requested"):
            try:
                panel_widget.fullscreen_requested.disconnect(self._handle_panel_fullscreen_request)
            except TypeError:
                pass

    def _handle_panel_fullscreen_request(self, _):
        self.fullscreen_toggle_requested.emit(self)

    def mousePressEvent(self, event):
        super().mousePressEvent(event)

    def set_panel(self, panel_key: str, state: dict = None):
        cls = PANEL_TYPES.get(panel_key, UnconfiguredPanel)
        if self.panel_widget:
            self._unbind_fullscreen_signal(self.panel_widget)
            self.body_layout.removeWidget(self.panel_widget)
            self.panel_widget.deleteLater()
        self.panel_widget = cls()
        self._bind_fullscreen_signal(self.panel_widget)
        self.body_layout.addWidget(self.panel_widget)
        self.current_panel_key = panel_key
        self.title_label.setText(panel_key)
        self.panel_state = state or {}
        
        # Apply state if available
        if state and panel_key == "Projection Views":
            if isinstance(self.panel_widget, ProjectionViewPanel) and "tab_index" in state:
                self.panel_widget.set_current_tab_index(state["tab_index"])

    def get_panel_state(self) -> dict:
        """Get current state of the panel for serialization."""
        state = {}
        if self.current_panel_key == "Projection Views" and isinstance(self.panel_widget, ProjectionViewPanel):
            state["tab_index"] = self.panel_widget.get_current_tab_index()
        return state

    def serialize(self) -> dict:
        """Serialize panel configuration including state."""
        return {
            "key": self.current_panel_key,
            "state": self.get_panel_state(),
        }


class PanelGrid(QWidget):
    def __init__(self, defaults=None):
        super().__init__()
        self.columns: list[dict] = []
        self.horizontal_splitter = QSplitter(Qt.Orientation.Horizontal)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.horizontal_splitter)
        self.fullscreen_container = QWidget()
        self.fullscreen_container.setVisible(False)
        self.fullscreen_layout = QVBoxLayout(self.fullscreen_container)
        self.fullscreen_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.fullscreen_container)
        self.fullscreen_info = None

        defaults = defaults or [
            ["Target Table", "Ion Parameters"],
            ["Plot Controls"],
            ["Calculation Parameters"],
        ]
        self._build_from_columns(defaults)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)

    def _clear_columns(self):
        while self.horizontal_splitter.count():
            self.horizontal_splitter.widget(0).setParent(None)
        for column in self.columns:
            column["splitter"].setParent(None)
            for panel in column["panels"]:
                panel.deleteLater()
        self.columns.clear()

    def _build_from_columns(self, columns_data):
        self._clear_columns()
        if not columns_data:
            columns_data = [["Unconfigured"]]
        normalized = []
        for entry in columns_data:
            if isinstance(entry, dict):
                normalized.append(entry)
            else:
                normalized.append({"panels": entry or ["Unconfigured"], "splitter": None})
        for col_def in normalized:
            panels_info = col_def.get("panels") or ["Unconfigured"]
            splitter_state = col_def.get("splitter")
            self._create_column(panels_info, splitter_state)

    def _wire_panel(self, panel: ConfigurablePanel):
        panel.fullscreen_toggle_requested.connect(self._toggle_fullscreen)

    def _create_column(self, panels_info, splitter_state: str | None = None):
        column_splitter = QSplitter(Qt.Orientation.Vertical)
        panels = []
        for panel_data in panels_info:
            # Handle both old format (string) and new format (dict)
            if isinstance(panel_data, str):
                key = panel_data
                state = None
            elif isinstance(panel_data, dict):
                key = panel_data.get("key", "Unconfigured")
                state = panel_data.get("state")
            else:
                key = "Unconfigured"
                state = None
            
            panel = ConfigurablePanel(key)
            if state:
                panel.set_panel(key, state)
            self._wire_panel(panel)
            panels.append(panel)
            column_splitter.addWidget(panel)
        self.columns.append({"splitter": column_splitter, "panels": panels})
        self.horizontal_splitter.addWidget(column_splitter)
        if splitter_state:
            column_splitter.restoreState(QByteArray.fromBase64(splitter_state.encode()))

    def _panel_position(self, pos: QPoint):
        widget = self.childAt(pos)
        if isinstance(widget, ConfigurablePanel):
            for c, col in enumerate(self.columns):
                for r, panel in enumerate(col["panels"]):
                    if panel is widget:
                        return r, c
        return 0, 0

    def _column_at_pos(self, pos: QPoint) -> int:
        if not self.columns:
            return 0
        for idx, column in enumerate(self.columns):
            if column["splitter"].geometry().contains(pos):
                return idx
        return max(0, len(self.columns) - 1)

    def _row_at_pos(self, column_idx: int, pos: QPoint) -> int:
        column = self.columns[column_idx]
        local = column["splitter"].mapFrom(self, pos)
        for idx, panel in enumerate(column["panels"]):
            if panel.geometry().contains(local):
                return idx
        return len(column["panels"]) - 1

    def _show_context_menu(self, pos):
        if not self.columns:
            col_idx = 0
            row_idx = 0
        else:
            col_idx = self._column_at_pos(pos)
            row_idx = self._row_at_pos(col_idx, pos)
        menu = QMenu(self)
        add_row_above = menu.addAction("Add Row Above")
        add_row_below = menu.addAction("Add Row Below")
        del_row = menu.addAction("Delete Row")
        add_col_before = menu.addAction("Add Column Before")
        add_col_after = menu.addAction("Add Column After")
        del_col = menu.addAction("Delete Column")
        action = menu.exec(self.mapToGlobal(pos))
        if action == add_row_above:
            self._insert_row(col_idx, row_idx)
        elif action == add_row_below:
            self._insert_row(col_idx, row_idx + 1)
        elif action == del_row:
            self._remove_row(col_idx, row_idx)
        elif action == add_col_before:
            template = self.columns[col_idx]["panels"] if self.columns else []
            keys = ["Unconfigured"] * len(template) or ["Unconfigured"]
            self._insert_col(col_idx, keys)
        elif action == add_col_after:
            template = self.columns[col_idx]["panels"] if self.columns else []
            keys = ["Unconfigured"] * len(template) or ["Unconfigured"]
            self._insert_col(col_idx + 1, keys)
        elif action == del_col:
            self._remove_col(col_idx)

    def _insert_row(self, column_idx, index):
        if not self.columns:
            self._insert_col(0, ["Unconfigured"])
        column_idx = max(0, min(column_idx, len(self.columns) - 1))
        column = self.columns[column_idx]
        panels = column["panels"]
        splitter = column["splitter"]
        idx = max(0, min(index, len(panels)))
        panel = ConfigurablePanel("Unconfigured")
        self._wire_panel(panel)
        panels.insert(idx, panel)
        splitter.insertWidget(idx, panel)

    def _remove_row(self, column_idx, index):
        if not self.columns:
            return
        column_idx = max(0, min(column_idx, len(self.columns) - 1))
        column = self.columns[column_idx]
        if len(column["panels"]) <= 1:
            return
        idx = max(0, min(index, len(column["panels"]) - 1))
        panel = column["panels"][idx]
        if self.fullscreen_info and self.fullscreen_info["panel"] is panel:
            self._exit_fullscreen()
        column["panels"].pop(idx)
        panel.deleteLater()

    def _insert_col(self, index, keys):
        column_splitter = QSplitter(Qt.Orientation.Vertical)
        panels = []
        for key in (keys or ["Unconfigured"]):
            panel = ConfigurablePanel(key)
            self._wire_panel(panel)
            panels.append(panel)
            column_splitter.addWidget(panel)
        idx = max(0, min(index, len(self.columns)))
        self.columns.insert(idx, {"splitter": column_splitter, "panels": panels})
        self.horizontal_splitter.insertWidget(idx, column_splitter)

    def _remove_col(self, index):
        if self.column_count() <= 1 or index >= len(self.columns):
            return
        column = self.columns[index]
        if self.fullscreen_info and self.fullscreen_info["column_idx"] == index:
            self._exit_fullscreen()
        self.columns.pop(index)
        column["splitter"].setParent(None)
        for panel in column["panels"]:
            panel.deleteLater()

    def _toggle_fullscreen(self, panel: ConfigurablePanel):
        if self.fullscreen_info and self.fullscreen_info["panel"] is panel:
            self._exit_fullscreen()
        else:
            self._enter_fullscreen(panel)

    def _find_panel_location(self, panel):
        for c, col in enumerate(self.columns):
            for r, p in enumerate(col["panels"]):
                if p is panel:
                    return c, r
        return None

    def _enter_fullscreen(self, panel: ConfigurablePanel):
        if self.fullscreen_info:
            self._exit_fullscreen()
        location = self._find_panel_location(panel)
        if location is None:
            return
        col_idx, row_idx = location
        column = self.columns[col_idx]
        column["panels"].pop(row_idx)
        self.fullscreen_info = {
            "panel": panel,
            "column_idx": col_idx,
            "row_idx": row_idx,
            "panel_key": panel.current_panel_key,
            "h_state": QByteArray(self.horizontal_splitter.saveState()),
            "v_states": [QByteArray(col["splitter"].saveState()) for col in self.columns],
        }
        self.fullscreen_layout.addWidget(panel)
        self.horizontal_splitter.setVisible(False)
        self.fullscreen_container.setVisible(True)

    def _exit_fullscreen(self):
        if not self.fullscreen_info:
            return
        info = self.fullscreen_info
        panel = info["panel"]
        column = self.columns[info["column_idx"]]
        insert_idx = max(0, min(info["row_idx"], len(column["panels"])))
        column["panels"].insert(insert_idx, panel)
        column["splitter"].insertWidget(insert_idx, panel)
        self.fullscreen_layout.removeWidget(panel)
        self.fullscreen_container.setVisible(False)
        self.horizontal_splitter.setVisible(True)
        if info.get("h_state"):
            self.horizontal_splitter.restoreState(info["h_state"])
        for col, state in zip(self.columns, info.get("v_states", [])):
            col["splitter"].restoreState(state)
        self.fullscreen_info = None

    def row_count(self):
        return max((len(col["panels"]) for col in self.columns), default=0)

    def column_count(self):
        return len(self.columns)

    def serialize(self):
        data_columns = []
        for idx, column in enumerate(self.columns):
            panels_data = []
            for panel in column["panels"]:
                panels_data.append(panel.serialize())
            
            # Handle fullscreen panel
            if self.fullscreen_info and self.fullscreen_info["column_idx"] == idx:
                insert_idx = max(0, min(self.fullscreen_info["row_idx"], len(panels_data)))
                fullscreen_panel = self.fullscreen_info["panel"]
                panels_data.insert(insert_idx, fullscreen_panel.serialize())
            
            data_columns.append({
                "panels": panels_data,
                "splitter": base64.b64encode(bytes(column["splitter"].saveState())).decode(),
            })
        return {
            "splitter": base64.b64encode(bytes(self.horizontal_splitter.saveState())).decode(),
            "columns": data_columns,
        }

    def apply_layout(self, data: dict):
        if self.fullscreen_info:
            self._exit_fullscreen()
        columns_data = data.get("columns")
        if columns_data is not None:
            self._build_from_columns(columns_data)
            state = data.get("splitter")
            if state:
                self.horizontal_splitter.restoreState(QByteArray.fromBase64(state.encode()))
            return
        legacy = data.get("panels", [])
        rebuilt = list(zip(*legacy)) if legacy else [["Unconfigured"]]
        self._build_from_columns([list(filter(None, col)) or ["Unconfigured"] for col in rebuilt])
        if "h_splitter" in data:
            self.horizontal_splitter.restoreState(QByteArray.fromBase64(data["h_splitter"].encode()))
        for column, state in zip(self.columns, data.get("v_splitters", [])):
            column["splitter"].restoreState(QByteArray.fromBase64(state.encode()))

# ------------------------------ MC Results Widget (Embeddable) ------------------------------ #

class MCResultsWidget(QWidget):
    """
    Embeddable widget for MC Results that can be used in other applications.
    Contains a PanelGrid with toolbar for view management.
    """
    def __init__(self, parent=None, show_toolbar=True):
        super().__init__(parent)
        self.layout_store = LayoutStore(VIEWS_DIR)
        self.show_toolbar = show_toolbar
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(4)
        
        # Panel Grid
        self.panel_grid = PanelGrid()
        main_layout.addWidget(self.panel_grid)
        
        # Load default view if exists
        if "Default" in self.layout_store.views:
            self._apply_view("Default")

    def get_view_names(self) -> List[str]:
        """Return list of available view names."""
        return self.layout_store.names()

    def save_layout(self, name: str):
        """Save current layout with given name."""
        payload = {
            "panel_grid": self.panel_grid.serialize(),
        }
        self.layout_store.save(name.strip(), payload)

    def _apply_view(self, name: str):
        data = self.layout_store.get(name)
        if not data:
            return False
        try:
            if "panel_grid" in data:
                self.panel_grid.apply_layout(data["panel_grid"])
            return True
        except Exception:
            return False

    def apply_view(self, name: str) -> bool:
        """Apply a saved view by name. Returns True on success."""
        return self._apply_view(name)

    def open_manage_views_dialog(self, parent=None):
        """Open the manage views dialog."""
        dialog = ManageViewsDialog(self.layout_store.names(), parent or self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.layout_store.apply_management(dialog.current_order(), dialog.deleted_names())
            return True
        return False


# ------------------------------ Main Window ------------------------------ #

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SRIM UI – Modern Version")

        # Use a normal size of 1366x768, start maximized, stay resizable.
        screen = self.screen() or (self.windowHandle().screen() if self.windowHandle() else None)
        avail = screen.availableGeometry() if screen else None
        base_w, base_h = 1366, 768
        min_w = min(base_w, avail.width()) if avail else base_w
        min_h = min(base_h, avail.height()) if avail else base_h
        self.setMinimumSize(min_w, min_h)
        self.resize(base_w, base_h)
        QTimer.singleShot(0, self.showMaximized)

        self.layout_store = LayoutStore(VIEWS_DIR)
        self.panel_grid = PanelGrid()
        self.setCentralWidget(self.panel_grid)
        self.create_menu()
        self.create_toolbar()
        if "Default" in self.layout_store.views:
            self.apply_view("Default")

    # -------------------- Menü & Toolbar --------------------

    def create_menu(self):
        menubar = self.menuBar()
        self.view_menu = menubar.addMenu("View")
        self.save_view_act = QAction("Save Layout…", self)
        self.save_view_act.triggered.connect(self.save_layout)
        self.manage_views_act = QAction("Manage Views…", self)
        self.manage_views_act.triggered.connect(self._open_manage_views)
        self.refresh_view_menu()

    def refresh_view_menu(self):
        self.view_menu.clear()
        self.view_menu.addAction(self.save_view_act)
        self.view_menu.addSeparator()
        for name in self.layout_store.names():
            act = QAction(name, self)
            act.triggered.connect(lambda _, n=name: self.apply_view(n))
            self.view_menu.addAction(act)
        self.view_menu.addSeparator()
        self.view_menu.addAction(self.manage_views_act)

    def create_toolbar(self):
        tb = QToolBar("Main")
        tb.setMovable(False)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, tb)

        for text in ["Animate", "Continue", "Change TRIM"]:
            btn = QPushButton(text)
            tb.addWidget(btn)

    # -------------------- Layouts speichern/laden --------------------

    def save_layout(self):
        name, ok = QInputDialog.getText(self, "Save Layout", "Layout name:")
        if not ok or not name.strip():
            return
        payload = {
            "geometry": base64.b64encode(bytes(self.saveGeometry())).decode(),
            "panel_grid": self.panel_grid.serialize(),
        }
        self.layout_store.save(name.strip(), payload)
        self.refresh_view_menu()

    def apply_view(self, name):
        data = self.layout_store.get(name)
        if not data:
            QMessageBox.warning(self, "Layout", f"View '{name}' not found.")
            return
        try:
            if "geometry" in data:
                self.restoreGeometry(QByteArray.fromBase64(data["geometry"].encode()))
            if "panel_grid" in data:
                self.panel_grid.apply_layout(data["panel_grid"])
        except Exception:
            QMessageBox.warning(self, "Layout", f"View '{name}' could not be loaded.")

    def _open_manage_views(self):
        dialog = ManageViewsDialog(self.layout_store.names(), self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.layout_store.apply_management(dialog.current_order(), dialog.deleted_names())
            self.refresh_view_menu()


# ------------------------------ Manage Views Dialog ------------------------------ #

class ManageViewsDialog(QDialog):
    def __init__(self, names: List[str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Manage Views")
        self.deleted = set()

        layout = QVBoxLayout(self)
        self.list_widget = QListWidget()
        self.list_widget.addItems(names)
        layout.addWidget(self.list_widget)

        controls = QHBoxLayout()
        self.up_btn = QPushButton("Move Up")
        self.down_btn = QPushButton("Move Down")
        self.del_btn = QPushButton("Delete")
        controls.addWidget(self.up_btn)
        controls.addWidget(self.down_btn)
        controls.addWidget(self.del_btn)
        layout.addLayout(controls)

        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        layout.addWidget(button_box)

        self.up_btn.clicked.connect(self._move_up)
        self.down_btn.clicked.connect(self._move_down)
        self.del_btn.clicked.connect(self._delete_selected)
        self.list_widget.currentRowChanged.connect(self._update_buttons)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        self._update_buttons()

    def _move_up(self):
        row = self.list_widget.currentRow()
        if row > 0:
            item = self.list_widget.takeItem(row)
            self.list_widget.insertItem(row - 1, item)
            self.list_widget.setCurrentRow(row - 1)

    def _move_down(self):
        row = self.list_widget.currentRow()
        if 0 <= row < self.list_widget.count() - 1:
            item = self.list_widget.takeItem(row)
            self.list_widget.insertItem(row + 1, item)
            self.list_widget.setCurrentRow(row + 1)

    def _delete_selected(self):
        row = self.list_widget.currentRow()
        if row >= 0:
            item = self.list_widget.takeItem(row)
            self.deleted.add(item.text())
            self.list_widget.setCurrentRow(min(row, self.list_widget.count() - 1))
        self._update_buttons()

    def _update_buttons(self):
        row = self.list_widget.currentRow()
        count = self.list_widget.count()
        self.up_btn.setEnabled(row > 0)
        self.down_btn.setEnabled(0 <= row < count - 1)
        self.del_btn.setEnabled(row >= 0)

    def current_order(self):
        return [self.list_widget.item(i).text() for i in range(self.list_widget.count())]

    def deleted_names(self):
        return list(self.deleted)


# ------------------------------ Main ------------------------------

def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    pal = app.palette()
    pal.setColor(pal.ColorRole.Window, Qt.GlobalColor.white)
    pal.setColor(pal.ColorRole.Base, Qt.GlobalColor.white)
    pal.setColor(pal.ColorRole.Text, Qt.GlobalColor.black)
    pal.setColor(pal.ColorRole.WindowText, Qt.GlobalColor.black)
    app.setPalette(pal)

    w = MainWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
