from __future__ import annotations

import json
import os
import sys
from typing import Optional, Callable

from PyQt6.QtCore import Qt, QTimer, QObject, QEvent, pyqtSignal, QLocale
from PyQt6.QtGui import QKeyEvent
from PyQt6.QtWidgets import (
    QMainWindow,
    QTabWidget,
    QFileDialog,
    QInputDialog,
    QMessageBox,
    QLineEdit,
    QDoubleSpinBox,
    QSpinBox,
    QAbstractSpinBox,
    QApplication,
)

# Support BOTH:
#   - project-root:  python3 -m app.main_window
#   - inside app/:   (cd app && python3 ./main_window.py)
try:
    from state import AppState
    from ui.pages.koral_page import KoralPage
    from ui.pages.mcsetup_page import MCSetupPage
    from ui.pages.mcresults_page import MCResultsPage
except ModuleNotFoundError:
    # Running from inside ./app -> add project root to sys.path and retry
    _project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if _project_root not in sys.path:
        sys.path.insert(0, _project_root)

    from state import AppState
    from ui.pages.koral_page import KoralPage
    from ui.pages.mcsetup_page import MCSetupPage
    from ui.pages.mcresults_page import MCResultsPage

try:
    from ui.logging import subscribe as subscribe_logs
    from ui.logging import unsubscribe as unsubscribe_logs
    from ui.logging import log as emit_log
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.logging import subscribe as subscribe_logs  # type: ignore
    from OpenSRIM.ui.logging import unsubscribe as unsubscribe_logs  # type: ignore
    from OpenSRIM.ui.logging import log as emit_log  # type: ignore


class _LogBridge(QObject):
    message = pyqtSignal(str)

    def __init__(self, handler: Callable[[str], None], parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self.message.connect(handler)


class _CommaToDotFilter(QObject):
    """Application-wide event filter that replaces comma key presses with dots
    in ALL widgets (QLineEdit, QSpinBox, QDoubleSpinBox, QTableWidget editors, etc.)."""

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Type.KeyPress:
            try:
                text = event.text()
            except Exception:
                text = ""
            if text == ",":
                # Forward a '.' key press instead of the comma
                dot_event = QKeyEvent(
                    QEvent.Type.KeyPress,
                    Qt.Key.Key_Period,
                    event.modifiers(),
                    ".",
                )
                QApplication.sendEvent(obj, dot_event)
                return True  # eat the comma
        return False  # don't call super — just pass through


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # Ensure C locale is active (dot decimal separator) even if
        # MainWindow is instantiated outside of main().
        QLocale.setDefault(QLocale(QLocale.Language.C))

        self.setWindowTitle("KORAL / MC Simulation")

        self.state = AppState()

        self._create_standard_menu_bar()

        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        self.koral_tab = KoralPage(self.state, on_log=emit_log)
        self.mc_setup_tab = MCSetupPage(
            state=self.state,
            on_log=emit_log,
        )
        self.mc_results_tab = MCResultsPage()

        self.tab_widget.addTab(self.koral_tab, "KORAL")
        self.tab_widget.addTab(self.mc_setup_tab, "MC Setup")
        self.tab_widget.addTab(self.mc_results_tab, "MC Results")

        self.tab_widget.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(self.tab_widget.currentIndex())

        self._apply_startup_geometry()

        # Install global comma-to-dot filter for ALL widgets
        self._comma_filter = _CommaToDotFilter(self)
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self._comma_filter)

        self._log_bridge = _LogBridge(self._on_page_log, parent=self)
        subscribe_logs(self._receive_external_log)

        if not self.state.log_entries:
            emit_log("Ready for simulation.")

    def closeEvent(self, event):
        try:
            unsubscribe_logs(self._receive_external_log)
        except Exception:
            pass
        super().closeEvent(event)

    def _receive_external_log(self, message: str) -> None:
        if hasattr(self, "_log_bridge") and self._log_bridge:
            self._log_bridge.message.emit(str(message))

    # --- sizing: start maximized, but keep a sensible normal size (1366x768) ---
    def _apply_startup_geometry(self):
        screen = self.screen() or (self.windowHandle().screen() if self.windowHandle() else None)
        avail = screen.availableGeometry() if screen else None

        base_w, base_h = 1280, 768
        # Minimum should never exceed available screen size (otherwise cannot fit on small displays)
        if avail:
            min_w = min(base_w, avail.width())
            min_h = min(base_h, avail.height())
        else:
            min_w, min_h = base_w, base_h

        self.setMinimumSize(min_w, min_h)
        self.resize(base_w, base_h)

        # show maximized after the window is created (keeps "normal" size for restore)
        # QTimer.singleShot(0, self.showMaximized)  # <-- remove this line to start in base_w/base_h

    # --- menu bar ---
    def _create_standard_menu_bar(self):
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("File")
        save_act = file_menu.addAction("Save Configuration")
        save_act.triggered.connect(self._handle_save_configuration)
        save_as_act = file_menu.addAction("Save Configuration As…")
        save_as_act.triggered.connect(self._handle_save_configuration_as)
        file_menu.addSeparator()
        load_act = file_menu.addAction("Load Configuration…")
        load_act.triggered.connect(self._handle_load_configuration)

        self.view_menu = menu_bar.addMenu("View")
        save_view_act = self.view_menu.addAction("Save Current Layout…")
        save_view_act.triggered.connect(self._prompt_save_results_layout)
        manage_views_act = self.view_menu.addAction("Manage Views…")
        manage_views_act.triggered.connect(self._manage_results_views)

        self.export_menu = menu_bar.addMenu("Export")
        for label, kind in [("Export Data (CSV)", "csv"), ("Export Plot (PNG)", "png"), ("Export Report (PDF)", "pdf")]:
            action = self.export_menu.addAction(label)
            action.triggered.connect(lambda _, k=kind: self._handle_export(k))

        self.view_menu.menuAction().setVisible(False)
        self.export_menu.menuAction().setVisible(False)

    def _on_tab_changed(self, index: int):
        is_mc_results_tab = self.tab_widget.tabText(index) == "MC Results"
        self.view_menu.menuAction().setVisible(is_mc_results_tab)
        self.export_menu.menuAction().setVisible(is_mc_results_tab)

    # --- logging bridge ---
    def _on_page_log(self, message: str):
        entry = self.state.add_log(message)
        self.mc_setup_tab.update_latest_log(entry)
        try:
            self.koral_tab.update_latest_log(entry)
        except Exception:
            pass

    # --- config save/load ---
    def _handle_save_configuration(self):
        if self.state.current_config_path:
            self._save_configuration_to_path(self.state.current_config_path)
        else:
            self._handle_save_configuration_as()

    def _handle_save_configuration_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Configuration", "", "JSON Files (*.json);;All Files (*)")
        if path:
            self._save_configuration_to_path(path)

    def _save_configuration_to_path(self, path: str):
        try:
            payload = self.mc_setup_tab.collect_simulation_config()
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
            self.state.current_config_path = path
            emit_log(f"Configuration saved to {os.path.basename(path)}.")
        except OSError as exc:
            QMessageBox.warning(self, "Save Configuration", f"Unable to save file:\n{exc}")

    def _handle_load_configuration(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Configuration", "", "JSON Files (*.json);;All Files (*)")
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)
            self.mc_setup_tab.apply_simulation_config(payload)
            self.state.current_config_path = path
            emit_log(f"Configuration loaded from {os.path.basename(path)}.")
        except (OSError, json.JSONDecodeError) as exc:
            QMessageBox.warning(self, "Load Configuration", f"Unable to load file:\n{exc}")

    # --- MC Results menu actions ---
    def _prompt_save_results_layout(self):
        widget = self.mc_results_tab.get_results_widget()
        if not widget:
            return
        name, ok = QInputDialog.getText(self, "Save Layout", "Layout name:")
        if ok and name.strip():
            widget.save_layout(name.strip())
            widget.layout_store.reload()
            emit_log(f"Saved MC Results layout '{name.strip()}'.")

    def _manage_results_views(self):
        widget = self.mc_results_tab.get_results_widget()
        if widget and widget.open_manage_views_dialog(self):
            widget.layout_store.reload()
            emit_log("Updated MC Results view definitions.")

    def _handle_export(self, kind: str):
        label = {"csv": "data", "png": "plot", "pdf": "report"}.get(kind, kind)
        emit_log(f"Export requested for {label}.")
        QMessageBox.information(self, "Export", f"The {label} export will start shortly.")

def main():
    import sys
    from PyQt6.QtWidgets import QApplication
    from PyQt6.QtCore import QLocale

    # Force C locale globally so all QDoubleSpinBox / QSpinBox use dot
    # as decimal separator instead of the system locale (e.g. comma in DE)
    QLocale.setDefault(QLocale(QLocale.Language.C))

    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
