from __future__ import annotations

import json
import os
import sys
from typing import Optional, Callable

from PyQt6.QtCore import Qt, QTimer, QObject, QEvent, pyqtSignal, QLocale
from PyQt6.QtGui import QIcon, QKeyEvent
from PyQt6.QtWidgets import (
    QMainWindow,
    QTabWidget,
    QFileDialog,
    QInputDialog,
    QMessageBox,
    QApplication,
    QToolButton,
    QStyle,
)

# Support BOTH:
#   - project-root:  python3 -m app.main_window
#   - inside app/:   (cd app && python3 ./main_window.py)
try:
    from state import AppState
    from ui.pages.koral_page import KoralPage
    from ui.pages.mcsetup_page import MCSetupPage
    from ui.pages.mcresults_page import MCResultsPage
    from ui.pages.advanced_options_page import AdvancedOptionsPage
    from ui.pages.simulation.simulation_page import SinglePlotPage
except ModuleNotFoundError:
    # Running from inside ./app -> add project root to sys.path and retry
    _project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if _project_root not in sys.path:
        sys.path.insert(0, _project_root)

    from state import AppState
    from ui.pages.koral_page import KoralPage
    from ui.pages.mcsetup_page import MCSetupPage
    from ui.pages.mcresults_page import MCResultsPage
    from ui.pages.advanced_options_page import AdvancedOptionsPage
    from ui.pages.simulation.simulation_page import SinglePlotPage

try:
    from ui.logging import subscribe as subscribe_logs
    from ui.logging import unsubscribe as unsubscribe_logs
    from ui.logging import log as emit_log
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.logging import subscribe as subscribe_logs  # type: ignore
    from OpenSRIM.ui.logging import unsubscribe as unsubscribe_logs  # type: ignore
    from OpenSRIM.ui.logging import log as emit_log  # type: ignore


APP_ID = "opensrim"
APP_NAME = "OpenSRIM"
APP_DIR = os.path.dirname(os.path.abspath(__file__))
APP_ICON_PATH = os.path.join(APP_DIR, "design", "assets", "opensrim-icon-app-512.png")
APP_DESKTOP_PATH = os.path.join(
    os.environ.get("XDG_DATA_HOME", os.path.expanduser("~/.local/share")),
    "applications",
    f"{APP_ID}.desktop",
)


def load_app_icon() -> QIcon:
    return QIcon(APP_ICON_PATH) if os.path.exists(APP_ICON_PATH) else QIcon()


def _desktop_quoted(path: str) -> str:
    escaped = path.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def ensure_linux_desktop_entry() -> None:
    if not sys.platform.startswith("linux") or not os.path.exists(APP_ICON_PATH):
        return
    desktop_content = "\n".join(
        [
            "[Desktop Entry]",
            "Type=Application",
            f"Name={APP_NAME}",
            "Comment=Monte-Carlo ion transport simulation",
            f"Exec=python3 {_desktop_quoted(os.path.join(APP_DIR, 'main_window.py'))}",
            f"Icon={APP_ICON_PATH}",
            "Terminal=false",
            "Categories=Science;Education;",
            "StartupNotify=true",
            f"StartupWMClass={APP_NAME}",
            "X-OpenSRIM-Generated=true",
            "",
        ]
    )
    try:
        if os.path.exists(APP_DESKTOP_PATH):
            with open(APP_DESKTOP_PATH, "r", encoding="utf-8") as existing:
                current_content = existing.read()
            if current_content == desktop_content:
                return
            if "X-OpenSRIM-Generated=true" not in current_content:
                return
        os.makedirs(os.path.dirname(APP_DESKTOP_PATH), exist_ok=True)
        with open(APP_DESKTOP_PATH, "w", encoding="utf-8") as desktop_file:
            desktop_file.write(desktop_content)
        os.chmod(APP_DESKTOP_PATH, 0o644)
    except OSError:
        pass


def configure_application(app: QApplication) -> None:
    app.setApplicationName(APP_NAME)
    app.setApplicationDisplayName(APP_NAME)
    app.setDesktopFileName(APP_ID)
    app.setWindowIcon(load_app_icon())


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
                dot_event = QKeyEvent(
                    QEvent.Type.KeyPress,
                    Qt.Key.Key_Period,
                    event.modifiers(),
                    ".",
                )
                QApplication.sendEvent(obj, dot_event)
                return True
        return False


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # Force dot as decimal separator globally (independent of OS locale).
        QLocale.setDefault(QLocale(QLocale.Language.C))

        self.setWindowTitle("OpenSRIM")
        self.setWindowIcon(load_app_icon())

        self.state = AppState()

        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        self.koral_tab = KoralPage(self.state, on_log=emit_log)
        self.mc_setup_tab = MCSetupPage(
            state=self.state,
            on_log=emit_log,
        )
        # Track all MC Results pages (multi-tab support).
        self.mc_results_tabs: list = []
        # Keep the first one accessible for backwards compatibility.
        self.mc_results_tab = self._create_mc_results_tab()
        self.single_plot_page = SinglePlotPage()
        self.advanced_options_tab = AdvancedOptionsPage()
        self.advanced_options_tab.set_project_elements_provider(self.mc_setup_tab.get_project_element_symbols)
        self._advanced_previous_widget = None

        self.tab_widget.addTab(self.koral_tab, "KORAL")
        self.tab_widget.addTab(self.mc_setup_tab, "MC Setup")
        # Initial MC Results tab; more can be added via the "+" button or
        # automatically when a new simulation completes.
        self._mc_results_first_index = self.tab_widget.addTab(self.mc_results_tab, "MC Results")
        self.tab_widget.addTab(self.single_plot_page, "Single Plot")

        # "+" button on the tab bar to add additional MC Results tabs.
        self._add_tab_btn = QToolButton(self.tab_widget)
        self._add_tab_btn.setText("+")
        self._add_tab_btn.setToolTip("Open another MC Results tab")
        self._add_tab_btn.setAutoRaise(True)
        self._add_tab_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._add_tab_btn.clicked.connect(self._on_add_results_tab_clicked)
        self.tab_widget.setCornerWidget(self._add_tab_btn, Qt.Corner.TopRightCorner)

        # Allow closing additional MC Results tabs (not the first one).
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self._on_tab_close_requested)
        self._update_tab_close_buttons()

        self.mc_setup_tab.advanced_requested.connect(self._open_advanced_options)
        self.mc_results_tab.advanced_requested.connect(self._open_advanced_options)
        self.koral_tab.advanced_requested.connect(self._open_advanced_options)
        self.mc_setup_tab.save_requested.connect(self._handle_save_configuration)
        self.mc_setup_tab.load_requested.connect(self._handle_load_configuration)
        self.mc_setup_tab.simulation_finished.connect(self._on_simulation_finished)
        self.mc_setup_tab.results_update.connect(self._on_results_update)
        self.mc_setup_tab.advanced_simulation_settings_changed.connect(
            self.advanced_options_tab.apply_mc_setup_advanced_config
        )
        self.advanced_options_tab.mc_follow_recoils_changed.connect(self.mc_setup_tab.set_follow_recoils)
        self.advanced_options_tab.mc_replacement_collisions_changed.connect(self.mc_setup_tab.set_replacement_collisions)
        self.advanced_options_tab.mc_pmax_min_changed.connect(self.mc_setup_tab.set_pmax_min)
        self.advanced_options_tab.mc_pmax_max_changed.connect(self.mc_setup_tab.set_pmax_max)
        self.advanced_options_tab.mc_psi_min_changed.connect(self.mc_setup_tab.set_psi_min)
        self.advanced_options_tab.mc_de_min_changed.connect(self.mc_setup_tab.set_de_min)
        self.advanced_options_tab.mc_psi_min_surface_changed.connect(self.mc_setup_tab.set_psi_min_surface)
        self.advanced_options_tab.mc_de_min_surface_changed.connect(self.mc_setup_tab.set_de_min_surface)
        self.advanced_options_tab.mc_nthreads_changed.connect(self.mc_setup_tab.set_nthreads)
        self.advanced_options_tab.mc_rng_seed_changed.connect(self.mc_setup_tab.set_rng_seed)
        self.advanced_options_tab.mc_electronic_stopping_changed.connect(self.mc_setup_tab.set_electronic_stopping_model)
        self.advanced_options_tab.mc_scattering_algorithm_changed.connect(self.mc_setup_tab.set_scattering_algorithm)
        self.advanced_options_tab.mc_n_absc_changed.connect(self.mc_setup_tab.set_n_absc)
        self.advanced_options_tab.mc_lindhard_correction_changed.connect(self.mc_setup_tab.set_lindhard_correction)
        self.advanced_options_tab.atoms_columns_visibility_changed.connect(
            lambda disp, latt, surf: self.mc_setup_tab._set_atoms_energy_columns_visible(
                show_disp=disp, show_latt=latt, show_surf=surf
            )
        )
        self.advanced_options_tab.mc_ion_angle_changed.connect(self._apply_mc_ion_angle)

        # Display settings (Advanced Options → MC Results plot area).
        # Connected for the first tab; new tabs are wired in
        # _create_mc_results_tab().
        self._wire_display_settings_for_results_page(self.mc_results_tab)
        # Single plot tab also reacts to font-size changes.
        self.advanced_options_tab.plot_font_size_changed.connect(
            self.single_plot_page.set_font_size
        )
        self.advanced_options_tab.beam_from_top_changed.connect(
            self.single_plot_page.set_beam_from_top
        )
        self.advanced_options_tab.koral_solver_changed.connect(
            self.koral_tab.set_koral_solver_settings
        )
        self.advanced_options_tab.histogram_settings_changed.connect(
            self.mc_setup_tab.set_histogram_settings
        )
        self.mc_setup_tab.histogram_defaults_changed.connect(
            self.advanced_options_tab.apply_histogram_defaults
        )

        # initial sync
        try:
            self.advanced_options_tab.set_mc_ion_angle(self.mc_setup_tab.get_ion_angle())
        except Exception:
            pass
        try:
            self.advanced_options_tab.apply_mc_setup_advanced_config(
                self.mc_setup_tab.get_advanced_simulation_settings()
            )
        except Exception:
            pass
        try:
            self.advanced_options_tab.apply_histogram_defaults(
                self.mc_setup_tab.get_histogram_default_settings()
            )
        except Exception:
            pass
        try:
            persisted_display = self.state.get_persisted_display_settings()
            if persisted_display:
                self.advanced_options_tab.apply_config({"display_settings": persisted_display})
        except Exception:
            pass

        self.advanced_options_tab.toolbar_visibility_changed.connect(self._persist_display_settings)
        self.advanced_options_tab.borders_visibility_changed.connect(self._persist_display_settings)
        self.advanced_options_tab.beam_from_top_changed.connect(self._persist_display_settings)
        self.advanced_options_tab.plot_font_size_changed.connect(self._persist_display_settings)
        self.advanced_options_tab.columns_changed.connect(self._persist_display_settings)

        self.tab_widget.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(self.tab_widget.currentIndex())
        self.tab_widget.setCurrentWidget(self.mc_setup_tab)
        self._on_tab_changed(self.tab_widget.currentIndex())

        self._apply_startup_geometry()

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

        base_w, base_h = 1366, 768
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

    def _on_tab_changed(self, _index: int):
        pass

    def _create_mc_results_tab(self) -> MCResultsPage:
        page = MCResultsPage()
        # Forward double-click on a tile → load into Single Plot tab and switch.
        page.plot_open_in_single.connect(self._on_plot_open_in_single)
        # Each new tab also gets the advanced-options wiring.
        page.advanced_requested.connect(self._open_advanced_options)
        self.mc_results_tabs.append(page)
        # Wire display-settings only if advanced_options_tab already exists
        # (first call happens before it is created — that path handles its
        # wiring inline in __init__).
        if hasattr(self, "advanced_options_tab"):
            self._wire_display_settings_for_results_page(page)
        return page

    def _wire_display_settings_for_results_page(self, page: MCResultsPage) -> None:
        widget = page.get_results_widget()
        self.advanced_options_tab.toolbar_visibility_changed.connect(
            widget.set_plot_toolbar_visible
        )
        self.advanced_options_tab.columns_changed.connect(
            widget.set_plot_columns
        )
        self.advanced_options_tab.borders_visibility_changed.connect(
            widget.set_plot_borders_visible
        )
        self.advanced_options_tab.plot_font_size_changed.connect(
            widget.set_plot_font_size
        )
        self.advanced_options_tab.beam_from_top_changed.connect(
            widget.set_beam_from_top
        )

    def _is_mc_results_tab_index(self, index: int) -> bool:
        widget = self.tab_widget.widget(index)
        return isinstance(widget, MCResultsPage)

    def _is_advanced_tab(self, widget) -> bool:
        return isinstance(widget, AdvancedOptionsPage)

    def _update_tab_close_buttons(self) -> None:
        """Show close buttons on MC Results tabs (as long as one stays open) and Advanced tabs."""
        bar = self.tab_widget.tabBar()
        if bar is None:
            return
        total_results = sum(
            1 for i in range(self.tab_widget.count())
            if isinstance(self.tab_widget.widget(i), MCResultsPage)
        )
        for i in range(self.tab_widget.count()):
            widget = self.tab_widget.widget(i)
            closable = False
            if isinstance(widget, MCResultsPage):
                # Every results tab is closable as long as more than one exists,
                # so at least one results tab always remains open.
                closable = total_results > 1
            elif self._is_advanced_tab(widget):
                closable = True
            try:
                btn_right = bar.tabButton(i, bar.ButtonPosition.RightSide)
                btn_left = bar.tabButton(i, bar.ButtonPosition.LeftSide)
                if not closable:
                    if btn_right is not None:
                        btn_right.hide()
                    if btn_left is not None:
                        btn_left.hide()
                else:
                    if btn_right is not None:
                        btn_right.show()
                    if btn_left is not None:
                        btn_left.show()
            except Exception:
                pass

    def _on_add_results_tab_clicked(self) -> None:
        page = self._create_mc_results_tab()
        # Insert just after the last MC Results tab.
        insert_at = self.tab_widget.count()
        for i in range(self.tab_widget.count() - 1, -1, -1):
            if isinstance(self.tab_widget.widget(i), MCResultsPage):
                insert_at = i + 1
                break
        title = f"MC Results {len(self.mc_results_tabs)}"
        self.tab_widget.insertTab(insert_at, page, title)
        self.tab_widget.setCurrentWidget(page)
        self._update_tab_close_buttons()

    def _on_tab_close_requested(self, index: int) -> None:
        widget = self.tab_widget.widget(index)
        if self._is_advanced_tab(widget):
            self._close_advanced_options_tab(index)
            return
        if not isinstance(widget, MCResultsPage):
            return
        # Keep at least one MC Results tab open at all times.
        total_results = sum(
            1 for i in range(self.tab_widget.count())
            if isinstance(self.tab_widget.widget(i), MCResultsPage)
        )
        if total_results <= 1:
            return
        try:
            self.mc_results_tabs.remove(widget)
        except ValueError:
            pass
        if getattr(self, "_last_live_target", None) is widget:
            self._last_live_target = None
        self.tab_widget.removeTab(index)
        widget.deleteLater()
        # Keep the backwards-compat reference pointing at a live results tab
        # (at least one always remains, so the list is never empty here).
        if self.mc_results_tab is widget and self.mc_results_tabs:
            self.mc_results_tab = self.mc_results_tabs[0]
        self._update_tab_close_buttons()

    def _on_plot_open_in_single(self, plot_id: str, plot_info: object) -> None:
        try:
            self.single_plot_page.set_plot(str(plot_id), dict(plot_info) if isinstance(plot_info, dict) else {})
        except Exception:
            return
        self.tab_widget.setCurrentWidget(self.single_plot_page)

    def _open_advanced_options(self, section_id: str):
        current = self.tab_widget.currentWidget()
        if current is not self.advanced_options_tab:
            self._advanced_previous_widget = current
        # The page the Advanced Settings tab should sit next to.
        anchor = current if current is not self.advanced_options_tab else self._advanced_previous_widget

        # Always (re)position the Advanced Settings tab directly to the right of
        # the tab the user came from: remove it from its current spot (if open)
        # and re-insert next to the anchor page. This way clicking a gear on a
        # different page moves the tab next to that page instead of leaving it
        # where it was first opened.
        adv_index = self.tab_widget.indexOf(self.advanced_options_tab)
        if adv_index >= 0:
            self.tab_widget.removeTab(adv_index)
        anchor_index = self.tab_widget.indexOf(anchor) if anchor is not None else -1
        if anchor_index < 0:
            anchor_index = self.tab_widget.count() - 1
        self.tab_widget.insertTab(
            anchor_index + 1, self.advanced_options_tab, "Advanced Settings"
        )
        self._update_tab_close_buttons()

        # Show only the sections relevant to the originating page and label the
        # tab accordingly (MC Setup vs MC Results vs KORAL).
        context = self.advanced_options_tab.context_for_section(section_id)
        self.advanced_options_tab.set_visible_context(context)
        ctx_label = {
            "mc_setup": "MC Setup (Adv.)",
            "mc_results": "MC Results (Adv.)",
            "koral": "KORAL (Adv.)",
        }.get(context or "", "Advanced Settings")
        self.tab_widget.setTabText(
            self.tab_widget.indexOf(self.advanced_options_tab), ctx_label
        )

        self.tab_widget.setCurrentWidget(self.advanced_options_tab)
        # Sync the correct angle when opening the relevant section.
        if section_id == "ion_selection_mc":
            try:
                self.advanced_options_tab.set_mc_ion_angle(self.mc_setup_tab.get_ion_angle())
            except Exception:
                pass
        if section_id == "histogram_settings":
            try:
                self.advanced_options_tab.apply_histogram_defaults(
                    self.mc_setup_tab.get_histogram_default_settings()
                )
            except Exception:
                pass
        self.advanced_options_tab.open_section(section_id)

    def _close_advanced_options_tab(self, index: int) -> None:
        previous = self._advanced_previous_widget
        fallback = self.mc_setup_tab
        if previous is self.advanced_options_tab:
            previous = None
        if previous is not None and self.tab_widget.indexOf(previous) >= 0:
            self.tab_widget.setCurrentWidget(previous)
        else:
            self.tab_widget.setCurrentWidget(fallback)
        self.tab_widget.removeTab(index)
        self._advanced_previous_widget = None
        self._update_tab_close_buttons()

    def _apply_mc_ion_angle(self, angle: float) -> None:
        try:
            self.mc_setup_tab.set_ion_angle(angle)
        except Exception:
            pass

    # --- logging bridge ---
    def _on_page_log(self, message: str):
        entry = self.state.add_log(message)
        self.mc_setup_tab.update_latest_log(entry)
        try:
            self.koral_tab.update_latest_log(entry)
        except Exception:
            pass

    def _on_simulation_finished(self, results_dir: str) -> None:
        """Load results into the run's MC Results tab and switch to it.

        Live updates during the run already populated a results tab (tracked in
        ``_last_live_target``). Reuse that same tab on completion instead of
        spawning a new one — otherwise the first simulation would fill
        "MC Results" via live updates and then create a redundant
        "MC Results 2" here. A fresh tab is only created if no live update ran
        (e.g. update interval ≥ ion count). ``_last_live_target`` is reset so
        the *next* simulation starts in a new tab.
        """
        target = getattr(self, "_last_live_target", None)
        if target is None or target not in self.mc_results_tabs:
            target = self._get_or_create_results_tab_for_new_run()
        target.load_results_from_directory(results_dir)
        self.tab_widget.setCurrentWidget(target)
        emit_log(f"Results loaded from {results_dir}")
        # Next run should target a fresh tab.
        self._last_live_target = None

    def _on_results_update(self, results_dir: str) -> None:
        """Live-refresh the active MC Results page during simulation."""
        target = getattr(self, "_last_live_target", None)
        if target is None or target not in self.mc_results_tabs:
            target = self._get_or_create_results_tab_for_new_run()
            self._last_live_target = target
        target.load_results_from_directory(results_dir, silent=True)

    def _get_or_create_results_tab_for_new_run(self) -> MCResultsPage:
        """Pick the first 'unused' MC Results tab, or create a new one."""
        for page in self.mc_results_tabs:
            # An "unused" tab has no loaded results yet (empty plot list).
            widget = page.get_results_widget()
            try:
                used = widget.plot_area._tiles if hasattr(widget, "plot_area") else None
                if used is None or len(used) == 0:
                    return page
            except Exception:
                continue
        # All existing tabs have content — spawn a new one.
        self._on_add_results_tab_clicked()
        return self.mc_results_tabs[-1]

    def _persist_display_settings(self, *_args) -> None:
        try:
            payload = self.advanced_options_tab.collect_config()
        except Exception:
            return
        display = payload.get("display_settings")
        if isinstance(display, dict):
            self.state.set_persisted_display_settings(display)

    # --- config save/load ---
    def _handle_save_configuration(self):
        if self.state.current_config_path:
            self._save_configuration_to_path(self.state.current_config_path)
        else:
            self._handle_save_configuration_as()

    def _handle_save_configuration_as(self):
        start_path = self.state.current_config_path or self.state.get_dialog_start_directory()
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Configuration",
            start_path,
            "Config Files (*.config);;All Files (*)",
        )
        if path:
            if not path.lower().endswith(".config"):
                path = path + ".config"
            self._save_configuration_to_path(path)

    def _save_configuration_to_path(self, path: str):
        try:
            payload = {
                "format": "OpenSRIM.config",
                "version": 1,
                "mc_setup": self.mc_setup_tab.collect_simulation_config(),
                "koral": self.koral_tab.collect_config() if hasattr(self.koral_tab, "collect_config") else {},
                "advanced": self.advanced_options_tab.collect_config() if hasattr(self.advanced_options_tab, "collect_config") else {},
            }
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
            self.state.current_config_path = path
            self.state.remember_dialog_path(path)
            emit_log(f"Configuration saved to {os.path.basename(path)}.")
        except OSError as exc:
            QMessageBox.warning(self, "Save Configuration", f"Unable to save file:\n{exc}")

    def _handle_load_configuration(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Configuration",
            self.state.get_dialog_start_directory(),
            "Config Files (*.config);;JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                payload = json.load(fh)

            # Backwards compat: older configs were raw MC setup payloads.
            if isinstance(payload, dict) and "mc_setup" in payload:
                mc_payload = payload.get("mc_setup") or {}
                koral_payload = payload.get("koral") or {}
                adv_payload = payload.get("advanced") or {}
            else:
                mc_payload = payload
                koral_payload = {}
                adv_payload = {}

            self.mc_setup_tab.apply_simulation_config(mc_payload)

            if koral_payload and hasattr(self.koral_tab, "apply_config"):
                self.koral_tab.apply_config(koral_payload)

            if adv_payload and hasattr(self.advanced_options_tab, "apply_config"):
                self.advanced_options_tab.apply_config(adv_payload)

            # Ensure Advanced Options shows the angles actually applied.
            try:
                self.advanced_options_tab.set_mc_ion_angle(self.mc_setup_tab.get_ion_angle())
            except Exception:
                pass
            self.state.current_config_path = path
            self.state.remember_dialog_path(path)
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

    import signal

    QLocale.setDefault(QLocale(QLocale.Language.C))
    if os.environ.get("OPENSRIM_CREATE_DESKTOP_ENTRY") == "1":
        ensure_linux_desktop_entry()

    app = QApplication(sys.argv)
    configure_application(app)
    win = MainWindow()
    win.show()

    # Allow Ctrl+C from the terminal to quit the application.
    # Qt blocks Python's signal handling while the event loop runs; a periodic
    # timer gives Python a chance to process SIGINT.
    signal.signal(signal.SIGINT, lambda *_: app.quit())
    sigint_timer = QTimer()
    sigint_timer.setInterval(200)
    sigint_timer.timeout.connect(lambda: None)
    sigint_timer.start()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
