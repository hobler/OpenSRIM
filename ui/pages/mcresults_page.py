from __future__ import annotations

from PyQt6.QtWidgets import QWidget, QVBoxLayout

try:
    from ui.pages.simulation.simulation_page import MCResultsWidget
except ModuleNotFoundError:  # pragma: no cover
    from OpenSRIM.ui.pages.simulation.simulation_page import MCResultsWidget  # type: ignore


class MCResultsPage(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        self._results_widget = MCResultsWidget()
        layout.addWidget(self._results_widget)

    def get_results_widget(self):
        return self._results_widget
