from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QToolButton, QWidget


class AdvancedSettingsButton(QToolButton):
    """Small gear-icon button that opens a page's Advanced Settings section.

    Used consistently across KORAL / MC Setup / MC Results so every
    "open advanced options" affordance looks and behaves the same way,
    instead of each page building its own ad hoc button.
    """

    def __init__(self, tooltip: str = "", parent: QWidget | None = None):
        super().__init__(parent)
        self.setText("⚙")
        self.setToolTip(tooltip)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedSize(22, 22)
