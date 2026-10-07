from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt


class ToggleSwitch(QtWidgets.QWidget):
    """Ein einfacher, klickbarer Toggle-Switch ohne Animation."""
    toggled = QtCore.pyqtSignal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._checked = False
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    def sizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(50, 26)

    def minimumSizeHint(self) -> QtCore.QSize:
        return self.sizeHint()

    def isChecked(self) -> bool:
        return self._checked

    def setChecked(self, checked: bool):
        checked = bool(checked)
        if self._checked == checked:
            return
        self._checked = checked
        self.toggled.emit(self._checked)
        self.update()

    def toggle(self):
        if self.isEnabled():
            self.setChecked(not self._checked)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.toggle()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event):
        if event.key() in (Qt.Key.Key_Space, Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self.toggle()
            event.accept()
            return
        super().keyPressEvent(event)

    def paintEvent(self, _event):
        radius = self.height() / 2
        knob_radius = radius - 3
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)

        bg_color = QtGui.QColor(0, 150, 136) if self._checked else QtGui.QColor(160, 160, 160)
        p.setBrush(bg_color)
        p.setPen(QtCore.Qt.PenStyle.NoPen)

        rect = QtCore.QRect(0, 0, self.width(), self.height())
        p.drawRoundedRect(rect, radius, radius)

        x_left = 3
        x_right = self.width() - (knob_radius * 2) - 3
        x = x_right if self._checked else x_left
        knob_rect = QtCore.QRectF(x, 3, knob_radius * 2, knob_radius * 2)

        p.setBrush(QtGui.QColor(255, 255, 255))
        p.drawEllipse(knob_rect)
