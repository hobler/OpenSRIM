#!/usr/bin/env python3
import sys
from PyQt6 import QtCore, QtGui, QtWidgets


class AccordionItem(QtWidgets.QFrame):
    """
    Extrahierter Accordion-Code (leicht bereinigt), mit Height-Animation.
    Usage: AccordionItem("Title", content_widget, expanded=True/False)
    """
    toggled = QtCore.pyqtSignal(bool)

    def __init__(self, title: str, content: QtWidgets.QWidget, expanded: bool = True, parent=None):
        super().__init__(parent)
        self._expanded = expanded
        self._anim: QtCore.QPropertyAnimation | None = None

        self.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.setObjectName("accordion-item")
        self.setStyleSheet("""
        QFrame#accordion-item {
            border: 1px solid #e5e7eb;
            border-radius: 10px;
            background: #ffffff;
        }
        QPushButton#accordion-header {
            text-align: left;
            padding: 10px 14px;
            border: none;
            font-weight: 600;
            background: #f8fafc;
            border-top-left-radius: 10px;
            border-top-right-radius: 10px;
            border-bottom: 1px solid #e5e7eb;
        }
        QPushButton#accordion-header:hover { background: #f3f4f6; }
        QFrame#accordion-body {
            background: #ffffff;
            border-bottom-left-radius: 10px;
            border-bottom-right-radius: 10px;
        }
        """)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        self.headerBtn = QtWidgets.QPushButton(title)
        self.headerBtn.setObjectName("accordion-header")
        self.headerBtn.setCursor(QtGui.QCursor(QtCore.Qt.CursorShape.PointingHandCursor))
        self.headerBtn.setIcon(self._icon_for_state(expanded))
        self.headerBtn.clicked.connect(self.toggle)
        lay.addWidget(self.headerBtn)

        self.body = QtWidgets.QFrame()
        self.body.setObjectName("accordion-body")
        self.bodyLay = QtWidgets.QVBoxLayout(self.body)
        self.bodyLay.setContentsMargins(12, 12, 12, 12)
        self.bodyLay.setSpacing(8)
        self.bodyLay.addWidget(content)
        lay.addWidget(self.body)

        self._apply_state(expanded, animate=False)

    def _icon_for_state(self, expanded: bool) -> QtGui.QIcon:
        icon = QtGui.QIcon()
        if expanded:
            icon.addPixmap(self.style().standardPixmap(QtWidgets.QStyle.StandardPixmap.SP_ArrowDown))
        else:
            icon.addPixmap(self.style().standardPixmap(QtWidgets.QStyle.StandardPixmap.SP_ArrowRight))
        return icon

    def toggle(self):
        self._apply_state(not self._expanded, animate=True)

    def _apply_state(self, expanded: bool, animate: bool):
        self._expanded = expanded
        self.headerBtn.setIcon(self._icon_for_state(expanded))

        if animate:
            start = self.body.maximumHeight()

            if expanded:
                # kurz "auf" setzen, um die Zielhöhe zu berechnen
                self.body.setMaximumHeight(10**6)
                target = self.body.sizeHint().height()
                self.body.setMaximumHeight(start)
                self.body.setVisible(True)
            else:
                target = 0
                self.body.setVisible(True)

            anim = QtCore.QPropertyAnimation(self.body, b"maximumHeight", self)
            anim.setDuration(160)
            anim.setStartValue(max(start, 0))
            anim.setEndValue(target)
            anim.setEasingCurve(QtCore.QEasingCurve.Type.InOutCubic)
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

    def _on_anim_finished(self, expanded: bool):
        if expanded:
            self.body.setMaximumHeight(16777215)
            self.body.setVisible(True)
        else:
            self.body.setMaximumHeight(0)
            self.body.setVisible(False)


class DemoWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyQt6 Accordion Demo")
        self.resize(720, 600)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)

        root = QtWidgets.QVBoxLayout(central)
        root.setContentsMargins(14, 14, 14, 14)
        root.setSpacing(10)

        # 5 Sample Accordions untereinander
        for i in range(1, 6):
            content = self._make_sample_content(i)
            acc = AccordionItem(f"Accordion {i}", content, expanded=(i == 1))
            root.addWidget(acc)

        root.addStretch(1)

    def _make_sample_content(self, idx: int) -> QtWidgets.QWidget:
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        title = QtWidgets.QLabel(f"Inhalt für Accordion {idx}")
        title.setStyleSheet("font-weight: 600;")

        body = QtWidgets.QLabel(
            "Lorem ipsum dolor sit amet, consectetur adipiscing elit. "
            "Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua."
        )
        body.setWordWrap(True)

        form = QtWidgets.QFormLayout()
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)
        form.addRow("Feld A:", QtWidgets.QLineEdit(f"Value {idx}.A"))
        form.addRow("Feld B:", QtWidgets.QLineEdit(f"Value {idx}.B"))
        form.addRow("Option:", QtWidgets.QComboBox())
        form.itemAt(5).widget().addItems(["One", "Two", "Three"])  # ComboBox aus FormLayout ziehen

        btnrow = QtWidgets.QHBoxLayout()
        btnrow.addStretch(1)
        btnrow.addWidget(QtWidgets.QPushButton("Apply"))
        btnrow.addWidget(QtWidgets.QPushButton("Reset"))

        lay.addWidget(title)
        lay.addWidget(body)
        lay.addLayout(form)
        lay.addLayout(btnrow)
        return w


def main():
    app = QtWidgets.QApplication(sys.argv)
    w = DemoWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
