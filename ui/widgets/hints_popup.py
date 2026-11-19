from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import tomllib  # type: ignore
except ImportError:  # pragma: no cover
    import tomli as tomllib  # type: ignore

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QTabWidget,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
    QMainWindow,
)


@dataclass(frozen=True)
class HintItem:
    id: str
    title: str
    content_md: str


class HintRepository:
    """
    TOML schema (minimal):
    [[pages.KORAL]]
    id = "ion"
    title = "Ion Selection"
    content = '''# ...markdown...'''

    [[pages."MC Setup"]]
    ...
    """
    def __init__(self, toml_path: str | Path):
        self.toml_path = Path(toml_path)

    def load(self) -> Dict[str, List[HintItem]]:
        with open(self.toml_path, "rb") as fh:
            raw = tomllib.load(fh)
        pages = raw.get("pages", {})
        out: Dict[str, List[HintItem]] = {}
        for page_id, items in pages.items():
            if not isinstance(items, list):
                continue
            parsed: List[HintItem] = []
            for it in items:
                if not isinstance(it, dict):
                    continue
                hid = str(it.get("id", "")).strip()
                title = str(it.get("title", "")).strip()
                content = str(it.get("content", "")).rstrip()
                if not hid or not title:
                    continue
                parsed.append(HintItem(id=hid, title=title, content_md=content))
            out[str(page_id)] = parsed
        return out


class HintPopup(QDialog):
    def __init__(self, parent=None, *, repo: HintRepository):
        super().__init__(parent)
        self.setWindowTitle("Hints")
        self.resize(900, 600)

        self._repo = repo
        self._data_by_page: Dict[str, List[HintItem]] = {}
        self._current_page: Optional[str] = None
        self._current_hints: List[HintItem] = []

        root = QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(8)

        header = QHBoxLayout()
        self.page_label = QLabel("Page: -")
        self.page_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        header.addWidget(self.page_label, 1)

        root.addLayout(header)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        self.nav = QListWidget()
        self.nav.setMinimumWidth(240)
        self.nav.currentRowChanged.connect(self._on_nav_changed)

        self.viewer = QTextBrowser()
        self.viewer.setOpenExternalLinks(True)
        self.viewer.setMarkdown("Select a hint on the left.")

        splitter.addWidget(self.nav)
        splitter.addWidget(self.viewer)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        root.addWidget(splitter, 1)

        self.reload()

    def reload(self) -> None:
        self._data_by_page = self._repo.load()
        # re-apply current page filter
        if self._current_page:
            self.set_current_page(self._current_page)
        else:
            self._render_nav([])

    def set_current_page(self, page_id: str) -> None:
        self._current_page = str(page_id)
        self.page_label.setText(f"Page: {self._current_page}")
        hints = self._data_by_page.get(self._current_page, [])
        self._render_nav(hints)

    def show_hint(self, page_id: str, hint_id: str) -> None:
        self.set_current_page(page_id)
        idx = next((i for i, h in enumerate(self._current_hints) if h.id == hint_id), -1)
        if idx >= 0:
            self.nav.setCurrentRow(idx)
        else:
            self.viewer.setMarkdown(f"## Not found\nHint `{hint_id}` not found on page `{page_id}`.")
        self.show()
        self.raise_()
        self.activateWindow()

    def _render_nav(self, hints: List[HintItem]) -> None:
        self._current_hints = list(hints)
        self.nav.blockSignals(True)
        self.nav.clear()
        for h in self._current_hints:
            item = QListWidgetItem(h.title)
            item.setData(Qt.ItemDataRole.UserRole, h.id)
            self.nav.addItem(item)
        self.nav.blockSignals(False)

        if self._current_hints:
            self.nav.setCurrentRow(0)
        else:
            self.viewer.setMarkdown("No hints available for this page.")

    def _on_nav_changed(self, row: int) -> None:
        if not (0 <= row < len(self._current_hints)):
            return
        h = self._current_hints[row]
        self.viewer.setMarkdown(h.content_md)


class HintSystem:
    """Small facade that owns a single popup and returns '?' buttons bound to hints."""
    def __init__(self, *, repo_path: str | Path, parent=None):
        self.repo = HintRepository(repo_path)
        self.popup = HintPopup(parent, repo=self.repo)

    def set_current_page(self, page_id: str) -> None:
        self.popup.set_current_page(page_id)

    def make_hint_button(self, *, page_id: str, hint_id: str, parent=None) -> QPushButton:
        btn = QPushButton("?", parent)
        btn.setFixedSize(22, 22)
        btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        btn.setToolTip("Hint")
        btn.clicked.connect(lambda: self.popup.show_hint(page_id, hint_id))
        return btn


# ---------------- TEST APP ----------------
def _demo_page(title: str, hint_system: HintSystem) -> QWidget:
    w = QWidget()
    layout = QVBoxLayout(w)

    row1 = QHBoxLayout()
    row1.addWidget(QLabel(f"{title}: Ion Selection"))
    row1.addWidget(hint_system.make_hint_button(page_id=title, hint_id="ion"))
    row1.addStretch(1)
    layout.addLayout(row1)

    row2 = QHBoxLayout()
    row2.addWidget(QLabel(f"{title}: Elements Table"))
    row2.addWidget(hint_system.make_hint_button(page_id=title, hint_id="elements"))
    row2.addStretch(1)
    layout.addLayout(row2)

    layout.addStretch(1)
    return w


def main():
    app = QApplication([])

    # assumes this file lives at app/ui/hints/hints_popup.py
    root = Path(__file__).resolve().parents[2]  # .../app
    hints_path = "hints.toml"

    hs = HintSystem(repo_path=hints_path)

    win = QMainWindow()
    win.setWindowTitle("Hint Popup Demo")
    win.resize(1000, 650)

    tabs = QTabWidget()
    win.setCentralWidget(tabs)

    for page_id in ["KORAL", "MC Setup"]:
        tabs.addTab(_demo_page(page_id, hs), page_id)

    def _sync_page(_idx: int):
        hs.set_current_page(tabs.tabText(tabs.currentIndex()))

    tabs.currentChanged.connect(_sync_page)
    _sync_page(0)

    win.show()
    app.exec()


if __name__ == "__main__":
    main()
