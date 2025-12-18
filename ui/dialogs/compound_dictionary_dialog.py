from __future__ import annotations

import html
import json
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional, Union

from PyQt6.QtCore import Qt, pyqtSignal, QSize, QRectF, QPointF
from PyQt6.QtGui import QIcon, QPainter, QPixmap, QColor, QPen, QBrush
from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QListWidget,
    QListWidgetItem,
    QScrollArea,
    QLabel,
    QWidget,
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QDoubleSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QMessageBox,
    QToolButton,
    QApplication,
    QStyle,
    QComboBox,
    QInputDialog,
)

from app.ui.widgets.periodic_table_picker import PeriodicTableDialog


class _CompoundEditDialog(QDialog):
    """
    Small editor for one compound. Kept intentionally simple:
    - name_display
    - density_g_cm3
    - section
    - composition rows: Z, fraction
    """

    def __init__(self, parent=None, initial: Optional[dict] = None, sections: Optional[list[str]] = None):
        super().__init__(parent)
        self.setWindowTitle("Edit Compound" if initial else "Add Compound")
        self.setMinimumSize(600, 500)
        self._initial = dict(initial) if isinstance(initial, dict) else {}

        root = QVBoxLayout(self)

        # name_display
        root.addWidget(QLabel("Name"))
        self.ed_name = QLineEdit(self._initial.get("name_display") or self._initial.get("name") or "")
        root.addWidget(self.ed_name)

        # density
        root.addWidget(QLabel("Density (g/cm³)"))
        self.sp_density = QDoubleSpinBox()
        self.sp_density.setDecimals(6)
        self.sp_density.setRange(0.0, 10_000.0)
        dens = self._initial.get("density_g_cm3")
        try:
            self.sp_density.setValue(float(dens) if dens is not None else 0.0)
        except (TypeError, ValueError):
            self.sp_density.setValue(0.0)
        root.addWidget(self.sp_density)

        # section (dropdown + add category)
        root.addWidget(QLabel("Section"))
        sec_row = QHBoxLayout()

        self.cb_section = QComboBox()
        self.cb_section.setEditable(True)
        existing_sections = [s for s in (sections or []) if s]
        for s in sorted(set(existing_sections), key=str.lower):
            self.cb_section.addItem(s)

        initial_section = (self._initial.get("section") or "").strip()
        if initial_section and self.cb_section.findText(initial_section) < 0:
            self.cb_section.addItem(initial_section)
        self.cb_section.setCurrentText(initial_section)

        self.btn_add_section = QToolButton()
        # folder + plus (theme first, then Qt's "new folder" icon)
        self.btn_add_section.setIcon(
            self._theme_icon_multi(
                ["folder-new", "folder-add", "folder-new-symbolic"],
                QStyle.StandardPixmap.SP_FileDialogNewFolder,
            )
        )
        self.btn_add_section.setToolTip("Add category")
        self.btn_add_section.setAutoRaise(True)

        sec_row.addWidget(self.cb_section, 1)
        sec_row.addWidget(self.btn_add_section)
        root.addLayout(sec_row)

        self.btn_add_section.clicked.connect(self._add_section)

        # composition table
        root.addWidget(QLabel("Composition (Element / fraction)"))
        self.tbl = QTableWidget(0, 2)
        self.tbl.setHorizontalHeaderLabels(["Element", "fraction"])
        self.tbl.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl.verticalHeader().setVisible(False)
        root.addWidget(self.tbl)

        for part in (self._initial.get("composition") or []):
            if not isinstance(part, dict):
                continue
            z = part.get("Z")
            fr = part.get("fraction")
            self._append_row(z, fr)

        # row controls: only plus and trash
        row_btns = QHBoxLayout()

        # add row: plus icon
        self.btn_add_row = QToolButton()
        self.btn_add_row.setIcon(
            self._theme_icon_multi(
                ["list-add", "add", "plus", "document-new"],
                QStyle.StandardPixmap.SP_FileDialogNewFolder,
            )
        )
        self.btn_add_row.setToolTip("Add row")
        self.btn_add_row.setAutoRaise(True)

        # delete row: trash
        self.btn_del_row = QToolButton()
        self.btn_del_row.setIcon(
            self._theme_icon_multi(
                ["user-trash", "edit-delete", "trash"],
                QStyle.StandardPixmap.SP_TrashIcon,
            )
        )
        self.btn_del_row.setToolTip("Delete row")
        self.btn_del_row.setAutoRaise(True)

        row_btns.addWidget(self.btn_add_row)
        row_btns.addWidget(self.btn_del_row)
        row_btns.addStretch(1)
        root.addLayout(row_btns)

        # Button-Verhalten: Plus fügt Zeile hinzu und öffnet sofort Periodentabelle
        self.btn_add_row.clicked.connect(self._add_row_and_pick_element)
        self.btn_del_row.clicked.connect(self._delete_selected_rows)

        # ok/cancel
        actions = QHBoxLayout()
        actions.addStretch(1)
        self.btn_cancel = QPushButton("Cancel")
        self.btn_ok = QPushButton("OK")
        self.btn_ok.setDefault(True)
        actions.addWidget(self.btn_cancel)
        actions.addWidget(self.btn_ok)
        root.addLayout(actions)

        self.btn_cancel.clicked.connect(self.reject)
        self.btn_ok.clicked.connect(self.accept)

        self._z_density = self._load_element_densities()
        self._suppress_tbl_signals = False

        # composition table
        self.tbl.cellChanged.connect(self._on_tbl_changed)
        self.tbl.cellDoubleClicked.connect(self._on_tbl_cell_double_clicked)
        self._recompute_density()

    def _theme_icon(self, theme_name: str, fallback: QStyle.StandardPixmap) -> QIcon:
        ico = QIcon.fromTheme(theme_name)
        if not ico.isNull():
            return ico
        return self.style().standardIcon(fallback)

    def _theme_icon_multi(self, theme_names: list[str], fallback: QStyle.StandardPixmap) -> QIcon:
        for name in theme_names:
            ico = QIcon.fromTheme(name)
            if not ico.isNull():
                return ico
        return self.style().standardIcon(fallback)

    def _add_section(self):
        text, ok = QInputDialog.getText(self, "Add category", "Category name:")
        if not ok:
            return
        s = (text or "").strip()
        if not s:
            return
        if self.cb_section.findText(s) < 0:
            self.cb_section.addItem(s)
        self.cb_section.setCurrentText(s)

    def _load_element_densities(self) -> dict[int, float]:
        """
        Reads periodicTableJson from app/ui/widgets and returns {Z: density_g_cm3}.
        """
        widgets_dir = Path(__file__).resolve().parents[1] / "widgets"
        candidates = [
            widgets_dir / "periodicTableJson.json",
            widgets_dir / "periodic_table_json.json",
            widgets_dir / "periodicTable.json",
        ]
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            return {}

        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}

        # Support either {"elements":[...]} or just [...]
        elements = data.get("elements") if isinstance(data, dict) else data
        if not isinstance(elements, list):
            return {}

        out: dict[int, float] = {}
        for el in elements:
            if not isinstance(el, dict):
                continue
            try:
                z = int(el.get("number"))
            except (TypeError, ValueError):
                continue
            dens = el.get("density")
            if dens is None:
                continue
            try:
                out[z] = float(dens)
            except (TypeError, ValueError):
                continue
        return out

    def _on_tbl_changed(self, _row: int, _col: int):
        if self._suppress_tbl_signals:
            return
        self._recompute_density()

    def _recompute_density(self):
        # Weighted average using entered fractions, normalized by sum(fraction).
        total = 0.0
        weighted = 0.0

        for r in range(self.tbl.rowCount()):
            z_item = self.tbl.item(r, 0)
            f_item = self.tbl.item(r, 1)
            if not z_item or not f_item:
                continue

            z = z_item.data(Qt.ItemDataRole.UserRole)
            try:
                z = int(z) if z is not None else None
            except (TypeError, ValueError):
                z = None
            if z is None:
                continue

            try:
                fr = float((f_item.text() or "").strip())
            except ValueError:
                continue

            dens = self._z_density.get(z)
            if dens is None:
                continue

            total += fr
            weighted += fr * dens

        if total <= 0.0:
            return

        self._suppress_tbl_signals = True
        try:
            self.sp_density.setValue(weighted / total)
        finally:
            self._suppress_tbl_signals = False

    def _append_row(self, z: Any, fraction: Any):
        r = self.tbl.rowCount()
        self.tbl.insertRow(r)

        elem_item = QTableWidgetItem("")
        if z is not None:
            try:
                z_int = int(z)
                elem_item.setData(Qt.ItemDataRole.UserRole, z_int)
                elem_item.setText(f"Z={z_int}")
            except (TypeError, ValueError):
                elem_item.setText(str(z))

        self.tbl.setItem(r, 0, elem_item)
        self.tbl.setItem(r, 1, QTableWidgetItem("" if fraction is None else str(fraction)))

    def _add_row_and_pick_element(self):
        self._append_row(None, "")
        row = self.tbl.rowCount() - 1
        self.tbl.setCurrentCell(row, 0)
        self._pick_element_for_row(row)

    def _pick_element_for_row(self, row=None):
        if row is None:
            row = self.tbl.currentRow()
            if row < 0:
                return
        dlg = PeriodicTableDialog(self, compact=True, show_hover_info=True, bordered=True)

        def _apply(element: dict):
            try:
                z = int(element.get("number"))
            except (TypeError, ValueError):
                return
            symbol = element.get("symbol", "")
            name = element.get("name", "")
            txt = f"{symbol} (Z={z})" if symbol else f"Z={z}"
            if name:
                txt = f"{txt} — {name}"

            it = self.tbl.item(row, 0)
            if it is None:
                it = QTableWidgetItem()
                self.tbl.setItem(row, 0, it)  # only insert if not present (prevents ownership errors)

            it.setData(Qt.ItemDataRole.UserRole, z)
            it.setText(txt)
            self._recompute_density()

        dlg.element_selected.connect(_apply)
        dlg.exec()

    def result_compound(self) -> dict:
        name_display = self.ed_name.text().strip()
        section = self.cb_section.currentText().strip()
        density = float(self.sp_density.value())

        composition: list[dict] = []
        for r in range(self.tbl.rowCount()):
            z_item = self.tbl.item(r, 0)
            f_item = self.tbl.item(r, 1)
            if not z_item or not f_item:
                continue
            try:
                fr = float(f_item.text().strip())
            except ValueError:
                continue

            z = z_item.data(Qt.ItemDataRole.UserRole)
            if z is None:
                # fallback: try to parse "Z=.." from text
                txt = (z_item.text() or "").strip()
                if "Z=" in txt:
                    try:
                        z = int(txt.split("Z=", 1)[1].split(")", 1)[0].strip())
                    except ValueError:
                        z = None
            try:
                z = int(z) if z is not None else None
            except (TypeError, ValueError):
                z = None

            if z is None:
                continue
            composition.append({"Z": z, "fraction": fr})

        base = dict(self._initial)  # preserve unknown keys
        base["name_display"] = name_display or base.get("name_display") or base.get("name") or "Unnamed"
        base["name"] = base.get("name") or base["name_display"]
        base["density_g_cm3"] = density
        base["section"] = section
        base["composition"] = composition
        base["n_components"] = len(composition)
        return base

    def _delete_selected_rows(self):
        rows = sorted({i.row() for i in self.tbl.selectedIndexes()}, reverse=True)
        for r in rows:
            self.tbl.removeRow(r)
        self._recompute_density()

    def _on_tbl_cell_double_clicked(self, row: int, col: int):
        if col == 0:
            self._pick_element_for_row(row)

    def accept(self):
        # Prüfe, ob mindestens eine gültige Fraction eingegeben wurde
        has_fraction = False
        for r in range(self.tbl.rowCount()):
            f_item = self.tbl.item(r, 1)
            if f_item:
                try:
                    fr = float(f_item.text().strip())
                    if fr != 0.0:
                        has_fraction = True
                        break
                except Exception:
                    continue
        if not has_fraction:
            QMessageBox.warning(self, "Error", "You have to set fraction!")
            return
        super().accept()


class CompoundDictionaryDialog(QDialog):
    compound_selected = pyqtSignal(dict)

    def __init__(self, parent=None, editable: bool = True):
        super().__init__(parent)
        self._editable = editable

        self._compounds_path = Path(__file__).with_name("compounds.json")
        self._default_path = Path(__file__).with_name("compounds.default.json")

        self._ensure_default_backup()

        self.setWindowTitle("Compound Dictionary")
        self.resize(720, 520)
        self.setModal(True)

        self.compounds: list[dict] = self._load_compounds(self._compounds_path)
        self._collect_all_sections()
        self.current_index = None

        layout = QVBoxLayout(self)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.section_tree = QTreeWidget()
        self.section_tree.setHeaderHidden(True)
        self.section_tree.setSelectionMode(QTreeWidget.SelectionMode.SingleSelection)
        self.tabs.addTab(self.section_tree, "By Section")

        self.alpha_list = QListWidget()
        self.alpha_list.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.tabs.addTab(self.alpha_list, "Alphabetical")

        info_scroll = QScrollArea()
        info_scroll.setWidgetResizable(True)
        layout.addWidget(info_scroll)

        info_container = QWidget()
        info_layout = QVBoxLayout(info_container)
        info_layout.setContentsMargins(8, 8, 8, 8)
        info_scroll.setWidget(info_container)

        self.info_label = QLabel("Select a compound to view its details.")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self.info_label.setWordWrap(True)
        info_layout.addWidget(self.info_label)
        info_layout.addStretch()

        buttons = QHBoxLayout()
        buttons.addStretch()
        self.add_btn = QPushButton("Add to Current Layer")
        self.add_btn.setEnabled(False)
        buttons.addWidget(self.add_btn)
        close_btn = QPushButton("Close")
        buttons.addWidget(close_btn)
        layout.addLayout(buttons)

        self._populate_section_tree()
        self._populate_alpha_list()

        self.section_tree.itemSelectionChanged.connect(self._handle_section_selection)
        self.section_tree.itemDoubleClicked.connect(self._handle_section_double_click)
        self.alpha_list.currentRowChanged.connect(self._handle_alpha_selection)
        self.alpha_list.itemDoubleClicked.connect(self._handle_alpha_double_click)
        self.add_btn.clicked.connect(self._emit_selection)
        close_btn.clicked.connect(self.reject)

        # Append editing toolbar (keeps base layout intact)
        self._install_edit_toolbar()

    def _ensure_default_backup(self):
        # Create a copy of the shipped state once, so user can restore later.
        try:
            if self._default_path.exists():
                return
            if self._compounds_path.exists():
                shutil.copyfile(self._compounds_path, self._default_path)
            else:
                # create both as empty to avoid crashes
                self._default_path.write_text("[]\n", encoding="utf-8")
                self._compounds_path.write_text("[]\n", encoding="utf-8")
        except Exception:
            # don't hard-crash the UI
            pass

    def _load_compounds(self, path: Path) -> list[dict]:
        def _strip_json_comments_and_fix(s: str) -> str:
            # remove // comments
            s = re.sub(r"(?m)^\s*//.*?$", "", s)
            # remove /* ... */ comments
            s = re.sub(r"/\*.*?\*/", "", s, flags=re.S)
            # replace unicode ellipsis placeholders like {…} with {}
            s = s.replace("{…}", "{}").replace("{…}", "{}").replace("…", "null")
            # remove trailing commas before } or ]
            s = re.sub(r",(\s*[}\]])", r"\1", s)
            return s

        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            # fallback: tolerate comments / trailing commas / {…}
            try:
                txt = _strip_json_comments_and_fix(path.read_text(encoding="utf-8"))
                raw = json.loads(txt)
            except Exception as e:
                QMessageBox.warning(self, "Compound Dictionary", f"Could not load compounds:\n{e}")
                return []

        if not isinstance(raw, list):
            return []

        out: list[dict] = []
        for c in raw:
            if not isinstance(c, dict):
                continue
            c.setdefault("name_display", c.get("name", ""))
            c.setdefault("composition", [])
            out.append(c)
        return out

    def _collect_all_sections(self):
        # Initialisiere self._all_sections aus allen Compounds und ggf. gespeicherten leeren Kategorien
        self._all_sections = set()
        for c in self.compounds:
            s = (c.get("section") or "").strip()
            if s:
                self._all_sections.add(s)
        # Lade persistente leere Kategorien (falls gewünscht, hier im RAM)
        # Optional: persistente Speicherung in Datei möglich

    def _save_compounds(self):
        try:
            self._compounds_path.write_text(
                json.dumps(self.compounds, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
        except Exception as e:
            QMessageBox.critical(self, "Compound Dictionary", f"Could not save compounds:\n{e}")
        self._collect_all_sections()

    def _populate_section_tree(self):
        from collections import defaultdict
        sections = defaultdict(list)
        for idx, compound in enumerate(self.compounds):
            key = compound.get("section") or "Uncategorized"
            sections[key].append(idx)
        # Ergänze leere Kategorien
        for s in self._all_sections:
            if s not in sections:
                sections[s] = []
        for section, indices in sorted(sections.items(), key=lambda item: item[0].lower()):
            parent = QTreeWidgetItem([section])
            parent.setFlags(parent.flags() & ~Qt.ItemFlag.ItemIsSelectable)
            self.section_tree.addTopLevelItem(parent)
            for idx in sorted(indices, key=lambda i: self.compounds[i].get("name_display", "").lower()):
                child = QTreeWidgetItem([self.compounds[idx].get("name_display", self.compounds[idx].get("name", "Unnamed"))])
                child.setData(0, Qt.ItemDataRole.UserRole, idx)
                parent.addChild(child)
            parent.setExpanded(False)

    def _populate_alpha_list(self):
        self.alpha_indices = sorted(
            range(len(self.compounds)),
            key=lambda i: self.compounds[i].get("name_display", self.compounds[i].get("name", "")).lower(),
        )
        for idx in self.alpha_indices:
            name = self.compounds[idx].get("name_display", self.compounds[idx].get("name", "Unnamed"))
            item = QListWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, idx)
            self.alpha_list.addItem(item)

    def _handle_section_selection(self):
        item = self.section_tree.currentItem()
        if not item or item.childCount():
            self._set_current_index(None)
            return
        idx = item.data(0, Qt.ItemDataRole.UserRole)
        self._set_current_index(idx)

    def _handle_section_double_click(self, item, _col):
        if item and not item.childCount():
            self._emit_selection()

    def _handle_alpha_selection(self, row):
        if row < 0:
            self._set_current_index(None)
            return
        item = self.alpha_list.item(row)
        idx = item.data(Qt.ItemDataRole.UserRole)
        self._set_current_index(idx)

    def _handle_alpha_double_click(self, _item):
        self._emit_selection()

    def _set_current_index(self, idx):
        if idx is None:
            self.current_index = None
            self.add_btn.setEnabled(False)
            self.info_label.setText("Select a compound to view its details.")
            return
        self.current_index = idx
        compound = self.compounds[idx]
        self.info_label.setText(self._format_compound_info(compound))
        self.add_btn.setEnabled(True)

    def _format_compound_info(self, compound):
        name = html.escape(compound.get("name_display", compound.get("name", "Unnamed")))
        section = html.escape(compound.get("section", "Uncategorized"))
        density = compound.get("density_g_cm3")
        density_star = compound.get("density_g_cm3_star")
        density_text = []
        if density is not None:
            density_text.append(f"{density} g/cm³")
        if density_star is not None and density_star != density:
            density_text.append(f"(★ {density_star} g/cm³)")
        composition_lines = []
        total = 0.0
        parts = []
        for part in compound.get("composition", []):
            fraction = part.get("fraction")
            if fraction is None:
                continue
            parts.append(fraction)
            total += float(fraction)
        for part, fraction in zip(compound.get("composition", []), parts):
            z = part.get("Z")
            if z is None:
                continue
            percent = (float(fraction) / total * 100.0) if total else 0.0
            composition_lines.append(f"Z={z}: {fraction} ({percent:.2f}%)")
        notes = compound.get("notes") or []
        kv_items = compound.get("kv") or {}
        html_lines = [
            f"<h3>{name}</h3>",
            f"<p><b>Section:</b> {section}</p>",
        ]
        if density_text:
            html_lines.append(f"<p><b>Density:</b> {' '.join(html.escape(x) for x in density_text)}</p>")
        if composition_lines:
            html_lines.append("<p><b>Composition:</b><br>" + "<br>".join(html.escape(line) for line in composition_lines) + "</p>")
        if kv_items:
            html_lines.append("<p><b>Key Values:</b><br>" + "<br>".join(f"{html.escape(str(k))}: {html.escape(str(v))}" for k, v in kv_items.items()) + "</p>")
        if notes:
            html_lines.append("<p><b>Notes:</b><br>" + "<br>".join(html.escape(line) for line in notes if line) + "</p>")
        return "".join(html_lines)

    def _emit_selection(self):
        if self.current_index is None:
            return
        self.compound_selected.emit(self.compounds[self.current_index])
        self.accept()

    def _install_edit_toolbar(self):
        if not self._editable:
            return

        lay = self.layout()
        if lay is None:
            lay = QVBoxLayout(self)

        bar = QHBoxLayout()

        def _theme_icon_multi(names: list[str], fallback: QStyle.StandardPixmap) -> QIcon:
            for n in names:
                ico = QIcon.fromTheme(n)
                if not ico.isNull():
                    return ico
            return self.style().standardIcon(fallback)

        def _icon_btn(icon: QIcon, tooltip: str) -> QToolButton:
            b = QToolButton()
            b.setIcon(icon)
            b.setToolTip(tooltip)
            b.setAutoRaise(True)
            return b

        # add: plus
        self.btn_add = _icon_btn(
            _theme_icon_multi(["list-add", "add", "plus", "document-new"], QStyle.StandardPixmap.SP_FileDialogNewFolder),
            "Add compound",
        )
        # edit: pen
        self.btn_edit = _icon_btn(
            _theme_icon_multi(["document-edit", "edit-rename", "accessories-text-editor"], QStyle.StandardPixmap.SP_FileDialogDetailedView),
            "Edit selected compound",
        )
        # delete: trash can
        self.btn_delete = _icon_btn(
            _theme_icon_multi(["user-trash", "edit-delete", "trash"], QStyle.StandardPixmap.SP_TrashIcon),
            "Delete selected compound",
        )
        # add category: folder+plus
        self.btn_add_category = _icon_btn(
            self._folder_badge_icon(254, badge="plus", folder_fill="#289C0B", badge_border="#2D9C0B", badge_fill="#2D9C0B"), 
            "Add category",
        )
        # delete category: folder-minus icon preferred
        self.btn_delete_category = _icon_btn(
            self._folder_badge_icon(254, badge="minus", folder_fill="#AA2828"), 
            "Delete category (and all its compounds)",
        )

        self.btn_restore = QPushButton("Restore defaults")
        self.btn_restore.setMaximumWidth(self.btn_restore.sizeHint().width())

        bar.addWidget(self.btn_add)
        bar.addWidget(self.btn_edit)
        bar.addWidget(self.btn_delete)
        bar.addWidget(self.btn_add_category)
        bar.addWidget(self.btn_delete_category)
        bar.addSpacing(8)
        bar.addWidget(self.btn_restore)
        bar.addStretch(1)
        lay.addLayout(bar)

        self.btn_add.clicked.connect(self._add_compound)
        self.btn_edit.clicked.connect(self._edit_selected_compound)
        self.btn_delete.clicked.connect(self._delete_selected_compound)
        self.btn_restore.clicked.connect(self._restore_defaults)
        self.btn_add_category.clicked.connect(self._add_category_main)
        self.btn_delete_category.clicked.connect(self._delete_category_main)

    def _existing_sections(self) -> list[str]:
        # Gibt alle Kategorien zurück, auch leere
        return sorted(set(self._all_sections) | {"Custom"}, key=str.lower)

    def _add_category_main(self):
        text, ok = QInputDialog.getText(self, "Add category", "Category name:")
        if not ok:
            return
        s = (text or "").strip()
        if not s:
            return
        if s not in self._all_sections:
            self._all_sections.add(s)
            self._rebuild_compound_view()
        QMessageBox.information(self, "Kategorie hinzugefügt", f"Die Kategorie '{s}' ist nun auswählbar und sichtbar.")

    def _delete_category_main(self):
        # Zeige Auswahl-Dialog für Kategorie
        if not self._all_sections:
            QMessageBox.information(self, "Keine Kategorien", "Es sind keine Kategorien vorhanden.")
            return
        cats = sorted(self._all_sections)
        cat, ok = QInputDialog.getItem(self, "Kategorie löschen", "Kategorie auswählen:", cats, editable=False)
        if not ok or not cat:
            return
        if QMessageBox.question(self, "Kategorie löschen", f"Soll die Kategorie '{cat}' und alle darin enthaltenen Compounds gelöscht werden?") != QMessageBox.StandardButton.Yes:
            return
        # Lösche alle Compounds dieser Kategorie
        self.compounds = [c for c in self.compounds if (c.get("section") or "").strip() != cat]
        self._all_sections.discard(cat)
        self._save_compounds()
        self._rebuild_compound_view()

    def _add_compound(self):
        dlg = _CompoundEditDialog(self, initial={"section": "Custom"}, sections=self._existing_sections())
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        comp = dlg.result_compound()
        # Neue Kategorie ggf. merken
        sec = (comp.get("section") or "").strip()
        if sec:
            self._all_sections.add(sec)
        self.compounds.append(comp)
        self._save_compounds()
        self._rebuild_compound_view()

    def _edit_selected_compound(self):
        idx = self._selected_compound_index()
        if not (0 <= idx < len(self.compounds)):
            return
        dlg = _CompoundEditDialog(self, initial=self.compounds[idx], sections=self._existing_sections())
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        comp = dlg.result_compound()
        sec = (comp.get("section") or "").strip()
        if sec:
            self._all_sections.add(sec)
        self.compounds[idx] = comp
        self._save_compounds()
        self._rebuild_compound_view()

    def _delete_selected_compound(self):
        idx = self._selected_compound_index()
        if not (0 <= idx < len(self.compounds)):
            return
        name = self.compounds[idx].get("name_display") or self.compounds[idx].get("name") or "this compound"
        if QMessageBox.question(self, "Delete", f"Delete '{name}'?") != QMessageBox.StandardButton.Yes:
            return
        self.compounds.pop(idx)
        self._save_compounds()
        self._rebuild_compound_view()

    def _restore_defaults(self):
        if not self._default_path.exists():
            QMessageBox.warning(self, "Restore defaults", "No default backup found.")
            return
        if QMessageBox.question(self, "Restore defaults", "Restore original compounds and overwrite your edits?") != QMessageBox.StandardButton.Yes:
            return
        try:
            shutil.copyfile(self._default_path, self._compounds_path)
        except Exception as e:
            QMessageBox.critical(self, "Restore defaults", f"Failed:\n{e}")
            return
        self.compounds = self._load_compounds(self._compounds_path)
        self._rebuild_compound_view()

    def _selected_compound_index(self) -> int:
        return int(self.current_index) if self.current_index is not None else -1

    def _rebuild_compound_view(self):
        self.section_tree.clear()
        self.alpha_list.clear()
        self._populate_section_tree()
        self._populate_alpha_list()
        self._set_current_index(None)

    @staticmethod
    def _folder_badge_icon(
        size: int = 32,
        badge: str = "minus",                          # "minus" | "plus"
        folder_fill: Optional[Union[QColor, str]] = None,  # QColor oder "#RRGGBB"
        badge_fill: Union[QColor, str] = QColor(220, 0, 0),
        badge_border: Union[QColor, str] = QColor(140, 0, 0),
    ) -> QIcon:
        app = QApplication.instance()
        style = app.style()

        # Basis: symbolic folder aus Theme, sonst Qt-Fallback
        base = QIcon.fromTheme("folder-symbolic")
        if base.isNull():
            base = QIcon.fromTheme("inode-directory-symbolic")
        if base.isNull():
            base = style.standardIcon(QStyle.StandardPixmap.SP_DirIcon)

        pm = base.pixmap(QSize(size, size))
        out = QPixmap(pm.size())
        out.fill(Qt.GlobalColor.transparent)

        p = QPainter(out)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)

        # --- optional: folder tint (Füllfarbe) ---
        if folder_fill is not None:
            tint = folder_fill if isinstance(folder_fill, QColor) else QColor(folder_fill)
            # Copy, then tint alpha mask
            tinted = QPixmap(pm.size())
            tinted.fill(Qt.GlobalColor.transparent)

            pt = QPainter(tinted)
            pt.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            pt.drawPixmap(0, 0, pm)
            pt.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceIn)
            pt.fillRect(tinted.rect(), tint)
            pt.end()

            p.drawPixmap(0, 0, tinted)
        else:
            p.drawPixmap(0, 0, pm)

        s = float(size)

        # Badge Geometrie (unten rechts)
        badge_d = max(12.0, s * 0.5)
        margin  = max(1.0,  s * 0.06)
        badge_rect = QRectF(s - badge_d - margin, s - badge_d - margin, badge_d, badge_d)

        border_w = max(1.0, s * 0.04)
        bb = badge_border if isinstance(badge_border, QColor) else QColor(badge_border)
        bf = badge_fill if isinstance(badge_fill, QColor) else QColor(badge_fill)

        p.setPen(QPen(bb, border_w))
        p.setBrush(QBrush(bf))
        p.drawEllipse(badge_rect)

        # Symbol im Badge (weiß)
        cx = badge_rect.center().x()
        cy = badge_rect.center().y()
        inset = badge_d * 0.28
        x1 = cx - (badge_d / 2.0 - inset)
        x2 = cx + (badge_d / 2.0 - inset)

        stroke_w = max(2.0, s * 0.085)
        pen = QPen(Qt.GlobalColor.white, stroke_w)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        p.setPen(pen)

        # Horizontale Linie (Minus und Plus gemeinsam)
        p.drawLine(QPointF(x1, cy), QPointF(x2, cy))

        # Vertikale Linie nur für Plus
        if str(badge).lower() == "plus":
            y1 = cy - (badge_d / 2.0 - inset)
            y2 = cy + (badge_d / 2.0 - inset)
            p.drawLine(QPointF(cx, y1), QPointF(cx, y2))

        p.end()
        return QIcon(out)