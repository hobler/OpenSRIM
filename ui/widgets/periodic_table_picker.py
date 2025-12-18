import json
import os
from PyQt6.QtWidgets import (QDialog, QWidget, QGridLayout, QPushButton, 
                             QVBoxLayout, QLabel, QScrollArea, QHBoxLayout, QFrame,
                             QSizePolicy)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont, QColor, QPainter, QPen

class PeriodicTableDialog(QDialog):
    """Dialog for selecting elements from the periodic table"""
    element_selected = pyqtSignal(dict)
    
    def __init__(self, parent=None, compact=True, show_hover_info=True, bordered=False):
        super().__init__(parent)
        self.compact = compact
        self.show_hover_info = show_hover_info
        self.bordered = bordered
        self.setObjectName("periodicTableDialog")
        self.setWindowTitle("Periodic Table of Elements")
        self.setModal(True)
        
        if compact:
            # Calculate exact size based on content with 1px spacing
            width = 18 * 40 + 17 * 1 + 12 + 24 + 10  # tiles + spacing + label + margins
            height = (7 * 40 + 6 * 1 + 12) + 8 + (2 * 40 + 8) + 8 + 60 + 24  # adjusted for 2-row legend
            if show_hover_info:
                height += 120  # Extra space for hover info
            self.resize(width, height)
        else:
            # Start in fullscreen mode for expanded version
            self.showMaximized()
        
        base_style = "background-color: white;"
        border_style = "border: 2px solid #000;" if self.bordered else ""
        self.setStyleSheet(f"QDialog#periodicTableDialog {{ {base_style} {border_style} }}")
        
        # Load elements from JSON
        self.elements = self._load_elements()
        
        # Setup UI
        self._setup_ui()
        
    def _load_elements(self):
        """Load elements from JSON file"""
        json_path = os.path.join(os.path.dirname(__file__), 'PeriodicTableJSON.json')
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return {elem['number']: elem for elem in data['elements']}
    
    def _get_element_group(self, element):
        """Determine element group for color coding"""
        category = element['category'].lower()
        number = element['number']
        
        # Special case: Hydrogen is a nonmetal, not alkali metal
        if number == 1:
            return 'nonmetal'
        
        if 'alkali metal' in category and 'alkaline' not in category:
            return 'alkali_metal'
        elif 'alkaline earth' in category:
            return 'alkaline_earth_metal'
        elif 'lanthanide' in category:
            return 'lanthanide'
        elif 'actinide' in category:
            return 'actinide'
        elif 'transition metal' in category:
            return 'transition_metal'
        elif 'post-transition metal' in category:
            return 'post_transition_metal'
        elif 'metalloid' in category:
            return 'metalloid'
        elif 'noble gas' in category:
            return 'noble_gas'
        elif 'halogen' in category or number in [9, 17, 35, 53, 85]:
            return 'halogen'
        elif 'nonmetal' in category or 'diatomic nonmetal' in category:
            return 'nonmetal'
        else:
            return 'unknown'
    
    def _setup_ui(self):
        """Setup the user interface"""
        main_layout = QVBoxLayout(self)
        
        if self.compact:
            main_layout.setContentsMargins(12, 12, 12, 12)
        else:
            main_layout.setContentsMargins(20, 20, 20, 20)
        
        main_layout.setSpacing(0)
        
        # Main content (no scroll area for compact)
        if self.compact:
            content_widget = QWidget()
            content_layout = QVBoxLayout(content_widget)
            content_layout.setSpacing(8)
            content_layout.setContentsMargins(0, 0, 0, 0)
            
            # Main periodic table (with embedded hover info)
            main_table = self._create_main_table()
            content_layout.addWidget(main_table)
            
            # Lanthanides and Actinides section
            content_layout.addSpacing(8)
            series_section = self._create_series_section()
            content_layout.addWidget(series_section)
            
            # Legend at bottom
            content_layout.addSpacing(8)
            legend = self._create_legend()
            content_layout.addWidget(legend, alignment=Qt.AlignmentFlag.AlignCenter)
            
            main_layout.addWidget(content_widget)
        else:
            # Main content in scroll area for expanded
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            scroll.setStyleSheet("QScrollArea { border: none; background-color: white; }")
            
            content_widget = QWidget()
            content_layout = QVBoxLayout(content_widget)
            content_layout.setSpacing(16)
            content_layout.setContentsMargins(0, 0, 0, 0)
            
            # Main periodic table (with embedded hover info)
            main_table = self._create_main_table()
            content_layout.addWidget(main_table)
            
            # Lanthanides and Actinides section
            content_layout.addSpacing(16)
            series_section = self._create_series_section()
            content_layout.addWidget(series_section)
            
            # Legend at bottom
            content_layout.addSpacing(12)
            legend = self._create_legend()
            content_layout.addWidget(legend, alignment=Qt.AlignmentFlag.AlignCenter)
            
            content_layout.addStretch()
            
            scroll.setWidget(content_widget)
            main_layout.addWidget(scroll)
    
    def _create_hover_info(self):
        """Create hover information display panel embedded in grid"""
        info_panel = QFrame()
        info_panel.setObjectName("hoverInfoPanel")
        info_panel.setStyleSheet("""
            QFrame#hoverInfoPanel {
                background-color: transparent;
                border: none;
                border-radius: 6px;
                padding: 8px;
            }
        """)
        info_panel.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        info_panel.setMinimumHeight(120 if self.compact else 170)
        
        layout = QVBoxLayout(info_panel)
        layout.setContentsMargins(10, 8, 10, 8)
        layout.setSpacing(8 if self.compact else 12)
        
        font_size_main = 12 if self.compact else 18
        font_size_symbol = font_size_main + 4 if self.compact else font_size_main + 6
        font_size_detail = 9 if self.compact else 10
        
        main_layout = QHBoxLayout()
        main_layout.setSpacing(12 if self.compact else 20)
        
        self.info_number = QLabel("—")
        self.info_number.setFont(QFont("Inter", font_size_main, QFont.Weight.Bold))
        self.info_number.setStyleSheet("color: #111827; background-color: transparent;")
        self.info_number.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_number.setMinimumWidth(46 if self.compact else 56)
        main_layout.addWidget(self.info_number)
        
        self.info_symbol = QLabel("—")
        self.info_symbol.setFont(QFont("Inter", font_size_symbol, QFont.Weight.Bold))
        self.info_symbol.setStyleSheet("color: #111827; background-color: transparent;")
        self.info_symbol.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_symbol.setMinimumWidth(72 if self.compact else 96)
        main_layout.addWidget(self.info_symbol)
        
        self.info_name = QLabel("—")
        self.info_name.setFont(QFont("Inter", font_size_main, QFont.Weight.Bold))
        self.info_name.setStyleSheet("color: #111827; background-color: transparent;")
        self.info_name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_name.setMinimumWidth(110 if self.compact else 160)
        main_layout.addWidget(self.info_name)
        
        main_layout.addStretch()
        layout.addLayout(main_layout)
        
        detail_grid = QGridLayout()
        detail_grid.setContentsMargins(0, 0, 0, 0)
        detail_grid.setHorizontalSpacing(10 if self.compact else 16)
        detail_grid.setVerticalSpacing(4 if self.compact else 6)
        
        detail_columns = 2 if self.compact else 3
        self.info_labels = {}
        info_items = [
            ("atomic_mass", "Mass:", "amu"),
            ("density", "Density:", "g/cm³"),
            ("melt", "M.P.:", "K"),
            ("boil", "B.P.:", "K"),
            ("electron_affinity", "E.Aff.:", "kJ/mol"),
            ("electronegativity_pauling", "E.neg.:", ""),
        ]
        
        row, col = 0, 0
        for key, label_text, unit in info_items:
            combined_label = QLabel(f"{label_text} —")
            combined_label.setFont(QFont("Inter", font_size_detail))
            combined_label.setStyleSheet("color: #111827; background-color: transparent;")
            combined_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            combined_label.setWordWrap(False)
            combined_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            combined_label.setMinimumHeight(18 if self.compact else 22)
            detail_grid.addWidget(combined_label, row, col)
            
            self.info_labels[key] = (combined_label, label_text, unit)
            
            col += 1
            if col >= detail_columns:
                col = 0
                row += 1
        
        layout.addLayout(detail_grid)
        return info_panel
    
    def _update_hover_info(self, element):
        """Update hover info panel with element data"""
        if not self.show_hover_info:
            return
        
        # Check if hover info widgets exist
        if not hasattr(self, 'info_number'):
            return
        if not hasattr(self, 'info_symbol'):
            return
        if not hasattr(self, 'info_name'):
            return
        if not hasattr(self, 'info_labels'):
            return
        
        # Update main info (no labels)
        self.info_number.setText(str(element.get('number', '—')))
        self.info_symbol.setText(element.get('symbol', '—'))
        self.info_name.setText(element.get('name', '—'))
        
        # Update additional info (format: "Label: Value Unit")
        for key, (label_widget, label_text, unit) in self.info_labels.items():
            value = element.get(key, "—")
            if value is None or value == "":
                text = f"{label_text} —"
            else:
                if isinstance(value, float):
                    value_str = f"{value:.3f}"
                else:
                    value_str = str(value)
                
                if unit:
                    text = f"{label_text} {value_str} {unit}"
                else:
                    text = f"{label_text} {value_str}"
            
            label_widget.setText(text)
    
    def _create_legend(self):
        """Create color legend (compact version at bottom)"""
        legend_widget = QFrame()
        legend_widget.setStyleSheet("""
            QFrame {
                background-color: #F8F9FA;
                border: 1px solid #DEE2E6;
                border-radius: 6px;
                padding: 8px;
            }
        """)
        
        group_colors = self._get_group_colors()
        group_names = {
            'alkali_metal': 'Alkali',
            'alkaline_earth_metal': 'Alk. Earth',
            'transition_metal': 'Transition',
            'post_transition_metal': 'Post-trans.',
            'metalloid': 'Metalloid',
            'nonmetal': 'Nonmetal',
            'halogen': 'Halogen',
            'noble_gas': 'Noble Gas',
            'lanthanide': 'Lanthanide',
            'actinide': 'Actinide'
        }
        
        if self.compact:
            # Grid layout for two rows in compact mode
            grid_layout = QGridLayout(legend_widget)
            grid_layout.setSpacing(8)
            grid_layout.setContentsMargins(8, 8, 8, 8)
            
            row = 0
            col = 0
            for group, name in group_names.items():
                item_layout = QHBoxLayout()
                item_layout.setSpacing(4)
                
                color_box = QLabel()
                color_box.setFixedSize(12, 12)
                color_box.setStyleSheet(f"""
                    background-color: {group_colors[group]};
                    border-radius: 2px;
                    border: 1px solid rgba(0, 0, 0, 0.1);
                """)
                
                name_label = QLabel(name)
                name_label.setFont(QFont("Inter", 8))
                name_label.setStyleSheet("color: #495057; background-color: transparent; border: none;")
                
                item_layout.addWidget(color_box)
                item_layout.addWidget(name_label)
                
                container = QWidget()
                container.setLayout(item_layout)
                container.setStyleSheet("background-color: transparent;")
                
                grid_layout.addWidget(container, row, col)
                
                col += 1
                if col >= 5:  # 5 items per row
                    col = 0
                    row += 1
        else:
            # Single row layout for expanded mode
            layout = QHBoxLayout(legend_widget)
            layout.setSpacing(16)
            layout.setContentsMargins(12, 8, 12, 8)
            
            for group, name in group_names.items():
                item_layout = QHBoxLayout()
                item_layout.setSpacing(6)
                
                color_box = QLabel()
                color_box.setFixedSize(16, 16)
                color_box.setStyleSheet(f"""
                    background-color: {group_colors[group]};
                    border-radius: 3px;
                    border: 1px solid rgba(0, 0, 0, 0.1);
                """)
                
                name_label = QLabel(name)
                name_label.setFont(QFont("Inter", 10))
                name_label.setStyleSheet("color: #495057; background-color: transparent; border: none;")
                
                item_layout.addWidget(color_box)
                item_layout.addWidget(name_label)
                
                container = QWidget()
                container.setLayout(item_layout)
                container.setStyleSheet("background-color: transparent;")
                
                layout.addWidget(container)
            
            layout.addStretch()
        
        return legend_widget
    
    def _get_group_colors(self):
        """Get ACS-style color mapping"""
        return {
            'alkali_metal': '#F39C12',
            'alkaline_earth_metal': '#F7DC6F',
            'transition_metal': '#2E86C1',
            'post_transition_metal': '#4D5B6A',
            'metalloid': '#9A8F56',
            'nonmetal': '#7ED3A6',
            'halogen': '#E6F28C',
            'noble_gas': '#8E6EC8',
            'lanthanide': '#D98AD5',
            'actinide': '#6CC3B1',
            'unknown': '#E0E0E0'
        }
    
    def _create_main_table(self):
        """Create main periodic table"""
        table_widget = QWidget()
        table_widget.setStyleSheet("background-color: white;")
        layout = QGridLayout(table_widget)
        spacing = 1 if self.compact else 6
        layout.setSpacing(spacing)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Set horizontal and vertical spacing separately for better control
        if self.compact:
            layout.setHorizontalSpacing(1)
            layout.setVerticalSpacing(1)
        else:
            layout.setHorizontalSpacing(6)
            layout.setVerticalSpacing(6)
        
        # Add group labels (1-18) at top
        label_height = 12 if self.compact else 18
        font_size = 7 if self.compact else 9
        for group in range(1, 19):
            label = QLabel(str(group))
            label.setFont(QFont("Inter", font_size, QFont.Weight.DemiBold))
            label.setStyleSheet("color: #6C757D;")
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            label.setFixedHeight(label_height)
            layout.addWidget(label, 0, group)
        
        # Add period labels (1-7) on left
        label_width = 12 if self.compact else 18
        for period in range(1, 8):
            label = QLabel(str(period))
            label.setFont(QFont("Inter", font_size, QFont.Weight.DemiBold))
            label.setStyleSheet("color: #6C757D;")
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            label.setFixedWidth(label_width)
            layout.addWidget(label, period, 0)
        
        # Add hover info panel at columns 3-12, rows 1-3 (if enabled)
        if self.show_hover_info:
            self.hover_info = self._create_hover_info()
            layout.addWidget(self.hover_info, 1, 3, 3, 10)  # span 3 rows, 10 columns
        
        # Create element tiles
        group_colors = self._get_group_colors()
        lanthanides = range(57, 72)
        actinides = range(89, 104)
        
        for number, element in self.elements.items():
            if number in lanthanides or number in actinides:
                continue
            
            # Create placeholder for lanthanides/actinides
            if element['period'] == 6 and element['group'] == 3:
                placeholder = self._create_placeholder("57-71", group_colors['lanthanide'])
                layout.addWidget(placeholder, 6, 3)
                continue
            elif element['period'] == 7 and element['group'] == 3:
                placeholder = self._create_placeholder("89-103", group_colors['actinide'])
                layout.addWidget(placeholder, 7, 3)
                continue
            
            tile = self._create_element_tile(element)
            layout.addWidget(tile, element['period'], element['group'])
        
        return table_widget
    
    def _create_placeholder(self, text, color):
        """Create placeholder for lanthanides/actinides"""
        size = 40 if self.compact else 80
        font_size = 7 if self.compact else 10
        
        placeholder = QLabel(text)
        placeholder.setFixedSize(size, size)
        placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        placeholder.setFont(QFont("Inter", font_size, QFont.Weight.Bold))
        placeholder.setStyleSheet(f"""
            background-color: {color};
            color: white;
            border-radius: {4 if self.compact else 8}px;
            border: 1px solid rgba(0, 0, 0, 0.1);
        """)
        return placeholder
    
    def _create_element_tile(self, element):
        """Create an element tile button"""
        group = self._get_element_group(element)
        color = self._get_group_colors()[group]
        
        # Determine text color based on background brightness
        text_color = self._get_text_color(color)
        
        btn = QPushButton()
        size = 40 if self.compact else 80
        btn.setFixedSize(size, size)
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        
        # Enable mouse tracking to ensure hover events work
        btn.setMouseTracking(True)
        
        border_radius = 4 if self.compact else 8
        btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {color};
                border-radius: {border_radius}px;
                border: 1px solid rgba(0, 0, 0, 0.1);
                padding: 0px;
            }}
            QPushButton:hover {{
                background-color: {self._lighten_color(color, 1.05)};
                border: 2px solid #0066CC;
            }}
        """)
        
        # Create layout for tile content
        layout = QVBoxLayout(btn)
        margin = 3 if self.compact else 5
        layout.setContentsMargins(margin, margin, margin, margin)
        layout.setSpacing(0)
        
        if self.compact:
            # Compact version - only symbol and number
            # Atomic number (top left, small)
            num_label = QLabel(str(element['number']))
            num_label.setFont(QFont("Inter", 6, QFont.Weight.DemiBold))
            num_label.setStyleSheet(f"color: {text_color}; background-color: transparent;")
            num_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
            layout.addWidget(num_label)
            
            # Symbol (dominant, centered)
            symbol_label = QLabel(element['symbol'])
            symbol_label.setFont(QFont("Inter", 16, QFont.Weight.Bold))
            symbol_label.setStyleSheet(f"color: {text_color}; background-color: transparent;")
            symbol_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(symbol_label, stretch=1, alignment=Qt.AlignmentFlag.AlignCenter)
            
        else:
            # Full version - number, symbol, name, mass
            # Atomic number
            num_label = QLabel(str(element['number']))
            num_label.setFont(QFont("Inter", 9, QFont.Weight.DemiBold))
            num_label.setStyleSheet(f"color: {text_color}; background-color: transparent;")
            num_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
            layout.addWidget(num_label)
            
            # Symbol (dominant)
            symbol_label = QLabel(element['symbol'])
            symbol_label.setFont(QFont("Inter", 28, QFont.Weight.Bold))
            symbol_label.setStyleSheet(f"color: {text_color}; background-color: transparent;")
            symbol_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(symbol_label, stretch=1)
            
            # Name
            name_label = QLabel(element['name'])
            name_label.setFont(QFont("Inter", 7))
            name_label.setStyleSheet(f"color: {text_color}; background-color: transparent;")
            name_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            name_label.setWordWrap(False)
            layout.addWidget(name_label)
            
            # Atomic mass
            mass_label = QLabel(str(element['atomic_mass']))
            mass_label.setFont(QFont("Inter", 7))
            mass_label.setStyleSheet(f"color: {text_color}; background-color: transparent; opacity: 0.9;")
            mass_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(mass_label)
        
        # Tooltip
        btn.setToolTip(f"{element['name']}\nAtomic Number: {element['number']}\n"
                      f"Atomic Mass: {element['atomic_mass']}\nCategory: {element['category']}")
        
        # Connect hover event - use a proper event override
        def on_enter(event, elem=element):
            self._update_hover_info(elem)
            event.accept()
        
        btn.enterEvent = on_enter
        btn.clicked.connect(lambda checked, e=element: self._on_element_clicked(e))
        
        return btn
    
    def _create_series_section(self):
        """Create lanthanides and actinides section"""
        series_widget = QWidget()
        series_widget.setStyleSheet("background-color: white;")
        layout = QVBoxLayout(series_widget)
        spacing = 4 if self.compact else 8
        layout.setSpacing(spacing)
        layout.setContentsMargins(0, 0, 0, 0)
        
        font_size = 8 if self.compact else 9
        label_width = 60 if self.compact else 100
        
        # Lanthanides
        lan_layout = QHBoxLayout()
        lan_layout.setSpacing(1 if self.compact else 6)
        
        lan_label = QLabel("Lanthanides")
        lan_label.setFont(QFont("Inter", font_size, QFont.Weight.Bold))
        lan_label.setStyleSheet("color: #495057;")
        lan_label.setFixedWidth(label_width)
        lan_layout.addWidget(lan_label)
        
        for number in range(57, 72):
            tile = self._create_element_tile(self.elements[number])
            lan_layout.addWidget(tile)
        
        lan_layout.addStretch()
        layout.addLayout(lan_layout)
        
        # Actinides
        act_layout = QHBoxLayout()
        act_layout.setSpacing(1 if self.compact else 6)
        
        act_label = QLabel("Actinides")
        act_label.setFont(QFont("Inter", font_size, QFont.Weight.Bold))
        act_label.setStyleSheet("color: #495057;")
        act_label.setFixedWidth(label_width)
        act_layout.addWidget(act_label)
        
        for number in range(89, 104):
            if number in self.elements:
                tile = self._create_element_tile(self.elements[number])
                act_layout.addWidget(tile)
        
        act_layout.addStretch()
        layout.addLayout(act_layout)
        
        return series_widget
    
    def _get_text_color(self, bg_color):
        """Determine text color based on background brightness"""
        color = QColor(bg_color)
        r, g, b = color.red() / 255, color.green() / 255, color.blue() / 255
        luminance = 0.299 * r + 0.587 * g + 0.114 * b
        return "white" if luminance < 0.6 else "#212529"
    
    def _lighten_color(self, color, factor):
        """Lighten a color by a factor"""
        qcolor = QColor(color)
        h, s, l, a = qcolor.getHslF()
        l = min(1.0, l * factor)
        qcolor.setHslF(h, s, l, a)
        return qcolor.name()
    
    def _on_element_clicked(self, element):
        """Handle element button click"""
        self.element_selected.emit(element)
        self.accept()
    
    def paintEvent(self, event):
        super().paintEvent(event)
        if self.bordered:
            painter = QPainter(self)
            pen = QPen(Qt.GlobalColor.black)
            pen.setWidth(2)
            painter.setPen(pen)
            painter.drawRect(self.rect().adjusted(1, 1, -1, -1))


class PeriodicTableButton(QPushButton):
    """Button that opens periodic table dialog when clicked"""
    element_selected = pyqtSignal(dict)
    
    def __init__(self, text="Select Element", parent=None, compact=True, show_hover_info=True, bordered=False, update_button_text=True):
        super().__init__(text, parent)
        self.selected_element = None
        self.compact = compact
        self.show_hover_info = show_hover_info
        self.bordered = bordered
        self.update_button_text = update_button_text
        self.clicked.connect(self._open_periodic_table)
        
    def _open_periodic_table(self):
        """Open periodic table dialog"""
        dialog = PeriodicTableDialog(self, compact=self.compact, show_hover_info=self.show_hover_info, bordered=self.bordered)
        dialog.element_selected.connect(self._on_element_selected)
        dialog.exec()
    
    def _on_element_selected(self, element):
        """Store and emit selected element"""
        self.selected_element = element
        if self.update_button_text:
            self.setText(f"{element['symbol']} ({element['name']})")
        self.element_selected.emit(element)
    
    def get_element(self):
        """Get the currently selected element"""
        return self.selected_element
