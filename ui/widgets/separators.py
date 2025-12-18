from PyQt6.QtWidgets import QFrame


class Line(QFrame):
    """A thin divider line for visual grouping."""
    def __init__(self):
        super().__init__()
        self.setFrameShape(QFrame.Shape.HLine)
        self.setFrameShadow(QFrame.Shadow.Sunken)
        self.setObjectName("divider")
