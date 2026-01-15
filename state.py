from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

import json
import sys


def _resolve_periodic_table_json() -> Optional[Path]:
    """
    Try multiple common locations.
    Returns a Path if found, otherwise None.
    """
    here = Path(__file__).resolve().parent          # .../app
    project_root = here.parent                      # .../SRIM UI

    candidates = [
        here / "ui" / "widgets" / "PeriodicTableJSON.json",  # user-provided location (app/ui/widgets)
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def load_elements_by_number() -> Dict[int, Dict[str, Any]]:
    path = _resolve_periodic_table_json()
    if not path:
        # Keep app usable even if the resource is missing.
        print(
            "WARNING: PeriodicTableJSON.json not found. "
            "Element picker will be empty until the file is restored.",
            file=sys.stderr,
        )
        return {}

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"WARNING: Failed to read/parse {path}: {exc}", file=sys.stderr)
        return {}

    elements = data.get("elements", [])
    out: Dict[int, Dict[str, Any]] = {}
    for entry in elements:
        try:
            out[int(entry["number"])] = entry
        except Exception:
            continue
    return out


@dataclass
class AppState:
    unit_options: List[str] = field(default_factory=lambda: ["Ång", "nm", "µm", "mm", "cm", "m", "km"])
    energy_defaults: Dict[str, str] = field(default_factory=lambda: {"damage": "25", "disp": "25", "latt": "3", "surf": "3"})
    elements_by_number: Dict[int, Dict[str, Any]] = field(default_factory=load_elements_by_number)

    log_entries: List[str] = field(default_factory=list)
    current_config_path: Optional[str] = None

    def add_log(self, message: str) -> str:
        timestamp = datetime.now().strftime("%H:%M:%S")
        entry = f"[{timestamp}] {message}"
        self.log_entries.append(entry)
        self.log_entries = self.log_entries[-100:]
        return entry

    def clear_logs(self) -> None:
        self.log_entries.clear()

    def get_unit_options(self) -> List[str]:
        return list(self.unit_options)
