from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

import json
import os
import sys


def _resolve_app_config_dir() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "opensrim"
    if sys.platform.startswith("win"):
        base = os.environ.get("APPDATA")
        if base:
            return Path(base) / "opensrim"
        return Path.home() / "AppData" / "Roaming" / "opensrim"
    base = os.environ.get("XDG_CONFIG_HOME")
    if base:
        return Path(base) / "opensrim"
    return Path.home() / ".config" / "opensrim"


def _resolve_app_config_path(scope: str = "global") -> Path:
    filename = "settings.toml" if scope == "global" else f"{scope}.toml"
    return _resolve_app_config_dir() / filename


def _load_toml_module():
    try:
        import tomllib  # type: ignore
        return tomllib
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore
            return tomllib
        except ImportError:
            return None


def _load_app_settings(scope: str = "global") -> Dict[str, Any]:
    path = _resolve_app_config_path(scope)
    if not path.is_file():
        return {}
    toml_module = _load_toml_module()
    if toml_module is None:
        return {}
    try:
        with open(path, "rb") as fh:
            data = toml_module.load(fh)
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return repr(float(value))
    return json.dumps(str(value), ensure_ascii=False)


def _write_toml_lines(prefix: str, payload: Dict[str, Any], lines: List[str]) -> None:
    scalar_items: list[tuple[str, Any]] = []
    nested_items: list[tuple[str, Dict[str, Any]]] = []
    for key, value in payload.items():
        if isinstance(value, dict):
            nested_items.append((str(key), value))
        elif isinstance(value, (str, int, float, bool)):
            scalar_items.append((str(key), value))

    if prefix:
        lines.append(f"[{prefix}]")
    for key, value in scalar_items:
        lines.append(f"{key} = {_toml_value(value)}")

    if scalar_items and nested_items:
        lines.append("")

    for index, (key, nested) in enumerate(nested_items):
        next_prefix = f"{prefix}.{key}" if prefix else key
        _write_toml_lines(next_prefix, nested, lines)
        if index < len(nested_items) - 1:
            lines.append("")


def _save_app_settings(payload: Dict[str, Any], scope: str = "global") -> None:
    path = _resolve_app_config_path(scope)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        lines: List[str] = []
        _write_toml_lines("", payload, lines)
        text = "\n".join(line for line in lines if line is not None).strip()
        path.write_text((text + "\n") if text else "", encoding="utf-8")
    except OSError:
        pass


def _coerce_directory(path_value: str | Path | None) -> Optional[str]:
    if not path_value:
        return None
    try:
        path = Path(path_value).expanduser()
    except Exception:
        return None
    directory = path if path.is_dir() else path.parent
    if not str(directory):
        return None
    if directory.exists():
        return str(directory)
    parent = directory.parent
    if parent and parent.exists():
        return str(parent)
    return None


def get_last_used_directory(fallback: str | Path | None = None) -> str:
    payload = _load_app_settings("global")
    stored = _coerce_directory(payload.get("last_directory"))
    if stored:
        return stored
    fallback_dir = _coerce_directory(fallback)
    if fallback_dir:
        return fallback_dir
    return str(Path.home())


def remember_last_used_path(path_value: str | Path | None) -> None:
    directory = _coerce_directory(path_value)
    if not directory:
        return
    payload = _load_app_settings("global")
    payload["last_directory"] = directory
    _save_app_settings(payload, "global")


def get_persisted_display_settings() -> Dict[str, Any]:
    payload = _load_app_settings("mc_results")
    display = payload.get("display_settings")
    return dict(display) if isinstance(display, dict) else {}


def set_persisted_display_settings(settings: Dict[str, Any]) -> None:
    payload = _load_app_settings("mc_results")
    payload["display_settings"] = {
        str(key): value
        for key, value in settings.items()
        if isinstance(value, (str, int, float, bool))
    }
    _save_app_settings(payload, "mc_results")


def load_scoped_settings(scope: str) -> Dict[str, Any]:
    """Read a page-scoped settings file (e.g. "koral", "mc_setup", "single_plot").

    Each page/simulator gets its own <scope>.toml under the app config dir,
    so a user (or an editor) can see at a glance which file holds which
    page's persisted preferences, instead of one shared settings.toml.
    """
    return _load_app_settings(scope)


def save_scoped_settings(payload: Dict[str, Any], scope: str) -> None:
    """Write a page-scoped settings file. See load_scoped_settings()."""
    _save_app_settings(payload, scope)


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
    unit_options: List[str] = field(default_factory=lambda: ["Å", "nm", "µm", "mm", "cm", "m", "km"])
    energy_defaults: Dict[str, str] = field(default_factory=lambda: {"damage": "25", "disp": "25", "latt": "3", "surf": "3"})
    elements_by_number: Dict[int, Dict[str, Any]] = field(default_factory=load_elements_by_number)

    log_entries: List[str] = field(default_factory=list)
    current_config_path: Optional[str] = None

    def add_log(self, message: str) -> str:
        timestamp = datetime.now().strftime("%H:%M:%S")
        entry = f"[{timestamp}] {message}"
        self.log_entries.append(entry)
        self.log_entries = self.log_entries[-1000:]
        return entry

    def clear_logs(self) -> None:
        self.log_entries.clear()

    def get_unit_options(self) -> List[str]:
        return list(self.unit_options)

    def get_dialog_start_directory(self, fallback: str | Path | None = None) -> str:
        if fallback is None and self.current_config_path:
            fallback = Path(self.current_config_path).parent
        return get_last_used_directory(fallback)

    def remember_dialog_path(self, path_value: str | Path | None) -> None:
        remember_last_used_path(path_value)

    def get_persisted_display_settings(self) -> Dict[str, Any]:
        return get_persisted_display_settings()

    def set_persisted_display_settings(self, settings: Dict[str, Any]) -> None:
        set_persisted_display_settings(settings)
