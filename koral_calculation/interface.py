from __future__ import annotations

from dataclasses import dataclass
import importlib
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    try:
        import tomli as tomllib  # type: ignore
    except ModuleNotFoundError:  # pragma: no cover
        tomllib = None  # type: ignore


@dataclass(frozen=True)
class DiscoveredModel:
    """A model discovered on disk.

    `id` is the folder name.
    A model is considered primary if a `.primarymodel` file exists in the folder.
    """

    id: str
    path: Path
    is_primary: bool


@runtime_checkable
class KoralCalculationModel(Protocol):
    """Interface contract for a KORAL calculation model.

    Implementation details will be added when the backend is embedded.
    """

    @property
    def id(self) -> str: ...

    def run(self, request: dict[str, Any]) -> dict[str, Any]: ...


def load_model(model_id: str) -> KoralCalculationModel:
    """Load a calculation model implementation.

    Convention:
    - module: koral_calculation.<model_id>.model
    - provides either:
      - MODEL: an instance implementing KoralCalculationModel, or
      - create_model(): function returning such an instance, or
      - Model: class implementing the protocol.
    """

    mod = importlib.import_module(f"koral_calculation.{model_id}.model")
    if hasattr(mod, "MODEL"):
        model = getattr(mod, "MODEL")
        return model
    if hasattr(mod, "create_model"):
        model = getattr(mod, "create_model")()
        return model
    if hasattr(mod, "Model"):
        model = getattr(mod, "Model")()
        return model
    raise ImportError(f"Model '{model_id}' does not export MODEL/create_model/Model")


def filter_outputs(result: dict[str, Any], requested_outputs: list[str]) -> dict[str, Any]:
    """Filter a model result to include only requested outputs.

    Expected result schema:
      {
        'model_id': str,
        'energies_keV': [..],
        'outputs': { output_id: [..values..], ... }
      }
    """

    outputs = result.get("outputs")
    if not isinstance(outputs, dict):
        return result
    filtered = {k: outputs[k] for k in requested_outputs if k in outputs}
    result = dict(result)
    result["outputs"] = filtered
    return result


def models_root() -> Path:
    return Path(__file__).resolve().parent


def discover_models(root: Path | None = None) -> list[DiscoveredModel]:
    """Discover available KORAL calculation models.

    Models are direct subdirectories of `koral_calculation/`.
    Hidden folders and __pycache__ are ignored.
    """

    base = root or models_root()
    out: list[DiscoveredModel] = []

    if not base.exists() or not base.is_dir():
        return out

    for child in base.iterdir():
        if not child.is_dir():
            continue
        name = child.name
        if name.startswith(".") or name == "__pycache__":
            continue
        if name in {"__init__.py"}:
            continue
        is_primary = (child / ".primarymodel").exists()
        out.append(DiscoveredModel(id=name, path=child, is_primary=is_primary))

    out.sort(key=lambda m: m.id.lower())
    return out


UI_PARAMS_FILENAME = "ui_params.toml"
UI_PARAMS_FORMAT_VERSION = 1


def load_ui_parameters(model_id: str, *, root: Path | None = None) -> dict[str, dict[str, Any]]:
    """Load optional per-model UI parameter metadata.

    A model may provide a `ui_params.toml` file next to its `model.py` to let the
    KORAL UI tune input widgets (min/max, step, decimals, defaults, ...).

    The file is optional. Returns an empty dict when missing or invalid.
    """

    if tomllib is None:
        return {}

    base = root or models_root()
    path = base / str(model_id) / UI_PARAMS_FILENAME
    if not path.exists():
        return {}

    try:
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}

    if not isinstance(raw, dict):
        return {}

    fmt = raw.get("format_version", UI_PARAMS_FORMAT_VERSION)
    try:
        fmt_i = int(fmt)
    except (TypeError, ValueError):
        return {}
    if fmt_i != UI_PARAMS_FORMAT_VERSION:
        return {}

    params = raw.get("parameters")
    if not isinstance(params, dict):
        return {}

    out: dict[str, dict[str, Any]] = {}
    for param_id, spec in params.items():
        if not isinstance(param_id, str) or not param_id.strip():
            continue
        if not isinstance(spec, dict):
            continue
        if "min" not in spec or "max" not in spec:
            continue
        try:
            vmin = float(spec["min"])
            vmax = float(spec["max"])
        except (TypeError, ValueError):
            continue

        normalized: dict[str, Any] = {"min": vmin, "max": vmax}

        for key in ("decimals", "step", "default"):
            if key not in spec:
                continue
            val = spec.get(key)
            if val is None:
                continue
            try:
                if key == "decimals":
                    normalized[key] = int(val)
                else:
                    normalized[key] = float(val)
            except (TypeError, ValueError):
                continue

        for key in ("label", "unit", "description"):
            val = spec.get(key)
            if isinstance(val, str) and val.strip():
                normalized[key] = val.strip()

        out[param_id] = normalized

    return out
