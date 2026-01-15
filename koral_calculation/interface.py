from __future__ import annotations

from dataclasses import dataclass
import importlib
from pathlib import Path
from typing import Any, Protocol, runtime_checkable


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
