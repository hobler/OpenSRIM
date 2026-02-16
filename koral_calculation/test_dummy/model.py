# koral_calculation/test_dummy/model.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping

try:
    from ui.logging import log as emit_log
except ModuleNotFoundError:  # pragma: no cover
    try:
        from OpenSRIM.ui.logging import log as emit_log  # type: ignore
    except ModuleNotFoundError:  # pragma: no cover
        def emit_log(message: str) -> None:  # type: ignore
            return


# Default-Outputs, falls die UI keine Liste liefert
DEFAULT_OUTPUTS = [
    "elec_stop",
    "nucl_stop",
    "range_csda",
    "prange",
    "long_strag",
    "lat_strag",
]

# Zuordnung zu empfohlenen SI-Basiseinheiten
BASE_UNITS = {
    "elec_stop": "J/m",
    "nucl_stop": "J/m",
    "range_csda": "m",
    "prange": "m",
    "long_strag": "m",
    "lat_strag": "m",
}


def _get_nested(d: Mapping[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = d
    for k in keys:
        if not isinstance(cur, Mapping) or k not in cur:
            return default
        cur = cur[k]
    return cur


def _ensure_float_list(x: Any) -> List[float]:
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        out: List[float] = []
        for v in x:
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                out.append(float("nan"))
        return out
    # single value
    try:
        return [float(x)]
    except (TypeError, ValueError):
        return [float("nan")]


@dataclass
class DummyModel:
    model_id: str = "test_dummy"

    def run(self, request: Dict[str, Any]) -> Dict[str, Any]:
        try:
            ion_sym = _get_nested(request, "ion", "symbol", default="") or ""
            ion_z = _get_nested(request, "ion", "Z", default=None)
            emit_log(f"KORAL[{self.model_id}]: run started (ion={ion_sym}, Z={ion_z})")
        except Exception:
            emit_log(f"KORAL[{self.model_id}]: run started")

        # 1) Energie-Grid bestimmen (bevorzugt energies_keV)
        energies_keV = _ensure_float_list(_get_nested(request, "energy", "energies_keV"))
        if not energies_keV:
            # Fallback: energies_J -> in keV umrechnen, wenn vorhanden
            energies_J = _ensure_float_list(_get_nested(request, "energy", "energies_J"))
            if energies_J:
                e_charge = 1.602176634e-19  # C
                # 1 eV = e Joule; 1 keV = 1e3 eV
                energies_keV = [(Ej / e_charge) / 1e3 for Ej in energies_J]
            else:
                # Letzter Fallback: kleines Default-Grid
                energies_keV = [10.0, 100.0, 1000.0]

        n = len(energies_keV)
        emit_log(f"KORAL[{self.model_id}]: energy grid prepared (n={n})")

        # 2) gewünschte Outputs lesen
        requested = _get_nested(request, "output", "requested", default=None)
        if not isinstance(requested, list) or not requested:
            requested = list(DEFAULT_OUTPUTS)
        emit_log(f"KORAL[{self.model_id}]: requested outputs = {', '.join(map(str, requested))}")

        # 3) Dummy-Berechnungen (deterministisch, glatt, immer Länge n)
        # Ziel: plausible Shapes + positive Größen, aber bewusst "fake".
        outputs: Dict[str, List[float]] = {}
        units: Dict[str, str] = {}

        for out_id in requested:
            out_id_str = str(out_id)

            if out_id_str in ("elec_stop", "nucl_stop"):
                # Beispiel: stopping ~ a * sqrt(E) in J/m
                # Skalen so wählen, dass Zahlen nicht absurd groß/klein sind.
                scale = 1e-12 if out_id_str == "elec_stop" else 3e-13
                vals = [scale * (max(e, 0.0) ** 0.5) for e in energies_keV]  # J/m (Dummy)
                outputs[out_id_str] = vals
                units[out_id_str] = BASE_UNITS.get(out_id_str, "J/m")

            else:
                # Beispiel: Längen ~ b * E^(1.2) in m (Dummy)
                scale = 1e-9
                vals = [scale * (max(e, 0.0) ** 1.2) for e in energies_keV]  # m (Dummy)
                outputs[out_id_str] = vals
                units[out_id_str] = BASE_UNITS.get(out_id_str, "m")

        # 4) Diagnostics: Request nicht mutieren, aber Meta-Infos liefern
        ion_Z = _get_nested(request, "ion", "Z", default=None)
        tgt_elems = _get_nested(request, "target", "elements", default=[])
        diagnostics = {
            "note": "Dummy/Testmodell: liefert synthetische Werte in SI-Basiseinheiten.",
            "requested_outputs": requested,
            "n_energies": n,
            "ion_Z": ion_Z,
            "target_elements_count": len(tgt_elems) if isinstance(tgt_elems, list) else None,
        }

        emit_log(f"KORAL[{self.model_id}]: run finished")

        return {
            "model_id": self.model_id,
            "energies_keV": energies_keV,
            "outputs": outputs,
            "units": units,
            "diagnostics": diagnostics,
        }


# Loader-Export: MODEL Instanz
MODEL = DummyModel()
