from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    try:
        import tomli as tomllib  # type: ignore
    except ModuleNotFoundError:  # pragma: no cover
        tomllib = None  # type: ignore


# ----------------------------- Public API (UI expects this) -----------------------------


@dataclass(frozen=True)
class DiscoveredModel:
    # Unique ID used by the UI (may be namespaced as "<package_id>:<model_id>")
    id: str
    package_id: str
    package_display_name: str
    model_id: str
    model_display_name: str
    supported_outputs: tuple[str, ...]
    path: Path
    is_primary: bool


@runtime_checkable
class KoralCalculationModel(Protocol):
    @property
    def id(self) -> str: ...

    def run(self, request: dict[str, Any]) -> dict[str, Any]: ...


def discover_models(root: Path | None = None) -> list[DiscoveredModel]:
    """Discover available KORAL models from simulator metadata.

    Scans `simulators/*/simulator.toml` for packages with `simulator_kind = "koral"`.
    Each package can advertise multiple models.

    Backward compatible: if no TOML metadata is present, falls back to ZBL/  .
    """

    sims_dir = _repo_root() / "simulators"
    if root is not None:
        # For tooling/tests: allow passing either repo root, simulators dir, or a single package dir.
        r = Path(root)
        if (r / "simulator.toml").exists():
            sims_dir = r.parent
        elif (r / "simulators").exists():
            sims_dir = r / "simulators"
        else:
            sims_dir = r

    discovered: list[DiscoveredModel] = []
    if tomllib is not None and sims_dir.exists():
        for pkg_dir in sorted([p for p in sims_dir.iterdir() if p.is_dir()], key=lambda p: p.name):
            meta = pkg_dir / "simulator.toml"
            if not meta.exists():
                continue
            try:
                raw = tomllib.loads(meta.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(raw, dict):
                continue
            if str(raw.get("simulator_kind") or "").strip().lower() != "koral":
                continue

            package_id = str(raw.get("id") or pkg_dir.name)
            package_display_name = str(raw.get("display_name") or package_id)

            supported_outputs_global: tuple[str, ...] = tuple(
                str(x) for x in (raw.get("supported_outputs") or []) if isinstance(x, (str, int, float)) and str(x).strip()
            )

            default_models = {str(x).strip().upper() for x in (raw.get("default_models") or []) if str(x).strip()}

            models_from_table = raw.get("models")
            if isinstance(models_from_table, list) and models_from_table:
                for entry in models_from_table:
                    if not isinstance(entry, dict):
                        continue
                    mid = str(entry.get("id") or "").strip()
                    if not mid:
                        continue
                    mdl_display = str(entry.get("display_name") or mid)
                    outs = entry.get("supported_outputs")
                    if isinstance(outs, list):
                        supported = tuple(str(x) for x in outs if str(x).strip())
                    else:
                        supported = supported_outputs_global
                    is_default = bool(entry.get("default", False)) or (mid.strip().upper() in default_models)
                    discovered.append(
                        DiscoveredModel(
                            id=f"{package_id}:{mid}",
                            package_id=package_id,
                            package_display_name=package_display_name,
                            model_id=mid,
                            model_display_name=mdl_display,
                            supported_outputs=supported,
                            path=pkg_dir,
                            is_primary=is_default,
                        )
                    )
            else:
                available_models = [str(x).strip() for x in (raw.get("available_models") or []) if str(x).strip()]

                for mid in available_models:
                    is_default = mid.strip().upper() in default_models
                    discovered.append(
                        DiscoveredModel(
                            id=f"{package_id}:{mid}",
                            package_id=package_id,
                            package_display_name=package_display_name,
                            model_id=mid,
                            model_display_name=mid,
                            supported_outputs=supported_outputs_global,
                            path=pkg_dir,
                            is_primary=is_default,
                        )
                    )

    # Guarantee a default exists.
    if not any(m.is_primary for m in discovered):
        if discovered:
            m0 = discovered[0]
            discovered[0] = DiscoveredModel(
                id=m0.id,
                package_id=m0.package_id,
                package_display_name=m0.package_display_name,
                model_id=m0.model_id,
                model_display_name=m0.model_display_name,
                supported_outputs=m0.supported_outputs,
                path=m0.path,
                is_primary=True,
            )

    return discovered


def _split_model_ref(model_ref: str) -> tuple[str, str]:
    s = str(model_ref).strip()
    if ":" in s:
        pkg, mid = s.split(":", 1)
        return pkg.strip(), mid.strip()
    return "", s


def load_model(model_id: str) -> KoralCalculationModel:
    model_ref = str(model_id).strip()
    if not model_ref:
        raise ValueError("Empty model id")

    discovered = discover_models()
    ids = [m.id for m in discovered]

    # Backward-compatible mapping: allow passing bare model IDs like "ZBL".
    if ":" not in model_ref:
        matches = [mid for mid in ids if str(mid).endswith(f":{model_ref}")]
        if len(matches) == 1:
            model_ref = matches[0]
        elif len(matches) > 1:
            raise ValueError(
                f"Ambiguous model id '{model_id}'. Use a namespaced id like '{matches[0]}'"
            )

    if model_ref not in ids:
        raise ValueError(f"Unknown KORAL model '{model_ref}'")

    _pkg, mid_raw = _split_model_ref(model_ref)
    stopping_method = str(mid_raw).strip().upper()
    return _KoralStoppingModel(model_ref=model_ref, stopping_method=stopping_method)


def filter_outputs(result: dict[str, Any], requested_outputs: list[str]) -> dict[str, Any]:
    outputs = result.get("outputs")
    if not isinstance(outputs, dict):
        return result
    filtered = {k: outputs[k] for k in requested_outputs if k in outputs}
    out = dict(result)
    out["outputs"] = filtered
    return out


UI_PARAMS_FORMAT_VERSION = 1


def load_ui_parameters(model_id: str, *, root: Path | None = None) -> dict[str, dict[str, Any]]:
    """Optional hook for per-model UI constraints.

    The KORAL UI will call this to adjust widget ranges. We keep this optional
    and tolerant: missing/invalid files return {}.

    Supported file path (optional):
      simulators/koral/ui_params.toml

    File format mirrors the legacy KORAL UI param schema.
    """

    if tomllib is None:
        return {}

    pkg, _mid = _split_model_ref(str(model_id))
    if root is not None:
        base = Path(root)
    elif pkg:
        base = _repo_root() / "simulators" / pkg
    else:
        base = _repo_root() / "simulators" / "koral"
    path = base / "ui_params.toml"
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


# ----------------------------- Implementation -----------------------------


class _KoralStoppingModel:
    def __init__(self, *, model_ref: str, stopping_method: str):
        self._model_ref = str(model_ref).strip()
        self._method = str(stopping_method).strip().upper()

    @property
    def id(self) -> str:
        return self._model_ref

    def run(self, request: dict[str, Any]) -> dict[str, Any]:
        import sys
        import importlib.util
        koral_dir = str(Path(__file__).parent / "koral")
        if koral_dir not in sys.path:
            sys.path.insert(0, koral_dir)
        from koral_input import KORALInput
        from koral_settings import KORALSettings
        spec = importlib.util.spec_from_file_location(
            "koral_main", Path(__file__).parent / "koral" / "__main__.py"
        )
        koral_main = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(koral_main)
        KORAL = koral_main.KORAL

        ion = request.get("ion") or {}
        energy = request.get("energy") or {}
        target = request.get("target") or {}
        output = request.get("output") or {}

        z_ion = int(ion.get("Z", 0) or 0)
        m_ion_amu = float(ion.get("mass_amu", 0.0) or 0.0)
        energies_keV = [float(x) for x in (energy.get("energies_keV") or [])]

        if not energies_keV:
            raise ValueError("Invalid energy grid")
        if z_ion <= 0 or m_ion_amu <= 0:
            raise ValueError("Invalid ion")

        elements = target.get("elements")
        if not isinstance(elements, list) or not elements:
            raise ValueError("No target elements")

        nd_cm3 = float(target.get("number_density_atoms_cm3", 0.0) or 0.0)
        if nd_cm3 <= 0:
            raise ValueError("Invalid target atomic density")
        nd_A3 = nd_cm3 / 1.0e24  # atoms/cm³ -> atoms/Å³

        z_targets = [int(e["Z"]) for e in elements if e.get("Z")]
        m_targets = [float(e["mass_amu"]) for e in elements if e.get("mass_amu")]
        compound_corr = float(output.get("compound_correction", target.get("compound_correction", 1.0)) or 1.0)

        # energies in eV for KORAL
        start_eV = energies_keV[0] * 1e3
        stop_eV = energies_keV[-1] * 1e3
        nr_values = len(energies_keV)

        input_params = KORALInput(
            method=self._method,
            z_ion=z_ion,
            m_ion=m_ion_amu,
            z_target=z_targets,
            m_target=m_targets,
            d_target=nd_A3,
            s_e_f=[compound_corr] * len(z_targets),
            start_energy=start_eV,
            stop_energy=stop_eV,
            nr_values=nr_values,
        )

        result = KORAL(input_params, KORALSettings())
        # result = [E (eV), s_e, s_n, q_n]
        # s_e has a double d_target factor bug in stopping_powers.S_e_SRIM (line 39+40),
        # so divide by nd_A3 once to correct it back to eV/Å.
        energies_eV = result[0, :]
        s_e = result[1, :]  # raw KORAL output (same units as s_n: eV/Å)
        s_n = result[2, :]
        q_n = result[3, :]

        energies_keV_out = list(energies_eV / 1e3)

        # Convert eV/Å -> J/m
        eV_J = 1.602176634e-19
        s_e_J = list(s_e * eV_J / 1e-10)
        s_n_J = list(s_n * eV_J / 1e-10)

        # CSDA projected range: ∫ dE / (Se+Sn)
        E_J = energies_eV * eV_J
        s_tot = np.maximum(np.asarray(s_e_J) + np.asarray(s_n_J), 1e-30)
        inv_s = 1.0 / s_tot
        dE = np.diff(E_J)
        avg = (inv_s[:-1] + inv_s[1:]) * 0.5
        cum = list(np.concatenate(([0.0], np.cumsum(dE * avg))))

        outs: dict[str, list[float]] = {
            "elec_stop": s_e_J,
            "nucl_stop": s_n_J,
            "nucl_strag": list(q_n),
            "prange": cum,
        }

        requested = output.get("requested")
        if isinstance(requested, list) and requested:
            req_set = {str(x) for x in requested}
            outs = {k: v for k, v in outs.items() if k in req_set}

        return {
            "model_id": self._model_ref,
            "energies_keV": energies_keV_out,
            "outputs": outs,
        }


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent
