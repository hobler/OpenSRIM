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
        import os
        import shutil
        import tempfile
        import importlib.util

        koral_dir = Path(__file__).parent / "koral"
        koral_dir_str = str(koral_dir)
        if koral_dir_str not in sys.path:
            sys.path.insert(0, koral_dir_str)

        spec = importlib.util.spec_from_file_location(
            "koral_main", koral_dir / "__main__.py"
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

        z_targets: list[int] = []
        m_targets: list[float] = []
        c_targets: list[float] = []
        for e in elements:
            if not isinstance(e, dict):
                continue
            try:
                z = int(e.get("Z", 0) or 0)
                m = float(e.get("mass_amu", 0.0) or 0.0)
            except (TypeError, ValueError):
                continue
            if z <= 0 or m <= 0:
                continue
            try:
                ratio = float(e.get("ratio", 0.0) or 0.0)
            except (TypeError, ValueError):
                ratio = 0.0
            if ratio <= 0:
                continue
            z_targets.append(z)
            m_targets.append(m)
            c_targets.append(ratio)

        if not z_targets:
            raise ValueError("No valid target elements")

        # The new KORAL main expects d_target as a list (one entry per element).
        # It represents the total atomic number density of the compound in
        # atoms/Å³; per-element weighting is done via c_target fractions.
        d_targets = [nd_A3] * len(z_targets)

        # energies in eV for KORAL
        start_eV = energies_keV[0] * 1e3
        stop_eV = energies_keV[-1] * 1e3
        nr_values = len(energies_keV)

        toml_text = _build_koral_input_toml(
            method=self._method,
            z_ion=z_ion,
            m_ion=m_ion_amu,
            z_targets=z_targets,
            m_targets=m_targets,
            d_targets=d_targets,
            c_targets=c_targets,
            start_eV=start_eV,
            stop_eV=stop_eV,
            nr_values=nr_values,
        )

        # Create a dedicated temp dir per run for KORAL's file-based I/O.
        tmp_dir = tempfile.mkdtemp(prefix="koral_run_")
        try:
            toml_path = Path(tmp_dir) / "koral_input.toml"
            toml_path.write_text(toml_text, encoding="utf-8")

            # KORAL loads SRIM stopping power tables via the relative path
            # './data/SRIM_setab/'. Temporarily chdir to the repo root so that
            # relative path resolves correctly.
            repo_root = _repo_root()
            prev_cwd = os.getcwd()
            try:
                os.chdir(str(repo_root))
                rc = KORAL(str(tmp_dir))
            finally:
                os.chdir(prev_cwd)

            if rc != 1:
                raise RuntimeError(f"KORAL returned error code {rc}")

            csv_path = Path(tmp_dir) / "koral.csv"
            if not csv_path.exists():
                raise RuntimeError("KORAL did not produce output koral.csv")

            data = np.loadtxt(str(csv_path), delimiter=";", comments="#")
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] < 7:
                raise RuntimeError(f"Unexpected koral.csv shape: {data.shape}")
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

        # koral.csv columns:
        #   0: E        [eV]
        #   1: R_p      [Å]
        #   2: sigma_x  [Å]   (longitudinal straggling)
        #   3: sigma_z  [Å]   (lateral straggling)
        #   4: S_e      [eV/Å]
        #   5: S_n      [eV/Å]
        #   6: Q_n      [eV²/Å]
        E_eV      = data[:, 0]
        R_p_A     = data[:, 1]
        sigma_x_A = data[:, 2]
        sigma_z_A = data[:, 3]
        S_e_evA   = data[:, 4]
        S_n_evA   = data[:, 5]
        Q_n_ev2A  = data[:, 6]

        energies_keV_out = [float(x) for x in (E_eV / 1e3)]

        # Convert eV/Å -> J/m
        eV_J = 1.602176634e-19
        A_to_m = 1e-10
        s_e_J = [float(x) for x in (S_e_evA * eV_J / A_to_m)]
        s_n_J = [float(x) for x in (S_n_evA * eV_J / A_to_m)]

        # Convert Å -> m
        prange_m  = [float(x) for x in (R_p_A     * A_to_m)]
        sigma_x_m = [float(x) for x in (sigma_x_A * A_to_m)]
        sigma_z_m = [float(x) for x in (sigma_z_A * A_to_m)]

        outs: dict[str, list[float]] = {
            "elec_stop":  s_e_J,
            "nucl_stop":  s_n_J,
            "nucl_strag": [float(x) for x in Q_n_ev2A],
            "prange":     prange_m,
            "long_strag": sigma_x_m,
            "lat_strag":  sigma_z_m,
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


def _build_koral_input_toml(
    *,
    method: str,
    z_ion: int,
    m_ion: float,
    z_targets: list[int],
    m_targets: list[float],
    d_targets: list[float],
    c_targets: list[float],
    start_eV: float,
    stop_eV: float,
    nr_values: int,
) -> str:
    """Build the TOML input file expected by simulators/koral/__main__.py.

    The new KORAL main function reads [settings] and [params] sections from a
    TOML file placed in its working directory.
    """

    def _fmt_float(v: float) -> str:
        return repr(float(v))

    def _fmt_int_list(xs: list[int]) -> str:
        return "[" + ", ".join(str(int(x)) for x in xs) + "]"

    def _fmt_float_list(xs: list[float]) -> str:
        return "[" + ", ".join(repr(float(x)) for x in xs) + "]"

    lines = [
        "[settings]",
        "nr_iterations = 1",
        'integration_method = "LSODA"',
        "rtol = 1e-6",
        "atol = 1e-6",
        "",
        "[params]",
        f'method = "{str(method).upper()}"',
        f"z_ion = {int(z_ion)}",
        f"m_ion = {_fmt_float(m_ion)}",
        f"z_target = {_fmt_int_list(z_targets)}",
        f"m_target = {_fmt_float_list(m_targets)}",
        f"d_target = {_fmt_float_list(d_targets)}",
        f"c_target = {_fmt_float_list(c_targets)}",
        f"start_energy = {_fmt_float(start_eV)}",
        f"stop_energy = {_fmt_float(stop_eV)}",
        f"nr_values = {int(nr_values)}",
        "",
    ]
    return "\n".join(lines)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent
