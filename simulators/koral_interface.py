from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
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
        ion = request.get("ion") or {}
        energy = request.get("energy") or {}
        target = request.get("target") or {}
        output = request.get("output") or {}

        z_ion = int(ion.get("Z", 0) or 0)
        m_ion_amu = float(ion.get("mass_amu", 0.0) or 0.0)
        energies_keV = [float(x) for x in (energy.get("energies_keV") or [])]
        energies_J = [float(x) for x in (energy.get("energies_J") or [])]

        if not energies_keV or not energies_J or len(energies_keV) != len(energies_J):
            raise ValueError("Invalid energy grid")
        if z_ion <= 0 or m_ion_amu <= 0:
            raise ValueError("Invalid ion")

        elements = target.get("elements")
        if not isinstance(elements, list) or not elements:
            raise ValueError("No target elements")

        compound_corr = float(output.get("compound_correction", target.get("compound_correction", 1.0)) or 1.0)

        # Total atomic number density (1/Å^3).
        nd_cm3 = float(target.get("number_density_atoms_cm3", 0.0) or 0.0)
        if nd_cm3 <= 0:
            raise ValueError("Invalid target atomic density")
        nd_A3 = nd_cm3 / 1.0e24

        # Atomic fractions from UI ratios.
        ratios: list[float] = []
        z_targets: list[int] = []
        m_targets_amu: list[float] = []
        for e in elements:
            if not isinstance(e, dict):
                continue
            z2 = int(e.get("Z", 0) or 0)
            if z2 <= 0:
                continue
            try:
                r = float(e.get("ratio", 0.0) or 0.0)
            except (TypeError, ValueError):
                r = 0.0
            if r <= 0:
                continue
            try:
                m2 = float(e.get("mass_amu", 0.0) or 0.0)
            except (TypeError, ValueError):
                m2 = 0.0
            if m2 <= 0:
                continue
            ratios.append(r)
            z_targets.append(z2)
            m_targets_amu.append(m2)

        if not ratios:
            raise ValueError("Target stoichiometry is zero")

        total_ratio = float(sum(ratios))
        fracs = [r / total_ratio for r in ratios]

        # Stopping powers (J/m)
        s_e_J_per_m = _electronic_stopping_J_per_m(
            z_ion=z_ion,
            z_targets=z_targets,
            fracs=fracs,
            nd_A3=nd_A3,
            energies_J=np.asarray(energies_J, dtype=float),
            compound_correction=compound_corr,
        )

        s_n_J_per_m = _nuclear_stopping_J_per_m(
            method=self._method,
            z_ion=z_ion,
            m_ion_amu=m_ion_amu,
            z_targets=z_targets,
            m_targets_amu=m_targets_amu,
            fracs=fracs,
            nd_A3=nd_A3,
            energies_J=np.asarray(energies_J, dtype=float),
        )

        # CSDA projected range (m): R(E) = ∫ dE / (Se+Sn)
        s_tot = np.maximum(s_e_J_per_m + s_n_J_per_m, 1e-30)
        E = np.asarray(energies_J, dtype=float)

        # Ensure monotonic increasing energies for integration.
        if np.any(np.diff(E) < 0):
            order = np.argsort(E)
            E = E[order]
            energies_keV_sorted = [energies_keV[i] for i in order]
            s_e_J_per_m = s_e_J_per_m[order]
            s_n_J_per_m = s_n_J_per_m[order]
            s_tot = s_tot[order]
        else:
            energies_keV_sorted = energies_keV

        inv_s = 1.0 / s_tot
        dE = np.diff(E)
        avg = (inv_s[:-1] + inv_s[1:]) * 0.5
        cum = np.concatenate(([0.0], np.cumsum(dE * avg)))

        zeros = np.zeros_like(cum)

        outs: dict[str, list[float]] = {
            "elec_stop": [float(x) for x in s_e_J_per_m],
            "nucl_stop": [float(x) for x in s_n_J_per_m],
            "prange": [float(x) for x in cum],
            "long_strag": [float(x) for x in zeros],
            "lat_strag": [float(x) for x in zeros],
        }

        requested = output.get("requested")
        if isinstance(requested, list) and requested:
            outs = {k: v for k, v in outs.items() if k in {str(x) for x in requested}}

        return {
            "model_id": self._model_ref,
            "energies_keV": list(energies_keV_sorted),
            "outputs": outs,
        }


def _repo_root() -> Path:
    # This file lives in <repo>/simulators/; root is one level above.
    return Path(__file__).resolve().parent.parent


def _data_dir() -> Path:
    return _repo_root() / "data"


@lru_cache(maxsize=128)
def _load_srim_setab(ion_z: int) -> tuple[np.ndarray, np.ndarray]:
    path = _data_dir() / "SRIM_setab" / f"SRIM2013-{int(ion_z):02d}.dat"
    if not path.exists():
        raise FileNotFoundError(f"Missing SRIM table: {path}")
    data = np.loadtxt(path, skiprows=6, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Invalid SRIM table format: {path}")
    energies_eV = data[:, 0]
    stop_eV_A2 = data[:, 1:]
    return energies_eV, stop_eV_A2


def _electronic_stopping_J_per_m(
    *,
    z_ion: int,
    z_targets: list[int],
    fracs: list[float],
    nd_A3: float,
    energies_J: np.ndarray,
    compound_correction: float,
) -> np.ndarray:
    """Electronic stopping (J/m) from SRIM-2013 stopping cross sections."""

    energies_eV_tab, stop_eV_A2_tab = _load_srim_setab(int(z_ion))

    eV_J = 1.602176634e-19

    energies_eV = energies_J / eV_J

    mix_cs = np.zeros_like(energies_eV, dtype=float)

    for z2, fi in zip(z_targets, fracs):
        if z2 <= 0 or fi <= 0:
            continue
        col = int(z2) - 1
        if col < 0 or col >= stop_eV_A2_tab.shape[1]:
            continue

        y = stop_eV_A2_tab[:, col]
        mix_cs += fi * _interp_positive(energies_eV_tab, y, energies_eV)

    linear_eV_per_A = mix_cs * float(nd_A3) * float(compound_correction)

    return linear_eV_per_A * (eV_J / 1.0e-10)


@lru_cache(maxsize=1)
def _nlh_params() -> dict[tuple[int, int], tuple[np.ndarray, np.ndarray]]:
    """Return {(Z1,Z2): (sn_params[a,b,c,d], qn_params[a,b,c,d])}."""

    sn_path = _data_dir() / "NLH" / "sn_fit_params.txt"
    qn_path = _data_dir() / "NLH" / "qn_fit_params.txt"

    def _read(path: Path) -> dict[tuple[int, int], np.ndarray]:
        out: dict[tuple[int, int], np.ndarray] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 6:
                continue
            try:
                z1 = int(parts[0])
                z2 = int(parts[1])
                a, b, c, d = (float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5]))
            except Exception:
                continue
            out[(z1, z2)] = np.asarray([a, b, c, d], dtype=float)
        return out

    sn = _read(sn_path) if sn_path.exists() else {}
    qn = _read(qn_path) if qn_path.exists() else {}

    merged: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
    for key, sn_params in sn.items():
        qn_params = qn.get(key)
        if qn_params is None:
            qn_params = np.asarray([0.0, 0.0, 0.0, 0.0], dtype=float)
        merged[key] = (sn_params, qn_params)
    return merged


def _nuclear_stopping_J_per_m(
    *,
    method: str,
    z_ion: int,
    m_ion_amu: float,
    z_targets: list[int],
    m_targets_amu: list[float],
    fracs: list[float],
    nd_A3: float,
    energies_J: np.ndarray,
) -> np.ndarray:
    """Nuclear stopping (J/m) using ZBL or NLH universal forms."""

    mid = str(method).strip().upper()
    if mid not in {"ZBL", "NLH"}:
        raise ValueError(f"Unsupported nuclear stopping method '{method}'")

    total = np.zeros_like(energies_J, dtype=float)

    eV_J = 1.602176634e-19

    energies_keV = energies_J / (eV_J * 1.0e3)

    for z2, m2, fi in zip(z_targets, m_targets_amu, fracs):
        if fi <= 0:
            continue
        nd_i = float(nd_A3) * float(fi)

        z1 = float(z_ion)
        z2f = float(z2)
        m1 = float(m_ion_amu)
        m2f = float(m2)
        denom = (z1 * z2f * ((z1 ** 0.23) + (z2f ** 0.23)))
        if denom <= 0:
            continue
        eps = 32.53 * (m2f / (m1 + m2f)) * energies_keV / denom

        if mid == "ZBL":
            sn_red = np.log(1.0 + 1.1383 * eps) / (2.0 * (eps + 0.01321 * (eps**0.21226) + 0.19593 * (eps**0.5)))
        else:
            params = _nlh_params().get((int(z_ion), int(z2)))
            sn_params = params[0] if params is not None else np.asarray([1.0, 0.0, 0.0, 0.0], dtype=float)
            a_p, b_p, c_p, d_p = [float(x) for x in np.asarray(sn_params, dtype=float).tolist()]
            sn_red = np.log(1.0 + a_p * eps) / (2.0 * (eps + b_p * (eps**c_p) + d_p * (eps**0.5)))

        pref = 8.462 * (z1 * z2f * m1) / ((m1 + m2f) * ((z1 ** 0.23) + (z2f ** 0.23)))

        s_cs_eV_per_1e15 = pref * sn_red
        s_cs_eV_A2 = s_cs_eV_per_1e15 * 10.0

        s_lin_eV_per_A = s_cs_eV_A2 * nd_i

        total += np.asarray(s_lin_eV_per_A, dtype=float) * (eV_J / 1.0e-10)

    return total


def _interp_positive(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    """Interpolate y(x) at xq; uses log-log when y>0 and x>0."""

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xq = np.asarray(xq, dtype=float)

    ok = (x > 0) & (y > 0)
    if np.count_nonzero(ok) >= 2 and np.all(xq > 0):
        lx = np.log(x[ok])
        ly = np.log(y[ok])
        lyq = np.interp(np.log(xq), lx, ly, left=ly[0], right=ly[-1])
        return np.exp(lyq)

    return np.interp(xq, x, y, left=float(y[0]), right=float(y[-1]))
