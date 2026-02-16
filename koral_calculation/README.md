# KORAL Calculation Models

This folder contains pluggable calculation models used by the KORAL page.

## Discovery

A model is a subfolder:

- `koral_calculation/<model_id>/`

The KORAL UI discovers models by listing these folders.

## UI parameter limits (optional)

The KORAL UI can apply model-specific **min/max ranges** (and a few optional
display hints) for its input widgets.

If a model folder contains:

- `koral_calculation/<model_id>/ui_params.toml`

then the UI will load it and apply the constraints to the corresponding
widgets. If multiple models are selected at the same time, the UI uses the
**intersection** of all selected models' ranges (max of mins / min of maxes). If
the intersection is empty, the UI falls back to its built-in defaults.

### File format (`ui_params.toml`)

Top-level object:

- `format_version` (int, required): currently `1`
- `parameters` (object, required): mapping of parameter-id → spec

Parameter spec:

- `min` (number, required)
- `max` (number, required)
- `decimals` (int, optional): for `QDoubleSpinBox`
- `step` (number, optional): for `QSpinBox`/`QDoubleSpinBox`
- `default` (number, optional)
- `label` (string, optional)
- `unit` (string, optional)
- `description` (string, optional)

Supported parameter IDs (current KORAL page):

- `ion_mass_amu`
- `energy_min_keV`
- `energy_max_keV`
- `compound_correction`

Example:

```toml
format_version = 1

[parameters.energy_min_keV]
min = 0.0
max = 1000000.0
decimals = 2

[parameters.compound_correction]
min = 0.0
max = 10.0
step = 0.01
default = 1.0
```

### Primary model

If a model folder contains a file named `.primarymodel`, it is treated as a default/preferred model (preselected in the UI).

## Required Python module

Each model folder must contain:

- `koral_calculation/<model_id>/model.py`

The loader supports any one of these exports:

- `MODEL`: an instance with a `run(request: dict) -> dict` method
- `create_model() -> object`: returns an instance with `run(request) -> dict`
- `Model`: a class that can be instantiated with no arguments and provides `run(request) -> dict`

See loader implementation in `koral_calculation/interface.py`.

## Request contract

The KORAL UI builds a `request` dictionary and passes it to `run()`.

Your model should:

- Treat unknown keys as optional
- Avoid mutating the request
- Validate required inputs and either raise a `ValueError` with a clear message or return zeros + diagnostics

### Common fields

```python
request = {
  "ion": {
    "Z": int,
    "mass_amu": float,
    # ... optional
  },
  "energy": {
    # Preferred: explicit energy grid
    "energies_keV": [float, ...],

    # Also provided for SI convenience
    "energies_J": [float, ...],
  },
  "target": {
    "elements": [
      {"Z": int, "mass_amu": float, "ratio": float},
      # ratio is stoichiometric; the model should normalize to atomic fractions if needed
    ],

    # Density information (both may be present)
    "number_density_atoms_cm3": float,
    "number_density_atoms_m3": float,
    "density_g_cm3": float,
    "density_kg_m3": float,

    # Optional tuning inputs
    "compound_correction": float,
    "elec_scale": float,
    "nucl_scale": float,
  },
  "output": {
    "requested": ["elec_stop", "nucl_stop", "range_csda", ...],
    "units": {
      # Display units requested by UI (models may ignore this; UI converts for display)
      "elec_stop": "keV/µm",
      "prange": "Ång",
    }
  }
}
```

## Result contract

`run()` must return a dictionary:

```python
result = {
  "model_id": "your_model_id",
  "energies_keV": [float, ...],
  "outputs": {
    "output_id": [float, ...],
    # each list must match the length/order of energies_keV
  },

  # Optional but recommended
  "units": {
    "output_id": "J/m" | "m" | "...",
  },
  "diagnostics": {"any": "json-serializable"},
}
```

### Standardized base units (important)

Models must compute and return **standardized SI base units**. The UI is responsible for converting to display/export units.

Recommended base units:

- Stopping powers (`elec_stop`, `nucl_stop`): **J/m**
- Length-like outputs (`range_csda`, `prange`, `long_strag`, `lat_strag`): **m**

## Output IDs

The UI only displays outputs whose IDs match the selection list on the KORAL page.

The built-in SRIM-like model uses (and the UI expects) these SRIM-style IDs:

- `elec_stop` (electronic stopping)
- `nucl_stop` (nuclear stopping)
- `prange` (range)
- `range_csda` (range)
- `long_strag` (longitudinal straggling)
- `lat_strag` (lateral straggling)

You can add more outputs, but they must be selectable in the UI to be shown.
