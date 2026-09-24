import tomllib
from pathlib import Path


def read_params(toml_path: str | Path = None) -> dict:
    """Read input parameters from a TOML file.

    Args:
        toml_path: Path to the TOML input file. Defaults to ``input.toml``
                   in the same directory as this script.

    Returns:
        dict: A dictionary containing the input parameters.
    """
    script_dir = Path(__file__).parent
    provided_path = Path(toml_path) if toml_path is not None else None

    if provided_path is not None and provided_path.is_file():
        toml_path = provided_path
        workdir = toml_path.parent
    else:
        print("No configuration specified; using defaults")
        toml_path = script_dir / "defaults.toml"
        workdir = "../../data/opentrim/results/default"

    with open(toml_path, "rb") as f:
        params = tomllib.load(f)

    # Convert limit lists to tuples to match the original interface.
    for dist_group in ("depth_distribution", "lateral_distribution"):
        for section in params.get("output", {}).get(dist_group, {}).values():
            if "limits" in section:
                section["limits"] = tuple(section["limits"])

    for dist_group in ("backscattered_atoms", "transmitted_atoms"):
        for section in params.get("output", {}).get(dist_group, {}).values():
            if "limits" in section:
                section["limits"] = tuple(section["limits"])

    for section in params.get("output", {}).get("distribution_2d", {}).values():
        if "nbins" in section:
            section["nbins"] = tuple(section["nbins"])
        if "limits" in section:
            section["limits"] = tuple(tuple(axis_limits) for axis_limits in section["limits"])

    params["simulation"]["workdir"] = str(Path(workdir).resolve())

    return params
