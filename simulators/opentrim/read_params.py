import tomllib
from pathlib import Path


def read_params(toml_path = None) -> dict:
    """Read input parameters from a TOML file.

    Args:
        toml_path: Path to the TOML input file. Defaults to ``input.toml``
                   in the same directory as this script.

    Returns:
        dict: A dictionary containing the input parameters.
    """
    if toml_path is None:
        toml_path = "defaults.toml"
    toml_path = Path(toml_path)
    toml_path_fallback = Path(__file__).parent / Path(toml_path)
    if not toml_path.is_file():
        toml_path = toml_path_fallback
    if not toml_path.is_file():
        raise FileNotFoundError("The provided config is not a file or doesn't exist")

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

    return params
