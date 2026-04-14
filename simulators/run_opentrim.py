"""Headless OpenTRIM runner for UI integration.

Runs the simulation and writes results to files (histograms, moments,
progress) without opening matplotlib windows or printing to stdout.

Usage:
    python -m simulators.run_opentrim <input.toml>
"""
import sys
from pathlib import Path

if __package__ is None:
    project_root = str(Path(__file__).resolve().parent.parent)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    __package__ = "simulators"

from simulators.opentrim.read_params import read_params
from simulators.opentrim.init_params import get_params
from simulators.opentrim.stats import init_stats, zero_stats
from simulators.opentrim.simulator import simulate_chunked


def main(toml_path: str | None = None) -> None:
    input_params = read_params(toml_path)
    _nelem_ion, nelem_target, params = get_params(input_params)
    stats = init_stats(_nelem_ion, nelem_target, input_params)

    zero_stats(stats)
    nions = input_params["simulation"]["nions"]
    chunk_size = input_params["simulation"]["nions_update"]
    simulate_chunked(chunk_size, nions, params, stats, input_params)


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(config_path)
