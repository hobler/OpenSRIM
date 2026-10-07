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

from simulators.opentrim import save_output
from simulators.opentrim.read_params import read_params
from simulators.opentrim.init_params import get_params
from simulators.opentrim.stats import init_stats, zero_stats
from simulators.opentrim.simulator import simulate_chunked


class _StopRequested(Exception):
    """Raised from the update callback to unwind out of the chunk loop."""


def _workdir(input_params) -> Path:
    """Resolve the run's working directory the same way process_data does."""
    workdir = input_params["simulation"]["workdir"]
    return Path(save_output.__file__).parent / workdir


def main(toml_path: str | None = None) -> None:
    input_params = read_params(toml_path)
    _nelem_ion, nelem_target, params = get_params(input_params)
    stats = init_stats(_nelem_ion, nelem_target, input_params)

    zero_stats(stats)
    nions = input_params["simulation"]["nions"]
    chunk_size = input_params["simulation"]["nions_update"]

    out_path = _workdir(input_params)
    out_path.mkdir(parents=True, exist_ok=True)

    # Advertise the live run so any OpenSRIM instance that selects this
    # directory can offer a Stop button (status: running/stopped/done/error).
    (out_path / "status").write_text("running")

    # Cooperative cancellation. The UI signals an abort by creating a
    # "stop_requested" file in the working directory. simulate_chunked already
    # writes stats and progress before invoking the callback, so raising here
    # aborts after the current chunk with valid partial results on disk.
    def _on_update(done, total, _stats):
        if (out_path / "stop_requested").exists():
            raise _StopRequested

    stopped = False
    try:
        simulate_chunked(chunk_size, nions, params, stats, input_params,
                         upd_callback=_on_update)
    except _StopRequested:
        stopped = True
    except BaseException:
        # Leave a terminal status behind so observers never wait on a
        # crashed run that still says "running".
        (out_path / "status").write_text("error")
        raise

    (out_path / "status").write_text("stopped" if stopped else "done")


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else None
    main(config_path)
