# OpenSRIM

OpenSRIM is a desktop application for simulating ion implantation into
solids — projected range, straggling, stopping powers, sputtering and
related quantities — in the spirit of SRIM/TRIM. It provides a PyQt6
graphical interface over two simulation backends:

- **KORAL** — a fast analytical solver for stopping powers, range and
  straggling as a function of ion energy.
- **OpenTRIM** — a Monte Carlo binary-collision simulator for full ion/recoil
  cascades, damage statistics and 1-D/2-D distributions.

A **Single Plot** view lets you combine, style and post-process (e.g.
Gaussian convolution) results from either backend for closer analysis or
publication-quality figures.

> **This project is under active development.** Behavior, file formats and
> the UI are still changing, and not everything is polished or fully
> validated yet. If something looks wrong, confusing or missing, please
> [open an issue](https://github.com/hobler/OpenSRIM/issues) — feedback at
> this stage is genuinely useful and appreciated.

## Installation

OpenSRIM targets Python 3.11+ and Linux (Ubuntu); other platforms may work
but are untested. This early on, installation is aimed at technically
comfortable users (venv + pip), not a packaged one-click installer.

```bash
git clone git@github.com:hobler/OpenSRIM.git
cd OpenSRIM
./setup.sh          # installs system (Qt/X11) packages + creates ./venv
source venv/bin/activate
python3 main_window.py
```

`setup.sh` installs the non-Python system libraries PyQt6 needs (via `apt`)
and then creates a local virtual environment from `requirements.txt`. If you
already have a suitable Python environment, you can skip `setup.sh` and just
run:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python3 main_window.py
```

### Notes on running OpenTRIM

OpenTRIM uses Numba for just-in-time compilation. As long as you don't modify
jit-compiled code, set `ENABLE_CACHING = getattr(sys, "frozen", True)` in `simulators/opentrim/config.py` to avoid recompilation at every run.

OpenTRIM may also be run from the command line using a hand-writted or 
otherwise generated `input.toml` file:

```bash
source /path/to/venv/bin/activate
python3 /path/to/simulators/opentrim/__main__.py /path/to/input.toml
```

where `/path/to` has to be replaced by the respective directories of the files.
The output files will be generated in the directory of `input.toml`.

## License

MIT — see [LICENSE](LICENSE).
