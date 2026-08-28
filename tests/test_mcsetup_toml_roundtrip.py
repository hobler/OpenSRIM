from __future__ import annotations

import os
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent
for candidate in (str(ROOT), str(PARENT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from PyQt6.QtWidgets import QApplication

from OpenSRIM.simulators.opentrim.read_params import read_params
from OpenSRIM.state import AppState
from OpenSRIM.ui.pages.mcsetup_page import MCSetupPage


ROUNDTRIP_TOML = """
[simulation]
follow_recoils = false
nions = 321
nions_update = 17
rng_seed = 6789
nthreads = 2
workdir = "legacy-results"

[beam]
symbol = "B"
name = "Boron"
Z = 5
M = 11.009
energy = 12.5
tilt = 7.5

[[layer]]
name = "Top"
width = 1234.0
density = 0.05
compound_correction = 0.77
gas = false

    [[layer.element]]
    symbol = "Si"
    name = "Silicon"
    Z = 14
    M = 28.086
    stoichiometry = 0.6
    displacement_energy = 19.0
    lattice_binding_energy = 8.25
    surface_binding_energy = 5.5

[[layer]]
name = "Bottom"
width = 4321.0
density = 0.06
compound_correction = 1.23
gas = false

    [[layer.element]]
    symbol = "O"
    name = "Oxygen"
    Z = 8
    M = 15.999
    stoichiometry = 0.4
    displacement_energy = 21.0
    lattice_binding_energy = 9.75
    surface_binding_energy = 6.5

[models]
potential = "NLHlin"
electronic_stopping = "Lindhard"

[models.scattering_integrals]
algorithm = "magic"
n_absc = 9

[models.lindhard_correction]
"B->Si" = 1.75

[cascade]
pmax_min = 0.1
pmax_max = 5.2
psi_min = 6.3
de_min = 16.4
psi_min_surface = 7.5
de_min_surface = 17.6
replacement_collisions = true

[output]

[output.trajectories]
start = true
collisions = true
stopped = true
backscattered = false
transmitted = false

[output.depth_distribution.ion_recoils]
score = true
nbins = 77
limits = [1.0, 5555.0]

[output.depth_distribution.nuclear_energy_deposition]
score = true
nbins = 77
limits = [1.0, 5555.0]

[output.depth_distribution.electronic_energy_deposition]
score = true
nbins = 77
limits = [1.0, 5555.0]

[output.lateral_distribution.ion_recoils]
score = true
nbins = 77
limits = [-22.0, 33.0]

[output.lateral_distribution.nuclear_energy_deposition]
score = true
nbins = 77
limits = [-22.0, 33.0]

[output.lateral_distribution.electronic_energy_deposition]
score = true
nbins = 77
limits = [-22.0, 33.0]

[output.distribution_2d.ion_recoils]
score = true
nbins = [77, 77]
limits = [[1.0, 5555.0], [-22.0, 33.0]]

[output.distribution_2d.nuclear_energy_deposition]
score = true
nbins = [77, 77]
limits = [[1.0, 5555.0], [-22.0, 33.0]]

[output.distribution_2d.electronic_energy_deposition]
score = true
nbins = [77, 77]
limits = [[1.0, 5555.0], [-22.0, 33.0]]

[output.backscattered_atoms.energy]
score = true
nbins = 77
limits = [0.0, 12500.0]

[output.backscattered_atoms.angle]
score = true
nbins = 77
limits = [-80.0, 80.0]

[output.transmitted_atoms.energy]
score = true
nbins = 77
limits = [0.0, 12500.0]

[output.transmitted_atoms.angle]
score = true
nbins = 77
limits = [-80.0, 80.0]
"""


class MCSetupTomlRoundtripTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])

    def test_load_save_preserves_nondefault_toml_fields(self) -> None:
        page = MCSetupPage(AppState())

        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.toml"
            input_path.write_text(ROUNDTRIP_TOML, encoding="utf-8")

            params = read_params(input_path)
            page.apply_simulation_config(page._raw_toml_to_payload(params))
            saved = tomllib.loads(page._build_input_toml(tmp))

        self.assertEqual(saved["simulation"]["follow_recoils"], False)
        self.assertEqual(saved["simulation"]["nions"], 321)
        self.assertEqual(saved["simulation"]["nions_update"], 17)
        self.assertEqual(saved["simulation"]["rng_seed"], 6789)
        self.assertEqual(saved["simulation"]["nthreads"], 2)

        self.assertEqual(saved["beam"]["symbol"], "B")
        self.assertEqual(saved["beam"]["Z"], 5)
        self.assertAlmostEqual(saved["beam"]["M"], 11.009)
        self.assertAlmostEqual(saved["beam"]["energy"], 12.5)
        self.assertAlmostEqual(saved["beam"]["tilt"], 7.5)

        self.assertEqual([layer["name"] for layer in saved["layer"]], ["Top", "Bottom"])
        self.assertAlmostEqual(saved["layer"][0]["density"], 0.05)
        self.assertAlmostEqual(saved["layer"][1]["density"], 0.06)
        self.assertAlmostEqual(saved["layer"][0]["compound_correction"], 0.77)
        self.assertAlmostEqual(saved["layer"][1]["compound_correction"], 1.23)

        elements = [layer["element"][0] for layer in saved["layer"]]
        self.assertEqual([element["symbol"] for element in elements], ["Si", "O"])
        self.assertEqual([element["Z"] for element in elements], [14, 8])
        self.assertAlmostEqual(elements[0]["M"], 28.086)
        self.assertAlmostEqual(elements[1]["M"], 15.999)
        self.assertAlmostEqual(elements[0]["displacement_energy"], 19.0)
        self.assertAlmostEqual(elements[1]["displacement_energy"], 21.0)
        self.assertAlmostEqual(elements[0]["lattice_binding_energy"], 8.25)
        self.assertAlmostEqual(elements[1]["lattice_binding_energy"], 9.75)
        self.assertAlmostEqual(elements[0]["surface_binding_energy"], 5.5)
        self.assertAlmostEqual(elements[1]["surface_binding_energy"], 6.5)

        self.assertEqual(saved["models"]["potential"], "NLHlin")
        self.assertEqual(saved["models"]["electronic_stopping"], "Lindhard")
        self.assertEqual(saved["models"]["scattering_integrals"]["algorithm"], "magic")
        self.assertEqual(saved["models"]["scattering_integrals"]["n_absc"], 9)
        self.assertAlmostEqual(saved["models"]["lindhard_correction"]["B->Si"], 1.75)

        self.assertAlmostEqual(saved["cascade"]["pmax_min"], 0.1)
        self.assertAlmostEqual(saved["cascade"]["pmax_max"], 5.2)
        self.assertAlmostEqual(saved["cascade"]["psi_min"], 6.3)
        self.assertAlmostEqual(saved["cascade"]["de_min"], 16.4)
        self.assertAlmostEqual(saved["cascade"]["psi_min_surface"], 7.5)
        self.assertAlmostEqual(saved["cascade"]["de_min_surface"], 17.6)
        self.assertEqual(saved["cascade"]["replacement_collisions"], True)

        output = saved["output"]
        self.assertEqual(output["trajectories"]["start"], True)
        self.assertEqual(output["depth_distribution"]["ion_recoils"]["nbins"], 77)
        self.assertEqual(output["depth_distribution"]["ion_recoils"]["limits"], [1.0, 5555.0])
        self.assertEqual(output["lateral_distribution"]["ion_recoils"]["limits"], [-22.0, 33.0])
        self.assertEqual(output["distribution_2d"]["ion_recoils"]["nbins"], [77, 77])
        self.assertEqual(output["distribution_2d"]["ion_recoils"]["limits"], [[1.0, 5555.0], [-22.0, 33.0]])
        self.assertEqual(output["backscattered_atoms"]["energy"]["limits"], [0.0, 12500.0])
        self.assertEqual(output["backscattered_atoms"]["angle"]["limits"], [-80.0, 80.0])


if __name__ == "__main__":
    unittest.main()
