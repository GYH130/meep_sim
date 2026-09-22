"""Numerical equation checks only; these tests never import or run Meep."""
import cmath
import copy
import json
import math
from pathlib import Path
import unittest

from diagnostic_tools.flat_discrete import discrete_epsilon, planar_tm_prediction

MODEL = json.loads((Path(__file__).resolve().parents[1] /
                    "assets/material/route_a_material.json").read_text())


class FlatDiscreteTests(unittest.TestCase):
    def predict(self, resolution=16, **kwargs):
        return planar_tm_prediction(MODEL, 10.5, 30, resolution, **kwargs)

    def test_frozen_reference_predictions_and_passive_depth(self):
        for res, expected in ((16, .9611893951182707), (20, .9599548292785641),
                              (24, .9591838911132063)):
            result = self.predict(res)
            self.assertAlmostEqual(result["R"], expected, places=12)
            self.assertLess(abs(result["depth_factor"]), 1)
            self.assertTrue(result["prediction_only"])
            self.assertFalse(result["fdtd_run"])
            self.assertFalse(result["convergence_certified"])

    def test_continuum_limit_and_no_contrast(self):
        eps = discrete_epsilon(MODEL, 1 / 10.5, 0)
        impedance = cmath.sqrt(eps - .25) / eps
        fresnel = abs((impedance - math.sqrt(.75)) / (impedance + math.sqrt(.75)))**2
        self.assertLess(abs(self.predict(100000)["R"] - fresnel), 1e-6)
        air = {"epsilon": 1, "susceptibilities": []}
        self.assertLess(planar_tm_prediction(air, 10.5, 30, 16)["R"], 1e-20)

    def test_timestep_limit_and_interface_override(self):
        continuum = discrete_epsilon(MODEL, 1 / 10.5, 0)
        coarse = abs(discrete_epsilon(MODEL, 1 / 10.5, .015625) - continuum)
        fine = abs(discrete_epsilon(MODEL, 1 / 10.5, .0078125) - continuum)
        self.assertAlmostEqual(coarse / fine, 4, places=4)
        bulk = self.predict(interface_mode="bulk")
        override = self.predict(interface_epsilon=bulk["epsilon_discrete"])
        self.assertEqual(bulk["R"], override["R"])
        self.assertLess(abs(self.predict(courant=.125)["R"] - self.predict()["R"]), 2e-7)

    def test_invalid_inputs(self):
        for kwargs in ({"courant": 0}, {"courant": .6}, {"interface_mode": "unknown"},
                       {"interface_epsilon": 1-1j}, {"interface_epsilon": complex("nan")}):
            with self.assertRaises(ValueError):
                self.predict(**kwargs)
        for wave, theta, res in ((0, 30, 16), (10.5, 90, 16), (10.5, 30, -1)):
            with self.assertRaises(ValueError):
                planar_tm_prediction(MODEL, wave, theta, res)
        active = copy.deepcopy(MODEL)
        active["susceptibilities"][0]["gamma"] = -1
        with self.assertRaises(ValueError):
            discrete_epsilon(active, 1 / 10.5, .015625)


if __name__ == "__main__":
    unittest.main()
