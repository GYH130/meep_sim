"""Geometry and frozen-material checks that never import or run Meep."""

import copy
import json
import math
from pathlib import Path
import unittest

import numpy as np

from ti2d.geometry import periodic_slot_mask, polygon_distance, validate_geometry
from ti2d.materials import audit_model, epsilon, load_model, validate_fit


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "assets/material/route_a_material.json"
REPORT_PATH = MODEL_PATH.with_name("fit_report.json")
RAW_PATH = MODEL_PATH.with_name("ti_ordal_raw.csv")


def geometry(**changes):
    return {"period_um": 50.0, "surface_fill": .6, "axis_length_um": 40.0,
            "tilt_deg": 30.0, "thickness_um": 60.0, "min_wall_um": 1.0,
            "residual_um": 10.0, **changes}


class GeometryTests(unittest.TestCase):
    def test_independent_halfspace_and_derived_quantities(self):
        g = geometry()
        result = validate_geometry(g)
        self.assertTrue(result["valid"], result)
        sine, cosine = .5, math.sqrt(3) / 2
        width = 50 * .6 * cosine
        self.assertAlmostEqual(result["width_normal_um"], width)
        self.assertAlmostEqual(result["surface_opening_um"], 30)
        self.assertAlmostEqual(result["vertical_center_depth_um"], 40 * cosine)
        self.assertAlmostEqual(result["maximum_depth_um"], 40 * cosine + width * sine / 2)
        vertices = np.asarray(result["vertices_x_depth"])
        axis = vertices[:, 0] * sine + vertices[:, 1] * cosine
        normal = vertices[:, 0] * cosine - vertices[:, 1] * sine
        self.assertTrue(np.all(axis <= 40 + 1e-10))
        self.assertTrue(np.all(np.abs(normal) <= width / 2 + 1e-10))
        self.assertTrue(np.all(vertices[:, 1] >= 0))
        top_x = vertices[np.isclose(vertices[:, 1], 0), 0]
        np.testing.assert_allclose(np.sort(top_x), [-15, 15])
        # A fake axis>=0 top cap would incorrectly remove this entrance point.
        self.assertTrue(periodic_slot_mask(-10, 0, g))

    def test_exact_distance_does_not_reject_overlapping_x_bounding_boxes(self):
        g = geometry()
        g.update(period_um=10, surface_fill=.4, tilt_deg=60)
        result = validate_geometry(g)
        self.assertTrue(result["valid"], result)
        self.assertGreater(result["x_span_um"], g["period_um"])
        self.assertGreaterEqual(result["min_wall_um"], 1)
        self.assertGreaterEqual(result["min_wall_um"], result["wall_lower_bound_um"] - 1e-10)
        self.assertLess(min(result["periodic_shifts"]), -1)
        vertices = np.asarray(result["vertices_x_depth"])
        brute = min(polygon_distance(vertices, vertices + [k * 10, 0])
                    for k in range(-12, 13) if k)
        self.assertAlmostEqual(brute, result["min_wall_um"])

    def test_mask_union_periodicity_mirror_and_bottom(self):
        g = geometry()
        g.update(period_um=10, surface_fill=.4, tilt_deg=60)
        x, depth = 20 * math.sin(math.pi / 3), 20 * math.cos(math.pi / 3)
        self.assertTrue(periodic_slot_mask(x, depth, g))
        self.assertTrue(periodic_slot_mask(x - 20, depth, g))
        self.assertFalse(periodic_slot_mask(41 * math.sin(math.pi / 3),
                                           41 * math.cos(math.pi / 3), g))
        rng = np.random.default_rng(5)
        xs, ds = rng.uniform(-50, 50, 1000), rng.uniform(-3, 60, 1000)
        mask = periodic_slot_mask(xs, ds, g)
        np.testing.assert_array_equal(mask, periodic_slot_mask(xs + 70, ds, g))
        mirror = {**g, "tilt_deg": -60}
        np.testing.assert_array_equal(mask, periodic_slot_mask(-xs, ds, mirror))
        self.assertFalse(np.any(mask[ds < 0]))

    def test_strict_residual_surface_bottom_and_input_gates(self):
        g = geometry()
        maximum = validate_geometry(g)["maximum_depth_um"]
        self.assertFalse(validate_geometry({**g, "thickness_um": maximum + 10})["valid"])
        self.assertTrue(validate_geometry({**g, "thickness_um": maximum + 10 + 1e-6})["valid"])
        short = validate_geometry({**g, "axis_length_um": 1})
        self.assertIn("SURFACE_BOTTOM_REGIME", short["status"])
        for change in ({"period_um": math.nan}, {"surface_fill": 1}, {"tilt_deg": 90},
                       {"min_wall_um": .9}, {"residual_um": 9.9}):
            self.assertFalse(validate_geometry({**g, **change})["valid"], change)
        narrow_wall = validate_geometry({**g, "surface_fill": .99})
        self.assertFalse(narrow_wall["valid"])
        self.assertIn("WALL_TOO_THIN", narrow_wall["status"])

    def test_convex_distance_overlap_containment_and_finite_endpoint(self):
        square = np.array([[0, 0], [2, 0], [2, 2], [0, 2]], float)
        self.assertEqual(polygon_distance(square, square + [2, 0]), 0)
        self.assertEqual(polygon_distance(square, square / 4 + [.5, .5]), 0)
        self.assertAlmostEqual(polygon_distance(square, square + [3, 4]), math.sqrt(5))


class MaterialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = load_model(MODEL_PATH)
        cls.report = json.loads(REPORT_PATH.read_text())

    def test_frozen_model_reproduces_source_fit_without_refitting(self):
        before = copy.deepcopy(self.model)
        result = validate_fit(self.model, RAW_PATH, self.report)
        self.assertTrue(result["valid"], result)
        self.assertEqual(result["original_nodes"], 6)
        self.assertEqual(result["interpolated_evaluation_points"], 301)
        self.assertAlmostEqual(result["normalization_scale"], 569.4522969223233)
        self.assertLess(result["normalized_rmse"], .02)
        self.assertLess(result["normalized_max_error"], .05)
        self.assertGreater(result["passivity_min_imag_epsilon"], 0)
        self.assertEqual(self.model, before)

    def test_epsilon_matches_independent_pole_formula(self):
        frequencies = np.array([1 / 6, 1 / 10.5, 1 / 16])
        expected = np.full(3, self.model["epsilon"], dtype=complex)
        for pole in self.model["susceptibilities"]:
            f0, gamma, sigma = (pole[key] for key in ("frequency", "gamma", "sigma"))
            denominator = -frequencies**2 - 1j * gamma * frequencies
            if pole["type"] == "Lorentz":
                denominator += f0**2
            expected += sigma * f0**2 / denominator
        np.testing.assert_allclose(epsilon(self.model, frequencies), expected)

    def test_audit_q_rejection_and_drude_not_resonance(self):
        result = audit_model(self.model, 4, .5, [6, 10.5, 16])
        highest = self.model["susceptibilities"][0]["frequency"]
        self.assertAlmostEqual(result["q"], math.pi * highest * .5 / 4)
        self.assertFalse(result["numerically_qualified"])
        resonance = copy.deepcopy(self.model)
        resonance["susceptibilities"][0]["frequency"] = 4 / math.pi
        with self.assertRaisesRegex(ValueError, "UNSTABLE"):
            audit_model(resonance, 2, .5, 10.5)
        drude = copy.deepcopy(self.model)
        drude["susceptibilities"][1]["frequency"] = 1000
        self.assertAlmostEqual(audit_model(drude, 4, .5, 10.5)["q"], result["q"])

    def test_audit_refuses_ranges_units_nonfinite_and_bad_poles(self):
        for wavelength in (5.99, 16.01, math.nan, -1):
            with self.assertRaises(ValueError):
                audit_model(self.model, 4, .5, wavelength)
        for resolution, courant in ((0, .5), (4, .51), (math.nan, .5), (4, math.inf)):
            with self.assertRaises(ValueError):
                audit_model(self.model, resolution, courant, 10.5)
        for change in ({"meep_length_unit_um": 2}, {"epsilon": math.nan}, {"fit_accepted": False}):
            with self.assertRaises(ValueError):
                audit_model({**self.model, **change}, 4, .5, 10.5)
        negative = copy.deepcopy(self.model)
        negative["susceptibilities"][0]["sigma"] = -1
        with self.assertRaises(ValueError):
            audit_model(negative, 4, .5, 10.5)

    def test_fit_refuses_stale_hash_metadata_and_changed_model(self):
        wrong_hash = {**self.model, "source_sha256": "0" * 64}
        self.assertFalse(validate_fit(wrong_hash, RAW_PATH, self.report)["valid"])
        stale_report = {**self.report, "normalization_scale": 1.0}
        self.assertFalse(validate_fit(self.model, RAW_PATH, stale_report)["valid"])
        changed = copy.deepcopy(self.model)
        changed["susceptibilities"][0]["sigma"] *= 2
        result = validate_fit(changed, RAW_PATH, self.report)
        self.assertFalse(result["valid"])
        self.assertIn("MATERIAL_FIT_RMSE_EXCEEDED", result["errors"])


if __name__ == "__main__":
    unittest.main()
