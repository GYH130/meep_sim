"""Fake Meep arrays exercise strict coordinate/shape and periodicity contracts."""

from types import SimpleNamespace
import unittest

import numpy as np

from ti2d.meep_arrays import native_volume_pairs, plane_fields


PERIOD, RESOLUTION, KX = 4.0, 4.0, .13


def corners(x, y):
    return (np.array([len(x), len(y), 1]),
            SimpleNamespace(x=x[0], y=y[0], z=0),
            SimpleNamespace(x=x[-1], y=y[-1], z=0))


class FakeSimulation:
    def __init__(self):
        self.monitor = SimpleNamespace(where=object())
        self.arrays, self.dimensions = {}, {}
        self.metadata_calls = 0
        self.dimension_calls = []
        self.x = np.arange(18) / RESOLUTION - 2.125
        self.weights = np.full(18, .25)
        self.weights[[0, -1]] = 0

    def add_pair(self, name, x, y):
        field = np.exp(2j * np.pi * KX * x[:, None]) * (1 + y[None, :] ** 2)
        self.arrays[name] = field
        self.arrays["D" + name[1:]] = (2 + 3j) * field
        self.dimensions[name] = corners(x, y)
        self.dimensions["D" + name[1:]] = corners(x, y)

    def get_dft_array(self, monitor, component, frequency):
        assert monitor is self.monitor
        assert frequency == 0
        return self.arrays[component]

    def get_array_slice_dimensions(self, component, vol):
        assert vol is self.monitor.where
        self.dimension_calls.append(component)
        return self.dimensions[component]

    def get_array_metadata(self, dft_cell):
        assert dft_cell is self.monitor
        self.metadata_calls += 1
        return self.x, np.array([2.0]), np.array([0.0]), self.weights


class NativeArrayTests(unittest.TestCase):
    def setUp(self):
        self.sim = FakeSimulation()
        self.sim.add_pair("Ex", np.arange(18) / 4 - 2.125, np.arange(17) / 4 - 3.5)
        self.sim.add_pair("Ey", np.arange(17) / 4 - 2, np.arange(18) / 4 - 3.625)
        self.pairs = {"Ex": ("Ex", "Dx"), "Ey": ("Ey", "Dy")}
        self.geometry = {"case": "flat", "thickness_um": 3}

    def read(self):
        return native_volume_pairs(self.sim, self.sim.monitor, self.pairs,
                                   RESOLUTION, PERIOD, KX, self.geometry, 0)

    def test_native_grids_remain_component_specific_and_complete_in_y(self):
        pairs, raw, metadata = self.read()
        self.assertEqual(self.sim.metadata_calls, 0)
        self.assertEqual(self.sim.dimension_calls, ["Ex", "Dx", "Ey", "Dy"])
        self.assertEqual(pairs[0]["E"].shape, (16, 17))
        self.assertEqual(pairs[1]["E"].shape, (16, 18))
        self.assertEqual(raw["Ex_E_full"].shape, (18, 17))
        self.assertEqual(raw["Ey_E_full"].shape, (17, 18))
        self.assertEqual(metadata["y_samples_removed"], 0)
        for pair in pairs:
            np.testing.assert_array_equal(pair["weights"], np.full(pair["E"].shape, 1 / 16))
            y = pair["coordinates"]["y"]
            expected = np.broadcast_to((y >= -3) & (y <= 0), pair["E"].shape)
            np.testing.assert_array_equal(pair["ti_mask"], expected)
        self.assertFalse(np.array_equal(pairs[0]["coordinates"]["x"], pairs[1]["coordinates"]["x"]))

    def test_slot_mask_on_native_coordinates(self):
        self.geometry = {"period_um": 4, "surface_fill": .5, "axis_length_um": 2,
                         "tilt_deg": 0, "thickness_um": 15}
        pairs, _, _ = self.read()
        for pair in pairs:
            x, y = np.meshgrid(pair["coordinates"]["x"], pair["coordinates"]["y"], indexing="ij")
            expected = (y >= -15) & (y <= 0) & ~((abs(x) <= 1) & (y >= -2) & (y <= 0))
            np.testing.assert_array_equal(pair["ti_mask"], expected)

    def test_no_shape_cropping_or_coordinate_repair(self):
        self.sim.arrays["Dx"] = self.sim.arrays["Dx"][:-1]
        with self.assertRaisesRegex(ValueError, "SHAPE_MISMATCH"):
            self.read()

    def test_reject_pair_on_distinct_grids_despite_equal_shapes(self):
        dims, lo, hi = self.sim.dimensions["Dx"]
        self.sim.dimensions["Dx"] = (dims, SimpleNamespace(x=lo.x + .125, y=lo.y, z=0),
                                      SimpleNamespace(x=hi.x + .125, y=hi.y, z=0))
        with self.assertRaisesRegex(ValueError, "NOT_COLLOCATED"):
            self.read()

    def test_reject_inconsistent_corner_and_bad_bloch_guard(self):
        dims, lo, hi = self.sim.dimensions["Ex"]
        original_hi = hi.x
        hi.x += .125
        with self.assertRaisesRegex(ValueError, "CORNER_SPACING"):
            self.read()
        hi.x = original_hi
        self.sim.arrays["Ex"][0, 0] += 1
        with self.assertRaisesRegex(ValueError, "BLOCH_GUARD_MISMATCH"):
            self.read()


class PlaneArrayTests(unittest.TestCase):
    def setUp(self):
        self.sim = FakeSimulation()
        field = np.exp(2j * np.pi * KX * self.sim.x)
        self.sim.arrays = {"Ez": field, "Hx": -2 * field}
        self.components = {"Ez": "Ez", "Hx": "Hx"}

    def read(self):
        return plane_fields(self.sim, self.sim.monitor, self.components, PERIOD, KX, RESOLUTION)

    def test_centered_plane_full_weights_and_canonical_selection_retained(self):
        x, arrays, raw = self.read()
        self.assertEqual(len(x), 16)
        self.assertEqual(arrays["Ez"].shape, (16,))
        self.assertEqual(raw["Ez_full"].shape, (18,))
        self.assertAlmostEqual(raw["weights_full"].sum(), 4)
        self.assertAlmostEqual(raw["weights_canonical_uniform"].sum(), 4)
        np.testing.assert_array_equal(raw["x_selected_indices"], np.arange(1, 17))
        self.assertEqual(self.sim.dimension_calls, [])

    def test_reject_nonvector_plane_and_weight_mismatch(self):
        self.sim.arrays["Ez"] = self.sim.arrays["Ez"][:, None]
        with self.assertRaisesRegex(ValueError, "SHAPE_MISMATCH"):
            self.read()
        self.sim.arrays["Ez"] = self.sim.arrays["Ez"][:, 0]
        self.sim.weights = self.sim.weights[:-1]
        with self.assertRaisesRegex(ValueError, "METADATA_SHAPE_MISMATCH"):
            self.read()

    def test_reject_plane_guard_instead_of_ignoring_it(self):
        self.sim.arrays["Hx"][-1] += .1
        with self.assertRaisesRegex(ValueError, "BLOCH_GUARD_MISMATCH"):
            self.read()


if __name__ == "__main__":
    unittest.main()
