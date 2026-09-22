"""Offline diagnostic contracts and overflow regressions; no Meep execution."""
import ast
import copy
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from diagnostic_tools import slot_stability_probe as probe
from ti2d.common import digest_file, digest_object


GEOMETRY = {"period_um": 20.0, "surface_fill": .45, "axis_length_um": 35.0,
            "tilt_deg": 30.0, "thickness_um": 60.0, "min_wall_um": 1.0, "residual_um": 10.0}


class SafeFieldArithmeticTests(unittest.TestCase):
    def fields(self, value=0j):
        return np.full((2, 3), value, dtype=complex), np.zeros((2, 3), dtype=complex)

    def summarize(self, ex, ey):
        return probe.field_summary(ex, ey, np.array([0.0, .03125]),
                                   np.array([-.03125, 0.0, .03125]), GEOMETRY, 0.0)[0]

    def test_original_1e302_field_preserved_without_power_overflow(self):
        ex, ey = self.fields(1e302 + 1e302j)
        ey[0, 1] = 1e302 - 1e302j
        with np.errstate(all="raise"):
            result = self.summarize(ex, ey)
        self.assertTrue(result["all_fields_finite"])
        self.assertIsNone(result["envelope_intensity"])
        self.assertAlmostEqual(result["log10_envelope_intensity"], 604 + math.log10(4))
        self.assertAlmostEqual(result["envelope_amplitude"] / 1e302, 2)
        self.assertEqual(result["raw_Ex_at_location"], {"real": 1e302, "imag": 1e302})
        self.assertEqual(result["index_xy"], [0, 1])
        json.dumps(result, allow_nan=False)

    def test_total_amplitude_can_exceed_binary64_without_nonfinite_input(self):
        ex, ey = self.fields(1.7e308 + 1.7e308j)
        ey[:] = ex
        with np.errstate(all="raise"):
            result = self.summarize(ex, ey)
        self.assertTrue(result["all_fields_finite"])
        self.assertIsNone(result["envelope_amplitude"])
        self.assertIsNone(result["envelope_intensity"])
        self.assertTrue(math.isfinite(result["log10_envelope_intensity"]))
        self.assertEqual(result["raw_Ex_at_location"]["real"], 1.7e308)
        json.dumps(result, allow_nan=False)

    def test_normal_values_match_direct_power_and_peak(self):
        rng = np.random.default_rng(4)
        ex, ey = [rng.normal(size=(2, 3)) + 1j * rng.normal(size=(2, 3)) for _ in range(2)]
        expected = np.abs(ex) ** 2 + np.abs(ey) ** 2
        result = self.summarize(ex, ey)
        self.assertAlmostEqual(result["envelope_intensity"], expected.max())
        self.assertEqual(result["index_xy"], list(np.unravel_index(expected.argmax(), expected.shape)))

    def test_zero_has_explicit_finite_encoding(self):
        result = self.summarize(*self.fields())
        self.assertEqual(result["envelope_intensity"], 0)
        self.assertEqual(result["envelope_amplitude"], 0)
        self.assertTrue(result["zero_field"])
        self.assertIsNone(result["log10_envelope_intensity"])
        self.assertIsNone(result["log10_intensity_over_reference_peak"])
        json.dumps(result, allow_nan=False)

    def test_tiny_nonzero_field_preserves_logarithm(self):
        result = self.summarize(*self.fields(1e-300))
        self.assertFalse(result["zero_field"])
        self.assertEqual(result["envelope_intensity"], 0)
        self.assertAlmostEqual(result["log10_envelope_intensity"], -600)
        self.assertIn("underflows", result["intensity_encoding_note"])

    def test_nonfinite_immediately_detected_and_json_safe(self):
        for value in (complex(math.inf, 0), complex(0, -math.inf), complex(math.nan, 0)):
            ex, ey = self.fields(2)
            ey[1, 2] = value
            result = self.summarize(ex, ey)
            self.assertFalse(result["all_fields_finite"])
            self.assertEqual(result["index_xy"], [1, 2])
            self.assertEqual(result["nonfinite_grid_points"], 1)
            self.assertEqual(probe.GrowthGuard().sample(result), "NONFINITE_FIELD")
            json.dumps(result, allow_nan=False)

    def test_point_fields_safe_at_extremes(self):
        finite = probe.point_fields_summary(1e302 + 1e302j, 1e302)
        self.assertTrue(finite["finite"])
        self.assertIsNone(finite["intensity"])
        self.assertAlmostEqual(finite["log10_intensity"], 604 + math.log10(3))
        nonfinite = probe.point_fields_summary(math.nan, 0)
        self.assertFalse(nonfinite["finite"])
        json.dumps([finite, nonfinite], allow_nan=False)

    def test_shape_mismatch_never_crops(self):
        for ex, ey in ((np.zeros((3, 4)), np.zeros((3, 5))),
                       (np.zeros(4), np.zeros(4)), (np.zeros((0, 4)), np.zeros((0, 4)))):
            with self.assertRaisesRegex(ValueError, "FIELD_SHAPE_MISMATCH"):
                probe.log_power_arrays(ex, ey)


class GrowthGuardTests(unittest.TestCase):
    @staticmethod
    def sample(ratio):
        return {"all_fields_finite": True, "log10_intensity_over_reference_peak": ratio}

    def test_requires_all_three_large_samples_and_both_strict_growth_jumps(self):
        guard = probe.GrowthGuard()
        observed = [guard.sample(self.sample(ratio)) for ratio in (10, 12.2, 13.3, 14.4)]
        self.assertEqual(observed, [None, None, None, "EARLY_GROWTH_DETECTED"])

    def test_large_plateau_or_single_spike_does_not_pass_growth_condition(self):
        for sequence in ((13, 13, 13, 13), (13, 15, 14), (13, 14, 15), (12, 14, 16), (13, None, 16)):
            guard = probe.GrowthGuard()
            self.assertTrue(all(guard.sample(self.sample(item)) is None for item in sequence))

    def test_resets_effectively_after_small_or_zero_sample(self):
        guard = probe.GrowthGuard()
        values = [guard.sample(self.sample(item)) for item in (13, 15, None, 17, 19, 21)]
        self.assertEqual(values, [None, None, None, None, None, "EARLY_GROWTH_DETECTED"])


class GridAndPatchTests(unittest.TestCase):
    def setUp(self):
        self.x = np.arange(25, dtype=float) / 32
        self.y = np.arange(30, dtype=float) / 32 - 1
        self.ex = np.zeros((25, 30), complex)
        self.ey = self.ex.copy()
        self.metadata = [self.x, self.y, np.array([0.0]), np.ones((25, 30)) / 32 ** 2]

    def test_exact_shape_and_spacing_required(self):
        probe.checked_grid(self.ex, self.ey, self.metadata, 32)
        variants = [(self.ex[:, :-1], self.ey, self.metadata),
                    (self.ex, self.ey[:-1], self.metadata),
                    (self.ex, self.ey, [self.x, self.y, np.array([0.0]), np.ones((24, 30))]),
                    (self.ex, self.ey, [self.x + np.arange(25) * .001, *self.metadata[1:]])]
        for ex, ey, metadata in variants:
            with self.assertRaises(ValueError):
                probe.checked_grid(ex, ey, metadata, 32)

    def test_nonfinite_metadata_refused_before_field_interpretation(self):
        self.metadata[0][4] = math.nan
        with self.assertRaisesRegex(ValueError, "GRID_NONFINITE"):
            probe.checked_grid(self.ex, self.ey, self.metadata, 32)

    def test_small_patch_keeps_original_values_and_coordinate_masks(self):
        self.ex[24, 29] = 1e302j
        saved = probe.local_patch(self.ex, self.ey, self.x, self.y, (24, 29), GEOMETRY, 0)
        self.assertEqual(saved["Ex_raw"].shape, (17, 17))
        self.assertEqual(saved["Ex_raw"][-1, -1], 1e302j)
        np.testing.assert_array_equal(saved["x_um"], self.x[-17:])
        np.testing.assert_array_equal(saved["y_um"], self.y[-17:])
        np.testing.assert_array_equal(saved["ideal_ti_mask"], saved["ideal_slab_mask"] & ~saved["ideal_slot_mask"])
        np.testing.assert_array_equal(saved["global_peak_index_xy"], [24, 29])


class CollectiveFailureTests(unittest.TestCase):
    def test_local_error_goes_through_allreduce_then_allgather(self):
        calls = []
        comm = SimpleNamespace(rank=2, allreduce=lambda value, op: calls.append(("reduce", value)) or value,
                               allgather=lambda value: calls.append(("gather", value)) or [None, None, value, None])
        def fail():
            raise ValueError("bad shape")
        with self.assertRaisesRegex(probe.CollectiveDiagnosticError, "bad shape"):
            probe.collective_call(comm, SimpleNamespace(MAX=1), "shape", fail)
        self.assertEqual([call[0] for call in calls], ["reduce", "gather"])

    def test_remote_error_is_propagated_even_after_local_success(self):
        comm = SimpleNamespace(rank=0, allreduce=lambda value, op: 1,
                               allgather=lambda value: [value, {"rank": 1, "error": "peer grid mismatch"}])
        with self.assertRaisesRegex(probe.CollectiveDiagnosticError, "peer grid mismatch"):
            probe.collective_call(comm, SimpleNamespace(MAX=1), "shape", lambda: 42)


class IsolationTests(unittest.TestCase):
    def test_import_does_not_import_meep_mpi_or_scipy(self):
        command = [sys.executable, "-B", "-c",
                   "import sys; import diagnostic_tools.slot_stability_probe; "
                   "assert not any(n.split('.')[0] in ('meep','mpi4py','scipy') for n in sys.modules)"]
        subprocess.run(command, cwd=Path(probe.__file__).resolve().parents[1], check=True, capture_output=True)

    def test_local_runtime_refused_before_meep_import(self):
        with patch.object(sys, "platform", "darwin"):
            with self.assertRaisesRegex(ValueError, "MEEP_SERVER_ONLY"):
                probe.main("/nonexistent/plan.json", "averaging_on", "/nonexistent")

    def test_simulation_arguments_match_production_except_explicit_averaging(self):
        root = Path(probe.__file__).resolve().parents[1]
        def simulation_keywords(path):
            tree = ast.parse(path.read_text())
            return [{item.arg: ast.unparse(item.value) for item in node.keywords}
                    for node in ast.walk(tree) if isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "mp" and node.func.attr == "Simulation"][0]
        production = simulation_keywords(root / "ti2d/solver.py")
        diagnostic = simulation_keywords(Path(probe.__file__))
        self.assertEqual(diagnostic.pop("eps_averaging"), "ARMS[case_id]")
        diagnostic = {key: value.replace("case[", "c[").replace("comm.size", "ranks") for key, value in diagnostic.items()}
        self.assertEqual(diagnostic, production)


class PlanValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.stage = Path(self.temp.name)
        module = self.stage / "diagnostic_tools/slot_stability_probe.py"
        module.parent.mkdir()
        module.write_text("offline probe fixture")
        material = self.stage / "model.json"
        material.write_text('{"fixture":true}')
        case = {"implementation_sha256": probe.FROZEN_CORE_SHA256, "resolution": 32,
                "courant": .25, "polarization": "p", "emission_theta_deg": 30,
                "case": "structure", "mpi_ranks": 4, "oxide_layer_in_model": False,
                "material_path": str(material), "material_sha256": digest_file(material)}
        source = self.stage / "source_case.json"
        source.write_text(json.dumps(case))
        self.run = self.stage / "repair_diagnostics/slot_stability_fixture"
        self.run.mkdir(parents=True)
        self.plan_path = self.run / "plan.json"
        self.plan = {"schema_version": 1, "run_id": self.run.name, "stage_root": str(self.stage),
                     "run_dir": str(self.run), "source_case_path": str(source),
                     "source_case_sha256": digest_file(source), "implementation_sha256": probe.FROZEN_CORE_SHA256,
                     "material_sha256": digest_file(material), "probe_sha256": digest_file(module),
                     "cases": [{"id": key, "eps_averaging": value} for key, value in probe.ARMS.items()],
                     "diagnostic": {"stop_time": 100.0, "sample_interval": 2.0,
                                    "reference_intensity_peak": probe.REFERENCE_INTENSITY_PEAK,
                                    "growth_ratio_limit": 1e12, "growth_step_factor": 10.0,
                                    "growth_consecutive_samples": 3, "patch_radius_cells": 8},
                     "budgets": {"mpi_ranks": 4}}
        for target, replacement in (("__file__", str(module)), ("FROZEN_CASE_SHA256", digest_file(source))):
            patcher = patch.object(probe, target, replacement)
            patcher.start()
            self.addCleanup(patcher.stop)
        patcher = patch.object(probe, "implementation_digest", return_value=probe.FROZEN_CORE_SHA256)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.write_plan()

    def write_plan(self):
        self.plan["plan_sha256"] = digest_object({key: value for key, value in self.plan.items() if key != "plan_sha256"})
        self.plan_path.write_text(json.dumps(self.plan))

    def inspect(self, arm="averaging_on", output=None):
        return probe.load_inspection(self.plan_path, arm, output or self.run / arm)

    def test_two_arms_same_immutable_physical_config_and_never_qualify(self):
        on, off = self.inspect()[2], self.inspect("averaging_off")[2]
        self.assertEqual(on["case_config"], off["case_config"])
        self.assertEqual(on["physical_case_sha256"], off["physical_case_sha256"])
        self.assertTrue(on["eps_averaging"])
        self.assertFalse(off["eps_averaging"])
        self.assertFalse(on["fdtd_started"])
        self.assertFalse(on["qualification"])
        self.assertTrue(on["causal_interpretation_requires_averaging_on_growth_reproduction"])

    def test_frozen_inputs_and_stop_limits_fail_closed(self):
        initial = copy.deepcopy(self.plan)
        for field, value in (("implementation_sha256", "0" * 64), ("probe_sha256", "0" * 64),
                             ("source_case_sha256", "0" * 64), ("material_sha256", "0" * 64),
                             ("cases", [{"id": "averaging_on", "eps_averaging": False}])):
            self.plan = copy.deepcopy(initial)
            self.plan[field] = value
            self.write_plan()
            with self.assertRaises(ValueError):
                self.inspect()
        self.plan = copy.deepcopy(initial)
        self.plan["diagnostic"]["stop_time"] = 526
        self.write_plan()
        with self.assertRaisesRegex(ValueError, "DIAGNOSTIC_LIMITS_CHANGED"):
            self.inspect()

    def test_plan_signature_and_output_scope_enforced(self):
        self.plan["plan_sha256"] = "0" * 64
        self.plan_path.write_text(json.dumps(self.plan))
        with self.assertRaisesRegex(ValueError, "PLAN_HASH_MISMATCH"):
            self.inspect()
        self.write_plan()
        with self.assertRaisesRegex(ValueError, "OUTPUT_OR_PLAN_PATH"):
            self.inspect(output=self.stage / "runs/production")

    def test_existing_evidence_blocks_retry(self):
        output = self.run / "averaging_on"
        output.mkdir()
        (output / "trace.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "RETRY_OR_OVERWRITE"):
            self.inspect()


if __name__ == "__main__":
    unittest.main()
