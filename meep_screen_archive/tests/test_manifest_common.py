"""Offline campaign, periodic-grid, and exact-cache contract checks.

All generated manifests and synthetic result evidence live in temporary
directories. The test suite never imports Meep or launches a solver.
"""

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from ti2d import manifest
from ti2d.common import atomic_json, atomic_npz, canonical_period, digest_file, digest_object, jsonable
from ti2d.queue import aggregate_gate, case_fingerprint, reusable_result, seal_result, verify_manifest
from ti2d.solver import layout, result_checks


class CanonicalPeriodTests(unittest.TestCase):
    def setUp(self):
        self.period = 4.0
        self.resolution = 4
        self.kx = -0.137
        self.x = np.arange(-8, 8) / self.resolution

    def fields(self, x):
        envelope = 1 + 0.4 * np.cos(2 * np.pi * x / self.period)
        carrier = np.exp(2j * np.pi * self.kx * x)
        return {"E": envelope * carrier, "H": (-0.3 + 0.2j) * envelope * carrier}

    def call(self, x=None, arrays=None, **kwargs):
        x = self.x if x is None else x
        arguments = dict(period=self.period, kx=self.kx, resolution=self.resolution)
        arguments.update(kwargs)
        return canonical_period(x, self.fields(x) if arrays is None else arrays, **arguments)

    def test_retains_full_bloch_phase_without_unrequested_shift(self):
        original = self.fields(self.x)
        x, arrays, keep = self.call(arrays=original)
        np.testing.assert_array_equal(x, self.x)
        self.assertTrue(np.all(keep))
        for name in original:
            np.testing.assert_array_equal(arrays[name], original[name])
        self.assertGreater(abs(arrays["E"][0].imag), 0.01)

    def test_midpoints_are_retained_with_exact_period_measure(self):
        midpoint = self.x + 0.5 / self.resolution
        x, arrays, keep = self.call(x=midpoint)
        np.testing.assert_array_equal(x, midpoint)
        self.assertEqual(len(x) / self.resolution, self.period)
        self.assertTrue(np.all(keep))

    def test_endpoint_duplicate_removed_once_with_bloch_audit(self):
        endpoints = np.arange(-8, 9) / self.resolution
        x, arrays, keep = self.call(x=endpoints)
        self.assertEqual(len(x), 16)
        self.assertEqual(int(keep.sum()), 16)
        self.assertFalse(keep[-1])
        np.testing.assert_allclose(arrays["E"], self.fields(self.x)["E"])
        self.assertAlmostEqual(np.sum(np.ones(len(x)) / self.resolution), self.period)

    def test_two_interpolation_guards_removed_without_double_weight(self):
        guarded = np.arange(-9, 9) / self.resolution
        x, arrays, keep = self.call(x=guarded)
        self.assertFalse(keep[0])
        self.assertFalse(keep[-1])
        self.assertEqual(int(keep.sum()), 16)
        np.testing.assert_array_equal(x, self.x)
        np.testing.assert_allclose(arrays["E"], self.fields(self.x)["E"])

    def test_multidimensional_fields_preserve_other_axes(self):
        guarded = np.arange(-9, 9) / self.resolution
        scalar = self.fields(guarded)["E"]
        values = scalar[:, None, None] * np.array([1, 2, 3])[None, :, None]
        x, arrays, keep = self.call(x=guarded, arrays={"E": values})
        self.assertEqual(arrays["E"].shape, (16, 3, 1))
        np.testing.assert_array_equal(arrays["E"], values[keep])

    def test_corrupt_guard_and_wrong_bloch_kx_fail(self):
        guarded = np.arange(-9, 9) / self.resolution
        arrays = self.fields(guarded)
        arrays["E"][0] += 0.001
        with self.assertRaisesRegex(ValueError, "BLOCH_GUARD_MISMATCH"):
            self.call(x=guarded, arrays=arrays)
        with self.assertRaisesRegex(ValueError, "BLOCH_GUARD_MISMATCH"):
            self.call(x=guarded, kx=-self.kx)

    def test_shape_mismatch_fails_without_slicing_or_broadcasting(self):
        with self.assertRaisesRegex(ValueError, "SHAPE_MISMATCH"):
            self.call(arrays={"E": np.ones(15)})
        with self.assertRaises(ValueError):
            self.call(arrays={"E": np.array(1.0)})

    def test_nonuniform_missing_and_noninteger_period_grids_fail(self):
        moved = self.x.copy()
        moved[3] += 0.01
        for x, kwargs in ((moved, {}), (self.x[:-1], {}),
                          (self.x, {"period": 4.1}), (self.x[::-1], {})):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.call(x=x, **kwargs)

    def test_nonfinite_and_complex_coordinates_rejected(self):
        for value in (np.nan, np.inf):
            x = self.x.copy()
            x[3] = value
            with self.assertRaises(ValueError):
                self.call(x=x, arrays={"E": np.ones(16)})
        with self.assertRaises(ValueError):
            self.call(x=self.x.astype(complex) + 1j, arrays={"E": np.ones(16)})

    def test_nonfinite_fields_fail_with_and_without_guards(self):
        for x in (self.x, np.arange(-9, 9) / self.resolution):
            for value in (np.nan, np.inf):
                arrays = self.fields(x)
                arrays["E"][0] = value
                with self.subTest(guards=len(x) > 16, value=value), self.assertRaises(ValueError):
                    self.call(x=x, arrays=arrays)

    def test_finite_positive_period_resolution_and_real_kx_required(self):
        for kwargs in ({"period": 0}, {"period": -1}, {"period": np.nan},
                       {"resolution": 0}, {"resolution": -1}, {"resolution": np.inf},
                       {"kx": np.nan}, {"kx": np.inf}, {"kx": 0.1 + 0.3j},
                       {"kx": np.complex128(0.1 + 0.3j)},
                       {"resolution": np.complex128(4 + 0.3j)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.call(**kwargs)


class CommonArtifactTests(unittest.TestCase):
    def test_canonical_digest_is_order_independent_and_numpy_aware(self):
        first = {"b": np.array([1, 2]), "a": np.float64(0.25), "z": 1 + 2j}
        second = {"z": {"real": 1.0, "imag": 2.0}, "a": 0.25, "b": [1, 2]}
        self.assertEqual(digest_object(first), digest_object(second))
        self.assertNotEqual(digest_object(first), digest_object({**first, "a": 0.26}))
        json.dumps(jsonable(first), allow_nan=False)

    def test_nonfinite_serialization_does_not_overwrite_existing_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.json"
            atomic_json(path, {"valid": True, "power": 0.2})
            original = path.read_bytes()
            for bad in (np.nan, np.inf, complex(1, np.nan)):
                with self.assertRaises(ValueError):
                    atomic_json(path, {"power": bad})
                self.assertEqual(path.read_bytes(), original)
            self.assertEqual([p.name for p in Path(directory).iterdir()], ["evidence.json"])

    def test_npz_preserves_complex_arrays_and_rank(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fields.npz"
            values = (np.arange(12).reshape(4, 3, 1) + 1j) * 0.25
            atomic_npz(path, E=values)
            with np.load(path, allow_pickle=False) as stored:
                np.testing.assert_array_equal(stored["E"], values)
            self.assertEqual([p.name for p in Path(directory).iterdir()], ["fields.npz"])


class ManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="ti2d-manifest-test-")
        cls.run_dir = (Path(cls.temporary.name) / "locked_run").resolve()
        root = Path(manifest.__file__).resolve().parents[1]
        cls.source_assets = {path: digest_file(path) for path in (root / "assets").rglob("*") if path.is_file()}
        cls.campaign = manifest.create(cls.run_dir)
        cls.tasks = {task["id"]: task for task in cls.campaign["tasks"]}
        cls.configs = {name: json.loads(Path(task["case_path"]).read_text()) for name, task in cls.tasks.items()}

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_exactly_twenty_design_conditions_and_thirty_three_validation_cases(self):
        validation = [task for task in self.tasks.values() if task["phase"] == "validation"]
        pilot = [task for task in self.tasks.values() if task["phase"] == "pilot"]
        self.assertEqual(len(validation), 33)
        self.assertEqual(len(pilot), 20)
        self.assertEqual(len(self.tasks), 53)
        self.assertEqual(set(self.campaign["production_gate_task_ids"]), {task["id"] for task in validation})
        expected = {(geometry, pol, angle) for geometry in
                    ("baseline", "period_high", "duty_high", "length_high", "tilt_high")
                    for pol in ("s", "p") for angle in (-30.0, 30.0)}
        actual = {(task["geometry_id"], task["polarization"], task["emission_theta_deg"]) for task in pilot}
        self.assertEqual(actual, expected)

    def test_every_case_and_source_asset_exists_inside_expected_paths(self):
        self.assertTrue(self.source_assets)
        for path, original_hash in self.source_assets.items():
            self.assertTrue(path.is_file())
            self.assertEqual(digest_file(path), original_hash)
        for name, task in self.tasks.items():
            case_path = Path(task["case_path"])
            self.assertTrue(case_path.is_file())
            self.assertTrue(case_path.is_relative_to(self.run_dir))
            material = Path(self.configs[name]["material_path"])
            self.assertTrue(material.is_file())
            self.assertEqual(digest_file(material), self.configs[name]["material_sha256"])
        self.assertTrue((self.run_dir / "material_recheck.json").is_file())
        self.assertTrue((self.run_dir / "geometries.json").is_file())

    def test_manifest_digest_and_existing_path_are_immutable(self):
        payload = deepcopy(self.campaign)
        saved_digest = payload.pop("manifest_sha256")
        self.assertEqual(saved_digest, digest_object(payload))
        path = self.run_dir / "manifest.json"
        originals = {file: digest_file(file) for file in self.run_dir.rglob("*") if file.is_file()}
        with self.assertRaises(FileExistsError):
            manifest.create(self.run_dir)
        self.assertEqual({file: digest_file(file) for file in self.run_dir.rglob("*") if file.is_file()}, originals)
        self.assertEqual(json.loads(path.read_text()), self.campaign)

    def test_runtime_admission_checks_manifest_and_implementation_before_first_lock(self):
        implementation = self.campaign["implementation_sha256"]
        verify_manifest(self.campaign, self.run_dir, implementation)
        changed = deepcopy(self.campaign)
        changed["budgets"]["single_solver_seconds"] += 1
        with self.assertRaisesRegex(ValueError, "MANIFEST_PAYLOAD_FINGERPRINT"):
            verify_manifest(changed, self.run_dir, implementation)
        with self.assertRaisesRegex(ValueError, "IMPLEMENTATION_CHANGED"):
            verify_manifest(self.campaign, self.run_dir, "different-implementation")

    def test_case_edit_after_manifest_creation_cannot_be_admitted_as_a_new_first_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            campaign = manifest.create(Path(directory) / "new_run")
            case_path = Path(campaign["tasks"][0]["case_path"])
            changed = json.loads(case_path.read_text())
            changed["wavelength_um"] = 9.0
            atomic_json(case_path, changed)
            with self.assertRaisesRegex(ValueError, "CASE_PAYLOAD_FINGERPRINT"):
                verify_manifest(campaign, case_path.parent.parent, campaign["implementation_sha256"])
            self.assertFalse((case_path.parent.parent / "lock_manifest.json").exists())

    def test_geometry_requirements_and_one_parameter_at_a_time(self):
        evidence = json.loads((self.run_dir / "geometries.json").read_text())
        definitions, audits = evidence["definitions"], evidence["audit"]
        baseline = definitions["baseline"]
        self.assertEqual(baseline["period_um"], 50)
        self.assertEqual(baseline["axis_length_um"], 40)
        self.assertEqual(baseline["tilt_deg"], 30)
        self.assertAlmostEqual(audits["baseline"]["width_normal_um"], 30, places=7)
        changed = {"period_high": ("period_um", 58), "duty_high": ("surface_fill", 0.8),
                   "length_high": ("axis_length_um", 48), "tilt_high": ("tilt_deg", 45),
                   "straight": ("tilt_deg", 0), "mirror": ("tilt_deg", -30)}
        for geometry, (key, value) in changed.items():
            self.assertEqual(definitions[geometry], {**baseline, key: value})
        for name, geometry in definitions.items():
            audit = audits[name]
            self.assertTrue(audit["valid"])
            self.assertEqual(geometry["thickness_um"], 60)
            self.assertEqual(geometry["min_wall_um"], 1)
            self.assertEqual(geometry["residual_um"], 10)
            self.assertGreaterEqual(audit["min_wall_um"], 1)
            self.assertGreater(audit["residual_um"], 10)
            self.assertTrue(audit["vertices_x_depth"])

    def test_locked_physics_and_stop_settings_are_consistent(self):
        for name, config in self.configs.items():
            self.assertEqual(config["id"], name)
            self.assertEqual(config["wavelength_um"], 10.5)
            self.assertIn(config["emission_theta_deg"], (-30.0, 30.0))
            self.assertIn(config["polarization"], ("s", "p"))
            self.assertEqual(config["courant"], 0.25)
            self.assertEqual(config["stop_window_um"], 20)
            self.assertEqual(config["stop_consecutive_windows"], 3)
            self.assertEqual(config["flux_drift_tolerance"], 0.001)
            self.assertEqual(config["intensity_decay"], 1e-6 if config["strict_stop"] else 1e-4)
            self.assertFalse(config["oxide_layer_in_model"])
            self.assertEqual(config["max_solver_seconds"], 21600)

    def test_production_reuse_is_exact_and_local_to_current_validation(self):
        aliases = [task for task in self.tasks.values() if "reused_from" in task]
        self.assertEqual(len(aliases), 8)
        for alias in aliases:
            source = self.tasks[alias["reused_from"]]
            self.assertEqual(alias["phase"], "pilot")
            self.assertEqual(source["phase"], "validation")
            self.assertEqual(alias["depends_on"], [source["id"]])
            self.assertEqual(alias["case_fingerprint"], source["case_fingerprint"])
            self.assertEqual(case_fingerprint(alias["case_path"], "fixed-test-implementation"),
                             case_fingerprint(source["case_path"], "fixed-test-implementation"))
            left = {key: value for key, value in self.configs[alias["id"]].items() if key != "id"}
            right = {key: value for key, value in self.configs[source["id"]].items() if key != "id"}
            self.assertEqual(left, right)
        for task in self.tasks.values():
            config = self.configs[task["id"]]
            self.assertEqual(task["case_fingerprint"], digest_object({key: value for key, value in config.items() if key != "id"}))

    def test_physics_or_strict_stop_change_invalidates_reuse(self):
        original_task = next(task for task in self.tasks.values()
                             if task["geometry_id"] == "baseline" and task["resolution"] == 20
                             and task["phase"] == "validation" and not task["strict_stop"])
        original = self.configs[original_task["id"]]
        baseline_hash = case_fingerprint(original_task["case_path"], "fixed-test-implementation")
        for key, value in (("resolution", 24), ("intensity_decay", 1e-6),
                           ("emission_theta_deg", -30), ("polarization", "s"),
                           ("courant", 0.2), ("wavelength_um", 10.4)):
            changed = {**original, key: value}
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "changed.json"
                atomic_json(path, changed)
                with self.subTest(key=key):
                    self.assertNotEqual(case_fingerprint(path, "fixed-test-implementation"), baseline_hash)
        self.assertNotEqual(case_fingerprint(original_task["case_path"], "different-implementation"), baseline_hash)

    def test_mesh_mirror_strict_and_boundary_comparisons_cover_validation(self):
        comparisons = self.campaign["gate_comparisons"]
        counts = {kind: sum(item["kind"] == kind for item in comparisons)
                  for kind in ("mesh", "mirror", "strict", "mesh_boundary")}
        self.assertEqual(counts, {"mesh": 4, "mirror": 8, "strict": 1, "mesh_boundary": 4})
        self.assertEqual({task_id for item in comparisons for task_id in item["ids"]},
                         set(self.campaign["production_gate_task_ids"]))
        for item in comparisons:
            self.assertEqual(item["metrics"], ["R", "T", "A_flux"])
            self.assertEqual(item["tolerance"], 0.001 if item["kind"] == "strict" else 0.01)

    def test_legacy_or_time_limited_results_never_qualify_production(self):
        self.assertIn("TIME_LIMIT_UNQUALIFIED", self.campaign["legacy_status"])
        self.assertIn("no legacy cache reuse", self.campaign["legacy_status"])
        synthetic_old = {task_id: {"status": "TIME_LIMIT_UNQUALIFIED", "metrics": {"R": 0.3, "T": 0.1, "A_flux": 0.6}}
                         for task_id in self.campaign["production_gate_task_ids"]}
        self.assertFalse(aggregate_gate(self.campaign, synthetic_old)["passed"])
        self.assertFalse(aggregate_gate(self.campaign, {})["passed"])
        for task in self.tasks.values():
            self.assertNotIn("D15", task["id"])
            self.assertNotIn("status", task)
            self.assertNotIn("qualified", task)

    def test_failed_material_or_geometry_audit_cannot_lock_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "bad_material"
            with patch.object(manifest, "validate_fit", return_value={"valid": False}):
                with self.assertRaisesRegex(ValueError, "MATERIAL_RECHECK_FAILED"):
                    manifest.create(run)
            self.assertFalse((run / "manifest.json").exists())
            self.assertFalse((run / "cases").exists())
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "bad_geometry"
            with patch.object(manifest, "validate_geometry", return_value={"valid": False}):
                with self.assertRaisesRegex(ValueError, "GEOMETRY_FAILED"):
                    manifest.create(run)
            self.assertFalse((run / "manifest.json").exists())
            self.assertFalse((run / "cases").exists())

    def test_declared_budgets_are_bounded(self):
        self.assertEqual(self.campaign["budgets"], {
            "total_active_seconds": 172800, "single_solver_seconds": 21600,
            "max_parallel": 4, "memory_fraction": 0.7, "mpi_ranks": 4})
        self.assertEqual(self.campaign["approved_max_targeted_repair_rounds"], 2)
        self.assertEqual(self.campaign["repairs_performed"], 0)


class CacheAdmissionTests(unittest.TestCase):
    def test_unsealed_or_time_limited_evidence_cannot_be_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            result = {"status": "TIME_LIMIT_UNQUALIFIED", "stop": {"converged": False}}
            atomic_json(path / "result.json", result)
            self.assertFalse(seal_result(path, "expected-case"))
            self.assertIsNone(reusable_result(path, "expected-case"))
            self.assertFalse((path / "completion.json").exists())

    def test_reuse_requires_exact_fingerprint_and_all_evidence_checksums(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            result = {"status": "QUALIFIED", "stop": {"converged": True, "stop_reason": "CONVERGED"},
                      "case_config": {"case": "structure"}, "metrics": {
                "R": 0.3, "T": 0.1, "A_flux": 0.6, "A_vol": 0.6,
                "order_R": 0.3, "order_T": 0.1, "reference_pseudo_reflection": 0.0}}
            atomic_json(path / "result.json", result)
            atomic_npz(path / "evidence.npz", E=np.ones((3, 4), dtype=complex))
            self.assertIsNone(reusable_result(path, "exact-case"))
            self.assertTrue(seal_result(path, "exact-case"))
            self.assertIsNotNone(reusable_result(path, "exact-case"))
            self.assertIsNone(reusable_result(path, "different-case"))
            atomic_npz(path / "evidence.npz", E=np.zeros((3, 4), dtype=complex))
            self.assertIsNone(reusable_result(path, "exact-case"))


class SolverGateTests(unittest.TestCase):
    def metrics(self):
        return {"R": 0.3, "T": 0.1, "A_flux": 0.6, "A_vol": 0.6,
                "order_R": 0.3, "order_T": 0.1, "reference_pseudo_reflection": 0.0}

    def test_layout_places_sources_and_monitors_outside_titanium_and_pml(self):
        config = {"geometry": {"period_um": 50, "thickness_um": 60},
                  "pml_um": 6, "air_above_um": 12, "air_below_um": 12}
        positions = layout(config)
        self.assertEqual(positions["height"], 96)
        self.assertEqual(positions["surface"] - positions["bottom"], 60)
        self.assertGreater(positions["top"], positions["source"])
        self.assertGreater(positions["source"], positions["refl"])
        self.assertGreater(positions["refl"], positions["surface"])
        self.assertGreater(positions["bottom"], positions["trans"])
        self.assertGreater(positions["trans"], positions["non_pml_bottom"])

    def test_independent_metrics_allow_consistent_converged_structure(self):
        checks = result_checks({"case": "structure"}, self.metrics(), True)
        self.assertTrue(all(checks.values()))

    def test_unconverged_timeout_cannot_pass_even_when_flux_numbers_match(self):
        checks = result_checks({"case": "structure"}, self.metrics(), False)
        self.assertFalse(checks["stopped"])
        self.assertFalse(all(checks.values()))

    def test_timeout_or_smoke_reason_overrides_final_collect_convergence(self):
        for reason in ("TIME_LIMIT_UNQUALIFIED", "SMOKE_ONLY", "CRASHED"):
            with self.subTest(reason=reason):
                checks = result_checks({"case": "structure"}, self.metrics(), True, stop_reason=reason)
                self.assertFalse(checks["stopped"])
                self.assertFalse(all(checks.values()))

    def test_flat_requires_independent_fresnel_evidence(self):
        checks = result_checks({"case": "flat"}, self.metrics(), True, analytical=None)
        self.assertFalse(all(checks.values()))
        consistent = result_checks({"case": "flat"}, self.metrics(), True,
                                   analytical={"R": 0.3, "T": 0.1, "A": 0.6})
        self.assertTrue(all(consistent.values()))

    def test_unknown_case_is_rejected_instead_of_becoming_an_unchecked_flat(self):
        checks = result_checks({"case": "strcuture"}, self.metrics(), True)
        self.assertFalse(checks["case_known"])
        self.assertFalse(all(checks.values()))

    def test_each_independent_physics_failure_blocks_qualification(self):
        changes = {"A_vol": 0.65, "order_R": 0.35, "order_T": 0.15,
                   "reference_pseudo_reflection": 0.002,
                   "R": -0.01, "T": 1.01, "A_flux": 1.01}
        for metric, value in changes.items():
            changed = {**self.metrics(), metric: value}
            with self.subTest(metric=metric):
                self.assertFalse(all(result_checks({"case": "structure"}, changed, True).values()))

    def test_nonfinite_metrics_block_qualification(self):
        for metric in self.metrics():
            for value in (float("nan"), float("inf")):
                with self.subTest(metric=metric, value=value):
                    checks = result_checks({"case": "structure"}, {**self.metrics(), metric: value}, True)
                    self.assertFalse(all(checks.values()))

    def test_planar_reflection_has_tighter_fresnel_gate(self):
        checks = result_checks({"case": "flat"}, self.metrics(), True,
                               analytical={"R": 0.303, "T": 0.1, "A": 0.597})
        self.assertFalse(checks["fresnel_R"])
        self.assertTrue(checks["fresnel_RTA"])


if __name__ == "__main__":
    unittest.main()
