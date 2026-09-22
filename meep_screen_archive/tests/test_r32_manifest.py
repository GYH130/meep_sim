"""Temporary synthetic artifacts only; never start Meep, MPI, or SSH."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import create_r32_run as r32
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest
from ti2d.queue import FileLock, SolverLedger, aggregate_gate, read_json, verify_manifest


class R32ManifestTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.stage = Path(temporary.name).resolve()
        self.source = self.stage / "runs" / r32.SOURCE_RUN_ID
        self.run = self.stage / "runs" / r32.RUN_ID
        self.case_path = self.source / "cases" / (r32.SOURCE_TASK_ID + ".json")
        self.result_path = self.source / "tasks" / r32.SOURCE_TASK_ID / "attempt_1" / "result.json"
        self.material = self.stage / "material.json"
        atomic_json(self.material, {"model": "synthetic frozen material"})
        self.implementation = implementation_digest()
        self.case = {"id": r32.SOURCE_TASK_ID, "case": "flat", "resolution": 16, "polarization": "p",
                     "emission_theta_deg": 30, "mpi_ranks": 4, "max_solver_seconds": 1800.0,
                     "geometry": {"period_um": 50, "thickness_um": 60}, "courant": .25,
                     "wavelength_um": 10.5, "pml_um": 6, "intensity_decay": 1e-4,
                     "stop_window_um": 20, "stop_consecutive_windows": 3, "flux_drift_tolerance": .001,
                     "material_path": str(self.material), "material_sha256": digest_file(self.material),
                     "implementation_sha256": self.implementation}
        self.manifest = {"schema_version": 1, "implementation_sha256": self.implementation,
                         "tasks": [{"id": r32.SOURCE_TASK_ID, "case_path": str(self.case_path),
                                    "phase": "validation", "depends_on": [], "geometry_id": "baseline",
                                    "resolution": 16, "gate_role": "flat"}],
                         "repair_budget": {"repair_active_solver_cap_seconds": 1800},
                         "repair_usage_at_creation": {"repair_active_solver_seconds": 17},
                         "repairs_performed": 1}
        self.result = {"status": "NUMERICAL_FAILED", "case_config": deepcopy(self.case),
                       "checks": {"stopped": True, "fresnel_R": False, "finite": True, "fresnel_RTA": True}}
        self.save_source()
        self.ledger = SolverLedger(self.stage / "budget_ledger.json")
        self.ledger.add_interval(100, 100 + r32.BASELINE_ACTIVE_SECONDS, 4, "historical", "NUMERICAL_FAILED")
        atomic_json(self.stage / "repair1.json", {"immutable_previous_budget": 1800})
        # Only the fixture's approval hashes differ from production; the validation path is unchanged.
        for name, value in (("APPROVED_SOURCE_CASE_SHA256", digest_file(self.case_path)),
                            ("APPROVED_LEDGER_SHA256", digest_file(self.ledger.path))):
            mock = patch.object(r32, name, value)
            mock.start()
            self.addCleanup(mock.stop)

    def save_source(self):
        atomic_json(self.case_path, self.case)
        self.manifest["tasks"][0]["case_fingerprint"] = digest_object(
            {key: value for key, value in self.case.items() if key != "id"})
        self.manifest.pop("manifest_sha256", None)
        self.manifest["manifest_sha256"] = digest_object(self.manifest)
        atomic_json(self.source / "manifest.json", self.manifest)
        atomic_json(self.result_path, self.result)

    def test_exact_diff_single_condition_and_separate_three_hour_shared_cap(self):
        originals = {path: path.read_bytes() for path in self.source.rglob("*") if path.is_file()}
        ledger_bytes = self.ledger.path.read_bytes()
        repair_bytes = (self.stage / "repair1.json").read_bytes()
        manifest = r32.create(self.stage)
        verify_manifest(manifest, self.run, self.implementation)
        self.assertEqual([task["id"] for task in manifest["tasks"]], [r32.TASK_ID])
        task = manifest["tasks"][0]
        self.assertEqual(task["phase"], "validation")
        self.assertEqual(task["depends_on"], [])
        self.assertEqual(task["resolution"], 32)
        self.assertNotIn("reused_from", task)
        self.assertNotIn("concurrency_probe", manifest)
        self.assertEqual(manifest["gate_comparisons"], [])
        self.assertEqual(manifest["pilot_design_conditions"], 0)
        self.assertEqual(manifest["budgets"]["max_parallel"], 1)
        self.assertEqual(manifest["budgets"]["single_solver_seconds"], 10800)
        self.assertEqual(manifest["budgets"]["stage_total_active_seconds"], 172800)
        self.assertEqual(manifest["budgets"]["total_active_seconds"], r32.BASELINE_ACTIVE_SECONDS + 10800)
        self.assertEqual(manifest["ledger_path"], str(self.ledger.path))
        for key in ("repair_budget", "repair_usage_at_creation", "repairs_performed"):
            self.assertNotIn(key, manifest)
        new_case = read_json(task["case_path"])
        self.assertEqual(new_case, dict(self.case, id=r32.TASK_ID, resolution=32, max_solver_seconds=10800.0))
        self.assertEqual(new_case["implementation_sha256"], self.implementation)
        diff = read_json(self.run / "config_diff.json")
        self.assertEqual(set(diff), {"id", "resolution", "max_solver_seconds"})
        authorization = read_json(self.run / "authorization.json")
        self.assertEqual(authorization["case_diff"], diff)
        self.assertEqual(authorization["ledger_snapshot"], self.ledger.data)
        self.assertEqual(authorization["source_case_sha256"], digest_file(self.case_path))
        self.assertEqual(authorization["source_manifest_sha256"], digest_file(self.source / "manifest.json"))
        self.assertEqual(manifest["authorization_sha256"], digest_object(authorization))
        self.assertFalse(authorization["automatic_retry_authorized"])
        self.assertEqual(originals, {path: path.read_bytes() for path in originals})
        self.assertEqual(ledger_bytes, self.ledger.path.read_bytes())
        self.assertEqual(repair_bytes, (self.stage / "repair1.json").read_bytes())
        self.assertEqual(len(list((self.run / "cases").glob("*.json"))), 1)
        self.assertFalse(aggregate_gate(manifest, {r32.TASK_ID: {"status": "QUALIFIED"}})["passed"])
        with self.assertRaises(FileExistsError):
            r32.create(self.stage)

    def test_missing_or_open_ledger_does_not_create_new_run(self):
        self.ledger.path.unlink()
        with self.assertRaisesRegex(RuntimeError, "EXISTING_SHARED_LEDGER_REQUIRED"):
            r32.create(self.stage)
        self.ledger.add_interval(30000, None, 4, "active", "RUNNING")
        with self.assertRaisesRegex(RuntimeError, "UNCLOSED_PRIOR_INTERVAL"):
            r32.create(self.stage)
        self.assertFalse(self.run.exists())

    def test_changed_implementation_and_source_manifest_are_rejected(self):
        with patch.object(r32, "implementation_digest", return_value="modified solver"):
            with self.assertRaisesRegex(ValueError, "IMPLEMENTATION_CHANGED"):
                r32.create(self.stage)
        self.manifest["scope"] = "unhashed mutation"
        atomic_json(self.source / "manifest.json", self.manifest)
        with self.assertRaisesRegex(ValueError, "MANIFEST_PAYLOAD_FINGERPRINT"):
            r32.create(self.stage)
        self.assertFalse(self.run.exists())

    def test_material_bytes_are_rechecked(self):
        atomic_json(self.material, {"model": "changed"})
        with self.assertRaisesRegex(ValueError, "FROZEN_MATERIAL_CHANGED"):
            r32.create(self.stage)
        self.assertFalse(self.run.exists())

    def test_resealed_unexpected_source_physics_are_rejected(self):
        for update in ({"courant": .5}, {"resolution": 20}, {"unexpected_physics_switch": True}):
            original = deepcopy(self.case)
            with self.subTest(update=update):
                self.case.update(update)
                self.result["case_config"] = deepcopy(self.case)
                self.save_source()
                with self.assertRaisesRegex(ValueError, "UNEXPECTED_SOURCE_PHYSICS_OR_NUMERICS"):
                    r32.create(self.stage)
                self.assertFalse(self.run.exists())
            self.case = original

    def test_source_failure_must_be_only_fresnel_and_stopped(self):
        for changes in ({"status": "QUALIFIED"}, {"checks": {"stopped": False, "fresnel_R": False}},
                        {"checks": {"stopped": True, "fresnel_R": False, "finite": False}},
                        {"case_config": {}}):
            with self.subTest(changes=changes):
                atomic_json(self.result_path, dict(self.result, **changes))
                with self.assertRaisesRegex(ValueError, "SOURCE_RESULT_NOT_APPROVED_FRESNEL_FAILURE"):
                    r32.create(self.stage)
                self.assertFalse(self.run.exists())

    def test_different_baseline_and_rewritten_ledger_cannot_receive_fresh_allowance(self):
        self.ledger.data["intervals"][0]["end"] += 1
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "APPROVED_BASELINE_MISMATCH"):
            r32.create(self.stage)
        self.ledger.data["intervals"][0]["end"] -= 1
        self.ledger.data["intervals"][0]["status"] = "rewritten"
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "APPROVED_LEDGER_CHANGED"):
            r32.create(self.stage)

    def test_frozen_full_history_allows_only_appends_and_counts_failed_attempts(self):
        r32.create(self.stage)
        authorization = read_json(self.run / "authorization.json")
        self.ledger.add_interval(30000, 30125, 4, r32.TASK_ID, "CRASHED")
        with FileLock(self.stage / "budget_ledger.lock"):
            verified = r32.verify_authorization(self.stage, authorization)
        self.assertEqual(verified.total_active_seconds() - r32.BASELINE_ACTIVE_SECONDS, 125)
        self.ledger.data["intervals"][0]["status"] = "changed history with identical duration"
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "HISTORICAL_LEDGER_CHANGED"):
            r32.verify_authorization(self.stage, authorization)
        authorization["ledger_snapshot"]["intervals"][0]["status"] = "tampered snapshot"
        with self.assertRaisesRegex(RuntimeError, "FROZEN_LEDGER_SNAPSHOT_CHANGED"):
            r32.verify_authorization(self.stage, authorization)


if __name__ == "__main__":
    unittest.main()
