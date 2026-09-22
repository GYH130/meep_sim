"""Synthetic temporary artifacts only; no Meep, MPI, SSH, or production writes."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import create_baseline_r32_run as baseline
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest
from ti2d.queue import FileLock, SolverLedger, aggregate_gate, read_json, verify_manifest


class BaselineR32ManifestTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.stage = Path(temporary.name).resolve()
        self.source = self.stage / "runs" / baseline.SOURCE_RUN_ID
        self.run = self.stage / "runs" / baseline.RUN_ID
        self.case_path = self.source / "cases" / (baseline.SOURCE_TASK_ID + ".json")
        self.result_path = self.source / "tasks" / baseline.SOURCE_TASK_ID / "attempt_1" / "result.json"
        self.material = self.stage / "material.json"
        atomic_json(self.material, {"model": "synthetic frozen material"})
        self.implementation = implementation_digest()
        self.case = {
            "id": baseline.SOURCE_TASK_ID, "case": "flat", "geometry_id": "baseline",
            "resolution": 32, "polarization": "p", "emission_theta_deg": 30.0,
            "mpi_ranks": 4, "max_solver_seconds": 10800.0, "strict_stop": False,
            "geometry": {"period_um": 50.0, "thickness_um": 60.0, "tilt_deg": 30.0},
            "courant": .25, "wavelength_um": 10.5, "pml_um": 6.0, "intensity_decay": 1e-4,
            "stop_window_um": 20.0, "stop_consecutive_windows": 3, "flux_drift_tolerance": .001,
            "reference_intensity_floor_fraction": 1.0,
            "stop_normalization_policy": "max(own_historic_peak, reference_source_on_non_pml_peak)",
            "mpi_partition_policy": "fixed_binary_Y0_X0_ranks0123",
            "material_path": str(self.material), "material_sha256": digest_file(self.material),
            "implementation_sha256": self.implementation,
        }
        self.manifest = {
            "schema_version": 1, "run_id": baseline.SOURCE_RUN_ID,
            "implementation_sha256": self.implementation, "material_sha256": self.case["material_sha256"],
            "tasks": [{"id": baseline.SOURCE_TASK_ID, "case_path": str(self.case_path),
                       "phase": "validation", "depends_on": [], "geometry_id": "baseline",
                       "resolution": 32, "gate_role": "flat"}],
        }
        evidence = {}
        for name in baseline.SOURCE_EVIDENCE_NAMES:
            path = self.result_path.parent / name
            atomic_json(path, {"synthetic_evidence": name})
            evidence[name] = digest_file(path)
        self.result = {
            "id": baseline.SOURCE_TASK_ID, "status": "QUALIFIED", "case_config": deepcopy(self.case),
            "implementation_sha256": self.implementation,
            "checks": {name: True for name in baseline.SOURCE_CHECK_NAMES},
            "stop": {"converged": True, "stop_reason": "CONVERGED"},
            "metrics": {"R": .96, "T": 0., "A_flux": .04, "A_vol": .04,
                        "order_R": .96, "order_T": 0., "reference_pseudo_reflection": 1e-8},
            "analytical": {"R": .96, "T": 0., "A": .04}, "evidence_sha256": evidence,
        }
        self.save_source()
        self.ledger = SolverLedger(self.stage / "budget_ledger.json")
        self.ledger.add_interval(100, 100 + baseline.BASELINE_ACTIVE_SECONDS, 4,
                                 "historical_including_failed_attempts", "QUALIFIED")
        for name, value in (("APPROVED_SOURCE_CASE_SHA256", digest_file(self.case_path)),
                            ("APPROVED_LEDGER_SHA256", digest_file(self.ledger.path))):
            mock = patch.object(baseline, name, value)
            mock.start()
            self.addCleanup(mock.stop)

    def save_source(self):
        atomic_json(self.case_path, self.case)
        self.manifest["tasks"][0]["case_fingerprint"] = digest_object(
            {key: value for key, value in self.case.items() if key != "id"})
        self.manifest.pop("manifest_sha256", None)
        self.manifest["manifest_sha256"] = digest_object(self.manifest)
        atomic_json(self.source / "manifest.json", self.manifest)
        self.result["case_sha256"] = digest_file(self.case_path)
        atomic_json(self.result_path, self.result)

    def test_exact_three_field_diff_single_baseline_and_six_hour_shared_cap(self):
        originals = {path: path.read_bytes() for path in self.source.rglob("*") if path.is_file()}
        ledger_bytes = self.ledger.path.read_bytes()
        manifest = baseline.create(self.stage)
        verify_manifest(manifest, self.run, self.implementation)
        self.assertEqual([task["id"] for task in manifest["tasks"]], [baseline.TASK_ID])
        task = manifest["tasks"][0]
        self.assertEqual(task["phase"], "validation")
        self.assertEqual(task["gate_role"], "baseline")
        self.assertEqual(task["depends_on"], [])
        self.assertEqual(task["resolution"], 32)
        self.assertNotIn("reused_from", task)
        self.assertNotIn("concurrency_probe", manifest)
        self.assertEqual(manifest["production_gate_task_ids"], [baseline.TASK_ID])
        self.assertEqual(manifest["gate_comparisons"], [])
        self.assertEqual(manifest["pilot_design_conditions"], 0)
        self.assertEqual(manifest["budgets"]["max_parallel"], 1)
        self.assertEqual(manifest["budgets"]["single_solver_seconds"], 21600)
        self.assertEqual(manifest["budgets"]["stage_total_active_seconds"], 172800)
        self.assertEqual(manifest["budgets"]["total_active_seconds"], 51721.36579871178)
        self.assertEqual(manifest["ledger_path"], str(self.ledger.path))
        new_case = read_json(task["case_path"])
        self.assertEqual(new_case, dict(self.case, id=baseline.TASK_ID, case="structure", max_solver_seconds=21600.0))
        diff = read_json(self.run / "config_diff.json")
        self.assertEqual(diff, {"id": {"before": baseline.SOURCE_TASK_ID, "after": baseline.TASK_ID},
                                "case": {"before": "flat", "after": "structure"},
                                "max_solver_seconds": {"before": 10800., "after": 21600.}})
        authorization = read_json(self.run / "authorization.json")
        self.assertEqual(authorization["case_diff"], diff)
        self.assertEqual(authorization["ledger_snapshot"], self.ledger.data)
        self.assertEqual(authorization["source_case_sha256"], digest_file(self.case_path))
        self.assertEqual(authorization["source_manifest_sha256"], digest_file(self.source / "manifest.json"))
        self.assertEqual(authorization["source_result_sha256"], digest_file(self.result_path))
        self.assertEqual(authorization["source_evidence_sha256"], self.result["evidence_sha256"])
        self.assertEqual(manifest["authorization_sha256"], digest_object(authorization))
        self.assertFalse(authorization["automatic_retry_authorized"])
        self.assertEqual(originals, {path: path.read_bytes() for path in originals})
        self.assertEqual(ledger_bytes, self.ledger.path.read_bytes())
        self.assertEqual(len(list((self.run / "cases").glob("*.json"))), 1)
        self.assertFalse(aggregate_gate(manifest, {baseline.TASK_ID: {"status": "QUALIFIED"}})["passed"])

    def test_repeat_creation_is_rejected_without_overwrite(self):
        baseline.create(self.stage)
        created = {path: path.read_bytes() for path in self.run.rglob("*") if path.is_file()}
        with self.assertRaisesRegex(FileExistsError, "BASELINE_R32_RUN_ALREADY_EXISTS"):
            baseline.create(self.stage)
        self.assertEqual(created, {path: path.read_bytes() for path in created})

    def test_missing_or_open_ledger_rejected(self):
        self.ledger.path.unlink()
        with self.assertRaisesRegex(RuntimeError, "EXISTING_SHARED_LEDGER_REQUIRED"):
            baseline.create(self.stage)
        self.ledger.add_interval(40000, None, 4, "active", "RUNNING")
        with self.assertRaisesRegex(RuntimeError, "UNCLOSED_PRIOR_INTERVAL"):
            baseline.create(self.stage)
        self.assertFalse(self.run.exists())

    def test_changed_baseline_or_rewritten_history_rejected(self):
        self.ledger.data["intervals"][0]["end"] += 1
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "APPROVED_BASELINE_MISMATCH"):
            baseline.create(self.stage)
        self.ledger.data["intervals"][0]["end"] -= 1
        self.ledger.data["intervals"][0]["status"] = "rewritten"
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "APPROVED_LEDGER_CHANGED"):
            baseline.create(self.stage)
        self.assertFalse(self.run.exists())

    def test_failed_missing_false_or_nonboolean_source_checks_rejected(self):
        changes = ({"status": "NUMERICAL_FAILED"}, {"checks": {}},
                   {"checks": dict(self.result["checks"], stopped=False)},
                   {"checks": dict(self.result["checks"], fresnel_R=1)},
                   {"checks": dict(self.result["checks"], unexpected=False)},
                   {"case_config": {}}, {"case_sha256": "0" * 64},
                   {"implementation_sha256": "0" * 64})
        for update in changes:
            with self.subTest(update=update):
                atomic_json(self.result_path, dict(self.result, **update))
                with self.assertRaisesRegex(ValueError, "SOURCE_RESULT_NOT_APPROVED_QUALIFIED"):
                    baseline.create(self.stage)
                self.assertFalse(self.run.exists())

    def test_source_qualification_is_recomputed_from_metrics_and_stop(self):
        for update in ({"stop": {"converged": False}},
                       {"metrics": dict(self.result["metrics"], A_vol=.5)}):
            with self.subTest(update=update):
                atomic_json(self.result_path, dict(self.result, **update))
                with self.assertRaisesRegex(ValueError, "QUALIFIED"):
                    baseline.create(self.stage)
                self.assertFalse(self.run.exists())

    def test_all_eight_evidence_entries_required_and_rehashed(self):
        for change in ("empty", "missing_entry", "extra_entry", "bad_hash", "missing_file", "changed_file", "symlink"):
            with self.subTest(change=change):
                result = deepcopy(self.result)
                path = self.result_path.parent / "orders.json"
                original = path.read_bytes()
                if change == "empty":
                    result["evidence_sha256"] = {}
                elif change == "missing_entry":
                    result["evidence_sha256"].pop("orders.json")
                elif change == "extra_entry":
                    result["evidence_sha256"]["../outside.json"] = "0" * 64
                elif change == "bad_hash":
                    result["evidence_sha256"]["orders.json"] = "0" * 64
                elif change == "missing_file":
                    path.unlink()
                elif change == "changed_file":
                    atomic_json(path, {"changed": True})
                else:
                    path.unlink()
                    path.symlink_to(self.material)
                atomic_json(self.result_path, result)
                with self.assertRaisesRegex(ValueError, "SOURCE_EVIDENCE"):
                    baseline.create(self.stage)
                self.assertFalse(self.run.exists())
                if path.is_symlink():
                    path.unlink()
                atomic_json(path, {"synthetic_evidence": "orders.json"})
                self.assertEqual(path.read_bytes(), original)

    def test_resealed_out_of_scope_parameters_rejected(self):
        for update in ({"courant": .5}, {"resolution": 20}, {"polarization": "s"},
                       {"emission_theta_deg": -30}, {"wavelength_um": 10.6},
                       {"mpi_partition_policy": "different"}, {"unexpected_physics_switch": True}):
            original = deepcopy(self.case)
            with self.subTest(update=update):
                self.case.update(update)
                self.result["case_config"] = deepcopy(self.case)
                self.save_source()
                with self.assertRaisesRegex(ValueError, "UNEXPECTED_SOURCE_PHYSICS_OR_NUMERICS"):
                    baseline.create(self.stage)
                self.assertFalse(self.run.exists())
            self.case = original

    def test_changed_implementation_material_and_manifest_rejected(self):
        with patch.object(baseline, "implementation_digest", return_value="modified implementation"):
            with self.assertRaisesRegex(ValueError, "APPROVED_IMPLEMENTATION_CHANGED"):
                baseline.create(self.stage)
        atomic_json(self.material, {"model": "changed"})
        with self.assertRaisesRegex(ValueError, "FROZEN_MATERIAL_CHANGED"):
            baseline.create(self.stage)
        self.manifest["scope"] = "unhashed mutation"
        atomic_json(self.source / "manifest.json", self.manifest)
        with self.assertRaisesRegex(ValueError, "MANIFEST_PAYLOAD_FINGERPRINT"):
            baseline.create(self.stage)
        self.assertFalse(self.run.exists())

    def test_full_history_allows_appends_and_counts_failed_attempts(self):
        baseline.create(self.stage)
        authorization = read_json(self.run / "authorization.json")
        self.ledger.add_interval(40000, 40125, 4, baseline.TASK_ID, "CRASHED")
        with FileLock(self.stage / "budget_ledger.lock"):
            verified = baseline.verify_authorization(self.stage, authorization)
        self.assertEqual(verified.total_active_seconds() - baseline.BASELINE_ACTIVE_SECONDS, 125)
        self.ledger.data["intervals"][0]["status"] = "same duration rewritten history"
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "HISTORICAL_LEDGER_CHANGED"):
            baseline.verify_authorization(self.stage, authorization)
        authorization["ledger_snapshot"]["intervals"][0]["status"] = "tampered snapshot"
        with self.assertRaisesRegex(RuntimeError, "FROZEN_LEDGER_SNAPSHOT_CHANGED"):
            baseline.verify_authorization(self.stage, authorization)

    def test_authorization_cannot_expand_budget_or_enable_retry(self):
        baseline.create(self.stage)
        authorization = read_json(self.run / "authorization.json")
        for update in ({"phase_active_solver_cap_seconds": 21601},
                       {"effective_total_active_seconds": 172800},
                       {"baseline_active_solver_seconds": 0}, {"automatic_retry_authorized": True}):
            with self.subTest(update=update):
                with self.assertRaisesRegex(RuntimeError, "BASELINE_R32_AUTHORIZATION_CHANGED"):
                    baseline.verify_authorization(self.stage, dict(authorization, **update))


if __name__ == "__main__":
    unittest.main()
