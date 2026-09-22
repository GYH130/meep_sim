"""Pure temporary-artifact checks; no Meep, MPI, or server connection."""
from copy import deepcopy
import contextlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import Mock, patch

from create_repair_run import create
from repair_budget import (BASELINE_ACTIVE_SECONDS, SNAPSHOT_RELATIVE,
                           locked_repair_budget, repair_usage)
import run_repair_probe
from ti2d.common import digest_file, implementation_digest
from ti2d.queue import FileLock, Queue, SolverLedger, aggregate_gate, atomic_json


class RepairBudgetTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.stage = Path(self.temporary.name)
        self.ledger = SolverLedger(self.stage / "budget_ledger.json")
        self.ledger.add_interval(100, 100 + BASELINE_ACTIVE_SECONDS, 4, "historical", "NUMERICAL_FAILED")
        self.original = deepcopy(self.ledger.data)
        snapshot = self.stage / SNAPSHOT_RELATIVE
        snapshot.parent.mkdir(parents=True)
        with tarfile.open(snapshot, "w:gz") as archive:
            archive.add(self.ledger.path, arcname="budget_ledger.json")

    def freeze(self):
        with FileLock(self.stage / "budget_ledger.lock"):
            return locked_repair_budget(self.stage, self.ledger)

    def test_diagnostic_usage_consumes_same_frozen_allowance_without_reset(self):
        self.ledger.add_interval(30000, 30150, 4, "repair_probe", "TIME_LIMIT_UNQUALIFIED")
        budget = self.freeze()
        usage = repair_usage(self.ledger, budget)
        self.assertEqual(budget["effective_total_active_seconds"], BASELINE_ACTIVE_SECONDS + 1800)
        self.assertEqual(usage["repair_active_solver_seconds"], 150)
        self.assertEqual(usage["repair_remaining_solver_seconds"], 1650)
        self.assertEqual(usage["stage_cap_seconds"], 172800)
        before = (self.stage / "repair1.json").read_bytes()
        self.freeze()
        self.assertEqual(before, (self.stage / "repair1.json").read_bytes())
        self.assertEqual(self.ledger.data["intervals"][:1], self.original["intervals"])

    def test_open_interval_blocks_another_launch(self):
        self.ledger.add_interval(30000, None, 4, "prior_probe", "RUNNING")
        with self.assertRaisesRegex(RuntimeError, "UNCLOSED_PRIOR"):
            self.freeze()

    def test_modified_history_and_contract_are_rejected(self):
        self.freeze()
        self.ledger.data["intervals"][0]["status"] = "changed"
        with self.assertRaisesRegex(RuntimeError, "HISTORICAL_LEDGER_CHANGED"):
            self.freeze()
        self.ledger.data = deepcopy(self.original)
        path = self.stage / "repair1.json"
        changed = json.loads(path.read_text())
        changed["repair_active_solver_cap_seconds"] = 3600
        atomic_json(path, changed)
        with self.assertRaisesRegex(RuntimeError, "FROZEN_REPAIR_BUDGET_CHANGED"):
            self.freeze()

    def test_missing_ledger_does_not_get_a_fresh_budget(self):
        missing = SolverLedger(self.stage / "missing" / "budget_ledger.json")
        with self.assertRaisesRegex(RuntimeError, "EXISTING_SHARED_LEDGER_REQUIRED"):
            locked_repair_budget(self.stage / "missing", missing)

    def source_run(self):
        source = self.stage / "runs" / "stage_20260921_B"
        source.mkdir(parents=True)
        material = self.stage / "material.json"
        atomic_json(material, {"label": "frozen"})
        task_id = "flat_r16_p_plus30"
        case = {"id": task_id, "case": "flat", "resolution": 16, "polarization": "p",
                "emission_theta_deg": 30, "mpi_ranks": 4, "max_solver_seconds": 21600,
                "geometry": {}, "material_path": str(material), "material_sha256": digest_file(material),
                "implementation_sha256": "old-code", "stop_consecutive_windows": 3}
        atomic_json(source / "case.json", case)
        task = {"id": task_id, "case_path": "case.json", "phase": "validation",
                "depends_on": [], "gate_role": "flat", "geometry_id": "baseline"}
        atomic_json(source / "manifest.json", {"tasks": [task, {"id": "unauthorized_extra"}],
                    "concurrency_probe": {"task_ids": [task_id, "unauthorized_extra"]},
                    "gate_comparisons": [{"kind": "mesh", "ids": [task_id, "unauthorized_extra"]}],
                    "manifest_sha256": "old-manifest"})
        return source

    def test_builder_contains_only_one_task_and_does_not_qualify_campaign(self):
        source = self.source_run()
        originals = {p.name: p.read_bytes() for p in source.iterdir()}
        self.ledger.add_interval(30000, 30150, 4, "repair_probe", "TIME_LIMIT_UNQUALIFIED")
        manifest = create(self.stage, overrides={"stop_post_source_rta_drift": True})
        self.assertEqual([t["id"] for t in manifest["tasks"]], ["flat_r16_p_plus30"])
        self.assertNotIn("concurrency_probe", manifest)
        self.assertEqual(manifest["pilot_design_conditions"], 0)
        self.assertEqual(manifest["gate_comparisons"], [])
        self.assertEqual(manifest["implementation_sha256"], implementation_digest())
        self.assertEqual(manifest["budgets"]["total_active_seconds"], BASELINE_ACTIVE_SECONDS + 1800)
        self.assertEqual(manifest["repair_usage_at_creation"]["repair_remaining_solver_seconds"], 1650)
        self.assertFalse(aggregate_gate(manifest, {"flat_r16_p_plus30": {"status": "QUALIFIED"}})["passed"])
        self.assertEqual(originals, {p.name: p.read_bytes() for p in source.iterdir()})
        queue = Queue(Path(manifest["run_dir"]) / "manifest.json", phases=("validation",))
        self.assertEqual(queue.ledger.path.resolve(), self.ledger.path.resolve())
        self.assertEqual(len(queue.manifest["tasks"]), 1)
        with self.assertRaises(FileExistsError):
            create(self.stage)

    def test_builder_refuses_physics_changes(self):
        self.source_run()
        with self.assertRaisesRegex(ValueError, "PRESERVE_PHYSICS"):
            create(self.stage, overrides={"resolution": 20})
        self.assertFalse((self.stage / "runs" / "repair_20260922_C").exists())

    def invoke_mock_probe(self, **mocks):
        atomic_json(self.stage / "repair_flux_probe.py", {"synthetic": True})
        resources = {"cpu_quota": 4, "memory_available_bytes": 4 * 1024**3}
        with contextlib.ExitStack() as stack:
            stack.enter_context(patch.object(run_repair_probe, "_assert_server"))
            stack.enter_context(patch.object(run_repair_probe, "cgroup_resources", return_value=resources))
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            for target, config in mocks.items():
                stack.enter_context(patch(target, **config))
            return run_repair_probe.run_probe(self.stage)

    def test_probe_spawn_failure_closes_ledger(self):
        code = self.invoke_mock_probe(**{
            "run_repair_probe.subprocess.Popen": {"side_effect": OSError("synthetic launch failure")}})
        self.assertEqual(code, 1)
        ledger = SolverLedger(self.ledger.path)
        self.assertEqual(ledger.data["intervals"][-1]["status"], "CRASHED")
        self.assertIsNotNone(ledger.data["intervals"][-1]["end"])

    def test_probe_timeout_cannot_be_overwritten_by_exit_zero(self):
        child = Mock(pid=4321, returncode=None)
        child.poll.side_effect = lambda: child.returncode
        def terminate(process, grace):
            self.assertIs(process, child)
            child.returncode = 0
        code = self.invoke_mock_probe(**{
            "run_repair_probe.subprocess.Popen": {"return_value": child},
            "run_repair_probe.process_tree_rss": {"return_value": 1024},
            "run_repair_probe.time.monotonic": {"side_effect": [0, 151, 151]},
            "run_repair_probe.terminate_child": {"side_effect": terminate}})
        self.assertEqual(code, 1)
        self.assertEqual(SolverLedger(self.ledger.path).data["intervals"][-1]["status"], "TIME_LIMIT_UNQUALIFIED")


if __name__ == "__main__":
    unittest.main()
