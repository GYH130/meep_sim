"""Frozen-F and off-only continuation tests; subprocesses are always mocked."""
from copy import deepcopy
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch

import run_slot_stability_diagnostic as runner
import run_slot_stability_off as continuation
import test_slot_stability_runner as runner_tests
from ti2d.common import atomic_json, digest_file, digest_object
from ti2d.queue import read_json


class SlotStabilityOffTests(unittest.TestCase):
    def setUp(self):
        # Reuse the existing synthetic frozen-physics fixture, not its solver.
        self.fixture = runner_tests.SlotStabilityRunnerTests(methodName="runTest")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.stage, self.old = self.fixture.stage, self.fixture.run
        self.old_plan = runner.prepare(self.stage)
        self.run = self.stage / "repair_diagnostics" / continuation.RUN_ID
        shutil.copyfile(continuation.__file__, self.stage / Path(continuation.__file__).name)
        old_out = self.old / "averaging_on"
        old_out.mkdir()
        result = runner_tests.synthetic_probe_result(old_out, self.old_plan, "averaging_on", "EARLY_GROWTH_DETECTED")
        times = list(range(0, 21, 2)) + [20 + .25 / 32]
        ratios = [-40, -35, -30, -25, -20, -15, -5, 9,
                  18.756162058621562, 27.911704157091428, 37.067246255553826, 37.103]
        samples = [{"simulation_time": t, "envelope": {"all_fields_finite": True,
                    "log10_intensity_over_reference_peak": ratio},
                    "early_guard_status": "EARLY_GROWTH_DETECTED" if i >= 10 else None}
                   for i, (t, ratio) in enumerate(zip(times, ratios))]
        trace = read_json(old_out / "trace.json")
        trace["samples"] = samples
        atomic_json(old_out / "trace.json", trace)
        result.update(t_reached=times[-1], trace_summary={"sample_count": len(samples), "last_sample": samples[-1]})
        result["evidence_sha256"]["trace.json"] = digest_file(old_out / "trace.json")
        atomic_json(old_out / "result.json", result)
        old_duration = continuation.APPROVED_CURRENT_ACTIVE_SECONDS - runner.BASELINE_ACTIVE_SECONDS
        ledger = runner.closed_ledger(self.stage)
        ledger.add_interval(40000, 40000 + old_duration, 4,
                            continuation.SOURCE_RUN_ID + ":averaging_on", "CRASHED", attempt=1)
        atomic_json(old_out / "process.json", {"case_id": "averaging_on", "status": "CRASHED",
                    "qualification": False, "exit_code": 0, "pid": 987654321, "elapsed_s": old_duration,
                    "error": "RuntimeError: GROWTH_EVIDENCE_MISSING", "evidence_validated": False})
        atomic_json(self.old / "status.json", {"status": "FINISHED_NOT_QUALIFICATION", "qualification": False})
        self.old_hashes = continuation._file_hashes(self.old)
        self.ledger_before = deepcopy(ledger.data)
        for name, value in (("APPROVED_F_PLAN_SHA256", digest_file(self.old / "plan.json")),
                            ("APPROVED_F_ON_RESULT_SHA256", digest_file(old_out / "result.json")),
                            ("APPROVED_LEDGER_SHA256", digest_file(ledger.path)),
                            ("SERVER_STAGE", self.stage)):
            mock = patch.object(continuation, name, value)
            mock.start()
            self.addCleanup(mock.stop)
        mock = patch.object(runner, "session_live_members", return_value=[])
        mock.start()
        self.addCleanup(mock.stop)

    def prepare(self):
        self.plan = continuation.prepare(self.stage)
        return self.plan

    def verify(self):
        return continuation.verify_prepared(self.stage, read_json(self.run / "plan.json"),
                 read_json(self.run / "authorization.json"), read_json(self.run / "offline_reaudit.json"),
                 runner.closed_ledger(self.stage))

    def execute(self, status="SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION", duration=1.0):
        clock, children = runner_tests.FakeClock(), []
        def popen(command, **kwargs):
            self.assertEqual(command[command.index("--case-id") + 1], "averaging_off")
            self.assertEqual(command[command.index("--plan") + 1], str(self.run / "plan.json"))
            self.assertTrue(kwargs["start_new_session"])
            self.assertLessEqual(float(kwargs["env"]["TI2D_DIAGNOSTIC_TIMEOUT_SECONDS"]), 1170)
            child = runner_tests.FakeProcess(clock, duration, self.run / "averaging_off", self.plan,
                                             "averaging_off", status)
            children.append(child)
            return child
        def members(session_id):
            if session_id == 987654321:
                return []
            child = children[0]
            return [] if child.returncode is not None else [{"pid": child.pid, "session": child.pid,
                "pgid": child.pid, "state": "R", "start_ticks": 123}]
        resources = {"cpu_quota": 8, "memory_available_bytes": 4 * 1024**3,
                     "memory_limit_bytes": 8 * 1024**3, "memory_current_bytes": 4 * 1024**3}
        with patch.object(continuation.sys, "platform", "linux"), \
                patch.object(runner.time, "time", clock.time), \
                patch.object(runner.time, "monotonic", clock.monotonic), \
                patch.object(runner.time, "sleep", clock.sleep), \
                patch.object(runner, "cgroup_resources", return_value=resources), \
                patch.object(runner, "process_tree_rss", return_value=1024), \
                patch.object(runner, "session_live_members", side_effect=members), \
                patch.object(runner.subprocess, "Popen", side_effect=popen) as launcher:
            report = continuation.run(self.stage)
        return report, children, launcher

    def test_prepare_reaudits_first_trigger_and_preserves_entire_f_and_ledger(self):
        with patch.object(runner.subprocess, "Popen", side_effect=AssertionError("prepare must not launch")):
            plan = self.prepare()
        audit = read_json(self.run / "offline_reaudit.json")
        self.assertEqual(audit["growth_trigger_evidence"]["witness_times"], [16, 18, 20])
        self.assertEqual(audit["offline_reaudit_status"], "EARLY_GROWTH_DETECTED")
        self.assertEqual(audit["original_ledger_status_preserved"], "CRASHED")
        self.assertEqual(plan["execute_case_ids"], ["averaging_off"])
        self.assertEqual(plan["cases"], self.old_plan["cases"])
        self.assertTrue(plan["no_new_on_run"])
        self.assertFalse(plan["qualification"])
        self.assertEqual(self.old_hashes, continuation._file_hashes(self.old))
        self.assertEqual(self.ledger_before, runner.closed_ledger(self.stage).data)
        self.verify()

    def test_original_phase_budget_is_not_reset(self):
        plan = self.prepare()
        usage = read_json(self.run / "status.json")["budget"]
        self.assertEqual(plan["budgets"]["baseline_active_seconds"], runner.BASELINE_ACTIVE_SECONDS)
        self.assertEqual(plan["budgets"]["phase_seconds"], 1800)
        self.assertEqual(plan["budgets"]["single_seconds"], 1200)
        self.assertAlmostEqual(usage["phase_active_solver_seconds"],
                               continuation.APPROVED_CURRENT_ACTIVE_SECONDS - runner.BASELINE_ACTIVE_SECONDS)
        self.assertLess(usage["phase_remaining_solver_seconds"], 1560)

    def test_other_closed_work_further_reduces_off_allowance(self):
        self.prepare()
        ledger = runner.closed_ledger(self.stage)
        ledger.add_interval(45000, 46500, 4, "other_authorized_work", "CRASHED")
        _, _, launcher = self.execute()
        env = launcher.call_args.kwargs["env"]
        remaining = 1800 - (continuation.APPROVED_CURRENT_ACTIVE_SECONDS - runner.BASELINE_ACTIVE_SECONDS) - 1500
        self.assertAlmostEqual(float(env["TI2D_DIAGNOSTIC_TIMEOUT_SECONDS"]), remaining - 30)
        self.assertLess(float(env["TI2D_DIAGNOSTIC_TIMEOUT_SECONDS"]), 30)

    def test_prior_g_attempt_blocks_launch_even_without_output_directory(self):
        self.prepare()
        ledger = runner.closed_ledger(self.stage)
        ledger.add_interval(45000, 45001, 4, continuation.RUN_ID + ":averaging_off", "CRASHED")
        with self.assertRaisesRegex(RuntimeError, "ALREADY_ATTEMPTED_NO_RESUME"):
            self.verify()

    def test_run_launches_only_off_once_and_joint_report_stays_unqualified(self):
        self.prepare()
        report, children, launcher = self.execute(duration=2)
        self.assertEqual(launcher.call_count, 1)
        self.assertEqual(len(children), 1)
        self.assertFalse((self.run / "averaging_on").exists())
        self.assertEqual(report["run_id"], continuation.RUN_ID)
        self.assertTrue(report["both_arms_complete_evidence"])
        self.assertEqual(report["interpretation"], "AVERAGING_OFF_SHORT_WINDOW_STABILITY_ONLY")
        self.assertTrue(report["no_new_on_run"])
        self.assertFalse(report["qualification"])
        self.assertFalse(report["formal_R_T_A_conclusion"])
        self.assertEqual(report["arms"][0]["source_run_id"], continuation.SOURCE_RUN_ID)
        self.assertEqual(self.old_hashes, continuation._file_hashes(self.old))
        current = runner.closed_ledger(self.stage).data
        self.assertEqual(current["intervals"][:-1], self.ledger_before["intervals"])
        self.assertEqual(current["intervals"][-1]["task_id"], continuation.RUN_ID + ":averaging_off")
        self.assertEqual(current["intervals"][-2]["status"], "CRASHED")
        with patch.object(continuation.sys, "platform", "linux"), \
                patch.object(runner.subprocess, "Popen", side_effect=AssertionError("no retry")):
            with self.assertRaisesRegex(RuntimeError, "ALREADY_ATTEMPTED|ALREADY_STARTED"):
                continuation.run(self.stage)

    def test_crashed_or_falsely_qualified_probe_cannot_be_promoted(self):
        self.prepare()
        report, _, launcher = self.execute(status="QUALIFIED")
        self.assertEqual(launcher.call_count, 1)
        self.assertFalse(report["both_arms_complete_evidence"])
        self.assertFalse(report["qualification"])
        self.assertEqual(report["arms"][1]["status"], "CRASHED")
        self.assertEqual(report["interpretation"], "INCOMPLETE_EVIDENCE_NO_ATTRIBUTION")

    def test_first_trigger_must_be_valid_even_if_later_window_grows(self):
        path = self.old / "averaging_on/trace.json"
        trace = read_json(path)
        trace["samples"][10]["envelope"]["log10_intensity_over_reference_peak"] = 28
        atomic_json(path, trace)
        with self.assertRaisesRegex(RuntimeError, "EVIDENCE_HASH_MISMATCH"):
            self.prepare()

    def test_old_off_directory_or_attempt_disallows_continuation(self):
        (self.old / "averaging_off").mkdir()
        with self.assertRaisesRegex(RuntimeError, "OFF_ALREADY_ATTEMPTED"):
            self.prepare()

    def test_live_original_session_disallows_prepare(self):
        with patch.object(runner, "session_live_members", return_value=[{"pid": 42}]):
            with self.assertRaisesRegex(RuntimeError, "SESSION_NOT_CONFIRMED_STOPPED"):
                self.prepare()

    def test_open_ledger_disallows_prepare(self):
        ledger = runner.closed_ledger(self.stage)
        ledger.add_interval(45000, None, 4, "unfinished", "RUNNING")
        with self.assertRaisesRegex(RuntimeError, "UNCLOSED_PRIOR_INTERVAL"):
            self.prepare()

    def test_historical_status_and_source_log_changes_fail_closed(self):
        self.prepare()
        ledger = runner.closed_ledger(self.stage)
        ledger.data["intervals"][-1]["status"] = "EARLY_GROWTH_DETECTED"
        atomic_json(ledger.path, ledger.data)
        with self.assertRaisesRegex(RuntimeError, "HISTORICAL_LEDGER_CHANGED"):
            self.verify()
        atomic_json(ledger.path, self.ledger_before)
        (self.old / "extra_unapproved.txt").write_text("changed after prepare")
        with self.assertRaisesRegex(RuntimeError, "F_EVIDENCE_OR_OFFLINE_REAUDIT_CHANGED"):
            self.verify()

    def test_no_on_authorization_and_no_budget_reset_even_with_rehashed_plan(self):
        self.prepare()
        original = read_json(self.run / "plan.json")
        for key, value in (("execute_case_ids", ["averaging_on", "averaging_off"]),
                           ("budgets", dict(original["budgets"], baseline_active_seconds=continuation.APPROVED_CURRENT_ACTIVE_SECONDS))):
            changed = deepcopy(original)
            changed[key] = value
            changed["plan_sha256"] = digest_object({name: item for name, item in changed.items() if name != "plan_sha256"})
            atomic_json(self.run / "plan.json", changed)
            with self.assertRaisesRegex(RuntimeError, "OFF_ONLY_PLAN_CONTRACT_CHANGED"):
                self.verify()

    def test_frozen_continuation_or_runner_change_disallows_run(self):
        self.prepare()
        (self.stage / Path(continuation.__file__).name).write_text("modified")
        with self.assertRaisesRegex(RuntimeError, "IMPLEMENTATION_CHANGED"):
            self.verify()

    def test_local_run_refused_before_launch(self):
        with patch.object(continuation.sys, "platform", "darwin"), \
                patch.object(runner.subprocess, "Popen", side_effect=AssertionError("no local FDTD")):
            with self.assertRaisesRegex(RuntimeError, "SERVER_ONLY"):
                continuation.run(self.stage)

    def test_real_f_snapshot_reaudits_without_meep_when_available(self):
        source = Path(runner.__file__).parent / "delivered" / continuation.SOURCE_RUN_ID
        if not (source / "plan.json").is_file():
            self.skipTest("real F read-only download not present")
        result = runner.validate_probe_result(source / "averaging_on", read_json(source / "plan.json"), "averaging_on")
        self.assertEqual(result["growth_trigger_evidence"]["witness_times"], [16.0, 18.0, 20.0])
        self.assertEqual(result["t_reached"], 20.0078125)
        self.assertFalse(result["qualification"])


if __name__ == "__main__":
    unittest.main()
