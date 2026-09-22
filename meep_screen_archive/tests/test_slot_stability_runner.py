"""Synthetic artifacts and mocked children only: no SSH, Meep or MPI launch."""
from copy import deepcopy
from pathlib import Path
import shutil
import signal
import tempfile
import unittest
from unittest.mock import patch

import run_slot_stability_diagnostic as runner
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest
from ti2d.queue import FileLock, SolverLedger, read_json


class FakeClock:
    def __init__(self):
        self.elapsed = 0.0

    def time(self):
        return 50000.0 + self.elapsed

    def monotonic(self):
        return self.elapsed

    def sleep(self, amount):
        self.elapsed += amount


def synthetic_probe_result(output, plan, case_id, status):
    if status == "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION":
        times, ratios = list(range(0, 101, 2)), [0.0] * 51
    elif status == "EARLY_GROWTH_DETECTED":
        times, ratios = [0, 2, 4], [13.0, 15.0, 17.0]
    else:
        times, ratios = [0, 2], [0.0, None]
    samples = [{"simulation_time": t, "envelope": {
        "all_fields_finite": not (status == "NONFINITE_FIELD" and i == len(times) - 1),
        "log10_intensity_over_reference_peak": ratio,
    }, "early_guard_status": "EARLY_GROWTH_DETECTED" if status == "EARLY_GROWTH_DETECTED" and i == 2 else None}
        for i, (t, ratio) in enumerate(zip(times, ratios))]
    atomic_json(output / "trace.json", {"case_id": case_id,
                "eps_averaging": runner.ARMS[case_id], "qualification": False, "samples": samples})
    atomic_json(output / "geometry_metadata.json", {"synthetic": True})
    atomic_json(output / "sample_grid.npz", {"synthetic_not_npz": True})
    result = {"case_id": case_id, "eps_averaging": runner.ARMS[case_id], "qualification": False,
              "status": status, "t_reached": times[-1], "plan_sha256": plan["plan_sha256"],
              "physical_case_sha256": plan["source_case_sha256"],
              "trace_summary": {"sample_count": len(samples), "last_sample": samples[-1]},
              "evidence_sha256": {name: digest_file(output / name) for name in
                                  ("trace.json", "geometry_metadata.json", "sample_grid.npz")}}
    for key in ("source_case_sha256", "implementation_sha256", "material_sha256", "probe_sha256"):
        result[key] = plan[key]
    atomic_json(output / "result.json", result)
    atomic_json(output / "heartbeat.json", {"t_reached": times[-1]})
    return result


class FakeProcess:
    def __init__(self, clock, duration, output, plan, case_id, status, returncode=0):
        self.clock, self.duration, self.start = clock, duration, clock.elapsed
        self.output, self.plan, self.case_id, self.status = output, plan, case_id, status
        self.returncode, self.final_returncode, self.pid = None, returncode, 1234567

    def poll(self):
        if self.returncode is None and self.clock.elapsed - self.start >= self.duration:
            self.returncode = self.final_returncode
            if self.status is not None:
                synthetic_probe_result(self.output, self.plan, self.case_id, self.status)
        return self.returncode


class SlotStabilityRunnerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.stage = Path(temporary.name).resolve()
        self.source = self.stage / "runs" / runner.SOURCE_RUN_ID
        self.run = self.stage / "repair_diagnostics" / runner.RUN_ID
        self.material = self.stage / "material.json"
        atomic_json(self.material, {"synthetic_frozen_material": True})
        self.case_path = self.source / "cases" / (runner.SOURCE_TASK_ID + ".json")
        self.failure_path = self.source / "tasks" / runner.SOURCE_TASK_ID / "attempt_1" / "result.json"
        self.case = {
            "id": runner.SOURCE_TASK_ID, "case": "structure", "geometry_id": "baseline",
            "resolution": 32, "courant": .25, "polarization": "p", "emission_theta_deg": 30.0,
            "wavelength_um": 10.5, "mpi_ranks": 4,
            "mpi_partition_policy": "fixed_binary_Y0_X0_ranks0123",
            "implementation_sha256": implementation_digest(), "geometry": {"tilt_deg": 30},
            "material_path": str(self.material), "material_sha256": digest_file(self.material),
        }
        atomic_json(self.case_path, self.case)
        manifest = {"run_id": runner.SOURCE_RUN_ID, "implementation_sha256": implementation_digest(),
                    "material_sha256": self.case["material_sha256"],
                    "tasks": [{"id": runner.SOURCE_TASK_ID, "case_path": str(self.case_path),
                               "case_fingerprint": digest_object({k: v for k, v in self.case.items() if k != "id"})}]}
        manifest["manifest_sha256"] = digest_object(manifest)
        atomic_json(self.source / "manifest.json", manifest)
        atomic_json(self.failure_path, {"status": "CRASHED", "error": "envelope_power must be finite"})
        self.ledger = SolverLedger(self.stage / "budget_ledger.json")
        self.ledger.add_interval(10, 10 + runner.BASELINE_ACTIVE_SECONDS, 4, "history", "CRASHED")
        (self.stage / "diagnostic_tools").mkdir()
        source_stage = Path(runner.__file__).parent
        shutil.copyfile(source_stage / "diagnostic_tools/slot_stability_probe.py",
                        self.stage / "diagnostic_tools/slot_stability_probe.py")
        shutil.copyfile(runner.__file__, self.stage / Path(runner.__file__).name)
        for name, value in (
                ("APPROVED_SOURCE_CASE_SHA256", digest_file(self.case_path)),
                ("APPROVED_FAILURE_SHA256", digest_file(self.failure_path)),
                ("APPROVED_LEDGER_SHA256", digest_file(self.ledger.path)),
                ("APPROVED_IMPLEMENTATION_SHA256", implementation_digest()),
                ("SERVER_STAGE", self.stage)):
            mock = patch.object(runner, name, value)
            mock.start()
            self.addCleanup(mock.stop)
        self.clock = FakeClock()

    def prepared(self):
        self.plan = runner.prepare(self.stage)
        return self.plan

    def check_frozen(self):
        return runner.verify_frozen(self.stage, read_json(self.run / "plan.json"),
                                   read_json(self.run / "authorization.json"), runner.closed_ledger(self.stage))

    def execute(self, statuses=None, *, duration=1.0, rss=1024, spawn_error=None,
                callback=None, returncode=0, resources=None, cleanup_error=None,
                orphan_after_exit=False, term_resistant=False, direct_arm=None):
        statuses = statuses or ["EARLY_GROWTH_DETECTED", "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION"]
        created, orphan_alive = [], {}

        def popen(command, **kwargs):
            self.assertTrue(kwargs["start_new_session"])
            self.assertEqual(kwargs["cwd"], self.stage)
            self.assertEqual(command[:6], [runner.MPIEXEC, "-n", "4", runner.PYTHON, "-B", "-m"])
            self.assertEqual(command[6], "diagnostic_tools.slot_stability_probe")
            self.assertTrue(all(kwargs["env"][key] == "1" for key in runner.THREAD_ENV))
            self.assertLessEqual(float(kwargs["env"]["TI2D_DIAGNOSTIC_TIMEOUT_SECONDS"]), 1170)
            if spawn_error:
                self.clock.sleep(.1)
                raise spawn_error
            case_id = command[command.index("--case-id") + 1]
            child = FakeProcess(self.clock, duration, self.run / case_id, self.plan, case_id,
                                statuses[len(created)], returncode=returncode)
            child.pid += len(created)
            created.append(child)
            orphan_alive[child.pid] = orphan_after_exit
            if callback:
                callback(child)
            return child

        def members(session_id):
            child = next(c for c in created if c.pid == session_id)
            leader = ([{"pid": child.pid, "session": child.pid, "pgid": child.pid,
                        "state": "R", "start_ticks": 123}] if child.returncode is None else [])
            orphan = ([{"pid": child.pid + 100, "session": child.pid, "pgid": child.pid + 100,
                        "state": "R", "start_ticks": 124}]
                      if child.returncode is not None and orphan_alive[child.pid] else [])
            return leader + orphan

        def terminate(member, session_id, sig):
            child = next(c for c in created if c.pid == session_id)
            self.assertEqual(member["session"], child.pid)
            self.assertIn(member["pid"], (child.pid, child.pid + 100))
            self.clock.sleep(.25)
            if cleanup_error:
                raise cleanup_error
            if term_resistant and sig == signal.SIGTERM:
                return
            if member["pid"] == child.pid:
                child.returncode = -int(sig)
            else:
                orphan_alive[child.pid] = False

        resources = resources or {"cpu_quota": 8, "memory_available_bytes": 4 * 1024**3,
                                  "memory_limit_bytes": 8 * 1024**3, "memory_current_bytes": 4 * 1024**3}
        with patch.object(runner.sys, "platform", "linux"), \
                patch.object(runner.time, "time", self.clock.time), \
                patch.object(runner.time, "monotonic", self.clock.monotonic), \
                patch.object(runner.time, "sleep", self.clock.sleep), \
                patch.object(runner, "cgroup_resources", return_value=resources), \
                patch.object(runner, "process_tree_rss", return_value=rss), \
                patch.object(runner.subprocess, "Popen", side_effect=popen) as popen_mock, \
                patch.object(runner, "session_live_members", side_effect=members), \
                patch.object(runner, "signal_session_member", side_effect=terminate) as terminate_mock:
            if direct_arm is None:
                report = runner.run(self.stage)
            else:
                ledger = runner.closed_ledger(self.stage)
                arm = runner.run_arm(self.stage, self.plan, ledger, direct_arm, [False])
                report = runner.report_for(self.plan, [arm], ledger)
        return report, created, popen_mock, terminate_mock

    def test_prepare_freezes_all_input_bytes_and_exact_two_arm_contract_without_solver(self):
        original = {p: p.read_bytes() for p in self.source.rglob("*") if p.is_file()}
        ledger_bytes = self.ledger.path.read_bytes()
        with patch.object(runner.subprocess, "Popen", side_effect=AssertionError("no launch")):
            plan = self.prepared()
        self.assertEqual(plan["cases"], [{"id": "averaging_on", "eps_averaging": True},
                                         {"id": "averaging_off", "eps_averaging": False}])
        self.assertEqual(plan["diagnostic"]["stop_time"], 100)
        self.assertEqual(plan["diagnostic"]["reference_intensity_peak"], .3334312033062042)
        self.assertEqual(plan["budgets"]["phase_seconds"], 1800)
        self.assertEqual(plan["budgets"]["single_seconds"], 1200)
        self.assertEqual(plan["budgets"]["stage_total_seconds"], 172800)
        self.assertFalse(plan["qualification"])
        self.assertEqual(ledger_bytes, self.ledger.path.read_bytes())
        self.assertEqual(original, {path: path.read_bytes() for path in original})
        auth = read_json(self.run / "authorization.json")
        self.assertEqual(auth["ledger_snapshot"], self.ledger.data)
        for entry in auth["inputs"].values():
            self.assertEqual(Path(entry["source_path"]).read_bytes(), Path(entry["snapshot_path"]).read_bytes())
        self.check_frozen()

    def test_prepare_rejects_overwrite_missing_open_or_changed_ledger(self):
        self.prepared()
        with self.assertRaisesRegex(FileExistsError, "ALREADY_EXISTS"):
            runner.prepare(self.stage)
        # Independent fixtures keep the original failed guard side-effect free.
        self.ledger.path.unlink()
        with self.assertRaisesRegex(RuntimeError, "EXISTING_SHARED_LEDGER"):
            runner.closed_ledger(self.stage)
        self.ledger.add_interval(40000, None, 4, "unfinished", "RUNNING")
        with self.assertRaisesRegex(RuntimeError, "UNCLOSED"):
            runner.closed_ledger(self.stage)

    def test_changed_approved_baseline_refused_before_preparation(self):
        self.ledger.data["intervals"][0]["end"] += 1
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "BASELINE_MISMATCH"):
            runner.prepare(self.stage)
        self.assertFalse(self.run.exists())

    def test_failure_reason_must_be_exact_even_if_fixture_hash_resealed(self):
        atomic_json(self.failure_path, {"status": "CRASHED", "error": "different"})
        with patch.object(runner, "APPROVED_FAILURE_SHA256", digest_file(self.failure_path)):
            with self.assertRaisesRegex(RuntimeError, "EXPECTED_SOURCE_NUMERICAL_FAILURE"):
                runner.prepare(self.stage)

    def test_hash_or_source_material_change_is_refused(self):
        with patch.object(runner, "implementation_digest", return_value="changed"):
            with self.assertRaisesRegex(RuntimeError, "IMPLEMENTATION_CHANGED"):
                runner.prepare(self.stage)
        atomic_json(self.material, {"changed": True})
        with self.assertRaisesRegex(RuntimeError, "MATERIAL_CHANGED"):
            runner.prepare(self.stage)

    def test_rewritten_history_rejected_but_closed_appends_reduce_allowance(self):
        self.prepared()
        self.ledger.add_interval(40000, 40500, 4, "another_closed_task", "CRASHED")
        self.check_frozen()
        self.assertAlmostEqual(runner.budget_usage(self.ledger)["effective_remaining_solver_seconds"], 1300)
        self.ledger.data["intervals"][0]["status"] = "rewritten"
        atomic_json(self.ledger.path, self.ledger.data)
        with self.assertRaisesRegex(RuntimeError, "HISTORICAL_LEDGER_CHANGED"):
            self.check_frozen()

    def test_changed_probe_and_third_arm_rejected(self):
        self.prepared()
        atomic_json(self.stage / "diagnostic_tools/slot_stability_probe.py", {"changed": True})
        with self.assertRaisesRegex(RuntimeError, "PROBE_CHANGED"):
            self.check_frozen()
        plan = read_json(self.run / "plan.json")
        plan["cases"].append({"id": "third", "eps_averaging": False})
        plan["plan_sha256"] = digest_object({k: v for k, v in plan.items() if k != "plan_sha256"})
        atomic_json(self.run / "plan.json", plan)
        with self.assertRaisesRegex(RuntimeError, "PLAN_CONTRACT"):
            self.check_frozen()

    def test_two_serial_arms_complete_and_charge_wall_and_rank_hours(self):
        self.prepared()
        report, children, popen, terminate = self.execute()
        self.assertEqual(popen.call_count, 2)
        self.assertEqual(len(children), 2)
        self.assertEqual(terminate.call_count, 0)
        self.assertEqual(children[1].start, children[0].start + 1)
        self.assertTrue(report["both_arms_complete_evidence"])
        self.assertEqual(report["interpretation"], "AVERAGING_OFF_SHORT_WINDOW_STABILITY_ONLY")
        self.assertFalse(report["qualification"])
        self.assertFalse(report["formal_R_T_A_conclusion"])
        ledger = runner.closed_ledger(self.stage)
        self.assertEqual([r["attempt"] for r in ledger.data["intervals"][-2:]], [1, 1])
        self.assertAlmostEqual(ledger.total_active_seconds() - runner.BASELINE_ACTIVE_SECONDS, 2)
        self.assertAlmostEqual(report["budget"]["rank_hours"], (runner.BASELINE_ACTIVE_SECONDS + 2) * 4 / 3600)
        self.assertEqual(read_json(self.run / "status.json")["status"], "FINISHED_NOT_QUALIFICATION")
        self.assertTrue((self.run / "report.md").is_file())

    def test_on_completion_without_failure_reproduction_forbids_attribution(self):
        self.prepared()
        report, *_ = self.execute(["SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION"] * 2)
        self.assertEqual(report["interpretation"], "ORIGINAL_FAILURE_NOT_REPRODUCED_NO_ATTRIBUTION")
        self.assertIn("cannot attribute", report["conclusion"])

    def test_nonfinite_is_diagnostic_and_allows_only_off_arm(self):
        self.prepared()
        report, _, popen, _ = self.execute(["NONFINITE_FIELD", "NONFINITE_FIELD"])
        self.assertEqual(popen.call_count, 2)
        self.assertEqual(report["interpretation"], "BOTH_ARMS_UNSTABLE_NO_FIX_ESTABLISHED")

    def test_spawn_failure_charged_and_stops_second_arm(self):
        self.prepared()
        report, _, popen, _ = self.execute(spawn_error=OSError("mock spawn failure"))
        self.assertEqual(popen.call_count, 1)
        self.assertEqual([arm["status"] for arm in report["arms"]], ["CRASHED", "NOT_RUN"])
        self.assertAlmostEqual(report["budget"]["phase_active_solver_seconds"], .1)
        self.assertFalse(report["both_arms_complete_evidence"])
        runner.closed_ledger(self.stage)

    def test_nonzero_exit_never_promoted_to_success(self):
        self.prepared()
        report, _, popen, _ = self.execute(returncode=9)
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(report["arms"][0]["status"], "CRASHED")
        self.assertFalse(report["arms"][0]["evidence_validated"])

    def test_memory_limit_terminates_owned_child_and_accounts_shutdown(self):
        self.prepared()
        report, children, popen, terminate = self.execute(duration=100, rss=4 * 1024**3)
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(terminate.call_count, 1)
        self.assertEqual(children[0].returncode, -15)
        self.assertEqual(report["arms"][0]["status"], "MEMORY_LIMIT")
        self.assertAlmostEqual(report["budget"]["phase_active_solver_seconds"], .45)
        runner.closed_ledger(self.stage)

    def test_cleanup_failure_never_falsely_closes_active_interval(self):
        self.prepared()
        report, _, popen, _ = self.execute(duration=100, rss=4 * 1024**3,
                                          cleanup_error=OSError("mock unkillable child"))
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(report["arms"][0]["status"], "CRASHED")
        self.assertTrue(report["arms"][0]["cleanup_unconfirmed"])
        self.assertEqual(read_json(self.run / "status.json")["status"], "CHILD_STOP_UNCONFIRMED")
        self.assertIsNone(read_json(self.ledger.path)["intervals"][-1]["end"])
        with self.assertRaisesRegex(RuntimeError, "UNCLOSED_PRIOR_INTERVAL"):
            runner.closed_ledger(self.stage)
        self.assertLessEqual(report["arms"][0]["elapsed_s"], 30)

    def test_exited_mpiexec_with_live_rank_cannot_close_or_start_next_arm_early(self):
        self.prepared()
        report, children, popen, signals = self.execute(orphan_after_exit=True)
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(report["arms"][0]["status"], "CRASHED")
        self.assertIn("LIVE_SESSION_MEMBERS", report["arms"][0]["error"])
        self.assertTrue(report["arms"][0]["session_cleanup"]["stopped"])
        self.assertEqual(signals.call_args.args[0]["pid"], children[0].pid + 100)
        self.assertGreater(report["arms"][0]["elapsed_s"], 1)
        runner.closed_ledger(self.stage)

    def test_term_resistant_rank_is_killed_within_single_shared_shutdown_deadline(self):
        self.prepared()
        report, _, popen, signals = self.execute(duration=100, rss=4 * 1024**3,
                                                term_resistant=True)
        self.assertEqual(popen.call_count, 1)
        self.assertEqual([call.args[2] for call in signals.call_args_list], [signal.SIGTERM, signal.SIGKILL])
        self.assertLessEqual(report["arms"][0]["elapsed_s"], 30)
        self.assertGreaterEqual(report["arms"][0]["elapsed_s"], 20)
        runner.closed_ledger(self.stage)

    def test_shared_ledger_and_run_locks_remain_owned_while_children_launch(self):
        self.prepared()

        def assert_locked(_child):
            for path in (self.stage / "budget_ledger.lock", self.run / "run.lock"):
                with self.assertRaisesRegex(RuntimeError, "Already locked"):
                    with FileLock(path):
                        self.fail("lock unexpectedly released")

        self.execute(callback=assert_locked)

    def test_signal_interrupt_stops_subsequent_arm_and_closes_ledger(self):
        self.prepared()

        def interrupt(_child):
            import signal
            signal.getsignal(signal.SIGTERM)(signal.SIGTERM, None)

        report, _, popen, terminate = self.execute(duration=100, callback=interrupt)
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(terminate.call_count, 1)
        self.assertEqual(report["arms"][0]["status"], "INTERRUPTED")
        runner.closed_ledger(self.stage)

    def test_phase_exhaustion_blocks_launch(self):
        self.prepared()
        self.ledger.add_interval(40000, 41790, 4, "other_accounted_work", "CRASHED")
        report, _, popen, _ = self.execute()
        self.assertEqual(popen.call_count, 0)
        self.assertEqual(report["arms"][0]["status"], "NOT_RUN")

    def test_timeout_reserves_shutdown_and_uses_reduced_phase_allowance(self):
        self.prepared()
        self.ledger.add_interval(40000, 41760, 4, "other_accounted_work", "CRASHED")
        report, _, popen, terminate = self.execute(duration=100)
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(terminate.call_count, 1)
        arm = report["arms"][0]
        self.assertEqual(arm["status"], "TIME_LIMIT_UNQUALIFIED")
        self.assertEqual(arm["allowance_seconds_including_shutdown"], 40)
        self.assertAlmostEqual(arm["elapsed_s"], 10.45)
        self.assertAlmostEqual(report["budget"]["phase_active_solver_seconds"], 1770.45)
        self.assertEqual(read_json(self.run / "averaging_on/invocation.json")["solver_timeout_seconds"], 10)

    def test_run_is_server_only_and_never_resumes_existing_attempt(self):
        self.prepared()
        with patch.object(runner.sys, "platform", "darwin"):
            with self.assertRaisesRegex(RuntimeError, "SERVER_ONLY"):
                runner.run(self.stage)
        self.execute()
        with patch.object(runner.sys, "platform", "linux"):
            with self.assertRaisesRegex(RuntimeError, "EXISTING_DIAGNOSTIC_ATTEMPT_NO_RESUME"):
                runner.run(self.stage)

    def test_evidence_hash_incomplete_trace_and_forged_success_rejected(self):
        plan = self.prepared()
        out = self.run / "averaging_on"
        out.mkdir()
        original = synthetic_probe_result(out, plan, "averaging_on", "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION")
        runner.validate_probe_result(out, plan, "averaging_on")
        for update in ({"evidence_sha256": {}}, {"t_reached": 98},
                       {"source_case_sha256": "bad"}, {"status": "QUALIFIED"}):
            with self.subTest(update=update):
                atomic_json(out / "result.json", dict(original, **update))
                with self.assertRaises(RuntimeError):
                    runner.validate_probe_result(out, plan, "averaging_on")
        atomic_json(out / "result.json", original)
        atomic_json(out / "sample_grid.npz", {"tampered": True})
        with self.assertRaisesRegex(RuntimeError, "EVIDENCE_HASH_MISMATCH"):
            runner.validate_probe_result(out, plan, "averaging_on")

    def test_report_never_attributes_incomplete_evidence(self):
        plan = self.prepared()
        arms = [{"case_id": "averaging_on", "status": "EARLY_GROWTH_DETECTED", "evidence_validated": True},
                {"case_id": "averaging_off", "status": "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION",
                 "evidence_validated": False}]
        report = runner.report_for(plan, arms, self.ledger)
        self.assertEqual(report["interpretation"], "INCOMPLETE_EVIDENCE_NO_ATTRIBUTION")
        self.assertFalse(report["both_arms_complete_evidence"])

    def test_caller_plan_run_id_controls_arm_status_accounting_and_report(self):
        self.prepared()
        original_id = runner.RUN_ID
        new_id = "slot_stability_20260922_G_off"
        self.plan = dict(self.plan, run_id=new_id,
                         run_dir=str(self.stage / "repair_diagnostics" / new_id))
        self.run = Path(self.plan["run_dir"])
        self.run.mkdir()
        report, _, popen, _ = self.execute(["SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION"], direct_arm="averaging_off")
        self.assertEqual(popen.call_count, 1)
        self.assertEqual(runner.RUN_ID, original_id)
        self.assertEqual(report["run_id"], new_id)
        self.assertEqual(read_json(self.run / "status.json")["run_id"], new_id)
        self.assertEqual(read_json(self.run / "averaging_off/process.json")["run_id"], new_id)
        self.assertEqual(read_json(self.ledger.path)["intervals"][-1]["task_id"], new_id + ":averaging_off")

    def test_actual_growth_trigger_survives_extra_dt_snapshot_without_rewriting_evidence(self):
        plan = self.prepared()
        out = self.run / "averaging_on"
        out.mkdir()
        result = synthetic_probe_result(out, plan, "averaging_on", "EARLY_GROWTH_DETECTED")
        trace = read_json(out / "trace.json")
        samples = [growth_sample(t, 0.0) for t in range(0, 16, 2)] + actual_growth_samples()
        trace["samples"] = samples
        atomic_json(out / "trace.json", trace)
        result["t_reached"] = 20.0078125
        result["trace_summary"] = {"sample_count": len(samples), "last_sample": samples[-1]}
        result["evidence_sha256"]["trace.json"] = digest_file(out / "trace.json")
        atomic_json(out / "result.json", result)
        before = {path: path.read_bytes() for path in out.iterdir()}
        validated = runner.validate_probe_result(out, plan, "averaging_on")
        self.assertEqual(validated["growth_trigger_evidence"]["trigger_index"], 10)
        self.assertEqual(validated["growth_trigger_evidence"]["witness_times"], [16, 18, 20])
        self.assertEqual(before, {path: path.read_bytes() for path in before})

    def test_nonfinite_point_probe_is_valid_even_when_full_grid_envelope_is_finite(self):
        plan = self.prepared()
        out = self.run / "averaging_on"
        out.mkdir()
        result = synthetic_probe_result(out, plan, "averaging_on", "NONFINITE_FIELD")
        trace = read_json(out / "trace.json")
        last = trace["samples"][-1]
        last["envelope"]["all_fields_finite"] = True
        last["probes"] = {"slot_bottom": {"finite": False}}
        atomic_json(out / "trace.json", trace)
        result["trace_summary"]["last_sample"] = last
        result["evidence_sha256"]["trace.json"] = digest_file(out / "trace.json")
        atomic_json(out / "result.json", result)
        self.assertEqual(runner.validate_probe_result(out, plan, "averaging_on")["status"], "NONFINITE_FIELD")
        last["probes"]["slot_bottom"]["finite"] = True
        atomic_json(out / "trace.json", trace)
        result["evidence_sha256"]["trace.json"] = digest_file(out / "trace.json")
        atomic_json(out / "result.json", result)
        with self.assertRaisesRegex(RuntimeError, "NONFINITE_EVIDENCE_MISSING"):
            runner.validate_probe_result(out, plan, "averaging_on")


def growth_sample(simulation_time, logratio, trigger=None):
    return {"simulation_time": simulation_time, "early_guard_status": trigger,
            "envelope": {"all_fields_finite": True, "log10_intensity_over_reference_peak": logratio},
            "probes": {"source": {"finite": True}}}


def actual_growth_samples():
    return [growth_sample(16, 18.756162), growth_sample(18, 27.911704),
            growth_sample(20, 37.067246, "EARLY_GROWTH_DETECTED"),
            growth_sample(20.0078125, 37.103010, "EARLY_GROWTH_DETECTED")]


class GrowthTriggerEvidenceTests(unittest.TestCase):
    def test_original_t16_t18_t20_witness_ignores_extra_final_dt_sample(self):
        samples = actual_growth_samples()
        original = deepcopy(samples)
        evidence = runner.growth_trigger_evidence(samples)
        self.assertEqual(evidence, {"trigger_index": 2, "trigger_time": 20,
                                   "witness_times": [16, 18, 20],
                                   "witness_logratios": [18.756162, 27.911704, 37.067246]})
        self.assertEqual(samples, original)

    def test_no_recorded_trigger_is_not_replaced_by_any_numerically_passing_window(self):
        samples = actual_growth_samples()
        for sample in samples:
            sample["early_guard_status"] = None
        with self.assertRaisesRegex(RuntimeError, "FIRST_TRIGGER_EVIDENCE_MISSING"):
            runner.growth_trigger_evidence(samples)
        with self.assertRaisesRegex(RuntimeError, "FIRST_TRIGGER_EVIDENCE_MISSING"):
            runner.growth_trigger_evidence([])

    def test_wrong_first_trigger_is_not_replaced_by_later_valid_trigger(self):
        samples = [growth_sample(0, 13), growth_sample(2, 13.5),
                   growth_sample(4, 14, "EARLY_GROWTH_DETECTED"),
                   growth_sample(6, 17), growth_sample(8, 20, "EARLY_GROWTH_DETECTED")]
        with self.assertRaisesRegex(RuntimeError, "FIRST_TRIGGER_THRESHOLDS_NOT_MET"):
            runner.growth_trigger_evidence(samples)
        samples = actual_growth_samples()
        samples[0]["early_guard_status"] = "EARLY_GROWTH_DETECTED"
        with self.assertRaisesRegex(RuntimeError, "FIRST_TRIGGER_EVIDENCE_MISSING"):
            runner.growth_trigger_evidence(samples)

    def test_irregular_or_shifted_sample_intervals_and_nonfinite_times_are_rejected(self):
        for times in ([16, 18.5, 20], [16.1, 18.1, 20.1], [16, 18, 20.0078125],
                      [16, 16, 20], [16, float("nan"), 20], [16, float("inf"), 20]):
            samples = actual_growth_samples()[:3]
            for sample, value in zip(samples, times):
                sample["simulation_time"] = value
            with self.subTest(times=times), self.assertRaisesRegex(RuntimeError, "SAMPLE_INTERVAL_INVALID"):
                runner.growth_trigger_evidence(samples)

    def test_thresholds_remain_strict_and_all_witness_fields_finite(self):
        for ratios in ([12, 15, 18], [13, 14, 16], [13, float("inf"), 18],
                       [13, float("nan"), 18], [13, None, 18], [13, True, 18]):
            samples = actual_growth_samples()[:3]
            for sample, value in zip(samples, ratios):
                sample["envelope"]["log10_intensity_over_reference_peak"] = value
            with self.subTest(ratios=ratios), self.assertRaisesRegex(RuntimeError, "THRESHOLDS_NOT_MET"):
                runner.growth_trigger_evidence(samples)
        samples = actual_growth_samples()
        samples[1]["probes"]["source"]["finite"] = False
        with self.assertRaisesRegex(RuntimeError, "THRESHOLDS_NOT_MET"):
            runner.growth_trigger_evidence(samples)
        with self.assertRaisesRegex(RuntimeError, "DIAGNOSTIC_CONTRACT_CHANGED"):
            runner.growth_trigger_evidence(actual_growth_samples(), dict(runner.DIAGNOSTIC, growth_ratio_limit=1e10))


class SessionAuditTests(unittest.TestCase):
    def test_orphan_ranks_in_separate_groups_count_but_unrelated_and_zombie_do_not(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for pid, session_id, group, state in ((31, 31, 31, "S"), (32, 31, 32, "R"),
                                                   (33, 33, 33, "S"), (34, 31, 34, "Z")):
                directory = root / str(pid)
                directory.mkdir()
                fields = [state, "1", str(group), str(session_id)] + ["0"] * 15 + [str(pid * 100)]
                (directory / "stat").write_text(f"{pid} (mpi rank ) worker) " + " ".join(fields))
            members = runner.session_live_members(31, root)
            self.assertEqual({entry["pid"] for entry in members}, {31, 32})
            self.assertEqual({entry["pgid"] for entry in members}, {31, 32})
            self.assertEqual({entry["start_ticks"] for entry in members}, {3100, 3200})

    def test_signal_rechecks_session_start_time_and_live_state_before_each_pid(self):
        member = {"pid": 32, "session": 31, "pgid": 32, "state": "R", "start_ticks": 3200}
        for replacement in (None, dict(member, session=90), dict(member, start_ticks=9900),
                            dict(member, state="Z")):
            with self.subTest(replacement=replacement), \
                    patch.object(runner, "_process_identity", return_value=replacement), \
                    patch.object(runner.os, "kill") as kill:
                runner.signal_session_member(member, 31, signal.SIGTERM)
                kill.assert_not_called()
        with patch.object(runner, "_process_identity", return_value=member), \
                patch.object(runner.os, "kill") as kill:
            runner.signal_session_member(member, 31, signal.SIGTERM)
            kill.assert_called_once_with(32, signal.SIGTERM)


if __name__ == "__main__":
    unittest.main()
