"""Run only the previously unexecuted averaging-off arm in a new directory.

The completed F averaging-on run and its CRASHED supervisor/ledger status stay
immutable. Its growth evidence is re-audited offline into G, without a new on
run, a field-state restart, a reset budget, or any production qualification.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import signal
import sys

import run_slot_stability_diagnostic as runner
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.queue import FileLock, read_json

RUN_ID = "slot_stability_20260922_G_off"
SOURCE_RUN_ID = "slot_stability_20260922_F"
SERVER_STAGE = Path("/root/autodl-tmp/meep_sim/meep_screen")
APPROVED_F_PLAN_SHA256 = "314bc57323b570e16bb16885460ab4489d29f316ebfab987aff30faa3ff2798b"
APPROVED_F_ON_RESULT_SHA256 = "47f7de031180048fbfd05d8c9e0f7a1439ff424ca7d0b788fa30336864939502"
APPROVED_LEDGER_SHA256 = "581d12f6fc7e4df8ab7ef60329686d564dc0c9bd2169270fd59076957f48affc"
APPROVED_CURRENT_ACTIVE_SECONDS = 31194.026368379593
EXECUTE_CASE_IDS = ["averaging_off"]


def _run_path(stage, run_id):
    path = stage / "repair_diagnostics" / run_id
    if path.resolve() != path:
        raise RuntimeError("DIAGNOSTIC_RUN_SYMLINK_PATH_PROHIBITED")
    return path


def _file_hashes(directory):
    """Freeze every existing F file, including original supervisor evidence."""
    result = {}
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise RuntimeError("SOURCE_DIAGNOSTIC_SYMLINK_PROHIBITED")
        if path.is_file():
            result[str(path.relative_to(directory))] = digest_file(path)
    return result


def _prefix_unchanged(ledger, snapshot):
    prefix = snapshot["intervals"]
    if (ledger.data["intervals"][:len(prefix)] != prefix or
            {key: value for key, value in ledger.data.items() if key != "intervals"} !=
            {key: value for key, value in snapshot.items() if key != "intervals"}):
        raise RuntimeError("HISTORICAL_LEDGER_CHANGED")


def audit_original_on(stage, ledger):
    """Read-only F evidence verification; never rewrite its result or ledger."""
    stage = Path(stage).resolve()
    source = _run_path(stage, SOURCE_RUN_ID)
    plan_path, result_path = source / "plan.json", source / "averaging_on/result.json"
    for path, expected in ((plan_path, APPROVED_F_PLAN_SHA256),
                           (result_path, APPROVED_F_ON_RESULT_SHA256)):
        runner._regular_file(path, "F_FROZEN_EVIDENCE")
        if digest_file(path) != expected:
            raise RuntimeError("APPROVED_F_EVIDENCE_CHANGED")
    plan = read_json(plan_path)
    if (plan.get("run_id") != SOURCE_RUN_ID or plan.get("run_dir") != str(source)
            or plan.get("stage_root") != str(stage) or plan.get("qualification") is not False
            or plan.get("diagnostic") != runner.DIAGNOSTIC
            or plan.get("budgets") != runner._budget_contract()
            or plan.get("cases") != [{"id": key, "eps_averaging": value} for key, value in runner.ARMS.items()]
            or plan.get("case_order") != list(runner.ARMS)
            or plan.get("automatic_retry_authorized") is not False
            or plan.get("implementation_sha256") != runner.APPROVED_IMPLEMENTATION_SHA256
            or implementation_digest() != runner.APPROVED_IMPLEMENTATION_SHA256):
        raise RuntimeError("F_FROZEN_PLAN_CONTRACT_CHANGED")
    if plan.get("plan_sha256") != digest_object({key: value for key, value in plan.items() if key != "plan_sha256"}):
        raise RuntimeError("F_PLAN_HASH_CHANGED")
    authorization = read_json(runner._regular_file(source / "authorization.json", "F_AUTHORIZATION"))
    if digest_object(authorization) != plan["authorization_sha256"]:
        raise RuntimeError("F_AUTHORIZATION_CHANGED")
    snapshot = authorization["ledger_snapshot"]
    if digest_object(snapshot) != authorization["ledger_snapshot_sha256"]:
        raise RuntimeError("F_LEDGER_SNAPSHOT_CHANGED")
    _prefix_unchanged(ledger, snapshot)
    paths, case = runner._source_inputs(stage)
    paths["ledger_snapshot"] = stage / "budget_ledger.json"
    if set(authorization["inputs"]) != set(paths):
        raise RuntimeError("F_FROZEN_INPUT_SET_CHANGED")
    for name, original in paths.items():
        entry = authorization["inputs"][name]
        frozen = source / "inputs" / (name + ".json")
        runner._regular_file(frozen, "F_FROZEN_INPUT")
        if (entry["source_path"] != str(original) or entry["snapshot_path"] != str(frozen)
                or digest_file(frozen) != entry["sha256"]):
            raise RuntimeError("F_FROZEN_INPUT_CHANGED")
        if name == "ledger_snapshot":
            if read_json(frozen) != snapshot or entry["sha256"] != runner.APPROVED_LEDGER_SHA256:
                raise RuntimeError("F_ORIGINAL_LEDGER_CHANGED")
        elif digest_file(original) != entry["sha256"]:
            raise RuntimeError("F_PHYSICAL_INPUT_CHANGED")
    for name in ("source_case", "source_manifest", "source_result", "material"):
        if plan.get(name + "_path") != str(paths[name]) or plan.get(name + "_sha256") != digest_file(paths[name]):
            raise RuntimeError("F_PLAN_SOURCE_CHANGED")
    if (plan["material_sha256"] != case["material_sha256"] or
            digest_file(stage / "diagnostic_tools/slot_stability_probe.py") != plan["probe_sha256"]):
        raise RuntimeError("F_MATERIAL_OR_PROBE_CHANGED")
    if (source / "averaging_off").exists() or (source / "averaging_off").is_symlink():
        raise RuntimeError("AVERAGING_OFF_ALREADY_ATTEMPTED")
    old_records = [record for record in ledger.data["intervals"]
                   if str(record.get("task_id", "")).startswith(SOURCE_RUN_ID + ":")]
    if (len(old_records) != 1 or old_records[0].get("task_id") != SOURCE_RUN_ID + ":averaging_on"
            or old_records[0].get("status") != "CRASHED" or old_records[0].get("attempt") != 1
            or old_records[0].get("ranks") != 4 or old_records[0].get("end") is None):
        raise RuntimeError("F_ON_ONLY_CLOSED_CRASHED_ATTEMPT_REQUIRED")
    process = read_json(runner._regular_file(source / "averaging_on/process.json", "F_PROCESS"))
    if (process.get("case_id") != "averaging_on" or process.get("status") != "CRASHED"
            or process.get("exit_code") != 0 or process.get("qualification") is not False
            or "GROWTH_EVIDENCE_MISSING" not in str(process.get("error", ""))
            or process.get("cleanup_unconfirmed") is True):
        raise RuntimeError("EXPECTED_F_SUPERVISOR_MISCLASSIFICATION_REQUIRED")
    pid = process.get("pid")
    if isinstance(pid, bool) or not isinstance(pid, int) or pid <= 0 or runner.session_live_members(pid):
        raise RuntimeError("F_ORIGINAL_SESSION_NOT_CONFIRMED_STOPPED")
    result = runner.validate_probe_result(source / "averaging_on", plan, "averaging_on")
    if result["status"] != "EARLY_GROWTH_DETECTED":
        raise RuntimeError("F_GROWTH_REPRODUCTION_REQUIRED")
    trigger = runner.growth_trigger_evidence(read_json(source / "averaging_on/trace.json")["samples"], plan["diagnostic"])
    arm = {"case_id": "averaging_on", "status": result["status"], "qualification": False,
           "evidence_validated": True, "probe_result": result, "result_sha256": digest_file(result_path),
           "elapsed_s": process.get("elapsed_s", old_records[0]["end"] - old_records[0]["start"]),
           "simulation_time": result["t_reached"], "reused_existing_evidence": True,
           "source_run_id": SOURCE_RUN_ID, "source_directory": str(source / "averaging_on"),
           "original_supervisor_status": "CRASHED", "new_solver_run": False}
    audit = {"schema_version": 1, "qualification": False, "source_run_id": SOURCE_RUN_ID,
             "source_plan_file_sha256": digest_file(plan_path),
             "source_result_file_sha256": digest_file(result_path),
             "original_supervisor_status": "CRASHED", "original_supervisor_error": process["error"],
             "original_ledger_status_preserved": "CRASHED", "offline_reaudit_status": "EARLY_GROWTH_DETECTED",
             "growth_trigger_evidence": trigger, "reused_on_arm": arm,
             "no_new_on_run": True, "field_state_resume": False,
             "explanation": "The first latched growth trigger has three qualifying fixed-interval samples. A later one-timestep terminal snapshot changed the final-three window and caused the original supervisor rejection. Original F artifacts and accounting are preserved."}
    return plan, audit


def prepare(stage):
    stage = Path(stage).resolve()
    run = _run_path(stage, RUN_ID)
    with FileLock(stage / "budget_ledger.lock"):
        if run.exists() or run.is_symlink():
            raise FileExistsError("OFF_CONTINUATION_ALREADY_EXISTS_NO_OVERWRITE")
        ledger = runner.closed_ledger(stage)
        if (digest_file(ledger.path) != APPROVED_LEDGER_SHA256 or
                abs(ledger.total_active_seconds() - APPROVED_CURRENT_ACTIVE_SECONDS) > 1e-6):
            raise RuntimeError("APPROVED_POST_F_LEDGER_CHANGED")
        old_plan, audit = audit_original_on(stage, ledger)
        source_hashes = _file_hashes(_run_path(stage, SOURCE_RUN_ID))
        authorization = {"schema_version": 1, "run_id": RUN_ID, "source_run_id": SOURCE_RUN_ID,
                         "execute_case_ids": EXECUTE_CASE_IDS, "no_new_on_run": True,
                         "qualification": False, "automatic_retry_authorized": False,
                         "ledger_path": str(ledger.path), "ledger_file_sha256": APPROVED_LEDGER_SHA256,
                         "ledger_snapshot": deepcopy(ledger.data), "ledger_snapshot_sha256": digest_object(ledger.data),
                         "source_directory_files_sha256": source_hashes,
                         "offline_reaudit_sha256": digest_object(audit), "budgets": runner._budget_contract()}
        plan = deepcopy(old_plan)
        plan.update(run_id=RUN_ID, run_dir=str(run), created_utc=utcnow(),
                    runner_sha256=digest_file(stage / Path(runner.__file__).name),
                    continuation_sha256=digest_file(stage / Path(__file__).name),
                    execute_case_ids=list(EXECUTE_CASE_IDS), no_new_on_run=True,
                    source_run_id=SOURCE_RUN_ID, source_plan_file_sha256=APPROVED_F_PLAN_SHA256,
                    source_on_result_file_sha256=APPROVED_F_ON_RESULT_SHA256,
                    offline_reaudit_sha256=digest_object(audit),
                    authorization_sha256=digest_object(authorization), qualification=False,
                    continuation_scope="Only averaging_off is executed from fresh fields; averaging_on is reused from immutable F evidence. Original 1800-second phase budget is retained.")
        plan["plan_sha256"] = digest_object({key: value for key, value in plan.items() if key != "plan_sha256"})
        run.mkdir(parents=True, exist_ok=False)
        atomic_json(run / "authorization.json", authorization)
        atomic_json(run / "offline_reaudit.json", audit)
        atomic_json(run / "plan.json", plan)
        atomic_json(run / "status.json", {"run_id": RUN_ID, "status": "PREPARED_OFF_ONLY_NOT_STARTED",
                                          "qualification": False, "no_new_on_run": True,
                                          "budget": runner.budget_usage(ledger)})
        return plan


def verify_prepared(stage, plan, authorization, audit, ledger):
    run = _run_path(stage, RUN_ID)
    expected = {"schema_version": 1, "run_id": RUN_ID, "run_dir": str(run), "stage_root": str(stage),
                "execute_case_ids": EXECUTE_CASE_IDS, "no_new_on_run": True,
                "source_run_id": SOURCE_RUN_ID, "qualification": False, "automatic_retry_authorized": False,
                "source_plan_file_sha256": APPROVED_F_PLAN_SHA256,
                "source_on_result_file_sha256": APPROVED_F_ON_RESULT_SHA256,
                "budgets": runner._budget_contract(), "diagnostic": runner.DIAGNOSTIC,
                "case_order": list(runner.ARMS),
                "cases": [{"id": key, "eps_averaging": value} for key, value in runner.ARMS.items()]}
    if any(plan.get(key) != value for key, value in expected.items()):
        raise RuntimeError("OFF_ONLY_PLAN_CONTRACT_CHANGED")
    if (plan.get("plan_sha256") != digest_object({key: value for key, value in plan.items() if key != "plan_sha256"})
            or plan.get("authorization_sha256") != digest_object(authorization)):
        raise RuntimeError("OFF_ONLY_PLAN_OR_AUTHORIZATION_CHANGED")
    if (authorization.get("run_id") != RUN_ID or authorization.get("execute_case_ids") != EXECUTE_CASE_IDS
            or authorization.get("budgets") != runner._budget_contract()
            or authorization.get("ledger_path") != str(ledger.path)
            or authorization.get("ledger_file_sha256") != APPROVED_LEDGER_SHA256
            or authorization.get("ledger_snapshot_sha256") != digest_object(authorization["ledger_snapshot"])):
        raise RuntimeError("OFF_ONLY_AUTHORIZATION_CONTRACT_CHANGED")
    _prefix_unchanged(ledger, authorization["ledger_snapshot"])
    if any(str(record.get("task_id", "")).startswith(RUN_ID + ":") for record in ledger.data["intervals"]):
        raise RuntimeError("OFF_CONTINUATION_ALREADY_ATTEMPTED_NO_RESUME")
    old_plan, fresh_audit = audit_original_on(stage, ledger)
    if (_file_hashes(_run_path(stage, SOURCE_RUN_ID)) != authorization["source_directory_files_sha256"]
            or audit != fresh_audit or digest_object(audit) != authorization["offline_reaudit_sha256"]
            or digest_object(audit) != plan["offline_reaudit_sha256"]):
        raise RuntimeError("F_EVIDENCE_OR_OFFLINE_REAUDIT_CHANGED")
    for key in ("source_case_path", "source_case_sha256", "source_manifest_path", "source_manifest_sha256",
                "source_result_path", "source_result_sha256", "material_path", "material_sha256",
                "implementation_sha256", "probe_sha256"):
        if plan.get(key) != old_plan.get(key):
            raise RuntimeError("OFF_CONTINUATION_PHYSICAL_INPUT_CHANGED")
    if (digest_file(stage / Path(runner.__file__).name) != plan["runner_sha256"]
            or digest_file(stage / Path(__file__).name) != plan["continuation_sha256"]):
        raise RuntimeError("OFF_CONTINUATION_IMPLEMENTATION_CHANGED")
    runner.budget_usage(ledger)


def run(stage, plan_path=None):
    stage = Path(stage).resolve()
    if sys.platform != "linux" or stage != SERVER_STAGE:
        raise RuntimeError("SERVER_ONLY: no local FDTD execution")
    run_dir = _run_path(stage, RUN_ID)
    if plan_path is not None and Path(plan_path).resolve() != run_dir / "plan.json":
        raise RuntimeError("FIXED_OFF_PLAN_PATH_REQUIRED")
    runner._regular_file(run_dir / "plan.json", "OFF_PREPARED_PLAN")
    with FileLock(stage / "budget_ledger.lock"), FileLock(run_dir / "run.lock"):
        ledger = runner.closed_ledger(stage)
        plan, authorization, audit = [read_json(runner._regular_file(run_dir / name, "OFF_FROZEN_INPUT"))
                                      for name in ("plan.json", "authorization.json", "offline_reaudit.json")]
        verify_prepared(stage, plan, authorization, audit, ledger)
        if (any((run_dir / name).exists() for name in ("started.json", "report.json", "averaging_on", "averaging_off"))):
            raise RuntimeError("OFF_CONTINUATION_ALREADY_STARTED_NO_RESUME")
        cancelled = [False]
        arms = [deepcopy(audit["reused_on_arm"])]
        def cancel(*unused):
            cancelled[0] = True
        previous = {sig: signal.signal(sig, cancel) for sig in (signal.SIGTERM, signal.SIGINT)}
        atomic_json(run_dir / "started.json", {"run_id": RUN_ID, "supervisor_pid": os.getpid(),
                    "started_utc": utcnow(), "plan_sha256": plan["plan_sha256"],
                    "execute_case_ids": EXECUTE_CASE_IDS, "no_new_on_run": True, "automatic_resume": False})
        try:
            # Deliberately one literal arm call: no loop over the two-arm schema.
            arms.append(runner.run_arm(stage, plan, ledger, "averaging_off", cancelled))
        except Exception as exc:
            arms.append({"case_id": "averaging_off", "status": "CRASHED", "qualification": False,
                         "error": f"{type(exc).__name__}: {exc}"})
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
            report = runner.report_for(plan, arms, ledger)
            report.update(no_new_on_run=True, source_run_id=SOURCE_RUN_ID,
                          execute_case_ids=EXECUTE_CASE_IDS, field_state_resume=False,
                          source_on_result_file_sha256=APPROVED_F_ON_RESULT_SHA256,
                          offline_reaudit_sha256=digest_object(audit),
                          original_supervisor_and_ledger_status_preserved="CRASHED")
            runner.write_report(run_dir, report)
            atomic_json(run_dir / "status.json", {"run_id": RUN_ID,
                        "status": "CHILD_STOP_UNCONFIRMED" if any(arm.get("cleanup_unconfirmed") for arm in arms)
                                  else "FINISHED_NOT_QUALIFICATION",
                        "qualification": False, "no_new_on_run": True,
                        "interpretation": report["interpretation"], "budget": report["budget"],
                        "arms": {arm["case_id"]: arm["status"] for arm in report["arms"]}})
        return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--stage-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--plan", type=Path)
    args = parser.parse_args(argv)
    if args.prepare and args.plan is not None:
        parser.error("--prepare uses the fixed plan path")
    result = prepare(args.stage_root) if args.prepare else run(args.stage_root, args.plan)
    print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
    return 0 if args.prepare or result["both_arms_complete_evidence"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
