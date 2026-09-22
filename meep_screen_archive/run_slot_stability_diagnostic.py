"""Freeze and supervise exactly two short slot-stability arms, never qualify.

``--prepare`` only validates and freezes inputs. ``--run`` is foreground: its
caller may detach it; this module never launches or resumes its own supervisor.
All solver imports and MPI execution remain in the isolated probe.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from ti2d.common import digest_file, digest_object, implementation_digest, utcnow
from ti2d.queue import (FileLock, MPIEXEC, PYTHON, SolverLedger, atomic_json,
                        cgroup_resources, process_tree_rss, read_json,
                        verify_manifest)

RUN_ID = "slot_stability_20260922_F"
SOURCE_RUN_ID = "baseline_r32_20260922_E"
SOURCE_TASK_ID = "baseline_r32_p_plus30"
SERVER_STAGE = Path("/root/autodl-tmp/meep_sim/meep_screen")
BASELINE_ACTIVE_SECONDS = 30952.50267100334
PHASE_CAP_SECONDS = 1800.0
STAGE_CAP_SECONDS = 172800.0
PER_CASE_SECONDS = 1200.0
SHUTDOWN_RESERVE_SECONDS = 30.0
APPROVED_SOURCE_CASE_SHA256 = "3e6b9631a47a3c19040e74fd3c39cfef8df228f45e7bf766782c3e1fb414a85c"
APPROVED_IMPLEMENTATION_SHA256 = "24d09c50dcf6c0cb9271947dbcb04f03b62d9d0500bd7bace34e6909d54834c1"
APPROVED_LEDGER_SHA256 = "74170c24b8e3bfc809631843ae49412dcaba094aa2a5ea5d2d0bd6569185089f"
APPROVED_FAILURE_SHA256 = "73aac73982bb49409d6c322b1cb9cfde1d9206fb8187dce09c3fc55bd984b2f3"
ARMS = {"averaging_on": True, "averaging_off": False}
DIAGNOSTIC = {
    "stop_time": 100.0, "sample_interval": 2.0,
    "reference_intensity_peak": 0.3334312033062042,
    "growth_ratio_limit": 1e12, "growth_consecutive_samples": 3,
    "growth_step_factor": 10.0, "patch_radius_cells": 8,
}
CONTINUE_STATUSES = {"EARLY_GROWTH_DETECTED", "NONFINITE_FIELD",
                     "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION"}
STOP_STATUSES = {"CRASHED", "MEMORY_LIMIT", "INTERRUPTED", "TIME_LIMIT_UNQUALIFIED"}
THREAD_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")


def _regular_file(path, label):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"{label}_REGULAR_FILE_REQUIRED")
    return path


def closed_ledger(stage):
    path = _regular_file(Path(stage) / "budget_ledger.json", "EXISTING_SHARED_LEDGER")
    ledger = SolverLedger(path)
    for record in ledger.data["intervals"]:
        if record.get("end") is None:
            raise RuntimeError("UNCLOSED_PRIOR_INTERVAL")
        if not all(isinstance(record.get(key), (int, float)) and
                   math.isfinite(record[key]) for key in ("start", "end")):
            raise RuntimeError("NONFINITE_LEDGER_INTERVAL")
        if record["end"] < record["start"]:
            raise RuntimeError("REVERSED_LEDGER_INTERVAL")
    return ledger


def budget_usage(ledger, now=None):
    total = ledger.total_active_seconds(now)
    if total < BASELINE_ACTIVE_SECONDS - 1e-6:
        raise RuntimeError("SHARED_LEDGER_BELOW_FROZEN_BASELINE")
    used = max(0.0, total - BASELINE_ACTIVE_SECONDS)
    return {
        "stage_active_solver_seconds": total, "stage_cap_seconds": STAGE_CAP_SECONDS,
        "stage_remaining_solver_seconds": max(0.0, STAGE_CAP_SECONDS - total),
        "phase_active_solver_seconds": used, "phase_cap_seconds": PHASE_CAP_SECONDS,
        "phase_remaining_solver_seconds": max(0.0, PHASE_CAP_SECONDS - used),
        "effective_remaining_solver_seconds": max(0.0, min(STAGE_CAP_SECONDS - total,
                                                            PHASE_CAP_SECONDS - used)),
        "rank_hours": ledger.rank_hours(now),
    }


def _source_inputs(stage):
    source = stage / "runs" / SOURCE_RUN_ID
    paths = {
        "source_manifest": source / "manifest.json",
        "source_case": source / "cases" / (SOURCE_TASK_ID + ".json"),
        "source_result": source / "tasks" / SOURCE_TASK_ID / "attempt_1" / "result.json",
    }
    for name, path in paths.items():
        _regular_file(path, name.upper())
    if digest_file(paths["source_case"]) != APPROVED_SOURCE_CASE_SHA256:
        raise RuntimeError("APPROVED_SOURCE_CASE_CHANGED")
    if digest_file(paths["source_result"]) != APPROVED_FAILURE_SHA256:
        raise RuntimeError("APPROVED_SOURCE_FAILURE_CHANGED")
    core = implementation_digest()
    if core != APPROVED_IMPLEMENTATION_SHA256:
        raise RuntimeError("APPROVED_IMPLEMENTATION_CHANGED")
    manifest, case, failure = (read_json(paths[key]) for key in
                               ("source_manifest", "source_case", "source_result"))
    verify_manifest(manifest, source, core)
    if (manifest.get("run_id") != SOURCE_RUN_ID or len(manifest["tasks"]) != 1 or
            manifest["tasks"][0]["id"] != SOURCE_TASK_ID):
        raise RuntimeError("EXACT_SOURCE_TASK_REQUIRED")
    manifest_case = Path(manifest["tasks"][0]["case_path"])
    if not manifest_case.is_absolute():
        manifest_case = source / manifest_case
    if manifest_case.resolve() != paths["source_case"].resolve():
        raise RuntimeError("SOURCE_MANIFEST_CASE_PATH_CHANGED")
    expected = {"id": SOURCE_TASK_ID, "case": "structure", "geometry_id": "baseline",
                "resolution": 32, "courant": .25, "polarization": "p",
                "emission_theta_deg": 30.0, "wavelength_um": 10.5, "mpi_ranks": 4,
                "mpi_partition_policy": "fixed_binary_Y0_X0_ranks0123"}
    if any(case.get(key) != value for key, value in expected.items()):
        raise RuntimeError("UNAPPROVED_SOURCE_CONFIGURATION")
    if (failure.get("status") != "CRASHED" or
            failure.get("error") != "envelope_power must be finite"):
        raise RuntimeError("EXPECTED_SOURCE_NUMERICAL_FAILURE_REQUIRED")
    material = Path(case["material_path"])
    paths["material"] = material if material.is_absolute() else stage / material
    _regular_file(paths["material"], "MATERIAL")
    if (digest_file(paths["material"]) != case["material_sha256"] or
            manifest.get("material_sha256") != case["material_sha256"]):
        raise RuntimeError("FROZEN_MATERIAL_CHANGED")
    return paths, case


def _budget_contract():
    return {
        "baseline_active_seconds": BASELINE_ACTIVE_SECONDS,
        "phase_seconds": PHASE_CAP_SECONDS, "stage_total_seconds": STAGE_CAP_SECONDS,
        "single_seconds": PER_CASE_SECONDS, "mpi_ranks": 4, "memory_fraction": .7,
    }


def prepare(stage):
    """Freeze all inputs under the shared lock; never import or launch a solver."""
    stage = Path(stage).resolve()
    run = stage / "repair_diagnostics" / RUN_ID
    if run.resolve() != run:
        raise RuntimeError("DIAGNOSTIC_RUN_SYMLINK_PATH_PROHIBITED")
    with FileLock(stage / "budget_ledger.lock"):
        if run.exists() or run.is_symlink():
            raise FileExistsError("DIAGNOSTIC_RUN_ALREADY_EXISTS: no overwrite or automatic resume")
        ledger = closed_ledger(stage)
        if abs(ledger.total_active_seconds() - BASELINE_ACTIVE_SECONDS) > 1e-6:
            raise RuntimeError("APPROVED_BASELINE_MISMATCH")
        if digest_file(ledger.path) != APPROVED_LEDGER_SHA256:
            raise RuntimeError("APPROVED_LEDGER_CHANGED")
        paths, case = _source_inputs(stage)
        paths["ledger_snapshot"] = ledger.path
        probe = _regular_file(stage / "diagnostic_tools" / "slot_stability_probe.py", "PROBE")
        runner = _regular_file(stage / Path(__file__).name, "RUNNER")
        input_hashes = {name: digest_file(path) for name, path in paths.items()}
        frozen_inputs = {name: {"source_path": str(path), "sha256": input_hashes[name],
                                "snapshot_path": str(run / "inputs" / (name + ".json"))}
                         for name, path in paths.items()}
        authorization = {
            "schema_version": 1, "run_id": RUN_ID, "created_utc": utcnow(),
            "scope": "Two serial short-window stability diagnostics only: eps_averaging on then off; no qualification, retry, third arm, scan, R/T/A inference, or production cache writes",
            "source_case_sha256": APPROVED_SOURCE_CASE_SHA256,
            "source_result_sha256": APPROVED_FAILURE_SHA256,
            "implementation_sha256": APPROVED_IMPLEMENTATION_SHA256,
            "ledger_file_sha256": APPROVED_LEDGER_SHA256,
            "ledger_path": str(ledger.path), "ledger_snapshot": deepcopy(ledger.data),
            "ledger_snapshot_sha256": digest_object(ledger.data),
            "inputs": frozen_inputs, "budgets": _budget_contract(),
            "max_parallel": 1, "threads_per_rank": 1,
            "shutdown_reserve_seconds": SHUTDOWN_RESERVE_SECONDS,
            "effective_total_active_seconds": min(STAGE_CAP_SECONDS,
                                                   BASELINE_ACTIVE_SECONDS + PHASE_CAP_SECONDS),
            "automatic_retry_authorized": False, "qualification": False,
        }
        plan = {
            "schema_version": 1, "run_id": RUN_ID, "stage_root": str(stage), "run_dir": str(run),
            "created_utc": authorization["created_utc"], "qualification": False,
            "source_case_path": str(paths["source_case"]),
            "source_case_sha256": APPROVED_SOURCE_CASE_SHA256,
            "source_manifest_path": str(paths["source_manifest"]),
            "source_manifest_sha256": input_hashes["source_manifest"],
            "source_result_path": str(paths["source_result"]),
            "source_result_sha256": APPROVED_FAILURE_SHA256,
            "material_path": str(paths["material"]), "material_sha256": case["material_sha256"],
            "implementation_sha256": APPROVED_IMPLEMENTATION_SHA256,
            "probe_sha256": digest_file(probe), "runner_sha256": digest_file(runner),
            "cases": [{"id": key, "eps_averaging": value} for key, value in ARMS.items()],
            "case_order": list(ARMS), "diagnostic": deepcopy(DIAGNOSTIC),
            "budgets": _budget_contract(), "authorization_sha256": digest_object(authorization),
            "automatic_retry_authorized": False,
        }
        plan["plan_sha256"] = digest_object(plan)
        run.mkdir(parents=True, exist_ok=False)
        (run / "inputs").mkdir()
        for name, path in paths.items():
            shutil.copyfile(path, frozen_inputs[name]["snapshot_path"])
        atomic_json(run / "authorization.json", authorization)
        atomic_json(run / "plan.json", plan)
        atomic_json(run / "status.json", {"status": "PREPARED_NOT_STARTED", "run_id": RUN_ID,
                                          "qualification": False, "budget": budget_usage(ledger)})
        return plan


def verify_frozen(stage, plan, authorization, ledger):
    """Must be called while owning budget_ledger.lock, before any new attempt."""
    stage = Path(stage).resolve()
    run = stage / "repair_diagnostics" / RUN_ID
    expected = {"schema_version": 1, "run_id": RUN_ID, "stage_root": str(stage),
                "run_dir": str(run), "qualification": False,
                "source_case_sha256": APPROVED_SOURCE_CASE_SHA256,
                "source_result_sha256": APPROVED_FAILURE_SHA256,
                "implementation_sha256": APPROVED_IMPLEMENTATION_SHA256,
                "case_order": list(ARMS),
                "cases": [{"id": key, "eps_averaging": value} for key, value in ARMS.items()],
                "diagnostic": DIAGNOSTIC, "budgets": _budget_contract(),
                "automatic_retry_authorized": False}
    if any(plan.get(key) != value for key, value in expected.items()):
        raise RuntimeError("FROZEN_PLAN_CONTRACT_CHANGED")
    if plan.get("plan_sha256") != digest_object({k: v for k, v in plan.items() if k != "plan_sha256"}):
        raise RuntimeError("FROZEN_PLAN_HASH_CHANGED")
    if digest_object(authorization) != plan.get("authorization_sha256"):
        raise RuntimeError("FROZEN_AUTHORIZATION_CHANGED")
    expected_auth = {"run_id": RUN_ID, "source_case_sha256": APPROVED_SOURCE_CASE_SHA256,
                     "source_result_sha256": APPROVED_FAILURE_SHA256,
                     "implementation_sha256": APPROVED_IMPLEMENTATION_SHA256,
                     "ledger_file_sha256": APPROVED_LEDGER_SHA256,
                     "ledger_path": str(stage / "budget_ledger.json"),
                     "budgets": _budget_contract(), "qualification": False,
                     "automatic_retry_authorized": False}
    if any(authorization.get(key) != value for key, value in expected_auth.items()):
        raise RuntimeError("FROZEN_AUTHORIZATION_CONTRACT_CHANGED")
    snapshot = authorization["ledger_snapshot"]
    if digest_object(snapshot) != authorization["ledger_snapshot_sha256"]:
        raise RuntimeError("FROZEN_LEDGER_SNAPSHOT_CHANGED")
    prefix = snapshot["intervals"]
    if (ledger.data["intervals"][:len(prefix)] != prefix or
            {k: v for k, v in ledger.data.items() if k != "intervals"} !=
            {k: v for k, v in snapshot.items() if k != "intervals"}):
        raise RuntimeError("HISTORICAL_LEDGER_CHANGED")
    if any(str(record.get("task_id", "")).startswith(RUN_ID + ":")
           for record in ledger.data["intervals"]):
        raise RuntimeError("EXISTING_DIAGNOSTIC_ATTEMPT_NO_RESUME")
    paths, case = _source_inputs(stage)
    paths["ledger_snapshot"] = stage / "budget_ledger.json"
    if set(authorization["inputs"]) != set(paths):
        raise RuntimeError("FROZEN_INPUT_SET_CHANGED")
    for name, source in paths.items():
        entry = authorization["inputs"][name]
        frozen = run / "inputs" / (name + ".json")
        if entry["snapshot_path"] != str(frozen) or entry["source_path"] != str(source):
            raise RuntimeError("FROZEN_INPUT_PATH_CHANGED")
        _regular_file(frozen, "FROZEN_INPUT")
        if digest_file(frozen) != entry["sha256"]:
            raise RuntimeError("FROZEN_INPUT_CHANGED")
        if name != "ledger_snapshot" and digest_file(source) != entry["sha256"]:
            raise RuntimeError("SOURCE_INPUT_CHANGED_AFTER_PREPARE")
    if (read_json(run / "inputs" / "ledger_snapshot.json") != snapshot or
            authorization["inputs"]["ledger_snapshot"]["sha256"] != APPROVED_LEDGER_SHA256):
        raise RuntimeError("FROZEN_LEDGER_SNAPSHOT_CHANGED")
    for name in ("source_case", "source_manifest", "source_result", "material"):
        if plan.get(name + "_path") != str(paths[name]) or plan.get(name + "_sha256") != digest_file(paths[name]):
            raise RuntimeError("PLAN_SOURCE_CHANGED")
    if case["material_sha256"] != plan["material_sha256"]:
        raise RuntimeError("PLAN_MATERIAL_CHANGED")
    if digest_file(stage / "diagnostic_tools" / "slot_stability_probe.py") != plan["probe_sha256"]:
        raise RuntimeError("PROBE_CHANGED_AFTER_PREPARE")
    if digest_file(stage / Path(__file__).name) != plan["runner_sha256"]:
        raise RuntimeError("RUNNER_CHANGED_AFTER_PREPARE")
    budget_usage(ledger)


def growth_trigger_evidence(samples, diagnostic=DIAGNOSTIC):
    """Validate the first recorded growth trigger, not a later terminal snapshot.

    Meep can append a final sample one dt after the stop callback has fired.
    That sample is useful terminal evidence but is not a new two-unit guard
    observation. The witness belongs to the first explicitly marked trigger.
    """
    fixed = {"sample_interval": 2.0, "growth_ratio_limit": 1e12,
             "growth_step_factor": 10.0, "growth_consecutive_samples": 3}
    if any(diagnostic.get(key) != value for key, value in fixed.items()):
        raise RuntimeError("GROWTH_DIAGNOSTIC_CONTRACT_CHANGED")
    trigger = next((index for index, sample in enumerate(samples)
                    if sample.get("early_guard_status") == "EARLY_GROWTH_DETECTED"), None)
    if trigger is None or trigger < 2:
        raise RuntimeError("GROWTH_FIRST_TRIGGER_EVIDENCE_MISSING")
    witness = samples[trigger - 2:trigger + 1]
    times = [sample.get("simulation_time") for sample in witness]
    ratios = [sample.get("envelope", {}).get("log10_intensity_over_reference_peak")
              for sample in witness]

    def finite_number(value):
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)

    interval = diagnostic["sample_interval"]
    if (any(not finite_number(value) or value < 0 for value in times) or
            any(not math.isclose(value / interval, round(value / interval), rel_tol=0.0, abs_tol=1e-9)
                for value in times) or
            any(not math.isclose(b - a, interval, rel_tol=0.0, abs_tol=1e-9)
                for a, b in zip(times, times[1:]))):
        raise RuntimeError("GROWTH_FIRST_TRIGGER_SAMPLE_INTERVAL_INVALID")
    if (any(sample.get("envelope", {}).get("all_fields_finite") is not True for sample in witness) or
            any(point.get("finite") is False for sample in witness
                for point in sample.get("probes", {}).values()) or
            any(not finite_number(value) or value <= 12 for value in ratios) or
            any(b - a <= 1 for a, b in zip(ratios, ratios[1:]))):
        raise RuntimeError("GROWTH_FIRST_TRIGGER_THRESHOLDS_NOT_MET")
    return {"trigger_index": trigger, "trigger_time": times[-1],
            "witness_times": times, "witness_logratios": ratios}


def validate_probe_result(out, plan, case_id):
    out = Path(out)
    result = read_json(_regular_file(out / "result.json", "PROBE_RESULT"))
    expected = {"case_id": case_id, "eps_averaging": ARMS[case_id],
                "plan_sha256": plan["plan_sha256"], "qualification": False}
    for key in ("source_case_sha256", "implementation_sha256", "material_sha256", "probe_sha256"):
        expected[key] = plan[key]
    expected["physical_case_sha256"] = plan["source_case_sha256"]
    if any(result.get(key) != value for key, value in expected.items()):
        raise RuntimeError("PROBE_RESULT_PROVENANCE_MISMATCH")
    if result.get("status") not in CONTINUE_STATUSES | STOP_STATUSES:
        raise RuntimeError("UNKNOWN_PROBE_STATUS")
    evidence = result.get("evidence_sha256")
    required = {"trace.json", "geometry_metadata.json", "sample_grid.npz"}
    allowed = required | {f"peak_patch_{i:02d}.npz" for i in range(4)}
    if (not isinstance(evidence, dict) or not set(evidence).issubset(allowed) or
            (result["status"] in CONTINUE_STATUSES and not required.issubset(evidence))):
        raise RuntimeError("PROBE_EVIDENCE_INCOMPLETE")
    for name, fingerprint in evidence.items():
        if Path(name).name != name or name in ("result.json", "process.json", "status.json"):
            raise RuntimeError("PROBE_EVIDENCE_PATH_INVALID")
        path = _regular_file(out / name, "PROBE_EVIDENCE")
        if digest_file(path) != fingerprint:
            raise RuntimeError("PROBE_EVIDENCE_HASH_MISMATCH")
    if result["status"] in CONTINUE_STATUSES:
        trace = read_json(out / "trace.json")
        samples = trace.get("samples", [])
        summary = result.get("trace_summary", {})
        if (trace.get("case_id") != case_id or trace.get("eps_averaging") is not ARMS[case_id] or
                trace.get("qualification") is not False or not samples or
                summary.get("sample_count") != len(samples) or summary.get("last_sample") != samples[-1]):
            raise RuntimeError("PROBE_TRACE_INCOMPLETE")
        times = [sample.get("simulation_time") for sample in samples]
        if (any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in times) or
                abs(times[0]) > .01 or any(b <= a or b - a > 2.01 for a, b in zip(times, times[1:])) or
                abs(result.get("t_reached", -1) - times[-1]) > .01):
            raise RuntimeError("PROBE_TRACE_TIMES_INVALID")
        if result["status"] == "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION":
            if times[-1] < 100 or len(samples) < 51 or any(
                    sample.get("envelope", {}).get("all_fields_finite") is not True or
                    any(point.get("finite") is False for point in sample.get("probes", {}).values())
                    for sample in samples):
                raise RuntimeError("SHORT_WINDOW_NOT_COMPLETE")
        elif result["status"] == "NONFINITE_FIELD":
            final = samples[-1]
            if (final.get("envelope", {}).get("all_fields_finite") is not False and
                    not any(point.get("finite") is False for point in final.get("probes", {}).values())):
                raise RuntimeError("NONFINITE_EVIDENCE_MISSING")
        else:
            result["growth_trigger_evidence"] = growth_trigger_evidence(samples, plan["diagnostic"])
    return result


def _heartbeat_time(out):
    try:
        return read_json(out / "heartbeat.json").get("t_reached")
    except (OSError, ValueError):
        return None


def _process_identity(pid, proc_root="/proc"):
    """Read Linux identity without assuming comm has no spaces or parentheses."""
    try:
        fields = (Path(proc_root) / str(pid) / "stat").read_text().rsplit(")", 1)[1].split()
    except (FileNotFoundError, ProcessLookupError):
        return None
    return {"pid": int(pid), "state": fields[0], "pgid": int(fields[2]),
            "session": int(fields[3]), "start_ticks": int(fields[19])}


def session_live_members(session_id, proc_root="/proc"):
    """Include orphan ranks and separate process groups in our owned session."""
    root = Path(proc_root)
    if not root.is_dir():
        raise RuntimeError("PROCESS_SESSION_AUDIT_UNAVAILABLE")
    members = []
    for path in root.iterdir():
        if not path.name.isdigit():
            continue
        record = _process_identity(int(path.name), root)
        if record and record["session"] == session_id and record["state"] not in ("Z", "X", "x"):
            members.append(record)
    return members


def signal_session_member(member, session_id, sig):
    """Recheck start time and session just before signaling a specific PID."""
    current = _process_identity(member["pid"])
    if (current is None or current["session"] != session_id or
            current["start_ticks"] != member["start_ticks"] or current["state"] in ("Z", "X", "x")):
        return
    try:
        os.kill(current["pid"], sig)
    except ProcessLookupError:
        pass


def stop_owned_session(child, deadline, publish):
    """One absolute shutdown deadline, including TERM, KILL, reap and audit.

    The production helper only waits for mpiexec. Here MPI ranks can live in
    different process groups, so all non-zombie members of its newly created
    session must exit before the solver interval can close. No blocking wait
    can restart or extend the 30-second reserve.
    """
    started, errors, sent = time.monotonic(), set(), set()
    kill_at = min(started + 20.0, max(started, deadline - 5.0))
    next_write = started
    while True:
        child.poll()
        try:
            members = session_live_members(child.pid)
        except Exception as exc:
            errors.add(f"SESSION_AUDIT_FAILED: {type(exc).__name__}: {exc}")
            return {"stopped": False, "live_members": None, "errors": sorted(errors)}
        if child.returncode is not None and not members:
            return {"stopped": True, "live_members": [], "errors": sorted(errors),
                    "elapsed_s": time.monotonic() - started}
        now = time.monotonic()
        if now >= deadline:
            return {"stopped": False, "live_members": members, "errors": sorted(errors),
                    "elapsed_s": now - started}
        sig = signal.SIGKILL if now >= kill_at else signal.SIGTERM
        for member in members:
            if time.monotonic() >= deadline:
                break
            identity = (member["pid"], member["start_ticks"], int(sig))
            if identity in sent:
                continue
            try:
                signal_session_member(member, child.pid, sig)
                sent.add(identity)
            except Exception as exc:
                errors.add(f"PID {member['pid']} SIGNAL {int(sig)}: {type(exc).__name__}: {exc}")
        if now >= next_write:
            publish()
            next_write = now + 5.0
        time.sleep(min(.2, max(0.0, deadline - time.monotonic())))


def run_arm(stage, plan, ledger, case_id, cancelled):
    """One owned process group and one interval; accounting includes shutdown."""
    out = Path(plan["run_dir"]) / case_id
    if out.exists() or out.is_symlink():
        raise RuntimeError("EXISTING_DIAGNOSTIC_ATTEMPT_NO_RESUME")
    usage = budget_usage(ledger)
    allowance = min(PER_CASE_SECONDS, usage["effective_remaining_solver_seconds"])
    if cancelled[0] or allowance <= SHUTDOWN_RESERVE_SECONDS:
        return {"case_id": case_id, "status": "NOT_RUN", "reason": "interrupted or no solver allowance",
                "budget": usage, "qualification": False}
    resources = cgroup_resources()
    memory_budget = int(.7 * resources["memory_available_bytes"])
    if resources["cpu_quota"] < 4 or memory_budget < 1024**3:
        return {"case_id": case_id, "status": "NOT_RUN", "reason": "INSUFFICIENT_CGROUP_RESOURCES",
                "resources": resources, "budget": usage, "qualification": False}
    out.mkdir(exist_ok=False)
    command = [MPIEXEC, "-n", "4", PYTHON, "-B", "-m", "diagnostic_tools.slot_stability_probe",
               "--plan", str(Path(plan["run_dir"]) / "plan.json"),
               "--case-id", case_id, "--output", str(out)]
    env = {**os.environ, **{name: "1" for name in THREAD_ENV},
           "PYTHONDONTWRITEBYTECODE": "1",
           "TI2D_DIAGNOSTIC_TIMEOUT_SECONDS": str(allowance - SHUTDOWN_RESERVE_SECONDS)}
    atomic_json(out / "invocation.json", {"command": command, "plan_sha256": plan["plan_sha256"],
                "resources_before": resources, "memory_budget_bytes": memory_budget,
                "timeout_seconds_including_shutdown": allowance,
                "shutdown_reserve_seconds": SHUTDOWN_RESERVE_SECONDS,
                "thread_environment": {key: env[key] for key in THREAD_ENV},
                "solver_timeout_seconds": allowance - SHUTDOWN_RESERVE_SECONDS,
                "budget_before": usage})
    started, monotonic_start = time.time(), time.monotonic()
    child, index, peak, next_write = None, None, 0, 0.0
    status, error, result = "CRASHED", None, None
    process = {}

    def write_process(current_status):
        process.update({"run_id": plan["run_id"], "case_id": case_id, "status": current_status, "qualification": False,
                        "pid": child.pid if child is not None else None, "supervisor_pid": os.getpid(),
                        "start_epoch": started, "updated_epoch": time.time(),
                        "elapsed_s": time.monotonic() - monotonic_start,
                        "simulation_time": _heartbeat_time(out), "mpi_ranks": 4, "threads_per_rank": 1,
                        "peak_process_tree_rss_bytes": peak, "memory_budget_bytes": memory_budget,
                        "allowance_seconds_including_shutdown": allowance,
                        "exit_code": child.returncode if child is not None else None,
                        "error": error, "budget": budget_usage(ledger)})
        atomic_json(out / "process.json", process)
        atomic_json(Path(plan["run_dir"]) / "process.json", process)
        atomic_json(Path(plan["run_dir"]) / "status.json", {
            "run_id": plan["run_id"], "status": current_status, "active_case_id": case_id,
            "qualification": False, "updated_epoch": process["updated_epoch"],
            "simulation_time": process["simulation_time"], "budget": process["budget"]})

    try:
        index = ledger.add_interval(started, None, 4, plan["run_id"] + ":" + case_id, "RUNNING", attempt=1)
        with (out / "stdout.log").open("xb") as log:
            if cancelled[0]:
                status = "INTERRUPTED"
            else:
                child = subprocess.Popen(command, cwd=stage, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, start_new_session=True)
                write_process("RUNNING")
                while child.poll() is None:
                    peak = max(peak, process_tree_rss(child.pid))
                    elapsed = time.monotonic() - monotonic_start
                    if cancelled[0]:
                        status = "INTERRUPTED"
                    elif elapsed >= allowance - SHUTDOWN_RESERVE_SECONDS:
                        status = "TIME_LIMIT_UNQUALIFIED"
                    elif peak > memory_budget:
                        status = "MEMORY_LIMIT"
                    else:
                        if elapsed >= next_write:
                            write_process("RUNNING")
                            next_write = elapsed + 5.0
                        time.sleep(min(.5, max(.01, allowance - SHUTDOWN_RESERVE_SECONDS - elapsed)))
                        continue
                    write_process("STOPPING_" + status)
                    break
                if cancelled[0] and status == "CRASHED":
                    status = "INTERRUPTED"
                if status == "CRASHED":
                    try:
                        result = validate_probe_result(out, plan, case_id)
                        if child.returncode == 0 and result["status"] in CONTINUE_STATUSES:
                            status = result["status"]
                        elif result["status"] in STOP_STATUSES:
                            status = result["status"]
                        else:
                            error = "PROBE_EXIT_STATUS_MISMATCH"
                    except Exception as exc:
                        error = f"{type(exc).__name__}: {exc}"
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    finally:
        stopped = child is None
        try:
            if child is not None:
                child.poll()
                live_members = session_live_members(child.pid)
                stopped = child.returncode is not None and not live_members
                if not stopped:
                    if child.returncode is not None and live_members:
                        status = "CRASHED"
                        error = f"{error or ''}; MPI_LAUNCHER_EXITED_WITH_LIVE_SESSION_MEMBERS"
                    deadline = min(monotonic_start + allowance,
                                   time.monotonic() + SHUTDOWN_RESERVE_SECONDS)
                    process["shutdown_deadline_monotonic"] = deadline
                    cleanup = stop_owned_session(child, deadline, lambda: write_process("STOPPING_" + status))
                    process["session_cleanup"] = cleanup
                    stopped = cleanup["stopped"]
        except Exception as exc:
            status = "CRASHED"
            error = f"{error or ''}; CLEANUP_ERROR: {type(exc).__name__}: {exc}"
        finally:
            if not stopped:
                status = "CRASHED"
                error = f"{error or ''}; CHILD_STOP_UNCONFIRMED: ledger interval remains open"
                process["cleanup_unconfirmed"] = True
            if index is not None and stopped:
                ledger.close_interval(index, time.time(), status)
            write_process(status)
    process["evidence_validated"] = result is not None and status in CONTINUE_STATUSES
    if result is not None:
        process["probe_result"] = result
        process["result_sha256"] = digest_file(out / "result.json")
    atomic_json(out / "process.json", process)
    return deepcopy(process)


def report_for(plan, arms, ledger):
    """A successful short window supports only a bounded diagnostic statement."""
    by_id = {entry["case_id"]: entry for entry in arms}
    for case_id in ARMS:
        by_id.setdefault(case_id, {"case_id": case_id, "status": "NOT_RUN",
                                   "reason": "earlier arm stopped or budget unavailable", "qualification": False})
    on, off = (by_id[key] for key in ARMS)
    both_complete = all(entry.get("evidence_validated") is True and entry["status"] in CONTINUE_STATUSES
                        for entry in (on, off))
    if not both_complete:
        interpretation = "INCOMPLETE_EVIDENCE_NO_ATTRIBUTION"
        conclusion = "The two-arm evidence is incomplete; no fix or causal attribution is established."
    elif on["status"] == "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION":
        interpretation = "ORIGINAL_FAILURE_NOT_REPRODUCED_NO_ATTRIBUTION"
        conclusion = "Averaging on reached t=100 without reproducing the original failure; this comparison cannot attribute that failure to averaging."
    elif off["status"] == "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION":
        interpretation = "AVERAGING_OFF_SHORT_WINDOW_STABILITY_ONLY"
        conclusion = "Averaging on showed instability and averaging off completed t=100. This supports an averaging-related short-window instability; it does not establish a production fix or convergence."
    else:
        interpretation = "BOTH_ARMS_UNSTABLE_NO_FIX_ESTABLISHED"
        conclusion = "Both arms showed instability; disabling averaging did not establish short-window stability."
    return {"schema_version": 1, "run_id": plan["run_id"], "created_utc": utcnow(),
            "plan_sha256": plan["plan_sha256"], "qualification": False,
            "formal_R_T_A_conclusion": False, "automatic_retry_authorized": False,
            "arms": list(by_id.values()), "both_arms_complete_evidence": both_complete,
            "interpretation": interpretation, "conclusion": conclusion,
            "budget": budget_usage(ledger), "ledger": ledger.as_dict(),
            "limitations": "Short t<=100 diagnostic only; source may still be on. No official R/T/A, qualification, full decay, convergence, or production repair is inferred."}


def write_report(run, report):
    atomic_json(run / "report.json", report)
    lines = ["# Slot stability short diagnostic", "", report["conclusion"], "",
             "| Arm | Status | Simulation time | Wall seconds | Evidence checked |",
             "| --- | --- | ---: | ---: | --- |"]
    for entry in report["arms"]:
        simulated = entry.get("probe_result", {}).get("t_reached", entry.get("simulation_time"))
        lines.append(f"| {entry['case_id']} | {entry['status']} | {simulated} | {entry.get('elapsed_s', 0):.3f} | {entry.get('evidence_validated', False)} |")
    budget = report["budget"]
    lines.extend(["", f"Phase active solver time: {budget['phase_active_solver_seconds']:.3f} / {PHASE_CAP_SECONDS:.0f} s. "
                  f"Stage: {budget['stage_active_solver_seconds']:.3f} / {STAGE_CAP_SECONDS:.0f} s. "
                  f"Total rank-hours: {budget['rank_hours']:.6f}.", "", report["limitations"], ""])
    # The run is new and single-use; no production document is overwritten.
    with (run / "report.md").open("x", encoding="utf-8") as stream:
        stream.write("\n".join(lines))


def run(stage, plan_path=None):
    stage = Path(stage).resolve()
    if sys.platform != "linux" or stage != SERVER_STAGE:
        raise RuntimeError("SERVER_ONLY: no local FDTD execution")
    run_dir = stage / "repair_diagnostics" / RUN_ID
    if run_dir.resolve() != run_dir:
        raise RuntimeError("DIAGNOSTIC_RUN_SYMLINK_PATH_PROHIBITED")
    expected_plan = run_dir / "plan.json"
    if plan_path is not None and Path(plan_path).resolve() != expected_plan:
        raise RuntimeError("FIXED_PLAN_PATH_REQUIRED")
    _regular_file(expected_plan, "PREPARED_PLAN")
    with FileLock(stage / "budget_ledger.lock"), FileLock(run_dir / "run.lock"):
        ledger = closed_ledger(stage)
        plan, authorization = read_json(expected_plan), read_json(run_dir / "authorization.json")
        verify_frozen(stage, plan, authorization, ledger)
        if ((run_dir / "started.json").exists() or (run_dir / "report.json").exists() or
                any((run_dir / case_id).exists() for case_id in ARMS)):
            raise RuntimeError("EXISTING_DIAGNOSTIC_ATTEMPT_NO_RESUME")
        cancelled, arms = [False], []

        def cancel(*unused):
            cancelled[0] = True

        previous = {sig: signal.signal(sig, cancel) for sig in (signal.SIGTERM, signal.SIGINT)}
        atomic_json(run_dir / "started.json", {"supervisor_pid": os.getpid(), "started_utc": utcnow(),
                                              "plan_sha256": plan["plan_sha256"], "automatic_resume": False})
        try:
            for case_id in ARMS:
                if cancelled[0]:
                    break
                arm = run_arm(stage, plan, ledger, case_id, cancelled)
                arms.append(arm)
                if arm["status"] not in CONTINUE_STATUSES:
                    break
        except Exception as exc:
            arms.append({"case_id": next((name for name in ARMS if name not in {a['case_id'] for a in arms}), "supervisor"),
                         "status": "CRASHED", "error": f"{type(exc).__name__}: {exc}", "qualification": False})
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
            report = report_for(plan, arms, ledger)
            write_report(run_dir, report)
            terminal = ("CHILD_STOP_UNCONFIRMED" if any(a.get("cleanup_unconfirmed") for a in arms)
                        else "FINISHED_NOT_QUALIFICATION")
            atomic_json(run_dir / "status.json", {"run_id": RUN_ID, "status": terminal,
                "interpretation": report["interpretation"], "qualification": False,
                "updated_utc": utcnow(), "budget": report["budget"],
                "arms": {entry["case_id"]: entry["status"] for entry in report["arms"]}})
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
        parser.error("--prepare uses the fixed plan path and does not accept --plan")
    result = prepare(args.stage_root) if args.prepare else run(args.stage_root, args.plan)
    print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
    return 0 if args.prepare or result["both_arms_complete_evidence"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
