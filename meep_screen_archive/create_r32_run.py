"""Freeze the one approved r32 flat test and its shared allowance; never solve."""
from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import sys

from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.queue import FileLock, SolverLedger, read_json, verify_manifest

RUN_ID = "flat_r32_20260922_D"
SOURCE_RUN_ID = "repair_20260922_C"
SOURCE_TASK_ID = "flat_r16_p_plus30"
TASK_ID = "flat_r32_p_plus30"
BASELINE_ACTIVE_SECONDS = 22581.670486450195
PHASE_CAP_SECONDS = 10800.0
STAGE_CAP_SECONDS = 172800.0
# Byte hash of the reviewed, unchanged repair r16 case; not a regenerated config.
APPROVED_SOURCE_CASE_SHA256 = "d71c30b3517d0edb068b1c5649671bf539c589b2aa3e6710883eee3e38e4a03d"
APPROVED_LEDGER_SHA256 = "bfc77312db537e2c4ee4639ac41f870eb60dae9672fe37cdf9208a383738185a"


def closed_ledger(stage):
    path = Path(stage).resolve() / "budget_ledger.json"
    if not path.is_file():
        raise RuntimeError("EXISTING_SHARED_LEDGER_REQUIRED")
    ledger = SolverLedger(path)
    for record in ledger.data["intervals"]:
        if record.get("end") is None:
            raise RuntimeError("UNCLOSED_PRIOR_INTERVAL")
        if not all(math.isfinite(record[key]) for key in ("start", "end")):
            raise RuntimeError("NONFINITE_LEDGER_INTERVAL")
    return ledger


def verify_authorization(stage, authorization):
    """Call under budget_ledger.lock; permit appends but never alter frozen history."""
    ledger = closed_ledger(stage)
    snapshot = authorization["ledger_snapshot"]
    expected_cap = min(STAGE_CAP_SECONDS, BASELINE_ACTIVE_SECONDS + PHASE_CAP_SECONDS)
    expected = {"baseline_active_solver_seconds": BASELINE_ACTIVE_SECONDS,
                "phase_active_solver_cap_seconds": PHASE_CAP_SECONDS,
                "stage_total_active_seconds": STAGE_CAP_SECONDS,
                "effective_total_active_seconds": expected_cap,
                "ledger_path": str(ledger.path), "run_id": RUN_ID, "task_id": TASK_ID}
    if any(authorization.get(key) != value for key, value in expected.items()):
        raise RuntimeError("R32_AUTHORIZATION_CHANGED")
    if digest_object(snapshot) != authorization["ledger_snapshot_sha256"]:
        raise RuntimeError("FROZEN_LEDGER_SNAPSHOT_CHANGED")
    prefix = snapshot["intervals"]
    if (ledger.data["intervals"][:len(prefix)] != prefix or
            {k: v for k, v in ledger.data.items() if k != "intervals"} !=
            {k: v for k, v in snapshot.items() if k != "intervals"}):
        raise RuntimeError("HISTORICAL_LEDGER_CHANGED")
    total = ledger.total_active_seconds()
    if total < BASELINE_ACTIVE_SECONDS - 1e-6:
        raise RuntimeError("SHARED_LEDGER_BELOW_FROZEN_BASELINE")
    return ledger


def create(stage):
    stage = Path(stage).resolve()
    run, source = stage / "runs" / RUN_ID, stage / "runs" / SOURCE_RUN_ID
    with FileLock(stage / "budget_ledger.lock"):
        if run.exists():
            raise FileExistsError("R32_RUN_ALREADY_EXISTS: no overwrite or automatic retry")
        ledger = closed_ledger(stage)
        if abs(ledger.total_active_seconds() - BASELINE_ACTIVE_SECONDS) > 1e-6:
            raise RuntimeError("APPROVED_BASELINE_MISMATCH")
        if digest_file(ledger.path) != APPROVED_LEDGER_SHA256:
            raise RuntimeError("APPROVED_LEDGER_CHANGED")
        manifest_path = source / "manifest.json"
        prior = read_json(manifest_path)
        implementation = implementation_digest()
        verify_manifest(prior, source, implementation)
        if len(prior["tasks"]) != 1 or prior["tasks"][0]["id"] != SOURCE_TASK_ID:
            raise ValueError("EXACT_SOURCE_TASK_REQUIRED")
        task = deepcopy(prior["tasks"][0])
        source_case_path = Path(task["case_path"])
        if not source_case_path.is_absolute():
            source_case_path = source / source_case_path
        original_case = read_json(source_case_path)
        if digest_file(source_case_path) != APPROVED_SOURCE_CASE_SHA256:
            raise ValueError("UNEXPECTED_SOURCE_PHYSICS_OR_NUMERICS")
        material = Path(original_case["material_path"])
        if not material.is_absolute():
            material = stage / material
        if digest_file(material) != original_case["material_sha256"]:
            raise ValueError("FROZEN_MATERIAL_CHANGED")
        result_path = source / "tasks" / SOURCE_TASK_ID / "attempt_1" / "result.json"
        result = read_json(result_path)
        checks = result.get("checks", {})
        if (result.get("status") != "NUMERICAL_FAILED" or checks.get("stopped") is not True or
                {key for key, value in checks.items() if value is not True} != {"fresnel_R"} or
                checks.get("fresnel_R") is not False or result.get("case_config") != original_case):
            raise ValueError("SOURCE_RESULT_NOT_APPROVED_FRESNEL_FAILURE")
        case = dict(original_case, id=TASK_ID, resolution=32, max_solver_seconds=PHASE_CAP_SECONDS)
        diff = {key: {"before": original_case[key], "after": case[key]}
                for key in original_case if original_case[key] != case[key]}
        if set(diff) != {"id", "resolution", "max_solver_seconds"}:
            raise ValueError("UNAUTHORIZED_CASE_DIFF")
        authorization = {
            "schema_version": 1, "run_id": RUN_ID, "task_id": TASK_ID, "created_utc": utcnow(),
            "baseline_active_solver_seconds": BASELINE_ACTIVE_SECONDS,
            "phase_active_solver_cap_seconds": PHASE_CAP_SECONDS,
            "stage_total_active_seconds": STAGE_CAP_SECONDS,
            "effective_total_active_seconds": min(STAGE_CAP_SECONDS, BASELINE_ACTIVE_SECONDS + PHASE_CAP_SECONDS),
            "ledger_path": str(ledger.path), "ledger_snapshot": deepcopy(ledger.data),
            "ledger_snapshot_sha256": digest_object(ledger.data), "ledger_file_sha256": digest_file(ledger.path),
            "source_manifest_path": str(manifest_path), "source_manifest_sha256": digest_file(manifest_path),
            "source_case_path": str(source_case_path), "source_case_sha256": digest_file(source_case_path),
            "source_result_path": str(result_path), "source_result_sha256": digest_file(result_path),
            "implementation_sha256": implementation, "material_sha256": original_case["material_sha256"],
            "case_diff": diff, "automatic_retry_authorized": False,
            "budget_meaning": "10800 s cumulative active solver wall time for reference plus target, including all failed or interrupted attempts; also inside original 172800 s shared stage cap",
            "scope": "One flat r32 p +30 condition only; no pilot, mesh sweep, or full-validation qualification",
        }
        verify_authorization(stage, authorization)
        run.mkdir(parents=True, exist_ok=False)
        case_path = run / "cases" / (TASK_ID + ".json")
        atomic_json(case_path, case)
        task.pop("reused_from", None)
        task.update(id=TASK_ID, phase="validation", depends_on=[], resolution=32, case_path=str(case_path),
                    case_fingerprint=digest_object({k: v for k, v in case.items() if k != "id"}))
        manifest = {"schema_version": 1, "run_id": RUN_ID, "created_utc": authorization["created_utc"],
                    "stage_root": str(stage), "run_dir": str(run), "ledger_path": str(ledger.path),
                    "implementation_sha256": implementation, "material_path": original_case["material_path"],
                    "material_sha256": original_case["material_sha256"], "tasks": [task],
                    "production_gate_task_ids": [TASK_ID], "gate_comparisons": [], "pilot_design_conditions": 0,
                    "scope": authorization["scope"], "phase_budget": {key: authorization[key] for key in (
                        "baseline_active_solver_seconds", "phase_active_solver_cap_seconds", "stage_total_active_seconds",
                        "effective_total_active_seconds", "budget_meaning")},
                    "authorization_path": str(run / "authorization.json"),
                    "authorization_sha256": digest_object(authorization),
                    "budgets": {"total_active_seconds": authorization["effective_total_active_seconds"],
                                "stage_total_active_seconds": STAGE_CAP_SECONDS,
                                "single_solver_seconds": PHASE_CAP_SECONDS, "max_parallel": 1,
                                "memory_fraction": .7, "mpi_ranks": 4}}
        manifest["manifest_sha256"] = digest_object(manifest)
        verify_manifest(manifest, run, implementation)
        atomic_json(run / "authorization.json", authorization)
        atomic_json(run / "config_diff.json", diff)
        atomic_json(run / "manifest.json", manifest)
        return manifest


def main():
    stage = Path(__file__).resolve().parent
    if sys.platform != "linux" or str(stage) != "/root/autodl-tmp/meep_sim/meep_screen":
        raise RuntimeError("SERVER_ONLY")
    manifest = create(stage)
    print(json.dumps({"run_id": RUN_ID, "task_ids": [TASK_ID], "phase_budget": manifest["phase_budget"],
                      "manifest": str(Path(manifest["run_dir"]) / "manifest.json"), "solver_launched": False}, indent=2))


if __name__ == "__main__":
    main()
