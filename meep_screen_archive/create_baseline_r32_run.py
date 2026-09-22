"""Freeze the one approved baseline r32 condition and shared budget; never solve."""
from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import sys

from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.queue import FileLock, SolverLedger, read_json, validate_result, verify_manifest

RUN_ID = "baseline_r32_20260922_E"
SOURCE_RUN_ID = "flat_r32_20260922_D"
SOURCE_TASK_ID = "flat_r32_p_plus30"
TASK_ID = "baseline_r32_p_plus30"
BASELINE_ACTIVE_SECONDS = 30121.365798711777
PHASE_CAP_SECONDS = 21600.0
STAGE_CAP_SECONDS = 172800.0
APPROVED_SOURCE_CASE_SHA256 = "d1b293573cf3ee9a28a24e00f57341bdcf432d120fec9e8f4aa645c38c54ccb2"
APPROVED_LEDGER_SHA256 = "9390fd237e0a5093d7c933539aa6a193f9a3134e96e9d5bd049b1080b7baa557"
APPROVED_IMPLEMENTATION_SHA256 = "24d09c50dcf6c0cb9271947dbcb04f03b62d9d0500bd7bace34e6909d54834c1"
SOURCE_EVIDENCE_NAMES = frozenset((
    "monitor_planes.npz", "volume_ED.npz", "volume_metadata.json", "orders.json",
    "structure_stop.json", "reference_stop.json", "reference_preload_audit.json",
    "structure_chunks.h5",
))
SOURCE_CHECK_NAMES = frozenset((
    "stopped", "case_known", "independent_absorption", "reflection_orders",
    "transmission_orders", "reference_pseudo_reflection", "finite",
    "passive_flux_within_tolerance", "fresnel_R", "fresnel_RTA",
))


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
    """Call under budget_ledger.lock; allow appends, never rewritten history."""
    ledger = closed_ledger(stage)
    expected = {
        "run_id": RUN_ID, "task_id": TASK_ID, "ledger_path": str(ledger.path),
        "baseline_active_solver_seconds": BASELINE_ACTIVE_SECONDS,
        "phase_active_solver_cap_seconds": PHASE_CAP_SECONDS,
        "stage_total_active_seconds": STAGE_CAP_SECONDS,
        "effective_total_active_seconds": min(STAGE_CAP_SECONDS, BASELINE_ACTIVE_SECONDS + PHASE_CAP_SECONDS),
        "ledger_file_sha256": APPROVED_LEDGER_SHA256,
        "implementation_sha256": APPROVED_IMPLEMENTATION_SHA256,
        "source_case_sha256": APPROVED_SOURCE_CASE_SHA256,
        "automatic_retry_authorized": False,
    }
    if any(authorization.get(key) != value for key, value in expected.items()):
        raise RuntimeError("BASELINE_R32_AUTHORIZATION_CHANGED")
    snapshot = authorization["ledger_snapshot"]
    if digest_object(snapshot) != authorization["ledger_snapshot_sha256"]:
        raise RuntimeError("FROZEN_LEDGER_SNAPSHOT_CHANGED")
    prefix = snapshot["intervals"]
    if (ledger.data["intervals"][:len(prefix)] != prefix or
            {k: v for k, v in ledger.data.items() if k != "intervals"} !=
            {k: v for k, v in snapshot.items() if k != "intervals"}):
        raise RuntimeError("HISTORICAL_LEDGER_CHANGED")
    if ledger.total_active_seconds() < BASELINE_ACTIVE_SECONDS - 1e-6:
        raise RuntimeError("SHARED_LEDGER_BELOW_FROZEN_BASELINE")
    return ledger


def verify_source_result(result, case, result_path):
    checks = result.get("checks")
    if (result.get("status") != "QUALIFIED" or not isinstance(checks, dict) or
            not SOURCE_CHECK_NAMES.issubset(checks) or
            any(value is not True for value in checks.values()) or
            result.get("id") != SOURCE_TASK_ID or result.get("case_config") != case or
            result.get("case_sha256") != APPROVED_SOURCE_CASE_SHA256 or
            result.get("implementation_sha256") != APPROVED_IMPLEMENTATION_SHA256):
        raise ValueError("SOURCE_RESULT_NOT_APPROVED_QUALIFIED")
    validate_result(result)
    evidence = result.get("evidence_sha256")
    if not isinstance(evidence, dict) or set(evidence) != SOURCE_EVIDENCE_NAMES:
        raise ValueError("SOURCE_EVIDENCE_SET_INCOMPLETE_OR_UNEXPECTED")
    for name, expected_hash in evidence.items():
        path = Path(result_path).parent / name
        if path.is_symlink() or not path.is_file() or digest_file(path) != expected_hash:
            raise ValueError(f"SOURCE_EVIDENCE_FINGERPRINT_MISMATCH:{name}")
    return deepcopy(evidence)


def create(stage):
    stage = Path(stage).resolve()
    run, source = stage / "runs" / RUN_ID, stage / "runs" / SOURCE_RUN_ID
    with FileLock(stage / "budget_ledger.lock"):
        if run.exists():
            raise FileExistsError("BASELINE_R32_RUN_ALREADY_EXISTS: no overwrite or automatic retry")
        ledger = closed_ledger(stage)
        if abs(ledger.total_active_seconds() - BASELINE_ACTIVE_SECONDS) > 1e-6:
            raise RuntimeError("APPROVED_BASELINE_MISMATCH")
        if digest_file(ledger.path) != APPROVED_LEDGER_SHA256:
            raise RuntimeError("APPROVED_LEDGER_CHANGED")
        implementation = implementation_digest()
        if implementation != APPROVED_IMPLEMENTATION_SHA256:
            raise ValueError("APPROVED_IMPLEMENTATION_CHANGED")
        manifest_path = source / "manifest.json"
        prior = read_json(manifest_path)
        verify_manifest(prior, source, implementation)
        if (prior.get("run_id") != SOURCE_RUN_ID or len(prior["tasks"]) != 1 or
                prior["tasks"][0]["id"] != SOURCE_TASK_ID):
            raise ValueError("EXACT_SOURCE_TASK_REQUIRED")
        source_task = prior["tasks"][0]
        source_case_path = Path(source_task["case_path"])
        if not source_case_path.is_absolute():
            source_case_path = source / source_case_path
        original_case = read_json(source_case_path)
        if digest_file(source_case_path) != APPROVED_SOURCE_CASE_SHA256:
            raise ValueError("UNEXPECTED_SOURCE_PHYSICS_OR_NUMERICS")
        required_case = {"id": SOURCE_TASK_ID, "case": "flat", "geometry_id": "baseline",
                         "resolution": 32, "polarization": "p", "emission_theta_deg": 30.0,
                         "wavelength_um": 10.5, "mpi_ranks": 4, "max_solver_seconds": 10800.0,
                         "strict_stop": False}
        if any(original_case.get(key) != value for key, value in required_case.items()):
            raise ValueError("UNEXPECTED_SOURCE_PHYSICS_OR_NUMERICS")
        material = Path(original_case["material_path"])
        if not material.is_absolute():
            material = stage / material
        if (digest_file(material) != original_case["material_sha256"] or
                prior.get("material_sha256") != original_case["material_sha256"]):
            raise ValueError("FROZEN_MATERIAL_CHANGED")
        result_path = source / "tasks" / SOURCE_TASK_ID / "attempt_1" / "result.json"
        result = read_json(result_path)
        source_evidence = verify_source_result(result, original_case, result_path)
        case = dict(deepcopy(original_case), id=TASK_ID, case="structure", max_solver_seconds=PHASE_CAP_SECONDS)
        diff = {key: {"before": original_case[key], "after": case[key]}
                for key in original_case if original_case[key] != case[key]}
        if set(diff) != {"id", "case", "max_solver_seconds"}:
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
            "source_evidence_sha256": source_evidence,
            "implementation_sha256": implementation, "material_sha256": original_case["material_sha256"],
            "case_diff": diff, "automatic_retry_authorized": False,
            "budget_meaning": "21600 s cumulative active solver wall time for reference plus target, including all failed or interrupted attempts; also inside original 172800 s shared stage cap",
            "scope": "One baseline slanted-slot r32 p +30 condition at 10.5 um only; no pilot, concurrency probe, mesh sweep, automatic retry, or full-validation qualification",
        }
        verify_authorization(stage, authorization)
        run.mkdir(parents=True, exist_ok=False)
        case_path = run / "cases" / (TASK_ID + ".json")
        atomic_json(case_path, case)
        task = {
            "id": TASK_ID, "phase": "validation", "depends_on": [], "gate_role": "baseline",
            "case_path": str(case_path), "geometry_id": case["geometry_id"], "resolution": case["resolution"],
            "polarization": case["polarization"], "emission_theta_deg": case["emission_theta_deg"],
            "tilt_deg": case["geometry"]["tilt_deg"], "strict_stop": case["strict_stop"],
            "case_fingerprint": digest_object({key: value for key, value in case.items() if key != "id"}),
        }
        manifest = {
            "schema_version": 1, "run_id": RUN_ID, "created_utc": authorization["created_utc"],
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
                        "memory_fraction": .7, "mpi_ranks": case["mpi_ranks"]},
        }
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
