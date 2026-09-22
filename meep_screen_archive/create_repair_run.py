"""Create one immutable flat-repair manifest; never launch a solver."""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
import re
import sys

from repair_budget import locked_repair_budget, repair_usage
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.queue import FileLock, SolverLedger, verify_manifest


def create(stage, run_id="repair_20260922_C", source_run="stage_20260921_B", overrides=None):
    stage = Path(stage).resolve()
    if not all(re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in (run_id, source_run)) or run_id == source_run:
        raise ValueError("NEW_RUN_ID_REQUIRED")
    run, source = stage / "runs" / run_id, stage / "runs" / source_run
    if run.exists():
        raise FileExistsError("REPAIR_RUN_ALREADY_EXISTS: do not overwrite or automatically retry")
    if not (stage / "budget_ledger.json").is_file():
        raise RuntimeError("EXISTING_SHARED_LEDGER_REQUIRED")
    with FileLock(stage / "budget_ledger.lock"):
        ledger = SolverLedger(stage / "budget_ledger.json")
        budget = locked_repair_budget(stage, ledger)
        usage = repair_usage(ledger, budget)
        if usage["effective_remaining_solver_seconds"] <= 30:
            raise RuntimeError("REPAIR_BUDGET_EXHAUSTED")
        prior = json.loads((source / "manifest.json").read_text())
        selected = [t for t in prior["tasks"] if t["id"] == budget["allowed_flat_task_id"]]
        if len(selected) != 1:
            raise ValueError("EXACT_FLAT_TASK_REQUIRED")
        task = deepcopy(selected[0])
        source_case = Path(task["case_path"])
        if not source_case.is_absolute():
            source_case = source / source_case
        case = json.loads(source_case.read_text())
        if not (case["case"] == "flat" and case["resolution"] == 16 and case["polarization"] == "p"
                and case["emission_theta_deg"] == 30 and case["mpi_ranks"] == 4):
            raise ValueError("FLAT_CASE_SCOPE_MISMATCH")
        overrides = dict(overrides or {})
        protected = {"id", "case", "geometry_id", "geometry", "resolution", "polarization",
                     "emission_theta_deg", "mpi_ranks", "material_path", "material_sha256",
                     "wavelength_um", "courant", "pml_um", "air_above_um", "air_below_um",
                     "source_fwidth_fraction", "source_cutoff", "oxide_layer_in_model"}
        if protected.intersection(overrides):
            raise ValueError("REPAIR_MUST_PRESERVE_PHYSICS")
        case.update(overrides)
        case.setdefault('reference_intensity_floor_fraction', 1.0)
        case.setdefault('stop_normalization_policy', 'max(own_historic_peak, reference_source_on_non_pml_peak)')
        case.setdefault('mpi_partition_policy', 'fixed_binary_Y0_X0_ranks0123')
        implementation = implementation_digest()
        case.update(implementation_sha256=implementation,
                    max_solver_seconds=min(float(case.get("max_solver_seconds", 21600)),
                                           budget["repair_active_solver_cap_seconds"]))
        if digest_file(case["material_path"]) != case["material_sha256"]:
            raise ValueError("FROZEN_MATERIAL_CHANGED")
        run.mkdir(parents=True, exist_ok=False)
        case_path = run / "cases" / (task["id"] + ".json")
        atomic_json(case_path, case)
        task.pop("reused_from", None)
        task.update(phase="validation", depends_on=[], case_path=str(case_path),
                    case_fingerprint=digest_object({k: v for k, v in case.items() if k != "id"}))
        manifest = deepcopy(prior)
        manifest.pop("manifest_sha256", None)
        manifest.pop("concurrency_probe", None)
        manifest.update(run_id=run_id, created_utc=utcnow(), run_dir=str(run), stage_root=str(stage),
                        implementation_sha256=implementation, ledger_path=str(ledger.path.resolve()),
                        tasks=[task], production_gate_task_ids=[task["id"]], gate_comparisons=[],
                        pilot_design_conditions=0, repairs_performed=1,
                        scope="Targeted r16 p +30 flat repair only; full validation and pilot remain unqualified",
                        repair_budget=budget, repair_usage_at_creation=usage,
                        source_manifest_path=str(source / "manifest.json"),
                        source_manifest_sha256=digest_file(source / "manifest.json"),
                        source_case_sha256=digest_file(source_case))
        manifest["budgets"] = {"total_active_seconds": budget["effective_total_active_seconds"],
                               "stage_total_active_seconds": budget["stage_total_active_seconds"],
                               "single_solver_seconds": budget["repair_active_solver_cap_seconds"],
                               "max_parallel": 1, "memory_fraction": .7, "mpi_ranks": 4}
        manifest["manifest_sha256"] = digest_object(manifest)
        verify_manifest(manifest, run, implementation)
        atomic_json(run / "manifest.json", manifest)
        return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="repair_20260922_C")
    parser.add_argument("--source-run", default="stage_20260921_B")
    parser.add_argument("--case-overrides-json", type=Path)
    args = parser.parse_args()
    stage = Path(__file__).resolve().parent
    if sys.platform != "linux" or str(stage) != "/root/autodl-tmp/meep_sim/meep_screen":
        raise RuntimeError("SERVER_ONLY")
    overrides = json.loads(args.case_overrides_json.read_text()) if args.case_overrides_json else None
    result = create(stage, args.run_id, args.source_run, overrides)
    print(json.dumps({"run_id": result["run_id"], "manifest": str(Path(result["run_dir"]) / "manifest.json"),
                      "task_ids": [t["id"] for t in result["tasks"]],
                      "repair_usage": result["repair_usage_at_creation"],
                      "full_validation_qualified": False, "solver_launched": False}, indent=2))


if __name__ == "__main__":
    main()
