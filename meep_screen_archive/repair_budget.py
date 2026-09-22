"""Frozen repair allowance within the unchanged cumulative stage ledger.

Call ``locked_repair_budget`` while holding ``budget_ledger.lock``.
No solver imports or import-time side effects.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import tarfile

from ti2d.queue import atomic_json, sha256, stable_hash

BASELINE_ACTIVE_SECONDS = 21608.23800969124
STAGE_CAP_SECONDS = 172800.0
REPAIR_CAP_SECONDS = 1800.0
SNAPSHOT_RELATIVE = Path("snapshots/20260922T015134Z_repair1/pre_repair_code_and_state.tar.gz")


def _closed_union(records):
    merged = []
    for record in records:
        if record.get("end") is None:
            raise RuntimeError("UNCLOSED_BASELINE_INTERVAL")
        if not all(math.isfinite(record[k]) for k in ("start", "end")):
            raise RuntimeError("NONFINITE_BASELINE_INTERVAL")
    for start, end in sorted((r["start"], r["end"]) for r in records):
        if end < start:
            raise RuntimeError("REVERSED_BASELINE_INTERVAL")
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return sum(b - a for a, b in merged)


def locked_repair_budget(stage, ledger):
    """Freeze/verify repair1.json and the exact historical ledger prefix."""
    stage = Path(stage).resolve()
    if ledger.path.resolve() != stage / "budget_ledger.json" or not ledger.path.is_file():
        raise RuntimeError("EXISTING_SHARED_LEDGER_REQUIRED")
    if any(r.get("end") is None for r in ledger.data["intervals"]):
        raise RuntimeError("UNCLOSED_PRIOR_INTERVAL: reconcile the owned process before launch")
    snapshot = stage / SNAPSHOT_RELATIVE
    if not snapshot.is_file():
        raise RuntimeError("PRE_REPAIR_SNAPSHOT_REQUIRED")
    with tarfile.open(snapshot, "r:gz") as archive:
        members = [m for m in archive.getmembers()
                   if m.isfile() and Path(m.name).name == "budget_ledger.json"]
        if len(members) != 1:
            raise RuntimeError("AMBIGUOUS_SNAPSHOT_LEDGER")
        with archive.extractfile(members[0]) as stream:
            original = json.load(stream)
    if original.get("version") != 1 or not isinstance(original.get("intervals"), list):
        raise RuntimeError("INVALID_SNAPSHOT_LEDGER")
    baseline = original["intervals"]
    if abs(_closed_union(baseline) - BASELINE_ACTIVE_SECONDS) > 1e-6:
        raise RuntimeError("PRE_REPAIR_BASELINE_MISMATCH")
    if ledger.data["intervals"][:len(baseline)] != baseline:
        raise RuntimeError("HISTORICAL_LEDGER_CHANGED")
    expected = {
        "schema_version": 1, "repair_id": "repair1", "repair_round": 1,
        "approved_max_targeted_repair_rounds": 2,
        "baseline_active_solver_seconds": BASELINE_ACTIVE_SECONDS,
        "stage_total_active_seconds": STAGE_CAP_SECONDS,
        "repair_active_solver_cap_seconds": REPAIR_CAP_SECONDS,
        "effective_total_active_seconds": min(STAGE_CAP_SECONDS, BASELINE_ACTIVE_SECONDS + REPAIR_CAP_SECONDS),
        "ledger_path": str(ledger.path.resolve()), "source_snapshot": str(snapshot),
        "source_snapshot_sha256": sha256(snapshot), "baseline_interval_count": len(baseline),
        "baseline_intervals_sha256": stable_hash(baseline), "allowed_flat_task_id": "flat_r16_p_plus30",
        "scope": "Nonzero-cache diagnostic and one r16 p +30 flat condition only; no full validation or pilot",
        "budget_meaning": "1800 s combined active solver wall time, including failed and interrupted attempts; no convergence promise",
    }
    expected["budget_sha256"] = stable_hash(expected)
    path = stage / "repair1.json"
    if path.exists():
        with path.open() as stream:
            if json.load(stream) != expected:
                raise RuntimeError("FROZEN_REPAIR_BUDGET_CHANGED")
    else:
        atomic_json(path, expected)
    return expected


def repair_usage(ledger, budget):
    total = ledger.total_active_seconds()
    baseline = budget["baseline_active_solver_seconds"]
    if total < baseline - 1e-6:
        raise RuntimeError("SHARED_LEDGER_BELOW_FROZEN_BASELINE")
    used = max(0.0, total - baseline)
    return {
        "stage_active_solver_seconds": total, "stage_cap_seconds": budget["stage_total_active_seconds"],
        "stage_remaining_solver_seconds": max(0.0, budget["stage_total_active_seconds"] - total),
        "repair_active_solver_seconds": used, "repair_cap_seconds": budget["repair_active_solver_cap_seconds"],
        "repair_remaining_solver_seconds": max(0.0, budget["repair_active_solver_cap_seconds"] - used),
        "effective_remaining_solver_seconds": max(0.0, budget["effective_total_active_seconds"] - total),
    }
