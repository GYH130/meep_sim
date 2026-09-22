"""Server-only bounded diagnostic, charged to the frozen shared repair budget."""
from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

from repair_budget import locked_repair_budget, repair_usage
from ti2d.queue import (FileLock, MPIEXEC, PYTHON, SolverLedger, atomic_json,
                        cgroup_resources, process_tree_rss, sha256, terminate_child)


def _assert_server(stage):
    if sys.platform != "linux" or str(stage) != "/root/autodl-tmp/meep_sim/meep_screen":
        raise RuntimeError("SERVER_ONLY")


def run_probe(stage):
    stage = Path(stage).resolve()
    _assert_server(stage)
    if not (stage / "budget_ledger.json").is_file():
        raise RuntimeError("EXISTING_SHARED_LEDGER_REQUIRED")
    with FileLock(stage / "budget_ledger.lock"):
        ledger = SolverLedger(stage / "budget_ledger.json")
        budget = locked_repair_budget(stage, ledger)
        resources = cgroup_resources()
        memory_budget = int(.7 * resources["memory_available_bytes"])
        allowance = min(180.0, repair_usage(ledger, budget)["effective_remaining_solver_seconds"])
        if resources["cpu_quota"] < 4 or memory_budget < 1024**3 or allowance <= 30:
            raise RuntimeError("INSUFFICIENT_REPAIR_BUDGET_OR_RESOURCES")
        out = stage / "repair_diagnostics" / (time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + "_" + uuid.uuid4().hex[:8])
        out.mkdir(parents=True, exist_ok=False)
        script = stage / "repair_flux_probe.py"
        command = [MPIEXEC, "-n", "4", PYTHON, "-B", str(script), "--output", str(out)]
        env = {**os.environ, **{name: "1" for name in
               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
               "PYTHONDONTWRITEBYTECODE": "1"}
        atomic_json(out / "invocation.json", {
            "command": command, "script_sha256": sha256(script),
            "wrapper_sha256": sha256(__file__), "repair_budget_sha256": budget["budget_sha256"],
            "resources_before": resources, "memory_budget_bytes": memory_budget,
            "timeout_seconds_including_shutdown": allowance, "shutdown_reserve_seconds": 30,
            "usage_before": repair_usage(ledger, budget),
        })
        cancelled = [False]
        def cancel(*unused):
            cancelled[0] = True
        prior_handlers = {sig: signal.signal(sig, cancel) for sig in (signal.SIGTERM, signal.SIGINT)}
        start, monotonic_start = time.time(), time.monotonic()
        index, child = None, None
        status, error, peak = "CRASHED", None, 0
        try:
            index = ledger.add_interval(start, None, 4, "repair1_nonzero_flux_probe", "RUNNING")
            with (out / "stdout.log").open("wb") as log:
                if cancelled[0]:
                    status = "INTERRUPTED_NO_FIELD_RESTART"
                else:
                    child = subprocess.Popen(command, cwd=stage, env=env, stdout=log,
                                             stderr=subprocess.STDOUT, start_new_session=True)
                    atomic_json(out / "process_start.json", {"pid": child.pid, "start_epoch": start})
                    while child.poll() is None:
                        peak = max(peak, process_tree_rss(child.pid))
                        elapsed = time.monotonic() - monotonic_start
                        if cancelled[0]:
                            status = "INTERRUPTED_NO_FIELD_RESTART"
                        elif elapsed >= allowance - 30:
                            status = "TIME_LIMIT_UNQUALIFIED"
                        elif peak > memory_budget:
                            status = "MEMORY_LIMIT"
                        else:
                            time.sleep(min(.5, max(.01, allowance - 30 - elapsed)))
                            continue
                        terminate_child(child, grace=20)
                        break
                    if status == "CRASHED" and child.returncode == 0:
                        status = "DIAGNOSTIC_COMPLETED_NOT_QUALIFICATION"
        except Exception as exc:
            error = repr(exc)
        finally:
            try:
                if child is not None and child.poll() is None:
                    terminate_child(child, grace=20)
            finally:
                if index is not None:
                    ledger.close_interval(index, time.time(), status)
                for sig, handler in prior_handlers.items():
                    signal.signal(sig, handler)
                result = {"status": status, "exit_code": child.returncode if child is not None else None,
                          "error": error, "path": str(out), "elapsed_s": time.monotonic() - monotonic_start,
                          "peak_process_tree_rss_bytes": peak, "mpi_ranks": 4,
                          "memory_budget_bytes": memory_budget, "solver_usage": ledger.as_dict(),
                          "repair_budget": repair_usage(ledger, budget)}
                atomic_json(out / "process.json", result)
                print(json.dumps(result), flush=True)
        return 0 if status == "DIAGNOSTIC_COMPLETED_NOT_QUALIFICATION" else 1


def main():
    return run_probe(Path(__file__).resolve().parent)


if __name__ == "__main__":
    raise SystemExit(main())
