"""Unified bounded-stage entry points. Run in the server meep_sim environment."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from .queue import Queue, SolverLedger, aggregate_gate, atomic_json, cgroup_resources, read_json
from .report import collect, generate_report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=["preflight","offline-test","validate","pilot","status","resume","report"])
    parser.add_argument("--manifest")
    parser.add_argument("--detach",action="store_true")
    args=parser.parse_args(argv)
    if args.command == "offline-test":
        # These tests use synthetic arrays only; no Meep solver is invoked.
        return subprocess.call([sys.executable,"-B","-m","unittest","discover","-s","tests","-v"])
    if not args.manifest:
        parser.error("--manifest is required")
    path=Path(args.manifest).resolve()
    if args.command == "preflight":
        queue=Queue(path)
        resources=cgroup_resources()
        print(json.dumps({"status":"READY_FOR_OFFLINE_AND_SHORT_TESTS","resources":resources,
                          "tasks":len(queue.manifest["tasks"]),"fingerprints":queue.locked,
                          "warning":"Preflight does not certify FDTD or authorize extra conditions."},indent=2))
        return 0
    if args.command == "status":
        status=read_json(path.parent/"status.json") if (path.parent/"status.json").exists() else {"state":"NOT_STARTED"}
        brief={k:status.get(k) for k in ("run_id","state","active_task","message","remaining_solver_seconds","usage")}
        brief["task_counts"]={s:sum(t.get("status")==s for t in status.get("tasks",{}).values())
                              for s in sorted({t.get("status") for t in status.get("tasks",{}).values()})}
        print(json.dumps(brief,indent=2))
        return 0
    if args.command == "report":
        result=generate_report(path)
        print(json.dumps({"state":result["state"],"report":str(path.parent/"reports"/"handoff.md")},indent=2))
        return 0
    if args.command == "pilot":
        manifest,_,_,results=collect(path)
        if not aggregate_gate(manifest,results)["passed"]:
            parser.error("Pilot forbidden: qualification gate not passed")
    phase="pilot" if args.command=="pilot" else "all"
    if args.detach:
        log_path=path.parent/"queue.log"
        command=[sys.executable,"-B","-m","ti2d.queue","--manifest",str(path),"--phase",phase]
        with log_path.open("ab") as log:
            child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                   start_new_session=True,env={**os.environ,"PYTHONDONTWRITEBYTECODE":"1"})
        print(json.dumps({"state":"DISPATCHED_NOT_YET_VERIFIED","pid":child.pid,
                          "log":str(log_path),"status":str(path.parent/"status.json")},indent=2))
        return 0
    Queue(path,("pilot",) if phase=="pilot" else ("validation","pilot")).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
