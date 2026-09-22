"""Bounded, restartable serial solver queue. No solver runs on import.

Resume means relaunching an incomplete condition, not resuming a time step.
Only checksum-complete qualified conditions can be reused. Background jobs do
not need an AI conversation or an open SSH session.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

PYTHON = "/root/miniconda3/envs/meep_sim/bin/python"
MPIEXEC = "/root/miniconda3/envs/meep_sim/bin/mpiexec"
TERMINAL_STATUSES = {"QUALIFIED", "NUMERICAL_FAILED", "TIME_LIMIT_UNQUALIFIED",
                     "CRASHED", "MEMORY_LIMIT", "MISSING_OUTPUT"}


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                   separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def ledger_path_for_manifest(manifest, run_dir):
    """A repair run may explicitly share its predecessor's cumulative ledger."""
    default = str(Path(manifest["stage_root"]) / "budget_ledger.json") if manifest.get("stage_root") else "ledger.json"
    path = Path(manifest.get("ledger_path", default))
    return path if path.is_absolute() else Path(run_dir) / path


class FileLock:
    def __init__(self, path):
        self.path = Path(path)
        self.stream = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.stream = self.path.open("a+")
        try:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            self.stream.close()
            raise RuntimeError(f"Already locked: {self.path}") from exc
        return self

    def __exit__(self, *unused):
        fcntl.flock(self.stream.fileno(), fcntl.LOCK_UN)
        self.stream.close()


class SolverLedger:
    """Wall-clock union and rank-hours, including failed/partial attempts."""
    def __init__(self, path):
        self._lock = threading.RLock()
        self.path = Path(path)
        self.data = read_json(path) if self.path.exists() else {"version": 1, "intervals": []}
        if self.data.get("version") != 1 or not isinstance(self.data.get("intervals"), list):
            raise ValueError("Unrecognized or corrupt solver ledger; refusing a fresh budget")
        for record in self.data["intervals"]:
            if record["ranks"] < 1 or not math.isfinite(record["start"]):
                raise ValueError("Invalid solver interval")
            if record.get("end") is not None and record["end"] < record["start"]:
                raise ValueError("Reversed solver interval")

    def add_interval(self, start_epoch, end_epoch, ranks, task_id, status, attempt=1):
        if end_epoch is not None and end_epoch < start_epoch:
            raise ValueError("Negative duration")
        with self._lock:
            self.data["intervals"].append(dict(start=float(start_epoch),
                end=None if end_epoch is None else float(end_epoch), ranks=int(ranks),
                task_id=str(task_id), status=status, attempt=int(attempt)))
            atomic_json(self.path, self.data)
            return len(self.data["intervals"]) - 1

    def close_interval(self, index, end_epoch, status):
        with self._lock:
            record = self.data["intervals"][index]
            record.update(end=max(float(end_epoch), record["start"]), status=status)
            atomic_json(self.path, self.data)

    def intervals(self, now=None):
        now = time.time() if now is None else now
        with self._lock:
            return [(r["start"], r["end"] if r.get("end") is not None else now)
                    for r in self.data["intervals"]]

    def total_active_seconds(self, now=None):
        merged = []
        for start, end in sorted(self.intervals(now)):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        return sum(max(0., b-a) for a, b in merged)

    def rank_hours(self, now=None):
        with self._lock:
            return sum(max(0., end-start)*r["ranks"] for r, (start, end)
                       in zip(self.data["intervals"], self.intervals(now))) / 3600

    def as_dict(self, now=None):
        intervals = self.intervals(now)
        return {"active_solver_seconds": self.total_active_seconds(now),
                "rank_hours": self.rank_hours(now),
                "calendar_span_seconds": max((b for a,b in intervals), default=0)
                    - min((a for a,b in intervals), default=0),
                "attempts": len(intervals)}


def cgroup_resources(root="/sys/fs/cgroup"):
    root = Path(root)
    quota_text = (root / "cpu.max").read_text().split()
    if quota_text[0] == "max":
        # A cpuset is a container limit, unlike host nproc/free output.
        ranges = (root / "cpuset.cpus.effective").read_text().strip().split(",")
        cpus = sum(int(p.split("-")[-1])-int(p.split("-")[0])+1 for p in ranges)
    else:
        cpus = int(quota_text[0]) / int(quota_text[1])
    memory_limit = (root / "memory.max").read_text().strip()
    if memory_limit == "max":
        raise RuntimeError("No finite container memory limit: supply an audited limit before running")
    limit, current = int(memory_limit), int((root / "memory.current").read_text())
    return {"cpu_quota": cpus, "memory_limit_bytes": limit,
            "memory_current_bytes": current, "memory_available_bytes": max(0, limit-current)}


def process_tree_rss(root_pid, proc_root="/proc"):
    records = {}
    page = os.sysconf("SC_PAGE_SIZE")
    for entry in Path(proc_root).glob("[0-9]*"):
        try:
            # comm can contain spaces and ')'; parse after its final closing bracket.
            stat = (entry / "stat").read_text().rsplit(")", 1)[1].split()
            ppid = int(stat[1])
            rss = int((entry / "statm").read_text().split()[1]) * page
            records[int(entry.name)] = (ppid, rss)
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError, IndexError):
            continue
    children, frontier = {root_pid}, [root_pid]
    while frontier:
        parent = frontier.pop()
        new = [pid for pid,(ppid,_) in records.items() if ppid == parent and pid not in children]
        children.update(new)
        frontier.extend(new)
    return sum(records.get(pid, (0,0))[1] for pid in children)


def pair_admission(resources, ranks, memory_fraction=.7):
    group = int(resources["memory_available_bytes"] * memory_fraction)
    if 2*ranks > resources["cpu_quota"]:
        return {"admitted":False,"reason":"Two jobs exceed container CPU quota"}
    if group < 2*1024**3:
        return {"admitted":False,"reason":"Two jobs lack 1 GiB each within group memory cap"}
    return {"admitted":True,"group_memory_budget_bytes":group,
            "memory_budget_bytes":group//2,"resources":resources,"parallel_jobs":2}


def concurrency_assessment(serial, concurrent):
    report={"policy":"return_to_serial","measured_slowdown":None,"compared":False}
    if serial.get("status") != "QUALIFIED" or concurrent.get("status") != "QUALIFIED":
        return {**report,"reason":"Probe or serial qualification failed"}
    if serial.get("reference_reused") is not False or concurrent.get("reference_reused") is not False:
        return {**report,"reason":"Reference cache reuse prevents matched timing comparison"}
    baseline=serial.get("queue_process",{}).get("elapsed_s")
    measured=concurrent.get("queue_process",{}).get("elapsed_s")
    if not all(isinstance(x,(int,float)) and math.isfinite(x) and x>0 for x in (baseline,measured)):
        return {**report,"reason":"Comparable elapsed measurements missing"}
    slowdown=measured/baseline
    return {**report,"compared":True,"measured_slowdown":slowdown,
            "reason":"Slowdown exceeds 30%; serial required" if slowdown>1.3 else
                     "Bounded two-job test passed timing criterion; remaining heterogeneous cases conservatively serial"}


def implementation_fingerprint(package_dir=None):
    if package_dir is None:
        from .common import implementation_digest
        return implementation_digest()
    package_dir = Path(package_dir or Path(__file__).parent)
    return stable_hash({str(p.relative_to(package_dir)): sha256(p)
                        for p in sorted(package_dir.rglob("*.py"))})


def verify_manifest(manifest, run_dir, implementation_sha):
    """Verify creator's hashes before writing the first runtime lock or solver."""
    from .common import digest_object
    expected = digest_object({k:v for k,v in manifest.items() if k != "manifest_sha256"})
    if manifest.get("manifest_sha256") != expected:
        raise ValueError("MANIFEST_PAYLOAD_FINGERPRINT_MISMATCH")
    if manifest.get("implementation_sha256") != implementation_sha:
        raise ValueError("IMPLEMENTATION_CHANGED_SINCE_MANIFEST")
    for task in manifest["tasks"]:
        path = Path(task["case_path"])
        path = path if path.is_absolute() else Path(run_dir)/path
        case = read_json(path)
        expected_case = digest_object({k:v for k,v in case.items() if k != "id"})
        if task.get("case_fingerprint") != expected_case:
            raise ValueError(f"CASE_PAYLOAD_FINGERPRINT_MISMATCH:{task['id']}")
        if case.get("implementation_sha256") != implementation_sha:
            raise ValueError(f"CASE_IMPLEMENTATION_MISMATCH:{task['id']}")


def case_fingerprint(case_path, implementation_sha=None):
    case_path = Path(case_path)
    case = copy.deepcopy(read_json(case_path))
    for key in ("id", "task_id", "phase", "output", "output_dir"):
        case.pop(key, None)
    material = case.get("material_path")
    material_sha = None
    if material:
        material_path = Path(material)
        if not material_path.is_absolute():
            material_path = case_path.parent / material_path
        if not material_path.exists():
            # Production cases normally use paths relative to the stage root.
            material_path = Path(material)
        material_sha = sha256(material_path)
    return stable_hash({"config":case, "material_sha256":material_sha,
                        "implementation_sha256":implementation_sha or implementation_fingerprint()})


def validate_result(result):
    if result.get("status") not in TERMINAL_STATUSES:
        raise ValueError("Solver result has missing/unknown terminal status")
    if result["status"] == "QUALIFIED":
        if result.get("stop", {}).get("converged") is not True:
            raise ValueError("QUALIFIED without converged stop evidence")
        if result.get("stop", {}).get("stop_reason") != "CONVERGED":
            raise ValueError("QUALIFIED with a non-convergence stop reason")
        for metric in ("R", "T", "A_flux", "A_vol", "order_R", "order_T",
                       "reference_pseudo_reflection"):
            value = result.get("metrics", {}).get(metric)
            if not isinstance(value, (int,float)) or isinstance(value,bool) or not math.isfinite(value):
                raise ValueError(f"Invalid qualified metric: {metric}")
        from .solver import result_checks
        checks = result_checks(result.get('case_config', {}), result['metrics'], True,
                               result.get('analytical'), result['stop']['stop_reason'])
        if not all(checks.values()):
            raise ValueError('QUALIFIED does not satisfy the fixed scientific gates')
    return result


def seal_result(task_dir, fingerprint):
    task_dir = Path(task_dir)
    result = validate_result(read_json(task_dir / "result.json"))
    if result["status"] != "QUALIFIED":
        return False
    files = []
    for path in sorted(task_dir.rglob("*")):
        if (not path.is_file() or path.name in {"completion.json", "task.lock", "heartbeat.json"}
                or ".tmp." in path.name):
            continue
        if path.is_symlink():
            raise ValueError("Evidence symlink not allowed")
        files.append({"path":str(path.relative_to(task_dir)), "size":path.stat().st_size,
                      "sha256":sha256(path)})
    atomic_json(task_dir / "completion.json", {"fingerprint":fingerprint,
                 "status":"QUALIFIED", "files":files, "sealed_epoch":time.time()})
    return True


def reusable_result(task_dir, fingerprint):
    task_dir = Path(task_dir)
    try:
        seal = read_json(task_dir / "completion.json")
        if seal.get("fingerprint") != fingerprint or seal.get("status") != "QUALIFIED":
            return None
        if "result.json" not in {item["path"] for item in seal.get("files",[])}:
            return None
        for record in seal["files"]:
            path = task_dir / record["path"]
            if not path.resolve().is_relative_to(task_dir.resolve()) or path.is_symlink():
                return None
            if path.stat().st_size != record["size"] or sha256(path) != record["sha256"]:
                return None
        return validate_result(read_json(task_dir / "result.json"))
    except (OSError, ValueError, KeyError, TypeError):
        return None


def aggregate_gate(manifest, results):
    required = manifest.get("production_gate_task_ids", [])
    failures = []
    if manifest.get('qualification_subset_only'):
        failures.append('Expiry batch is a validation subset, not full production qualification')
    if not required:
        failures.append("No locked production qualification tasks")
    for task_id in required:
        if results.get(task_id, {}).get("status") != "QUALIFIED":
            failures.append(f"{task_id}: not qualified")
    comparisons = []
    for comparison in manifest.get("gate_comparisons", []):
        ids = comparison["ids"]
        values = {}
        passed = all(results.get(task_id,{}).get("status") == "QUALIFIED" for task_id in ids)
        for metric in comparison.get("metrics", ["R", "T", "A_flux"]):
            measured = [results.get(task_id,{}).get("metrics",{}).get(metric) for task_id in ids]
            valid = all(isinstance(v,(int,float)) and math.isfinite(v) for v in measured)
            spread = max(measured)-min(measured) if valid and measured else None
            values[metric] = spread
            passed = passed and spread is not None and spread < comparison["tolerance"]
        comparisons.append({**comparison, "spread":values, "passed":bool(passed)})
        if not passed:
            failures.append(f"{comparison['kind']} {','.join(ids)}: failed/incomplete")
    kinds = {item.get("kind") for item in manifest.get("gate_comparisons", [])}
    for kind in ("mesh", "mirror", "strict"):
        if kind not in kinds:
            failures.append(f"Missing locked {kind} comparison")
    return {"passed":not failures, "failures":failures, "comparisons":comparisons,
            "scope":"specified Ti model, 2D, 10.5 um, +/-30 deg, s/p only"}


def terminate_child(process, grace=20.):
    """Only signal the session created for this queue's own child."""
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=5)
    except ProcessLookupError:
        pass


class Queue:
    def __init__(self, manifest_path, phases=("validation", "pilot")):
        self.manifest_path = Path(manifest_path).resolve()
        self.run_dir = self.manifest_path.parent
        self.manifest = read_json(self.manifest_path)
        self.phases = phases
        self.budget = self.manifest["budgets"]
        self.ledger = SolverLedger(ledger_path_for_manifest(self.manifest,self.run_dir))
        self.results = {}
        self.statuses = {}
        self.cancelled = False
        self.children = {}
        self._state_lock = threading.RLock()
        self.concurrency_probe = {"status":"NOT_RUN","policy":"serial_except_locked_two_job_probe"}
        self.impl_sha = implementation_fingerprint()
        verify_manifest(self.manifest,self.run_dir,self.impl_sha)
        self.manifest_sha = sha256(self.manifest_path)
        self.locked = {"manifest_sha256":self.manifest_sha,
                       "implementation_sha256":self.impl_sha,
                       "case_sha256":{t["id"]:case_fingerprint(self.case_path(t), self.impl_sha)
                                      for t in self.manifest["tasks"]}}
        existing = self.run_dir / "lock_manifest.json"
        if existing.exists() and read_json(existing) != self.locked:
            raise RuntimeError("Locked manifest/implementation changed; preserve this run and use an explicit new run")

    def case_path(self, task):
        path = Path(task["case_path"])
        return path if path.is_absolute() else self.run_dir / path

    def remaining_seconds(self):
        return max(0., min(self.budget['total_active_seconds'] - self.ledger.total_active_seconds(),
                           self.budget['hard_deadline_epoch'] - time.time()))

    def checkpoint(self):
        """Offline report/export only. No AI/API calls or solver launches."""
        from .report import generate_report
        generate_report(self.manifest_path)
        exporter = Path(__file__).resolve().parents[1] / 'export_archive.py'
        subprocess.run([PYTHON, '-B', str(exporter)], check=True, timeout=180)

    def cancel(self, signum, frame):
        self.cancelled = True

    def snapshot(self, state, active=None, message=None):
        with self._state_lock:
            return self._snapshot_locked(state,active,message)

    def _snapshot_locked(self, state, active=None, message=None):
        resources = cgroup_resources()
        summary = {"run_id":self.manifest["run_id"], "state":state,
            "updated_epoch":time.time(), "queue_pid":os.getpid(), "active_task":active,
            "active_tasks":sorted(self.children),
            "message":message, "resources":resources, "usage":self.ledger.as_dict(),
            "remaining_solver_seconds":self.remaining_seconds(),
            "hard_deadline_epoch":self.budget['hard_deadline_epoch'],
            "tasks":self.statuses, "gate":aggregate_gate(self.manifest,self.results),
            "parallel_jobs":len(self.children), "concurrency_probe":self.concurrency_probe,
            "concurrency_reason":"One locked, useful two-job comparison; other conditions stay serial unless controlled evidence justifies wider concurrency",
            "ai_usage":{"actual_usage_available":False, "rmb_cost":None,
                        "note":"The queue makes no AI/API calls. Server billing and this coding turn's AI usage are separate; neither RMB cost is inferred."}}
        atomic_json(self.run_dir / "status.json", summary)
        return summary

    def execute(self, task, attempt, admission=None):
        task_dir = self.run_dir / "tasks" / task["id"] / f"attempt_{attempt}"
        task_dir.mkdir(parents=True, exist_ok=False)
        resources = cgroup_resources()
        ranks = self.budget.get("mpi_ranks",4)
        if ranks > resources["cpu_quota"]:
            raise RuntimeError("Requested MPI ranks exceed measured container CPU quota")
        memory_budget = (admission["memory_budget_bytes"] if admission else
                         int(resources["memory_available_bytes"] * self.budget.get("memory_fraction",.7)))
        if memory_budget < 1024**3:
            raise RuntimeError("Less than 1 GiB admitted task-tree memory: no solver launched")
        remaining = self.remaining_seconds()
        timeout = min(self.budget["single_solver_seconds"], remaining)
        if timeout <= 30:
            return {"status":"TIME_LIMIT_UNQUALIFIED", "reason":"Insufficient total budget to start"}, task_dir
        # Reserve graceful shutdown within both the single-solve and union caps.
        effective_timeout = max(1., timeout - 30.)
        command = [MPIEXEC,"-n",str(ranks),PYTHON,"-B","-m","ti2d.solver",
                   "--case",str(self.case_path(task)),"--output",str(task_dir)]
        env = os.environ.copy()
        env.update({name:"1" for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS",
                   "MKL_NUM_THREADS","NUMEXPR_NUM_THREADS")})
        env.update(PYTHONDONTWRITEBYTECODE="1", TI2D_WALL_TIMEOUT_SECONDS=str(effective_timeout),
                   TI2D_MEMORY_BUDGET_BYTES=str(memory_budget))
        atomic_json(task_dir / "invocation.json", {"command":command, "memory_budget_bytes":memory_budget,
                   "timeout_seconds":timeout,"fingerprint":self.locked["case_sha256"][task["id"]],
                   "resources_before":resources,"attempt":attempt})
        start = time.time()
        cleanup_deadline = time.monotonic() + timeout
        record = self.ledger.add_interval(start,None,ranks,task["id"],"RUNNING",attempt)
        reason, peak, solver_end_epoch, final_status, process = None, 0, None, "CRASHED", None
        with FileLock(task_dir / "task.lock"), (task_dir / "stdout.log").open("wb") as log:
            try:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                           env=env, start_new_session=True)
                with self._state_lock:
                    self.children[task["id"]]=process.pid
                atomic_json(task_dir / "process.json", {"pid":process.pid,"start_epoch":start})
                while process.poll() is None:
                    elapsed = time.time()-start
                    peak = max(peak, process_tree_rss(process.pid))
                    if self.cancelled:
                        reason = "INTERRUPTED_NO_FIELD_RESTART"
                    elif elapsed >= effective_timeout:
                        reason = "TIME_LIMIT_UNQUALIFIED"
                    elif peak > memory_budget:
                        reason = "MEMORY_LIMIT"
                    if reason:
                        break
                    self.snapshot("RUNNING",task["id"])
                    time.sleep(min(5., max(.1,effective_timeout-elapsed)))
                from .session_guard import stop_owned_session
                cleanup = stop_owned_session(process, cleanup_deadline,
                    lambda: self.snapshot('STOPPING',task['id']))
                atomic_json(task_dir / 'session_cleanup.json', cleanup)
                if not cleanup['stopped']:
                    raise RuntimeError('OWNED_SESSION_STILL_LIVE_KEEP_LEDGER_OPEN')
                exit_code = process.returncode
                solver_end_epoch = time.time()
                if reason == "INTERRUPTED_NO_FIELD_RESTART":
                    status = "CRASHED"
                elif reason:
                    status = reason
                else:
                    status = None
                result_path = task_dir / "result.json"
                try:
                    result = validate_result(read_json(result_path))
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    result = {"status":"MISSING_OUTPUT","reason":str(exc)}
                if not reason and exit_code and result.get("status") in ("QUALIFIED","MISSING_OUTPUT"):
                    status = "CRASHED"
                if status:
                    result["solver_reported_status"] = result.get("status")
                    result.update(status=status, reason=reason or f"process exit code {exit_code}")
                final_status = result["status"]
                result.update(queue_process={"exit_code":exit_code,"elapsed_s":time.time()-start,
                    "peak_process_tree_rss_bytes":peak,"mpi_ranks":ranks,
                    "memory_budget_bytes":memory_budget,"attempt":attempt})
                atomic_json(result_path,result)
                if result["status"] == "QUALIFIED":
                    seal_result(task_dir,self.locked["case_sha256"][task["id"]])
                return result, task_dir
            finally:
                from .session_guard import stop_owned_session
                cleanup = (stop_owned_session(process, cleanup_deadline, lambda: None)
                           if process else {'stopped': True})
                with self._state_lock:
                    self.children.pop(task["id"],None)
                if cleanup['stopped']:
                    self.ledger.close_interval(record,time.time(),final_status)
                else:
                    raise RuntimeError('OWNED_SESSION_STILL_LIVE_KEEP_LEDGER_OPEN')

    def eligible_pair(self, task):
        probe=self.manifest.get("concurrency_probe")
        if not probe or self.concurrency_probe["status"] != "NOT_RUN":
            return None
        ids=probe.get("task_ids",[])
        if len(ids)!=2 or task["id"]!=ids[0]:
            return None
        by_id={t["id"]:t for t in self.manifest["tasks"]}
        if any(i not in by_id for i in ids):
            raise ValueError("Locked concurrency probe refers to unknown task")
        serial=self.results.get(probe.get("serial_task_id"),{})
        if serial.get("status")!="QUALIFIED":
            self.concurrency_probe={"status":"SKIPPED","reason":"Serial timing baseline unavailable"}
            return None
        paired=[by_id[i] for i in ids]
        serial_task=by_id.get(probe.get("serial_task_id"))
        serial_case=read_json(self.case_path(serial_task)) if serial_task else {}
        probe_case=read_json(self.case_path(paired[0]))
        stripped=lambda case:{k:v for k,v in case.items() if k not in ("id","emission_theta_deg")}
        if (not serial_case or stripped(serial_case)!=stripped(probe_case)
                or serial_case.get("case")!="flat"
                or serial_case.get("emission_theta_deg")!=-probe_case.get("emission_theta_deg",0)):
            self.concurrency_probe={"status":"SKIPPED","reason":"Locked timing cases are not opposite-angle same-layout flat controls"}
            return None
        if any(t["phase"]!="validation" or t.get("reused_from") or
               any(self.results.get(dep,{}).get("status")!="QUALIFIED" for dep in t.get("depends_on",[]))
               for t in paired):
            raise ValueError("Concurrency probe must contain two unblocked, required validation solves")
        if any((self.run_dir/"tasks"/t["id"]).exists() for t in paired):
            self.concurrency_probe={"status":"SKIPPED","reason":"Prior probe attempt exists; no duplicate launch"}
            return None
        admission=pair_admission(cgroup_resources(),self.budget.get("mpi_ranks",4),self.budget.get("memory_fraction",.7))
        if not admission["admitted"]:
            self.concurrency_probe={"status":"SKIPPED",**admission}
            return None
        return paired,serial,admission

    def execute_pair(self, paired, serial, admission):
        self.concurrency_probe={"status":"RUNNING","task_ids":[t["id"] for t in paired],**admission}
        outputs={}
        with ThreadPoolExecutor(max_workers=2,thread_name_prefix="ti2d-needed-probe") as pool:
            futures={pool.submit(self.execute,t,1,admission):t for t in paired}
            for future in as_completed(futures):
                task=futures[future]
                try:
                    result,task_dir=future.result()
                except Exception:
                    self.cancelled=True
                    raise
                outputs[task["id"]]=(result,task_dir)
                if result["status"]!="QUALIFIED":
                    # Stop this queue's peer condition, never unrelated user processes.
                    self.cancelled=True
        assessment=concurrency_assessment(serial,outputs[paired[0]["id"]][0])
        self.concurrency_probe={**self.concurrency_probe,**assessment,"status":"COMPLETE"}
        atomic_json(self.run_dir/"concurrency_probe.json",self.concurrency_probe)
        return outputs

    def completed_comparison_failed(self):
        gate=aggregate_gate(self.manifest,self.results)
        return any(not comparison["passed"] and
                   all(self.results.get(task_id,{}).get("status")=="QUALIFIED"
                       for task_id in comparison["ids"])
                   for comparison in gate["comparisons"])

    def run(self):
        with FileLock(self.run_dir / "queue.lock"), FileLock(self.ledger.path.with_suffix(".lock")):
            # Another run could have changed this shared ledger after construction.
            self.ledger = SolverLedger(self.ledger.path)
            atomic_json(self.run_dir / "lock_manifest.json",self.locked)
            if any(r.get("end") is None for r in self.ledger.data["intervals"]):
                raise RuntimeError("Unclosed prior solver interval: inspect surviving process and reconcile ledger before resume; no duplicate launch")
            signal.signal(signal.SIGTERM,self.cancel)
            signal.signal(signal.SIGINT,self.cancel)
            for task in self.manifest["tasks"]:
                self.statuses[task["id"]] = {"status":"NOT_RUN"}
            failed_designs = set()
            final = "COMPLETE"
            for task in self.manifest["tasks"]:
                task_id = task["id"]
                if task_id in self.results:
                    continue
                fingerprint = self.locked["case_sha256"][task_id]
                task_root = self.run_dir / "tasks" / task_id
                attempts = sorted(task_root.glob("attempt_*")) if task_root.exists() else []
                result = None
                for existing in reversed(attempts):
                    result = reusable_result(existing,fingerprint)
                    if result:
                        self.results[task_id] = result
                        self.statuses[task_id] = {"status":"QUALIFIED","reused":True,"result_path":str(existing / "result.json")}
                        break
                if result:
                    continue
                if task.get("reused_from"):
                    source = task["reused_from"]
                    if (self.locked["case_sha256"].get(source) == fingerprint
                            and self.results.get(source,{}).get("status") == "QUALIFIED"):
                        self.results[task_id] = self.results[source]
                        self.statuses[task_id] = {"status":"QUALIFIED","reused_from":source,
                             "result_path":self.statuses[source]["result_path"]}
                        continue
                    self.statuses[task_id] = {"status":"BLOCKED","reason":"Exact qualified reuse source unavailable"}
                    final = "BLOCKED"
                    break
                if task["phase"] not in self.phases:
                    continue
                if self.cancelled:
                    final = "INTERRUPTED"
                    break
                if self.remaining_seconds() <= 30:
                    final = "BUDGET_EXHAUSTED_PARTIAL"
                    break
                if task["phase"] == "pilot" and not aggregate_gate(self.manifest,self.results)["passed"]:
                    final = "VALIDATION_FAILED_OR_INCOMPLETE"
                    break
                if any(self.results.get(t,{}).get("status") != "QUALIFIED" for t in task.get("depends_on",[])):
                    self.statuses[task_id] = {"status":"BLOCKED","reason":"Dependency not qualified"}
                    final = "PARTIAL"
                    continue
                if task.get("geometry_id") in failed_designs:
                    self.statuses[task_id] = {"status":"ISOLATED_DESIGN","reason":"Another excitation of this design failed"}
                    continue
                if attempts:
                    # No automatic numerical/crash retries, including across process restarts.
                    self.statuses[task_id] = {"status":"NEEDS_REVIEW","reason":"Prior incomplete/unqualified attempt; no field checkpoint continuation or automatic retry"}
                    final = "NEEDS_REVIEW"
                    break
                pair=self.eligible_pair(task)
                if pair is not None:
                    outputs=self.execute_pair(*pair)
                    with self._state_lock:
                        for paired_id,(paired_result,paired_dir) in outputs.items():
                            self.results[paired_id]=paired_result
                            self.statuses[paired_id]={"status":paired_result["status"],
                                "result_path":str(paired_dir/"result.json")}
                    if (any(r["status"]!="QUALIFIED" for r,d in outputs.values())
                            or self.completed_comparison_failed()):
                        final="VALIDATION_FAILED"
                        break
                    self.snapshot("RUNNING")
                    continue
                result, task_dir = self.execute(task,1)
                self.results[task_id] = result
                self.statuses[task_id] = {"status":result["status"],"result_path":str(task_dir / "result.json")}
                self.snapshot('CHECKPOINT')
                self.checkpoint()
                if result["status"] != "QUALIFIED":
                    if task["phase"] == "validation":
                        final = "VALIDATION_FAILED"
                        break
                    failed_designs.add(task.get("geometry_id"))
                    final = "PARTIAL"
                if task["phase"]=="validation" and self.completed_comparison_failed():
                    final="VALIDATION_FAILED"
                    break
                self.snapshot("RUNNING")
            if final=="COMPLETE" and not aggregate_gate(self.manifest,self.results)["passed"]:
                final="VALIDATION_FAILED_OR_INCOMPLETE"
            summary = self.snapshot(final)
            self.checkpoint()
            return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--phase",choices=["all","validation","pilot"],default="all")
    args = parser.parse_args(argv)
    phases = ("validation","pilot") if args.phase == "all" else (args.phase,)
    try:
        Queue(args.manifest,phases).run()
    except Exception as exc:
        run_dir = Path(args.manifest).resolve().parent
        atomic_json(run_dir / "queue_error.json",
                    {"status":"QUEUE_ERROR","error":str(exc),"epoch":time.time()})
        if "Already locked:" not in str(exc):
            status_path = run_dir / "status.json"
            try:
                status = read_json(status_path) if status_path.exists() else {}
            except (OSError,ValueError):
                status = {}
            status.update(state="QUEUE_ERROR",active_task=None,message=str(exc),updated_epoch=time.time())
            atomic_json(status_path,status)
            try:
                from .report import generate_report
                generate_report(args.manifest)
            except Exception:
                pass
        raise


if __name__ == "__main__":
    main()
