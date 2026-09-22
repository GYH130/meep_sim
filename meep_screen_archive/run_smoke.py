"""Server wrapper: bounded MPI smoke, charged to the shared stage budget."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from ti2d.queue import FileLock, SolverLedger, process_tree_rss, cgroup_resources, terminate_child


stage = Path(__file__).resolve().parent
if str(stage) != '/root/autodl-tmp/meep_sim/meep_screen' or sys.platform != 'linux':
    raise RuntimeError('Server only')
attempt = stage / 'smoke' / time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())
attempt.mkdir(parents=True, exist_ok=False)
ledger = SolverLedger(stage / 'budget_ledger.json')
with FileLock(stage / 'budget_ledger.lock'):
    resources = cgroup_resources()
    if resources['cpu_quota'] < 4 or ledger.total_active_seconds() + 180 > 172800:
        raise RuntimeError('Insufficient resource or solver budget')
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
           'MKL_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1'}
    command = ['/root/miniconda3/envs/meep_sim/bin/mpiexec', '-n', '4', sys.executable,
               '-B', str(stage / 'solver_smoke.py'), '--output', str(attempt)]
    start = time.time()
    index = ledger.add_interval(start, None, 4, 'api_smoke', 'RUNNING')
    status = 'CRASHED'
    with (attempt / 'stdout.log').open('wb') as log:
        child = subprocess.Popen(command, cwd=stage, env=env, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        try:
            while child.poll() is None:
                if time.time() - start > 150:
                    status = 'TIME_LIMIT_UNQUALIFIED'
                    terminate_child(child, grace=20)
                    break
                if process_tree_rss(child.pid) > .7 * resources['memory_available_bytes']:
                    status = 'MEMORY_LIMIT'
                    terminate_child(child, grace=20)
                    break
                time.sleep(.5)
            if child.returncode == 0:
                status = 'SMOKE_PASSED_NOT_PHYSICALLY_QUALIFIED'
        finally:
            if child.poll() is None:
                terminate_child(child, grace=20)
            ledger.close_interval(index, time.time(), status)
    print(json.dumps({'status': status, 'path': str(attempt), 'solver_usage': ledger.as_dict()}))
    raise SystemExit(0 if child.returncode == 0 else 1)
