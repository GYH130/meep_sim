"""J: new bounded matched-pair diagnostic; preserve I and all qualification gates."""
from pathlib import Path
from copy import deepcopy
import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.queue import FileLock, SolverLedger, cgroup_resources, process_tree_rss, PYTHON, MPIEXEC
from ti2d.session_guard import session_live_members, stop_owned_session
from ti2d.solver import layout
from recover_tail_analysis import recover

ROOT=Path(__file__).resolve().parent
SERVER=Path('/root/autodl-tmp/meep_sim/meep_screen/expiry_campaign')
RUN=ROOT/'runs/tail_20260923_J'
SCRIPTS=('tail_probe_v2.py','tail_analysis.py','recover_tail_analysis.py','run_tail_j.py',
         'test_recover_tail_analysis.py','test_tail_j.py')
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')


def allowance(plan,ledger,now,single_cap):
    return max(0.,min(single_cap,21600.,172800.-ledger.total_active_seconds(now),
                     plan['ledger_active_before']+32400.-ledger.total_active_seconds(now)))


def prepare():
    with FileLock(ROOT.parent/'budget_ledger.lock'):
        if RUN.exists(): raise FileExistsError('SINGLE_USE_RUN_ALREADY_EXISTS')
        ledger=SolverLedger(ROOT.parent/'budget_ledger.json')
        if any(r.get('end') is None for r in ledger.data['intervals']):
            raise RuntimeError('PRIOR_SOLVER_INTERVAL_STILL_OPEN')
        h=ROOT/'runs/expiry_20260922_H'
        c=json.loads((h/'cases/slot32p.json').read_text())
        old=h/'tasks/slot32p/attempt_1'
        if implementation_digest()!=c['implementation_sha256']:
            raise RuntimeError('H_CORE_FINGERPRINT_CHANGED')
        if json.loads((old/'result.json').read_text())['status']!='TIME_LIMIT_UNQUALIFIED':
            raise RuntimeError('UNEXPECTED_H_STATE')
        if digest_file(c['material_path'])!=c['material_sha256']:
            raise RuntimeError('MATERIAL_CHANGED')
        override=json.loads((ROOT/'RENTAL_DEADLINE_OVERRIDE.json').read_text())
        if override['hard_deadline_epoch'] is not None:
            raise RuntimeError('RENTAL_POLICY_CHANGED_REQUIRES_NEW_PLAN')
        if 172800.-ledger.total_active_seconds()<32400.:
            raise RuntimeError('INSUFFICIENT_TOTAL_SOLVER_BUDGET')
        previous=ROOT/'runs/tail_20260923_I'
        recovered=recover(previous)
        if recovered['arms']['pml6_control']['status']!='TIME_LIMIT_UNQUALIFIED':
            raise RuntimeError('UNEXPECTED_I_SOLVER_STATUS')
        (RUN/'inputs').mkdir(parents=True)
        source={'H_case.json':h/'cases/slot32p.json','H_result.json':old/'result.json',
                'H_stop.json':old/'structure_stop.json','ledger_before.json':ledger.path,
                'rental_override.json':ROOT/'RENTAL_DEADLINE_OVERRIDE.json',
                'I_recovered_summary.json':previous/'reports_recovered/diagnostic_summary.json',
                'I_process.json':previous/'pml6_control/process.json'}
        source.update({p.name:p for p in Path(c['material_path']).parent.iterdir() if p.is_file()})
        for n,p in source.items(): shutil.copyfile(p,RUN/'inputs'/n)
        shutil.copytree(ROOT/'ti2d',RUN/'inputs/frozen_core',ignore=shutil.ignore_patterns('__pycache__'))
        repo=ROOT.parents[1]
        status=subprocess.run(['git','-C',str(repo),'status','--porcelain'],capture_output=True,check=True).stdout
        (RUN/'inputs/git_status_before.txt').write_bytes(status)
        diff=subprocess.run(['git','-C',str(repo),'diff','--binary','HEAD'],capture_output=True,check=True).stdout
        (RUN/'inputs/protected_changes.patch').write_bytes(diff)
        atomic_json(RUN/'inputs/directory_presence.json',{str(p):'present' if p.exists() else 'not found'
                    for p in (repo/'sources',repo/'review',repo/'src')})
        cases=[]
        for arm,pml,cap in (('pml6_control',6.,14400.),('pml12_test',12.,18000.)):
            cc=deepcopy(c); cc.update(id=arm,pml_um=pml,diagnostic_stop_time=1200.,max_solver_seconds=cap)
            cc['diagnostic_only']=True
            if any(layout(cc)[k]!=layout(c)[k] for k in ('source','refl','trans','surface','bottom','top','non_pml_bottom','non_pml_height')):
                raise RuntimeError('PHYSICAL_MONITOR_OR_STRUCTURE_POSITION_CHANGED')
            path=RUN/'cases'/(arm+'.json'); atomic_json(path,cc)
            cases.append(dict(id=arm,path=str(path),sha256=digest_file(path),single_cap=cap))
        plan=dict(run_id=RUN.name,created_utc=utcnow(),tasks=cases,qualification=False,
            source_case_sha256=digest_file(source['H_case.json']),implementation_sha256=implementation_digest(),
            inputs={str(p.relative_to(RUN)):digest_file(p) for p in (RUN/'inputs').rglob('*') if p.is_file()},
            scripts={n:digest_file(ROOT/n) for n in SCRIPTS},ledger_active_before=ledger.total_active_seconds(),
            budgets=dict(phase_solver_seconds=32400.,single_solver_seconds={'pml6_control':14400.,'pml12_test':18000.},
                         total_solver_seconds=172800.,hard_deadline_epoch=None,mpi_ranks=4,max_parallel=1),
            objective='Replay H tail with complex waveforms; vary only PML thickness 6 to 12 um; identify spatial/spectral and boundary sensitivity',
            restart='Fresh zero-field runs; no timestep resume; I remains timed-out and counted',
            measured_cost_basis=dict(I_seconds=8955.034698180854,I_simulation_time=914.625,
                 target_simulation_time=1200.,control_expected_hours=3.27,test_estimated_hours=3.68,
                 estimate_not_guarantee=True),
            production_retries=0,scope='diagnostic only; no qualification, cache, production sweep or automatic retry')
        plan['sha256']=digest_object(plan); atomic_json(RUN/'plan.json',plan)
        atomic_json(RUN/'status.json',dict(state='PREPARED',qualification=False,plan_sha256=plan['sha256']))
        return plan


def package():
    exports=ROOT/'exports'; exports.mkdir(exist_ok=True)
    temporary=exports/'J_raw_evidence.partial.tar.gz'
    with tarfile.open(temporary,'w:gz') as tf:
        tf.add(RUN,arcname=RUN.name)
        for n in SCRIPTS: tf.add(ROOT/n,arcname='code/'+n)
    path=exports/'J_raw_evidence.tar.gz'; os.replace(temporary,path)
    atomic_json(exports/'J_export.json',dict(file=path.name,sha256=digest_file(path),bytes=path.stat().st_size,
                                           completed_utc=utcnow(),qualification=False))


def run():
    plan=json.loads((RUN/'plan.json').read_text())
    if digest_object({k:v for k,v in plan.items() if k!='sha256'})!=plan['sha256']:
        raise RuntimeError('PLAN_CHANGED')
    for name,expected in plan['scripts'].items():
        if digest_file(ROOT/name)!=expected: raise RuntimeError('DIAGNOSTIC_CODE_CHANGED')
    for name,expected in plan['inputs'].items():
        if digest_file(RUN/name)!=expected: raise RuntimeError('FROZEN_INPUT_CHANGED')
    if implementation_digest()!=plan['implementation_sha256']: raise RuntimeError('CORE_CHANGED')
    cancelled=[False]
    for sig in (signal.SIGTERM,signal.SIGINT):
        signal.signal(sig,lambda *_:cancelled.__setitem__(0,True))
    records=[]
    with FileLock(RUN/'run.lock'),FileLock(ROOT.parent/'budget_ledger.lock'):
        ledger=SolverLedger(ROOT.parent/'budget_ledger.json')
        if any(r.get('end') is None for r in ledger.data['intervals']): raise RuntimeError('OPEN_SOLVER_INTERVAL')
        try:
            for task in plan['tasks']:
                if cancelled[0]: break
                limit=allowance(plan,ledger,time.time(),task['single_cap'])
                if limit<300: break
                out=RUN/task['id']; out.mkdir(exist_ok=False)
                if digest_file(task['path'])!=task['sha256']: raise RuntimeError('CASE_CHANGED')
                resources=cgroup_resources(); memory=int(.7*resources['memory_available_bytes'])
                if resources['cpu_quota']<4 or memory<2*1024**3: raise RuntimeError('RESOURCES_UNAVAILABLE')
                env=dict(os.environ,**{n:'1' for n in THREADS},PYTHONDONTWRITEBYTECODE='1',TAIL_WALL_SECONDS=str(limit-45))
                command=[MPIEXEC,'-n','4',PYTHON,'-B',str(ROOT/'tail_probe_v2.py'),'--case',task['path'],'--output',str(out)]
                atomic_json(out/'invocation.json',dict(command=command,allowance_s=limit,resources=resources,
                    memory_budget_bytes=memory,threads_per_rank=1))
                child=None; began=time.monotonic(); peak=0; state='CRASHED'; cleanup=None
                index=ledger.add_interval(time.time(),None,4,RUN.name+':'+task['id'],'RUNNING')
                def publish():
                    atomic_json(RUN/'status.json',dict(state='RUNNING' if state=='CRASHED' else state,
                        active_task=task['id'],supervisor_pid=os.getpid(),pid=child.pid if child else None,
                        updated_utc=utcnow(),elapsed_s=time.monotonic()-began,peak_tree_rss_bytes=peak,
                        memory_budget_bytes=memory,allowance_s=limit,ledger=ledger.as_dict(),qualification=False))
                try:
                    with (out/'stdout.log').open('xb') as log:
                        child=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                        while child.poll() is None:
                            peak=max(peak,process_tree_rss(child.pid)); publish()
                            if cancelled[0]: state='INTERRUPTED'; break
                            if peak>memory: state='MEMORY_LIMIT'; break
                            if time.monotonic()-began>=limit-30: state='TIME_LIMIT_UNQUALIFIED'; break
                            time.sleep(5)
                    if child.poll() is not None and (out/'result.json').exists():
                        candidate=json.loads((out/'result.json').read_text())['status']
                        state=candidate if child.returncode==0 or candidate!='DIAGNOSTIC_COMPLETE' else 'CRASHED'
                finally:
                    if child is None: cleanup={'stopped':True}
                    else: cleanup=stop_owned_session(child,min(began+limit,time.monotonic()+30),publish)
                    atomic_json(out/'process.json',dict(status=state,cleanup=cleanup,exit_code=child.poll() if child else None,
                        elapsed_s=time.monotonic()-began,peak_tree_rss_bytes=peak))
                    if cleanup['stopped']: ledger.close_interval(index,time.time(),state)
                    else: raise RuntimeError('OWNED_MPI_SESSION_NOT_CONFIRMED_STOPPED')
                records.append(dict(id=task['id'],status=state))
                # Store actual solver state before independent report generation.
                atomic_json(RUN/'status.json',dict(state='ANALYZING',tasks=records,qualification=False,updated_utc=utcnow()))
                try:
                    analysis=recover(RUN)
                except Exception as exc:
                    atomic_json(RUN/'analysis_error.json',dict(error=str(exc),qualification=False,
                                solver_status=state,updated_utc=utcnow()))
                    raise RuntimeError('POSTPROCESSING_FAILED:'+str(exc)) from exc
                if state!='DIAGNOSTIC_COMPLETE': break
                if task['id']=='pml6_control' and not analysis['control_reproduced']: break
            atomic_json(RUN/'status.json',dict(state='DIAGNOSTIC_COMPLETE' if len(records)==2 and all(r['status']=='DIAGNOSTIC_COMPLETE' for r in records) else 'DIAGNOSTIC_INCOMPLETE',
                active_task=None,tasks=records,qualification=False,ledger=ledger.as_dict(),updated_utc=utcnow()))
        except Exception as exc:
            atomic_json(RUN/'status.json',dict(state='FAILED',error=str(exc),tasks=records,qualification=False,updated_utc=utcnow()))
            raise
        finally:
            try:
                recover(RUN)
            except Exception as exc:
                atomic_json(RUN/'analysis_error.json',dict(error=str(exc),qualification=False))
            finally:
                package()


if __name__=='__main__':
    if sys.platform!='linux' or ROOT!=SERVER: raise RuntimeError('SERVER_ONLY')
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['prepare','run']); a=p.parse_args()
    if a.action=='prepare': print(json.dumps(prepare(),indent=2))
    else: run()
