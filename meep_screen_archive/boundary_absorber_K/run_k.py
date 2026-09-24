"""K: two-case, dependency-gated absorber verification; server only, no AI polling."""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.geometry import validate_geometry
from ti2d.materials import load_model, audit_model
from ti2d.queue import FileLock, Queue, SolverLedger, verify_manifest

ROOT = Path(__file__).resolve().parent
SERVER = Path('/root/autodl-tmp/meep_sim/meep_screen/boundary_absorber_K')
RUN = ROOT / 'runs/boundary_20260924_K'
SOLVER_DEADLINE = datetime(2026, 9, 24, 12, 30, tzinfo=timezone.utc).timestamp()
RENTAL_EXPIRY = datetime(2026, 9, 24, 13, 53, tzinfo=timezone.utc).timestamp()
H_CORE = '67f6d37394b51edaa5562fe27058723da3eae1e642d4d43ff4fe956aea6a568e'
MATERIAL = '4791b3868810fb4a3fa93c8469b488650270c0f47a18804c9b26bd7354153b5e'


def cap_seconds(task, remaining):
    return max(0., min(float(task['single_cap_seconds']), 21600., remaining))


def prepare():
    with FileLock(ROOT.parent / 'budget_ledger.lock'):
        if RUN.exists():
            raise FileExistsError('SINGLE_USE_RUN_ALREADY_EXISTS')
        ledger = SolverLedger(ROOT.parent / 'budget_ledger.json')
        if any(r.get('end') is None for r in ledger.data['intervals']):
            raise RuntimeError('OPEN_PRIOR_SOLVER_INTERVAL')
        if min(SOLVER_DEADLINE-time.time(), 172800-ledger.total_active_seconds()) < 3600:
            raise RuntimeError('INSUFFICIENT_WINDOW_NO_SOLVER_LAUNCHED')
        history = ROOT.parent / 'expiry_campaign/runs/expiry_20260922_H'
        base = json.loads((history / 'cases/slot32p.json').read_text())
        if base['implementation_sha256'] != H_CORE or base['material_sha256'] != MATERIAL:
            raise RuntimeError('FROZEN_H_INPUT_CHANGED')
        if digest_file(base['material_path']) != MATERIAL:
            raise RuntimeError('MATERIAL_CHANGED')
        (RUN / 'inputs').mkdir(parents=True)
        shutil.copytree(Path(base['material_path']).parent, ROOT / 'assets/material')
        base['material_path'] = str(ROOT / 'assets/material/route_a_material.json')
        for name in ('flat32p', 'slot32p'):
            shutil.copyfile(history / 'cases' / (name+'.json'), RUN / 'inputs' / ('H_'+name+'_case.json'))
            for filename in ('result.json', 'structure_stop.json'):
                source = history / 'tasks' / name / 'attempt_1' / filename
                if source.exists():
                    shutil.copyfile(source, RUN / 'inputs' / ('H_'+name+'_'+filename))
        shutil.copyfile(ledger.path, RUN / 'inputs/ledger_before.json')
        repo = ROOT.parents[1]
        (RUN / 'inputs/git_status_before.txt').write_bytes(subprocess.run(
            ['git','-C',str(repo),'status','--porcelain'], check=True, capture_output=True).stdout)
        (RUN / 'inputs/protected_changes.patch').write_bytes(subprocess.run(
            ['git','-C',str(repo),'diff','--binary','HEAD'], check=True, capture_output=True).stdout)
        tracked = subprocess.run(['git','-C',str(repo),'diff','--name-only','HEAD'],
                                 check=True, capture_output=True, text=True).stdout.splitlines()
        atomic_json(RUN / 'inputs/protected_fingerprints.json', {
            name: digest_file(repo/name) for name in tracked if (repo/name).is_file()})
        atomic_json(RUN / 'inputs/directory_presence.json', {
            str(p): 'present' if p.exists() else 'not found' for p in
            (repo/'sources', repo/'review', repo/'src')})
        implementation = implementation_digest()
        model = load_model(base['material_path'])
        tasks, differences = [], {}
        for name, kind, cap in (('K_flat32p','flat',14400.), ('K_slot32p','structure',21600.)):
            case = deepcopy(base)
            case.update(id=name, case=kind, pml_um=12., boundary_kind='absorber',
                        boundary_profile={'name':'quadratic','R_asymptotic':1e-15,'mean_stretch':1.0},
                        implementation_sha256=implementation,
                        max_solver_seconds=cap-120.)
            if not validate_geometry(case['geometry'])['valid']:
                raise RuntimeError('INVALID_GEOMETRY')
            audit_model(model, case['resolution'], case['courant'], case['wavelength_um'])
            path = RUN / 'cases' / (name+'.json')
            atomic_json(path, case)
            tasks.append(dict(id=name, phase='validation', case_path=str(path),
                depends_on=[tasks[-1]['id']] if tasks else [], single_cap_seconds=cap,
                geometry_id='baseline', gate_role='flat' if kind=='flat' else 'baseline',
                resolution=32, polarization='p', emission_theta_deg=30, tilt_deg=30,
                strict_stop=False, case_fingerprint=digest_object({k:v for k,v in case.items() if k!='id'})))
            differences[name] = {k: {'before':base.get(k),'after':v} for k,v in case.items() if v!=base.get(k)}
        atomic_json(RUN / 'config_diff.json', differences)
        manifest = dict(schema_version=1, run_id=RUN.name, created_utc=utcnow(),
            stage_root=str(ROOT), run_dir=str(RUN), ledger_path=str(ledger.path),
            implementation_sha256=implementation, material_path=base['material_path'],
            material_sha256=MATERIAL, tasks=tasks,
            production_gate_task_ids=[t['id'] for t in tasks], gate_comparisons=[],
            qualification_subset_only=True, pilot_design_conditions=0,
            scope='Absorber12 versus prior PML12 boundary implementation; full-observable flat then baseline slot p30 at r32; no production sweep',
            missing_full_qualification=['mesh convergence','strict stop','mirror and straight-slot symmetry','s polarization','L48 boundary geometry'],
            budgets=dict(total_active_seconds=min(172800., ledger.total_active_seconds()+36000.),
                stage_total_active_seconds=172800., phase_solver_seconds=36000.,
                single_solver_seconds=21600., hard_deadline_epoch=SOLVER_DEADLINE,
                hard_deadline_utc='2026-09-24T12:30:00Z', max_parallel=1,
                mpi_ranks=4, memory_fraction=.7),
            rental_expiry_user_confirmed_epoch=RENTAL_EXPIRY,
            rental_note='User confirmed approximately 2026-09-24 21:53 Asia/Shanghai; stop solvers by 20:30, reserve 83 minutes for backup',
            automatic_retry_authorized=False,
            scripts_sha256={p.name:digest_file(p) for p in sorted(ROOT.glob('*.py'))},
            inputs_sha256={str(p.relative_to(RUN)):digest_file(p) for p in (RUN/'inputs').rglob('*') if p.is_file()})
        manifest['manifest_sha256'] = digest_object(manifest)
        verify_manifest(manifest, RUN, implementation)
        atomic_json(RUN/'manifest.json', manifest)
        atomic_json(RUN/'status.json', dict(state='PREPARED', production_qualified=False,
                    updated_utc=utcnow(), hard_deadline_epoch=SOLVER_DEADLINE))
        return {'manifest':str(RUN/'manifest.json'), 'tasks':len(tasks),
                'remaining_hours':(SOLVER_DEADLINE-time.time())/3600, 'solver_launched':False}


class KQueue(Queue):
    def __init__(self):
        super().__init__(RUN/'manifest.json', ('validation',))
        for name, expected in self.manifest['scripts_sha256'].items():
            if digest_file(ROOT/name) != expected:
                raise RuntimeError('K_SCRIPT_CHANGED:'+name)
        for name, expected in self.manifest['inputs_sha256'].items():
            if digest_file(RUN/name) != expected:
                raise RuntimeError('FROZEN_INPUT_CHANGED:'+name)

    def execute(self, task, attempt, admission=None):
        original = self.budget['single_solver_seconds']
        self.budget['single_solver_seconds'] = cap_seconds(task, self.remaining_seconds())
        try:
            return super().execute(task, attempt, admission)
        finally:
            self.budget['single_solver_seconds'] = original

    def checkpoint(self):
        from ti2d.report import generate_report
        from export_k import main as export
        generate_report(self.manifest_path)
        lines = ['# K absorber verification — two cases only', '',
                 'This is not the 20-condition size study. No production gate is unlocked.', '',
                 'Only outer boundary implementation differs from PML12; frozen material and scientific tolerances are unchanged.', '',
                 '| Condition | Status |', '|---|---|']
        for task in self.manifest['tasks']:
            lines.append('| '+task['id']+' | '+self.statuses.get(task['id'],{}).get('status','NOT_RUN')+' |')
        lines += ['', 'Full observables and fixed acceptance checks: per-task result.json.',
                  'Missing production evidence: mesh, stricter stop, mirror/straight-slot, s polarization and L48 checks.',
                  'Solver cutoff: 2026-09-24 20:30 China time. Rental expires approximately 21:53.',
                  '4 MPI ranks, serial only. K cap 10 h within original 48 h; no automatic retries.',
                  'Raw fields: separate checksum bundle, not GitHub. AI/RMB cost: unavailable, not inferred.']
        (RUN/'reports/K_HANDOFF.md').write_text('\n'.join(lines)+'\n')
        export(final=False)


def run():
    os.chdir(ROOT)
    os.environ['MPLCONFIGDIR'] = str(ROOT/'cache/matplotlib')
    os.environ['XDG_CACHE_HOME'] = str(ROOT/'cache')
    os.environ['TMPDIR'] = str(ROOT/'tmp')
    (ROOT/'tmp').mkdir(exist_ok=True)
    try:
        KQueue().run()
    except Exception as exc:
        atomic_json(RUN/'queue_error.json', dict(state='QUEUE_ERROR',error=str(exc),utc=utcnow()))
        raise
    finally:
        from export_k import main as export
        export(final=True)


if __name__ == '__main__':
    if sys.platform != 'linux' or ROOT != SERVER:
        raise RuntimeError('SERVER_ONLY_NO_LOCAL_FDTD')
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('prepare','run'))
    args = parser.parse_args()
    if args.action == 'prepare':
        print(json.dumps(prepare()))
    else:
        run()
