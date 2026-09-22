"""Freeze a finite expiry validation queue; no Meep import or execution."""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import sys
from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.geometry import validate_geometry
from ti2d.materials import load_model, audit_model
from ti2d.queue import FileLock, SolverLedger, verify_manifest

ROOT = Path(__file__).resolve().parent
STAGE = ROOT.parent
RUN = ROOT / 'runs' / 'expiry_20260922_H'
DEADLINE = datetime(2026, 9, 23, 8, 52, 44, tzinfo=timezone.utc).timestamp()
SOURCE_SHA = '3e6b9631a47a3c19040e74fd3c39cfef8df228f45e7bf766782c3e1fb414a85c'


def create():
    if sys.platform != 'linux' or str(STAGE) != '/root/autodl-tmp/meep_sim/meep_screen':
        raise RuntimeError('SERVER_ONLY')
    import json, time
    with FileLock(STAGE / 'budget_ledger.lock'):
        if RUN.exists():
            raise FileExistsError('NO_OVERWRITE_OR_AUTOMATIC_RETRY')
        if time.time() >= DEADLINE - 3600:
            raise RuntimeError('INSUFFICIENT_RENTAL_WINDOW')
        ledger = SolverLedger(STAGE / 'budget_ledger.json')
        if any(x.get('end') is None for x in ledger.data['intervals']):
            raise RuntimeError('OPEN_PRIOR_SOLVER_INTERVAL')
        original = STAGE / 'runs/baseline_r32_20260922_E/cases/baseline_r32_p_plus30.json'
        if digest_file(original) != SOURCE_SHA:
            raise RuntimeError('FROZEN_SOURCE_CASE_CHANGED')
        base = json.loads(original.read_text())
        if digest_file(base['material_path']) != base['material_sha256']:
            raise RuntimeError('FROZEN_MATERIAL_CHANGED')
        implementation = implementation_digest()
        # Sequential, dependency-gated validation only. No automatic size scan.
        specs = [
            ('flat32p', 'flat', 32, 'p', 30, 30, False),
            ('slot32p', 'structure', 32, 'p', 30, 30, False),
            ('strict32p', 'structure', 32, 'p', 30, 30, True),
            ('mirror32p', 'structure', 32, 'p', -30, -30, False),
            ('flat32s', 'flat', 32, 's', 30, 30, False),
            ('slot32s', 'structure', 32, 's', 30, 30, False),
            ('flat24p', 'flat', 24, 'p', 30, 30, False),
            ('slot24p', 'structure', 24, 'p', 30, 30, False),
            ('flat20p', 'flat', 20, 'p', 30, 30, False),
            ('slot20p', 'structure', 20, 'p', 30, 30, False),
        ]
        tasks, diffs = [], {}
        model = load_model(base['material_path'])
        for name, kind, resolution, pol, theta, tilt, strict in specs:
            c = deepcopy(base)
            c.update(id=name, case=kind, resolution=resolution, polarization=pol,
                     emission_theta_deg=theta, strict_stop=strict,
                     intensity_decay=1e-6 if strict else 1e-4, eps_averaging=False,
                     implementation_sha256=implementation, max_solver_seconds=21600.)
            c['geometry']['tilt_deg'] = tilt
            if not validate_geometry(c['geometry'])['valid']:
                raise RuntimeError('INVALID_GEOMETRY:' + name)
            audit_model(model, resolution, c['courant'], c['wavelength_um'])
            path = RUN / 'cases' / (name + '.json')
            atomic_json(path, c)
            tasks.append(dict(id=name, phase='validation', case_path=str(path),
                depends_on=[tasks[-1]['id']] if tasks else [], geometry_id='baseline',
                resolution=resolution, polarization=pol, emission_theta_deg=theta,
                tilt_deg=tilt, strict_stop=strict, gate_role=kind,
                case_fingerprint=digest_object({k:v for k,v in c.items() if k!='id'})))
            diffs[name] = {k:dict(before=base.get(k),after=v) for k,v in c.items() if base.get(k)!=v}
        comparisons = [
            dict(kind='strict', ids=['slot32p','strict32p'], tolerance=.001),
            dict(kind='mirror', ids=['slot32p','mirror32p'], tolerance=.01),
            dict(kind='mesh', ids=['flat32p','flat24p'], tolerance=.01),
            dict(kind='mesh', ids=['slot32p','slot24p'], tolerance=.01),
            dict(kind='mesh', ids=['flat32p','flat24p','flat20p'], tolerance=.01),
            dict(kind='mesh', ids=['slot32p','slot24p','slot20p'], tolerance=.01),
        ]
        manifest = dict(schema_version=1,run_id=RUN.name,created_utc=utcnow(),stage_root=str(ROOT),
            run_dir=str(RUN),ledger_path=str(ledger.path),implementation_sha256=implementation,
            material_path=base['material_path'],material_sha256=base['material_sha256'],
            tasks=tasks,production_gate_task_ids=[t['id'] for t in tasks],gate_comparisons=comparisons,
            qualification_subset_only=True,pilot_design_conditions=0,
            scope='Finite eps_averaging=False validation subset; no qualified size sensitivity yet',
            missing_full_qualification=['straight-slot +/- angles','s mirror','boundary L48 r24','complete production-grid qualification'],
            budgets=dict(total_active_seconds=min(172800.,ledger.total_active_seconds()+72000.),
                stage_total_active_seconds=172800.,single_solver_seconds=21600.,
                hard_deadline_epoch=DEADLINE,hard_deadline_utc='2026-09-23T08:52:44Z',
                max_parallel=1,mpi_ranks=4,memory_fraction=.7),
            expiry_assumption='24h from 2026-09-22T12:52:44Z user report; 4h reserved; not provider-verified',
            history_ledger_sha256=digest_file(ledger.path),automatic_retry_authorized=False)
        manifest['manifest_sha256'] = digest_object(manifest)
        verify_manifest(manifest,RUN,implementation)
        atomic_json(RUN/'ledger_before.json',ledger.data)
        atomic_json(RUN/'config_diff.json',diffs)
        atomic_json(RUN/'manifest.json',manifest)
        print(json.dumps({'manifest':str(RUN/'manifest.json'),'tasks':len(tasks),'hard_deadline_utc':manifest['budgets']['hard_deadline_utc'],'solver_launched':False}))


if __name__ == '__main__':
    create()
