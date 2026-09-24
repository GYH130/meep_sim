"""Freeze L boundary-strength screening inputs; never launch Meep here.

Narrow cases are uniform planar diagnostics, NOT a 2 um scientific period
study. The same unmodified K worker and acceptance tolerances are used.
"""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import json
import shutil
import subprocess
import sys
import time

from ti2d.common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from ti2d.geometry import validate_geometry
from ti2d.materials import load_model, audit_model
from ti2d.queue import FileLock, SolverLedger, verify_manifest

ROOT = Path(__file__).resolve().parent
SERVER = Path('/root/autodl-tmp/meep_sim/meep_screen/boundary_screen_L')
RUN = ROOT / 'runs/boundary_20260924_L'
SOLVER_DEADLINE = datetime(2026, 9, 24, 12, 30, tzinfo=timezone.utc).timestamp()
RENTAL_EXPIRY = datetime(2026, 9, 24, 13, 53, tzinfo=timezone.utc).timestamp()
K_CORE = '9c2b217a329f816195b201d9bae76a6bcc5d3144ded21e031610560e55ad0578'
MATERIAL = '4791b3868810fb4a3fa93c8469b488650270c0f47a18804c9b26bd7354153b5e'
K_PSEUDO_REFLECTION = 0.0017166551733512714
K_INPUT_POWER = 36607.7710324763
CONTROL_ID = 'L_narrow_control_1e15'
CANDIDATE_IDS = ['L_narrow_candidate_1e8', 'L_narrow_candidate_1e6']
FULL_IDS = ['L_full_flat_1e8', 'L_full_flat_1e6']
NARROW_GEOMETRY = dict(period_um=2., surface_fill=.2, axis_length_um=.5,
                       tilt_deg=30, thickness_um=60., min_wall_um=1., residual_um=10.)


def build_cases(base, material_path, implementation):
    """Return every possible condition before selection; no adaptive additions."""
    if not validate_geometry(NARROW_GEOMETRY)['valid']:
        raise ValueError('INVALID_NARROW_PLACEHOLDER_GEOMETRY')
    specifications = [(CONTROL_ID, 1e-15, True, 1200.),
                      (CANDIDATE_IDS[0], 1e-8, True, 1200.),
                      (CANDIDATE_IDS[1], 1e-6, True, 1200.),
                      (FULL_IDS[0], 1e-8, False, 14400.),
                      (FULL_IDS[1], 1e-6, False, 14400.)]
    built = []
    for name, reflectivity, narrow, cap in specifications:
        case = deepcopy(base)
        case.update(id=name, case='flat', material_path=str(material_path),
                    implementation_sha256=implementation,
                    max_solver_seconds=cap - 120.,
                    geometry_id='uniform_planar_diagnostic_width_2um' if narrow else 'baseline',
                    boundary_kind='absorber', pml_um=12.,
                    boundary_profile=dict(name='quadratic', R_asymptotic=reflectivity, mean_stretch=1.0))
        if narrow:
            case['geometry'] = deepcopy(NARROW_GEOMETRY)
        case['diagnostic_scope'] = (
            'Uniform planar computational width only; placeholder slot is not instantiated; '
            'not a scientific period or geometry-sensitivity condition.' if narrow else
            'Full-width 50 um planar qualification only; no slot or production scan authorized.')
        if not validate_geometry(case['geometry'])['valid']:
            raise ValueError('INVALID_CASE_GEOMETRY:' + name)
        built.append((name, case, cap))
    return built


def control_proof(result, case, source_path):
    """Small exact-scalar baseline, linked to original immutable result hash."""
    ref = result['reference']
    required = ('stop', 'downward_flux', 'direction', 'plane_closure', 'empty_transmission')
    if result['status'] != 'NUMERICAL_FAILED' or result.get('failure_stage') != 'reference':
        raise ValueError('K_CONTROL_NOT_EXPECTED_REFERENCE_FAILURE')
    if any(ref['checks'].get(k) is not True for k in required):
        raise ValueError('K_CONTROL_OTHER_CHECK_FAILED')
    if ref['checks'].get('empty_pseudo_reflection') is not False:
        raise ValueError('K_CONTROL_PSEUDO_GATE_NOT_FAILED')
    if abs(ref['pseudo_reflection'] - K_PSEUDO_REFLECTION) > 1e-14:
        raise ValueError('K_PSEUDO_REFERENCE_CHANGED')
    if abs(ref['input_power'] - K_INPUT_POWER) > 1e-8:
        raise ValueError('K_INPUT_POWER_CHANGED')
    period = float(case['geometry']['period_um'])
    if period != 50.:
        raise ValueError('K_REFERENCE_WIDTH_CHANGED')
    return dict(source_result_path=str(source_path), source_result_sha256=digest_file(source_path),
                status=result['status'], failure_stage=result['failure_stage'],
                checks=ref['checks'], reference_status=ref['status'],
                identity=ref['identity'], pseudo_reflection=ref['pseudo_reflection'],
                input_power=ref['input_power'], computational_width_um=period,
                input_power_per_width=ref['input_power'] / period,
                front_projection=ref['front_projection'], back_projection=ref['back_projection'],
                evidence_sha256=ref.get('evidence_sha256', {}),
                control_comparison=dict(pseudo_reflection_absolute_tolerance=1e-5,
                                        input_power_per_width_relative_tolerance=1e-3,
                                        required_true_checks=list(required),
                                        expected_status='NUMERICAL_FAILED',
                                        diagnostic_reproduction_is_qualification=False))


def prepare():
    if sys.platform != 'linux' or ROOT != SERVER:
        raise RuntimeError('SERVER_ONLY_PREPARATION')
    with FileLock(ROOT.parent / 'budget_ledger.lock'):
        if RUN.exists() or (ROOT / 'assets/material').exists():
            raise FileExistsError('SINGLE_USE_L_INPUTS_ALREADY_EXIST')
        ledger = SolverLedger(ROOT.parent / 'budget_ledger.json')
        if any(record.get('end') is None for record in ledger.data['intervals']):
            raise RuntimeError('OPEN_PRIOR_SOLVER_INTERVAL')
        if min(SOLVER_DEADLINE-time.time(), 172800-ledger.total_active_seconds()) < 1200:
            raise RuntimeError('INSUFFICIENT_WINDOW_NO_SOLVER_LAUNCHED')
        history = ROOT.parent / 'boundary_absorber_K/runs/boundary_20260924_K'
        case_path = history / 'cases/K_flat32p.json'
        result_path = history / 'tasks/K_flat32p/attempt_1/result.json'
        base = json.loads(case_path.read_text())
        result = json.loads(result_path.read_text())
        implementation = implementation_digest()
        if implementation != K_CORE or base['implementation_sha256'] != K_CORE:
            raise RuntimeError('UNMODIFIED_K_CORE_REQUIRED')
        if base['material_sha256'] != MATERIAL or digest_file(base['material_path']) != MATERIAL:
            raise RuntimeError('MATERIAL_CHANGED')
        proof = control_proof(result, base, result_path)
        (RUN / 'inputs').mkdir(parents=True)
        shutil.copytree(Path(base['material_path']).parent, ROOT / 'assets/material')
        material_path = ROOT / 'assets/material/route_a_material.json'
        shutil.copyfile(case_path, RUN / 'inputs/K_flat32p_case.json')
        atomic_json(RUN / 'inputs/K_reference_scalars.json', proof)
        for name in ('H_flat32p_case.json', 'H_slot32p_case.json',
                     'H_flat32p_result.json', 'H_slot32p_result.json'):
            src = history / 'inputs' / name
            if src.is_file():
                shutil.copyfile(src, RUN / 'inputs' / name)
        shutil.copyfile(ledger.path, RUN / 'inputs/ledger_before.json')
        repo = ROOT.parents[1]
        for filename, args in (
            ('git_status_before.txt', ['status', '--porcelain']),
            ('protected_changes.patch', ['diff', '--binary', 'HEAD'])):
            (RUN / 'inputs' / filename).write_bytes(subprocess.run(
                ['git', '-C', str(repo)] + args, check=True, capture_output=True).stdout)
        tracked = subprocess.run(['git', '-C', str(repo), 'diff', '--name-only', 'HEAD'],
                                 check=True, capture_output=True, text=True).stdout.splitlines()
        atomic_json(RUN / 'inputs/protected_fingerprints.json', {
            name: digest_file(repo/name) for name in tracked if (repo/name).is_file()})
        atomic_json(RUN / 'inputs/directory_presence.json', {
            str(p): 'present' if p.exists() else 'not found'
            for p in (repo/'sources', repo/'review', repo/'src')})
        model = load_model(material_path)
        tasks, differences = [], {}
        for name, case, cap in build_cases(base, material_path, implementation):
            audit_model(model, case['resolution'], case['courant'], case['wavelength_um'])
            path = RUN / ('configs/future' if name in FULL_IDS else 'cases') / (name + '.json')
            atomic_json(path, case)
            task = dict(id=name, phase='validation', case_path=str(path),
                              depends_on=[], single_cap_seconds=cap,
                              geometry_id=case['geometry_id'], gate_role='flat',
                              resolution=32, polarization='p', emission_theta_deg=30,
                              tilt_deg=30, strict_stop=False,
                              case_fingerprint=digest_object({k:v for k,v in case.items() if k!='id'}))
            if name not in FULL_IDS:
                tasks.append(task)
            differences[name] = {k: dict(before=base.get(k), after=v)
                                 for k,v in case.items() if v != base.get(k)}
        atomic_json(RUN / 'config_diff.json', differences)
        manifest = dict(
            schema_version=1, run_id=RUN.name, created_utc=utcnow(), stage_root=str(ROOT),
            run_dir=str(RUN), ledger_path=str(ledger.path), implementation_sha256=implementation,
            material_path=str(material_path), material_sha256=MATERIAL, tasks=tasks,
            production_gate_task_ids=[], gate_comparisons=[],
            qualification_subset_only=True, pilot_design_conditions=0,
            scope='Fixed two-strength Absorber12 narrow planar boundary screening with replay control only; full-width and slot deferred because rental window is insufficient',
            l_selection=dict(control_id=CONTROL_ID, candidate_ids=CANDIDATE_IDS,
                             full_ids=[], future_full_ids=FULL_IDS,
                             full_status='DEFERRED_RENTAL_WINDOW',
                             first_qualified_candidate_only=True,
                             candidate_reference_pseudo_reflection_max=5e-4,
                             abort_after_timeout_or_crash=True,
                             reference_proof_path=str(RUN/'inputs/K_reference_scalars.json'),
                             control_status='DIAGNOSTIC_REPRODUCED_NOT_QUALIFIED',
                             skip_second_candidate_after_first_pass=True,
                             stop_on_control_mismatch=True, no_slot=True),
            missing_full_qualification=['baseline slot stop','mesh convergence','strict stop',
                'mirror and straight-slot symmetry','s polarization','L48 boundary geometry'],
            budgets=dict(total_active_seconds=min(172800., ledger.total_active_seconds()+3600.),
                         stage_total_active_seconds=172800., phase_solver_seconds=3600.,
                         single_solver_seconds=21600., hard_deadline_epoch=SOLVER_DEADLINE,
                         hard_deadline_utc='2026-09-24T12:30:00Z', max_parallel=1,
                         mpi_ranks=4, memory_fraction=.7),
            rental_expiry_user_confirmed_epoch=RENTAL_EXPIRY,
            rental_note='Expiry approximately 2026-09-24 21:53 China time; solver cutoff20:30; 83 minutes reserved for backup.',
            automatic_retry_authorized=False,
            scripts_sha256={p.name:digest_file(p) for p in sorted(ROOT.glob('*.py'))},
            inputs_sha256={str(p.relative_to(RUN)):digest_file(p)
                           for p in (RUN/'inputs').rglob('*') if p.is_file()},
            future_configs_sha256={str(p.relative_to(RUN)):digest_file(p)
                                   for p in (RUN/'configs/future').glob('*.json')},
            material_files_sha256={str(p.relative_to(ROOT)):digest_file(p)
                                   for p in (ROOT/'assets/material').rglob('*') if p.is_file()})
        manifest['manifest_sha256'] = digest_object(manifest)
        verify_manifest(manifest, RUN, implementation)
        atomic_json(RUN/'manifest.json', manifest)
        atomic_json(RUN/'status.json', dict(state='PREPARED', production_qualified=False,
                    updated_utc=utcnow(), hard_deadline_epoch=SOLVER_DEADLINE))
        return dict(manifest=str(RUN/'manifest.json'), tasks=len(tasks),
                    remaining_hours=(SOLVER_DEADLINE-time.time())/3600, solver_launched=False)


if __name__ == '__main__':
    print(json.dumps(prepare()))
