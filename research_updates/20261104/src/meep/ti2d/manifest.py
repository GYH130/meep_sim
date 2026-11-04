"""Build the locked finite campaign, with no runtime expansion or refitting."""
from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path

from .common import atomic_json, digest_file, digest_object, implementation_digest, utcnow
from .geometry import validate_geometry
from .materials import load_model, validate_fit


def create(run_dir):
    root = Path(__file__).resolve().parents[1]
    run_dir = Path(run_dir).resolve()
    if (run_dir / 'manifest.json').exists():
        raise FileExistsError('Manifest already locked; use status/resume, not regenerate')
    run_dir.mkdir(parents=True, exist_ok=True)
    material = root / 'assets/material/route_a_material.json'
    model = load_model(material)
    fit = validate_fit(model, material.with_name('ti_ordal_raw.csv'), material.with_name('fit_report.json'))
    atomic_json(run_dir / 'material_recheck.json', fit)
    if not fit['valid']:
        raise ValueError('FROZEN_MATERIAL_RECHECK_FAILED')
    base = dict(period_um=50., surface_fill=.692820323, axis_length_um=40.,
                tilt_deg=30., thickness_um=60., min_wall_um=1., residual_um=10.)
    geometries = {'baseline': base, 'period_high': dict(base, period_um=58.),
                  'duty_high': dict(base, surface_fill=.8),
                  'length_high': dict(base, axis_length_um=48.),
                  'tilt_high': dict(base, tilt_deg=45.),
                  'straight': dict(base, tilt_deg=0.), 'mirror': dict(base, tilt_deg=-30.)}
    qualified_geometry = {name: validate_geometry(g) for name, g in geometries.items()}
    atomic_json(run_dir / 'geometries.json', {'definitions': geometries, 'audit': qualified_geometry})
    if not all(a['valid'] for a in qualified_geometry.values()):
        raise ValueError('MANIFEST_GEOMETRY_FAILED')
    tasks, comparisons = [], []
    code_hash = implementation_digest()
    seen = {}

    def add(geometry_id, resolution, pol, theta, role, phase='validation', strict=False, flat=False):
        label = f'{"flat" if flat else geometry_id}_r{resolution}_{pol}_{"plus" if theta > 0 else "minus"}30'
        label += '_strict' if strict else ''
        if phase == 'pilot':
            label = 'pilot_' + label
        cfg = {'schema_version': 1, 'id': label, 'geometry_id': geometry_id,
               'geometry': deepcopy(geometries[geometry_id]), 'case': 'flat' if flat else 'structure',
               'material_path': str(material), 'material_sha256': digest_file(material),
               'material_label': 'Ordal_Route_A_frozen', 'implementation_sha256': code_hash,
               'resolution': resolution, 'courant': .25, 'wavelength_um': 10.5,
               'polarization': pol, 'emission_theta_deg': theta,
               'oxide_layer_in_model': False, 'mpi_ranks': 4,
               'source_fwidth_fraction': .2, 'source_cutoff': 5.,
               'pml_um': 6., 'air_above_um': 12., 'air_below_um': 12.,
               'sample_interval_um': 2., 'stop_window_um': 20.,
               'stop_consecutive_windows': 3, 'intensity_decay': 1e-6 if strict else 1e-4,
               'flux_drift_tolerance': .001,
               'reference_intensity_floor_fraction': 1.0,
               'stop_normalization_policy': 'max(own_historic_peak, reference_source_on_non_pml_peak)',
               'mpi_partition_policy': 'fixed_binary_Y0_X0_ranks0123',
               'max_solver_seconds': 21600., 'strict_stop': strict,
               'qualification_scope': '10.5um; emission ±30deg; s/p; listed 2D geometries only'}
        config_path = run_dir / 'cases' / (label + '.json')
        atomic_json(config_path, cfg)
        fingerprint = digest_object({k: v for k, v in cfg.items() if k != 'id'})
        task = {'id': label, 'phase': phase, 'case_path': str(config_path), 'depends_on': [],
                'gate_role': role, 'geometry_id': geometry_id, 'resolution': resolution,
                'polarization': pol, 'emission_theta_deg': theta,
                'tilt_deg': geometries[geometry_id]['tilt_deg'], 'strict_stop': strict,
                'case_fingerprint': fingerprint}
        if fingerprint in seen:
            task['reused_from'] = seen[fingerprint]
            task['depends_on'] = [seen[fingerprint]]
        else:
            seen[fingerprint] = label
        tasks.append(task)
        return label

    baseline_ids, length_ids = {}, {}
    # Front-load cheapest diagnostic planar cases. No D00/D04 reproduction.
    for pol in ('p', 's'):
        ids = [add('baseline', r, pol, 30., 'flat', flat=True) for r in (16, 20, 24)]
        comparisons.append(dict(kind='mesh', ids=ids, tolerance=.01, metrics=['R', 'T', 'A_flux']))
        neg = add('baseline', 20, pol, -30., 'flat', flat=True)
        comparisons.append(dict(kind='mirror', ids=[ids[1], neg], tolerance=.01, metrics=['R', 'T', 'A_flux']))
    for pol in ('p', 's'):
        ids = [add('baseline', r, pol, 30., 'baseline') for r in (16, 20, 24)]
        comparisons.append(dict(kind='mesh', ids=ids, tolerance=.01, metrics=['R', 'T', 'A_flux']))
        baseline_ids[(pol, 30.)] = ids[1]
        baseline_ids[(pol, -30.)] = add('baseline', 20, pol, -30., 'baseline')
    strict_id = add('baseline', 20, 'p', 30., 'baseline', strict=True)
    comparisons.append(dict(kind='strict', ids=[baseline_ids[('p', 30.)], strict_id],
                            tolerance=.001, metrics=['R', 'T', 'A_flux']))
    for pol in ('p', 's'):
        ids = [add('straight', 20, pol, angle, 'straight') for angle in (30., -30.)]
        comparisons.append(dict(kind='mirror', ids=ids, tolerance=.01, metrics=['R', 'T', 'A_flux']))
        for theta in (30., -30.):
            mirror = add('mirror', 20, pol, -theta, 'mirror')
            comparisons.append(dict(kind='mirror', ids=[baseline_ids[(pol, theta)], mirror],
                                    tolerance=.01, metrics=['R', 'T', 'A_flux']))
    for pol in ('p', 's'):
        for theta in (30., -30.):
            ids = [add('length_high', r, pol, theta, 'length_boundary') for r in (20, 24)]
            length_ids[(pol, theta)] = ids[0]
            comparisons.append(dict(kind='mesh_boundary', ids=ids, tolerance=.01, metrics=['R', 'T', 'A_flux']))
    validation_ids = [task['id'] for task in tasks]
    for geo in ('baseline', 'period_high', 'duty_high', 'length_high', 'tilt_high'):
        for pol in ('p', 's'):
            for theta in (30., -30.):
                add(geo, 20, pol, theta, 'pilot', phase='pilot')
    manifest = {'schema_version': 1, 'run_id': run_dir.name, 'created_utc': utcnow(),
                'stage_root': str(root), 'run_dir': str(run_dir), 'implementation_sha256': code_hash,
                'ledger_path': str(root / 'budget_ledger.json'),
                'material_path': str(material), 'material_sha256': digest_file(material),
                'budgets': {'total_active_seconds': 172800, 'single_solver_seconds': 21600,
                            'max_parallel': 4, 'memory_fraction': .7, 'mpi_ranks': 4},
                'tasks': tasks, 'production_gate_task_ids': validation_ids,
                'concurrency_probe': {'serial_task_id': 'flat_r20_p_plus30',
                                      'task_ids': ['flat_r20_p_minus30', 'flat_r16_s_plus30']},
                'gate_comparisons': comparisons, 'pilot_design_conditions': 20,
                'legacy_status': 'Three D15c long cases remain TIME_LIMIT_UNQUALIFIED; no legacy cache reuse.',
                'scope': 'specified Ti optical model; 2D slanted slots, not 3D holes or measured high-temperature samples',
                'interpretation': 'Local one-sided machining-parameter sensitivities; fixed F changes width with P; fixed opening changes normal width with tilt.',
                'approved_max_targeted_repair_rounds': 2, 'repairs_performed': 0}
    manifest['manifest_sha256'] = digest_object(manifest)
    atomic_json(run_dir / 'manifest.json', manifest)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('run_dir')
    args = parser.parse_args()
    m = create(args.run_dir)
    print(f"Locked {len(m['production_gate_task_ids'])} validation tasks and 20 pilot conditions; exact-case reuse enabled.")
