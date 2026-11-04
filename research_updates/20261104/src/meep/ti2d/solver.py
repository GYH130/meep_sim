"""Server-only MPI Meep driver. Every exit is explicit; a timeout is never pass.

Physical vertical z is mapped to Meep y. Emission direction is stored separately
from reciprocal incidence k=(-sin(theta),-cos(theta))*f. All monitor extraction
is collective; ordinary global evidence is written by rank zero only.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
import signal
import time
import traceback
import uuid
from pathlib import Path

import numpy as np

from .common import atomic_json, atomic_npz, digest_file, digest_object, implementation_digest, utcnow
from .diagnostics import StopTracker, fresnel_slab, integrate_absorption, periodic_projection
from .geometry import validate_geometry
from .materials import MaterialStabilityError, audit_model, build_medium, epsilon, load_model
from .flux_reference import make_partition, partition_spec, collective_layout_signature, assert_loaded_reference
from .reference_cache import validated_reference
from .boundary import boundary_spec, build_boundary_layers


class FieldNumericalError(FloatingPointError):
    """Nonfinite or overflow-risk fields: abort, never clip or normalize away."""

    def __init__(self, code, audit):
        self.code = code
        self.audit = audit
        super().__init__(code + ": " + str(audit.get('component', 'unknown')))


def checked_field_power(values, component, contribution_count=1):
    """Check E before squaring, preserving exact finite powers on valid input.

    The bound reserves room for summing the electric components and is solely
    an overflow guard, not a convergence or physical acceptance threshold.
    Zero/subnormal source tails remain valid. Inputs are never clipped.
    """
    if not isinstance(contribution_count, int) or contribution_count < 1:
        raise ValueError('FIELD_POWER_CONTRIBUTION_COUNT_INVALID')
    array = np.asarray(values, dtype=np.complex128)
    audit = {'component': str(component), 'shape': list(array.shape),
             'contribution_count': contribution_count,
             'policy': 'abort_without_clipping', 'physical_acceptance_check': False}
    if array.size == 0:
        raise FieldNumericalError('FIELD_ARRAY_EMPTY', audit)
    if not np.isfinite(array).all():
        audit['nonfinite_count'] = int(np.count_nonzero(~np.isfinite(array)))
        raise FieldNumericalError('FIELD_NONFINITE', audit)
    safe_magnitude = math.sqrt(np.finfo(np.float64).max / contribution_count) / 2
    component_limit = safe_magnitude / math.sqrt(2)
    real_max, imag_max = float(np.abs(array.real).max()), float(np.abs(array.imag).max())
    if real_max > component_limit or imag_max > component_limit:
        audit.update(max_abs_real=real_max, max_abs_imag=imag_max,
                     safe_magnitude_limit=safe_magnitude,
                     per_real_imag_component_limit=component_limit)
        raise FieldNumericalError('FIELD_POWER_OVERFLOW_RISK', audit)
    with np.errstate(over='raise', invalid='raise', under='ignore'):
        try:
            power = np.abs(array) ** 2
        except FloatingPointError as exc:
            raise FieldNumericalError('FIELD_POWER_OVERFLOW', audit) from exc
    if not np.isfinite(power).all():
        raise FieldNumericalError('FIELD_POWER_NONFINITE', audit)
    return power


def layout(c):
    g = c['geometry']
    pml, above, below, thick = c['pml_um'], c['air_above_um'], c['air_below_um'], g['thickness_um']
    height = 2 * pml + above + below + thick
    top = height / 2 - pml
    surface = top - above
    bottom = surface - thick
    source = top - .20 * above
    refl = top - .65 * above
    if 'reflection_monitor_above_surface_um' in c:
        offset = c['reflection_monitor_above_surface_um']
        if (isinstance(offset, bool) or not isinstance(offset, (int, float))
                or not math.isfinite(offset) or offset <= 0):
            raise ValueError('REFLECTION_MONITOR_OFFSET_INVALID')
        refl = surface + offset
        # A reflected-flux monitor must remain in source-free, non-PML air.
        # The layout is part of the reference identity; shifted monitors can
        # never reuse the reference of the historical monitor position.
        tolerance = 1e-12 * max(1., abs(surface), abs(source))
        if not surface + tolerance < refl < source - tolerance:
            raise ValueError('REFLECTION_MONITOR_MUST_BE_BETWEEN_SURFACE_AND_SOURCE')
    return dict(period=g['period_um'], height=height, top=top, surface=surface, bottom=bottom,
                source=source, refl=refl,
                trans=bottom - .5 * below, non_pml_bottom=-height / 2 + pml,
                non_pml_height=height - 2 * pml)


def stopping_layout(c):
    """Keep convergence sampling independent of diagnostic reflection planes."""
    return layout({k: v for k, v in c.items()
                   if k != 'reflection_monitor_above_surface_um'})


def resolve_source_identity(c):
    """Resolve a same-polarization source without changing the default drive.

    An explicit source is a flat-control diagnostic only.  In this 2D x/y
    geometry, s fields are Ez/Hx/Hy and p fields are Hz/Ex/Ey; the permitted
    alternatives use the tangential component at the horizontal source plane.
    Each choice needs its own equally driven empty-cell reference.  No claim
    of equal source amplitude or directional purity follows from this switch.
    """
    pol = c.get('polarization')
    if pol not in ('s', 'p'):
        raise ValueError('POLARIZATION_INVALID')
    default = 'Ez' if pol == 's' else 'Hz'
    explicit = 'diagnostic_source_component' in c
    if explicit and c.get('case') != 'flat':
        raise ValueError('DIAGNOSTIC_SOURCE_REQUIRES_FLAT_CASE')
    component = c['diagnostic_source_component'] if explicit else default
    allowed = ('Ez', 'Hx') if pol == 's' else ('Hz', 'Ex')
    if not isinstance(component, str) or component not in allowed:
        raise ValueError('DIAGNOSTIC_SOURCE_POLARIZATION_MISMATCH')
    return {'polarization': pol, 'component': component,
            'default_component': default,
            'diagnostic_source_component': component if explicit else None,
            'diagnostic_override': explicit}


def result_checks(c, metrics, stopped, analytical=None, stop_reason='CONVERGED'):
    checks = {'stopped': bool(stopped and stop_reason == 'CONVERGED'),
              'case_known': c.get('case') in ('flat', 'structure'),
              'independent_absorption': abs(metrics['A_vol'] - metrics['A_flux']) < .01,
              'reflection_orders': abs(metrics['order_R'] - metrics['R']) < .01,
              'transmission_orders': abs(metrics['order_T'] - metrics['T']) < .01,
              'reference_pseudo_reflection': metrics['reference_pseudo_reflection'] < .001,
              'finite': all(math.isfinite(float(v)) for v in metrics.values()),
              'passive_flux_within_tolerance': all(-.001 < metrics[k] < 1.001 for k in ('R', 'T', 'A_flux'))}
    if c.get('case') == 'flat' and analytical is None:
        checks['fresnel_evidence_present'] = False
    if analytical is not None:
        checks['fresnel_R'] = abs(metrics['R'] - analytical['R']) < .0025
        checks['fresnel_RTA'] = max(abs(metrics[k] - analytical[a]) for k, a in
                                    (('R', 'R'), ('T', 'T'), ('A_flux', 'A'))) < .01
    return checks


def stop_threshold(c, reference):
    return c['intensity_decay'] if reference else c.get('exploratory_structure_intensity_decay', c['intensity_decay'])


def resolve_minimum_post_source_time(c):
    """Optional lower bound on simulated tail time, never a forced stop.

    This applies equally to fresh reference and structure calculations.  A
    missing field preserves the original strict-convergence stopping behavior.
    JSON strings, null, booleans and nonfinite/negative values fail closed.
    """
    value = c.get('minimum_post_source_time_um', 0.0)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError('MINIMUM_POST_SOURCE_TIME_MUST_BE_FINITE_NONNEGATIVE')
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError('MINIMUM_POST_SOURCE_TIME_MUST_BE_FINITE_NONNEGATIVE') from exc
    if not math.isfinite(value) or value < 0:
        raise ValueError('MINIMUM_POST_SOURCE_TIME_MUST_BE_FINITE_NONNEGATIVE')
    return value


def convergence_stop_requested(converged, simulation_time, source_end, minimum_post_source_time):
    """Require both numerical convergence and the matched minimum tail time."""
    return bool(converged and simulation_time >= source_end + minimum_post_source_time)


def classify_result(checks, reason, exploratory=False):
    if reason == 'CONVERGED' and all(checks.values()):
        return 'EXPLORATORY_CONVERGED_NOT_QUALIFIED' if exploratory else 'QUALIFIED'
    return 'NUMERICAL_FAILED' if reason in ('CONVERGED', 'NUMERICAL_FAILED_EARLY') else reason


NO_PER_CASE_CAP_POLICY = 'NO_PER_CASE_CAP_TOTAL_BUDGET_SUPERVISED'
NO_TIME_CAP_POLICY = 'NO_PER_CASE_OR_CUMULATIVE_TIME_CAP_SUPERVISED'
LEGACY_WALLTIME_POLICY = 'PER_CASE_CAPPED_AT_6H'


def resolve_wall_limit(c, environ=None):
    """Resolve only elapsed wall time, never numerical stopping/qualification.

    Unlimited *per-case* time requires an explicit case policy, JSON null cap,
    and matching controller environment. The controller remains responsible for
    the selected external accounting/resource policy. The distinct no-time-cap
    policy removes cumulative enforcement but still requires process supervision.
    Legacy inputs retain their six-hour ceiling. Invalid limits fail closed.
    """
    environ = os.environ if environ is None else environ
    policy = c.get('walltime_policy', LEGACY_WALLTIME_POLICY)
    env_limit = environ.get('TI2D_WALL_TIMEOUT_SECONDS', 21600.)
    if policy in (NO_PER_CASE_CAP_POLICY, NO_TIME_CAP_POLICY):
        if ('max_solver_seconds' not in c or c['max_solver_seconds'] is not None
                or env_limit != 'none'):
            raise ValueError('NO_PER_CASE_CAP_REQUIRES_EXPLICIT_NULL_AND_CONTROLLER_NONE')
        return None
    if policy != LEGACY_WALLTIME_POLICY:
        raise ValueError('UNKNOWN_WALLTIME_POLICY')

    def finite_positive(value):
        if isinstance(value, bool):
            raise ValueError('WALLTIME_LIMIT_MUST_BE_FINITE_POSITIVE')
        try:
            seconds = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError('WALLTIME_LIMIT_MUST_BE_FINITE_POSITIVE') from exc
        if not math.isfinite(seconds) or seconds <= 0:
            raise ValueError('WALLTIME_LIMIT_MUST_BE_FINITE_POSITIVE')
        return seconds

    return min(finite_positive(env_limit),
               finite_positive(c.get('max_solver_seconds', 21600.)), 21600.)


def walltime_stop_requested(terminated, elapsed_seconds, wall_limit):
    return bool(terminated or (wall_limit is not None and elapsed_seconds >= wall_limit))


def main(case_path, output):
    # Import only on the server's numerical entry point.
    import meep as mp
    from mpi4py import MPI
    from .meep_arrays import native_volume_pairs, plane_fields

    comm = MPI.COMM_WORLD
    rank, ranks = comm.rank, comm.size
    master = rank == 0
    mp.verbosity(0)
    c = json.loads(Path(case_path).read_text())
    out = Path(output).resolve()
    stage = Path(__file__).resolve().parents[1]
    if stage not in out.parents:
        raise ValueError('All solver outputs must remain under meep_screen')
    if master:
        out.mkdir(parents=True, exist_ok=True)
    comm.Barrier()
    started = time.time()
    start_clock = time.monotonic()
    wall_limit = resolve_wall_limit(c)
    terminated = [False]
    signal.signal(signal.SIGTERM, lambda *_: terminated.__setitem__(0, True))
    implementation = implementation_digest()
    base_result = {'schema_version': 1, 'id': c['id'], 'started_utc': utcnow(),
                   'case_config': c, 'case_sha256': digest_file(case_path),
                   'implementation_sha256': implementation, 'mpi_ranks': ranks,
                   'material_label': 'Ordal_Route_A_frozen', 'oxide_layer_in_model': False,
                   'walltime_policy': c.get('walltime_policy', LEGACY_WALLTIME_POLICY),
                   'per_case_wall_limit_seconds': wall_limit,
                   'external_total_budget_supervision_required': wall_limit is None and c.get('walltime_policy') != NO_TIME_CAP_POLICY,
                   'external_process_supervision_required': wall_limit is None,
                   'cumulative_time_cap_seconds': None if c.get('walltime_policy') == NO_TIME_CAP_POLICY else 'external_policy',
                   'restart_semantics': 'New run or qualified reference reuse; no field-state timestep resume',
                   'scope': 'specified Ti optical model; 2D only; real oxidation/sample temperature unknown'}
    active_sim = None
    cache_lock = None
    try:
        if ranks != c['mpi_ranks']:
            raise ValueError('MPI_RANKS_MISMATCH')
        minimum_post_source_time = resolve_minimum_post_source_time(c)
        base_result['minimum_post_source_time_um'] = minimum_post_source_time
        if c.get('stop_normalization_policy') != 'strict_historical_peak':
            raise ValueError('THIS_ISOLATED_STAGE_REQUIRES_STRICT_HISTORICAL_PEAK')
        if c.get('exploratory_structure_intensity_decay') is not None and not (
                c['exploratory_structure_intensity_decay'] == 1e-3 and c['intensity_decay'] == 1e-4
                and c.get('qualification_scope') == 'EXPLORATORY_ONLY_NOT_FORMAL_QUALIFICATION'
                and c.get('case') == 'structure' and ranks == 16):
            raise ValueError('EXPLORATORY_THRESHOLD_SCOPE_MISMATCH')
        if c.get('implementation_sha256') != implementation:
            raise ValueError('IMPLEMENTATION_CHANGED_SINCE_MANIFEST')
        material_path = Path(c['material_path'])
        if digest_file(material_path) != c['material_sha256']:
            raise ValueError('MATERIAL_FINGERPRINT_MISMATCH')
        model = load_model(material_path)
        stability = audit_model(model, c['resolution'], c['courant'], c['wavelength_um'])
        geometry_report = validate_geometry(c['geometry'])
        if not geometry_report['valid']:
            raise ValueError('GEOMETRY_NOT_QUALIFIED')
        if c['oxide_layer_in_model']:
            raise ValueError('OXIDE_OUTSIDE_AUTHORIZED_SCOPE')
        if c.get('case') not in ('flat', 'structure'):
            raise ValueError('CASE_INVALID')
        boundary = boundary_spec(c)
        base_result.update(stability=stability, geometry_audit=geometry_report,
                           boundary=boundary,
                           boundary_qualification='Individual-case checks only; boundary independence not assumed')
        ly = layout(c)
        f = 1 / c['wavelength_um']
        theta = math.radians(c['emission_theta_deg'])
        kx = -f * math.sin(theta)
        incident = (kx, -f * math.cos(theta), 0.)
        pol = c['polarization']
        if pol not in ('s', 'p'):
            raise ValueError('POLARIZATION_INVALID')
        components = {'Ex': mp.Ex, 'Hz': mp.Hz} if pol == 'p' else {'Ez': mp.Ez, 'Hx': mp.Hx}
        ecomps = (mp.Ex, mp.Ey) if pol == 'p' else (mp.Ez,)
        edpairs = {'Ex': (mp.Ex, mp.Dx), 'Ey': (mp.Ey, mp.Dy)} if pol == 'p' else {'Ez': (mp.Ez, mp.Dz)}
        source_identity = resolve_source_identity(c)
        source_component = getattr(mp, source_identity['component'])
        pulse = mp.GaussianSource(frequency=f, fwidth=c['source_fwidth_fraction'] * f,
                                  cutoff=c['source_cutoff'], is_integrated=True)
        source_end = float(pulse.start_time + 2 * pulse.width * pulse.cutoff)
        sources = [mp.Source(pulse, source_component, center=mp.Vector3(0, ly['source']),
                             size=mp.Vector3(ly['period'], 0),
                             amp_func=lambda r: np.exp(2j * np.pi * kx * r.x))]
        base_result.update(layout=ly, source_end_time=source_end,
                           source_identity=source_identity,
                           fixed_mpi_partition=partition_spec(ranks, ly['period']),
                           emission_direction=[math.sin(theta), 0., math.cos(theta)],
                           reciprocal_incident_k_simulation_xy=incident,
                           physical_to_simulation_axes='physical (x,z) -> Meep (x,y)')

        if c.get('eps_averaging') is not False:
            raise ValueError('EXPIRY_CAMPAIGN_REQUIRES_EXPLICIT_EPS_AVERAGING_FALSE')

        def simulation(reference):
            objects = []
            if not reference:
                medium = build_medium(model)
                objects = [mp.Block(center=mp.Vector3(0, .5 * (ly['surface'] + ly['bottom'])),
                                    size=mp.Vector3(ly['period'], c['geometry']['thickness_um'], mp.inf), material=medium)]
                if c['case'] == 'structure':
                    for shift in geometry_report['periodic_shifts']:
                        vertices = [mp.Vector3(x + shift * ly['period'], ly['surface'] - d)
                                    for x, d in geometry_report['vertices_x_depth']]
                        objects.append(mp.Prism(vertices=vertices, height=mp.inf,
                                                axis=mp.Vector3(0, 0, 1), material=mp.air))
            return mp.Simulation(cell_size=mp.Vector3(ly['period'], ly['height']),
                                 boundary_layers=build_boundary_layers(mp, c),
                                 geometry=objects, sources=sources, resolution=c['resolution'],
                                 Courant=c['courant'], dimensions=2, k_point=mp.Vector3(kx, 0),
                                 force_complex_fields=True, split_chunks_evenly=True,
                                 eps_averaging=c['eps_averaging'],
                                 chunk_layout=make_partition(mp, ranks, ly['period']))

        def monitors(sim, include_volume):
            result = {}
            for name in ('refl', 'trans'):
                volume = mp.Volume(center=mp.Vector3(0, ly[name]), size=mp.Vector3(ly['period'], 0))
                result[name] = sim.add_flux(f, 0, 1, mp.FluxRegion(center=volume.center, size=volume.size, direction=mp.Y))
                result[name + '_dft'] = sim.add_dft_fields(list(components.values()), [f], where=volume,
                                                          yee_grid=False, decimation_factor=1)
            if include_volume:
                # A small air halo contains both native interface grids. No PML/source inside.
                volume = mp.Volume(center=mp.Vector3(0, .5 * (ly['surface'] + ly['bottom'])),
                                   size=mp.Vector3(ly['period'], c['geometry']['thickness_um'] + 2 / c['resolution']))
                result['volume'] = sim.add_dft_fields([v for pair in edpairs.values() for v in pair],
                                                       [f], where=volume, yee_grid=True, decimation_factor=1)
            sim.init_sim()
            sim.run(until=0)
            return result

        # Probes cover both vacuum propagation and slot cavity; coordinates are folded periodically.
        # Moving a measurement plane must not also move its stopping probe.
        stop_layout = stopping_layout(c)
        probes = {'front_air': (0., stop_layout['refl']), 'back_air': (0., ly['trans'])}
        if c['case'] == 'structure':
            alpha = math.radians(c['geometry']['tilt_deg'])
            for fraction in (.25, .50, .75):
                distance = c['geometry']['axis_length_um'] * fraction
                px = (distance * math.sin(alpha) + ly['period'] / 2) % ly['period'] - ly['period'] / 2
                probes['cavity_' + str(fraction)] = (px, ly['surface'] - distance * math.cos(alpha))
        probe_vectors = {key: mp.Vector3(*xy) for key, xy in probes.items()}
        non_pml = mp.Volume(center=mp.Vector3(), size=mp.Vector3(ly['period'], ly['non_pml_height']))

        def run_until_stopped(sim, mons, reference, incident_power, prefix, reference_stop=None):
            reference_peak = (None if reference else
                              reference_stop['source_on_non_pml_envelope_peak'])
            tracker = StopTracker(source_end, window=c['stop_window_um'],
                                  consecutive_windows=c['stop_consecutive_windows'],
                                  decay_threshold=stop_threshold(c, reference), rt_tolerance=c['flux_drift_tolerance'],
                                  reference_intensity_peak=reference_peak,
                                  reference_provenance=(None if reference else 'qualified_reference_source_on_non_pml_peak'),
                                  reference_floor_fraction=c.get('reference_intensity_floor_fraction', 1.0),
                                  normalization_policy=c['stop_normalization_policy'])
            phase_start = time.monotonic()
            peak_flux = [0.]
            last_sample = [-1.]
            last_write = [0.]
            reason = [None]

            def collect(s):
                now = float(s.meep_time())
                if now <= last_sample[0] + 1e-8:
                    return
                last_sample[0] = now
                # Collective array extraction occurs on every MPI rank, not inside rank-zero guards.
                intensity = None
                field_failure = None
                for component in ecomps:
                    values = s.get_array(component=component, vol=non_pml, cmplx=True)
                    try:
                        power = checked_field_power(values, str(component), len(ecomps))
                        intensity = power if intensity is None else intensity + power
                    except FieldNumericalError as exc:
                        field_failure = field_failure or {'code': exc.code, 'audit': exc.audit}
                # Finish every collective field extraction before synchronizing an
                # error, so one rank cannot leave others inside a Meep reduction.
                probe_power = {}
                for name, vec in probe_vectors.items():
                    probe_values = [s.get_field_point(component, vec) for component in ecomps]
                    try:
                        probe_power[name] = float(checked_field_power(
                            probe_values, 'probe_' + name, len(ecomps)).sum())
                    except FieldNumericalError as exc:
                        field_failure = field_failure or {'code': exc.code, 'audit': exc.audit}
                failures = comm.allgather(field_failure)
                failure = next((item for item in failures if item is not None), None)
                if failure is not None:
                    failure['audit'].update(simulation_time=now, phase=prefix,
                                            mpi_synchronized=True)
                    raise FieldNumericalError(failure['code'], failure['audit'])
                envelope = float(np.max(intensity))
                raw_r = float(mp.get_fluxes(mons['refl'])[0])
                raw_t = float(mp.get_fluxes(mons['trans'])[0])
                peak_flux[0] = max(peak_flux[0], abs(raw_r), abs(raw_t), 1e-280)
                normalizer = peak_flux[0] if reference else incident_power
                rval = -raw_r / normalizer if reference else raw_r / normalizer
                tval = -raw_t / normalizer
                tracker.sample(now, probe_power, envelope, rval, tval)
                tracker.trace[-1]['phase_wall_elapsed_s'] = time.monotonic() - phase_start
                # Only assess gross failures after the source and three complete stable
                # flux windows. A preloaded negative incident flux before arrival is normal.
                windows = tracker.summary()['windows'][-3:]
                if (not reference and now >= source_end + 3 * c['stop_window_um'] and
                        len(windows) == 3 and all(w['covered'] and w['rt_stable'] for w in windows) and
                        (rval < -.01 or rval > 1.01 or tval < -.01 or tval > 1.01)):
                    reason[0] = 'NUMERICAL_FAILED_EARLY'
                if master and time.monotonic() - last_write[0] >= 30:
                    atomic_json(out / 'heartbeat.json', dict(phase=prefix, simulation_time=now,
                                source_end=source_end, wall_elapsed_s=time.monotonic() - start_clock,
                                stop_converged=tracker.converged,
                                minimum_post_source_time_um=minimum_post_source_time,
                                eligible_to_stop=convergence_stop_requested(
                                    tracker.converged, now, source_end, minimum_post_source_time),
                                last_sample=tracker.trace[-1]))
                    atomic_json(out / (prefix + '_stop_partial.json'), tracker.summary())
                    last_write[0] = time.monotonic()

            def stop(s):
                flag = comm.allreduce(int(walltime_stop_requested(
                    terminated[0], time.monotonic() - start_clock, wall_limit)), op=MPI.MAX)
                if flag:
                    reason[0] = 'TIME_LIMIT_UNQUALIFIED'
                    return True
                if reason[0] == 'NUMERICAL_FAILED_EARLY':
                    return True
                target = c.get('benchmark_structure_time') if not reference else None
                if target is not None:
                    if not (math.isfinite(target) and target > source_end):
                        raise ValueError('BENCHMARK_MUST_EXTEND_BEYOND_SOURCE_END')
                    if s.meep_time() >= target:
                        reason[0] = 'PERFORMANCE_DIAGNOSTIC_ONLY'
                        return True
                    return False
                if c.get('smoke_max_time') is not None and s.meep_time() >= c['smoke_max_time']:
                    reason[0] = 'SMOKE_ONLY'
                    return True
                return convergence_stop_requested(
                    tracker.converged, float(s.meep_time()), source_end, minimum_post_source_time)

            sim.run(mp.at_every(c['sample_interval_um'], collect), until=stop)
            collect(sim)
            status = reason[0] or ('CONVERGED' if tracker.converged else 'TIME_LIMIT_UNQUALIFIED')
            summary = tracker.summary()
            summary.update(trace=tracker.trace, stop_reason=status, elapsed_s=time.monotonic() - phase_start,
                           minimum_post_source_time_um=minimum_post_source_time,
                           minimum_eligible_stop_time=source_end + minimum_post_source_time,
                           minimum_post_source_time_reached=float(sim.meep_time()) >= source_end + minimum_post_source_time,
                           intensity_definition='sum of |E_component|^2; strict own historical peak, not current tail',
                           flux_trace_semantics=('reference front/back downward flux divided by historical flux peak; R trace is not reflectance'
                                                 if reference else 'R=scattered upward flux/input; T=downward transmitted flux/input'))
            if master:
                atomic_json(out / (prefix + '_stop.json'), summary)
            return summary, status

        # A cache key covers every reference-affecting condition, build and MPI partition.
        reference_identity = {k: c[k] for k in ('resolution', 'courant', 'wavelength_um', 'polarization',
                                'emission_theta_deg', 'source_fwidth_fraction', 'source_cutoff', 'pml_um',
                                'air_above_um', 'air_below_um', 'intensity_decay', 'sample_interval_um',
                                'stop_window_um', 'stop_consecutive_windows', 'flux_drift_tolerance')}
        reference_identity.update(layout=ly, boundary=boundary, implementation_sha256=implementation,
                                  source_identity=source_identity,
                                  eps_averaging=c['eps_averaging'],
                                  material_sha256=c['material_sha256'], meep_version=mp.__version__, mpi_ranks=ranks,
                                  fixed_partition=partition_spec(ranks, ly['period']),
                                  minimum_post_source_time_um=minimum_post_source_time,
                                  stop_normalization_policy=c['stop_normalization_policy'],
                                  reference_intensity_floor_fraction=c.get('reference_intensity_floor_fraction', 1.0))
        refkey = digest_object(reference_identity)
        cache_dir = stage / 'cache' / 'references' / refkey
        if master:
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_lock = (cache_dir / '.lock').open('a+')
            fcntl.flock(cache_lock, fcntl.LOCK_EX)
        comm.Barrier()
        cached = None
        if master and not c.get('smoke_max_time'):
            try:
                pointer = json.loads((cache_dir / 'qualified.json').read_text())
                candidate = json.loads((cache_dir / pointer['manifest']).read_text())
                if validated_reference(candidate, cache_dir, reference_identity, ranks):
                    cached = candidate
            except (OSError, ValueError, KeyError, TypeError, AttributeError):
                cached = None
        cached = comm.bcast(cached, root=0)
        ref_manifest = cached
        if cached is None:
            attempt_name = comm.bcast(uuid.uuid4().hex if master else None, root=0)
            attempt_dir = cache_dir / attempt_name
            if master:
                attempt_dir.mkdir()
            comm.Barrier()
            active_sim = simulation(True)
            refmons = monitors(active_sim, False)
            chunk_file = attempt_dir / 'chunks.h5'
            active_sim.dump_chunk_layout(str(chunk_file))  # Collective Meep-owned HDF5 operation.
            comm.Barrier()
            ref_layout_signature = collective_layout_signature(active_sim, chunk_file, comm)
            ref_stop, ref_reason = run_until_stopped(active_sim, refmons, True, None, 'reference')
            xref, fields_ref, arrays_ref = plane_fields(active_sim, refmons['refl_dft'], components, ly['period'], kx, c['resolution'])
            xback, fields_back, arrays_back = plane_fields(active_sim, refmons['trans_dft'], components, ly['period'], kx, c['resolution'])
            front = periodic_projection(xref, fields_ref, ly['period'], kx, f, pol)
            back = periodic_projection(xback, fields_back, ly['period'], kx, f, pol)
            raw_incident = float(mp.get_fluxes(refmons['refl'])[0])
            input_power = -raw_incident
            # The reference should travel downward below the source. Opposite power is the background artifact.
            ref_checks = {'stop': bool(ref_stop['converged'] and ref_reason == 'CONVERGED'), 'downward_flux': input_power > 0,
                          'direction': front['down_power'] > front['up_power'],
                          'plane_closure': abs(front['net_up_power'] - raw_incident) / max(input_power, 1e-280) < .01,
                          'empty_pseudo_reflection': front['up_power'] / max(front['down_power'], 1e-280) < .001,
                          'empty_transmission': abs(-float(mp.get_fluxes(refmons['trans'])[0]) / max(input_power, 1e-280) - 1) < .01}
            ref_data = active_sim.get_flux_data(refmons['refl'])
            # Opaque flux buffers are rank-local; each rank owns exactly its named file.
            atomic_npz(attempt_dir / f'flux_rank{rank:03d}.npz', E=ref_data.E, H=ref_data.H)
            if master:
                atomic_npz(attempt_dir / 'reference_planes.npz', x=xref,
                           source_component=np.asarray(source_identity['component']),
                           source_polarization=np.asarray(pol),
                           source_diagnostic_override=np.asarray(source_identity['diagnostic_override']),
                           **{'front_' + key: value for key, value in fields_ref.items()},
                           **{'front_raw_' + key: value for key, value in arrays_ref.items()},
                           **{'back_raw_' + key: value for key, value in arrays_back.items()})
            comm.Barrier()
            ref_manifest = {'status': 'QUALIFIED' if all(ref_checks.values()) else ref_reason if ref_reason != 'CONVERGED' else 'NUMERICAL_FAILED',
                            'identity': reference_identity, 'checks': ref_checks, 'stop': ref_stop,
                            'chunk_layout_signature': ref_layout_signature,
                            'input_power': input_power, 'front_projection': front, 'back_projection': back,
                            'pseudo_reflection': front['up_power'] / max(front['down_power'], 1e-280),
                            'attempt': attempt_name, 'reference_cached': False}
            if master:
                ref_manifest['evidence_sha256'] = {str(p.relative_to(cache_dir)): digest_file(p) for p in attempt_dir.iterdir() if p.is_file()}
                atomic_json(attempt_dir / 'manifest.json', ref_manifest)
                if ref_manifest['status'] == 'QUALIFIED':
                    atomic_json(cache_dir / 'qualified.json', {'manifest': attempt_name + '/manifest.json'})
            active_sim.reset_meep()
            active_sim = None
            if ref_manifest['status'] != 'QUALIFIED':
                final = dict(base_result, status=ref_manifest['status'], stop=ref_stop, reference=ref_manifest,
                             failure_stage='reference', checks=ref_checks, elapsed_s=time.monotonic() - start_clock)
                if master:
                    atomic_json(out / 'result.json', final)
                return 2
        else:
            input_power = float(cached['input_power'])
            attempt_dir = cache_dir / cached['attempt']
            chunk_file = attempt_dir / 'chunks.h5'
        comm.Barrier()
        if master and cache_lock is not None:
            fcntl.flock(cache_lock, fcntl.LOCK_UN)
            cache_lock.close()
            cache_lock = None
        with np.load(attempt_dir / 'reference_planes.npz', allow_pickle=False) as data:
            xref = data['x'].copy()
            fields_ref = {name: data['front_' + name].copy() for name in components}
        with np.load(attempt_dir / f'flux_rank{rank:03d}.npz', allow_pickle=False) as data:
            # Keep the Meep namedtuple type, without depending on its export location.
            flux_buffers = (data['E'].copy(), data['H'].copy())
        active_sim = simulation(False)
        mons = monitors(active_sim, True)
        structure_chunk_file = out / 'structure_chunks.h5'
        active_sim.dump_chunk_layout(str(structure_chunk_file))
        comm.Barrier()
        structure_layout_signature = collective_layout_signature(active_sim, structure_chunk_file, comm)
        if structure_layout_signature != ref_manifest['chunk_layout_signature']:
            raise ValueError('REFERENCE_SPATIAL_CHUNK_OR_PROCESS_OWNERSHIP_MISMATCH')
        empty_data = active_sim.get_flux_data(mons['refl'])
        if comm.allreduce(int(empty_data.E.shape != flux_buffers[0].shape or
                              empty_data.H.shape != flux_buffers[1].shape), op=MPI.MAX):
            raise ValueError('REFERENCE_CHUNK_OR_FLUX_SHAPE_MISMATCH')
        active_sim.load_minus_flux_data(mons['refl'], type(empty_data)(*flux_buffers))
        preload_audit = assert_loaded_reference(active_sim, mons['refl'], flux_buffers, input_power, mp, comm)
        active_sim.run(until=0)
        preload_audit_after_run = assert_loaded_reference(active_sim, mons['refl'], flux_buffers, input_power, mp, comm)
        if master:
            atomic_json(out / 'reference_preload_audit.json', {
                'after_load': preload_audit, 'after_zero_step_run': preload_audit_after_run,
                'chunk_layout_signature': structure_layout_signature})
        structure_stop, reason = run_until_stopped(active_sim, mons, False, input_power, 'structure', ref_manifest['stop'])
        xr, fr, rawr = plane_fields(active_sim, mons['refl_dft'], components, ly['period'], kx, c['resolution'])
        xt, ft, rawt = plane_fields(active_sim, mons['trans_dft'], components, ly['period'], kx, c['resolution'])
        if not np.array_equal(xr, xref):
            raise ValueError('REFERENCE_DFT_COORDINATE_MISMATCH')
        scattered = {name: fr[name] - fields_ref[name] for name in components}
        projection_r = periodic_projection(xr, scattered, ly['period'], kx, f, pol)
        projection_t = periodic_projection(xt, ft, ly['period'], kx, f, pol)
        geom_for_mask = dict(c['geometry'], case=c['case'])
        pairs, raw_volume, volume_meta = native_volume_pairs(active_sim, mons['volume'], edpairs,
                                                           c['resolution'], ly['period'], kx,
                                                           geom_for_mask, ly['surface'])
        absorption = integrate_absorption(pairs, f, input_power)
        R = float(mp.get_fluxes(mons['refl'])[0]) / input_power
        T = -float(mp.get_fluxes(mons['trans'])[0]) / input_power
        metrics = {'R': R, 'T': T, 'A_flux': 1 - R - T, 'A_vol': absorption['A'],
                   'order_R': projection_r['up_power'] / input_power,
                   'order_T': projection_t['down_power'] / input_power,
                   'reference_pseudo_reflection': float(ref_manifest['pseudo_reflection'])}
        analytical = fresnel_slab(complex(epsilon(model, f)), c['wavelength_um'], c['emission_theta_deg'], pol,
                                 c['geometry']['thickness_um']) if c['case'] == 'flat' else None
        checks = result_checks(c, metrics, structure_stop['converged'], analytical, reason)
        exploratory = c.get('exploratory_structure_intensity_decay') is not None
        status = classify_result(checks, reason, exploratory)
        if master:
            atomic_npz(out / 'monitor_planes.npz', x_refl=xr, x_trans=xt,
                       source_component=np.asarray(source_identity['component']),
                       source_polarization=np.asarray(pol),
                       source_diagnostic_override=np.asarray(source_identity['diagnostic_override']),
                       **{'total_' + name: arr for name, arr in fr.items()},
                       **{'incident_' + name: arr for name, arr in fields_ref.items()},
                       **{'scattered_' + name: arr for name, arr in scattered.items()},
                       **{'trans_' + name: arr for name, arr in ft.items()},
                       **{'raw_refl_' + name: arr for name, arr in rawr.items()},
                       **{'raw_trans_' + name: arr for name, arr in rawt.items()})
            atomic_npz(out / 'volume_ED.npz', **raw_volume)
            atomic_json(out / 'volume_metadata.json', volume_meta)
            atomic_json(out / 'orders.json', {'reflection': projection_r, 'transmission': projection_t})
            immutable_evidence = ('monitor_planes.npz', 'volume_ED.npz', 'volume_metadata.json',
                                  'orders.json', 'structure_stop.json', 'reference_stop.json',
                                  'reference_preload_audit.json', 'structure_chunks.h5')
            evidence = {name: digest_file(out / name) for name in immutable_evidence if (out / name).is_file()}
            final = dict(base_result, status=status, metrics=metrics, checks=checks, stop=structure_stop,
                         analytical=analytical, absorption=absorption, reference_cache_key=refkey,
                         reference_attempt=str(attempt_dir), reference_reused=cached is not None,
                         reference_checks=ref_manifest['checks'], reference_stop=ref_manifest['stop'],
                         reference_preload_audit=preload_audit_after_run,
                         numerical_error_bounds={'validated': False, 'reason': 'aggregate mesh/strict-stop evidence required'},
                         formal_qualification=False if exploratory else None,
                         checks_scope='exploratory stop plus unchanged physical checks' if exploratory else 'strict stop and physical checks',
                         evidence_sha256=evidence, elapsed_s=time.monotonic() - start_clock, completed_utc=utcnow(),
                         sum_identity_not_conservation_proof=True, unclipped=True)
            atomic_json(out / 'result.json', final)
        return 0 if status == 'QUALIFIED' else 2
    except Exception as exc:
        if master:
            if hasattr(exc, 'audit'):
                failure_name = ('numerical_guard_failure.json'
                                if isinstance(exc, (MaterialStabilityError, FieldNumericalError))
                                else 'reference_preload_failure.json')
                atomic_json(out / failure_name, exc.audit)
            atomic_json(out / 'result.json', dict(base_result, status='CRASHED', error=str(exc),
                        failure_code=getattr(exc, 'code', 'UNCLASSIFIED_SOLVER_EXCEPTION'),
                        failure_stage=('material_preflight' if isinstance(exc, MaterialStabilityError)
                                       else 'field_sampling' if isinstance(exc, FieldNumericalError)
                                       else 'solver'),
                        traceback=traceback.format_exc(), elapsed_s=time.monotonic() - start_clock,
                        stop={'converged': False}, completed_utc=utcnow()))
        return 3
    finally:
        if active_sim is not None:
            active_sim.reset_meep()
        if cache_lock is not None:
            cache_lock.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    raise SystemExit(main(args.case, args.output))
