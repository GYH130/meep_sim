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
from .materials import audit_model, build_medium, epsilon, load_model
from .flux_reference import make_partition, partition_spec, collective_layout_signature, assert_loaded_reference
from .reference_cache import validated_reference
from .boundary import boundary_spec, build_boundary_layers


def layout(c):
    g = c['geometry']
    pml, above, below, thick = c['pml_um'], c['air_above_um'], c['air_below_um'], g['thickness_um']
    height = 2 * pml + above + below + thick
    top = height / 2 - pml
    surface = top - above
    bottom = surface - thick
    return dict(period=g['period_um'], height=height, top=top, surface=surface, bottom=bottom,
                source=top - .20 * above, refl=top - .65 * above,
                trans=bottom - .5 * below, non_pml_bottom=-height / 2 + pml,
                non_pml_height=height - 2 * pml)


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
    wall_limit = min(float(os.environ.get('TI2D_WALL_TIMEOUT_SECONDS', 21600)),
                     float(c.get('max_solver_seconds', 21600)), 21600.)
    terminated = [False]
    signal.signal(signal.SIGTERM, lambda *_: terminated.__setitem__(0, True))
    implementation = implementation_digest()
    base_result = {'schema_version': 1, 'id': c['id'], 'started_utc': utcnow(),
                   'case_config': c, 'case_sha256': digest_file(case_path),
                   'implementation_sha256': implementation, 'mpi_ranks': ranks,
                   'material_label': 'Ordal_Route_A_frozen', 'oxide_layer_in_model': False,
                   'restart_semantics': 'New run or qualified reference reuse; no field-state timestep resume',
                   'scope': 'specified Ti optical model; 2D only; real oxidation/sample temperature unknown'}
    active_sim = None
    cache_lock = None
    try:
        if ranks != c['mpi_ranks']:
            raise ValueError('MPI_RANKS_MISMATCH')
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
        source_component = mp.Hz if pol == 'p' else mp.Ez
        pulse = mp.GaussianSource(frequency=f, fwidth=c['source_fwidth_fraction'] * f,
                                  cutoff=c['source_cutoff'], is_integrated=True)
        source_end = float(pulse.start_time + 2 * pulse.width * pulse.cutoff)
        sources = [mp.Source(pulse, source_component, center=mp.Vector3(0, ly['source']),
                             size=mp.Vector3(ly['period'], 0),
                             amp_func=lambda r: np.exp(2j * np.pi * kx * r.x))]
        base_result.update(layout=ly, source_end_time=source_end,
                           fixed_mpi_partition=partition_spec(ranks),
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
                                 chunk_layout=make_partition(mp, ranks))

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
        probes = {'front_air': (0., ly['refl']), 'back_air': (0., ly['trans'])}
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
                                  decay_threshold=c['intensity_decay'], rt_tolerance=c['flux_drift_tolerance'],
                                  reference_intensity_peak=reference_peak,
                                  reference_provenance=(None if reference else 'qualified_reference_source_on_non_pml_peak'),
                                  reference_floor_fraction=c.get('reference_intensity_floor_fraction', 1.0))
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
                for component in ecomps:
                    values = s.get_array(component=component, vol=non_pml, cmplx=True)
                    power = np.abs(values) ** 2
                    intensity = power if intensity is None else intensity + power
                envelope = float(np.max(intensity))
                probe_power = {name: float(sum(abs(s.get_field_point(component, vec)) ** 2 for component in ecomps))
                               for name, vec in probe_vectors.items()}
                raw_r = float(mp.get_fluxes(mons['refl'])[0])
                raw_t = float(mp.get_fluxes(mons['trans'])[0])
                peak_flux[0] = max(peak_flux[0], abs(raw_r), abs(raw_t), 1e-280)
                normalizer = peak_flux[0] if reference else incident_power
                rval = -raw_r / normalizer if reference else raw_r / normalizer
                tval = -raw_t / normalizer
                tracker.sample(now, probe_power, envelope, rval, tval)
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
                                stop_converged=tracker.converged, last_sample=tracker.trace[-1]))
                    atomic_json(out / (prefix + '_stop_partial.json'), tracker.summary())
                    last_write[0] = time.monotonic()

            def stop(s):
                flag = comm.allreduce(int(terminated[0] or time.monotonic() - start_clock >= wall_limit), op=MPI.MAX)
                if flag:
                    reason[0] = 'TIME_LIMIT_UNQUALIFIED'
                    return True
                if reason[0] == 'NUMERICAL_FAILED_EARLY':
                    return True
                if c.get('smoke_max_time') is not None and s.meep_time() >= c['smoke_max_time']:
                    reason[0] = 'SMOKE_ONLY'
                    return True
                return tracker.converged

            sim.run(mp.at_every(c['sample_interval_um'], collect), until=stop)
            collect(sim)
            status = reason[0] or ('CONVERGED' if tracker.converged else 'TIME_LIMIT_UNQUALIFIED')
            summary = tracker.summary()
            summary.update(trace=tracker.trace, stop_reason=status, elapsed_s=time.monotonic() - phase_start,
                           intensity_definition='sum of |E_component|^2; reference-excitation-anchored historic-peak floor',
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
                                  eps_averaging=c['eps_averaging'],
                                  material_sha256=c['material_sha256'], meep_version=mp.__version__, mpi_ranks=ranks,
                                  fixed_partition=partition_spec(ranks),
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
        status = 'QUALIFIED' if all(checks.values()) else (
            'NUMERICAL_FAILED' if reason in ('CONVERGED', 'NUMERICAL_FAILED_EARLY') else reason)
        if master:
            atomic_npz(out / 'monitor_planes.npz', x_refl=xr, x_trans=xt,
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
                         evidence_sha256=evidence, elapsed_s=time.monotonic() - start_clock, completed_utc=utcnow(),
                         sum_identity_not_conservation_proof=True, unclipped=True)
            atomic_json(out / 'result.json', final)
        return 0 if status == 'QUALIFIED' else 2
    except Exception as exc:
        if master:
            if hasattr(exc, 'audit'):
                atomic_json(out / 'reference_preload_failure.json', exc.audit)
            atomic_json(out / 'result.json', dict(base_result, status='CRASHED', error=str(exc),
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
