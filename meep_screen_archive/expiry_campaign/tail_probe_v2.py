"""J diagnostic probe: frozen I physics, explicit per-case cap and separate terminal sample."""
from pathlib import Path
import argparse
import json
import math
import os
import signal
import sys
import time
import traceback

import numpy as np
from ti2d.common import atomic_json, atomic_npz, digest_file, implementation_digest
from ti2d.flux_reference import make_partition
from ti2d.geometry import validate_geometry
from ti2d.materials import audit_model, build_medium, load_model
from ti2d.solver import layout

ROOT = Path(__file__).resolve().parent
SERVER = Path('/root/autodl-tmp/meep_sim/meep_screen/expiry_campaign')


def probe_positions(c):
    ly = layout(c)
    points = {'front_air': (0., ly['refl']), 'back_air': (0., ly['trans'])}
    a = math.radians(c['geometry']['tilt_deg'])
    for frac in (.25, .5, .75):
        d = c['geometry']['axis_length_um'] * frac
        x = (d * math.sin(a) + ly['period']/2) % ly['period'] - ly['period']/2
        points['cavity_' + str(frac)] = (x, ly['surface'] - d * math.cos(a))
    points.update(front_near_pml=(0., ly['top']-1.),
                  back_near_slab=(0., ly['bottom']-1.),
                  back_near_pml=(0., ly['non_pml_bottom']+1.),
                  back_off_axis=(ly['period']/4, ly['trans']))
    return points


def main(case_path, output):
    if sys.platform != 'linux' or ROOT != SERVER:
        raise RuntimeError('FDTD_SERVER_ONLY')
    import meep as mp
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    master = comm.rank == 0
    mp.verbosity(0)
    c = json.loads(Path(case_path).read_text())
    out = Path(output).resolve()
    if ROOT not in out.parents:
        raise ValueError('OUTPUT_OUTSIDE_DATA_PROJECT')
    if comm.size != 4 or c['mpi_ranks'] != 4:
        raise ValueError('EXACT_FOUR_MPI_RANKS_REQUIRED')
    if c['implementation_sha256'] != implementation_digest():
        raise ValueError('FROZEN_CORE_CHANGED')
    if digest_file(c['material_path']) != c['material_sha256']:
        raise ValueError('FROZEN_MATERIAL_CHANGED')
    if c['eps_averaging'] is not False or c['case'] != 'structure' or c['polarization'] != 'p':
        raise ValueError('DIAGNOSTIC_CONTRACT_MISMATCH')
    model = load_model(c['material_path'])
    stability = audit_model(model, c['resolution'], c['courant'], c['wavelength_um'])
    geom = validate_geometry(c['geometry'])
    if not geom['valid']:
        raise ValueError('INVALID_GEOMETRY')
    ly = layout(c)
    f = 1 / c['wavelength_um']
    kx = -f * math.sin(math.radians(c['emission_theta_deg']))
    pulse = mp.GaussianSource(frequency=f, fwidth=c['source_fwidth_fraction']*f,
                             cutoff=c['source_cutoff'], is_integrated=True)
    source_end = float(pulse.start_time + 2*pulse.width*pulse.cutoff)
    objects = [mp.Block(center=mp.Vector3(0, (ly['surface']+ly['bottom'])/2),
                        size=mp.Vector3(ly['period'], c['geometry']['thickness_um'], mp.inf),
                        material=build_medium(model))]
    for shift in geom['periodic_shifts']:
        vertices = [mp.Vector3(x+shift*ly['period'], ly['surface']-d)
                    for x, d in geom['vertices_x_depth']]
        objects.append(mp.Prism(vertices=vertices, height=mp.inf,
                               axis=mp.Vector3(0, 0, 1), material=mp.air))
    sim = mp.Simulation(cell_size=mp.Vector3(ly['period'], ly['height']),
        boundary_layers=[mp.PML(c['pml_um'], direction=mp.Y)], geometry=objects,
        sources=[mp.Source(pulse, mp.Hz, center=mp.Vector3(0, ly['source']),
                           size=mp.Vector3(ly['period'], 0),
                           amp_func=lambda r: np.exp(2j*np.pi*kx*r.x))],
        resolution=c['resolution'], Courant=c['courant'], dimensions=2,
        k_point=mp.Vector3(kx, 0), force_complex_fields=True,
        split_chunks_evenly=True, eps_averaging=False, chunk_layout=make_partition(mp, 4))
    points = probe_positions(c)
    names = list(points)
    components = {'Ex': mp.Ex, 'Ey': mp.Ey, 'Hz': mp.Hz}
    line_y = {'front': ly['refl'], 'back': ly['trans'], 'back_near_pml': ly['non_pml_bottom']+1.}
    line_vols = {n: mp.Volume(center=mp.Vector3(0, y), size=mp.Vector3(ly['period'], 0))
                 for n, y in line_y.items()}
    full_vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(ly['period'], ly['non_pml_height']))
    interrupted = [False]
    signal.signal(signal.SIGTERM, lambda *_: interrupted.__setitem__(0, True))
    began = time.monotonic()
    wall = min(21600., float(c['max_solver_seconds']), float(os.environ['TAIL_WALL_SECONDS']))
    if not 60. <= wall <= 21600.:
        raise ValueError('INVALID_DIAGNOSTIC_WALL_CAP')
    trace, points_t, point_values, lines_t, line_values = [], [], [], [], {}
    chunks, spatial = [], []
    count, last_time, last_flush, last_write = [0], [-1.], [0.], [0.]
    xline = {}
    reason = [None]
    if master:
        out.mkdir(parents=True, exist_ok=True)
        atomic_json(out/'metadata.json', dict(case=c, layout=ly, probe_names=names,
            probe_xy=points, line_y=line_y, point_component_names=list(components),
            source_end=source_end, geometry_audit=geom, stability=stability,
            point_dt=.25, line_dt=2., spatial_dt=20., line_grid_stride=4,
            snapshot_grid_stride=8, intensity_definition='abs(Ex)^2+abs(Ey)^2',
            diagnostic_only=True, qualification=False,
            passive_monitors_omitted='Production flux/DFT/volume accumulators; fields/source/material/grid unchanged',
            restart='New zero-field simulation; no field-state restart', meep_version=mp.__version__))
    comm.Barrier()

    def flush():
        if not master or not points_t:
            return
        path = out / ('waveforms_%04d.npz' % len(chunks))
        atomic_npz(path, point_t=np.asarray(points_t), point_fields=np.asarray(point_values),
                   line_t=np.asarray(lines_t),
                   **{'x_'+n: x for n, x in xline.items()},
                   **{n: np.asarray(a) for n, a in line_values.items()})
        chunks.append(dict(file=path.name, sha256=digest_file(path),
                           first_time=points_t[0], last_time=points_t[-1]))
        points_t.clear(); point_values.clear(); lines_t.clear(); line_values.clear()
        atomic_json(out/'waveform_manifest.json', chunks)
        atomic_json(out/'power_trace.json', trace)
        atomic_json(out/'spatial_summary.json', spatial)

    def collect(s):
        now = float(s.meep_time())
        if now <= last_time[0]+1e-8:
            return
        last_time[0] = now
        count[0] += 1
        values = np.asarray([[s.get_field_point(comp, mp.Vector3(*points[n]))
                              for comp in components.values()] for n in names], complex)
        if not np.isfinite(values).all():
            reason[0] = 'NONFINITE_FIELD'
            return
        if master:
            points_t.append(now); point_values.append(values)
        if count[0] % 8 == 0:
            powers = np.sum(np.abs(values[:, :2])**2, axis=1)
            if not np.isfinite(powers).all():
                reason[0] = 'NONFINITE_FIELD'; return
            if master:
                trace.append(dict(time=now, probe_power=dict(zip(names, powers.tolist()))))
                lines_t.append(now)
            for n, vol in line_vols.items():
                metadata = s.get_array_metadata(vol=vol)
                x = np.asarray(metadata[0])
                keep = np.flatnonzero((x >= -ly['period']/2-1e-9) & (x < ly['period']/2-1e-9))
                if len(keep) != round(ly['period']*c['resolution']):
                    raise ValueError('LINE_CANONICAL_PERIOD_MISMATCH')
                selection = keep[::4]
                for cn, comp in components.items():
                    a = np.asarray(s.get_array(component=comp, vol=vol, cmplx=True)).squeeze()
                    if a.shape != x.shape or not np.isfinite(a).all():
                        raise ValueError('LINE_ARRAY_METADATA_MISMATCH_OR_NONFINITE')
                    if master:
                        xline[n] = x[selection]
                        line_values.setdefault(n+'_'+cn, []).append(a[selection])
        if count[0] % 80 == 0:
            ex = np.asarray(s.get_array(component=mp.Ex, vol=full_vol, cmplx=True))
            ey = np.asarray(s.get_array(component=mp.Ey, vol=full_vol, cmplx=True))
            x, y, z, w = s.get_array_metadata(vol=full_vol)
            x, y = np.asarray(x), np.asarray(y)
            if ex.shape != (len(x), len(y)) or ey.shape != ex.shape:
                raise ValueError('FULL_FIELD_METADATA_MISMATCH')
            power = np.abs(ex)**2 + np.abs(ey)**2
            if not np.isfinite(power).all():
                reason[0] = 'NONFINITE_FIELD'; return
            if master:
                ix, iy = np.unravel_index(np.argmax(power), power.shape)
                spatial.append(dict(time=now, max_power=float(power[ix, iy]),
                    peak_xy=[float(x[ix]), float(y[iy])],
                    front_air_peak=float(power[:, y>ly['surface']].max()),
                    back_air_peak=float(power[:, y<ly['bottom']].max())))
                if count[0] % 400 == 0:
                    atomic_npz(out/('intensity_t%06d.npz' % round(now)),
                               x=x[::8], y=y[::8], intensity=power[::8, ::8],
                               time=np.asarray(now))
        if master and (now-last_flush[0] >= 100 or time.monotonic()-last_write[0] >= 180):
            flush(); last_flush[0] = now; last_write[0] = time.monotonic()
        if master and count[0] % 80 == 0:
            atomic_json(out/'heartbeat.json', dict(simulation_time=now, source_end=source_end,
                wall_elapsed_s=time.monotonic()-began, status='RUNNING_DIAGNOSTIC',
                latest_power=trace[-1] if trace else None))

    def stop(s):
        if comm.allreduce(int(interrupted[0] or time.monotonic()-began>=wall), op=MPI.MAX):
            reason[0] = 'TIME_LIMIT_UNQUALIFIED'; return True
        if reason[0]:
            return True
        return s.meep_time() >= c['diagnostic_stop_time']

    error = None
    try:
        sim.init_sim(); sim.run(until=0)
        sim.run(mp.at_every(.25, collect), until=stop)
        # Do not append an off-grid termination time to the uniform FFT stream.
        terminal = np.asarray([[sim.get_field_point(comp, mp.Vector3(*points[n]))
                                for comp in components.values()] for n in names], complex)
        if not np.isfinite(terminal).all():
            reason[0] = 'NONFINITE_FIELD'
        if master and np.isfinite(terminal).all():
            atomic_npz(out/'terminal_sample.npz', point_t=np.asarray([float(sim.meep_time())]),
                       point_fields=terminal[None, :, :])
            atomic_json(out/'terminal_sample.json', dict(file='terminal_sample.npz',
                sha256=digest_file(out/'terminal_sample.npz'), time=float(sim.meep_time())))
        state = reason[0] or 'DIAGNOSTIC_COMPLETE'
    except Exception:
        error = traceback.format_exc()
        state = 'CRASHED'
    finally:
        if master:
            flush()
            atomic_json(out/'result.json', dict(status=state, qualification=False,
                simulation_time=float(sim.meep_time()), source_end=source_end,
                elapsed_s=time.monotonic()-began, waveform_chunks=len(chunks),
                case_sha256=digest_file(case_path), error=error))
        sim.reset_meep()
    return 0 if state=='DIAGNOSTIC_COMPLETE' else 2


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--case', required=True); p.add_argument('--output', required=True)
    a = p.parse_args()
    raise SystemExit(main(a.case, a.output))
