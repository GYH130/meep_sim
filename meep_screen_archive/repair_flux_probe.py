"""Server-only, nonzero cached-field lifecycle diagnostic; never qualification."""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

from ti2d.common import atomic_json, digest_file
from ti2d.materials import audit_model, build_medium, load_model
from ti2d.solver import layout
from ti2d.flux_reference import make_partition, collective_layout_signature, assert_loaded_reference


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    stage = Path(__file__).resolve().parent
    if sys.platform != 'linux' or str(stage) != '/root/autodl-tmp/meep_sim/meep_screen':
        raise RuntimeError('SERVER_ONLY')
    import meep as mp
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
    if comm.size != 4:
        raise RuntimeError('FOUR_RANKS_REQUIRED')
    mp.verbosity(0)
    out = Path(args.output).resolve()
    if not out.is_relative_to(stage / 'repair_diagnostics'):
        raise ValueError('OUTPUT_SCOPE')
    prior = json.loads((stage / 'runs/stage_20260921_B/tasks/flat_r16_p_plus30/attempt_1/result.json').read_text())
    c = prior['case_config']
    refdir = Path(prior['reference_attempt'])
    ref = json.loads((refdir / 'manifest.json').read_text())
    # Read historical data as diagnostic evidence only, never as a new qualified cache.
    for name, sha in ref['evidence_sha256'].items():
        if digest_file(refdir.parent / name) != sha:
            raise ValueError('OLD_REFERENCE_EVIDENCE_CHANGED')
    with np.load(refdir / f'flux_rank{rank:03d}.npz', allow_pickle=False) as d:
        expected = {k: d[k].copy() for k in ('E', 'H')}
    model = load_model(c['material_path'])
    audit_model(model, c['resolution'], c['courant'], c['wavelength_um'])
    ly = layout(c)
    f = 1 / c['wavelength_um']
    kx = -f * math.sin(math.radians(c['emission_theta_deg']))
    pulse = mp.GaussianSource(frequency=f, fwidth=c['source_fwidth_fraction'] * f,
                              cutoff=c['source_cutoff'], is_integrated=True)
    source = mp.Source(pulse, mp.Hz, center=mp.Vector3(0, ly['source']),
                       size=mp.Vector3(ly['period'], 0),
                       amp_func=lambda r: np.exp(2j * np.pi * kx * r.x))
    sim = mp.Simulation(cell_size=mp.Vector3(ly['period'], ly['height']),
                        boundary_layers=[mp.PML(c['pml_um'], direction=mp.Y)],
                        geometry=[mp.Block(center=mp.Vector3(0, 0),
                                           size=mp.Vector3(ly['period'], c['geometry']['thickness_um'], mp.inf),
                                           material=build_medium(model))],
                        sources=[source], resolution=c['resolution'], Courant=c['courant'],
                        dimensions=2, k_point=mp.Vector3(kx, 0), force_complex_fields=True,
                        split_chunks_evenly=True, chunk_layout=str(refdir / 'chunks.h5'))
    monitors = {}
    for name in ('refl', 'trans'):
        vol = mp.Volume(center=mp.Vector3(0, ly[name]), size=mp.Vector3(ly['period'], 0))
        monitors[name] = sim.add_flux(f, 0, 1, mp.FluxRegion(center=vol.center, size=vol.size, direction=mp.Y))
        sim.add_dft_fields([mp.Ex, mp.Hz], [f], where=vol, yee_grid=False, decimation_factor=1)
    sim.add_dft_fields([mp.Ex, mp.Dx, mp.Ey, mp.Dy], [f],
                       where=mp.Volume(center=mp.Vector3(), size=mp.Vector3(ly['period'], c['geometry']['thickness_um'] + 2 / c['resolution'])),
                       yee_grid=True, decimation_factor=1)
    sim.init_sim()
    sim.run(until=0)
    flux = monitors['refl']
    points = []

    def checkpoint(label):
        data = sim.get_flux_data(flux)
        values = {}
        for name in ('E', 'H'):
            a = np.asarray(getattr(data, name))
            b = expected[name]
            values[name] = {'shape': list(a.shape), 'expected_norm': float(np.linalg.norm(b)),
                            'actual_norm': float(np.linalg.norm(a)),
                            'expected_nonzero_range': ([int(np.flatnonzero(b)[0]), int(np.flatnonzero(b)[-1])] if np.any(b) else []),
                            'actual_nonzero_range': ([int(np.flatnonzero(a)[0]), int(np.flatnonzero(a)[-1])] if np.any(a) else []),
                            'max_error_vs_negative': float(np.max(abs(a + b))) if a.shape == b.shape else None}
        record = {'checkpoint': label, 'time': float(sim.meep_time()),
                  'normalized_flux': float(mp.get_fluxes(flux)[0]) / ref['input_power'],
                  'ranks': comm.gather(values, root=0)}
        if rank == 0:
            points.append(record)
            atomic_json(out / 'probe.json', {'status': 'DIAGNOSTIC_NOT_QUALIFICATION',
                        'reference_dir': str(refdir), 'input_power': ref['input_power'],
                        'checkpoints': points})
            print(json.dumps(record), flush=True)

    checkpoint('before_load')
    empty = sim.get_flux_data(flux)
    sim.load_minus_flux_data(flux, type(empty)(E=expected['E'], H=expected['H']))
    checkpoint('after_load')
    sim.run(until=0)
    checkpoint('after_second_zero_step_run')
    audit_model(model, c['resolution'], c['courant'], c['wavelength_um'])
    sim.run(until=2 * c['courant'] / c['resolution'])
    checkpoint('after_two_steps')
    # Nonphysical sentinel identifies each current rank's opaque array index ownership.
    sim.load_flux_data(flux, type(empty)(E=np.ones_like(expected['E']), H=np.ones_like(expected['H'])))
    checkpoint('nonphysical_unit_sentinel_ownership')
    complete = {name: comm.allreduce(a, op=MPI.SUM) for name, a in expected.items()}
    sim.load_minus_flux_data(flux, type(empty)(E=complete['E'], H=complete['H']))
    checkpoint('diagnostic_only_global_composite_load')
    sim.reset_meep()

    def fixed_sim(reference):
        objects = [] if reference else [mp.Block(center=mp.Vector3(),
                         size=mp.Vector3(ly['period'], c['geometry']['thickness_um'], mp.inf),
                         material=build_medium(model))]
        result = mp.Simulation(cell_size=mp.Vector3(ly['period'], ly['height']),
                    boundary_layers=[mp.PML(c['pml_um'], direction=mp.Y)],
                    geometry=objects, sources=[source], resolution=c['resolution'], Courant=c['courant'],
                    dimensions=2, k_point=mp.Vector3(kx, 0), force_complex_fields=True,
                    split_chunks_evenly=True, chunk_layout=make_partition(mp, 4))
        fluxes = {}
        for name in ('refl', 'trans'):
            vol = mp.Volume(center=mp.Vector3(0, ly[name]), size=mp.Vector3(ly['period'], 0))
            fluxes[name] = result.add_flux(f, 0, 1, mp.FluxRegion(center=vol.center, size=vol.size, direction=mp.Y))
            result.add_dft_fields([mp.Ex, mp.Hz], [f], where=vol, yee_grid=False, decimation_factor=1)
        if not reference:
            result.add_dft_fields([mp.Ex, mp.Dx, mp.Ey, mp.Dy], [f],
                where=mp.Volume(center=mp.Vector3(), size=mp.Vector3(ly['period'], c['geometry']['thickness_um'] + 2 / c['resolution'])),
                yee_grid=True, decimation_factor=1)
        result.init_sim()
        result.run(until=0)
        chunkpath = out / ('fixed_reference_chunks.h5' if reference else 'fixed_structure_chunks.h5')
        result.dump_chunk_layout(str(chunkpath))
        comm.Barrier()
        signature = collective_layout_signature(result, chunkpath, comm)
        return result, fluxes['refl'], signature

    # Seed a fresh fixed-layout reference monitor with the observed NONZERO historical
    # data. This is a deterministic serialization regression, not a new physical run
    # or a qualified-cache conversion. A real reference is recomputed in the flat trial.
    sim, flux, fixed_reference_layout = fixed_sim(True)
    empty = sim.get_flux_data(flux)
    sim.load_flux_data(flux, type(empty)(E=complete['E'], H=complete['H']))
    d = sim.get_flux_data(flux)
    expected = {'E': d.E.copy(), 'H': d.H.copy()}
    checkpoint('fixed_reference_nonzero_seed')
    sim.reset_meep()
    sim, flux, fixed_structure_layout = fixed_sim(False)
    if fixed_reference_layout != fixed_structure_layout:
        raise ValueError('FIXED_LAYOUT_REGRESSION_FAILED')
    empty = sim.get_flux_data(flux)
    sim.load_minus_flux_data(flux, type(empty)(E=expected['E'], H=expected['H']))
    audits = {'after_load': assert_loaded_reference(sim, flux, expected, ref['input_power'], mp, comm)}
    checkpoint('fixed_after_load')
    sim.run(until=0)
    audits['after_zero_step_run'] = assert_loaded_reference(sim, flux, expected, ref['input_power'], mp, comm)
    checkpoint('fixed_after_zero_step_run')
    audit_model(model, c['resolution'], c['courant'], c['wavelength_um'])
    sim.run(until=2 * c['courant'] / c['resolution'])
    audits['after_two_steps'] = assert_loaded_reference(sim, flux, expected, ref['input_power'], mp, comm)
    checkpoint('fixed_after_two_steps')
    if rank == 0:
        atomic_json(out / 'fixed_regression.json', {'status': 'PASSED_NONZERO_LIFECYCLE_NOT_PHYSICS',
                    'layout': fixed_reference_layout, 'audits': audits,
                    'source': 'historical field values injected solely as deterministic nonzero test data',
                    'total_timesteps': 4})
    sim.reset_meep()


if __name__ == '__main__':
    main()
