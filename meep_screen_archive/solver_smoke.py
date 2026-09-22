"""Server-only four-rank API smoke; four simulations, exactly two steps each.

This is never a physical qualification, a cache entry, or a production result.
The launcher must wrap this process in the campaign's cumulative SolverLedger.
No Meep import or solver initialization occurs when this module is imported.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np

from ti2d.common import atomic_json, atomic_npz, digest_file, utcnow
from ti2d.geometry import validate_geometry
from ti2d.materials import audit_model, build_medium, epsilon, load_model
from ti2d.meep_arrays import native_volume_pairs, plane_fields
from ti2d.solver import layout


STATUS = "SMOKE_ONLY_NEVER_QUALIFIED"


def _assert_server(stage):
    if sys.platform != "linux":
        raise RuntimeError("SERVER_ONLY: this smoke must not run on the local machine")
    if not Path(stage).resolve().is_relative_to(Path("/root/autodl-tmp/meep_sim")):
        raise RuntimeError("SERVER_ONLY: expected the authorized server repository")


def main(output=None):
    stage = Path(__file__).resolve().parent
    _assert_server(stage)
    # Meep and MPI are available only at this explicitly server-guarded entry.
    import meep as mp
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank, ranks = comm.rank, comm.size
    master = rank == 0
    if ranks != 4:
        raise RuntimeError("SMOKE_MPI_RANKS_MISMATCH: launch with exactly four MPI ranks")
    mp.verbosity(0)
    output = (stage / "smoke" if output is None else Path(output)).resolve()
    if output != stage / "smoke" and not output.is_relative_to(stage / "smoke"):
        raise ValueError("SMOKE_OUTPUT_OUTSIDE_STAGE_SMOKE")
    if master:
        output.mkdir(parents=True, exist_ok=True)
    comm.Barrier()
    interval_start = time.time()
    started_utc = utcnow()
    model_path = stage / "assets/material/route_a_material.json"
    model = load_model(model_path)
    geometry = dict(period_um=50., surface_fill=.692820323, axis_length_um=40.,
                    tilt_deg=30., thickness_um=60., min_wall_um=1., residual_um=10.)
    geometry_report = validate_geometry(geometry)
    if not geometry_report["valid"]:
        raise ValueError("SMOKE_BASELINE_GEOMETRY_INVALID")
    config = dict(geometry=geometry, pml_um=6., air_above_um=12., air_below_um=12.,
                  resolution=16., courant=.25, wavelength_um=10.5)
    box = layout(config)
    resolution, courant = config["resolution"], config["courant"]
    wavelength = config["wavelength_um"]
    f, dt = 1 / wavelength, courant / resolution
    kx = -f * math.sin(math.radians(30.))
    summary = dict(status=STATUS, physical_qualification=False,
                   smoke_passed=False, started_utc=started_utc,
                   interval_start_epoch=interval_start, mpi_ranks=ranks,
                   material_sha256=digest_file(model_path), geometry=geometry_report,
                   configuration=config, layout=box, polarizations=[], timed_entries=[],
                   ledger_accounting="Launcher must record this process interval in the shared SolverLedger",
                   scope="API, memory layout, arrays, periodic masks and subtraction round trip only")

    def admitted(label, action):
        # Every initialization, zero-step monitor setup, and timed run is gated.
        audit = audit_model(model, resolution, courant, wavelength)
        began = time.time()
        result = action()
        summary["timed_entries"].append(dict(label=label, start_epoch=began,
                                             end_epoch=time.time(), ranks=ranks, q=audit["q"]))
        return result

    try:
        audit_model(model, resolution, courant, wavelength)
        medium = build_medium(model)
        analytical = complex(epsilon(model, f))
        meep_epsilon = np.asarray(medium.epsilon(f), complex)
        target = np.eye(3, dtype=complex) * analytical
        if meep_epsilon.shape != (3, 3) or not np.allclose(meep_epsilon, target, rtol=1e-12, atol=1e-12):
            raise ValueError("SMOKE_MEEP_ANALYTICAL_EPSILON_MISMATCH")
        summary["material_epsilon_check"] = dict(analytical_epsilon=analytical,
                                                meep_epsilon=meep_epsilon,
                                                maximum_absolute_error=float(np.max(abs(meep_epsilon - target))))
        if master:
            atomic_json(output / "smoke_summary.json", summary)

        for pol in ("p", "s"):
            components = {"Ex": mp.Ex, "Hz": mp.Hz} if pol == "p" else {"Ez": mp.Ez, "Hx": mp.Hx}
            edpairs = {"Ex": (mp.Ex, mp.Dx), "Ey": (mp.Ey, mp.Dy)} if pol == "p" else {"Ez": (mp.Ez, mp.Dz)}
            source_component = mp.Hz if pol == "p" else mp.Ez
            pulse = mp.GaussianSource(frequency=f, fwidth=.2 * f, cutoff=5., is_integrated=True)
            source = mp.Source(pulse, source_component, center=mp.Vector3(0, box["source"]),
                               size=mp.Vector3(box["period"], 0),
                               amp_func=lambda r: np.exp(2j * np.pi * kx * r.x))
            chunk_path = output / (pol + "_chunks.h5")
            reference_data = None
            pol_record = {"polarization": pol, "status": STATUS, "simulations": []}
            for case in ("vacuum", "structure"):
                objects = []
                if case == "structure":
                    objects.append(mp.Block(center=mp.Vector3(0, .5 * (box["surface"] + box["bottom"])),
                                            size=mp.Vector3(box["period"], geometry["thickness_um"], mp.inf),
                                            material=medium))
                    for shift in geometry_report["periodic_shifts"]:
                        vertices = [mp.Vector3(x + shift * box["period"], box["surface"] - depth)
                                    for x, depth in geometry_report["vertices_x_depth"]]
                        objects.append(mp.Prism(vertices=vertices, height=mp.inf,
                                                axis=mp.Vector3(0, 0, 1), material=mp.air))
                audit_model(model, resolution, courant, wavelength)
                sim = mp.Simulation(cell_size=mp.Vector3(box["period"], box["height"]),
                                    boundary_layers=[mp.PML(config["pml_um"], direction=mp.Y)],
                                    geometry=objects, sources=[source], resolution=resolution,
                                    Courant=courant, dimensions=2, k_point=mp.Vector3(kx, 0),
                                    force_complex_fields=True, split_chunks_evenly=True,
                                    chunk_layout=str(chunk_path) if case == "structure" else None)
                plane = mp.Volume(center=mp.Vector3(0, box["refl"]), size=mp.Vector3(box["period"], 0))
                flux = sim.add_flux(f, 0, 1, mp.FluxRegion(center=plane.center, size=plane.size, direction=mp.Y))
                plane_monitor = sim.add_dft_fields(list(components.values()), [f], where=plane,
                                                  yee_grid=False, decimation_factor=1)
                slab_volume = mp.Volume(center=mp.Vector3(0, .5 * (box["surface"] + box["bottom"])),
                                        size=mp.Vector3(box["period"], geometry["thickness_um"] + 2 / resolution))
                volume_monitor = sim.add_dft_fields([component for pair in edpairs.values() for component in pair],
                                                   [f], where=slab_volume, yee_grid=True, decimation_factor=1)
                label = pol + "_" + case
                admitted(label + ":init", sim.init_sim)
                admitted(label + ":zero_step_monitor_initialization", lambda: sim.run(until=0))
                if case == "vacuum":
                    # Meep owns this collective HDF5 write, like the production solver.
                    sim.dump_chunk_layout(str(chunk_path))
                else:
                    sim.load_minus_flux_data(flux, reference_data)
                    # Check the rank-local subtraction buffers before adding any fields.
                    loaded = sim.get_flux_data(flux)
                    for name in ("E", "H"):
                        if not np.array_equal(np.asarray(getattr(loaded, name)), -np.asarray(getattr(reference_data, name))):
                            raise ValueError("SMOKE_FLUX_SUBTRACTION_ROUNDTRIP_MISMATCH:" + name)
                start_time = float(sim.meep_time())
                admitted(label + ":two_timesteps", lambda: sim.run(until=2 * dt))
                end_time = float(sim.meep_time())
                if not np.isclose(end_time - start_time, 2 * dt, rtol=0, atol=1e-12):
                    raise ValueError("SMOKE_TIMESTEP_COUNT_MISMATCH")
                plane_x, plane_values, raw_plane = plane_fields(sim, plane_monitor, components,
                                                                box["period"], kx, resolution)
                pair_records, raw_volume, grid_metadata = native_volume_pairs(
                    sim, volume_monitor, edpairs, resolution, box["period"], kx,
                    geometry, box["surface"])
                flux_data = sim.get_flux_data(flux)
                rank_buffers = comm.gather({"E": np.array(flux_data.E, copy=True),
                                            "H": np.array(flux_data.H, copy=True)}, root=0)
                record = dict(case=case, status=STATUS, simulation_start_time=start_time,
                              simulation_end_time=end_time, timestep_count=2,
                              flux=float(mp.get_fluxes(flux)[0]), plane_shape=list(plane_x.shape),
                              plane_components={name: list(value.shape) for name, value in plane_values.items()},
                              volume_grid=grid_metadata,
                              subtraction_loaded=(case == "structure"),
                              mask_interpretation="Baseline Ti geometry mask; vacuum case exercises coordinate API only")
                if master:
                    atomic_npz(output / (label + "_planes.npz"), **raw_plane)
                    atomic_npz(output / (label + "_volume.npz"), **raw_volume)
                    for buffer_rank, buffer in enumerate(rank_buffers):
                        atomic_npz(output / f"{label}_flux_rank{buffer_rank:03d}.npz", **buffer)
                    record["flux_buffer_shapes_by_rank"] = [
                        {name: list(value.shape) for name, value in buffer.items()} for buffer in rank_buffers]
                    atomic_json(output / (label + "_probe.json"), record)
                if case == "vacuum":
                    # Own independent arrays after resetting the reference simulation.
                    reference_data = type(flux_data)(E=np.array(flux_data.E, copy=True),
                                                      H=np.array(flux_data.H, copy=True))
                pol_record["simulations"].append(record)
                sim.reset_meep()
                comm.Barrier()
            if master:
                pol_record["chunk_layout_sha256"] = digest_file(chunk_path)
            summary["polarizations"].append(pol_record)
            if master:
                atomic_json(output / "smoke_summary.json", summary)
        summary["smoke_passed"] = True
        summary["simulation_count"] = 4
        summary["total_timesteps"] = 8
        summary["status"] = STATUS
        summary["interval_end_epoch"] = time.time()
        summary["completed_utc"] = utcnow()
        summary["wall_seconds"] = summary["interval_end_epoch"] - interval_start
        summary["rank_hours"] = ranks * summary["wall_seconds"] / 3600
        if master:
            atomic_json(output / "smoke_summary.json", summary)
            print(f"{STATUS}: API smoke passed; 4 simulations, 8 timesteps; no physical qualification.")
        return 0
    except Exception:
        summary.update(smoke_passed=False, status=STATUS, failure=traceback.format_exc(),
                       interval_end_epoch=time.time(), failed_rank=rank)
        if master:
            atomic_json(output / "smoke_summary.json", summary)
        # A rank-local failure must not leave other ranks stuck in collective I/O.
        comm.Abort(2)
        return 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output))
