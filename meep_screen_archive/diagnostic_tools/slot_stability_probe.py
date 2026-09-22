"""Isolated, server-only short-window slot diagnostic; never qualification.

The two arms differ only in ``eps_averaging``. Both omit passive production
DFT/flux/volume monitors and incident-reference preload. Thus the averaging-on
arm must reproduce rapid growth before an arm comparison supports causality.
No production result, check, cache, queue, or solver source is modified.

Importing this module and --inspect-only never import Meep, MPI, or SciPy.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import signal
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from ti2d.common import atomic_json, atomic_npz, digest_file, digest_object, implementation_digest, utcnow
from ti2d.geometry import periodic_slot_mask, validate_geometry


FROZEN_CORE_SHA256 = "24d09c50dcf6c0cb9271947dbcb04f03b62d9d0500bd7bace34e6909d54834c1"
FROZEN_CASE_SHA256 = "3e6b9631a47a3c19040e74fd3c39cfef8df228f45e7bf766782c3e1fb414a85c"
SERVER_STAGE = Path("/root/autodl-tmp/meep_sim/meep_screen")
REFERENCE_INTENSITY_PEAK = 0.3334312033062042
T_STOP = 100.0
SAMPLE_INTERVAL = 2.0
PATCH_WIDTH = 17
MAX_PATCHES = 4
ARMS = {"averaging_on": True, "averaging_off": False}
FINAL_STATUSES = {"EARLY_GROWTH_DETECTED", "NONFINITE_FIELD",
                  "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION", "TIME_LIMIT_UNQUALIFIED",
                  "INTERRUPTED", "CRASHED"}
MAX_LOG10_FLOAT = math.log10(np.finfo(float).max)


class CollectiveDiagnosticError(RuntimeError):
    """A Python-side failure observed and reported consistently by every rank."""


def collective_call(comm, mpi, name, function):
    """Synchronize Python failures after a same-order operation on every rank.

    Native MPI/Meep crashes cannot be recovered here; the external runner must
    retain its timeout and process-group cleanup. In particular this does not
    claim to rescue a rank that crashes inside a native collective operation.
    """
    value, failure = None, None
    try:
        value = function()
    except Exception as exc:
        failure = {"rank": int(comm.rank), "operation": name,
                   "error": f"{type(exc).__name__}: {exc}"}
    failed = comm.allreduce(int(failure is not None), op=mpi.MAX)
    if failed:
        errors = [entry for entry in comm.allgather(failure) if entry is not None]
        raise CollectiveDiagnosticError(json.dumps(errors, sort_keys=True))
    return value


def scalar_record(value):
    """JSON-safe exact finite values, with explicit nonfinite representations."""
    value = float(value)
    if math.isfinite(value):
        return value
    return {"value": None, "nonfinite": "nan" if math.isnan(value) else "inf" if value > 0 else "-inf"}


def complex_record(value):
    value = complex(value)
    return {"real": scalar_record(value.real), "imag": scalar_record(value.imag)}


def representable_power(log10_power):
    """Never square a large field or hide overflow by clipping it."""
    if log10_power is None:
        return None, "nonfinite input field"
    if log10_power == -math.inf:
        return 0.0, "exactly zero field"
    if log10_power > MAX_LOG10_FLOAT:
        return None, "intensity exceeds binary64; finite log10_intensity retained"
    try:
        power = 10.0 ** log10_power
    except OverflowError:
        return None, "intensity exceeds binary64; finite log10_intensity retained"
    return power, "underflows binary64; log10_intensity retained" if power == 0 else None


def log_power_arrays(ex, ey):
    """Compute log10(|Ex|²+|Ey|²) by scaling the four real components.

    Finite flags are taken from original fields before magnitude arithmetic.
    Only dimensionless values bounded by one are squared. Zero is represented
    internally by -inf, then encoded as null plus an explicit zero flag in JSON.
    """
    ex, ey = np.asarray(ex), np.asarray(ey)
    if ex.ndim != 2 or ex.shape != ey.shape or not ex.size:
        raise ValueError(f"FIELD_SHAPE_MISMATCH: Ex={ex.shape}, Ey={ey.shape}; no cropping allowed")
    if ex.dtype.kind not in "fc" or ey.dtype.kind not in "fc":
        raise ValueError("FIELD_DTYPE_INVALID")
    finite = np.isfinite(ex) & np.isfinite(ey)
    output = np.full(ex.shape, np.nan, dtype=float)
    parts = [ex.real, ex.imag, ey.real, ey.imag]
    # Avoid stacking four full-grid absolute arrays on every MPI rank.
    scale = np.abs(parts[0])
    for part in parts[1:]:
        np.maximum(scale, np.abs(part), out=scale)
    zero = finite & (scale == 0)
    valid = finite & (scale > 0)
    output[zero] = -math.inf
    if np.any(valid):
        normalized_power = sum((part[valid] / scale[valid]) ** 2 for part in parts)
        output[valid] = 2 * np.log10(scale[valid]) + np.log10(normalized_power)
    return output, finite


def checked_grid(ex, ey, metadata, resolution):
    """Match the full centered arrays to metadata without cropping or padding."""
    if len(metadata) != 4:
        raise ValueError("GRID_METADATA_COUNT_INVALID")
    x, y, z, weights = [np.asarray(value) for value in metadata]
    if x.ndim != 1 or y.ndim != 1 or z.shape != (1,) or min(x.size, y.size) < 2:
        raise ValueError("GRID_COORDINATE_SHAPE_INVALID")
    shape = (x.size, y.size)
    if np.shape(ex) != shape or np.shape(ey) != shape or weights.shape != shape:
        raise ValueError(f"GRID_FIELD_SHAPE_MISMATCH: expected {shape}; Ex {np.shape(ex)}; Ey {np.shape(ey)}; weights {weights.shape}")
    if any(value.dtype.kind not in "fiu" or not np.isfinite(value).all() for value in (x, y, z, weights)):
        raise ValueError("GRID_NONFINITE_OR_NONREAL_METADATA")
    if not all(np.allclose(np.diff(axis), 1 / resolution, atol=1e-9, rtol=1e-9) for axis in (x, y)):
        raise ValueError("GRID_SPACING_MISMATCH")
    if np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError("GRID_WEIGHTS_INVALID")
    return x, y, z, weights


def point_region(x, y, geometry, surface):
    depth = surface - y
    slot = bool(periodic_slot_mask(x, depth, geometry))
    slab = bool(0 <= depth <= geometry["thickness_um"])
    return {"x_um": float(x), "y_um": float(y), "depth_um": float(depth),
            "ideal_geometry_region": "slot_air" if slot and slab else "titanium" if slab else "exterior_air",
            "inside_ideal_slot": slot, "inside_slab_extent": slab,
            "classification_note": "analytic geometry, not an effective-material tensor or a Yee-cell volume fraction"}


def field_summary(ex, ey, x, y, geometry, surface):
    log_power, finite = log_power_arrays(ex, ey)
    all_finite = bool(finite.all())
    if all_finite:
        index = np.unravel_index(int(np.argmax(log_power)), log_power.shape)
        maximum = float(log_power[index])
        zero = maximum == -math.inf
        intensity, note = representable_power(maximum)
        log_json = None if zero else maximum
        amplitude = representable_power(maximum / 2)[0]
    else:
        index = tuple(int(value) for value in np.argwhere(~finite)[0])
        maximum, zero, intensity, amplitude, log_json = None, False, None, None, None
        note = "nonfinite input field; location is first nonfinite sample, not a finite peak"
    ix, iy = [int(value) for value in index]
    record = {"all_fields_finite": all_finite, "nonfinite_grid_points": int((~finite).sum()),
              "envelope_amplitude": amplitude, "envelope_intensity": intensity,
              "log10_envelope_intensity": log_json, "zero_field": zero,
              "intensity_encoding_note": note,
              "log10_intensity_over_reference_peak": (log_json - math.log10(REFERENCE_INTENSITY_PEAK)
                                                       if log_json is not None else None),
              "index_xy": [ix, iy], "raw_Ex_at_location": complex_record(ex[ix, iy]),
              "raw_Ey_at_location": complex_record(ey[ix, iy]),
              "location": point_region(x[ix], y[iy], geometry, surface)}
    return record, index


def point_fields_summary(ex, ey):
    log_power, finite = log_power_arrays(np.array([[ex]], complex), np.array([[ey]], complex))
    logarithm = float(log_power[0, 0]) if finite[0, 0] else None
    power, note = representable_power(logarithm)
    return {"raw_Ex": complex_record(ex), "raw_Ey": complex_record(ey),
            "finite": bool(finite[0, 0]), "intensity": power,
            "log10_intensity": logarithm if logarithm is not None and math.isfinite(logarithm) else None,
            "zero_field": logarithm == -math.inf, "intensity_encoding_note": note}


class GrowthGuard:
    """Require three adjacent above-threshold samples and both >10x jumps."""

    def __init__(self):
        self.history = []

    def sample(self, summary):
        if not summary["all_fields_finite"]:
            return "NONFINITE_FIELD"
        ratio = summary["log10_intensity_over_reference_peak"]
        self.history.append(ratio)
        self.history = self.history[-3:]
        if (len(self.history) == 3 and all(item is not None and item > 12 for item in self.history)
                and all(b - a > 1 for a, b in zip(self.history, self.history[1:]))):
            return "EARLY_GROWTH_DETECTED"
        return None


def local_patch(ex, ey, x, y, index, geometry, surface):
    """At most 17x17 samples; raw complex values, exact coordinates and masks."""
    bounds = []
    for center, extent in zip(index, np.shape(ex)):
        start = max(0, min(int(center) - PATCH_WIDTH // 2, extent - PATCH_WIDTH))
        bounds.append(slice(start, min(extent, start + PATCH_WIDTH)))
    xs, ys = bounds
    xx, yy = np.meshgrid(x[xs], y[ys], indexing="ij")
    depth = surface - yy
    slab = (depth >= 0) & (depth <= geometry["thickness_um"])
    slot = periodic_slot_mask(xx, depth, geometry)
    return {"Ex_raw": np.asarray(ex)[xs, ys], "Ey_raw": np.asarray(ey)[xs, ys],
            "x_um": np.asarray(x)[xs], "y_um": np.asarray(y)[ys], "depth_um": depth,
            "ideal_slab_mask": slab, "ideal_slot_mask": slot, "ideal_ti_mask": slab & ~slot,
            "global_peak_index_xy": np.asarray(index, dtype=int),
            "patch_start_index_xy": np.array([xs.start, ys.start], dtype=int)}


def load_inspection(plan_path, case_id, output, *, server=False):
    """Validate the immutable plan, source bytes, production core and arm."""
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    if plan.get("schema_version") != 1:
        raise ValueError("PLAN_SCHEMA_INVALID")
    supplied_hash = plan.get("plan_sha256")
    if supplied_hash != digest_object({key: value for key, value in plan.items() if key != "plan_sha256"}):
        raise ValueError("PLAN_HASH_MISMATCH")
    stage = Path(plan["stage_root"]).resolve()
    run_dir = Path(plan["run_dir"]).resolve()
    if stage != Path(__file__).resolve().parents[1]:
        raise ValueError("PLAN_STAGE_ROOT_MISMATCH")
    if run_dir.parent != stage / "repair_diagnostics" or run_dir.name != plan["run_id"]:
        raise ValueError("DIAGNOSTIC_RUN_ROOT_INVALID")
    if output != run_dir / case_id or plan_path != run_dir / "plan.json":
        raise ValueError("DIAGNOSTIC_OUTPUT_OR_PLAN_PATH_INVALID")
    if server and (sys.platform != "linux" or stage != SERVER_STAGE):
        raise ValueError("MEEP_SERVER_ONLY: local FDTD is prohibited")
    if plan.get("cases") != [{"id": key, "eps_averaging": value} for key, value in ARMS.items()]:
        raise ValueError("DIAGNOSTIC_ARMS_CHANGED")
    if case_id not in ARMS:
        raise ValueError("DIAGNOSTIC_ARM_INVALID")
    diagnostic = plan["diagnostic"]
    expected = {"stop_time": T_STOP, "sample_interval": SAMPLE_INTERVAL,
                "reference_intensity_peak": REFERENCE_INTENSITY_PEAK,
                "growth_ratio_limit": 1e12, "growth_step_factor": 10.0,
                "growth_consecutive_samples": 3, "patch_radius_cells": 8}
    if diagnostic != expected:
        raise ValueError("DIAGNOSTIC_LIMITS_CHANGED")
    if plan.get("implementation_sha256") != FROZEN_CORE_SHA256 or implementation_digest() != FROZEN_CORE_SHA256:
        raise ValueError("FROZEN_PRODUCTION_CORE_MISMATCH")
    if plan.get("probe_sha256") != digest_file(__file__):
        raise ValueError("DIAGNOSTIC_PROBE_HASH_MISMATCH")
    case_path = Path(plan["source_case_path"]).resolve()
    if plan.get("source_case_sha256") != FROZEN_CASE_SHA256 or digest_file(case_path) != FROZEN_CASE_SHA256:
        raise ValueError("FROZEN_PHYSICAL_CASE_BYTES_MISMATCH")
    case = json.loads(case_path.read_text())
    if (case.get("implementation_sha256") != FROZEN_CORE_SHA256 or case.get("resolution") != 32
            or case.get("courant") != .25 or case.get("polarization") != "p"
            or case.get("emission_theta_deg") != 30 or case.get("case") != "structure"
            or case.get("mpi_ranks") != 4 or case.get("oxide_layer_in_model") is not False):
        raise ValueError("FROZEN_PHYSICAL_CASE_PARAMETERS_MISMATCH")
    if (plan.get("material_sha256") != case.get("material_sha256")
            or digest_file(case["material_path"]) != case["material_sha256"]):
        raise ValueError("FROZEN_ORDAL_MATERIAL_MISMATCH")
    if plan["budgets"].get("mpi_ranks") != 4:
        raise ValueError("MPI_RANKS_MISMATCH")
    if server and case_path != SERVER_STAGE / "runs/baseline_r32_20260922_E/cases/baseline_r32_p_plus30.json":
        raise ValueError("FROZEN_SOURCE_CASE_PATH_MISMATCH")
    if any((output / name).exists() for name in ("result.json", "trace.json", "heartbeat.json")):
        raise ValueError("DIAGNOSTIC_RETRY_OR_OVERWRITE_PROHIBITED")
    report = {"schema_version": 1, "case_id": case_id, "eps_averaging": ARMS[case_id],
              "qualification": False, "fdtd_started": False,
              "physical_case_sha256": FROZEN_CASE_SHA256,
              "source_case_sha256": FROZEN_CASE_SHA256,
              "implementation_sha256": FROZEN_CORE_SHA256,
              "material_sha256": plan["material_sha256"],
              "probe_sha256": digest_file(__file__), "plan_sha256": supplied_hash,
              "case_config": case, "diagnostic": diagnostic, "output": str(output),
              "only_simulation_parameter_difference": "eps_averaging",
              "passive_monitors_omitted": ["DFT", "flux", "volume_absorption", "reference_preload"],
              "causal_interpretation_requires_averaging_on_growth_reproduction": True,
              "scope": "short-window field stability diagnostic; no production qualification or checks"}
    return plan, case, report


def _probe_coordinates(case, ly):
    probes = {"front_air": (0.0, ly["refl"]), "back_air": (0.0, ly["trans"])}
    angle = math.radians(case["geometry"]["tilt_deg"])
    for fraction in (.25, .50, .75):
        distance = case["geometry"]["axis_length_um"] * fraction
        px = (distance * math.sin(angle) + ly["period"] / 2) % ly["period"] - ly["period"] / 2
        probes["cavity_" + str(fraction)] = (px, ly["surface"] - distance * math.cos(angle))
    return probes


def main(plan_path, case_id, output):
    """The only numerical entry point: exclusively on the approved server."""
    stage = Path(__file__).resolve().parents[1]
    if sys.platform != "linux" or stage != SERVER_STAGE:
        raise ValueError("MEEP_SERVER_ONLY: local FDTD is prohibited")
    # No import of either runtime is allowed above this server-only guard.
    import meep as mp
    from mpi4py import MPI
    from ti2d.flux_reference import make_partition, partition_spec
    from ti2d.materials import audit_model, build_medium, load_model
    from ti2d.solver import layout

    comm = MPI.COMM_WORLD
    mp.verbosity(0)
    master = comm.rank == 0
    started = time.monotonic()
    terminated = [False]
    signal.signal(signal.SIGTERM, lambda *_: terminated.__setitem__(0, True))
    signal.signal(signal.SIGINT, lambda *_: terminated.__setitem__(0, True))
    out = Path(output).resolve()
    base, trace, patches, active_sim = None, [], [], None
    status, failure, t_reached, source_end = "CRASHED", None, 0.0, None
    plan, case, inspection = collective_call(comm, MPI, "validate_plan", lambda: load_inspection(plan_path, case_id, output, server=True))
    if comm.size != case["mpi_ranks"]:
        raise ValueError("MPI_RANKS_MISMATCH")
    wall_limit = min(float(os.environ["TI2D_DIAGNOSTIC_TIMEOUT_SECONDS"]), float(plan["budgets"]["single_seconds"]))
    if not math.isfinite(wall_limit) or wall_limit <= 0:
        raise ValueError("DIAGNOSTIC_TIMEOUT_INVALID")
    collective_call(comm, MPI, "create_output", lambda: out.mkdir(parents=True, exist_ok=True) if master else None)
    base = {key: value for key, value in inspection.items() if key != "fdtd_started"}
    base.update(started_utc=utcnow(), mpi_ranks=int(comm.size), qualification=False,
                no_dft=True, no_reference_run=True, no_reference_preload=True,
                production_cache_read=False, production_cache_write=False,
                production_checks_modified=False, timeout_seconds=wall_limit,
                source_end_time=None, t_reached=0.0,
                numeric_encoding="Original complex fields are retained at peaks/probes and in small patches. No clipping. Null intensity beyond binary64 retains finite log10 intensity.",
                growth_guard_semantics="Three adjacent samples must each have intensity/reference_peak > 1e12, and both adjacent intensity growth factors must be > 10.")
    collective_call(comm, MPI, "write_initial_heartbeat", lambda: atomic_json(out / "heartbeat.json", {
        "phase": "INITIALIZING", "case_id": case_id, "qualification": False,
        "eps_averaging": ARMS[case_id], "t_reached": 0.0, "wall_elapsed_seconds": time.monotonic() - started,
        "updated_utc": utcnow()}) if master else None)
    try:
        model = collective_call(comm, MPI, "load_material", lambda: load_model(case["material_path"]))
        stability = collective_call(comm, MPI, "audit_material", lambda: audit_model(model, case["resolution"], case["courant"], case["wavelength_um"]))
        geometry_report = validate_geometry(case["geometry"])
        if not geometry_report["valid"]:
            raise ValueError("GEOMETRY_INVALID")
        ly = layout(case)
        f = 1 / case["wavelength_um"]
        theta = math.radians(case["emission_theta_deg"])
        kx = -f * math.sin(theta)
        pulse = mp.GaussianSource(frequency=f, fwidth=case["source_fwidth_fraction"] * f,
                                  cutoff=case["source_cutoff"], is_integrated=True)
        source_end = float(pulse.start_time + 2 * pulse.width * pulse.cutoff)
        sources = [mp.Source(pulse, mp.Hz, center=mp.Vector3(0, ly["source"]),
                             size=mp.Vector3(ly["period"], 0),
                             amp_func=lambda r: np.exp(2j * np.pi * kx * r.x))]
        medium = build_medium(model)
        objects = [mp.Block(center=mp.Vector3(0, .5 * (ly["surface"] + ly["bottom"])),
                            size=mp.Vector3(ly["period"], case["geometry"]["thickness_um"], mp.inf), material=medium)]
        for shift in geometry_report["periodic_shifts"]:
            vertices = [mp.Vector3(x + shift * ly["period"], ly["surface"] - depth)
                        for x, depth in geometry_report["vertices_x_depth"]]
            objects.append(mp.Prism(vertices=vertices, height=mp.inf, axis=mp.Vector3(0, 0, 1), material=mp.air))
        active_sim = collective_call(comm, MPI, "create_simulation", lambda: mp.Simulation(
            cell_size=mp.Vector3(ly["period"], ly["height"]),
            boundary_layers=[mp.PML(case["pml_um"], direction=mp.Y)], geometry=objects, sources=sources,
            resolution=case["resolution"], Courant=case["courant"], dimensions=2, k_point=mp.Vector3(kx, 0),
            force_complex_fields=True, split_chunks_evenly=True, chunk_layout=make_partition(mp, comm.size),
            eps_averaging=ARMS[case_id]))
        probes = _probe_coordinates(case, ly)
        probe_vectors = {key: mp.Vector3(*xy) for key, xy in probes.items()}
        non_pml = mp.Volume(center=mp.Vector3(), size=mp.Vector3(ly["period"], ly["non_pml_height"]))
        base.update(source_end_time=source_end, layout=ly, geometry_audit=geometry_report,
                    stability_rough_screen=stability, fixed_mpi_partition=partition_spec(comm.size),
                    probe_coordinates_xy=probes,
                    reciprocal_incident_k_simulation_xy=[kx, -f * math.cos(theta), 0.0],
                    physical_to_simulation_axes="physical (x,z) -> Meep (x,y)",
                    field_grid="centered get_array Ex/Ey as used by production non-PML envelope; not native staggered arrays",
                    passive_omission_note="DFT, flux and volume monitors and reference preload are passive. These omissions still require averaging-on reproduction of rapid growth for causal comparison.")
        collective_call(comm, MPI, "initialize_simulation", active_sim.init_sim)
        collective_call(comm, MPI, "initialize_zero_time", lambda: active_sim.run(until=0))
        metadata = collective_call(comm, MPI, "get_centered_metadata", lambda: active_sim.get_array_metadata(vol=non_pml))
        last_sample = [-1.0]
        reason = [None]
        guard = GrowthGuard()
        grid_saved = [False]
        patch_levels = set()

        def collect(sim):
            nonlocal t_reached
            now = float(sim.meep_time())
            if now <= last_sample[0] + 1e-8:
                return
            # Every rank calls these in the same sequence, even if local
            # validation will fail. Shape validation has its own allreduce.
            ex = collective_call(comm, MPI, "get_Ex", lambda: sim.get_array(component=mp.Ex, vol=non_pml, cmplx=True))
            ey = collective_call(comm, MPI, "get_Ey", lambda: sim.get_array(component=mp.Ey, vol=non_pml, cmplx=True))
            x, y, _, weights = collective_call(comm, MPI, "validate_centered_grid", lambda: checked_grid(ex, ey, metadata, case["resolution"]))
            summary, peak_index = collective_call(comm, MPI, "summarize_envelope", lambda: field_summary(ex, ey, x, y, case["geometry"], ly["surface"]))
            sampled_probes = {}
            for name, vector in probe_vectors.items():
                point_ex = collective_call(comm, MPI, "get_probe_Ex_" + name, lambda vector=vector: sim.get_field_point(mp.Ex, vector))
                point_ey = collective_call(comm, MPI, "get_probe_Ey_" + name, lambda vector=vector: sim.get_field_point(mp.Ey, vector))
                sampled_probes[name] = collective_call(comm, MPI, "summarize_probe_" + name, lambda: point_fields_summary(point_ex, point_ey))
            local_reason = guard.sample(summary)
            if not all(item["finite"] for item in sampled_probes.values()):
                local_reason = "NONFINITE_FIELD"
            reason_code = comm.allreduce({None: 0, "EARLY_GROWTH_DETECTED": 1, "NONFINITE_FIELD": 2}[local_reason], op=MPI.MAX)
            reason[0] = {0: reason[0], 1: "EARLY_GROWTH_DETECTED", 2: "NONFINITE_FIELD"}[reason_code]
            sample = {"simulation_time": now, "wall_elapsed_seconds": time.monotonic() - started,
                      "envelope": summary, "probes": sampled_probes,
                      "early_guard_status": reason[0]}
            trace.append(sample)
            t_reached, last_sample[0] = now, now

            def persist():
                if not master:
                    return
                if not grid_saved[0]:
                    atomic_npz(out / "sample_grid.npz", x_um=x, y_um=y, z_um=np.asarray(metadata[2]))
                    atomic_json(out / "geometry_metadata.json", {
                        "geometry": case["geometry"], "geometry_audit": geometry_report, "layout": ly,
                        "array_shape": list(np.shape(ex)), "metadata_weight_shape": list(weights.shape),
                        "metadata_weight_sum": float(np.sum(weights)),
                        "coordinates_source": "get_array_metadata(vol=non_pml)",
                        "shape_policy": "exact Ex/Ey/weights/coordinate match; no cropping",
                        "ideal_mask_interpretation": "analytic region only, not averaged dielectric-tensor or volume-fraction evidence"})
                    grid_saved[0] = True
                ratio = summary["log10_intensity_over_reference_peak"]
                terminal = reason[0] is not None or now >= T_STOP or terminated[0]
                level = ("terminal" if terminal else "ratio_1e12" if ratio is not None and ratio > 12
                         else "ratio_1e6" if ratio is not None and ratio > 6 else "first_sample" if now >= SAMPLE_INTERVAL else None)
                if level is not None and level not in patch_levels and len(patches) < MAX_PATCHES:
                    patch_name = f"peak_patch_{len(patches):02d}.npz"
                    atomic_npz(out / patch_name, simulation_time=np.array(now),
                               **local_patch(ex, ey, x, y, peak_index, case["geometry"], ly["surface"]))
                    patches.append({"file": patch_name, "simulation_time": now, "trigger": level,
                                    "sha256": digest_file(out / patch_name)})
                    patch_levels.add(level)
                atomic_json(out / "trace.json", {"schema_version": 1, "qualification": False,
                            "case_id": case_id, "eps_averaging": ARMS[case_id], "samples": trace,
                            "reference_intensity_peak": REFERENCE_INTENSITY_PEAK, "patches": patches})
                atomic_json(out / "heartbeat.json", {"phase": "SAMPLING", "case_id": case_id,
                            "eps_averaging": ARMS[case_id], "qualification": False,
                            "t_reached": now, "source_end_time": source_end, "early_guard_status": reason[0],
                            "wall_elapsed_seconds": time.monotonic() - started,
                            "updated_utc": utcnow(), "last_sample": sample})
            collective_call(comm, MPI, "persist_sample", persist)

        def stop(sim):
            interruption = comm.allreduce(int(terminated[0]), op=MPI.MAX)
            timeout = comm.allreduce(int(time.monotonic() - started >= wall_limit), op=MPI.MAX)
            if reason[0] is not None:
                return True
            if interruption:
                reason[0] = "INTERRUPTED"
            elif timeout:
                reason[0] = "TIME_LIMIT_UNQUALIFIED"
            elif float(sim.meep_time()) >= T_STOP:
                reason[0] = "SHORT_WINDOW_COMPLETED_NOT_QUALIFICATION"
            return reason[0] is not None

        collect(active_sim)
        active_sim.run(mp.at_every(SAMPLE_INTERVAL, collect), until=stop)
        collect(active_sim)
        status = reason[0] or "CRASHED"
    except Exception as exc:
        failure = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
        status = "CRASHED"
        if active_sim is not None:
            try:
                t_reached = float(active_sim.meep_time())
            except Exception:
                pass
    # Sampling/setup Python exceptions are synchronized by collective_call.
    # Native process failures are classified and sealed by the external runner.
    def finalize():
        if not master:
            return
        evidence = {path.name: digest_file(path) for path in sorted(out.iterdir())
                    if path.is_file() and (path.name in {"trace.json", "geometry_metadata.json", "sample_grid.npz"}
                                           or path.name.startswith("peak_patch_"))}
        summary = {"sample_count": len(trace), "first_time": trace[0]["simulation_time"] if trace else None,
                   "last_time": trace[-1]["simulation_time"] if trace else None,
                   "max_log10_intensity": max((sample["envelope"]["log10_envelope_intensity"]
                                              for sample in trace if sample["envelope"]["log10_envelope_intensity"] is not None), default=None),
                   "last_sample": trace[-1] if trace else None}
        result = dict(base, status=status, qualification=False, finished_utc=utcnow(),
                      t_reached=t_reached, source_end_time=source_end,
                      elapsed_seconds=time.monotonic() - started, trace_summary=summary,
                      evidence_sha256=evidence, patches=patches, error=failure,
                      interpretation="A completed t<=100 window is not source-end decay or stability qualification. Averaging-on must reproduce growth before drawing causal conclusions from the comparison.")
        atomic_json(out / "result.json", result)
        atomic_json(out / "heartbeat.json", {"phase": "FINISHED", "status": status, "case_id": case_id,
                    "eps_averaging": ARMS[case_id], "qualification": False, "t_reached": t_reached,
                    "wall_elapsed_seconds": time.monotonic() - started, "updated_utc": utcnow(),
                    "result_sha256": digest_file(out / "result.json")})
    collective_call(comm, MPI, "finalize_evidence", finalize)
    return 2 if status in {"CRASHED", "INTERRUPTED", "TIME_LIMIT_UNQUALIFIED"} else 0


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--case-id", required=True, choices=ARMS)
    parser.add_argument("--output", required=True)
    parser.add_argument("--inspect-only", action="store_true", help="Validate configuration/hashes only; no MPI or Meep initialization")
    args = parser.parse_args()
    if args.inspect_only:
        _, _, report = load_inspection(args.plan, args.case_id, args.output)
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
        return 0
    return main(args.plan, args.case_id, args.output)


if __name__ == "__main__":
    sys.exit(cli())
