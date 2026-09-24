"""Deterministic rank ownership and readback checks for cached incident fields.

This module imports no Meep or MPI runtime. Meep's chunk HDF5 file records grid
volumes, but does not preserve MPI owners; both simulations must receive the
same explicit BinaryPartition and compare the resulting layout signatures.
See https://meep.readthedocs.io/en/latest/Parallel_Meep/#user-specified-cell-partition
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

from .common import digest_object


def collective_layout_signature(sim, path, comm):
    """Broadcast read failures too, so peers never hang awaiting rank zero."""
    payload = None
    if comm.rank == 0:
        try:
            payload = {'signature': layout_signature(path, sim.structure.get_chunk_owners()), 'error': None}
        except Exception as exc:
            payload = {'signature': None, 'error': f'{type(exc).__name__}: {exc}'}
    payload = comm.bcast(payload, root=0)
    if payload['error']:
        raise ValueError('CHUNK_LAYOUT_READ_FAILED: ' + payload['error'])
    return payload['signature']


def partition_spec(ranks):
    """Return the JSON-portable partition used by the qualified four-rank path."""
    if isinstance(ranks, (bool, np.bool_)) or not isinstance(ranks, (int, np.integer)) or ranks != 4:
        raise ValueError("REFERENCE_PARTITION_REQUIRES_FOUR_RANKS")
    return [["Y", 0.0], [["X", 0.0], 0, 1], [["X", 0.0], 2, 3]]


def make_partition(mp, ranks):
    """Construct a fresh Meep partition with explicit leaf process IDs."""
    def convert(node):
        if isinstance(node, int):
            return node
        (axis, position), left, right = node
        return [(getattr(mp, axis), position), convert(left), convert(right)]

    return mp.BinaryPartition(data=convert(partition_spec(ranks)))


def layout_signature(chunk_file, owners):
    """Hash the logical chunk coordinates/sizes AND ordered MPI owners.

    HDF5 encoding/compression/byte order are deliberately excluded. Meep writes
    one flattened triple per chunk into each of gv_nums and gv_origins.
    """
    import h5py

    owners = np.asarray(owners)
    if owners.ndim != 1 or owners.size == 0 or owners.dtype.kind not in "iu" or np.any(owners < 0):
        raise ValueError("REFERENCE_CHUNK_OWNERS_INVALID")
    datasets = {}
    with h5py.File(str(chunk_file), "r") as source:
        for name in ("gv_nums", "gv_origins"):
            if name not in source:
                raise ValueError(f"REFERENCE_CHUNK_DATASET_MISSING:{name}")
            values = np.asarray(source[name][()])
            if (values.ndim != 1 or values.size != 3 * owners.size or
                    values.dtype.kind not in "iuf" or not np.isfinite(values).all()):
                raise ValueError(f"REFERENCE_CHUNK_DATASET_INVALID:{name}")
            if name == "gv_nums":
                if np.any(values < 0) or np.any(values != np.floor(values)):
                    raise ValueError("REFERENCE_CHUNK_GRID_COUNTS_INVALID")
                logical = [int(value) for value in values]
            else:
                logical = [float(value) for value in values]
            datasets[name] = {"shape": list(values.shape), "values": logical}
    payload = {"schema_version": 1, "datasets": datasets,
               "owners": [int(value) for value in owners]}
    return dict(payload, sha256=digest_object(payload))


def _buffers(data):
    if isinstance(data, Mapping):
        return {name: np.asarray(data[name]) for name in ("E", "H")}
    if hasattr(data, "E") and hasattr(data, "H"):
        return {name: np.asarray(getattr(data, name)) for name in ("E", "H")}
    if isinstance(data, (tuple, list)) and len(data) == 2:
        return {name: np.asarray(value) for name, value in zip(("E", "H"), data)}
    raise ValueError("expected E/H arrays or a FluxData pair")


def _array_stats(array):
    numeric = array.dtype.kind in "iufc"
    finite = bool(numeric and np.isfinite(array).all())
    nonzero = int(np.count_nonzero(array)) if numeric else None
    with np.errstate(over="ignore", invalid="ignore"):
        magnitudes = np.abs(array.astype(np.complex128)) if finite else None
        maximum = float(np.max(magnitudes)) if finite and array.size else 0.0 if finite else None
        power = float(np.sum(magnitudes ** 2, dtype=float)) if finite else None
    if maximum is not None and not math.isfinite(maximum):
        maximum = None
    if power is not None and not math.isfinite(power):
        power = None
    return {"shape": list(array.shape), "finite": finite, "nonzero": nonzero,
            "max_abs": maximum, "power": power}


def buffer_check(expected, loaded):
    """Compare same-sign rank-local E/H buffers, retaining zero-owner ranks.

    The per-component tolerance is atol=1e-12*max(abs(expected),abs(loaded))
    plus rtol=1e-12*abs(expected). No dimensionful floor can hide a tiny signal.
    Local zero buffers are allowed; callers must enforce global nonzero input.
    """
    report = {"valid": False, "shape_match": False, "finite": False,
              "expected_nonzero": 0, "loaded_nonzero": 0,
              "expected_power": None, "loaded_power": None,
              "components": {}, "errors": []}
    try:
        expected, loaded = _buffers(expected), _buffers(loaded)
    except (KeyError, TypeError, ValueError) as exc:
        report["errors"].append(f"BUFFER_INPUT_INVALID:{exc}")
        return report
    for name in ("E", "H"):
        wanted, actual = expected[name], loaded[name]
        a, b = _array_stats(wanted), _array_stats(actual)
        same_shape = wanted.shape == actual.shape
        finite = (a["finite"] and b["finite"] and a["max_abs"] is not None and
                  b["max_abs"] is not None and a["power"] is not None and b["power"] is not None)
        item = {"expected": a, "loaded": b, "shape_match": same_shape,
                "finite": finite, "max_abs_error": None, "atol": None,
                "rtol": 1e-12, "valid": False}
        if not same_shape:
            report["errors"].append(f"BUFFER_SHAPE_MISMATCH:{name}")
        if not finite:
            report["errors"].append(f"BUFFER_NONFINITE_OR_NONNUMERIC:{name}")
        if same_shape and finite:
            atol = 1e-12 * max(a["max_abs"], b["max_abs"])
            with np.errstate(over="ignore", invalid="ignore"):
                error = np.abs(actual.astype(np.complex128) - wanted.astype(np.complex128))
                valid = bool(np.all(error <= atol + 1e-12 * np.abs(wanted.astype(np.complex128))))
            max_error = float(np.max(error)) if error.size else 0.0
            item.update(atol=atol, max_abs_error=max_error if math.isfinite(max_error) else None,
                        valid=valid)
            if not valid:
                report["errors"].append(f"BUFFER_VALUE_MISMATCH:{name}")
        report["components"][name] = item
    report["shape_match"] = all(item["shape_match"] for item in report["components"].values())
    report["finite"] = all(item["finite"] for item in report["components"].values())
    for field in ("expected", "loaded"):
        stats = [item[field] for item in report["components"].values()]
        report[field + "_nonzero"] = sum(item["nonzero"] or 0 for item in stats)
        if all(item["power"] is not None for item in stats):
            power = sum(item["power"] for item in stats)
            report[field + "_power"] = power if math.isfinite(power) else None
    report["valid"] = bool(all(item["valid"] for item in report["components"].values()) and
                           report["expected_power"] is not None and report["loaded_power"] is not None)
    return report


def assert_loaded_reference(sim, flux, expected, input_power, mp, comm):
    """Collectively verify load_minus_flux_data before starting the target run.

    ``expected`` is the ORIGINAL cached incident E/H data. Readback must equal
    its negative. Every rank calls get_flux_data, get_fluxes, and allgather in
    the same order, including when local input validation fails. A failed audit
    raises ValueError on every rank. This function never advances the simulation.
    """
    local = {"rank": int(comm.rank), "errors": [], "input_power": None,
             "normalized_flux": None, "buffer_check": None}
    try:
        if (isinstance(input_power, (bool, np.bool_)) or not np.isscalar(input_power) or
                np.iscomplexobj(input_power)):
            raise ValueError("not a real scalar")
        power = float(input_power)
        if not math.isfinite(power) or power <= 0:
            raise ValueError("not finite and positive")
        local["input_power"] = power
    except (TypeError, ValueError, OverflowError):
        local["errors"].append("REFERENCE_INPUT_POWER_INVALID")
    try:
        loaded = sim.get_flux_data(flux)
        original = _buffers(expected)
        if any(value.dtype.kind not in "iufc" for value in original.values()):
            raise ValueError("cached buffers must be numeric E/H arrays")
        negative = {name: -value.astype(np.complex128) for name, value in original.items()}
        local["buffer_check"] = buffer_check(negative, loaded)
        if not local["buffer_check"]["valid"]:
            local["errors"].append("REFERENCE_LOADED_BUFFER_MISMATCH")
    except Exception as exc:
        # A rank-local Python error must not strand peers in the flux collective.
        local["errors"].append(f"REFERENCE_BUFFER_READBACK_FAILED:{exc}")
    try:
        values = np.asarray(mp.get_fluxes(flux))
        if values.shape != (1,) or values.dtype.kind not in "iuf" or not np.isfinite(values).all():
            raise ValueError("requires one finite real flux")
        if local["input_power"] is not None:
            normalized = float(values[0]) / local["input_power"]
            if not math.isfinite(normalized):
                raise ValueError("nonfinite normalized flux")
            local["normalized_flux"] = normalized
            if abs(normalized + 1.0) > 1e-9:
                local["errors"].append("REFERENCE_LOADED_FLUX_NOT_MINUS_ONE")
    except Exception as exc:
        local["errors"].append(f"REFERENCE_FLUX_READBACK_FAILED:{exc}")
    ranks = comm.allgather(local)
    errors = sorted({error for report in ranks for error in report["errors"]})
    if len(ranks) != comm.size or sorted(report["rank"] for report in ranks) != list(range(comm.size)):
        errors.append("REFERENCE_AUDIT_RANKS_INVALID")
    global_powers = {}
    for name in ("expected_power", "loaded_power"):
        values = [(report.get("buffer_check") or {}).get(name) for report in ranks]
        total = sum(values) if all(value is not None for value in values) else None
        if total is None or not math.isfinite(total) or total <= 0:
            errors.append("REFERENCE_GLOBAL_" + name.upper() + "_NOT_POSITIVE")
            total = None if total is None or not math.isfinite(total) else total
        global_powers[name] = total
    inputs = [report["input_power"] for report in ranks]
    if all(value is not None for value in inputs) and inputs:
        if max(inputs) - min(inputs) > 1e-12 * max(inputs):
            errors.append("REFERENCE_INPUT_POWER_DIFFERS_ACROSS_RANKS")
    audit = {"passed": not errors, "ranks": ranks, "errors": errors,
             "global_buffer_power": global_powers, "normalized_flux_target": -1.0,
             "normalized_flux_tolerance": 1e-9, "comparison": "loaded buffers == -cached incident buffers"}
    if errors:
        error = ValueError("REFERENCE_LOADED_AUDIT_FAILED:" + ";".join(errors))
        error.audit = audit
        raise error
    return audit
