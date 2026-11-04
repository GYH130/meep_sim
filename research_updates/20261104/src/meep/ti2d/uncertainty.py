"""Observed numerical error evidence, without claiming rigorous error bounds.

This module is a pure post-processing function: it does not read files, run a
solver, infer unmeasured errors, or transfer evidence between physical cases.
Only qualified results in the supplied manifest are considered. Mesh evidence
varies resolution alone; strict-stop evidence varies only the two declared
stopping controls. Every other saved configuration value must match exactly.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from numbers import Real


_METRICS = ("R", "T", "A_flux")
_IDENTITY_FIELDS = ("geometry", "case", "material_sha256", "wavelength_um",
                    "polarization", "emission_theta_deg", "resolution", "courant",
                    "implementation_sha256", "strict_stop", "intensity_decay")
_LABEL_FIELDS = ("geometry_id", "resolution", "polarization", "emission_theta_deg", "strict_stop")
_NON_PHYSICAL_FIELDS = {"id", "task_id", "phase", "output", "output_dir"}


def _finite(value):
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(float(value))


def _configuration(result, task, manifest):
    config = result.get("case_config")
    if not isinstance(config, Mapping) or any(name not in config for name in _IDENTITY_FIELDS):
        return None, "complete_saved_case_config_required"
    if (not isinstance(config["geometry"], Mapping) or not config["geometry"]
            or config["case"] not in ("flat", "structure")
            or config["polarization"] not in ("s", "p")
            or not isinstance(config["strict_stop"], bool)):
        return None, "invalid_saved_case_identity"
    if any(not _finite(config[name]) for name in
           ("wavelength_um", "emission_theta_deg", "resolution", "courant", "intensity_decay")):
        return None, "nonfinite_saved_case_identity"
    if min(config[name] for name in ("wavelength_um", "resolution", "courant", "intensity_decay")) <= 0:
        return None, "nonpositive_saved_case_control"
    if any(not isinstance(config[name], str) or not config[name]
           for name in ("material_sha256", "implementation_sha256")):
        return None, "missing_material_or_implementation_identity"
    for name in _LABEL_FIELDS:
        if name in task and config.get(name) != task[name]:
            return None, "saved_case_disagrees_with_manifest:" + name
    if manifest.get("material_sha256") not in (None, config["material_sha256"]):
        return None, "saved_material_disagrees_with_manifest"
    if manifest.get("implementation_sha256") not in (None, config["implementation_sha256"]):
        return None, "saved_implementation_disagrees_with_manifest"
    return dict(config), None


def _same_config(left, right, varying=()):
    excluded = _NON_PHYSICAL_FIELDS | set(varying)
    return ({name: value for name, value in left.items() if name not in excluded}
            == {name: value for name, value in right.items() if name not in excluded})


def _spans(records):
    answer = {}
    for metric in _METRICS:
        values = [record["result"].get("metrics", {}).get(metric) for record in records]
        answer[metric] = float(max(values) - min(values)) if all(_finite(value) for value in values) else None
    return answer


def _matching_comparisons(manifest, records, target, kind):
    varying = ("resolution",) if kind == "mesh" else ("strict_stop", "intensity_decay")
    allowed_kinds = ("mesh", "mesh_boundary") if kind == "mesh" else ("strict",)
    candidates = []
    for comparison in manifest.get("gate_comparisons", []):
        if comparison.get("kind") not in allowed_kinds:
            continue
        ids = comparison.get("ids", [])
        if len(ids) < 2 or len(set(ids)) != len(ids) or any(name not in records for name in ids):
            continue
        members = [records[name] for name in ids]
        if any(member["config"] is None for member in members):
            continue
        # The target itself must have an exact sampled counterpart; being an
        # unsampled intermediate resolution is not matching mesh evidence.
        if not any(_same_config(target, member["config"]) for member in members):
            continue
        if not all(_same_config(target, member["config"], varying) for member in members):
            continue
        if kind == "mesh":
            if len({member["config"]["resolution"] for member in members}) < 2:
                continue
        else:
            if len(members) != 2 or {member["config"]["strict_stop"] for member in members} != {False, True}:
                continue
            strict = next(member["config"] for member in members if member["config"]["strict_stop"])
            ordinary = next(member["config"] for member in members if not member["config"]["strict_stop"])
            if strict["intensity_decay"] >= ordinary["intensity_decay"]:
                continue
        candidates.append({"ids": list(ids), "kind": comparison["kind"], "values": _spans(members)})
    values = {metric: max((item["values"][metric] for item in candidates
                           if item["values"][metric] is not None), default=None)
              for metric in _METRICS}
    return values, [{"ids": item["ids"], "kind": item["kind"]} for item in candidates]


def _tail_ranges(result):
    """Full cross-window range, never the maximum of three within-window ranges."""
    stop = result.get("stop", {})
    windows, trace = stop.get("windows"), stop.get("trace")
    if not isinstance(windows, list) or len(windows) < 3 or not isinstance(trace, list):
        return None
    windows = windows[-3:]
    if any(not isinstance(window, Mapping) or not _finite(window.get("start"))
           or not _finite(window.get("end")) or window["end"] <= window["start"] for window in windows):
        return None
    if any(not math.isclose(first["end"], second["start"], abs_tol=1e-9, rel_tol=1e-12)
           for first, second in zip(windows, windows[1:])):
        return None
    if any(not isinstance(item, Mapping) or not _finite(item.get("time")) for item in trace):
        return None
    start, end = windows[0]["start"], windows[-1]["end"]
    samples = [item for item in trace if start < item["time"] <= end + 1e-10]
    if any(sum(window["start"] < item["time"] <= window["end"] + 1e-10 for item in samples) < 2
           for window in windows):
        return None
    if any(not _finite(item.get(metric)) for item in samples for metric in ("R", "T")):
        return None
    return {"R": float(max(item["R"] for item in samples) - min(item["R"] for item in samples)),
            "T": float(max(item["T"] for item in samples) - min(item["T"] for item in samples)),
            "start": float(start), "end": float(end), "samples": len(samples)}


def derive_error_evidence(manifest, results):
    """Return flat, JSON-safe diagnostic rows for qualified manifest results.

    ``A_error_estimate`` sums the observed mesh A span, strict-stop A
    difference, and |A_vol-A_flux| only when all three have exact-case
    coverage. Available finite-film Fresnel A differences and the last three
    complete windows' R-range + T-range are added. Missing optional evidence
    remains explicit; it is never silently represented as a measured zero.

    This sum is an empirical diagnostic, not a confidence interval, proven
    upper bound, extrapolated continuum error, or validated formal bound.
    A high-resolution member cannot inherit strict-stop evidence from a
    different resolution, angle, geometry, polarization, or material.
    """
    if not isinstance(manifest, Mapping) or not isinstance(results, Mapping):
        raise ValueError("manifest and results must be mappings")
    records = {}
    for task in manifest.get("tasks", []):
        if not isinstance(task, Mapping) or not isinstance(task.get("id"), str):
            raise ValueError("manifest tasks require string ids")
        name = task["id"]
        if name in records:
            raise ValueError("duplicate qualified task id")
        result = results.get(name)
        if not isinstance(result, Mapping) or result.get("status") != "QUALIFIED":
            continue
        config, identity_error = _configuration(result, task, manifest)
        records[name] = {"task": task, "result": result, "config": config, "identity_error": identity_error}
    rows = []
    for name, record in records.items():
        task, result, config = record["task"], record["result"], record["config"]
        identity = config or task
        row = {"id": name, "phase": task.get("phase"), "status": "QUALIFIED",
               "geometry_id": identity.get("geometry_id"), "case": identity.get("case"),
               "resolution": identity.get("resolution"), "polarization": identity.get("polarization"),
               "emission_theta_deg": identity.get("emission_theta_deg"),
               "wavelength_um": identity.get("wavelength_um"),
               "material_sha256": identity.get("material_sha256"),
               "implementation_sha256": identity.get("implementation_sha256"),
               "identity_verified": config is not None, "identity_error": record["identity_error"],
               "estimate_kind": "empirical_not_rigorous", "validated": False,
               "formal_error_bound_available": False, "A_error_estimate": None,
               "estimate_components": {}, "missing_coverage": []}
        empty = {metric: None for metric in _METRICS}
        mesh, mesh_sources = _matching_comparisons(manifest, records, config, "mesh") if config else (empty, [])
        strict, strict_sources = _matching_comparisons(manifest, records, config, "strict") if config else (empty, [])
        row.update(mesh_comparisons=mesh_sources, strict_stop_comparisons=strict_sources)
        for metric in _METRICS:
            row["mesh_span_" + metric] = mesh[metric]
            row["strict_stop_delta_" + metric] = strict[metric]
            row["flat_fresnel_delta_" + metric] = None
        if any(mesh[metric] is None for metric in _METRICS):
            row["missing_coverage"].append("mesh")
        if any(strict[metric] is None for metric in _METRICS):
            row["missing_coverage"].append("strict_stop")
        metrics = result.get("metrics", {})
        a_flux, a_volume = metrics.get("A_flux"), metrics.get("A_vol")
        absorption = abs(float(a_volume) - float(a_flux)) if _finite(a_flux) and _finite(a_volume) else None
        row["absorption_delta_A"] = absorption
        if absorption is None:
            row["missing_coverage"].append("absorption")
        analytical = result.get("analytical")
        if config is not None and config["case"] == "flat":
            for metric, analytical_name in (("R", "R"), ("T", "T"), ("A_flux", "A")):
                value = analytical.get(analytical_name) if isinstance(analytical, Mapping) else None
                if _finite(value) and _finite(metrics.get(metric)):
                    row["flat_fresnel_delta_" + metric] = abs(float(metrics[metric]) - float(value))
            if any(row["flat_fresnel_delta_" + metric] is None for metric in _METRICS):
                row["missing_coverage"].append("flat_fresnel")
        tail = _tail_ranges(result)
        row.update(tail_R_range=tail["R"] if tail else None,
                   tail_T_range=tail["T"] if tail else None,
                   tail_A_estimate=tail["R"] + tail["T"] if tail else None,
                   tail_start=tail["start"] if tail else None,
                   tail_end=tail["end"] if tail else None,
                   tail_samples=tail["samples"] if tail else 0)
        if tail is None:
            row["missing_coverage"].append("tail")
        if config is None:
            row["missing_coverage"].append("exact_case_identity")
        if all(value is not None for value in (mesh["A_flux"], strict["A_flux"], absorption)):
            components = {"mesh_span_A_flux": mesh["A_flux"], "strict_stop_delta_A_flux": strict["A_flux"],
                          "absorption_delta_A": absorption}
            if row["flat_fresnel_delta_A_flux"] is not None:
                components["flat_fresnel_delta_A_flux"] = row["flat_fresnel_delta_A_flux"]
            if row["tail_A_estimate"] is not None:
                components["tail_A_estimate"] = row["tail_A_estimate"]
            row["estimate_components"] = components
            row["A_error_estimate"] = float(sum(components.values()))
        row["interpretation"] = ("Observed finite numerical differences only; no rigorous bound or confidence level. "
                                 "Missing coverage is unmeasured, and no formal error bound is validated.")
        rows.append(row)
    return rows
