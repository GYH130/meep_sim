#!/usr/bin/env python3
"""Curate scientific observations from local evidence snapshots, without solver calls.

Usage: python extract_scientific_data.py --workspace PATH [--output PATH]
The private evidence snapshots are not distributed. Public JSON records preserve
their relative identifiers and SHA256 hashes, but omit operational identities.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

MATERIAL = "283197fec119be9ca9c2827179ade13c566e240df2ef787959522fcf159f7c09"
CURRENT_IMPL = "bb222ce35efe30bb664c1e59a10a3e6861f838c521c99f30b95c32a21832fead"
SOURCES = {
    "fill_final": "meep_screen/lab_fill_bridge_20261013/PROGRESS_CHECK_20261014_1345.json",
    "fill_preflight": "meep_screen/lab_fill_bridge_20261013/PREFLIGHT_VERIFIED.json",
    "baseline": "meep_screen/lab_fill_bridge_20261013/BASELINE_F065_VERIFIED.json",
    "predecessor": "meep_screen/lab_fill_bridge_20261013/PREDECESSOR_VERIFIED.json",
    "endpoints": "meep_screen/lab_parallel_exploration_restart_20261012/PROGRESS_CHECK_20261013_1537.json",
    "night": "meep_screen/lab_validation_night_20261011/PROGRESS_CHECK_20261012_0845.json",
    "source_controls": "meep_screen/lab_parallel_exploration_20261012/PROGRESS_CHECK_20261012_1547.json",
    "field_audit": "meep_screen/lab_monitor_polarization_diagnostic_20261012_v2/FIELD_DIAGNOSIS_20261012.json",
    "material_model": "meep_screen/lab_fill_bridge_20261013/assets/material/route_a_material.json",
    "material_fit": "meep_screen/lab_fill_bridge_20261013/assets/material/fit_report.json",
    "stability": "meep_screen/stability_diagnostic_20261011/diagnosis.json",
}

def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace", type=Path, required=True)
    p.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    args = p.parse_args()
    src, hashes = {}, {}
    for name, rel in SOURCES.items():
        raw = (args.workspace / rel).read_bytes()
        src[name] = json.loads(raw)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    manifest = [{"id": k, "relative_path": v, "sha256": hashes[k], "distributed": k == "material_model"} for k, v in SOURCES.items()]

    def provenance(*ids):
        return {"source_ids": list(ids), "source_sha256": {i: hashes[i] for i in ids}}

    def record(rid, stage, fill, wavelength, angle, pol, resolution, metrics, status, **extra):
        return dict(id=rid, stage=stage, fill=fill, wavelength_um=wavelength,
                    theta_deg=angle, polarization=pol, resolution=resolution,
                    metrics=metrics, checks=None, status=status, stop=None,
                    material_sha256=MATERIAL, result_sha256=None, case_sha256=None,
                    implementation_sha256=None, **extra)

    records = []
    night_by_slug = {r["slug"]: r for r in src["night"]["rows"]}
    for b in src["baseline"]["cases"]:
        c = b["case"]
        n = night_by_slug[b["slug"]]
        assert c["material_sha256"] == MATERIAL
        assert b["metrics"] == n["metrics"] and b["result_sha"] == n["result_sha256"]
        r = record(b["slug"], "lab_validation_night_20261011", c["geometry"]["surface_fill"], c["wavelength_um"],
                   c["emission_theta_deg"], c["polarization"], c["resolution"], b["metrics"], n["status"], **provenance("baseline", "night"))
        r.update(result_sha256=b["result_sha"], case_sha256=b["case_sha"], implementation_sha256=c["implementation_sha256"],
                 failed_checks=n["failed_checks"], stop={"converged": True, "decay_threshold": c["exploratory_structure_intensity_decay"],
                 "evidence_scope": "aggregate completion receipt confirms two converged structures; detailed stop samples not embedded"})
        records.append(r)

    predecessors = {r["slug"]: r for r in src["predecessor"]["tasks"]}
    endpoint_evidence = {r["slug"]: r for r in src["endpoints"]["final_evidence_verification"]["tasks"]}
    for e in src["endpoints"]["results"]:
        spec = predecessors[e["slug"]]
        ev = endpoint_evidence[e["slug"]]
        assert all(ev["evidence_hashes_match"].values())
        r = record(e["slug"], "lab_parallel_exploration_restart_20261012", spec["surface_fill"], spec["wavelength_um"],
                   spec["angle"], spec["polarization"], spec["resolution"], e["metrics"], e["status"], **provenance("endpoints", "predecessor"))
        r.update(checks=e["checks"], stop=ev["final_stop"], result_sha256=e["result_sha256"],
                 case_sha256=spec["case_sha256"], implementation_sha256=spec["implementation_sha256"],
                 reference_reused=e["reference_reused"], raw_evidence_hashes_verified=all(ev["evidence_hashes_match"].values()),
                 last_windows=ev["last_windows"], finished_utc=e["completed_utc"])
        records.append(r)

    fills = {r["slug"]: r for r in src["fill_preflight"]["verified"]["tasks"]}
    for f in src["fill_final"]["live_snapshot"]["tasks"]:
        spec = fills[f["slug"]]
        assert f["result_hash_matches"] and f["case_matches"] and f["stop_matches"]
        r = record(f["slug"], "lab_fill_bridge_20261013", spec["surface_fill"], spec["wavelength_um"], spec["angle"],
                   spec["polarization"], 20, f["result"]["metrics"], f["result"]["status"], **provenance("fill_final", "fill_preflight"))
        r.update(checks=f["result"]["checks"], stop=f["stop"], case_sha256=spec["case_sha256"], implementation_sha256=CURRENT_IMPL,
                 finished_utc=f["finished_utc"], result_hash_verified=True,
                 result_hash_note="Completion receipt verifies hash match but does not embed its value; value intentionally null.")
        records.append(r)
    records.sort(key=lambda r: (r["fill"], r["wavelength_um"], -r["theta_deg"]))
    assert len(records) == 10 and len({r["id"] for r in records}) == 10
    for r in records:
        m = r["metrics"]
        assert abs(m["R"] + m["T"] + m["A_flux"] - 1) < 1e-12
        r["derived"] = {"flux_minus_volume_pp": 100 * (m["A_flux"] - m["A_vol"]), "signed_T_percent": 100 * m["T"]}
    scope = "Exploratory 2D blind-slot p polarization r20; discrete wavelengths/angles, not formal qualification or a continuous band spectrum."
    dump(args.output / "current_structured_results.json", dict(schema_version=1, scope=scope, material_sha256=MATERIAL,
         value_units="All R/T/A metrics are dimensionless fractions; pp means percentage points.", records=records, source_manifest=manifest))

    flats = []
    for n in src["night"]["rows"]:
        if n["case_kind"] != "flat":
            continue
        r = record(n["slug"], "lab_validation_night_20261011", None, n["wavelength_um"], n["angle"], n["polarization"],
                   n["resolution"], n["metrics"], n["status"], **provenance("night"))
        r.update(result_sha256=n["result_sha256"], failed_checks=n["failed_checks"], analytic=n["analytical"], R_error_pp=n["R_error_percentage_points"])
        flats.append(r)
    assert len(flats) == 12
    dump(args.output / "flat_validation_results.json", dict(schema_version=1, scope="Flat controls; not structured full-width validation.",
         material_sha256=MATERIAL, counts={"passed": 4, "accuracy_failed": 8}, fresnel_R_absolute_tolerance=0.0025,
         records=flats, source_manifest=[m for m in manifest if m["id"] == "night"]))

    controls = []
    for d in src["source_controls"]["completed_diagnostics"]:
        component = d["slug"].split("_")[3]
        r = record(d["slug"], "lab_parallel_exploration_20261012", None, 10.5, 0,
                   d["slug"].split("_")[1][0], 20, d["metrics"], d["state"], **provenance("source_controls"))
        r.update(checks=d["checks"], result_sha256=d["result_sha256"], source_component=component,
                 source_type="electric" if component.startswith("E") else "magnetic", stop={"reason": d["reason"]})
        controls.append(r)
    audit = []
    for a in src["field_audit"]["saved_field_audits"]:
        audit.append({"id": a["slug"], "power_balance": a["analysis"]["power_balance"], "front_cross": a["analysis"]["front_cross"],
                      "result_sha256": a["input_hashes"]["result.json"], "raw_input_sha256": a["input_hashes"], "input_unchanged": a["input_unchanged"]})
    dump(args.output / "source_type_diagnostic.json", dict(schema_version=1, scope="Normal-incidence 10.5 um r20 flat controls only; not a physical structured s/p contrast.",
         records=controls, comparison=src["source_controls"]["diagnostic_interpretation"], saved_field_audits=audit,
         limitations=src["field_audit"]["limitations"], source_manifest=[m for m in manifest if m["id"] in ("source_controls", "field_audit")]))

    geometries = {g["surface_fill"]: g for g in src["fill_preflight"]["prepare"]["geometry_checks"]}
    base = src["baseline"]["cases"][0]["case"]
    # Same finite tilted-slot half-plane formula documented in ti2d/geometry.py.
    g = base["geometry"]
    c, s = math.cos(math.radians(g["tilt_deg"])), math.sin(math.radians(g["tilt_deg"]))
    width = g["period_um"] * g["surface_fill"] * c
    half_open = g["period_um"] * g["surface_fill"] / 2
    depth = g["axis_length_um"] * c + width * abs(s) / 2
    geometries[.65] = dict(g, surface_opening_um=2*half_open, width_normal_um=width,
        maximum_depth_um=depth, min_wall_um=g["period_um"]*c-width, residual_um=g["thickness_um"]-depth,
        vertices_x_depth=[[-half_open,0],[half_open,0],[g["axis_length_um"]*s+width*c/2,g["axis_length_um"]*c-width*s/2],
                          [g["axis_length_um"]*s-width*c/2,depth]],
        derivation="Calculated from preserved geometry parameters using the finite tilted-slot half-plane model; not a new solve.")
    settings_keys = ("resolution", "courant", "polarization", "mpi_ranks", "source_fwidth_fraction", "source_cutoff", "pml_um", "air_above_um", "air_below_um", "sample_interval_um", "stop_window_um", "stop_consecutive_windows", "intensity_decay", "flux_drift_tolerance", "strict_stop", "stop_normalization_policy", "eps_averaging", "boundary_kind", "boundary_profile", "exploratory_structure_intensity_decay", "minimum_post_source_time_um", "oxide_layer_in_model")
    settings = {k: base[k] for k in settings_keys}
    dump(args.output / "geometry_and_settings.json", dict(schema_version=1, material_sha256=MATERIAL, geometry=[geometries[f] for f in sorted(geometries)],
         settings=settings, fill_definition="F = air-slot surface opening width / period; not metal volume or surface fraction.",
         geometry_scope="2D periodic finite tilted blind slot; not a 3D drilled hole. Changing F also changes maximum cavity depth and residual thickness.",
         scope_note="Common settings apply to the 10 structured points; flat diagnostic recipes have their recorded wavelengths, polarizations and resolution.",
         qualification="Exploratory only. A_flux versus A_vol discrepancy is a diagnostic, not a complete uncertainty interval.",
         settings_origin="baseline case plus fill bridge preflight; current implementation wrappers differ from the night baseline and are tracked per record.",
         source_manifest=[m for m in manifest if m["id"] in ("baseline", "fill_preflight")]))

    pairs = []
    for fill, wl in sorted({(r["fill"],r["wavelength_um"]) for r in records}):
        group = {r["theta_deg"]: r for r in records if r["fill"] == fill and r["wavelength_um"] == wl}
        a, b = group[30]["metrics"], group[-30]["metrics"]
        pairs.append(dict(fill=fill,wavelength_um=wl,polarization="p",resolution=20,
            plus_id=group[30]["id"],minus_id=group[-30]["id"],A_flux_plus=a["A_flux"],A_flux_minus=b["A_flux"],
            A_vol_plus=a["A_vol"],A_vol_minus=b["A_vol"],difference_flux_pp=100*(a["A_flux"]-b["A_flux"]),
            difference_vol_pp=100*(a["A_vol"]-b["A_vol"]),ratio_flux=a["A_flux"]/b["A_flux"],ratio_vol=a["A_vol"]/b["A_vol"]))
    dump(args.output / "direction_pairs.json", dict(schema_version=1,definition="Difference and ratio between the two sampled signed angles +30 and -30 degrees, not an angular directivity factor or bandwidth.",records=pairs))
    assert hashes["material_model"] == MATERIAL
    # Byte-for-byte model asset, no source-table redistribution or new fitting.
    shutil.copyfile(args.workspace / SOURCES["material_model"], args.output / "current_model.json")
    fit_keys = ("accepted", "canonical_accepted", "local_accepted", "band_um", "required_band_um", "material_sha256", "source_sha256",
                "local_metrics", "local_thresholds", "normalized_max_error", "normalized_rmse", "original_node_max_error", "original_node_rmse",
                "interpolation", "normalization", "normalization_scale", "qualified_fdtd", "physical_sample_validated", "timestep_audits",
                "selection_reason", "interpretation", "lorentz_terms", "selected_candidate_index", "passivity_min_imag_epsilon")
    fit = {k: src["material_fit"][k] for k in fit_keys}
    fit.update(schema_version=1, source_manifest=[m for m in manifest if m["id"] in ("material_fit", "material_model")],
               distinction="Frequency-domain fit error against the transcribed optical table is separate from same-model FDTD-versus-Fresnel error and from fabricated-sample uncertainty.")
    dump(args.output / "material_fit_summary.json", fit)
    stability_keys = ("name", "scope", "polarization", "control", "material_sha256", "candidate_index", "original_epsilon_inf", "effective_epsilon_inf",
                      "courant", "resolution", "no_pml", "no_interfaces", "geometry", "wavelength_um", "theta_deg", "target_simulation_time", "status",
                      "end_simulation_time", "maximum_recorded_abs_e", "source_note")
    s = src["stability"]
    stability = dict(schema_version=1, diagnostic_date_utc=s["created_utc"], evidence_level=s["evidence_level"],
                     test_count=len(s["controlled_results"]), judgement=s["judgement"], caveats=s["caveats"],
                     theoretical_risk_estimates=s["theoretical_risk_estimates"],
                     records=[{k:r[k] for k in stability_keys if k in r} for r in s["controlled_results"]],
                     conclusions_scope="Sixteen short independent tests isolate the early divergence mechanism of the rejected wideband fit; they are not an accuracy or long-time qualification of the later stable production material.",
                     subsequent_context="The current production material is separately identified by current_model.json. Its later ten structured results are exploratory and still retain negative signed transmission and absorption-bookkeeping discrepancies.",
                     source_manifest=[m for m in manifest if m["id"] == "stability"])
    assert stability["test_count"] == 16
    dump(args.output / "stability_diagnosis_summary.json", stability)
    # Automated fail-closed leak check on curated JSON, not the private snapshots.
    forbidden = ("/home/", "/Users/", "192.168.", "100.125.", "Administrator", "DESKTOP-", "password", "private_key")
    for target in args.output.glob("*.json"):
        body = target.read_text(encoding="utf-8")
        assert not any(t in body for t in forbidden), f"Unexpected operational data in {target.name}"
    print(json.dumps({"structured":len(records),"flat_controls":len(flats),"source_controls":len(controls),"pairs":len(pairs),"output":str(args.output)}))

if __name__ == "__main__":
    main()
