#!/usr/bin/env python3
"""Assess one hash-bound six-model saved-field export without relabeling history."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import re
from pathlib import Path

import floquet_reconstruct
from floquet_reconstruct import FRACTIONS, analyze_case, closure_diagnostic, finite, parse_export

IDS = (
    "F_RM1_A3H80_t3_l10p5_th15_phi30_p_cap22",
    "F_RM2_A3H80_t3_l10p5_th60_phi15_p_cap22",
    "AP2_A3H80_t3_l10p5_th15_phi30_p",
    "AP4_A3H80_t3_l10p5_th60_phi15_p",
    "RM1_A3H80_t3_l10p5_th15_phi30_p_cap22",
    "RM2_A3H80_t3_l10p5_th60_phi15_p_cap22",
)
FLUX_WARNINGS = {"R_orders_vs_net_flux", "flux_vs_total_loss_absolute", "flux_vs_total_loss_relative"}
TRUNCATION_SENSITIVITY_GATE = 0.0002
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def markers(text, name):
    prefix = "CODEX_FLOQUET_" + name + "_JSON"
    out = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(prefix) and line[len(prefix):].lstrip().startswith(("=", ":", "{")):
            out.append(json.loads(line[len(prefix):].lstrip(" :=\t")))
    return out


def unique_by_id(rows, label):
    require(isinstance(rows, list) and len(rows) == 6, f"{label}: exactly six records required")
    out = {}
    for r in rows:
        require(isinstance(r, dict), f"{label}: records must be objects")
        cid = r.get("case_id")
        require(cid in IDS and cid not in out, f"{label}: unknown or duplicate case {cid}")
        out[cid] = r
    require(set(out) == set(IDS), f"{label}: closed six IDs differ")
    return out


def verify_native(bundle):
    require(isinstance(bundle, dict), "completed export object required")
    text = bundle.get("native_stdout")
    require(isinstance(text, str), "native stdout text required")
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    require(digest == bundle.get("native_stdout_sha256"), "native stdout UTF-8 SHA256 mismatch")
    status, integrity = bundle.get("status", {}), bundle.get("integrity", {})
    require(status.get("state") == "COMPLETE_SAVED_FLOQUET_EXPORT_AWAITING_INDEPENDENT_ANALYSIS", "native controller not complete")
    require(type(status.get("native_exit_code")) is int and status["native_exit_code"] == 0, "native exit code not zero")
    require(status.get("resource_guard_tripped") is False, "resource guard status not explicitly clear")
    require(status.get("solver_executed") is False and status.get("solver_runs") == 0 and status.get("model_saved") is False, "controller is not explicitly read-only")
    require(integrity.get("source_and_copy_sha256_verified") is True and integrity.get("cases") == 6 and integrity.get("planes") == 30, "source/copy integrity not verified for six models")
    require(integrity.get("solver_executed") is False and integrity.get("model_saved") is False, "integrity receipt not read-only")
    begin, done = markers(text, "BEGIN"), markers(text, "DONE")
    require(len(begin) == len(done) == 1, "exactly one native BEGIN and DONE required")
    require(not markers(text, "ERROR"), "native ERROR marker present")
    begin, done = begin[0], done[0]
    require(begin.get("expected_cases") == 6 and begin.get("expected_planes") == 30 and begin.get("solver_runs") == 0 and begin.get("model_saved") is False, "BEGIN scope differs")
    require(done.get("state") == "COMPLETE_EXPORT_NOT_PHYSICAL_ACCEPTANCE" and done.get("completed") is True and done.get("planned_cases") == 6 and done.get("successful_cases") == 6 and done.get("successful_planes") == 30, "DONE is not complete six-case thirty-plane export")
    require(done.get("solver_runs") == 0 and done.get("solver_executed") is False and done.get("model_saved") is False and done.get("physical_acceptance") is False, "DONE unsafe execution or relabeling")
    require(status.get("native_done") == done, "controller native_done differs from stdout DONE")
    require(begin.get("runtime_plan_sha256") == integrity.get("runtime_plan_sha256") and HEX64.fullmatch(str(begin.get("runtime_plan_sha256"))), "runtime plan SHA mismatch")
    cb = unique_by_id(markers(text, "CASE_BEGIN"), "CASE_BEGIN")
    ce = unique_by_id(markers(text, "CASE_END"), "CASE_END")
    external_integrity = any(e.get("model_integrity_check_location") == "external_controller" for e in ce.values())
    external_hashes = None
    if external_integrity:
        require(all(e.get("model_integrity_check_location") == "external_controller" for e in ce.values()), "mixed native/external integrity protocols forbidden")
        external_hashes = unique_by_id(integrity.get("input_hashes"), "external controller input hashes")
        for cid, proof in external_hashes.items():
            require(proof.get("original_and_copy_actual_SHA_verified") is True, f"{cid}: actual original/copy hashes not verified externally")
            require(HEX64.fullmatch(str(proof.get("model_sha256"))) is not None and HEX64.fullmatch(str(proof.get("metrics_sha256"))) is not None, f"{cid}: missing external model/metrics SHA256")
    else:
        require(all("model_integrity_check_location" not in e for e in ce.values()), "unknown integrity protocol cannot default to native verification")
    cases = unique_by_id(parse_export(text), "parsed cases")
    planes = markers(text, "PLANE")
    require(len(planes) == 30, "exactly thirty PLANE markers required")
    for cid in IDS:
        b, e, c = cb[cid], ce[cid], cases[cid]
        require(b.get("solver_executed") is False and b.get("model_saved") is False and b.get("physical_acceptance") is False, f"{cid}: CASE_BEGIN unsafe flags")
        sha = b.get("model_sha256_expected") if external_integrity else b.get("model_sha256")
        require(HEX64.fullmatch(str(sha)) is not None, f"{cid}: model SHA256 expected binding missing")
        if external_integrity:
            require(b.get("model_integrity_check_location") == "external_controller" and "model_sha256" not in b, f"{cid}: native CASE_BEGIN must label expected SHA and external verification explicitly")
            require(e.get("model_sha256_expected") == sha == external_hashes[cid]["model_sha256"], f"{cid}: native embedded SHA and external actual SHA differ")
            require(all(k not in e for k in ("model_sha256_before", "model_sha256_after", "model_integrity_unchanged")), f"{cid}: external native marker must not falsely claim filesystem verification")
        else:
            require("model_sha256_expected" not in b and "model_integrity_check_location" not in b, f"{cid}: mixed CASE_BEGIN integrity protocol forbidden")
            require(e.get("model_sha256_before") == e.get("model_sha256_after") == sha and e.get("model_integrity_unchanged") is True, f"{cid}: model SHA differs before/after or native integrity not proved")
        # This is a separately normalized internal copy, not a fabricated native
        # marker. The immutable stdout and cb/ce dictionaries remain untouched.
        cases[cid] = {**c, "model_sha256": sha, "model_sha256_verification_source": "external_controller" if external_integrity else "native_before_after_and_external_controller", "model_sha256_native_field": "model_sha256_expected" if external_integrity else "model_sha256"}
        c = cases[cid]
        require(e.get("completed") is True and e.get("solver_runs") == 0 and e.get("model_saved") is False and e.get("physical_acceptance") is False, f"{cid}: CASE_END incomplete or unsafe")
        require(len(c["planes"]) == 5 and {p.get("plane_index") for p in c["planes"]} == set(range(5)), f"{cid}: five unique plane indices required")
        c["planes"].sort(key=lambda p: p["plane_index"])
        for j, p in enumerate(c["planes"]):
            require(p.get("cap_fraction") == FRACTIONS[j] and p.get("area_passed") is True and p.get("no_order_omitted") is True and p.get("no_missing_value_replaced_with_zero") is True, f"{cid}: fixed full-plane proof differs")
            require(p.get("physical_acceptance") is False and p.get("quadrature_order") == 8, f"{cid}: plane evaluation contract differs")
            require(p.get("order_count") == b.get("order_count") == len(p.get("orders", [])), f"{cid}: order count mismatch")
            coupons = p.get("fixed_reciprocal_orthogonality_coupons")
            require(isinstance(coupons, list), f"{cid}: orthogonality coupon inventory missing")
            expected = 0 if cid.startswith("F_") else 6
            require(len(coupons) == expected, f"{cid}: fixed orthogonality coupon inventory differs")
            for q in coupons:
                require(q.get("passed") is True and q.get("tolerance") == 1e-6 and finite(q.get("absolute_complex_defect"), "orthogonality defect") <= 1e-6, f"{cid}: native orthogonality coupon failed")
    original = unique_by_id(bundle.get("original_metrics"), "original metrics")
    if "input_ids" in status:
        require(set(status["input_ids"]) == set(IDS) and len(status["input_ids"]) == 6, "controller input IDs differ")
    return cases, original, {"native_stdout_sha256_verified": digest, "runtime_plan_sha256": begin["runtime_plan_sha256"], "native_complete": True, "case_count": 6, "plane_count": 30, "solver_runs": 0, "models_saved": False, "model_integrity_verification_mode": "external_controller_actual_original_and_copy_SHA_verified" if external_integrity else "native_before_after_and_controller_copy_SHA_verified", "external_input_hashes_verified": list(external_hashes.values()) if external_integrity else None, "native_done": done}


def minmax(values):
    return [float(min(values)), float(max(values))]


def original_gate_summary(metrics):
    checks = metrics.get("checks")
    require(isinstance(checks, dict) and checks, "original checks dictionary required")
    require(all(type(v) is bool for v in checks.values()), "original checks must remain explicit booleans")
    nonflux = {k: v for k, v in checks.items() if k not in FLUX_WARNINGS}
    source_names = ("requested_incident_wavevector", "source_mode_purity", "independent_incident_power", "shared_source_fit", "requested_global_incidence", "incident_wavevector_points_down", "side_floquet_source_binding", "explicit_global_source_binding", "normalized_source_mode_global_direction", "internal_E_fit_condition", "internal_E_fit_relative_L2_residual", "internal_E_fit_all_LOO_Pin", "internal_E_fit_cut_areas")
    source = {k: checks[k] for k in source_names if k in checks}
    require(source, "original source checks absent")
    return {"original_state": metrics.get("state"), "original_checks_unchanged": checks, "original_source_checks_separate": source, "original_source_checks_passed": all(source.values()), "original_non_rawflux_hard_checks_passed": all(nonflux.values()), "original_raw_flux_flags_unchanged": {k: checks[k] for k in sorted(FLUX_WARNINGS) if k in checks}, "original_top_raw_H_power_diagnostic_unchanged": metrics.get("original_top_raw_H_power_diagnostic"), "historical_accepted_value_unchanged": metrics.get("accepted"), "historical_results_not_reclassified": True}


def summarize_model(cid, case, metrics, analysis):
    flat = cid.startswith("F_")
    loss_key = "A_Ti_H" if flat else "A_loss_total"
    loss = finite(metrics.get(loss_key), loss_key)
    old = original_gate_summary(metrics)
    closure = closure_diagnostic(analysis, loss)
    summarized = []
    analytic_a = 1-finite(metrics.get("analytic_IBC_R"), "analytic_IBC_R") if flat else None
    allpass = analysis["independent_inventory_and_area_proof"]["verified"] and old["original_source_checks_passed"] and old["original_non_rawflux_hard_checks_passed"]
    for shell in analysis["fixed_nested_cutoffs"]:
        cutoff = shell["fixed_transverse_cutoff_over_k0"]
        row = {"cutoff_k0": cutoff, "mode_count": shell["mode_count"], "complete_no_mode_silently_removed": shell["complete_no_mode_silently_removed"]}
        if not shell["complete_no_mode_silently_removed"]:
            row["method_candidate_checks_passed"] = False
            summarized.append(row)
            allpass = False
            continue
        rows = [r for r in closure["rows"] if r.get("cutoff") == cutoff]
        row.update({"field_prediction_checks_passed": shell["field_prediction_checks_passed"], "global_relative_L2_E_residual": shell["full_fit_relative_L2_E_residual"], "max_LOO_heldout_relative_L2_E_residual": max(shell["LOO_held_out_relative_L2_E_residual"]), "max_full_or_LOO_design_condition": max(max([m["design_condition_2"], *[q["design_condition_2"] for q in m["leave_one_out"]]]) for m in analysis["mode_results"] if cutoff == 1.2 or m["padding11"]), "zeroth_Pin_W": shell["zeroth_order_incident_W"], "all_LOO_zeroth_Pin_W": shell["all_LOO_zeroth_order_incident_W"], "max_zeroth_full_or_LOO_Pin_error_W": max(abs(v-1) for v in [shell["zeroth_order_incident_W"], *shell["all_LOO_zeroth_order_incident_W"]]), "source_Pin_and_LOO_check_passed": shell["zeroth_and_all_LOO_specified_1W_check_passed"] and shell["specified_1W_incident_power_check_passed"], "total_propagating_incoming_W": shell["total_propagating_incoming_W"], "nonzero_order_incoming_W_preserved": shell["nonzero_order_incoming_W_preserved"], "evanescent_up_down_cross_W_preserved": shell["evanescent_up_down_cross_W"], "closure_absolute_threshold_unchanged": .005, "closure_relative_threshold_unchanged": .02, "measurements": []})
        for kind, key in (("measured_E_fitted_H", "all_five_A_measured_E_fitted_H_over_specified_Pin"), ("LOO_measured_E_predicted_H", "all_five_A_LOO_over_specified_Pin")):
            comparison = next(r for r in rows if r["kind"] == kind)
            values = shell[key]
            measure = {"method": kind, "all_five_A": values, "A_minmax": minmax(values), "max_absolute_closure_error": max(comparison["all_five_absolute_mismatch"]), "max_relative_closure_error": max(comparison["all_five_relative_mismatch"]), "all_five_preserved_closure_gates_passed": comparison["all_five_preserved_closure_gates_passed"]}
            if flat:
                measure.update({"analytic_flat_A": analytic_a, "max_absolute_error_vs_1_minus_analytic_R": max(abs(a-analytic_a) for a in values), "analytic_absolute_threshold": .002, "all_five_analytic_flat_check_passed": all(abs(a-analytic_a) <= .002 for a in values)})
            row["measurements"].append(measure)
        row["method_candidate_checks_passed"] = row["field_prediction_checks_passed"] and row["source_Pin_and_LOO_check_passed"] and all(m["all_five_preserved_closure_gates_passed"] and (not flat or m["all_five_analytic_flat_check_passed"]) for m in row["measurements"])
        allpass = allpass and row["method_candidate_checks_passed"]
        summarized.append(row)
    cutoff = analysis["fixed_cutoff_truncation_comparison"]
    truncation_max = None
    if cutoff["complete"]:
        two = analysis["fixed_nested_cutoffs"]
        truncation_max = max(abs(b-a) for key in ("all_five_A_measured_E_fitted_H_over_specified_Pin", "all_five_A_LOO_over_specified_Pin") for a,b in zip(two[0][key], two[1][key]))
    stable = truncation_max is not None and truncation_max <= TRUNCATION_SENSITIVITY_GATE
    return {"case_id": cid, "role": "analytic_flat_control" if flat else "coarse_cavity" if cid.startswith("AP") else "refined_cavity", "model_sha256": case["model_sha256"], "model_sha256_verification_source": case["model_sha256_verification_source"], "model_sha256_native_field": case["model_sha256_native_field"], "comparison_loss_key": loss_key, "comparison_loss": loss, "original_A_port_unchanged": metrics.get("A_port"), "original_A_flux_unchanged": metrics.get("A_flux", metrics.get("A_net")), "original_gates_separate_and_unchanged": old, "two_fixed_cutoff_summaries": summarized, "truncation_sensitivity": {"new_sensitivity_gate_not_old_closure_gate": TRUNCATION_SENSITIVITY_GATE, "max_all_plane_and_LOO_cutoff12_minus11_A_abs": truncation_max, "passed": stable}, "saved_point_method_candidate_checks_passed": bool(allpass and stable), "formal_qualified": False, "historical_raw_H_flags_changed": False, "fitself_conservation_used_as_validation": False, "full_reconstruction_diagnostic": analysis}


def mesh_comparisons(models):
    out=[]
    for coarse, fine in ((IDS[2],IDS[4]),(IDS[3],IDS[5])):
        rows=[]
        for old,new in zip(models[coarse]["two_fixed_cutoff_summaries"], models[fine]["two_fixed_cutoff_summaries"]):
            for a,b in zip(old.get("measurements",[]),new.get("measurements",[])):
                delta=[v-u for u,v in zip(a["all_five_A"],b["all_five_A"])]
                rows.append({"cutoff_k0":old["cutoff_k0"],"method":a["method"],"all_five_fine_minus_coarse_A":delta,"all_five_fine_minus_coarse_percentage_points":[100*v for v in delta],"max_absolute_fine_minus_coarse_A":max(abs(v) for v in delta)})
        out.append({"coarse_case_id":coarse,"fine_case_id":fine,"aligned_all_five_plane_comparisons":rows,"single_level_sensitivity_not_global_convergence":True})
    return out


def assess(bundle):
    cases, original, verification = verify_native(bundle)
    models={cid:summarize_model(cid,cases[cid],original[cid],analyze_case(cases[cid])) for cid in IDS}
    supported=all(m["saved_point_method_candidate_checks_passed"] for m in models.values())
    source_hashes={"floquet_reconstruct.py":hashlib.sha256(Path(floquet_reconstruct.__file__).read_bytes()).hexdigest(),"assess_export.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    return {"schema":"SAVED_FLOQUET_METHOD_REPAIR_ASSESSMENT_V1","created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"analysis_source_sha256":source_hashes,"state":"METHOD_REPAIR_SUPPORTED_ON_SAVED_POINTS_ONLY" if supported else "METHOD_REPAIR_NOT_YET_SUPPORTED_RETAIN_DIAGNOSTICS","method_repair_supported_on_saved_points":supported,"formal_qualified":False,"band_directionality_validated":False,"historical_accepted_or_raw_H_flags_changed":False,"source_material_geometry_or_solver_modified":False,"solver_executed":False,"native_evidence_verification":verification,"retained_original_closure_thresholds":{"absolute":.005,"relative":.02},"new_truncation_sensitivity_gate_separate":TRUNCATION_SENSITIVITY_GATE,"validation_policy":"All six saved models; both predetermined cutoffs; every fixed plane plus every held-out prediction; no favorable selection; fitted self-conservation is not validation","models":[models[cid] for cid in IDS],"aligned_coarse_fine_comparisons":mesh_comparisons(models)}


def compact_table(result):
    lines=["case | cutoff | A measured min–max % | A LOO min–max % | max closure pp / relative % | method checks"]
    for m in result["models"]:
        for s in m["two_fixed_cutoff_summaries"]:
            if not s.get("measurements"):
                lines.append(f"{m['case_id']} | {s['cutoff_k0']} | INCOMPLETE")
                continue
            a,b=s["measurements"]
            maxabs=max(a["max_absolute_closure_error"],b["max_absolute_closure_error"])
            maxrel=max(a["max_relative_closure_error"],b["max_relative_closure_error"])
            lines.append(f"{m['case_id']} | {s['cutoff_k0']:.1f} | {100*a['A_minmax'][0]:.5f}–{100*a['A_minmax'][1]:.5f} | {100*b['A_minmax'][0]:.5f}–{100*b['A_minmax'][1]:.5f} | {100*maxabs:.5f} / {100*maxrel:.4f} | {s['method_candidate_checks_passed']}")
    lines.append(result["state"]+"; formal_qualified=false; historical flags unchanged")
    return "\n".join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=assess(json.loads(args.input.read_text(encoding="utf-8-sig")))
    result["completed_export_file_sha256"]=hashlib.sha256(args.input.read_bytes()).hexdigest()
    with args.output.open("x",encoding="utf-8") as f:
        json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
        f.write("\n")
    print(compact_table(result))


if __name__=="__main__":
    main()

