#!/usr/bin/env python3
"""Saved-field E-only Floquet/Maxwell reconstruction, never solver acceptance.

Convention exp(+i omega t), up exp(-i kz dz), Re(kz)>=0, Im(kz)<=0.
No S parameter or dissipation value participates in any fit or normalization.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

FRACTIONS = (0.1, 0.25, 0.5, 0.75, 0.9)
ABS_GATE = 0.005
REL_GATE = 0.02
CONDITION_LIMIT = 10.0
RESIDUAL_LIMIT = 0.005
SCHEMA = "E_ONLY_FLOQUET_MAXWELL_DIAGNOSTIC_V1"


def finite(value: Any, name: str, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name}: finite number required")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name}: finite {'positive ' if positive else ''}number required")
    return value


def cpair(value: Any, name: str) -> complex:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{name}: [real, imaginary] required")
    return complex(finite(value[0], name), finite(value[1], name))


def pair(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def parse_export(text: str) -> list[dict]:
    """Accept a single case, {cases:[...]}, or native JSON marker stream."""
    stripped = text.lstrip("\ufeff \t\r\n")
    if stripped.startswith("{"):
        data = json.loads(stripped)
        cases = data.get("cases", [data])
        if not isinstance(cases, list) or not cases:
            raise ValueError("nonempty cases array required")
        return cases
    cases: dict[str, dict] = {}
    for line in text.splitlines():
        for marker in ("CODEX_FLOQUET_CASE_BEGIN_JSON", "CODEX_FLOQUET_PLANE_JSON"):
            if marker not in line:
                continue
            rest = line.split(marker, 1)[1].lstrip(" :=\t")
            row = json.loads(rest)
            case_id = row.get("case_id")
            if not isinstance(case_id, str) or not case_id:
                raise ValueError("marker case_id required")
            if marker.endswith("CASE_BEGIN_JSON"):
                if case_id in cases:
                    raise ValueError(f"duplicate case begin: {case_id}")
                cases[case_id] = {**row, "planes": []}
            else:
                if case_id not in cases:
                    raise ValueError(f"plane before case begin: {case_id}")
                cases[case_id]["planes"].append(row)
    if not cases:
        raise ValueError("no Floquet cases found")
    return list(cases.values())


def validate_case(case: dict) -> tuple[np.ndarray, list[tuple[int, int]], dict]:
    if not isinstance(case, dict) or not isinstance(case.get("case_id"), str):
        raise ValueError("case_id and object required")
    lam = finite(case.get("lambda_m"), "lambda_m", True)
    finite(case.get("Z0"), "Z0", True)
    pin = finite(case.get("nominal_Pin_W"), "nominal_Pin_W", True)
    if not math.isclose(pin, 1.0, rel_tol=0, abs_tol=1e-14):
        raise ValueError("specified 1 W source normalization must remain unchanged")
    height = finite(case.get("height_m"), "height_m")
    cap = finite(case.get("cap_m"), "cap_m", True)
    if tuple(case.get("cap_fractions", ())) != FRACTIONS:
        raise ValueError("the five predetermined cap fractions are required")
    planes = case.get("planes")
    if not isinstance(planes, list) or len(planes) != 5:
        raise ValueError("exactly five saved-field planes required")
    z = np.array([finite(p.get("z_m"), "z_m") for p in planes])
    if not np.allclose(z, height + cap * np.array(FRACTIONS), atol=1e-14, rtol=1e-10):
        raise ValueError("plane heights disagree with fixed cap fractions")
    areas = np.array([finite(p.get("area_m2"), "area_m2", True) for p in planes])
    if not np.allclose(areas, areas[0], rtol=1e-6, atol=0):
        raise ValueError("full-cell areas differ among planes")
    k0 = 2 * math.pi / lam
    inventory: dict[tuple[int, int], dict] = {}
    keys: list[tuple[int, int]] = []
    for j, plane in enumerate(planes):
        orders = plane.get("orders")
        if not isinstance(orders, list) or not orders:
            raise ValueError("nonempty complete order inventory required")
        seen = set()
        for row in orders:
            if any(type(row.get(k)) is not int for k in ("n1", "n2")):
                raise ValueError("integer reciprocal order indices required")
            key = row["n1"], row["n2"]
            if key in seen:
                raise ValueError(f"duplicate order {key}")
            seen.add(key)
            qx = finite(row.get("qx_per_m"), "qx_per_m")
            qy = finite(row.get("qy_per_m"), "qy_per_m")
            kz2 = finite(row.get("kz2_per_m2"), "kz2_per_m2")
            if not math.isclose(kz2, k0*k0-qx*qx-qy*qy, rel_tol=1e-9, abs_tol=1e-9*k0*k0):
                raise ValueError(f"dispersion mismatch for order {key}")
            rho = math.hypot(qx, qy) / k0
            if rho > 1.2 + 1e-9:
                raise ValueError("order outside prescribed fixed-kin 1.20 cutoff")
            for flag, expected in (("padding11", rho <= 1.1+1e-10), ("padding12", True), ("prop", kz2 > 0)):
                if flag in row and (type(row[flag]) is not bool or row[flag] != expected):
                    raise ValueError(f"{flag} inconsistent with physical fixed-kin dispersion")
            cpair(row.get("Ex"), "Ex")
            cpair(row.get("Ey"), "Ey")
            desc = dict(qx=qx, qy=qy, kz2=kz2, padding11=rho <= 1.1+1e-10)
            if j == 0:
                inventory[key] = desc
                keys.append(key)
            elif key not in inventory or any(not math.isclose(desc[k], inventory[key][k], rel_tol=1e-12, abs_tol=1e-8) for k in ("qx", "qy", "kz2")):
                raise ValueError("order identity/wavevector changed across planes")
        if seen != set(inventory):
            raise ValueError("order inventories differ across planes")
    if (0, 0) not in inventory:
        raise ValueError("incident zeroth order missing")
    if "kin_rad_per_m" in case:
        kin = case["kin_rad_per_m"]
        if not isinstance(kin, list) or len(kin) != 3 or finite(kin[2], "kin_z") >= 0:
            raise ValueError("specified downward kin required")
        if not all(math.isclose(finite(kin[i], "kin"), inventory[(0, 0)][k], abs_tol=.1, rel_tol=1e-10) for i,k in enumerate(("qx","qy"))):
            raise ValueError("zeroth Fourier wavevector is not the specified Bloch incidence")
        kv=np.array([finite(x,"kin") for x in kin])
        if not math.isclose(float(np.linalg.norm(kv)),k0,rel_tol=1e-9):
            raise ValueError("specified incident wavevector is not on vacuum dispersion")
    return z, keys, inventory


def independent_inventory_proof(case: dict, inventory: dict) -> dict:
    """Check every fixed-kin reciprocal order independently, not just common lists."""
    if "reciprocal_basis_per_m" not in case or "kin_rad_per_m" not in case:
        return {"verified":False,"state":"MISSING_RECIPROCAL_OR_INCIDENT_METADATA_NOT_COMPLETE"}
    basis=case["reciprocal_basis_per_m"]
    if not isinstance(basis,list) or len(basis)!=2 or any(not isinstance(b,list) or len(b)!=2 for b in basis):
        raise ValueError("reciprocal_basis_per_m must contain two 2D reciprocal vectors")
    B=np.array([[finite(x,"reciprocal basis") for x in b] for b in basis])
    if abs(float(np.linalg.det(B)))<1e-20 or np.linalg.cond(B)>1e6:
        raise ValueError("singular or ill-conditioned reciprocal basis")
    k0=2*math.pi/case["lambda_m"]
    kin=np.array(case["kin_rad_per_m"][:2],dtype=float)
    inv=np.linalg.inv(B)
    center=kin@inv
    radius=1.2*k0*np.linalg.norm(inv,axis=0)
    lower=np.ceil(-center-radius-1e-9).astype(int)
    upper=np.floor(-center+radius+1e-9).astype(int)
    if np.prod(upper-lower+1)>1_000_000:
        raise ValueError("reciprocal inventory exceeds bounded diagnostic enumeration")
    expected={}
    for n1 in range(int(lower[0]),int(upper[0])+1):
        for n2 in range(int(lower[1]),int(upper[1])+1):
            q=kin+np.array([n1,n2])@B
            if float(np.linalg.norm(q))<=1.2*k0*(1+1e-10):
                expected[(n1,n2)]=q
    if set(expected)!=set(inventory):
        raise ValueError(f"incomplete independent fixed-kin reciprocal inventory: missing={sorted(set(expected)-set(inventory))}, extra={sorted(set(inventory)-set(expected))}")
    for key,q in expected.items():
        if not np.allclose(q,[inventory[key]["qx"],inventory[key]["qy"]],rtol=1e-10,atol=.01):
            raise ValueError(f"fixed-kin reciprocal wavevector mismatch: {key}")
    expected_area=(2*math.pi)**2/abs(float(np.linalg.det(B)))
    if not all(math.isclose(p["area_m2"],expected_area,rel_tol=1e-6) for p in case["planes"]):
        raise ValueError("cut plane area is not the full reciprocal primitive-cell area")
    if "expected_cell_area_m2" in case and not math.isclose(finite(case["expected_cell_area_m2"],"expected_cell_area_m2",True),expected_area,rel_tol=1e-6):
        raise ValueError("geometric cell area disagrees with reciprocal primitive-cell area")
    return {"verified":True,"state":"INDEPENDENT_FIXED_KIN_RECIPROCAL_SET_AND_FULL_CELL_AREA_VERIFIED","expected_total_orders":len(expected),"expected_propagating_orders":sum(d["kz2"]>0 for d in inventory.values()),"expected_cell_area_m2":expected_area,"basis_layout":"two row vectors: q=kin_parallel+[n1,n2]@B"}


def design(z: np.ndarray, kz: complex, zref: float) -> np.ndarray:
    return np.column_stack((np.exp(-1j*kz*(z-zref)), np.exp(1j*kz*(z-zref))))


def maxwell_fields(et: np.ndarray, qx: float, qy: float, signed_kz: complex, k0: float, Z0: float) -> tuple[np.ndarray, np.ndarray]:
    e = np.array([et[0], et[1], -(qx*et[0]+qy*et[1])/signed_kz], dtype=complex)
    h = np.cross(np.array([qx,qy,signed_kz], dtype=complex), e) / (k0*Z0)
    return e, h


def poynting_z(e: np.ndarray, h: np.ndarray, area: float) -> float:
    return float(area*.5*np.real(e[0]*np.conj(h[1])-e[1]*np.conj(h[0])))


def fit_mode(z: np.ndarray, values: np.ndarray, desc: dict, k0: float, Z0: float, area: float) -> dict:
    kz2 = desc["kz2"]
    if abs(kz2) <= 1e-14*k0*k0:
        return {"state":"UNRESOLVED_GRAZING_MODE_NOT_DROPPED", "kz_re_im":[0.,0.], "near_grazing":True}
    kz = complex(math.sqrt(kz2), 0) if kz2 > 0 else complex(0, -math.sqrt(-kz2))
    zref = float(z[2])
    D = design(z, kz, zref)
    coeff, _, rank, _ = np.linalg.lstsq(D, values, rcond=None)
    condition = float(np.linalg.cond(D))
    prediction = D @ coeff
    error2 = float(np.sum(np.abs(prediction-values)**2))
    signal2 = float(np.sum(np.abs(values)**2))
    planes = []
    for j in range(len(z)):
        eup,hup = maxwell_fields(D[j,0]*coeff[0], desc["qx"], desc["qy"], kz, k0, Z0)
        edn,hdn = maxwell_fields(D[j,1]*coeff[1], desc["qx"], desc["qy"], -kz, k0, Z0)
        e,h = eup+edn,hup+hdn
        net = poynting_z(e,h,area)
        up,dn = poynting_z(eup,hup,area),poynting_z(edn,hdn,area)
        measured = poynting_z(values[j],h,area)
        planes.append({"plane_index":j,"fitted_EH_net_W":net,"measured_E_fitted_H_net_W":measured,"up_self_W":up,"down_self_W":dn,"up_down_cross_W":net-up-dn})
    loo = []
    for held in range(len(z)):
        keep = np.arange(len(z)) != held
        c,_,r,_ = np.linalg.lstsq(D[keep],values[keep],rcond=None)
        pred = D[held]@c
        eup,hup = maxwell_fields(D[held,0]*c[0], desc["qx"], desc["qy"], kz, k0, Z0)
        edn,hdn = maxwell_fields(D[held,1]*c[1], desc["qx"], desc["qy"], -kz, k0, Z0)
        loo.append({"held_out_plane_index":held,"design_condition_2":float(np.linalg.cond(D[keep])),"rank":int(r),"held_out_E_error2":float(np.sum(np.abs(pred-values[held])**2)),"held_out_E_signal2":float(np.sum(np.abs(values[held])**2)),"held_out_measured_E_predicted_H_net_W":poynting_z(values[held],hup+hdn,area),"held_out_fitted_EH_net_W":poynting_z(eup+edn,hup+hdn,area),"down_self_W":poynting_z(edn,hdn,area)})
    return {"state":"FITTED_NOT_ACCEPTED","kz_re_im":pair(kz),"reference_z_m":zref,"near_grazing":abs(kz)/k0<.01,"design_condition_2":condition,"rank":int(rank),"full_fit_E_error2":error2,"full_fit_E_signal2":signal2,"relative_L2_E_residual":math.sqrt(error2/max(signal2,1e-300)),"up_Ex_Ey_re_im":[pair(x) for x in coeff[0]],"down_Ex_Ey_re_im":[pair(x) for x in coeff[1]],"planes":planes,"leave_one_out":loo,"condition_passed":bool(condition<=CONDITION_LIMIT and all(x["design_condition_2"]<=CONDITION_LIMIT and x["rank"]==2 for x in loo) and rank==2)}


def analyze_case(case: dict) -> dict:
    z, keys, inventory = validate_case(case)
    inventory_proof=independent_inventory_proof(case,inventory)
    k0 = 2*math.pi/case["lambda_m"]
    area = case["planes"][0]["area_m2"]
    by_plane = [{(r["n1"],r["n2"]):r for r in p["orders"]} for p in case["planes"]]
    modes=[]
    for key in keys:
        vals=np.array([[cpair(p[key]["Ex"],"Ex"),cpair(p[key]["Ey"],"Ey")] for p in by_plane])
        result=fit_mode(z,vals,inventory[key],k0,case["Z0"],area)
        modes.append({"n1":key[0],"n2":key[1],**inventory[key],**result})
    unresolved=[(m["n1"],m["n2"]) for m in modes if "planes" not in m]
    shells=[]
    for cutoff in (1.1,1.2):
        selected=[m for m in modes if cutoff==1.2 or m["padding11"]]
        complete=all("planes" in m for m in selected)
        shell={"fixed_transverse_cutoff_over_k0":cutoff,"mode_count":len(selected),"complete_exported_inventory_fit":complete,"complete_no_mode_silently_removed":complete and inventory_proof["verified"]}
        if complete:
            err=sum(m["full_fit_E_error2"] for m in selected)
            signal=sum(m["full_fit_E_signal2"] for m in selected)
            full_res=math.sqrt(err/max(signal,1e-300))
            net=[sum(m["planes"][j]["fitted_EH_net_W"] for m in selected) for j in range(5)]
            measured=[sum(m["planes"][j]["measured_E_fitted_H_net_W"] for m in selected) for j in range(5)]
            held=[sum(m["leave_one_out"][j]["held_out_measured_E_predicted_H_net_W"] for m in selected) for j in range(5)]
            loo_res=[math.sqrt(sum(m["leave_one_out"][j]["held_out_E_error2"] for m in selected)/max(sum(m["leave_one_out"][j]["held_out_E_signal2"] for m in selected),1e-300)) for j in range(5)]
            incoming=sum(-m["planes"][2]["down_self_W"] for m in selected if m["kz2"]>0)
            source=next(m for m in selected if (m["n1"],m["n2"])==(0,0))
            zero_incoming=-source["planes"][2]["down_self_W"]
            loo_zero_incoming=[-r["down_self_W"] for r in source["leave_one_out"]]
            nonzero_incoming=sum(-m["planes"][2]["down_self_W"] for m in selected if m["kz2"]>0 and (m["n1"],m["n2"])!=(0,0))
            out=sum(m["planes"][2]["up_self_W"] for m in selected if m["kz2"]>0)
            shell.update({"full_fit_relative_L2_E_residual":full_res,"LOO_held_out_relative_L2_E_residual":loo_res,"all_modes_condition_passed":all(m["condition_passed"] for m in selected),"field_prediction_checks_passed":full_res<=RESIDUAL_LIMIT and max(loo_res)<=RESIDUAL_LIMIT and all(m["condition_passed"] for m in selected),"all_five_fitted_EH_net_W":net,"all_five_measured_E_fitted_H_net_W":measured,"all_five_LOO_measured_E_predicted_H_net_W":held,"all_five_A_measured_E_fitted_H_over_specified_Pin":[-v/case["nominal_Pin_W"] for v in measured],"all_five_A_LOO_over_specified_Pin":[-v/case["nominal_Pin_W"] for v in held],"total_propagating_incoming_W":incoming,"nonzero_order_incoming_W_preserved":nonzero_incoming,"nonzero_incoming_is_diagnostic_not_forced_zero":True,"propagating_outgoing_W":out,"incoming_error_vs_specified_Pin_W":abs(incoming-case["nominal_Pin_W"]),"specified_1W_incident_power_check_passed":abs(incoming-case["nominal_Pin_W"])<=ABS_GATE,"evanescent_up_down_cross_W":sum(m["planes"][2]["up_down_cross_W"] for m in selected if m["kz2"]<0)})
            shell.update({"zeroth_order_incident_W":zero_incoming,"all_LOO_zeroth_order_incident_W":loo_zero_incoming,"zeroth_and_all_LOO_specified_1W_check_passed":abs(zero_incoming-case["nominal_Pin_W"])<=ABS_GATE and max(abs(v-case["nominal_Pin_W"]) for v in loo_zero_incoming)<=ABS_GATE,"tiny_mode_relative_residuals_are_diagnostics_not_individual_hard_gates":True})
        shells.append(shell)
    cutoff_comparison={"complete":all("all_five_fitted_EH_net_W" in s for s in shells)}
    if cutoff_comparison["complete"]:
        cutoff_comparison.update({"cutoff12_minus_cutoff11_fitted_EH_W":[b-a for a,b in zip(shells[0]["all_five_fitted_EH_net_W"],shells[1]["all_five_fitted_EH_net_W"])],"cutoff12_minus_cutoff11_measured_E_fitted_H_W":[b-a for a,b in zip(shells[0]["all_five_measured_E_fitted_H_net_W"],shells[1]["all_five_measured_E_fitted_H_net_W"])],"no_favorable_cutoff_selected":True})
    return {"schema":SCHEMA,"case_id":case["case_id"],"state":"SAVED_FIELD_RECONSTRUCTION_DIAGNOSTIC_ONLY" if not unresolved and inventory_proof["verified"] else "INCOMPLETE_INVENTORY_GRAZING_OR_UNRESOLVED_DIAGNOSTIC_ONLY","solver_executed":False,"formal_qualified":False,"raw_H_historical_flags_changed":False,"fit_uses_S_parameters_or_losses":False,"normalization":"specified_1W_unchanged_not_fitted_Pin","nominal_Pin_W":case["nominal_Pin_W"],"independent_inventory_and_area_proof":inventory_proof,"phase_convention":"exp(+i*omega*t);up=exp(-i*kz*dz);Re(kz)>=0,Im(kz)<=0","cap_fractions":list(FRACTIONS),"unresolved_orders_not_dropped":unresolved,"preserved_absolute_gate":ABS_GATE,"preserved_relative_gate":REL_GATE,"diagnostic_field_residual_limit":RESIDUAL_LIMIT,"condition_limit":CONDITION_LIMIT,"reconstructed_plane_invariance_is_algebraic_not_validation":True,"native_completeness_receipt_must_be_verified_separately_by_controller":True,"selection_policy":"all_five_planes_equal_weights_and_both_fixed_cutoffs_no_favorable_plane_or_shell_selection","mode_results":modes,"fixed_nested_cutoffs":shells,"fixed_cutoff_truncation_comparison":cutoff_comparison}


def closure_diagnostic(result: dict, total_loss: float) -> dict:
    """Optional after-fit comparison; cannot alter input, fit, or normalization."""
    loss=finite(total_loss,"postfit total_loss")
    rows=[]
    for shell in result["fixed_nested_cutoffs"]:
        if not shell["complete_no_mode_silently_removed"]:
            rows.append({"cutoff":shell["fixed_transverse_cutoff_over_k0"],"comparison_complete":False})
            continue
        for kind,key in (("measured_E_fitted_H","all_five_A_measured_E_fitted_H_over_specified_Pin"),("LOO_measured_E_predicted_H","all_five_A_LOO_over_specified_Pin")):
            differences=[abs(a-loss) for a in shell[key]]
            relative=[d/max(abs(a),1e-12) for d,a in zip(differences,shell[key])]
            rows.append({"cutoff":shell["fixed_transverse_cutoff_over_k0"],"kind":kind,"all_five_absolute_mismatch":differences,"all_five_relative_mismatch":relative,"all_five_preserved_closure_gates_passed":max(differences)<=ABS_GATE and max(relative)<=REL_GATE,"formal_qualified":False})
    return {"comparison_after_fit_only":True,"loss_used_for_fit_or_normalization":False,"rows":rows}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    results={"schema":SCHEMA,"cases":[analyze_case(c) for c in parse_export(args.input.read_text(encoding="utf-8-sig"))]}
    content=json.dumps(results,indent=2,allow_nan=False)+"\n"
    if args.output:
        with args.output.open("x",encoding="utf-8") as f:
            f.write(content)
    else:
        print(content,end="")


if __name__=="__main__":
    main()

