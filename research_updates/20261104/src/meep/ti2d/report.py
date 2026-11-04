"""Offline reports and figures from saved solver evidence; never launches Meep."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from .queue import (SolverLedger, aggregate_gate, atomic_json, read_json,
                    ledger_path_for_manifest, reusable_result)


def _number(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


def collect(manifest_path):
    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent
    manifest = read_json(manifest_path)
    status = read_json(root / "status.json") if (root / "status.json").exists() else {}
    lock = read_json(root / "lock_manifest.json") if (root / "lock_manifest.json").exists() else {}
    rows, results = [], {}
    for task in manifest["tasks"]:
        task_status = status.get("tasks",{}).get(task["id"],{})
        result_path = task_status.get("result_path")
        if not result_path:
            attempts = sorted((root / "tasks" / task["id"]).glob("attempt_*/result.json"))
            result_path = str(attempts[-1]) if attempts else None
        result = {}
        if result_path:
            try:
                result = read_json(result_path)
            except (OSError,ValueError):
                result = {"status":"CORRUPT_OUTPUT"}
            if result.get("status") == "QUALIFIED":
                fingerprint = lock.get("case_sha256",{}).get(task["id"])
                verified = reusable_result(Path(result_path).parent,fingerprint) if fingerprint else None
                if verified is None:
                    result = {**result,"status":"UNVERIFIED_EVIDENCE",
                              "reason":"Missing seal, fingerprint mismatch, or corrupt evidence"}
                    task_status = {**task_status,"status":"UNVERIFIED_EVIDENCE"}
        results[task["id"]] = result
        metrics = result.get("metrics",{})
        case = read_json(root / task["case_path"]) if not Path(task["case_path"]).is_absolute() else read_json(task["case_path"])
        geometry = case.get("geometry",{})
        row = {"id":task["id"],"phase":task["phase"],"geometry_id":task.get("geometry_id"),
               "gate_role":task.get("gate_role"),
               "resolution":task.get("resolution"),"polarization":task.get("polarization"),
               "emission_theta_deg":task.get("emission_theta_deg"),
               "status":task_status.get("status",result.get("status","NOT_RUN")),
               "result_path":result_path or "", "material":case.get("material_path",""),
               "oxide_layer_in_model":False,
               "elapsed_s":result.get("queue_process",{}).get("elapsed_s",result.get("elapsed_s")),
               "mpi_ranks":result.get("queue_process",{}).get("mpi_ranks"),
               "peak_process_tree_rss_bytes":result.get("queue_process",{}).get("peak_process_tree_rss_bytes")}
        for key in ("R","T","A_flux","A_vol","order_R","order_T","reference_pseudo_reflection"):
            row[key] = metrics.get(key)
        for key in ("period_um","surface_fill","axis_length_um","tilt_deg","thickness_um"):
            row[key] = geometry.get(key,case.get(key))
        bounds = result.get("numerical_error_bounds",result.get("numeric_error_bounds",{}))
        row["A_error_bound"] = bounds.get("A") if bounds.get("validated") else None
        rows.append(row)
    return manifest,status,rows,results


def write_csv(path, rows, default_fields=("status",)):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row)) or list(default_fields)
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer = csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def paired_metrics(rows, gate_passed):
    pairs = []
    for geometry_id in dict.fromkeys(r["geometry_id"] for r in rows if r["phase"] == "pilot"):
        for polarization in ("s","p"):
            pair = {r["emission_theta_deg"]:r for r in rows if r["phase"] == "pilot"
                    and r["geometry_id"] == geometry_id and r["polarization"] == polarization}
            minus,plus = pair.get(-30),pair.get(30)
            row = {"geometry_id":geometry_id,"polarization":polarization,
                   "status":"UNVERIFIED","emissivity_minus30":None,"emissivity_plus30":None,
                   "angle_difference":None,"paired_error_bound":None,"normalized_contrast":None,
                   "interpretation":"Not calculated or not qualified"}
            if gate_passed and minus and plus and all(r["status"] == "QUALIFIED" for r in (minus,plus)):
                a,b = minus["A_flux"],plus["A_flux"]
                row.update(emissivity_minus30=a,emissivity_plus30=b,angle_difference=b-a,
                           status="EMISSIVITY_QUALIFIED_ERROR_BOUND_UNVERIFIED")
                errors = [minus["A_error_bound"],plus["A_error_bound"]]
                if all(_number(e) and e >= 0 for e in errors):
                    bound = sum(errors)
                    row.update(paired_error_bound=bound,status="QUALIFIED",
                        interpretation="Resolved" if abs(b-a)>bound else "Current precision cannot resolve the difference")
                    if abs(a+b) > 5*bound and a+b > 0:
                        row["normalized_contrast"] = (b-a)/(b+a)
                    else:
                        row["interpretation"] += "; contrast denominator below 5x paired error bound"
                else:
                    row["interpretation"] = "No validated numerical error bound; contrast unverified"
            pairs.append(row)
    return pairs


def sensitivities(rows,gate_passed):
    parameter = {"period_high":"period_um","duty_high":"surface_fill",
                 "length_high":"axis_length_um","tilt_high":"tilt_deg"}
    answer = []
    baseline = {(r["polarization"],r["emission_theta_deg"]):r for r in rows
                if r["phase"] == "pilot" and r["geometry_id"] == "baseline"}
    for row in rows:
        if row["phase"] != "pilot" or row["geometry_id"] not in parameter:
            continue
        base = baseline.get((row["polarization"],row["emission_theta_deg"]))
        p = parameter[row["geometry_id"]]
        item = {"geometry_id":row["geometry_id"],"parameter":p,"polarization":row["polarization"],
                "emission_theta_deg":row["emission_theta_deg"],"status":"UNVERIFIED",
                "delta_parameter":None,"delta_emissivity":None,"finite_difference":None,
                "effect_error_bound":None,"sensitivity_error_bound":None,"resolved":None}
        if (gate_passed and base and row["status"] == base["status"] == "QUALIFIED"
                and _number(row[p]) and _number(base[p]) and row[p] != base[p]):
            dx,dy = row[p]-base[p],row["A_flux"]-base["A_flux"]
            item.update(delta_parameter=dx,delta_emissivity=dy,finite_difference=dy/dx,
                        status="LOCAL_DIFFERENCE_ERROR_BOUND_UNVERIFIED")
            if all(_number(r["A_error_bound"]) and r["A_error_bound"]>=0 for r in (row,base)):
                error = row["A_error_bound"]+base["A_error_bound"]
                item.update(effect_error_bound=error,sensitivity_error_bound=error/abs(dx),
                            resolved=abs(dy)>error,status="QUALIFIED")
        answer.append(item)
    return answer


def _placeholder(ax,title,reason="Not run / unqualified"):
    ax.set_title(title)
    ax.text(.5,.5,reason,ha="center",va="center",transform=ax.transAxes,wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])


def make_plots(output,manifest,rows,results,pairs,sensitivity,gate):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True,exist_ok=True)
    paths = []
    # Geometry is sourced from validated case definitions; no inferred 3D shape.
    fig,axes = plt.subplots(1,5,figsize=(17,4),constrained_layout=True)
    geometries = []
    for task in manifest["tasks"]:
        if task["phase"] == "pilot" and task.get("geometry_id") not in [x[0] for x in geometries]:
            geometries.append((task.get("geometry_id"),task))
    from .geometry import validate_geometry
    for ax,item in zip(axes,geometries):
        name,task = item
        path = Path(task["case_path"])
        if not path.is_absolute():
            path = Path(manifest["_run_dir"]) / path
        case = read_json(path)
        try:
            geom = validate_geometry(case["geometry"])
            vertices = geom["vertices_x_depth"]
            period = case["geometry"]["period_um"]
            for shift in (-period,0,period):
                xs = [v[0]+shift for v in vertices]+[vertices[0][0]+shift]
                ys = [v[1] for v in vertices]+[vertices[0][1]]
                ax.fill(xs,ys,color="white",edgecolor="black")
            ax.set_facecolor("#b2b6bf")
            ax.set_xlim(-period/2,period/2)
            ax.set_ylim(60,0)
            ax.set_aspect("equal")
            ax.set_xlabel("x (um)")
            ax.set_ylabel("Depth (um)")
            ax.set_title(name)
        except (KeyError,ValueError,TypeError) as exc:
            _placeholder(ax,name,str(exc))
    path = output / "geometry.png"; fig.savefig(path,dpi=150); plt.close(fig); paths.append(str(path))
    fig,axes = plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    plotted = False
    for role in ("flat","baseline"):
        for pol in ("s","p"):
            selected = [r for r in rows if r["phase"]=="validation" and r["gate_role"]==role
                        and r["polarization"]==pol and r["emission_theta_deg"]==30 and _number(r["R"])]
            selected.sort(key=lambda r:r["resolution"])
            if selected:
                axes[0].plot([r["resolution"] for r in selected],[r["R"] for r in selected],"o-",label=f"{role} {pol}")
                plotted=True
    if plotted:
        axes[0].set(xlabel="Resolution (pixels/um)",ylabel="R (unclipped)",title="Grid evidence, not an assumed error bound")
        axes[0].legend(fontsize=8)
    else:
        _placeholder(axes[0],"Grid convergence")
    selected = [r for r in rows if _number(r["A_vol"]) and _number(r["A_flux"])]
    if selected:
        axes[1].plot(range(len(selected)),[abs(r["A_vol"]-r["A_flux"]) for r in selected],"o",label="|Avol-Aflux|")
        axes[1].axhline(.01,color="red",linestyle="--",label="Gate 0.01")
        axes[1].set(xlabel="Saved condition index",ylabel="Absolute error",title="Independent absorption check")
        axes[1].legend()
    else:
        _placeholder(axes[1],"Energy diagnostic errors")
    path=output/"convergence_energy.png";fig.savefig(path,dpi=150);plt.close(fig);paths.append(str(path))
    fig,axes=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for ax,pol in zip(axes,("s","p")):
        selected=[p for p in pairs if p["polarization"]==pol and _number(p["emissivity_plus30"])]
        if selected:
            x=list(range(len(selected)))
            ax.plot(x,[p["emissivity_minus30"] for p in selected],"o-",label="-30 deg")
            ax.plot(x,[p["emissivity_plus30"] for p in selected],"s-",label="+30 deg")
            ax.set_xticks(x,[p["geometry_id"] for p in selected],rotation=25,ha="right")
            ax.set(ylabel="Directional emissivity = qualified absorptivity",title=f"{pol} polarization, specified Ti model")
            ax.legend()
        else:
            _placeholder(ax,f"Paired emissivity: {pol}","Qualification gate has not passed" if not gate["passed"] else "Pilot not run")
    path=output/"paired_emissivity.png";fig.savefig(path,dpi=150);plt.close(fig);paths.append(str(path))
    fig,ax=plt.subplots(figsize=(10,4),constrained_layout=True)
    selected=[s for s in sensitivity if _number(s["delta_emissivity"])]
    if selected:
        ax.bar(range(len(selected)),[s["delta_emissivity"] for s in selected])
        ax.set_xticks(range(len(selected)),[f"{s['geometry_id']} {s['polarization']} {s['emission_theta_deg']}" for s in selected],rotation=55,ha="right")
        ax.set(ylabel="Emissivity change from baseline",title="Local machining-parameter differences (coupled definitions)")
    else:
        _placeholder(ax,"Local parameter sensitivity")
    path=output/"parameter_sensitivity.png";fig.savefig(path,dpi=150);plt.close(fig);paths.append(str(path))
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    plotted=False
    for task_id,result in results.items():
        trace=result.get("stop",{}).get("trace",[])
        if not isinstance(trace,list) or not trace:
            continue
        valid=[s for s in trace if _number(s.get("time"))]
        if not valid:
            continue
        intensity_points=[]
        drift_points=[]
        for sample in valid:
            ratios=[sample.get("intensity_ratio"),sample.get("envelope_ratio")]
            probes=sample.get("probe_ratio",{})
            if isinstance(probes,dict): ratios.extend(probes.values())
            ratios=[r for r in ratios if _number(r)]
            if ratios and max(ratios)>0:
                intensity_points.append((sample["time"],max(ratios)))
            drift=[sample.get("drift"),sample.get("window_R_range"),sample.get("window_T_range")]
            drift=[d for d in drift if _number(d)]
            if drift and max(drift)>0:
                drift_points.append((sample["time"],max(drift)))
        for ax,points in zip(axes,(intensity_points,drift_points)):
            if points:
                ax.semilogy([p[0] for p in points],[p[1] for p in points],label=task_id)
                plotted=True
        windows = result.get("stop",{}).get("windows",[])
        points = [(w.get("end"),max(w["R_range"],w["T_range"])) for w in windows
                  if _number(w.get("end")) and _number(w.get("R_range")) and _number(w.get("T_range"))]
        points = [(t,d) for t,d in points if d>0]
        if points:
            axes[1].semilogy([p[0] for p in points],[p[1] for p in points],label=task_id)
            plotted=True
    if not plotted:
        for ax,title in zip(axes,("Historical-peak normalized |E|^2 decay","R/T window absolute drift")):
            _placeholder(ax,title,"No compatible saved stop trace yet; inspect task stop evidence")
    else:
        axes[0].axhline(1e-4,color="red",linestyle="--")
        axes[1].axhline(1e-3,color="red",linestyle="--")
        for ax in axes:
            ax.set_xlabel("Meep time")
            if ax.lines: ax.legend(fontsize=5)
    path=output/"stopping.png";fig.savefig(path,dpi=150);plt.close(fig);paths.append(str(path))
    return paths


def generate_report(manifest_path):
    manifest_path=Path(manifest_path).resolve()
    root=manifest_path.parent
    manifest,status,rows,results=collect(manifest_path)
    manifest["_run_dir"]=str(root)
    output=root/"reports"
    output.mkdir(parents=True,exist_ok=True)
    gate=aggregate_gate(manifest,results)
    pairs=paired_metrics(rows,gate["passed"])
    sensitivity=sensitivities(rows,gate["passed"])
    from .uncertainty import derive_error_evidence
    error_evidence=derive_error_evidence(manifest,results)
    write_csv(output/"conditions.csv",rows)
    write_csv(output/"paired_emissivity.csv",pairs)
    write_csv(output/"local_sensitivity.csv",sensitivity)
    write_csv(output/"numerical_error_evidence.csv",error_evidence)
    usage=SolverLedger(ledger_path_for_manifest(manifest,root)).as_dict()
    report={"state":status.get("state","NOT_STARTED"),"scope":"Specified Ti optical model; 2D slanted slots, not 3D circular holes or high-temperature measured samples",
            "oxide_layer_in_model":False,"sample_oxidation_state":"unknown",
            "gate":gate,"usage":usage,"conditions":rows,"pairs":pairs,"sensitivities":sensitivity,
            "concurrency_probe":status.get("concurrency_probe",{"status":"NOT_RUN"}),
            "numerical_error_evidence":error_evidence,
            "ai_usage":status.get("ai_usage",{"actual_usage_available":False,"rmb_cost":None}),
            "numerical_error_note":"Grid differences are observations, not automatic rigorous error bounds. Contrast and resolved-effect claims require explicit validated error bounds."}
    try:
        report["figures"]=make_plots(output/"figures",manifest,rows,results,pairs,sensitivity,gate)
    except Exception as exc:
        # Preserve solver evidence and CSV even if the optional plotting dependency fails.
        report["figure_error"]=str(exc)
    atomic_json(output/"summary.json",report)
    qualified=sum(r["status"]=="QUALIFIED" for r in rows if r["phase"]=="pilot")
    failure_lines = ['- ' + x for x in gate['failures'][:5]]
    if len(gate['failures']) > 5:
        failure_lines.append(f"- {len(gate['failures'])-5} additional failed/not-run checks; full list in summary.json.")
    content=f"""# Bounded Ti 2D stage handoff

State: **{report['state']}**. Pilot conditions qualified: **{qualified}/20**.
Qualification gate: **{'PASS' if gate['passed'] else 'NOT PASSED'}**.

Scope: specified Ordal Route A Ti optical model, 2D slanted slots, 10.5 um, +/-30 deg, s/p only. No oxide layer; real-sample oxidation and high-temperature validity remain unknown.

## Resource use

- Active-solver interval union: {usage['active_solver_seconds']/3600:.4f} h / 48 h.
- MPI rank-hours: {usage['rank_hours']:.4f}; calendar span: {usage['calendar_span_seconds']/3600:.4f} h.
- No AI calls are made by this queue. AI usage and RMB costs are unavailable unless an actual account bill is supplied; no price is invented.
- Concurrency: one locked two-job test uses required validation conditions, each limited to half the group cap of 70% container-available memory. More than 30% slowdown or incomparable cache state returns to serial; heterogeneous later cases conservatively stay serial. Probe status: {report['concurrency_probe'].get('status','NOT_RUN')}. Restarting an interrupted condition is not time-step checkpoint recovery.

## Interpretation and unfinished work

{chr(10).join(failure_lines) or '- Bounded qualification passed; consult per-condition status for remaining pilot work.'}

- Preserve unclipped R, T, A_flux, A_vol. Algebraic R+T+A_flux=1 is not independent conservation evidence.
- A numeric value in conditions.csv does not imply qualification; use its status and this gate.
- Fixed surface fill with changed period changes opening width; fixed opening with changed tilt changes normal width. Only local finite differences are reported.
- A missing validated error bound means contrast is unverified, not zero uncertainty. Effects within paired error bounds are unresolved at current precision.
- No formal peak angle, angular spectrum, FWHM, optimal design, 3D or high-temperature conclusion.

## Evidence and next action

- `conditions.csv`, `paired_emissivity.csv`, `local_sensitivity.csv`, `numerical_error_evidence.csv`, `summary.json`, `figures/`. Observed grid/stop spreads are evidence, not a rigorous universal error bound.
- Original arrays and stop curves remain in task attempt directories; previous files are not overwritten.
- Inspect a failed condition's result.json and short log tail before authorizing any targeted repair or rerun. No automatic scan expansion.
- For status use `python -B -m ti2d.cli status --manifest {manifest_path}`; report regeneration is offline.
- Resume accepts exact fingerprint-matched qualified results. Unclosed intervals, corrupt evidence, and earlier unqualified attempts require review.
"""
    (output/"handoff.md").write_text(content,encoding="utf-8")
    allowed=["ti2d/*.py","tests/*.py",str(manifest_path),str(output/"handoff.md"),
             str(output/"summary.json"),str(output/"conditions.csv"),str(output/"paired_emissivity.csv"),
             str(output/"local_sensitivity.csv"),str(output/"numerical_error_evidence.csv"),str(output/"figures/*.png")]
    # This is an upload candidate list, not a claim a secret scanner has run.
    (output/"github_upload_candidates.txt").write_text("# Review and scan for secrets before upload; never upload credentials, sessions, caches, or large field files.\n"+"\n".join(allowed)+"\n",encoding="utf-8")
    return report
