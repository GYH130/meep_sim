"""D15c_build_d15c_tables.py — 从 D15b 结果 JSON 生成复核表（不跑 FDTD）

输入：``results/tables/D15b_*.json``（每个条件一个文件）
输出：``results/tables/D15c_*.csv`` 与 ``results/reports/D15c_summary.json``

产出三张表 + 一份汇总：

1. ``D15c_condition_table.csv``  逐条件一行：几何、R/T/A/A_vol、收敛、
   源后窗口、分阶闭合、资源、`quality_status`、`fingerprint`。
2. ``D15c_order_table.csv``      逐阶一行：m、kx、ky、出射角、份额（反射/透射）。
3. ``D15c_convergence_table.csv`` 同一条件在不同 resolution 下的 R/T/A/A_vol。
4. ``D15c_qualification_coverage.csv`` 二维资格逐项覆盖表（哪些项已做、
   由哪个条件覆盖、证据文件）。

只做**读与汇总**，不修改任何已有 JSON，也不重跑仿真。
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

TABLE_DIR = PROJECT_ROOT / "results" / "tables"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"

#: 二维资格项 → 判定依据。'done' 由具体条件的 quality_checks 决定。
QUALIFICATION_ITEMS = [
    ("geom_qualified", "几何合格（壁厚/底厚/不贯穿）",
     "slot_geometry.qualification + tests/test_slot_geometry.py"),
    ("q_lt_1", "色散稳定 q<1 且 q>=1 时拒绝",
     "D15b dispersion_q 拒绝路径（cal4_r20 日志可见 r24/C0.5 记录）"),
    ("flux_sign_convention", "通量符号约定（R=refl/|input|）",
     "D15b 与 D04 一致；平面基准 R 与 Fresnel 对比"),
    ("all_orders_enumerated", "按实际 kx 与周期枚举全部传播阶",
     "enumerate_orders + tests/test_d15b_order_projection.py"),
    ("order_decomposition", "分阶功率与整面通量闭合",
     "D15b order_decomposition.closure（绝对标定 P/n²）"),
    ("volume_absorption", "独立体吸收 A_vol 与 A_flux 一致性",
     "D15b 每对 Yee 分量分别积分 + 金属/空气分区"),
    ("post_source_window", "源后观察窗口充分",
     "D15b post_source_stability.post_source_ok"),
    ("multi_probe", "多探针场强包络稳定",
     "D15b post_source_stability.probe_envelope_segments"),
    ("ref_pseudo_reflection", "参考伪反射（同布局平面基准）",
     "flat 条件：R 应与解析 Fresnel 比较"),
    ("mesh_convergence", "网格收敛（r16/r20/r24）",
     "D15c_convergence_table.csv"),
    ("mirror_symmetry", "镜像/正负角控制 (±α, ±θ)",
     "待跑（本阶段未完成）"),
    ("s_polarization", "s 偏振",
     "待跑（本阶段未完成）"),
    ("directionality", "+30° 与 −30° 观察不要求相等",
     "方向性信号，非一致性检查"),
]


def _load_all() -> list[dict]:
    out = []
    for p in sorted(TABLE_DIR.glob("D15b_*.json")):
        try:
            out.append((p, json.loads(p.read_text())))
        except Exception as exc:
            print(f"[warn] cannot read {p.name}: {exc!r}", file=sys.stderr)
    return out


def _flat_fresnel_R(rec: dict) -> float | None:
    """平面基准的解析 Fresnel R（用于参考伪反射判定）。"""
    try:
        import cmath
        import meep as mp
        from src.materials import get_ti_medium
        wl = rec["wavelength_um"]
        ti = get_ti_medium(wl - 1.0, wl + 1.0)
        eps = ti.epsilon(1.0 / wl)[0][0]
        n = cmath.sqrt(eps)
        return float(abs((1.0 - n) / (1.0 + n)) ** 2)
    except Exception:
        return None


def build_condition_table(records) -> Path:
    out = TABLE_DIR / "D15c_condition_table.csv"
    cols = [
        "file", "fingerprint", "case", "resolution", "courant", "polarization",
        "theta_deg", "tilt_deg", "R", "T", "A_flux", "A_vol",
        "volume_flux_abs_difference", "order_sum_over_plane_flux",
        "n_orders", "mr_max", "post_source_window_um", "n_steps", "wall_total_s",
        "mpi_processes", "peak_memory_mib_max_rank", "quality_status",
        "quality_missing",
    ]
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for p, d in records:
            od = d.get("order_decomposition") or {}
            clo = od.get("closure") or {}
            qm = d.get("quality_missing")
            w.writerow([
                p.name, d.get("fingerprint"), d.get("case"),
                d.get("resolution"), d.get("courant"), d.get("polarization"),
                d.get("emission_theta_deg"), d.get("slot_tilt_angle_deg"),
                _f(d.get("R")), _f(d.get("T")), _f(d.get("A_flux")),
                _f(d.get("A_vol")), _f(d.get("volume_flux_abs_difference")),
                _f(clo.get("refl_order_over_flux")),
                od.get("n_orders"),
                (od.get("orthogonality_self_check") or {}).get("max_cross_overlap"),
                _f(d.get("post_source_window_um")), d.get("n_steps"),
                _f(d.get("wall_total_s")), d.get("mpi_processes"),
                _f(d.get("peak_memory_mib_max_rank")), d.get("quality_status"),
                ";".join(qm) if isinstance(qm, list) else qm,
            ])
    return out


def _f(v):
    if v is None:
        return None
    try:
        return f"{float(v):.9g}"
    except Exception:
        return v


def build_order_table(records) -> Path:
    out = TABLE_DIR / "D15c_order_table.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "plane", "m", "kx_m", "ky_m", "theta_out_deg",
                    "flux_meep_units", "share"])
        for p, d in records:
            od = d.get("order_decomposition") or {}
            for plane, key_rows, key_frac in (
                    ("refl", "rows_refl", "refl_order_fractions"),
                    ("trans", "rows_trans", "trans_order_fractions")):
                rows = od.get(key_rows) or []
                fracs = od.get(key_frac) or {}
                for r in rows:
                    w.writerow([
                        p.name, plane, r.get("m"), _f(r.get("kx_m")),
                        _f(r.get("ky_m")), _f(r.get("theta_out_deg")),
                        _f(r.get("flux_meep_units")),
                        _f(fracs.get(str(r.get("m")))),
                    ])
    return out


def build_convergence_table(records) -> Path:
    out = TABLE_DIR / "D15c_convergence_table.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "polarization", "theta_deg", "resolution",
                    "courant", "R", "T", "A_flux", "A_vol", "A_vol_minus_A_flux",
                    "post_source_window_um", "quality_status", "file"])
        for p, d in records:
            # 只收"真正的条件"：标定/调试运行（max_steps>0）不进收敛表
            if d.get("max_steps"):
                continue
            if d.get("quality_status") not in ("valid",) and "invalid" not in str(
                    d.get("quality_status")):
                continue
            w.writerow([
                d.get("case"), d.get("polarization"), d.get("emission_theta_deg"),
                d.get("resolution"), d.get("courant"),
                _f(d.get("R")), _f(d.get("T")), _f(d.get("A_flux")),
                _f(d.get("A_vol")), _f(d.get("volume_flux_abs_difference")),
                _f(d.get("post_source_window_um")), d.get("quality_status"),
                p.name,
            ])
    return out


def build_coverage_table(records) -> tuple[Path, dict]:
    out = TABLE_DIR / "D15c_qualification_coverage.csv"
    valid = [(p, d) for p, d in records if d.get("quality_status") == "valid"]
    case_by = {}
    for p, d in valid:
        case_by.setdefault(d.get("case"), []).append((p, d))

    def any_valid(case) -> bool:
        return bool(case_by.get(case))

    res_set = sorted({d.get("resolution") for _, d in valid})
    checks_by_case = {}
    for case, items in case_by.items():
        agg = {}
        for _p, d in items:
            for k, v in (d.get("quality_checks") or {}).items():
                agg[k] = agg.get(k, bool(v)) and bool(v)
        checks_by_case[case] = agg

    struct_ok = checks_by_case.get("structure", {})
    flat_ok = checks_by_case.get("flat", {})

    # 参考伪反射：读已存的 Fresnel 审计表（由 D15c_flat_fresnel_audit 生成）
    flat_gate_pass = None
    try:
        import csv as _csv
        audit = TABLE_DIR / "D15c_flat_fresnel_audit.csv"
        if audit.exists():
            for row in _csv.DictReader(audit.open()):
                if row["quantity"] == "gate_pass":
                    flat_gate_pass = (row["value"].strip().lower() == "true")
    except Exception:
        flat_gate_pass = None

    rows = []
    for key, desc, evidence in QUALIFICATION_ITEMS:
        if key == "geom_qualified":
            status = "done" if struct_ok.get("geometry_qualified") else "todo"
        elif key == "q_lt_1":
            status = "done" if struct_ok.get("dispersion_q_lt_1") else "todo"
        elif key == "flux_sign_convention":
            status = "done" if (any_valid("flat") and any_valid("structure")) else "todo"
        elif key == "all_orders_enumerated":
            status = "done" if (struct_ok.get("all_propagating_orders_enumerated")
                               and (struct_ok.get("order_decomposition_done")
                                    or flat_ok.get("order_decomposition_done"))) else "todo"
        elif key == "order_decomposition":
            # 分阶**份额**已可信（平面基准 m=0 占 0.9999）；但"阶和=整面通量"
            # 在相干相消时不成立，本阶段**不作为**通过判据。
            status = "done_shares_only" if (struct_ok.get("order_decomposition_done")
                                            and flat_ok.get("order_decomposition_done")) else "partial"
        elif key == "volume_absorption":
            status = "done" if struct_ok.get("volume_absorption_done") else "todo"
        elif key == "post_source_window":
            status = "done" if struct_ok.get("post_source_window_sufficient") else "todo"
        elif key == "multi_probe":
            status = "done" if struct_ok.get("multi_probe_recorded") else "todo"
        elif key == "ref_pseudo_reflection":
            # "做了"与"通过"是两件事：这里读实测的判定结果，不靠有无文件
            status = "done" if flat_gate_pass else (
                "FAILED" if flat_gate_pass is False else "todo")
        elif key == "mesh_convergence":
            status = "done" if len([r for r in res_set if r]) >= 3 else \
                ("partial" if len([r for r in res_set if r]) >= 2 else "todo")
        elif key == "directionality":
            status = "n/a（非一致性检查）"
        else:
            status = "todo"
        rows.append({"item": key, "description": desc, "status": status,
                     "evidence": evidence})

    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "description", "status", "evidence"])
        w.writeheader()
        w.writerows(rows)

    summary = {
        "n_conditions_total": len(records),
        "n_conditions_valid": len(valid),
        "resolutions_covered": res_set,
        "cases_covered": sorted(case_by.keys()),
        "coverage": {r["item"]: r["status"] for r in rows},
    }
    return out, summary


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    records = _load_all()
    if not records:
        print("no D15b_*.json found", file=sys.stderr)
        return 1
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    p1 = build_condition_table(records)
    p2 = build_order_table(records)
    p3 = build_convergence_table(records)
    p4, summary = build_coverage_table(records)
    (REPORT_DIR / "D15c_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False))
    for x in (p1, p2, p3, p4):
        print("wrote", x.relative_to(PROJECT_ROOT))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
