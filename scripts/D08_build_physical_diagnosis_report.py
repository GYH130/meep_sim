"""Build the final physical diagnosis report from D00-D07 outputs.

This is a pure post-processing script.  It does not import Meep and does not
run any simulations.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path("/private/tmp") / "matplotlib-codex-cache"),
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.io_utils import ensure_dir, project_path, save_figure, setup_logger


STATUS_ORDER = {"PASS": 3, "WARNING": 2, "FAIL": 1, "NOT_RUN": 0}
STATUS_COLORS = {
    "PASS": "#54A24B",
    "WARNING": "#F58518",
    "FAIL": "#E45756",
    "NOT_RUN": "#9D9D9D",
}


@dataclass
class StageResult:
    key: str
    title: str
    status: str
    summary: str
    metrics: dict


def _paths() -> dict[str, Path]:
    return {
        "d00_lossless": project_path(
            "results", "diagnostics", "tables",
            "D00_lossless_interface_validation.csv",
        ),
        "d00_flat": project_path(
            "results", "diagnostics", "tables",
            "D00_flat_ti_fresnel_validation.csv",
        ),
        "d01": project_path(
            "results", "diagnostics", "tables", "D01_solver_regression_metrics.csv",
        ),
        "d02": project_path(
            "results", "diagnostics", "tables", "D02_polarization_metrics.csv",
        ),
        "d03": project_path(
            "results", "diagnostics", "tables", "D03_resolution_convergence.csv",
        ),
        "d04": project_path(
            "results", "diagnostics", "tables", "D04_absorbed_power_integrals.csv",
        ),
        "d04_hotspot": project_path(
            "results", "diagnostics", "tables", "D04_hotspot_metrics.csv",
        ),
        "d05": project_path(
            "results", "diagnostics", "tables", "D05_directionality_metrics.csv",
        ),
        "d06": project_path(
            "results", "diagnostics", "tables", "D06_energy_channel_summary.csv",
        ),
        "d07": project_path(
            "results", "diagnostics", "tables", "D07_oxide_enhancement_metrics.csv",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D08_final_physical_diagnosis_report.md",
        ),
        "dashboard": project_path(
            "results", "diagnostics", "figures", "D08_diagnosis_dashboard.png",
        ),
    }


def _read_csv(path: Path, logger) -> pd.DataFrame | None:
    if not path.exists():
        logger.warning("Input not found: %s", path)
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:
        logger.warning("Could not read %s: %s", path, exc)
        return None


def _bool_series_all_true(series: pd.Series) -> bool:
    if series.empty:
        return False
    return series.astype(str).str.upper().isin(["TRUE", "PASS"]).all()


def _status_from_pass_column(df: pd.DataFrame | None) -> str:
    if df is None:
        return "NOT_RUN"
    if "pass_or_fail" not in df.columns:
        return "WARNING"
    vals = df["pass_or_fail"].astype(str).str.upper()
    if vals.eq("FAIL").any():
        return "FAIL"
    if vals.eq("WARNING").any():
        return "WARNING"
    if vals.eq("PASS").all():
        return "PASS"
    return "WARNING"


def _diagnose_d00(paths: dict[str, Path], logger) -> StageResult:
    lossless = _read_csv(paths["d00_lossless"], logger)
    flat = _read_csv(paths["d00_flat"], logger)
    if lossless is None or flat is None:
        return StageResult("D00", "Flux Sign And Fresnel Validation", "NOT_RUN",
                           "D00 输入文件不完整。", {})

    max_lossless_r = float(lossless.get("abs_R_error", pd.Series([np.nan])).max())
    max_lossless_t = float(lossless.get("abs_T_error", pd.Series([np.nan])).max())
    max_flat_a = float(flat.get("abs_A_error", pd.Series([np.nan])).max())
    selected_defs = set(lossless.get("selected_T_definition", pd.Series(dtype=str)).astype(str))
    signed_selected = selected_defs == {"signed"} or "signed" in selected_defs
    status = "PASS"
    if not signed_selected or max_lossless_r > 0.02 or max_lossless_t > 0.02:
        status = "FAIL"
    elif max_flat_a > 0.05:
        status = "WARNING"
    summary = (
        f"透射符号选择为 {sorted(selected_defs)}；无损界面最大 R/T 误差 "
        f"{max_lossless_r:.3g}/{max_lossless_t:.3g}；平面 Ti 最大 A 误差 "
        f"{max_flat_a:.3g}。"
    )
    return StageResult("D00", "Flux Sign And Fresnel Validation", status, summary, {
        "lossless_R_error": max_lossless_r,
        "lossless_T_error": max_lossless_t,
        "flat_Ti_A_error": max_flat_a,
        "signed_T_selected": signed_selected,
    })


def _diagnose_d01(paths: dict[str, Path], logger) -> StageResult:
    df = _read_csv(paths["d01"], logger)
    status = _status_from_pass_column(df)
    if df is None:
        return StageResult("D01", "Solver Regression", status, "D01 未运行。", {})
    worst_delta = float(df.get("mean_A_delta_after_flux_fix", pd.Series([0])).abs().max())
    summary = f"回归表状态 {status}；通量修正前后平均 A 最大变化 {worst_delta:.3g}。"
    return StageResult("D01", "Solver Regression", status, summary, {
        "max_mean_A_delta_after_flux_fix": worst_delta,
    })


def _diagnose_d02(paths: dict[str, Path], logger, hz_tol: float) -> StageResult:
    df = _read_csv(paths["d02"], logger)
    status = _status_from_pass_column(df)
    if df is None:
        return StageResult("D02", "Polarization Diagnostic", status, "D02 未运行。", {})
    max_delta = float(df.get("delta_mean_A_Hz_minus_Ez", pd.Series([0])).max())
    min_delta = float(df.get("delta_mean_A_Hz_minus_Ez", pd.Series([0])).min())
    hz_flag = bool(
        df.get("hz_mean_enhancement_flag", pd.Series([False])).astype(str)
        .str.upper().eq("TRUE").any()
    ) or max_delta >= hz_tol
    if status == "PASS" and hz_flag:
        status = "WARNING"
    summary = (
        f"Hz-Ez 的 8-13 um 平均 A 差值范围 [{min_delta:.3g}, {max_delta:.3g}]；"
        f"{'Hz 显著高于 Ez，固定 Ez 结论需谨慎。' if hz_flag else '未发现 Hz 明显高于 Ez。'}"
    )
    return StageResult("D02", "Polarization Diagnostic", status, summary, {
        "max_delta_mean_A_Hz_minus_Ez": max_delta,
        "hz_significantly_higher": hz_flag,
    })


def _diagnose_d03(paths: dict[str, Path], logger) -> StageResult:
    df = _read_csv(paths["d03"], logger)
    status = _status_from_pass_column(df)
    if df is None:
        return StageResult("D03", "Numerical Convergence", status, "D03 未运行。", {})
    max_diff = float(df.get(
        "max_abs_difference_from_highest_accuracy_case",
        pd.Series([np.nan]),
    ).max())
    failed = []
    if "pass_or_fail" in df.columns:
        failed = sorted(df.loc[
            df["pass_or_fail"].astype(str).str.upper() == "FAIL",
            "case_name",
        ].astype(str).unique().tolist())
    summary = (
        f"最大谱差 {max_diff:.3g}；失败对象 {failed if failed else '无'}。"
    )
    return StageResult("D03", "Numerical Convergence", status, summary, {
        "max_abs_difference_from_highest_accuracy_case": max_diff,
        "failed_cases": failed,
    })


def _diagnose_d04(paths: dict[str, Path], logger, local_e2_tol: float) -> StageResult:
    df = _read_csv(paths["d04"], logger)
    status = _status_from_pass_column(df)
    if df is None:
        return StageResult("D04", "Field And Absorbed Power", status, "D04 未运行。", {})
    max_volume_diff = float(df.get("volume_flux_abs_difference", pd.Series([np.nan])).max())
    hotspot = _read_csv(paths["d04_hotspot"], logger)
    max_metal_e2 = np.nan
    max_wall_fraction = np.nan
    weak_localization = False
    if hotspot is not None:
        metal_col = (
            "max_metal_e2_enhancement"
            if "max_metal_e2_enhancement" in hotspot.columns
            else "max_E2_metal_over_air_median"
        )
        wall_col = (
            "mean_wall_fraction"
            if "mean_wall_fraction" in hotspot.columns
            else "wall_power_fraction"
        )
        max_metal_e2 = float(hotspot.get(metal_col, pd.Series([np.nan])).max())
        max_wall_fraction = float(hotspot.get(wall_col, pd.Series([np.nan])).max())
        weak_localization = np.isfinite(max_metal_e2) and max_metal_e2 < local_e2_tol
    if status == "PASS" and weak_localization:
        status = "WARNING"
    summary = (
        f"Flux/volume 最大差 {max_volume_diff:.3g}；最大 metal |E|^2 增强 "
        f"{max_metal_e2:.3g}，最大槽壁耗散占比 {max_wall_fraction:.3g}。"
    )
    return StageResult("D04", "Field And Absorbed Power", status, summary, {
        "max_volume_flux_abs_difference": max_volume_diff,
        "max_metal_e2_enhancement": max_metal_e2,
        "max_wall_absorbed_fraction": max_wall_fraction,
        "weak_localized_loss_channel": weak_localization,
    })


def _diagnose_d05(paths: dict[str, Path], logger,
                  low_a_threshold: float, direction_ratio_tol: float) -> StageResult:
    df = _read_csv(paths["d05"], logger)
    if df is None:
        return StageResult("D05", "Angle-Resolved Absorptance", "NOT_RUN",
                           "D05 未运行。", {})
    band = df[df.get("metric_scope", pd.Series(dtype=str)).astype(str) == "band_average"]
    if band.empty:
        band = df.copy()
    mean_a = float(band.get("mean_A_theta", band.get("theta0_absorptance", pd.Series([np.nan]))).max())
    max_ratio = float(band.get("mean_directionality_ratio", band.get("ratio_at_theta0", pd.Series([1]))).max())
    max_asym = float(band.get("asymmetry_at_theta0", pd.Series([0])).abs().max())
    directional = max_ratio >= direction_ratio_tol or max_asym >= 0.05
    low_a = np.isfinite(mean_a) and mean_a < low_a_threshold
    status = "WARNING" if directional else "PASS"
    summary = (
        f"最大角平均 A 约 {mean_a:.3g}；最大方向比 {max_ratio:.3g}；"
        f"最大 A(+θ)-A(-θ) {max_asym:.3g}。"
    )
    return StageResult("D05", "Angle-Resolved Absorptance", status, summary, {
        "max_mean_A_theta": mean_a,
        "max_directionality_ratio": max_ratio,
        "max_asymmetry_at_theta0": max_asym,
        "low_absorption": low_a,
        "directional": directional,
    })


def _diagnose_d06(paths: dict[str, Path], logger) -> StageResult:
    df = _read_csv(paths["d06"], logger)
    status = _status_from_pass_column(df)
    if df is None:
        return StageResult("D06", "Diffraction Energy Channels", status, "D06 未运行。", {})
    max_resid = float(df.get("order_sum_residual", pd.Series([0])).abs().max())
    mean_reflection = float(df.get("R_flux_monitor", pd.Series([np.nan])).mean())
    mean_absorption = float(df.get("absorptance", pd.Series([np.nan])).mean())
    nonspec = float(df.get("nonspecular_reflection", pd.Series([np.nan])).mean())
    summary = (
        f"阶功率求和残差最大 {max_resid:.3g}；平均 R={mean_reflection:.3g}，"
        f"A={mean_absorption:.3g}，非镜面反射={nonspec:.3g}。"
    )
    return StageResult("D06", "Diffraction Energy Channels", status, summary, {
        "max_order_sum_residual": max_resid,
        "mean_reflection": mean_reflection,
        "mean_absorption": mean_absorption,
        "mean_nonspecular_reflection": nonspec,
    })


def _diagnose_d07(paths: dict[str, Path], logger, oxide_tol: float) -> StageResult:
    df = _read_csv(paths["d07"], logger)
    status = _status_from_pass_column(df)
    if df is None:
        return StageResult("D07", "Oxide And Multiscale Sensitivity", status,
                           "D07 未运行。", {})
    max_enh = float(df.get("enhancement_over_bare", pd.Series([0])).max())
    sig = bool(df.get("significant_enhancement", pd.Series([False])).astype(str)
               .str.upper().eq("TRUE").any()) or max_enh >= oxide_tol
    placeholder = True
    if status == "PASS" and sig:
        status = "WARNING"
    summary = (
        f"最大 placeholder 氧化层增强 {max_enh:.3g}；"
        f"{'超过显著阈值。' if sig else '未超过显著阈值。'} "
        "TiO2 仍为 placeholder，不能作实验定量预测。"
    )
    return StageResult("D07", "Oxide And Multiscale Sensitivity", status, summary, {
        "max_oxide_enhancement_over_bare": max_enh,
        "oxide_significant": sig,
        "uses_placeholder_tio2": placeholder,
    })


def _make_dashboard(stages: list[StageResult], paths: dict[str, Path]) -> None:
    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], width_ratios=[1.2, 1])
    ax_status = fig.add_subplot(gs[0, :])
    ax_abs = fig.add_subplot(gs[1, 0])
    ax_notes = fig.add_subplot(gs[1, 1])

    names = [s.key for s in stages]
    values = [STATUS_ORDER[s.status] for s in stages]
    colors = [STATUS_COLORS[s.status] for s in stages]
    ax_status.bar(names, values, color=colors)
    ax_status.set_ylim(0, 3.4)
    ax_status.set_yticks([0, 1, 2, 3])
    ax_status.set_yticklabels(["NOT_RUN", "FAIL", "WARNING", "PASS"])
    ax_status.set_title("D00-D07 Diagnostic Status")
    ax_status.grid(axis="y", alpha=0.25)
    for i, s in enumerate(stages):
        ax_status.text(i, values[i] + 0.08, s.status, ha="center", fontsize=8)

    metric_labels = []
    metric_values = []
    lookup = {s.key: s for s in stages}
    for key, metric, label in [
        ("D02", "max_delta_mean_A_Hz_minus_Ez", "Hz-Ez mean A"),
        ("D03", "max_abs_difference_from_highest_accuracy_case", "D03 max delta A"),
        ("D04", "max_metal_e2_enhancement", "Max metal E2 enh."),
        ("D05", "max_directionality_ratio", "D05 dir. ratio"),
        ("D06", "mean_nonspecular_reflection", "D06 nonspec R"),
        ("D07", "max_oxide_enhancement_over_bare", "D07 oxide enh."),
    ]:
        val = lookup.get(key, StageResult(key, "", "NOT_RUN", "", {})).metrics.get(metric, np.nan)
        if np.isfinite(val):
            metric_labels.append(label)
            metric_values.append(float(val))
    if metric_values:
        ax_abs.barh(metric_labels, metric_values, color="#4C78A8")
        ax_abs.set_title("Key Diagnostic Metrics")
        ax_abs.grid(axis="x", alpha=0.25)
    else:
        ax_abs.text(0.5, 0.5, "No metrics available", ha="center", va="center")
        ax_abs.set_axis_off()

    ax_notes.set_axis_off()
    notes = [
        "Decision Logic",
        "D00 fail: stop physical interpretation",
        "D03 fail: defer quantitative structure claims",
        "D02 warning: fixed Ez may miss polarization response",
        "D04 warning: localized loss channel may be weak",
        "D05/D06: directionality vs total absorption",
        "D07: oxide trends only if material model is valid",
    ]
    ax_notes.text(0, 1, "\n".join(notes), va="top", fontsize=10)
    fig.tight_layout()
    save_figure(fig, paths["dashboard"])
    plt.close(fig)


def _stage_table(stages: list[StageResult]) -> str:
    lines = ["| Stage | Status | Summary |", "|---|---|---|"]
    for s in stages:
        lines.append(f"| {s.key} {s.title} | **{s.status}** | {s.summary} |")
    return "\n".join(lines)


def _write_report(stages: list[StageResult], paths: dict[str, Path],
                  final_logic: dict, args: argparse.Namespace) -> None:
    ensure_dir(paths["report"].parent)
    confirmed = final_logic["confirmed"]
    hypotheses = final_logic["hypotheses"]
    ml_ready = final_logic["ml_ready"]
    top_compute = final_logic["top_compute"]
    top_experiment = final_logic["top_experiment"]

    report = f"""# D08 Final Physical Diagnosis Report

## Executive Diagnosis

{final_logic['headline']}

## Stage Status

{_stage_table(stages)}

## Confirmed Conclusions

{chr(10).join(f'- {item}' for item in confirmed)}

## Unconfirmed Hypotheses

{chr(10).join(f'- {item}' for item in hypotheses)}

## Next Computation To Spend Resources On

{chr(10).join(f'- {item}' for item in top_compute)}

## Next Experimental Data To Collect

{chr(10).join(f'- {item}' for item in top_experiment)}

## Machine-Learning Dataset Readiness

**{'值得进入' if ml_ready else '不建议进入'}机器学习数据集生成阶段。**

{final_logic['ml_reason']}

## Important Caveats

- D08 不重新运行 Meep，只汇总已存在 CSV；若上游 D00-D07 是 smoke test，本报告也只是 smoke 级别。
- 如果 D00 或 D03 未通过，不应把当前结构吸收率作为机器学习标签。
- D07 当前 TiO2 是 `placeholder_demo`，不能作为实验定量预测。

## Outputs

- Dashboard: `{paths['dashboard']}`
- Report: `{paths['report']}`

## Thresholds Used

```json
{json.dumps(vars(args), indent=2, ensure_ascii=False)}
```
"""
    paths["report"].write_text(report, encoding="utf-8")


def _final_logic(stages: list[StageResult]) -> dict:
    lookup = {s.key: s for s in stages}
    confirmed: list[str] = []
    hypotheses: list[str] = []
    top_compute: list[str] = []
    top_experiment: list[str] = []

    d00 = lookup["D00"]
    d03 = lookup["D03"]
    d02 = lookup["D02"]
    d04 = lookup["D04"]
    d05 = lookup["D05"]
    d06 = lookup["D06"]
    d07 = lookup["D07"]

    if d00.status == "FAIL":
        headline = "当前吸收率计算尚不可信，禁止解释结构物理。"
        ml_ready = False
        ml_reason = "D00 未通过，R/T/A 标签的符号或归一化仍不可靠。"
    elif d03.status == "FAIL":
        headline = "D00 已通过，但 D03 显示数值未完全收敛，结构定量结论暂缓。"
        ml_ready = False
        ml_reason = "D03 存在失败项，当前不应将这些结果作为机器学习训练标签。"
    else:
        headline = "基础求解器可信，当前可做物理趋势判断，但仍需注意模型保真度。"
        ml_ready = True
        ml_reason = "D00/D03 未阻塞；仍建议只把已收敛参数下的结果纳入数据集。"

    if d00.status == "PASS":
        confirmed.append("D00 确认了透射通量应使用 signed transmittance，平面 Ti 与 Fresnel 基准基本一致。")
    if d02.metrics.get("hz_significantly_higher", False):
        confirmed.append("D02 显示 Hz 相对 Ez 有明显增强，过去只看 Ez 可能遗漏偏振响应。")
    else:
        confirmed.append("D02 未显示 Hz 相对 Ez 的强烈、稳健跃迁，固定 Ez 结论不是唯一问题来源。")
    if d04.metrics.get("weak_localized_loss_channel", False):
        confirmed.append("D04 支持当前裸 Ti 几何没有形成强局域场增强。")
    if d05.metrics.get("low_absorption", False) and d05.metrics.get("directional", False):
        confirmed.append("D05 支持结构更像方向调制单元而非高发射背景层。")
    elif d05.metrics:
        confirmed.append("D05 当前结果未显示强方向性；角分辨结论仍受 smoke/采样范围限制。")
    if d06.status in {"PASS", "WARNING"} and d06.metrics:
        confirmed.append("D06 的衍射阶功率和总反射 flux 闭合，说明反射能量通道分析可用。")
    if d07.metrics.get("oxide_significant", False):
        confirmed.append("D07 placeholder 扫描显示氧化层可能显著增强吸收，但不能定量外推到真实 TiO2。")
    elif d07.metrics:
        confirmed.append("D07 placeholder smoke 未显示超过阈值的氧化层增强。")

    hypotheses.extend([
        "裸 Ti 斜槽提升有限可能来自缺少有效局域损耗模式，也可能来自 2D 理想几何过于简单。",
        "真实飞秒加工样品的高发射若存在，可能依赖氧化层、粗糙度、再凝固层或多尺度结构。",
        "斜槽几何更可能负责方向性/衍射通道调制，而不是单独提供宽带高吸收。",
    ])

    if d03.status == "FAIL":
        top_compute.append("优先重跑 D03 中失败的 slanted_Hz resolution 收敛，用 resolution>=48/64 确定统一基准。")
    top_compute.extend([
        "用 D03 推荐的收敛设置重跑 D02/D05/D06 的关键点，避免 smoke 参数主导结论。",
        "替换 D07 的 TiO2 placeholder，加入可靠 mid-IR TiO2/TiOx n,k 或 Drude-Lorentz 拟合。",
        "在氧化层模型通过后，再做几何参数扫描和数据集生成。",
    ])
    top_experiment.extend([
        "氧化层厚度和空间分布，尤其槽壁/槽底覆盖。",
        "TiO2/TiOx 成分、相态和 8-13 um 光学常数。",
        "截面几何：槽深、开口、侧壁角、圆角、周期分布。",
        "多尺度粗糙度和再凝固颗粒，用于决定是否需要 3D/有效介质模型。",
    ])

    return dict(
        headline=headline,
        confirmed=confirmed,
        hypotheses=hypotheses,
        top_compute=top_compute,
        top_experiment=top_experiment,
        ml_ready=ml_ready,
        ml_reason=ml_reason,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build final D00-D07 physical diagnosis report.")
    p.add_argument("--hz_delta_tol", type=float, default=0.02)
    p.add_argument("--local_e2_tol", type=float, default=1.0)
    p.add_argument("--low_absorption_threshold", type=float, default=0.2)
    p.add_argument("--directionality_ratio_tol", type=float, default=1.2)
    p.add_argument("--oxide_enhancement_tol", type=float, default=0.03)
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    logger = setup_logger("D08_build_physical_diagnosis_report")
    paths = _paths()
    ensure_dir(paths["report"].parent)
    ensure_dir(paths["dashboard"].parent)

    stages = [
        _diagnose_d00(paths, logger),
        _diagnose_d01(paths, logger),
        _diagnose_d02(paths, logger, args.hz_delta_tol),
        _diagnose_d03(paths, logger),
        _diagnose_d04(paths, logger, args.local_e2_tol),
        _diagnose_d05(paths, logger, args.low_absorption_threshold, args.directionality_ratio_tol),
        _diagnose_d06(paths, logger),
        _diagnose_d07(paths, logger, args.oxide_enhancement_tol),
    ]
    logic = _final_logic(stages)
    _make_dashboard(stages, paths)
    _write_report(stages, paths, logic, args)
    logger.info("Wrote %s", paths["report"])
    logger.info("Wrote %s", paths["dashboard"])


if __name__ == "__main__":
    main()
