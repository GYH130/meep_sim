"""D00_flux_sign_and_fresnel_validation.py — 通量符号与 Fresnel 基准诊断

研究目的
--------
在继续解释 Ti 微结构发射率之前，独立验证当前项目中 R、T、A 的通量符号、
归一化方式和半无限平面界面 Fresnel 基准是否一致。

诊断任务
--------
A. 无损空气 / n=1.5 半无限界面：
   - 保持 `run_periodic_2d_metal_spectrum()` 的源、反射 monitor、透射 monitor
     布局；
   - 让介质从界面向下延伸穿过 bottom PML，避免二次介质/空气界面；
   - 同时计算 T_raw = trans_flux / |input_flux| 与
     T_signed = -trans_flux / |input_flux|，用 Fresnel 解析解自动判定符号。

B. 平面 Ti 半无限界面近似：
   - 使用 `get_ti_medium()` 和同一个公共仿真函数；
   - 用 Ti.epsilon(f) 计算半无限空气/Ti Fresnel 反射率；
   - 仅在 Rakić Ti 标定有效范围内定量判定，超出波长标记 warning。

物理假设
--------
1. 所有长度单位为 μm，Meep 频率 f = 1 / wavelength_um；
2. 正入射、2D、Ez 偏振；对各向同性平面界面，正入射下 TE/TM 等价；
3. 无损介质延伸穿过下方 PML，界面没有吸收，故 R + T = 1；
4. Ti 基底在中红外近似不透明，半无限 Fresnel 基准取 A = 1 - R；
5. Ti Rakić 模型超过 12.398 μm 的部分只作为外推，不作严格通过/失败判定。

通过/失败判据
-------------
1. 无损界面：选定 T 定义后 max(|R-R_th|, |T-T_th|, |A|) <= lossless_tol；
2. 平面 Ti：有效波段内 max(|R_meep-R_Fresnel|, |A_meep-A_Fresnel|)
   <= ti_tol，且 max(|T_selected|) <= ti_trans_tol；
3. 公共函数返回的 T 字段应与无损界面选出的 T 定义一致。
"""

from __future__ import annotations

import argparse
import cmath
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meep as mp
import numpy as np
import pandas as pd

from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium
from src.simulation import run_periodic_2d_metal_spectrum


DEFAULTS = dict(
    wavelength_min_um=8.0,
    wavelength_max_um=13.0,
    nfreq=21,
    resolution=32,
    period_um=1.0,
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    decay_db=40.0,
    lossless_n=1.5,
    lossless_tol=3e-2,
    ti_tol=5e-2,
    ti_trans_tol=2e-2,
)


def _diagnostic_paths() -> dict[str, Path]:
    return {
        "lossless_csv": project_path(
            "results", "diagnostics", "tables",
            "D00_lossless_interface_validation.csv",
        ),
        "ti_csv": project_path(
            "results", "diagnostics", "tables",
            "D00_flat_ti_fresnel_validation.csv",
        ),
        "lossless_png": project_path(
            "results", "diagnostics", "figures",
            "D00_lossless_RT_validation.png",
        ),
        "ti_png": project_path(
            "results", "diagnostics", "figures",
            "D00_flat_ti_fresnel_comparison.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D00_flux_sign_validation_report.md",
        ),
    }


def _write_csv(df: pd.DataFrame, path: Path) -> Path:
    ensure_dir(path.parent)
    df.to_csv(path, index=False)
    return path


def _lossless_theory(n: float) -> tuple[float, float]:
    r = ((1.0 - n) / (1.0 + n)) ** 2
    t = 4.0 * n / (1.0 + n) ** 2
    return float(r), float(t)


def _ti_fresnel(wavelengths_um: np.ndarray, ti_medium) -> tuple[np.ndarray, np.ndarray]:
    r_values = []
    lo_valid, hi_valid = TI_RAKIC_VALID_LAMBDA_UM
    for wl in wavelengths_um:
        if wl < lo_valid or wl > hi_valid:
            r_values.append(np.nan)
            continue
        eps = ti_medium.epsilon(1.0 / float(wl))[0][0]
        n_ti = cmath.sqrt(eps)
        r_values.append(abs((1.0 - n_ti) / (1.0 + n_ti)) ** 2)
    r = np.asarray(r_values, dtype=float)
    return r, 1.0 - r


def _select_transmission_definition(lossless: pd.DataFrame) -> tuple[str, float, float]:
    raw_err = float(np.max(np.abs(lossless["T_raw"] - lossless["theory_T_or_A"])))
    signed_err = float(np.max(np.abs(lossless["T_signed"] - lossless["theory_T_or_A"])))
    selected = "raw" if raw_err <= signed_err else "signed"
    return selected, raw_err, signed_err


def _with_selected_transmission(df: pd.DataFrame, selected: str) -> pd.DataFrame:
    df = df.copy()
    if selected == "raw":
        df["T_selected"] = df["T_raw"]
    elif selected == "signed":
        df["T_selected"] = df["T_signed"]
    else:
        raise ValueError(f"unknown selected T definition: {selected}")
    df["A_selected"] = 1.0 - df["R"] - df["T_selected"]
    return df


def _base_rows(result: dict) -> pd.DataFrame:
    wl = result["wavelength_um"]
    input_flux = result.get("input_flux_raw")
    refl_flux = result.get("reflection_flux_raw")
    trans_flux = result.get("transmission_flux_raw", result.get("raw_trans_flux"))
    if input_flux is None or refl_flux is None or trans_flux is None:
        raise RuntimeError(
            "run_periodic_2d_metal_spectrum() 未返回 raw flux 字段；"
            "请先使用包含 D00 诊断字段的 src/simulation.py。"
        )
    df = pd.DataFrame({
        "wavelength_um": wl,
        "input_flux_raw": input_flux,
        "reflection_flux_raw": refl_flux,
        "transmission_flux_raw": trans_flux,
        "R": refl_flux / np.abs(input_flux),
        "T_raw": trans_flux / np.abs(input_flux),
        "T_signed": -trans_flux / np.abs(input_flux),
        "T_public": result["T"],
        "A_public": result["A"],
    })
    return df


def run_lossless_interface(args: argparse.Namespace, logger) -> tuple[pd.DataFrame, dict]:
    n = args.lossless_n
    medium = mp.Medium(index=n)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        # run_periodic_2d_metal_spectrum() places the transmission monitor between
        # y_surface-substrate_thickness and bottom PML.  Extending the dielectric
        # through the bottom PML avoids a spurious dielectric/air interface at the
        # PML inner boundary and keeps the transmission monitor inside the medium.
        thickness = substrate_thickness_um + 2.0 * args.pml_thickness_um
        return [
            mp.Block(
                material=medium,
                center=mp.Vector3(0, y_surface_um - 0.5 * thickness, 0),
                size=mp.Vector3(args.period_um, thickness, mp.inf),
            )
        ]

    logger.info("running lossless interface validation: n=%.4g", n)
    result = run_periodic_2d_metal_spectrum(
        geometry_factory=factory,
        period_um=args.period_um,
        wavelength_min_um=args.wavelength_min_um,
        wavelength_max_um=args.wavelength_max_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        nfreq=args.nfreq,
        decay_db=args.decay_db,
        source_component="Ez",
        logger=logger,
    )
    theory_r, theory_t = _lossless_theory(n)
    df = _base_rows(result)
    df["theory_R"] = theory_r
    df["theory_T_or_A"] = theory_t
    selected, raw_err, signed_err = _select_transmission_definition(df)
    df = _with_selected_transmission(df, selected)
    df["residual"] = df["R"] + df["T_selected"] - 1.0
    df["abs_R_error"] = np.abs(df["R"] - df["theory_R"])
    df["abs_T_error"] = np.abs(df["T_selected"] - df["theory_T_or_A"])
    df["selected_T_definition"] = selected
    df["normalization_note"] = (
        "R = reflection_flux_raw / abs(input_flux_raw); "
        "T_raw = transmission_flux_raw / abs(input_flux_raw); "
        "T_signed = -transmission_flux_raw / abs(input_flux_raw)"
    )
    summary = {
        "selected_T_definition": selected,
        "raw_T_max_abs_error": raw_err,
        "signed_T_max_abs_error": signed_err,
        "max_abs_R_error": float(df["abs_R_error"].max()),
        "max_abs_T_error": float(df["abs_T_error"].max()),
        "max_abs_A_selected": float(np.abs(df["A_selected"]).max()),
        "max_abs_residual": float(np.abs(df["residual"]).max()),
        "public_T_matches_selected": bool(
            np.allclose(df["T_public"], df["T_selected"], atol=1e-10, rtol=1e-10)
        ),
        "walltime_s": float(result["walltime_s"]),
        "geometry_y": result.get("geometry_y", {}),
    }
    summary["passed"] = bool(
        max(
            summary["max_abs_R_error"],
            summary["max_abs_T_error"],
            summary["max_abs_A_selected"],
        ) <= args.lossless_tol
    )
    logger.info("lossless selected T definition: %s", selected)
    logger.info("lossless summary: %s", json.dumps(summary, ensure_ascii=False, default=str))
    return df, summary


def run_flat_ti_validation(
    args: argparse.Namespace,
    selected_t: str,
    logger,
) -> tuple[pd.DataFrame, dict]:
    ti = get_ti_medium(
        lambda_min_um=args.wavelength_min_um,
        lambda_max_um=args.wavelength_max_um,
    )

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        thickness = substrate_thickness_um + 2.0 * args.pml_thickness_um
        return [
            mp.Block(
                material=ti,
                center=mp.Vector3(0, y_surface_um - 0.5 * thickness, 0),
                size=mp.Vector3(args.period_um, thickness, mp.inf),
            )
        ]

    logger.info("running flat Ti Fresnel validation")
    result = run_periodic_2d_metal_spectrum(
        geometry_factory=factory,
        period_um=args.period_um,
        wavelength_min_um=args.wavelength_min_um,
        wavelength_max_um=args.wavelength_max_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        nfreq=args.nfreq,
        decay_db=args.decay_db,
        source_component="Ez",
        logger=logger,
    )
    df = _base_rows(result)
    df = _with_selected_transmission(df, selected_t)
    theory_r, theory_a = _ti_fresnel(df["wavelength_um"].to_numpy(), ti)
    df["theory_R"] = theory_r
    df["theory_T_or_A"] = theory_a
    df["residual"] = df["A_selected"] - df["theory_T_or_A"]
    df["abs_R_error"] = np.abs(df["R"] - df["theory_R"])
    df["abs_A_error"] = np.abs(df["A_selected"] - df["theory_T_or_A"])
    df["selected_T_definition"] = selected_t
    lo_valid, hi_valid = TI_RAKIC_VALID_LAMBDA_UM
    valid = (df["wavelength_um"] >= lo_valid) & (df["wavelength_um"] <= hi_valid)
    df["ti_model_valid"] = valid
    df["validity_note"] = np.where(
        valid,
        "within Ti Rakić validity range",
        "warning: outside Ti Rakić validity range; excluded from pass/fail",
    )
    df["normalization_note"] = (
        "R = reflection_flux_raw / abs(input_flux_raw); "
        "T_selected follows lossless-interface sign validation; "
        "A_selected = 1 - R - T_selected"
    )

    valid_df = df[df["ti_model_valid"]]
    if valid_df.empty:
        raise RuntimeError("Ti 有效波段内没有采样点，无法定量验证 Fresnel 基准。")
    summary = {
        "max_abs_R_error_valid": float(valid_df["abs_R_error"].max()),
        "max_abs_A_error_valid": float(valid_df["abs_A_error"].max()),
        "max_abs_T_selected_valid": float(np.abs(valid_df["T_selected"]).max()),
        "mean_A_selected_8_13um": float(df["A_selected"].mean()),
        "mean_A_fresnel_8_13um": float(df["theory_T_or_A"].mean()),
        "warning_wavelength_count": int((~df["ti_model_valid"]).sum()),
        "public_T_matches_selected": bool(
            np.allclose(df["T_public"], df["T_selected"], atol=1e-10, rtol=1e-10)
        ),
        "walltime_s": float(result["walltime_s"]),
        "geometry_y": result.get("geometry_y", {}),
    }
    summary["passed"] = bool(
        max(summary["max_abs_R_error_valid"], summary["max_abs_A_error_valid"])
        <= args.ti_tol
        and summary["max_abs_T_selected_valid"] <= args.ti_trans_tol
    )
    logger.info("flat Ti summary: %s", json.dumps(summary, ensure_ascii=False, default=str))
    return df, summary


def plot_lossless(df: pd.DataFrame, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(df["wavelength_um"], df["R"], "o-", label="Meep R")
    ax.plot(df["wavelength_um"], df["T_raw"], "s--", label="T raw")
    ax.plot(df["wavelength_um"], df["T_signed"], "^--", label="T signed")
    ax.plot(df["wavelength_um"], df["T_selected"], "k-", lw=2.0, label="T selected")
    ax.axhline(float(df["theory_R"].iloc[0]), color="C0", ls=":", label="Fresnel R")
    ax.axhline(float(df["theory_T_or_A"].iloc[0]), color="C3", ls=":", label="Fresnel T")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Power fraction")
    ax.set_title("D00 lossless air / n interface flux-sign validation")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(ncol=2, fontsize=8)
    return save_figure(fig, out_path)


def plot_ti(df: pd.DataFrame, out_path: Path) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)
    ax0, ax1 = axes
    ax0.plot(df["wavelength_um"], df["R"], "o-", label="Meep R")
    ax0.plot(df["wavelength_um"], df["theory_R"], "k--", label="Fresnel R")
    ax0.plot(df["wavelength_um"], df["A_selected"], "s-", label="Meep A selected")
    ax0.plot(df["wavelength_um"], df["theory_T_or_A"], "k:", label="Fresnel A")
    ax0.axvspan(TI_RAKIC_VALID_LAMBDA_UM[1], df["wavelength_um"].max(),
                color="orange", alpha=0.10, label="Ti extrapolation")
    ax0.set_ylabel("Power fraction")
    ax0.set_title("D00 flat Ti Meep vs Fresnel")
    ax0.grid(True, ls=":", alpha=0.6)
    ax0.legend(ncol=2, fontsize=8)

    ax1.plot(df["wavelength_um"], df["R"] - df["theory_R"], "o-", label="R error")
    ax1.plot(df["wavelength_um"], df["A_selected"] - df["theory_T_or_A"],
             "s-", label="A error")
    ax1.plot(df["wavelength_um"], df["T_selected"], "^-", label="T selected")
    ax1.axhline(0, color="k", lw=0.8)
    ax1.set_xlabel("Wavelength (um)")
    ax1.set_ylabel("Difference")
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend(fontsize=8)
    return save_figure(fig, out_path)


def _status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def write_report(
    args: argparse.Namespace,
    lossless_summary: dict,
    ti_summary: dict,
    paths: dict[str, Path],
) -> Path:
    selected = lossless_summary["selected_T_definition"]
    public_matches = (
        lossless_summary["public_T_matches_selected"]
        and ti_summary["public_T_matches_selected"]
    )
    src_must_change = not public_matches
    if selected == "signed":
        old_recompute = (
            "本次验证选择 signed T。凡是在本次修正前按 raw T 写出的旧结果都应重新计算，"
            "包括 flat_ti_spectrum.csv、periodic_groove_spectrum.csv、"
            "slanted_groove_spectra.csv 以及依赖这些表的图和报告。"
        )
    elif src_must_change:
        old_recompute = (
            "公共函数当前 T 字段与选定定义不一致；依赖 run_periodic_2d_metal_spectrum() "
            "输出 T/A 的旧结果均需重新计算。"
        )
    else:
        old_recompute = (
            "本次 D00 运行时公共函数的 T 字段已与选定符号一致；若旧表缺少 raw flux 与 "
            "signed_transmittance 元数据，建议重新运行以便追溯。"
        )
    lines = [
        "# D00 Flux Sign And Fresnel Validation Report",
        "",
        "## Run Configuration",
        "",
        f"- wavelength range: {args.wavelength_min_um:g}-{args.wavelength_max_um:g} um",
        f"- nfreq: {args.nfreq}",
        f"- resolution: {args.resolution} pixels/um",
        f"- pml_thickness_um: {args.pml_thickness_um:g}",
        f"- substrate_thickness_um: {args.substrate_thickness_um:g}",
        f"- air_buffer_um: {args.air_buffer_um:g}",
        f"- source/monitor convention: source propagates toward -y; raw flux columns are preserved",
        "",
        "## Physical Assumptions",
        "",
        "- Length unit is um and Meep frequency is f = 1 / wavelength_um.",
        "- Normal incidence, 2D Ez polarization; flat isotropic interfaces are polarization independent at normal incidence.",
        "- The lossless dielectric extends through the lower PML to approximate a semi-infinite transmitted medium.",
        "- Flat Ti is compared with a semi-infinite Fresnel interface; wavelengths above the Ti Rakić upper validity bound are warnings, not pass/fail samples.",
        "",
        "## Pass / Fail Criteria",
        "",
        f"- lossless interface tolerance: {args.lossless_tol:g}",
        f"- flat Ti Fresnel tolerance: {args.ti_tol:g}",
        f"- flat Ti transmission tolerance: {args.ti_trans_tol:g}",
        "",
        "## Answers",
        "",
        f"1. Current transmission formula sign: selected definition is `{selected}`. "
        f"Lossless max error for raw T = {lossless_summary['raw_T_max_abs_error']:.3e}; "
        f"for signed T = {lossless_summary['signed_T_max_abs_error']:.3e}.",
        f"2. Flat Ti Meep vs Fresnel: `{_status(ti_summary['passed'])}` over valid Ti wavelengths. "
        f"max |R-R_Fresnel| = {ti_summary['max_abs_R_error_valid']:.3e}; "
        f"max |A-A_Fresnel| = {ti_summary['max_abs_A_error_valid']:.3e}; "
        f"max |T_selected| = {ti_summary['max_abs_T_selected_valid']:.3e}.",
        f"3. Must src/simulation.py be modified before microstructure diagnostics? "
        f"{'Yes' if src_must_change else 'No further modification is required'} based on this run. "
        f"Public T matches selected definition: {public_matches}.",
        f"4. Old results to recompute: {old_recompute}",
        "",
        "## Numerical Verification",
        "",
        f"- lossless interface: `{_status(lossless_summary['passed'])}`; "
        f"max |R-R_th| = {lossless_summary['max_abs_R_error']:.3e}, "
        f"max |T-T_th| = {lossless_summary['max_abs_T_error']:.3e}, "
        f"max |A_selected| = {lossless_summary['max_abs_A_selected']:.3e}.",
        f"- flat Ti warning samples outside material validity: {ti_summary['warning_wavelength_count']}.",
        "",
        "## Verified Conclusions",
        "",
        f"- The lossless n={args.lossless_n:g} interface selects `{selected}` transmission normalization.",
        "- The CSV files retain `input_flux_raw`, `reflection_flux_raw`, `transmission_flux_raw`, `T_raw`, `T_signed`, `T_selected`, and `A_selected`.",
        f"- The flat Ti comparison is {_status(ti_summary['passed'])} within the Rakić validity range sampled by this run.",
        "",
        "## Hypotheses Still Not Proven",
        "",
        "- Whether slanted grooves fail because they lack an absorbing cavity mode is not answered by D00; D00 only validates the measurement pipeline.",
        "- Whether the structure mainly redistributes directionality rather than total emissivity requires a separate angular-incidence absorptance diagnostic.",
        "- Whether oxidation or multiscale roughness dominates real femtosecond-laser Ti emissivity remains a modeling hypothesis.",
        "",
        "## Needs Experiment Or Higher-Fidelity Model",
        "",
        "- Ti optical constants near and beyond 12.398 um should be checked against FTIR or ellipsometry for the processed sample.",
        "- Oxide layers, rounded sidewalls, 3D finite grooves, and polarization averaging require follow-up models.",
        "",
        "## Output Files",
        "",
        f"- lossless CSV: `{paths['lossless_csv']}`",
        f"- flat Ti CSV: `{paths['ti_csv']}`",
        f"- lossless PNG: `{paths['lossless_png']}`",
        f"- flat Ti PNG: `{paths['ti_png']}`",
        f"- report: `{paths['report']}`",
    ]
    ensure_dir(paths["report"].parent)
    paths["report"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths["report"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--wavelength_min_um", type=float, default=DEFAULTS["wavelength_min_um"])
    p.add_argument("--wavelength_max_um", type=float, default=DEFAULTS["wavelength_max_um"])
    p.add_argument("--nfreq", type=int, default=DEFAULTS["nfreq"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--period_um", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--lossless_n", type=float, default=DEFAULTS["lossless_n"])
    p.add_argument("--lossless_tol", type=float, default=DEFAULTS["lossless_tol"])
    p.add_argument("--ti_tol", type=float, default=DEFAULTS["ti_tol"])
    p.add_argument("--ti_trans_tol", type=float, default=DEFAULTS["ti_trans_tol"])
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.wavelength_min_um <= 0 or args.wavelength_min_um >= args.wavelength_max_um:
        raise ValueError("wavelength_min_um must be positive and less than wavelength_max_um")
    if args.nfreq <= 1:
        raise ValueError("nfreq must be > 1")

    logger = setup_logger("D00_flux_sign_and_fresnel_validation")
    logger.info("=== D00_flux_sign_and_fresnel_validation ===")
    logger.info("args = %s", vars(args))
    mp.verbosity(1)

    paths = _diagnostic_paths()
    lossless_df, lossless_summary = run_lossless_interface(args, logger)
    selected_t = lossless_summary["selected_T_definition"]
    ti_df, ti_summary = run_flat_ti_validation(args, selected_t, logger)

    _write_csv(lossless_df, paths["lossless_csv"])
    _write_csv(ti_df, paths["ti_csv"])
    plot_lossless(lossless_df, paths["lossless_png"])
    plot_ti(ti_df, paths["ti_png"])
    write_report(args, lossless_summary, ti_summary, paths)

    logger.info("lossless CSV → %s", paths["lossless_csv"])
    logger.info("flat Ti CSV → %s", paths["ti_csv"])
    logger.info("lossless PNG → %s", paths["lossless_png"])
    logger.info("flat Ti PNG → %s", paths["ti_png"])
    logger.info("report → %s", paths["report"])
    logger.info("overall status: lossless=%s, flat_ti=%s",
                _status(lossless_summary["passed"]), _status(ti_summary["passed"]))
    return 0 if lossless_summary["passed"] and ti_summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
