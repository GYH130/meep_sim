"""D01_solver_regression_after_flux_fix.py — D00 通量修正后的求解器回归测试

研究目的
--------
在 D00 完成通量符号验证并修正透射率定义后，对平面 Ti、对称矩形槽、
退化斜槽和非对称斜槽做统一回归，确认修正没有引入新的求解器错误。

测试对象
--------
1. flat_ti：平面 Ti；
2. rectangular_groove：P=10 um, W=4 um, D=3 um；
3. slanted_tilt0：top=bottom=4 um, D=3 um, tilt=0 deg；
4. slanted_tilt20：top=bottom=4 um, D=3 um, tilt=20 deg。

物理假设
--------
1. 所有长度单位为 um，Meep 频率 f = 1 / wavelength_um；
2. 正入射、2D、Ez 偏振；不透明 Ti 基底下 emissivity_proxy = absorptance；
3. Ti 使用 Meep 内置 Rakić Drude-Lorentz 模型，超过 12.398 um 的点标记为
   材料外推，不作为严格 Fresnel 通过/失败判据；
4. 矩形槽和 tilt=0 斜槽几何等价，但 mp.Block 与 mp.Prism 的离散化可能带来
   O(1/resolution) 差异，默认用 5% 容差判定退化一致性。

通过/失败判据
-------------
1. 每个 case 的 corrected T 不应明显为负；
2. corrected A 应主要落在 [0, 1]；
3. R/T/A 和 raw flux 不应出现 NaN/Inf；
4. 输出必须保留 input_flux_raw、reflection_flux_raw、transmission_flux_raw；
5. flat_ti 在 Ti 有效波段内应接近 Fresnel；
6. rectangular_groove 与 slanted_tilt0 应通过退化一致性检查。
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

from src.geometry import (
    build_rectangular_groove_geometry,
    build_slanted_groove_geometry,
)
from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium
from src.simulation import run_periodic_2d_metal_spectrum


DEFAULTS = dict(
    wavelength_min_um=8.0,
    wavelength_max_um=13.0,
    nfreq=21,
    resolution=32,
    period_um=10.0,
    width_um=4.0,
    depth_um=3.0,
    tilt_deg=20.0,
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    decay_db=40.0,
    flat_fresnel_tol=5e-2,
    degenerate_tol=5e-2,
    range_tol=1e-3,
    negative_t_tol=1e-3,
    mean_change_tol=1e-3,
)


def _paths() -> dict[str, Path]:
    return {
        "spectra_csv": project_path(
            "results", "diagnostics", "tables",
            "D01_solver_regression_spectra.csv",
        ),
        "metrics_csv": project_path(
            "results", "diagnostics", "tables",
            "D01_solver_regression_metrics.csv",
        ),
        "figure": project_path(
            "results", "diagnostics", "figures",
            "D01_flat_rect_slanted_comparison.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D01_solver_regression_report.md",
        ),
    }


def _ti_fresnel(wavelengths_um: np.ndarray, ti_medium) -> tuple[np.ndarray, np.ndarray]:
    lo_valid, hi_valid = TI_RAKIC_VALID_LAMBDA_UM
    r_values = []
    for wl in wavelengths_um:
        if wl < lo_valid or wl > hi_valid:
            r_values.append(np.nan)
            continue
        eps = ti_medium.epsilon(1.0 / float(wl))[0][0]
        n_ti = cmath.sqrt(eps)
        r_values.append(abs((1.0 - n_ti) / (1.0 + n_ti)) ** 2)
    r = np.asarray(r_values, dtype=float)
    return r, 1.0 - r


def _result_to_df(result: dict, *, case_name: str, reference_kind: str,
                  note: str) -> pd.DataFrame:
    required = ["input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"]
    missing = [k for k in required if k not in result]
    if missing:
        raise RuntimeError(f"{case_name}: run result missing raw flux fields: {missing}")
    df = pd.DataFrame({
        "case_name": case_name,
        "wavelength_um": result["wavelength_um"],
        "R": result["R"],
        "T": result["T"],
        "A": result["A"],
        "emissivity_proxy": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "pre_fix_A_from_raw_T": 1.0 - result["R"] - result["raw_transmittance"],
        "reference_kind": reference_kind,
        "normalization_note": (
            "R = reflection_flux_raw / abs(input_flux_raw); "
            "corrected T = -transmission_flux_raw / abs(input_flux_raw); "
            "raw_transmittance is retained for D00/D01 provenance"
        ),
        "physical_note": note,
    })
    return df


def _run_case(case_name: str, args: argparse.Namespace, logger) -> pd.DataFrame:
    ti = get_ti_medium(
        lambda_min_um=args.wavelength_min_um,
        lambda_max_um=args.wavelength_max_um,
    )

    if case_name == "flat_ti":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
        reference_kind = "fresnel_flat_ti"
        note = "flat Ti baseline; compare R/A against semi-infinite Fresnel"
    elif case_name == "rectangular_groove":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_rectangular_groove_geometry(
                period_x_um=args.period_um,
                groove_width_um=args.width_um,
                groove_depth_um=args.depth_um,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
            )
        reference_kind = "flat_ti"
        note = "symmetric rectangular groove"
    elif case_name == "slanted_tilt0":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.width_um,
                bottom_width_um=args.width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=0.0,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
            )
        reference_kind = "rectangular_groove"
        note = "degenerate slanted groove; should match rectangular groove"
    elif case_name == "slanted_tilt20":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.width_um,
                bottom_width_um=args.width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=args.tilt_deg,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
            )
        reference_kind = "slanted_tilt0"
        note = "asymmetric slanted groove"
    else:
        raise ValueError(f"unknown case_name: {case_name}")

    logger.info(">>> running case %s", case_name)
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
    df = _result_to_df(result, case_name=case_name,
                       reference_kind=reference_kind, note=note)
    df["period_um"] = args.period_um
    df["width_um"] = args.width_um
    df["depth_um"] = args.depth_um if case_name != "flat_ti" else 0.0
    df["tilt_deg"] = (
        0.0 if case_name == "slanted_tilt0"
        else args.tilt_deg if case_name == "slanted_tilt20"
        else np.nan
    )
    df["walltime_s"] = float(result["walltime_s"])
    return df


def _band_mean(df: pd.DataFrame, col: str, lo: float = 8.0, hi: float = 13.0) -> float:
    band = df[(df["wavelength_um"] >= lo) & (df["wavelength_um"] <= hi)]
    if band.empty:
        return float("nan")
    return float(band[col].mean())


def _max_abs_delta(df: pd.DataFrame, ref: pd.DataFrame, col: str) -> float:
    if not np.allclose(df["wavelength_um"].to_numpy(), ref["wavelength_um"].to_numpy()):
        raise RuntimeError("reference wavelength grid mismatch")
    return float(np.max(np.abs(df[col].to_numpy() - ref[col].to_numpy())))


def _case_metrics(
    spectra: dict[str, pd.DataFrame],
    args: argparse.Namespace,
    ti_medium,
) -> tuple[pd.DataFrame, dict]:
    rows = []
    fresnel_r, fresnel_a = _ti_fresnel(spectra["flat_ti"]["wavelength_um"].to_numpy(),
                                       ti_medium)
    fresnel = pd.DataFrame({
        "wavelength_um": spectra["flat_ti"]["wavelength_um"],
        "R": fresnel_r,
        "A": fresnel_a,
    })
    valid_fresnel = np.isfinite(fresnel["R"].to_numpy())

    summary = {}
    for name, df in spectra.items():
        finite = bool(np.all(np.isfinite(df[[
            "R", "T", "A", "input_flux_raw", "reflection_flux_raw",
            "transmission_flux_raw", "raw_transmittance", "signed_transmittance",
        ]].to_numpy())))
        has_raw_flux = all(c in df.columns for c in (
            "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw",
        ))
        max_negative_t = float(max(0.0, -df["T"].min()))
        min_a = float(df["A"].min())
        max_a = float(df["A"].max())
        mean_a = _band_mean(df, "A")
        mean_a_pre = _band_mean(df, "pre_fix_A_from_raw_T")
        mean_delta_after_fix = mean_a - mean_a_pre

        if name == "flat_ti":
            valid = valid_fresnel
            max_dR = float(np.nanmax(np.abs(
                df["R"].to_numpy()[valid] - fresnel["R"].to_numpy()[valid]
            )))
            max_dA = float(np.nanmax(np.abs(
                df["A"].to_numpy()[valid] - fresnel["A"].to_numpy()[valid]
            )))
            reference_case = "fresnel_flat_ti"
            reference_tol = args.flat_fresnel_tol
        elif name == "rectangular_groove":
            ref = spectra["flat_ti"]
            max_dR = _max_abs_delta(df, ref, "R")
            max_dA = _max_abs_delta(df, ref, "A")
            reference_case = "flat_ti"
            reference_tol = float("nan")
        elif name == "slanted_tilt0":
            ref = spectra["rectangular_groove"]
            max_dR = _max_abs_delta(df, ref, "R")
            max_dA = _max_abs_delta(df, ref, "A")
            reference_case = "rectangular_groove"
            reference_tol = args.degenerate_tol
        else:
            ref = spectra["slanted_tilt0"]
            max_dR = _max_abs_delta(df, ref, "R")
            max_dA = _max_abs_delta(df, ref, "A")
            reference_case = "slanted_tilt0"
            reference_tol = float("nan")

        base_checks_pass = (
            finite
            and has_raw_flux
            and max_negative_t <= args.negative_t_tol
            and min_a >= -args.range_tol
            and max_a <= 1.0 + args.range_tol
        )
        if name == "flat_ti":
            reference_pass = max(max_dR, max_dA) <= args.flat_fresnel_tol
        elif name == "slanted_tilt0":
            reference_pass = max(max_dR, max_dA) <= args.degenerate_tol
        else:
            reference_pass = True
        passed = bool(base_checks_pass and reference_pass)

        row = {
            "case_name": name,
            "reference_case": reference_case,
            "max_negative_T": max_negative_t,
            "min_A": min_a,
            "max_A": max_a,
            "mean_A_8_13um": mean_a,
            "mean_A_pre_fix_8_13um": mean_a_pre,
            "mean_A_delta_after_flux_fix": mean_delta_after_fix,
            "max_abs_delta_R_vs_reference": max_dR,
            "max_abs_delta_A_vs_reference": max_dA,
            "reference_tol": reference_tol,
            "has_raw_flux": has_raw_flux,
            "all_finite": finite,
            "base_checks_pass": base_checks_pass,
            "reference_check_pass": reference_pass,
            "pass_or_fail": "PASS" if passed else "FAIL",
        }
        rows.append(row)
        summary[name] = row
    return pd.DataFrame(rows), summary


def _plot(spectra: dict[str, pd.DataFrame], metrics: pd.DataFrame,
          out_path: Path, ti_medium) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)
    ax_r, ax_a = axes
    labels = {
        "flat_ti": "flat Ti",
        "rectangular_groove": "rectangular groove",
        "slanted_tilt0": "slanted tilt=0",
        "slanted_tilt20": "slanted tilt=20",
    }
    for name, df in spectra.items():
        ax_r.plot(df["wavelength_um"], df["R"], lw=1.6, label=labels[name])
        ax_a.plot(df["wavelength_um"], df["A"], lw=1.6, label=labels[name])

    fresnel_r, fresnel_a = _ti_fresnel(spectra["flat_ti"]["wavelength_um"].to_numpy(),
                                       ti_medium)
    ax_r.plot(spectra["flat_ti"]["wavelength_um"], fresnel_r,
              "k--", lw=1.2, label="flat Fresnel")
    ax_a.plot(spectra["flat_ti"]["wavelength_um"], fresnel_a,
              "k--", lw=1.2, label="flat Fresnel")
    for ax in axes:
        ax.axvspan(TI_RAKIC_VALID_LAMBDA_UM[1],
                   max(df["wavelength_um"].max() for df in spectra.values()),
                   color="orange", alpha=0.10)
        ax.grid(True, ls=":", alpha=0.6)
        ax.legend(ncol=2, fontsize=8)
    ax_r.set_ylabel("Corrected R")
    ax_a.set_ylabel("Corrected A / emissivity proxy")
    ax_a.set_xlabel("Wavelength (um)")
    ax_r.set_title("D01 solver regression after flux-sign fix")
    return save_figure(fig, out_path)


def _write_report(
    args: argparse.Namespace,
    metrics: pd.DataFrame,
    paths: dict[str, Path],
) -> Path:
    row = {r["case_name"]: r for _, r in metrics.iterrows()}
    max_mean_change = float(np.max(np.abs(metrics["mean_A_delta_after_flux_fix"])))
    significant = max_mean_change > args.mean_change_tol
    deg = row["slanted_tilt0"]
    deg_pass = bool(deg["reference_check_pass"])
    old03_recompute = (
        "Yes. 旧的 03 表是在符号修正前生成的，应废弃并重算；"
        "本 D01 也用 pre_fix_A_from_raw_T 量化了差异。"
    )
    lines = [
        "# D01 Solver Regression After Flux Fix",
        "",
        "## Run Configuration",
        "",
        f"- wavelength range: {args.wavelength_min_um:g}-{args.wavelength_max_um:g} um",
        f"- nfreq: {args.nfreq}",
        f"- resolution: {args.resolution} pixels/um",
        f"- period / width / depth: {args.period_um:g} / {args.width_um:g} / {args.depth_um:g} um",
        f"- asymmetric tilt: {args.tilt_deg:g} deg",
        f"- pml / substrate / air buffer: {args.pml_thickness_um:g} / {args.substrate_thickness_um:g} / {args.air_buffer_um:g} um",
        "",
        "## Pass / Fail Criteria",
        "",
        f"- negative T tolerance: {args.negative_t_tol:g}",
        f"- A range tolerance around [0, 1]: {args.range_tol:g}",
        f"- flat Fresnel tolerance: {args.flat_fresnel_tol:g}",
        f"- rectangular vs tilt=0 tolerance: {args.degenerate_tol:g}",
        "",
        "## Required Answers",
        "",
        f"1. Did the 8-13 um mean absorptance change significantly after the flux fix? "
        f"{'Yes' if significant else 'No'}; max |mean_A_new - mean_A_old_rawT| = {max_mean_change:.3e}.",
        f"2. Should the old 03 slanted-groove conclusion be discarded/recomputed? {old03_recompute}",
        f"3. Rectangular groove vs tilt=0 slanted groove degeneracy check: "
        f"{'PASS' if deg_pass else 'FAIL'}; max |dR| = {deg['max_abs_delta_R_vs_reference']:.3e}, "
        f"max |dA| = {deg['max_abs_delta_A_vs_reference']:.3e}.",
        "",
        "## Verified Numerical Conclusions",
        "",
    ]
    for _, r in metrics.iterrows():
        lines.append(
            f"- {r['case_name']}: {r['pass_or_fail']}; "
            f"<A>_8-13 = {r['mean_A_8_13um']:.6f}; "
            f"max_negative_T = {r['max_negative_T']:.3e}; "
            f"A range = [{r['min_A']:.3e}, {r['max_A']:.3e}]; "
            f"reference = {r['reference_case']}."
        )
    lines.extend([
        "",
        "## Hypotheses Still Not Proven",
        "",
        "- D01 does not prove whether the groove lacks a real absorbing cavity mode; it only verifies solver consistency after D00.",
        "- D01 does not determine polarization dependence; Ez-only spectra remain a modeling limitation.",
        "- D01 does not prove whether directionality changes without total emissivity gain; angular-incidence absorptance is still needed.",
        "",
        "## Needs Experiment Or Higher-Fidelity Model",
        "",
        "- Ti optical constants near 13 um and processed-surface oxidation still need experiment or improved material models.",
        "- Rounded sidewalls, oxide layers, 3D finite grooves, and rough multiscale structures are not covered by this regression.",
        "",
        "## Output Files",
        "",
        f"- spectra CSV: `{paths['spectra_csv']}`",
        f"- metrics CSV: `{paths['metrics_csv']}`",
        f"- figure: `{paths['figure']}`",
        f"- report: `{paths['report']}`",
    ])
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
    p.add_argument("--width_um", type=float, default=DEFAULTS["width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--tilt_deg", type=float, default=DEFAULTS["tilt_deg"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--flat_fresnel_tol", type=float, default=DEFAULTS["flat_fresnel_tol"])
    p.add_argument("--degenerate_tol", type=float, default=DEFAULTS["degenerate_tol"])
    p.add_argument("--range_tol", type=float, default=DEFAULTS["range_tol"])
    p.add_argument("--negative_t_tol", type=float, default=DEFAULTS["negative_t_tol"])
    p.add_argument("--mean_change_tol", type=float, default=DEFAULTS["mean_change_tol"])
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.wavelength_min_um <= 0 or args.wavelength_min_um >= args.wavelength_max_um:
        raise ValueError("wavelength_min_um must be positive and < wavelength_max_um")
    if args.nfreq <= 1:
        raise ValueError("nfreq must be > 1")
    if args.depth_um >= args.substrate_thickness_um:
        raise ValueError("depth_um must be < substrate_thickness_um")

    logger = setup_logger("D01_solver_regression_after_flux_fix")
    logger.info("=== D01_solver_regression_after_flux_fix ===")
    logger.info("args = %s", vars(args))
    mp.verbosity(1)

    paths = _paths()
    case_names = ["flat_ti", "rectangular_groove", "slanted_tilt0", "slanted_tilt20"]
    spectra = {name: _run_case(name, args, logger) for name in case_names}
    spectra_df = pd.concat(spectra.values(), ignore_index=True)

    ti = get_ti_medium(lambda_min_um=args.wavelength_min_um,
                       lambda_max_um=args.wavelength_max_um)
    metrics, summary = _case_metrics(spectra, args, ti)

    ensure_dir(paths["spectra_csv"].parent)
    spectra_df.to_csv(paths["spectra_csv"], index=False)
    metrics.to_csv(paths["metrics_csv"], index=False)
    _plot(spectra, metrics, paths["figure"], ti)
    _write_report(args, metrics, paths)

    logger.info("spectra CSV → %s", paths["spectra_csv"])
    logger.info("metrics CSV → %s", paths["metrics_csv"])
    logger.info("figure → %s", paths["figure"])
    logger.info("report → %s", paths["report"])
    logger.info("metrics summary: %s", json.dumps(summary, ensure_ascii=False, default=str))

    ok = bool((metrics["pass_or_fail"] == "PASS").all())
    logger.info("overall status: %s", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
