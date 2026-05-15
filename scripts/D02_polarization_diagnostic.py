"""D02_polarization_diagnostic.py — Ez/Hz 偏振吸收诊断

研究目的
--------
比较二维周期 Ti 槽结构在 Ez 与 Hz 两种独立偏振下的光谱吸收响应，判断此前
仅使用 Ez 是否遗漏了更强的金属槽吸收模式。

物理假设
--------
1. 所有长度单位为 um，Meep 频率 f = 1 / wavelength_um；
2. 正入射、2D 周期模型；Ez 和 Hz 是二维模型中的两种独立偏振；
3. Ti 使用 Meep 内置 Rakić Drude-Lorentz 模型，超过 12.398 um 的点属于外推；
4. 不透明 Ti 基底下 emissivity_proxy = absorptance；
5. 当前仍是 2D 模型，不能把 Hz 更高自动解释为真实 3D 非偏振样品发射率更高。

通过/失败判据
-------------
1. 每个 case/polarization 的 R/T/A 和 raw flux 必须有限；
2. 输出必须保留 input_flux_raw、reflection_flux_raw、transmission_flux_raw；
3. 平面 Ti 在正入射下 Ez 与 Hz 应一致，默认要求 8-13 um 平均 A 差异
   <= flat_pol_tol；
4. 微结构偏振差异只作诊断结论，不作为求解器失败条件。
"""

from __future__ import annotations

import argparse
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

from src.geometry import build_slanted_groove_geometry
from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium
from src.simulation import run_periodic_2d_metal_spectrum


DEFAULTS = dict(
    wavelength_min_um=5.0,
    wavelength_max_um=15.0,
    nfreq=121,
    resolution=32,
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    decay_db=40.0,
    flat_pol_tol=2e-2,
    hz_enhancement_tol=2e-2,
)


def _paths() -> dict[str, Path]:
    return {
        "spectra_csv": project_path(
            "results", "diagnostics", "tables",
            "D02_polarization_spectra.csv",
        ),
        "metrics_csv": project_path(
            "results", "diagnostics", "tables",
            "D02_polarization_metrics.csv",
        ),
        "flat_png": project_path(
            "results", "diagnostics", "figures",
            "D02_flat_Ti_Ez_Hz_control.png",
        ),
        "spectra_png": project_path(
            "results", "diagnostics", "figures",
            "D02_slanted_Ez_Hz_spectra.png",
        ),
        "bar_png": project_path(
            "results", "diagnostics", "figures",
            "D02_mean_A_polarization_comparison.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D02_polarization_diagnostic_report.md",
        ),
    }


def _case_specs(include_tilt30: bool) -> list[dict]:
    specs = [
        dict(case_name="flat_ti", kind="flat", tilt_deg=np.nan,
             label="Flat Ti"),
        dict(case_name="symmetric_groove_tilt0", kind="slanted", tilt_deg=0.0,
             label="Symmetric groove, tilt=0"),
        dict(case_name="slanted_groove_tilt20", kind="slanted", tilt_deg=20.0,
             label="Slanted groove, tilt=20"),
    ]
    if include_tilt30:
        specs.append(dict(case_name="slanted_groove_tilt30", kind="slanted",
                          tilt_deg=30.0, label="Slanted groove, tilt=30"))
    return specs


def _run_one(spec: dict, polarization: str, args: argparse.Namespace, logger) -> pd.DataFrame:
    ti = get_ti_medium(
        lambda_min_um=args.wavelength_min_um,
        lambda_max_um=args.wavelength_max_um,
    )

    if spec["kind"] == "flat":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
    else:
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.top_width_um,
                bottom_width_um=args.bottom_width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=float(spec["tilt_deg"]),
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
            )

    logger.info(">>> case=%s polarization=%s", spec["case_name"], polarization)
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
        source_component=polarization,
        logger=logger,
    )

    required = ["input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"]
    missing = [k for k in required if k not in result]
    if missing:
        raise RuntimeError(f"{spec['case_name']} {polarization}: missing raw flux {missing}")

    df = pd.DataFrame({
        "case_name": spec["case_name"],
        "case_label": spec["label"],
        "polarization": polarization,
        "wavelength_um": result["wavelength_um"],
        "reflectance": result["R"],
        "transmittance": result["T"],
        "absorptance": result["A"],
        "emissivity_proxy": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "normalization_note": (
            "R = reflection_flux_raw / abs(input_flux_raw); "
            "T = -transmission_flux_raw / abs(input_flux_raw); "
            "A = 1 - R - T"
        ),
        "physics_note": (
            "2D normal-incidence polarization diagnostic; "
            "non-polarized real emission needs TE/TM averaging and 3D validation"
        ),
        "period_um": args.period_um,
        "top_width_um": 0.0 if spec["kind"] == "flat" else args.top_width_um,
        "bottom_width_um": 0.0 if spec["kind"] == "flat" else args.bottom_width_um,
        "depth_um": 0.0 if spec["kind"] == "flat" else args.depth_um,
        "tilt_deg": spec["tilt_deg"],
        "walltime_s": float(result["walltime_s"]),
    })
    return df


def _band(df: pd.DataFrame, lo: float = 8.0, hi: float = 13.0) -> pd.DataFrame:
    return df[(df["wavelength_um"] >= lo) & (df["wavelength_um"] <= hi)]


def _metrics(spectra: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    grouped = spectra.groupby(["case_name", "case_label", "polarization"], sort=False)
    base = {}
    for (case_name, case_label, pol), df in grouped:
        band = _band(df)
        if band.empty:
            raise RuntimeError(f"{case_name} {pol}: no samples in 8-13 um band")
        peak_idx = band["absorptance"].idxmax()
        finite = bool(np.all(np.isfinite(df[[
            "reflectance", "transmittance", "absorptance",
            "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw",
        ]].to_numpy())))
        has_raw_flux = all(c in df.columns for c in (
            "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw",
        ))
        base[(case_name, pol)] = dict(
            case_name=case_name,
            case_label=case_label,
            polarization=pol,
            mean_A_8_13um=float(band["absorptance"].mean()),
            peak_A_8_13um=float(band.loc[peak_idx, "absorptance"]),
            peak_wavelength_um=float(band.loc[peak_idx, "wavelength_um"]),
            min_A=float(df["absorptance"].min()),
            max_A=float(df["absorptance"].max()),
            max_abs_T=float(np.max(np.abs(df["transmittance"]))),
            has_raw_flux=has_raw_flux,
            all_finite=finite,
        )

    for case_name in spectra["case_name"].drop_duplicates():
        ez = base[(case_name, "Ez")]
        hz = base[(case_name, "Hz")]
        delta_mean = hz["mean_A_8_13um"] - ez["mean_A_8_13um"]
        delta_peak = hz["peak_A_8_13um"] - ez["peak_A_8_13um"]
        for pol in ("Ez", "Hz"):
            row = dict(base[(case_name, pol)])
            row["delta_mean_A_Hz_minus_Ez"] = delta_mean
            row["delta_peak_A_Hz_minus_Ez"] = delta_peak
            row["hz_mean_enhancement_flag"] = abs(delta_mean) > args.hz_enhancement_tol
            row["pass_or_fail"] = "PASS" if row["has_raw_flux"] and row["all_finite"] else "FAIL"
            rows.append(row)
    return pd.DataFrame(rows)


def _pivot_metric(metrics: pd.DataFrame, value: str) -> pd.DataFrame:
    return metrics.pivot_table(index=["case_name", "case_label"],
                               columns="polarization", values=value, aggfunc="first")


def _plot_flat(spectra: pd.DataFrame, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    flat = spectra[spectra["case_name"] == "flat_ti"]
    for pol in ("Ez", "Hz"):
        df = flat[flat["polarization"] == pol]
        ax.plot(df["wavelength_um"], df["absorptance"], lw=1.8, label=f"{pol} A")
    ax.axvspan(8, 13, color="orange", alpha=0.10, label="8-13 um")
    ax.axvline(TI_RAKIC_VALID_LAMBDA_UM[1], color="0.4", ls=":", lw=1.0,
               label="Ti validity upper bound")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Absorptance")
    ax.set_title("D02 flat Ti polarization control")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend()
    return save_figure(fig, out_path)


def _plot_structures(spectra: pd.DataFrame, out_path: Path) -> Path:
    cases = [c for c in spectra["case_name"].drop_duplicates() if c != "flat_ti"]
    fig, axes = plt.subplots(len(cases), 1, figsize=(7.6, 3.4 * len(cases)),
                             sharex=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        subset = spectra[spectra["case_name"] == case]
        label = str(subset["case_label"].iloc[0])
        for pol in ("Ez", "Hz"):
            df = subset[subset["polarization"] == pol]
            ax.plot(df["wavelength_um"], df["absorptance"], lw=1.6, label=pol)
        ax.axvspan(8, 13, color="orange", alpha=0.10)
        ax.axvline(TI_RAKIC_VALID_LAMBDA_UM[1], color="0.4", ls=":", lw=1.0)
        ax.set_ylabel("Absorptance")
        ax.set_title(label)
        ax.grid(True, ls=":", alpha=0.6)
        ax.legend()
    axes[-1].set_xlabel("Wavelength (um)")
    fig.suptitle("D02 structured Ti Ez/Hz spectra", y=0.995)
    return save_figure(fig, out_path)


def _plot_bar(metrics: pd.DataFrame, out_path: Path) -> Path:
    pivot = _pivot_metric(metrics, "mean_A_8_13um").reset_index()
    x = np.arange(len(pivot))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.bar(x - width / 2, pivot["Ez"], width, label="Ez")
    ax.bar(x + width / 2, pivot["Hz"], width, label="Hz")
    ax.set_xticks(x)
    ax.set_xticklabels(pivot["case_label"], rotation=18, ha="right")
    ax.set_ylabel("Mean absorptance, 8-13 um")
    ax.set_title("D02 polarization comparison")
    ax.grid(True, axis="y", ls=":", alpha=0.6)
    ax.legend()
    return save_figure(fig, out_path)


def _write_report(
    metrics: pd.DataFrame,
    args: argparse.Namespace,
    paths: dict[str, Path],
) -> Path:
    pivot_mean = _pivot_metric(metrics, "mean_A_8_13um")
    pivot_peak = _pivot_metric(metrics, "peak_A_8_13um")
    flat_delta = float(
        pivot_mean.loc[("flat_ti", "Flat Ti"), "Hz"]
        - pivot_mean.loc[("flat_ti", "Flat Ti"), "Ez"]
    )
    flat_ok = abs(flat_delta) <= args.flat_pol_tol

    structured = metrics[metrics["case_name"] != "flat_ti"].drop_duplicates("case_name")
    hz_stronger_cases = structured[
        structured["delta_mean_A_Hz_minus_Ez"] > args.hz_enhancement_tol
    ]["case_name"].tolist()
    ez_conclusion_still = len(hz_stronger_cases) == 0

    lines = [
        "# D02 Polarization Diagnostic Report",
        "",
        "## Run Configuration",
        "",
        f"- wavelength range: {args.wavelength_min_um:g}-{args.wavelength_max_um:g} um",
        f"- nfreq: {args.nfreq}",
        f"- resolution: {args.resolution} pixels/um",
        f"- geometry: P={args.period_um:g}, top={args.top_width_um:g}, bottom={args.bottom_width_um:g}, depth={args.depth_um:g} um",
        f"- include tilt=30: {args.include_tilt30}",
        "",
        "## Pass / Fail Criteria",
        "",
        "- all R/T/A and raw flux columns finite",
        "- raw flux columns retained for every wavelength sample",
        f"- flat Ti polarization consistency: |mean_A_Hz - mean_A_Ez| <= {args.flat_pol_tol:g}",
        f"- microstructure Hz enhancement flag: mean_A_Hz - mean_A_Ez > {args.hz_enhancement_tol:g}",
        "",
        "## Required Answers",
        "",
        f"1. Flat Ti Ez/Hz consistency: {'PASS' if flat_ok else 'FAIL'}; "
        f"delta mean A = {flat_delta:.3e}.",
        f"2. Does Hz show stronger microstructure absorption than Ez? "
        f"{'Yes for ' + ', '.join(hz_stronger_cases) if hz_stronger_cases else 'No clear Hz enhancement by the configured threshold'}.",
        f"3. Does the previous Ez-only 'structure enhancement is not obvious' conclusion still hold? "
        f"{'Yes under this 2D diagnostic' if ez_conclusion_still else 'Not fully; Hz should be included before making that conclusion'}.",
        "",
        "## Verified Numerical Conclusions",
        "",
    ]
    for (case_name, case_label), row in pivot_mean.iterrows():
        peak_row = pivot_peak.loc[(case_name, case_label)]
        delta = float(row["Hz"] - row["Ez"])
        lines.append(
            f"- {case_label}: mean A Ez={row['Ez']:.6f}, Hz={row['Hz']:.6f}, "
            f"delta={delta:.3e}; peak A Ez={peak_row['Ez']:.6f}, Hz={peak_row['Hz']:.6f}."
        )

    lines.extend([
        "",
        "## Hypotheses Still Not Proven",
        "",
        "- A higher Hz response in 2D would not automatically imply higher real non-polarized emission.",
        "- This diagnostic does not include 3D finite grooves, roughness, oxide layers, or angular incidence.",
        "- The result does not identify the microscopic absorption mode; field maps or mode diagnostics are needed.",
        "",
        "## Needs Experiment Or Higher-Fidelity Model",
        "",
        "- Non-polarized real samples require TE/TM averaging and 3D geometry validation.",
        "- Ti optical constants beyond the Rakić range and laser-induced oxidation need experimental confirmation.",
        "",
        "## Output Files",
        "",
        f"- spectra CSV: `{paths['spectra_csv']}`",
        f"- metrics CSV: `{paths['metrics_csv']}`",
        f"- flat control PNG: `{paths['flat_png']}`",
        f"- structured spectra PNG: `{paths['spectra_png']}`",
        f"- mean-A bar PNG: `{paths['bar_png']}`",
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
    p.add_argument("--top_width_um", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width_um", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--flat_pol_tol", type=float, default=DEFAULTS["flat_pol_tol"])
    p.add_argument("--hz_enhancement_tol", type=float,
                   default=DEFAULTS["hz_enhancement_tol"])
    p.add_argument("--include_tilt30", action="store_true",
                   help="also run optional slanted groove with tilt=30 deg")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.wavelength_min_um <= 0 or args.wavelength_min_um >= args.wavelength_max_um:
        raise ValueError("wavelength_min_um must be positive and < wavelength_max_um")
    if args.depth_um >= args.substrate_thickness_um:
        raise ValueError("depth_um must be < substrate_thickness_um")
    if args.nfreq <= 1:
        raise ValueError("nfreq must be > 1")

    logger = setup_logger("D02_polarization_diagnostic")
    logger.info("=== D02_polarization_diagnostic ===")
    logger.info("args = %s", vars(args))
    mp.verbosity(1)

    paths = _paths()
    frames = []
    for spec in _case_specs(args.include_tilt30):
        for pol in ("Ez", "Hz"):
            frames.append(_run_one(spec, pol, args, logger))
    spectra = pd.concat(frames, ignore_index=True)
    metrics = _metrics(spectra, args)

    ensure_dir(paths["spectra_csv"].parent)
    spectra.to_csv(paths["spectra_csv"], index=False)
    metrics.to_csv(paths["metrics_csv"], index=False)
    _plot_flat(spectra, paths["flat_png"])
    _plot_structures(spectra, paths["spectra_png"])
    _plot_bar(metrics, paths["bar_png"])
    _write_report(metrics, args, paths)

    flat_metrics = metrics[metrics["case_name"] == "flat_ti"].drop_duplicates("case_name")
    flat_delta = float(flat_metrics["delta_mean_A_Hz_minus_Ez"].iloc[0])
    flat_ok = abs(flat_delta) <= args.flat_pol_tol
    all_ok = bool((metrics["pass_or_fail"] == "PASS").all() and flat_ok)

    logger.info("spectra CSV → %s", paths["spectra_csv"])
    logger.info("metrics CSV → %s", paths["metrics_csv"])
    logger.info("flat PNG → %s", paths["flat_png"])
    logger.info("spectra PNG → %s", paths["spectra_png"])
    logger.info("bar PNG → %s", paths["bar_png"])
    logger.info("report → %s", paths["report"])
    logger.info("metrics summary: %s", json.dumps(
        metrics.to_dict(orient="records"), ensure_ascii=False, default=str))
    logger.info("overall status: %s", "PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
