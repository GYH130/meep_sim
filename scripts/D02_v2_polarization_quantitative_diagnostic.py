"""D02_v2 polarization diagnostic with single-wavelength quantitative metrics.

This diagnostics_v2 script avoids broadband-source endpoint artifacts by
running one Meep simulation per wavelength/polarization/case.  Quantitative
band metrics use wavelength-integrated averages over valid 8.0-12.25 um points
only.  12.5 and 13.0 um are retained as observation points and marked as
material extrapolation for the current Ti Rakić model.
"""

from __future__ import annotations

import argparse
import json
import logging
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
import meep as mp
import numpy as np
import pandas as pd

from src.geometry import build_slanted_groove_geometry
from src.io_utils import ensure_dir, project_path, save_figure
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium
from src.postprocess import (
    assess_quantitative_validity,
    material_validity_mask,
    opaque_substrate_transmission_check,
    wavelength_integrated_average,
)
from src.simulation import run_periodic_2d_metal_single_wavelength


DEFAULTS = dict(
    wavelengths_quantitative_um=[
        8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.25,
    ],
    wavelengths_observation_um=[12.5, 13.0],
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    resolution=48,
    pml_thickness_um=4.0,
    substrate_thickness_um=8.0,
    air_buffer_um=8.0,
    decay_db=60.0,
    fwidth_fraction=0.06,
    hz_enhancement_tol=0.02,
    flat_pol_tol=0.01,
    mirror_abs_tol=0.01,
    slanted_difference_tol=0.02,
)


@dataclass(frozen=True)
class CaseSpec:
    case_name: str
    case_label: str
    kind: str
    tilt_deg: float | None


def _paths() -> dict[str, Path]:
    return {
        "spectra": project_path(
            "results", "diagnostics_v2", "tables",
            "D02_v2_polarization_spectra.csv",
        ),
        "metrics": project_path(
            "results", "diagnostics_v2", "tables",
            "D02_v2_polarization_metrics.csv",
        ),
        "flat": project_path(
            "results", "diagnostics_v2", "figures", "D02_v2_flat_control.png",
        ),
        "spectra_fig": project_path(
            "results", "diagnostics_v2", "figures", "D02_v2_Ez_Hz_spectra.png",
        ),
        "unpolarized": project_path(
            "results", "diagnostics_v2", "figures",
            "D02_v2_unpolarized_proxy.png",
        ),
        "report": project_path(
            "results", "diagnostics_v2", "reports",
            "D02_v2_polarization_report.md",
        ),
        "log": project_path(
            "logs", "diagnostics_v2", "D02_v2_polarization.log",
        ),
    }


def _setup_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger("D02_v2_polarization")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler(path, encoding="utf-8")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info("Log file: %s", path)
    return logger


def _case_specs() -> dict[str, CaseSpec]:
    return {
        "flat_ti": CaseSpec("flat_ti", "Flat Ti", "flat", None),
        "symmetric_groove_tilt0": CaseSpec(
            "symmetric_groove_tilt0", "Symmetric groove, tilt=0", "groove", 0.0,
        ),
        "slanted_groove_tilt20": CaseSpec(
            "slanted_groove_tilt20", "Slanted groove, tilt=+20", "groove", 20.0,
        ),
        "mirrored_slanted_groove_tilt_minus20": CaseSpec(
            "mirrored_slanted_groove_tilt_minus20",
            "Mirrored slanted groove, tilt=-20",
            "groove",
            -20.0,
        ),
    }


def _parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return list(default)
    out = []
    for item in values:
        for token in item.replace(",", " ").split():
            out.append(float(token))
    return out


def _build_factory(spec: CaseSpec, ti, args: argparse.Namespace):
    if spec.kind == "flat":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
        return factory

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        return build_slanted_groove_geometry(
            period_x_um=args.period_um,
            top_width_um=args.top_width_um,
            bottom_width_um=args.bottom_width_um,
            depth_um=args.depth_um,
            tilt_angle_deg=float(spec.tilt_deg),
            substrate_thickness_um=substrate_thickness_um,
            y_surface=y_surface_um,
            medium_substrate=ti,
        )
    return factory


def _result_level(material_flag: str, trans_flag: str, quality_flag: str) -> str:
    if trans_flag == "FAIL" or quality_flag == "TRANSMISSION_FAILURE":
        return "FAIL"
    if material_flag != "VALID":
        return "NOT_QUANTITATIVE"
    if trans_flag != "NUMERICAL_PASS" or quality_flag != "VALID":
        return "WARNING"
    return "NUMERICAL_PASS"


def _run_one(spec: CaseSpec, polarization: str, wavelength_um: float,
             args: argparse.Namespace, logger) -> dict:
    ti = get_ti_medium(lambda_min_um=wavelength_um, lambda_max_um=wavelength_um)
    result = run_periodic_2d_metal_single_wavelength(
        geometry_factory=_build_factory(spec, ti, args),
        period_um=args.period_um,
        wavelength_um=wavelength_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        decay_db=args.decay_db,
        source_component=polarization,
        fwidth_fraction=args.fwidth_fraction,
        logger=logger,
    )
    mat_mask = material_validity_mask(
        np.array([wavelength_um]),
        TI_RAKIC_VALID_LAMBDA_UM[1],
    )
    quality = assess_quantitative_validity(
        np.array([wavelength_um]),
        np.array([result["T"]]),
        np.array([result["input_flux_raw"]]),
        mat_mask,
        result["source_mode"],
    )
    trans = opaque_substrate_transmission_check(np.array([result["T"]]))
    material_flag = str(quality["material_validity_flag"][0])
    trans_flag = str(trans["transmission_quality_flag"][0])
    quality_flag = str(quality["quality_flag"][0])
    in_valid_band = 8.0 <= wavelength_um <= 12.25
    valid_for_metric = (
        in_valid_band
        and material_flag == "VALID"
        and trans_flag == "NUMERICAL_PASS"
        and quality_flag == "VALID"
    )
    return {
        "case_name": spec.case_name,
        "case_label": spec.case_label,
        "polarization": polarization,
        "wavelength_um": float(wavelength_um),
        "R": result["R"],
        "T": result["T"],
        "A": result["A"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "material_validity_flag": material_flag,
        "transmission_quality_flag": trans_flag,
        "numerical_quality_flag": quality_flag,
        "result_level": _result_level(material_flag, trans_flag, quality_flag),
        "valid_for_quantitative_metric": bool(valid_for_metric),
        "solver_version": result["solver_version"],
        "source_mode": result["source_mode"],
        "resolution": result["resolution"],
        "pml_thickness_um": result["pml_thickness_um"],
        "substrate_thickness_um": result["substrate_thickness_um"],
        "air_buffer_um": result["air_buffer_um"],
        "decay_db": result["decay_db"],
        "fwidth_fraction": result["fwidth_fraction"],
        "material_model": "Ti_Rakic_meep_builtin",
        "period_um": args.period_um,
        "top_width_um": 0.0 if spec.kind == "flat" else args.top_width_um,
        "bottom_width_um": 0.0 if spec.kind == "flat" else args.bottom_width_um,
        "depth_um": 0.0 if spec.kind == "flat" else args.depth_um,
        "tilt_deg": np.nan if spec.kind == "flat" else spec.tilt_deg,
        "normalization_note": (
            "D00 convention: R=reflection_flux_raw/abs(input_flux_raw); "
            "T=-transmission_flux_raw/abs(input_flux_raw); A=1-R-T"
        ),
        "physics_note": (
            "2D normal-incidence polarization diagnostic. "
            "A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2 is not a 3D sample emissivity."
        ),
        "walltime_s": result["walltime_s"],
    }


def _metric_for_group(df: pd.DataFrame, args: argparse.Namespace) -> dict:
    valid = df[df["valid_for_quantitative_metric"]].copy()
    mean_a = wavelength_integrated_average(
        df["A"].to_numpy(),
        df["wavelength_um"].to_numpy(),
        8.0,
        12.25,
        valid_mask=df["valid_for_quantitative_metric"].to_numpy(dtype=bool),
    )
    avg_status = wavelength_integrated_average.last_status
    if valid.empty:
        peak_a = np.nan
        peak_wl = np.nan
    else:
        idx = valid["A"].idxmax()
        peak_a = float(valid.loc[idx, "A"])
        peak_wl = float(valid.loc[idx, "wavelength_um"])
    failed = bool((df["result_level"] == "FAIL").any())
    if failed:
        level = "FAIL"
    elif avg_status != "valid":
        level = "WARNING"
    elif len(valid) < args.min_valid_samples:
        level = "WARNING"
    elif args.run_scope == "full":
        level = "NUMERICAL_PASS"
    else:
        level = "CODE_PASS"
    return {
        "mean_A_8_12p25_valid": mean_a,
        "peak_A_valid_band": peak_a,
        "peak_wavelength_valid_band_um": peak_wl,
        "valid_sample_count": int(len(valid)),
        "invalid_sample_count": int(len(df) - len(valid)),
        "band_average_status": avg_status,
        "max_abs_T": float(np.nanmax(np.abs(df["T"].to_numpy()))),
        "has_material_extrapolation_points": bool(
            (df["material_validity_flag"] != "VALID").any()
        ),
        "has_transmission_failure": bool(
            (df["transmission_quality_flag"] == "FAIL").any()
        ),
        "result_level": level,
    }


def _compute_metrics(spectra: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    base = {}
    for (case_name, case_label, pol), df in spectra.groupby(
        ["case_name", "case_label", "polarization"],
        sort=False,
    ):
        metric = _metric_for_group(df, args)
        metric.update({
            "case_name": case_name,
            "case_label": case_label,
            "polarization": pol,
            "solver_version": str(df["solver_version"].iloc[0]),
            "source_mode": str(df["source_mode"].iloc[0]),
            "resolution": int(df["resolution"].iloc[0]),
            "pml_thickness_um": float(df["pml_thickness_um"].iloc[0]),
            "substrate_thickness_um": float(df["substrate_thickness_um"].iloc[0]),
            "air_buffer_um": float(df["air_buffer_um"].iloc[0]),
            "decay_db": float(df["decay_db"].iloc[0]),
            "material_model": str(df["material_model"].iloc[0]),
            "material_validity_flag": (
                "VALID"
                if not metric["has_material_extrapolation_points"]
                else "HAS_OBSERVATION_EXTRAPOLATION"
            ),
            "numerical_quality_flag": (
                "VALID"
                if metric["result_level"] in {"NUMERICAL_PASS", "CODE_PASS"}
                else metric["result_level"]
            ),
        })
        base[(case_name, pol)] = metric

    for case_name in spectra["case_name"].drop_duplicates():
        if (case_name, "Ez") not in base or (case_name, "Hz") not in base:
            for pol in ("Ez", "Hz"):
                if (case_name, pol) not in base:
                    continue
                row = dict(base[(case_name, pol)])
                row["delta_mean_A_Hz_minus_Ez"] = np.nan
                row["delta_peak_A_Hz_minus_Ez"] = np.nan
                row["hz_enhancement_flag"] = False
                rows.append(row)
            continue
        ez = base[(case_name, "Ez")]
        hz = base[(case_name, "Hz")]
        delta_mean = hz["mean_A_8_12p25_valid"] - ez["mean_A_8_12p25_valid"]
        delta_peak = hz["peak_A_valid_band"] - ez["peak_A_valid_band"]
        hz_flag = bool(delta_mean > args.hz_enhancement_tol)
        for pol in ("Ez", "Hz"):
            row = dict(base[(case_name, pol)])
            row["delta_mean_A_Hz_minus_Ez"] = delta_mean
            row["delta_peak_A_Hz_minus_Ez"] = delta_peak
            row["hz_enhancement_flag"] = hz_flag
            rows.append(row)
    metrics = pd.DataFrame(rows)

    flat = metrics[metrics["case_name"] == "flat_ti"]
    if set(flat["polarization"]) >= {"Ez", "Hz"}:
        ez = float(flat[flat["polarization"] == "Ez"]["mean_A_8_12p25_valid"].iloc[0])
        hz = float(flat[flat["polarization"] == "Hz"]["mean_A_8_12p25_valid"].iloc[0])
        metrics["flat_Ez_Hz_abs_delta"] = abs(hz - ez)
        metrics["flat_polarization_consistency_flag"] = np.where(
            abs(hz - ez) <= args.flat_pol_tol,
            "NUMERICAL_PASS",
            "WARNING",
        )
    else:
        metrics["flat_Ez_Hz_abs_delta"] = np.nan
        metrics["flat_polarization_consistency_flag"] = "NOT_RUN"

    mirror_delta = np.nan
    mirror_rows = []
    for pol in ("Ez", "Hz"):
        plus = metrics[
            (metrics["case_name"] == "slanted_groove_tilt20")
            & (metrics["polarization"] == pol)
        ]
        minus = metrics[
            (metrics["case_name"] == "mirrored_slanted_groove_tilt_minus20")
            & (metrics["polarization"] == pol)
        ]
        if not plus.empty and not minus.empty:
            delta = abs(
                float(plus["mean_A_8_12p25_valid"].iloc[0])
                - float(minus["mean_A_8_12p25_valid"].iloc[0])
            )
            mirror_rows.append(delta)
    if mirror_rows:
        mirror_delta = float(max(mirror_rows))
    metrics["slanted_mirror_mean_A_abs_delta"] = mirror_delta
    metrics["slanted_mirror_consistency_flag"] = (
        "NOT_RUN"
        if not np.isfinite(mirror_delta)
        else ("NUMERICAL_PASS" if mirror_delta <= args.mirror_abs_tol else "WARNING")
    )
    return metrics


def _compute_unpolarized_proxy(spectra: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "case_name", "case_label", "wavelength_um", "resolution",
        "pml_thickness_um", "substrate_thickness_um", "air_buffer_um", "decay_db",
        "material_model", "material_validity_flag", "source_mode", "solver_version",
    ]
    pivot = spectra.pivot_table(
        index=keys,
        columns="polarization",
        values=["A", "valid_for_quantitative_metric"],
        aggfunc="first",
    ).reset_index()
    rows = []
    for _, row in pivot.iterrows():
        try:
            ez_a = float(row[("A", "Ez")])
            hz_a = float(row[("A", "Hz")])
            ez_valid = bool(row[("valid_for_quantitative_metric", "Ez")])
            hz_valid = bool(row[("valid_for_quantitative_metric", "Hz")])
        except Exception:
            continue
        out = {key: row[(key, "")] if (key, "") in row.index else row[key] for key in keys}
        out["A_Ez"] = ez_a
        out["A_Hz"] = hz_a
        out["A_unpolarized_2D_proxy"] = 0.5 * (ez_a + hz_a)
        out["valid_for_quantitative_metric"] = ez_valid and hz_valid
        out["physics_note"] = (
            "A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2; not a true 3D sample emissivity."
        )
        rows.append(out)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def _plot_flat(spectra: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    flat = spectra[spectra["case_name"] == "flat_ti"]
    for pol, color in [("Ez", "#4C78A8"), ("Hz", "#F58518")]:
        df = flat[flat["polarization"] == pol].sort_values("wavelength_um")
        if df.empty:
            continue
        ax.plot(df["wavelength_um"], df["A"], marker="o", color=color, label=pol)
        invalid = df[~df["valid_for_quantitative_metric"]]
        if not invalid.empty:
            ax.scatter(invalid["wavelength_um"], invalid["A"], s=65,
                       facecolors="none", edgecolors=color)
    ax.axvspan(8, 12.25, color="#E8EEF7", alpha=0.6, label="quantitative band")
    ax.axvline(TI_RAKIC_VALID_LAMBDA_UM[1], color="0.3", ls=":", lw=1,
               label="Ti validity limit")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("A")
    ax.set_title("D02_v2 flat Ti Ez/Hz control")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_spectra(spectra: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.8), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, (case, df_case) in zip(axes, spectra.groupby("case_name", sort=False)):
        for pol, color in [("Ez", "#4C78A8"), ("Hz", "#F58518")]:
            df = df_case[df_case["polarization"] == pol].sort_values("wavelength_um")
            if df.empty:
                continue
            ax.plot(df["wavelength_um"], df["A"], marker="o", lw=1.4,
                    color=color, label=pol)
            invalid = df[~df["valid_for_quantitative_metric"]]
            if not invalid.empty:
                ax.scatter(invalid["wavelength_um"], invalid["A"], s=55,
                           facecolors="none", edgecolors=color)
        ax.axvspan(8, 12.25, color="#E8EEF7", alpha=0.45)
        ax.axvline(TI_RAKIC_VALID_LAMBDA_UM[1], color="0.3", ls=":", lw=1)
        ax.set_title(case, fontsize=9)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    for ax in axes[-2:]:
        ax.set_xlabel("Wavelength (um)")
    for ax in axes[::2]:
        ax.set_ylabel("A")
    fig.suptitle("D02_v2 Ez/Hz spectra; hollow markers excluded from quantitative metrics")
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_unpolarized(proxy: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    if proxy.empty:
        ax.text(0.5, 0.5, "No paired Ez/Hz data", ha="center", va="center")
        ax.set_axis_off()
    else:
        for case, df in proxy.groupby("case_name", sort=False):
            df = df.sort_values("wavelength_um")
            ax.plot(df["wavelength_um"], df["A_unpolarized_2D_proxy"],
                    marker="o", lw=1.5, label=case)
            invalid = df[~df["valid_for_quantitative_metric"]]
            if not invalid.empty:
                ax.scatter(invalid["wavelength_um"], invalid["A_unpolarized_2D_proxy"],
                           s=55, facecolors="none", edgecolors="black")
        ax.axvspan(8, 12.25, color="#E8EEF7", alpha=0.45)
        ax.axvline(TI_RAKIC_VALID_LAMBDA_UM[1], color="0.3", ls=":", lw=1)
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("A_unpolarized_2D_proxy")
        ax.set_title("D02_v2 non-polarized 2D proxy")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    save_figure(fig, out_path)
    plt.close(fig)


def _unpolarized_metric(proxy: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    if proxy.empty:
        return pd.DataFrame()
    rows = []
    for (case, label), df in proxy.groupby(["case_name", "case_label"], sort=False):
        avg = wavelength_integrated_average(
            df["A_unpolarized_2D_proxy"].to_numpy(),
            df["wavelength_um"].to_numpy(),
            8.0,
            12.25,
            valid_mask=df["valid_for_quantitative_metric"].to_numpy(dtype=bool),
        )
        rows.append({
            "case_name": case,
            "case_label": label,
            "mean_A_unpolarized_2D_proxy_8_12p25_valid": avg,
            "band_average_status": wavelength_integrated_average.last_status,
            "valid_sample_count": int(df["valid_for_quantitative_metric"].sum()),
        })
    return pd.DataFrame(rows)


def _write_report(paths: dict[str, Path], spectra: pd.DataFrame,
                  metrics: pd.DataFrame, proxy_metrics: pd.DataFrame,
                  args: argparse.Namespace) -> None:
    ensure_dir(paths["report"].parent)
    flat_delta = (
        float(metrics["flat_Ez_Hz_abs_delta"].dropna().iloc[0])
        if "flat_Ez_Hz_abs_delta" in metrics and metrics["flat_Ez_Hz_abs_delta"].notna().any()
        else np.nan
    )
    mirror_delta = (
        float(metrics["slanted_mirror_mean_A_abs_delta"].dropna().iloc[0])
        if "slanted_mirror_mean_A_abs_delta" in metrics
        and metrics["slanted_mirror_mean_A_abs_delta"].notna().any()
        else np.nan
    )
    hz_rows = metrics[metrics["polarization"] == "Ez"]
    hz_enhanced_cases = hz_rows.loc[
        hz_rows["hz_enhancement_flag"].astype(bool), "case_name"
    ].tolist()
    slanted = metrics[
        metrics["case_name"].isin([
            "symmetric_groove_tilt0", "slanted_groove_tilt20",
            "mirrored_slanted_groove_tilt_minus20",
        ])
    ]
    if not slanted.empty:
        spread = float(slanted["mean_A_8_12p25_valid"].max()
                       - slanted["mean_A_8_12p25_valid"].min())
    else:
        spread = np.nan

    if proxy_metrics.empty:
        proxy_answer = "未生成 Ez/Hz 配对的非偏振二维代理值。"
    else:
        lines = []
        for _, row in proxy_metrics.iterrows():
            lines.append(
                f"{row['case_name']}: "
                f"{row['mean_A_unpolarized_2D_proxy_8_12p25_valid']:.4f}"
            )
        proxy_answer = "; ".join(lines)

    code_pass = not spectra.empty and not metrics.empty
    numerical_pass = bool(
        code_pass
        and args.run_scope == "full"
        and not (spectra["result_level"] == "FAIL").any()
        and (metrics["result_level"] == "NUMERICAL_PASS").all()
        and (metrics["valid_sample_count"] >= args.min_valid_samples).all()
    )
    has_failure = bool(code_pass and (spectra["result_level"] == "FAIL").any())
    has_warning = bool(code_pass and (metrics["result_level"] != "NUMERICAL_PASS").any())
    if has_failure:
        result_level = "FAIL"
    elif args.run_scope == "smoke":
        result_level = "CODE_PASS"
    elif numerical_pass:
        result_level = "NUMERICAL_PASS"
    elif has_warning:
        result_level = "WARNING"
    else:
        result_level = "CODE_PASS"
    if args.run_scope == "smoke":
        physics_statement = (
            "本次 run_scope=smoke，只能说明脚本、字段和质量标记可运行；"
            "不能作为物理结论通过。"
        )
    else:
        physics_statement = (
            "本脚本提供 NUMERICAL_PASS 级偏振诊断；是否 PHYSICS_READY 仍需 D03_v2 "
            "收敛验证。"
        )
    metrics_columns = [
        "case_name",
        "polarization",
        "mean_A_8_12p25_valid",
        "peak_A_valid_band",
        "peak_wavelength_valid_band_um",
        "delta_mean_A_Hz_minus_Ez",
        "delta_peak_A_Hz_minus_Ez",
        "hz_enhancement_flag",
        "valid_sample_count",
        "result_level",
    ]
    table_df = metrics[[c for c in metrics_columns if c in metrics.columns]].copy()
    metrics_table = _dataframe_to_markdown(table_df)

    report = f"""# D02_v2 Polarization Quantitative Diagnostic

Overall result level: **{result_level}**

## Purpose

Re-run the Ez/Hz polarization diagnostic using single-wavelength simulations
and diagnostics_v2 quality filters.  This avoids broadband endpoint source
artifacts and excludes Ti Rakić extrapolation points from quantitative metrics.

## Physical Assumptions

- Length unit is um; Meep frequency is `f = 1 / wavelength_um`.
- Source propagates from +y to -y.
- D00 flux convention is used: `R=refl/abs(input)`, `T=-trans/abs(input)`,
  `A=1-R-T`.
- Ti uses Meep built-in Rakić Drude-Lorentz model; wavelengths above
  {TI_RAKIC_VALID_LAMBDA_UM[1]:.3f} um are marked `MATERIAL_EXTRAPOLATION`.
- `A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2` is not a true 3D sample emissivity.

## Pass / Warning Criteria

- Quantitative metrics only use 8.0-12.25 um, material-valid,
  `TRANSMISSION NUMERICAL_PASS`, and `numerical_quality_flag=VALID` points.
- 12.5 and 13.0 um are observation points only.
- Hz enhancement flag uses `delta_mean_A_Hz_minus_Ez > {args.hz_enhancement_tol:g}`;
  it does not use `abs(delta)`.
- smoke runs never become `PHYSICS_READY`.

{physics_statement}

## Required Answers

1. In wavelength-integrated valid-band terms, is Hz truly higher than Ez?
   {'Yes for: ' + ', '.join(hz_enhanced_cases) if hz_enhanced_cases else 'No case exceeded the Hz enhancement tolerance in this run.'}

2. Are symmetric and slanted grooves clearly different?
   Valid-band mean-A spread across groove cases/polarizations is {spread:.4g}.

3. What is the current bare-Ti 2D non-polarized proxy?
   {proxy_answer}

4. Is the result worth entering D03 convergence validation?
   {'Yes, as a numerical target set, but not as PHYSICS_READY.' if numerical_pass else 'Not yet; first resolve failed/warning quality rows, insufficient samples, or run only smoke scope.'}

## Controls

- Flat Ti Ez/Hz mean-A difference: {flat_delta:.4g}; tolerance {args.flat_pol_tol:g}.
- Slanted +20/-20 normal-incidence mean-A max difference: {mirror_delta:.4g};
  tolerance {args.mirror_abs_tol:g}.

## Metrics Table

{metrics_table}

## Outputs

- Spectra CSV: `{paths['spectra']}`
- Metrics CSV: `{paths['metrics']}`
- Flat control: `{paths['flat']}`
- Ez/Hz spectra: `{paths['spectra_fig']}`
- Unpolarized proxy: `{paths['unpolarized']}`
- Log: `{paths['log']}`

## Run Configuration

```json
{json.dumps(vars(args), indent=2, ensure_ascii=False)}
```
"""
    paths["report"].write_text(report, encoding="utf-8")


def _format_markdown_value(value) -> str:
    if pd.isna(value):
        return "nan"
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.6g}"
    return str(value)


def _dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No metrics rows generated._"
    headers = list(df.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append(
            "| " + " | ".join(_format_markdown_value(row[col]) for col in headers) + " |"
        )
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="D02_v2 single-wavelength polarization diagnostic.")
    p.add_argument(
        "--cases", nargs="+",
        default=[
            "flat_ti",
            "symmetric_groove_tilt0",
            "slanted_groove_tilt20",
            "mirrored_slanted_groove_tilt_minus20",
        ],
        choices=sorted(_case_specs().keys()),
    )
    p.add_argument("--polarizations", nargs="+", default=["Ez", "Hz"], choices=["Ez", "Hz"])
    p.add_argument("--wavelengths_um", nargs="*", default=None)
    p.add_argument("--include_observation_points", action="store_true", default=True)
    p.add_argument("--no_observation_points", dest="include_observation_points", action="store_false")
    p.add_argument("--period_um", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--top_width_um", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width_um", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float, default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--fwidth_fraction", type=float, default=DEFAULTS["fwidth_fraction"])
    p.add_argument("--hz_enhancement_tol", type=float, default=DEFAULTS["hz_enhancement_tol"])
    p.add_argument("--flat_pol_tol", type=float, default=DEFAULTS["flat_pol_tol"])
    p.add_argument("--mirror_abs_tol", type=float, default=DEFAULTS["mirror_abs_tol"])
    p.add_argument("--slanted_difference_tol", type=float, default=DEFAULTS["slanted_difference_tol"])
    p.add_argument("--min_valid_samples", type=int, default=3)
    p.add_argument("--run_scope", choices=["smoke", "full"], default="full")
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)
    logger = _setup_logger(paths["log"])

    wavelengths = _parse_float_list(
        args.wavelengths_um,
        DEFAULTS["wavelengths_quantitative_um"],
    )
    if args.include_observation_points and args.wavelengths_um is None:
        wavelengths += DEFAULTS["wavelengths_observation_um"]
    wavelengths = sorted(set(float(w) for w in wavelengths))

    case_map = _case_specs()
    rows = []
    for case_name in args.cases:
        spec = case_map[case_name]
        for pol in args.polarizations:
            for wl in wavelengths:
                logger.info(">>> case=%s pol=%s wl=%.6g", case_name, pol, wl)
                rows.append(_run_one(spec, pol, wl, args, logger))

    spectra = pd.DataFrame(rows)
    metrics = _compute_metrics(spectra, args)
    proxy = _compute_unpolarized_proxy(spectra)
    proxy_metrics = _unpolarized_metric(proxy, args)

    spectra.to_csv(paths["spectra"], index=False)
    metrics.to_csv(paths["metrics"], index=False)
    _plot_flat(spectra, paths["flat"])
    _plot_spectra(spectra, paths["spectra_fig"])
    _plot_unpolarized(proxy, paths["unpolarized"])
    _write_report(paths, spectra, metrics, proxy_metrics, args)

    logger.info("Wrote %s", paths["spectra"])
    logger.info("Wrote %s", paths["metrics"])
    logger.info("Wrote %s", paths["report"])


if __name__ == "__main__":
    main()
