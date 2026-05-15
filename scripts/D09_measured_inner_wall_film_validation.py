"""D09 measured inner-wall lossy-film validation.

This script uses one narrowband Meep run per wavelength and applies the measured
n,k film only to the inside sidewalls/bottom of a 2D slanted Ti groove.  The
main quantitative observable is flux-derived total absorptance A.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/private/tmp") / "matplotlib-codex-cache"))

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.geometry import (
    build_inner_wall_film_slanted_groove_geometry,
    build_slanted_groove_geometry,
    slanted_groove_vertices,
)
from src.io_utils import ensure_dir, project_path, save_figure
from src.materials import (
    TI_RAKIC_VALID_LAMBDA_UM,
    get_measured_lossy_wall_film_medium_single_wavelength,
    get_ti_medium,
    load_measured_nk_table,
)
from src.postprocess import opaque_substrate_transmission_check, wavelength_integrated_average
from src.simulation import run_periodic_2d_metal_single_wavelength


BAND_LO = 8.1014
BAND_HI = 12.962
PERIOD_UM = 10.0
TOP_WIDTH_UM = 4.0
BOTTOM_WIDTH_UM = 4.0
DEPTH_UM = 3.0
TILT_DEG = 20.0


def _paths() -> dict[str, Path]:
    base_fig = project_path("results", "diagnostics_v2", "figures")
    return {
        "spectra": project_path(
            "results", "diagnostics_v2", "tables",
            "D09_measured_inner_wall_film_spectra.csv",
        ),
        "metrics": project_path(
            "results", "diagnostics_v2", "tables",
            "D09_measured_inner_wall_film_metrics.csv",
        ),
        "resolution": project_path(
            "results", "diagnostics_v2", "tables", "D09_resolution_check.csv",
        ),
        "report": project_path(
            "results", "diagnostics_v2", "reports",
            "D09_measured_inner_wall_film_validation_report.md",
        ),
        "log": project_path(
            "logs", "diagnostics_v2",
            "D09_measured_inner_wall_film_validation.log",
        ),
        "geom_bare": base_fig / "D09_inner_wall_film_geometry_bare.png",
        "geom_200": base_fig / "D09_inner_wall_film_geometry_200nm.png",
        "geom_250": base_fig / "D09_inner_wall_film_geometry_250nm.png",
        "geom_300": base_fig / "D09_inner_wall_film_geometry_300nm.png",
        "ez": base_fig / "D09_Ez_spectral_comparison.png",
        "hz": base_fig / "D09_Hz_spectral_comparison.png",
        "proxy": base_fig / "D09_unpolarized_proxy_comparison.png",
        "mean_vs_t": base_fig / "D09_mean_A_vs_film_thickness.png",
        "enhancement": base_fig / "D09_enhancement_over_bare.png",
        "best_e2": base_fig / "D09_best_case_E2.png",
        "best_h2": base_fig / "D09_best_case_H2.png",
        "best_overlay": base_fig / "D09_best_case_geometry_overlay.png",
    }


def _setup_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger("D09_measured_inner_wall_film")
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
    return logger


def _mode_defaults(mode: str) -> dict:
    if mode == "smoke":
        return {
            "wavelengths_um": [9.0, 10.0, 12.0],
            "polarizations": ["Ez"],
            "film_thicknesses_um": [0.0, 0.25],
            "resolution_values": [48],
            "pml_thickness_um": 4.0,
            "substrate_thickness_um": 8.0,
            "air_buffer_um": 8.0,
            "decay_db": 20.0,
        }
    if mode == "screen":
        return {
            "wavelengths_um": [
                8.1014, 8.5, 9.0, 9.5, 10.0, 10.5,
                11.0, 11.5, 12.0, 12.5, 12.962,
            ],
            "polarizations": ["Ez", "Hz"],
            "film_thicknesses_um": [0.0, 0.20, 0.25, 0.30],
            "resolution_values": [64],
            "pml_thickness_um": 4.0,
            "substrate_thickness_um": 8.0,
            "air_buffer_um": 8.0,
            "decay_db": 60.0,
        }
    if mode == "refine":
        return {
            "wavelengths_um": None,
            "polarizations": ["Ez", "Hz"],
            "film_thicknesses_um": None,
            "resolution_values": [64, 80],
            "pml_thickness_um": 4.0,
            "substrate_thickness_um": 8.0,
            "air_buffer_um": 8.0,
            "decay_db": 60.0,
        }
    raise ValueError(f"Unknown mode: {mode}")


def _parse_float_list(values: list[str] | None) -> list[float] | None:
    if values is None:
        return None
    out = []
    for item in values:
        for token in item.replace(",", " ").split():
            out.append(float(token))
    return out


def _parse_str_list(values: list[str] | None) -> list[str] | None:
    if values is None:
        return None
    out = []
    for item in values:
        out.extend(item.replace(",", " ").split())
    return out


def _select_refine_wavelengths(nk_table: pd.DataFrame, max_points: int = 31) -> list[float]:
    band = nk_table[
        (nk_table["wavelength_um"] >= BAND_LO)
        & (nk_table["wavelength_um"] <= BAND_HI)
    ].copy()
    if len(band) <= max_points:
        return band["wavelength_um"].astype(float).tolist()
    idx = np.linspace(0, len(band) - 1, max_points).round().astype(int)
    return band.iloc[idx]["wavelength_um"].astype(float).tolist()


def _film_metadata_nan(wavelength_um: float) -> dict:
    return {
        "n_film": np.nan,
        "k_film": np.nan,
        "epsilon_real_film": np.nan,
        "epsilon_imag_film": np.nan,
        "D_conductivity_film": np.nan,
        "film_interpolation_flag": "NOT_APPLICABLE",
        "film_material_validity_flag": "NOT_APPLICABLE",
        "film_model_mode": "none",
        "film_data_range_min_um": np.nan,
        "film_data_range_max_um": np.nan,
        "film_warning": "",
    }


def _make_geometry_factory(case_name: str, thickness_um: float, film_medium, args):
    import meep as mp

    ti = get_ti_medium()
    air = mp.Medium(epsilon=1.0)
    if case_name == "flat_Ti":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
        return factory
    if case_name == "bare_slanted_groove":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.top_width_um,
                bottom_width_um=args.bottom_width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=args.tilt_angle_deg,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
                medium_groove=air,
            )
        return factory
    if case_name in {
        "inner_wall_film_slanted_groove",
        "sidewalls_only_inner_wall_film_slanted_groove",
    }:
        mode = (
            "sidewalls_only"
            if case_name == "sidewalls_only_inner_wall_film_slanted_groove"
            else "sidewalls_and_bottom"
        )
        if film_medium is None:
            raise ValueError("film_medium is required for coated groove cases")

        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_inner_wall_film_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.top_width_um,
                bottom_width_um=args.bottom_width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=args.tilt_angle_deg,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                film_thickness_um=thickness_um,
                medium_substrate=ti,
                medium_film=film_medium,
                medium_groove=air,
                coating_mode=mode,
            )
        return factory
    raise ValueError(f"Unknown case_name: {case_name}")


def _plot_geometry(thickness_um: float, out_path: Path, args, mode: str = "sidewalls_and_bottom") -> None:
    y_surface = 0.0
    outer = slanted_groove_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=y_surface,
    )
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    half_p = 0.5 * args.period_um
    substrate = plt.Rectangle(
        (-half_p, y_surface - args.substrate_thickness_um),
        args.period_um,
        args.substrate_thickness_um,
        facecolor="#7F7F7F",
        edgecolor="black",
        lw=1.0,
        label="Ti substrate",
    )
    ax.add_patch(substrate)
    if thickness_um <= 0:
        ax.fill(*zip(*(outer + [outer[0]])), facecolor="white", edgecolor="#1F77B4",
                lw=1.5, label="air groove")
    else:
        ax.fill(*zip(*(outer + [outer[0]])), facecolor="#D62728", alpha=0.72,
                edgecolor="black", lw=1.1, label="measured lossy wall film")
        inner_depth = args.depth_um - (thickness_um if mode == "sidewalls_and_bottom" else 0.0)
        inner = slanted_groove_vertices(
            top_width_um=args.top_width_um - 2.0 * thickness_um,
            bottom_width_um=args.bottom_width_um - 2.0 * thickness_um,
            depth_um=inner_depth,
            tilt_angle_deg=args.tilt_angle_deg,
            y_surface=y_surface,
        )
        ax.fill(*zip(*(inner + [inner[0]])), facecolor="white",
                edgecolor="#1F77B4", lw=1.5, label="air core")
    ax.axhline(y_surface, color="black", lw=1.0)
    ax.text(-half_p + 0.1, y_surface + 0.08, "top Ti land is uncoated", fontsize=8)
    ax.set_xlim(-half_p, half_p)
    ax.set_ylim(y_surface - args.depth_um - 1.0, y_surface + 0.6)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"Inner-wall film geometry, t={thickness_um:g} um")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.2)
    save_figure(fig, out_path)
    plt.close(fig)


def _generate_geometry_plots(paths: dict[str, Path], args) -> None:
    _plot_geometry(0.0, paths["geom_bare"], args)
    _plot_geometry(0.20, paths["geom_200"], args)
    _plot_geometry(0.25, paths["geom_250"], args)
    _plot_geometry(0.30, paths["geom_300"], args)


def _run_one(case_name: str, thickness_um: float, pol: str, wl: float,
             resolution: int, args, logger) -> dict:
    film_medium = None
    film_meta = _film_metadata_nan(wl)
    if case_name in {
        "inner_wall_film_slanted_groove",
        "sidewalls_only_inner_wall_film_slanted_groove",
    }:
        film_medium, meta = get_measured_lossy_wall_film_medium_single_wavelength(
            wl,
            args.nk_csv,
            allow_extrapolation=False,
        )
        film_meta = {
            "n_film": meta["n"],
            "k_film": meta["k"],
            "epsilon_real_film": meta["epsilon_real"],
            "epsilon_imag_film": meta["epsilon_imag"],
            "D_conductivity_film": meta["D_conductivity"],
            "film_interpolation_flag": meta["interpolation_flag"],
            "film_material_validity_flag": "VALID",
            "film_model_mode": meta["model_mode"],
            "film_data_range_min_um": meta["data_range_um"][0],
            "film_data_range_max_um": meta["data_range_um"][1],
            "film_warning": meta["warning"],
        }

    result = run_periodic_2d_metal_single_wavelength(
        geometry_factory=_make_geometry_factory(case_name, thickness_um, film_medium, args),
        period_um=args.period_um,
        wavelength_um=wl,
        resolution=resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        decay_db=args.decay_db,
        source_component=pol,
        fwidth_fraction=args.fwidth_fraction,
        solver_version="inner_wall_film_single_wavelength_v1",
        source_mode="single_wavelength_narrowband",
        logger=logger,
    )
    trans = opaque_substrate_transmission_check(np.array([result["T"]]))
    trans_flag = str(trans["transmission_quality_flag"][0])
    ti_flag = (
        "VALID"
        if wl <= TI_RAKIC_VALID_LAMBDA_UM[1]
        else "TI_RAKIC_EXTRAPOLATION"
    )
    material_flag = (
        "VALID"
        if ti_flag == "VALID" and film_meta["film_material_validity_flag"] in {"VALID", "NOT_APPLICABLE"}
        else ";".join(
            flag for flag in [ti_flag, film_meta["film_material_validity_flag"]]
            if flag not in {"VALID", "NOT_APPLICABLE"}
        )
    )
    numerical_flag = trans_flag
    valid_metric = trans_flag == "NUMERICAL_PASS" and film_meta["film_material_validity_flag"] != "FAIL"
    row = {
        "mode": args.mode,
        "case_name": case_name,
        "polarization": pol,
        "film_thickness_um": float(thickness_um),
        "coating_mode": (
            "sidewalls_only"
            if case_name == "sidewalls_only_inner_wall_film_slanted_groove"
            else ("none" if thickness_um == 0 else "sidewalls_and_bottom")
        ),
        "wavelength_um": wl,
        "R": result["R"],
        "T": result["T"],
        "A": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "transmission_quality_flag": trans_flag,
        "numerical_quality_flag": numerical_flag,
        "material_validity_flag": material_flag,
        "ti_material_validity_flag": ti_flag,
        "valid_for_band_metric": bool(valid_metric),
        "solver_version": result["solver_version"],
        "source_mode": result["source_mode"],
        "resolution": result["resolution"],
        "pml_thickness_um": result["pml_thickness_um"],
        "substrate_thickness_um": result["substrate_thickness_um"],
        "air_buffer_um": result["air_buffer_um"],
        "decay_db": result["decay_db"],
        "fwidth_fraction": result["fwidth_fraction"],
        "material_model": "Ti_Rakic + measured_lossy_wall_film_conductivity",
        "period_um": args.period_um,
        "top_width_um": args.top_width_um,
        "bottom_width_um": args.bottom_width_um,
        "depth_um": args.depth_um,
        "tilt_angle_deg": args.tilt_angle_deg,
        "film_grid_points": float(thickness_um * resolution),
        "normalization_note": (
            "D00 convention: R=reflection_flux_raw/abs(input_flux_raw); "
            "T=-transmission_flux_raw/abs(input_flux_raw); A=1-R-T"
        ),
        "walltime_s": result["walltime_s"],
    }
    row.update(film_meta)
    return row


def _compute_metrics(spectra: pd.DataFrame, args) -> pd.DataFrame:
    if spectra.empty:
        return pd.DataFrame()
    rows = []
    groups = ["case_name", "polarization", "film_thickness_um", "coating_mode", "resolution"]
    for keys, df in spectra.groupby(groups, sort=False):
        case, pol, thick, coating, resolution = keys
        valid = df["valid_for_band_metric"].to_numpy(dtype=bool)
        avg = wavelength_integrated_average(
            df["A"].to_numpy(),
            df["wavelength_um"].to_numpy(),
            BAND_LO,
            BAND_HI,
            valid_mask=valid,
        )
        status = wavelength_integrated_average.last_status
        valid_df = df[df["valid_for_band_metric"]]
        if valid_df.empty:
            peak_a = np.nan
            peak_wl = np.nan
        else:
            idx = valid_df["A"].idxmax()
            peak_a = float(valid_df.loc[idx, "A"])
            peak_wl = float(valid_df.loc[idx, "wavelength_um"])
        rows.append({
            "case_name": case,
            "polarization": pol,
            "film_thickness_um": thick,
            "coating_mode": coating,
            "resolution": resolution,
            "mean_A_data_covered_band": avg,
            "peak_A": peak_a,
            "peak_wavelength_um": peak_wl,
            "valid_sample_count": int(valid.sum()),
            "invalid_sample_count": int((~valid).sum()),
            "band_average_status": status,
            "min_film_grid_points": float(df["film_grid_points"].min()),
            "max_abs_T": float(np.nanmax(np.abs(df["T"].to_numpy()))),
            "any_transmission_failure": bool((df["transmission_quality_flag"] == "FAIL").any()),
            "any_ti_extrapolation": bool((df["ti_material_validity_flag"] != "VALID").any()),
        })
    metrics = pd.DataFrame(rows)
    bare = metrics[metrics["case_name"] == "bare_slanted_groove"]
    for idx, row in metrics.iterrows():
        ref = bare[
            (bare["polarization"] == row["polarization"])
            & (bare["resolution"] == row["resolution"])
        ]
        if ref.empty or pd.isna(row["mean_A_data_covered_band"]):
            abs_enh = np.nan
            rel_enh = np.nan
        else:
            ref_a = float(ref["mean_A_data_covered_band"].iloc[0])
            abs_enh = float(row["mean_A_data_covered_band"] - ref_a)
            rel_enh = float(abs_enh / ref_a) if ref_a != 0 else np.nan
        metrics.loc[idx, "enhancement_over_bare_absolute"] = abs_enh
        metrics.loc[idx, "enhancement_over_bare_relative"] = rel_enh
        if pd.isna(row["mean_A_data_covered_band"]):
            candidate = "NOT_QUANTITATIVE"
        elif row["mean_A_data_covered_band"] >= 0.50:
            candidate = "HIGH_EMISSIVITY_CANDIDATE"
        elif row["mean_A_data_covered_band"] < 0.20:
            candidate = "LOW_ABSORPTION"
        else:
            candidate = "MODERATE_ABSORPTION"
        metrics.loc[idx, "emissivity_candidate_class"] = candidate
        if pd.isna(abs_enh):
            enh_class = "NOT_QUANTITATIVE"
        elif abs_enh < 0.03:
            enh_class = "LIMITED_ENHANCEMENT"
        elif abs_enh <= 0.10:
            enh_class = "WORTH_OPTIMIZING"
        else:
            enh_class = "SIGNIFICANT_CANDIDATE"
        metrics.loc[idx, "enhancement_class"] = enh_class
    proxy_rows = []
    for (case, thick, coating, resolution), df in spectra.groupby(
        ["case_name", "film_thickness_um", "coating_mode", "resolution"],
        sort=False,
    ):
        pivot = df.pivot_table(
            index="wavelength_um",
            columns="polarization",
            values=["A", "valid_for_band_metric"],
            aggfunc="first",
        )
        if not {"Ez", "Hz"}.issubset(set(pivot["A"].columns)):
            continue
        wavelengths = pivot.index.to_numpy(dtype=float)
        proxy_a = 0.5 * (
            pivot[("A", "Ez")].to_numpy(dtype=float)
            + pivot[("A", "Hz")].to_numpy(dtype=float)
        )
        valid = (
            pivot[("valid_for_band_metric", "Ez")].to_numpy(dtype=bool)
            & pivot[("valid_for_band_metric", "Hz")].to_numpy(dtype=bool)
        )
        avg = wavelength_integrated_average(proxy_a, wavelengths, BAND_LO, BAND_HI, valid)
        proxy_rows.append({
            "case_name": case,
            "polarization": "unpolarized_2D_proxy",
            "film_thickness_um": thick,
            "coating_mode": coating,
            "resolution": resolution,
            "mean_A_data_covered_band": avg,
            "peak_A": float(np.nanmax(proxy_a)) if proxy_a.size else np.nan,
            "peak_wavelength_um": float(wavelengths[np.nanargmax(proxy_a)]) if proxy_a.size else np.nan,
            "valid_sample_count": int(valid.sum()),
            "invalid_sample_count": int((~valid).sum()),
            "band_average_status": wavelength_integrated_average.last_status,
            "min_film_grid_points": float(df["film_grid_points"].min()),
            "max_abs_T": float(np.nanmax(np.abs(df["T"].to_numpy()))),
            "any_transmission_failure": bool((df["transmission_quality_flag"] == "FAIL").any()),
            "any_ti_extrapolation": bool((df["ti_material_validity_flag"] != "VALID").any()),
            "physics_note": (
                "A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2; not a true 3D sample emissivity."
            ),
        })
    if proxy_rows:
        metrics = pd.concat([metrics, pd.DataFrame(proxy_rows)], ignore_index=True)
        bare = metrics[metrics["case_name"] == "bare_slanted_groove"]
        for idx, row in metrics[metrics["polarization"] == "unpolarized_2D_proxy"].iterrows():
            ref = bare[
                (bare["polarization"] == "unpolarized_2D_proxy")
                & (bare["resolution"] == row["resolution"])
            ]
            if ref.empty or pd.isna(row["mean_A_data_covered_band"]):
                continue
            ref_a = float(ref["mean_A_data_covered_band"].iloc[0])
            abs_enh = float(row["mean_A_data_covered_band"] - ref_a)
            metrics.loc[idx, "enhancement_over_bare_absolute"] = abs_enh
            metrics.loc[idx, "enhancement_over_bare_relative"] = (
                abs_enh / ref_a if ref_a != 0 else np.nan
            )
            metrics.loc[idx, "emissivity_candidate_class"] = (
                "HIGH_EMISSIVITY_CANDIDATE"
                if row["mean_A_data_covered_band"] >= 0.50
                else ("LOW_ABSORPTION" if row["mean_A_data_covered_band"] < 0.20 else "MODERATE_ABSORPTION")
            )
    return metrics


def _resolution_check(metrics: pd.DataFrame, args) -> pd.DataFrame:
    if args.mode != "refine" or metrics.empty:
        return pd.DataFrame([{
            "mode": args.mode,
            "resolution_check_status": "NOT_RUN",
            "note": "Resolution comparison is required only for refine mode.",
        }])
    rows = []
    for keys, df in metrics.groupby(["case_name", "polarization", "film_thickness_um"], sort=False):
        if not {64, 80}.issubset(set(df["resolution"].astype(int))):
            continue
        a64 = float(df[df["resolution"] == 64]["mean_A_data_covered_band"].iloc[0])
        a80 = float(df[df["resolution"] == 80]["mean_A_data_covered_band"].iloc[0])
        rows.append({
            "case_name": keys[0],
            "polarization": keys[1],
            "film_thickness_um": keys[2],
            "mean_A_res64": a64,
            "mean_A_res80": a80,
            "abs_delta_mean_A": abs(a80 - a64),
            "resolution_check_status": "NUMERICAL_PASS" if abs(a80 - a64) < 0.01 else "WARNING",
        })
    return pd.DataFrame(rows)


def _plot_spectra(spectra: pd.DataFrame, pol: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    df_pol = spectra[spectra["polarization"] == pol]
    if df_pol.empty:
        ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
        ax.set_axis_off()
    else:
        for (case, thick), df in df_pol.groupby(["case_name", "film_thickness_um"], sort=False):
            label = f"{case}, t={thick:g} um"
            ax.plot(df["wavelength_um"], df["A"], marker="o", lw=1.3, label=label)
        ax.axvspan(BAND_LO, BAND_HI, color="#E8EEF7", alpha=0.45)
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("A")
        ax.set_title(f"D09 {pol} spectral comparison")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_proxy(spectra: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    pivot = spectra.pivot_table(
        index=["case_name", "film_thickness_um", "wavelength_um"],
        columns="polarization",
        values="A",
        aggfunc="first",
    ).reset_index()
    if not {"Ez", "Hz"}.issubset(set(pivot.columns)):
        ax.text(0.5, 0.5, "No paired Ez/Hz data", ha="center", va="center")
        ax.set_axis_off()
    else:
        pivot["A_unpolarized_2D_proxy"] = 0.5 * (pivot["Ez"] + pivot["Hz"])
        for (case, thick), df in pivot.groupby(["case_name", "film_thickness_um"], sort=False):
            ax.plot(
                df["wavelength_um"],
                df["A_unpolarized_2D_proxy"],
                marker="o",
                lw=1.3,
                label=f"{case}, t={thick:g} um",
            )
        ax.axvspan(BAND_LO, BAND_HI, color="#E8EEF7", alpha=0.45)
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("A_unpolarized_2D_proxy")
        ax.set_title("D09 2D non-polarized proxy")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_metric_bars(metrics: pd.DataFrame, out_path: Path, value_col: str,
                      title: str, ylabel: str) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    df = metrics[
        metrics["case_name"].isin([
            "bare_slanted_groove", "inner_wall_film_slanted_groove",
        ])
    ].copy()
    if df.empty or value_col not in df:
        ax.text(0.5, 0.5, "No metric data", ha="center", va="center")
        ax.set_axis_off()
    else:
        labels = [
            f"{r.case_name}\n{r.polarization}\nt={r.film_thickness_um:g}"
            for r in df.itertuples()
        ]
        ax.bar(np.arange(len(df)), df[value_col].astype(float), color="#4C78A8")
        ax.set_xticks(np.arange(len(df)))
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25)
    save_figure(fig, out_path)
    plt.close(fig)


def _field_components(pol: str):
    import meep as mp

    if pol == "Ez":
        return [mp.Ez], [mp.Hx, mp.Hy], mp.Ez
    if pol == "Hz":
        return [mp.Ex, mp.Ey], [mp.Hz], mp.Hz
    raise ValueError(pol)


def _sum_abs2(arrays: list[np.ndarray]) -> np.ndarray:
    min_shape = tuple(min(a.shape[i] for a in arrays) for i in range(arrays[0].ndim))
    slices = tuple(slice(0, n) for n in min_shape)
    out = np.zeros(min_shape, dtype=float)
    for arr in arrays:
        out += np.abs(np.asarray(arr)[slices]) ** 2
    return out


def _plot_field_map(data: np.ndarray, title: str, path: Path, extent) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    vmax = np.nanpercentile(data, 99.5) if np.isfinite(data).any() else None
    im = ax.imshow(data.T, origin="lower", aspect="auto", extent=extent,
                   cmap="inferno", vmin=0, vmax=vmax)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    save_figure(fig, path)
    plt.close(fig)


def _write_placeholder_field_figures(paths: dict[str, Path], reason: str, args) -> dict:
    for key, title in [
        ("best_e2", "D09 |E|^2 field snapshot not run"),
        ("best_h2", "D09 |H|^2 field snapshot not run"),
    ]:
        fig, ax = plt.subplots(figsize=(7.0, 4.0))
        ax.text(0.5, 0.5, reason, ha="center", va="center", wrap=True)
        ax.set_title(title)
        ax.set_axis_off()
        save_figure(fig, paths[key])
        plt.close(fig)
    _plot_geometry(0.25, paths["best_overlay"], args)
    return {"field_snapshot_status": "NOT_RUN", "field_snapshot_note": reason}


def _run_best_field_snapshot(spectra: pd.DataFrame, metrics: pd.DataFrame,
                             paths: dict[str, Path], args, logger) -> dict:
    if args.mode == "smoke" or args.skip_fields:
        return _write_placeholder_field_figures(
            paths,
            "Field snapshots are skipped for smoke runs or when --skip-fields is set.",
            args,
        )
    candidates = metrics[
        (metrics["case_name"] == "inner_wall_film_slanted_groove")
        & metrics["enhancement_over_bare_absolute"].notna()
    ].sort_values("enhancement_over_bare_absolute", ascending=False)
    if candidates.empty:
        return _write_placeholder_field_figures(paths, "No coated best case available.", args)
    best = candidates.iloc[0]
    rows = spectra[
        (spectra["case_name"] == best["case_name"])
        & (spectra["polarization"] == best["polarization"])
        & (spectra["film_thickness_um"] == best["film_thickness_um"])
        & (spectra["resolution"] == best["resolution"])
    ].copy()
    rows = rows.sort_values("A", ascending=False)
    if rows.empty:
        return _write_placeholder_field_figures(paths, "No spectral row for best case.", args)
    row = rows.iloc[0]

    import meep as mp

    wl = float(row["wavelength_um"])
    pol = str(row["polarization"])
    thickness = float(row["film_thickness_um"])
    resolution = int(row["resolution"])
    ti = get_ti_medium()
    film, _ = get_measured_lossy_wall_film_medium_single_wavelength(wl, args.nk_csv)
    bottom_buffer = args.pml_thickness_um
    cell_y = (
        2 * args.pml_thickness_um
        + args.air_buffer_um
        + args.substrate_thickness_um
        + bottom_buffer
    )
    y_top_edge = 0.5 * cell_y
    y_top_pml_inner = y_top_edge - args.pml_thickness_um
    y_surface = y_top_pml_inner - args.air_buffer_um
    y_src = y_top_pml_inner - 0.25 * args.air_buffer_um
    y_refl = y_surface + 0.5 * args.air_buffer_um
    y_dft_min = y_surface - args.depth_um - 0.5
    y_dft_max = y_surface + 0.8
    y_center = 0.5 * (y_dft_min + y_dft_max)
    y_size = y_dft_max - y_dft_min
    import meep as mp
    air = mp.Medium(epsilon=1.0)
    geometry = build_inner_wall_film_slanted_groove_geometry(
        period_x_um=args.period_um,
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        substrate_thickness_um=args.substrate_thickness_um,
        y_surface=y_surface,
        film_thickness_um=thickness,
        medium_substrate=ti,
        medium_film=film,
        medium_groove=air,
    )
    e_comps, h_comps, src_c = _field_components(pol)
    fcen = 1.0 / wl
    sim = mp.Simulation(
        cell_size=mp.Vector3(args.period_um, cell_y, 0),
        boundary_layers=[mp.PML(thickness=args.pml_thickness_um, direction=mp.Y)],
        sources=[
            mp.Source(
                mp.GaussianSource(frequency=fcen, fwidth=args.fwidth_fraction * fcen, is_integrated=True),
                component=src_c,
                center=mp.Vector3(0, y_src, 0),
                size=mp.Vector3(args.period_um, 0, 0),
            )
        ],
        resolution=resolution,
        k_point=mp.Vector3(),
        geometry=geometry,
        dimensions=2,
    )
    dft = sim.add_dft_fields(
        e_comps + h_comps,
        fcen,
        0,
        1,
        center=mp.Vector3(0, y_center, 0),
        size=mp.Vector3(args.period_um, y_size, 0),
        yee_grid=True,
    )
    logger.info(">>> D09 field snapshot pol=%s wl=%.6g t=%.4g", pol, wl, thickness)
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    e2 = _sum_abs2([sim.get_dft_array(dft, comp, 0) for comp in e_comps])
    h2 = _sum_abs2([sim.get_dft_array(dft, comp, 0) for comp in h_comps])
    extent = (-0.5 * args.period_um, 0.5 * args.period_um, y_dft_min, y_dft_max)
    _plot_field_map(e2, f"|E|^2 best D09 case, {pol}, {wl:g} um", paths["best_e2"], extent)
    _plot_field_map(h2, f"|H|^2 best D09 case, {pol}, {wl:g} um", paths["best_h2"], extent)
    _plot_geometry(thickness, paths["best_overlay"], args)
    return {
        "field_snapshot_status": "COMPLETED",
        "field_snapshot_case": row["case_name"],
        "field_snapshot_polarization": pol,
        "field_snapshot_wavelength_um": wl,
        "field_snapshot_film_thickness_um": thickness,
    }


def _write_failure_report(paths: dict[str, Path], exc: Exception, args) -> None:
    ensure_dir(paths["report"].parent)
    paths["report"].write_text(
        "# D09 Measured Inner-Wall Film Validation\n\n"
        "Overall result level: **FAIL**\n\n"
        f"Failure: `{type(exc).__name__}: {exc}`\n\n"
        "No physical conclusion should be drawn from this failed run.\n",
        encoding="utf-8",
    )


def _write_report(paths: dict[str, Path], spectra: pd.DataFrame,
                  metrics: pd.DataFrame, resolution: pd.DataFrame,
                  field_info: dict, args) -> None:
    ensure_dir(paths["report"].parent)
    if args.mode == "smoke":
        level = "CODE_PASS" if not spectra.empty else "FAIL"
    elif args.mode == "screen":
        screen_quality = (
            not spectra.empty
            and not metrics.empty
            and int(spectra["resolution"].min()) >= 64
            and float(spectra["decay_db"].min()) >= 60.0
        )
        level = "NUMERICAL_SCREENING" if screen_quality else "WARNING"
    else:
        pass_resolution = (
            not resolution.empty
            and (resolution.get("resolution_check_status") == "NUMERICAL_PASS").all()
        )
        no_fail = not spectra.empty and not (spectra["transmission_quality_flag"] == "FAIL").any()
        level = "NUMERICAL_PASS" if pass_resolution and no_fail else "WARNING"

    best = metrics[
        metrics["case_name"] == "inner_wall_film_slanted_groove"
    ].sort_values("enhancement_over_bare_absolute", ascending=False).head(1)
    if best.empty:
        best_text = "No quantitative coated-groove metric is available."
    else:
        row = best.iloc[0]
        best_text = (
            f"Best coated case: t={row['film_thickness_um']:g} um, "
            f"{row['polarization']}, mean_A={row['mean_A_data_covered_band']:.4g}, "
            f"absolute enhancement={row['enhancement_over_bare_absolute']:.4g}, "
            f"relative enhancement={row['enhancement_over_bare_relative']:.4g}."
        )
    high = metrics[
        metrics["mean_A_data_covered_band"].fillna(-1) >= 0.50
    ] if not metrics.empty else pd.DataFrame()
    high_text = "yes" if not high.empty else "no"
    table = _markdown_table(metrics.head(30))
    paths["report"].write_text(
        "# D09 Measured Inner-Wall Film Validation\n\n"
        f"Overall result level: **{level}**\n\n"
        "## Scope\n\n"
        "- 2D periodic slanted Ti groove model only; not a final 3D slanted-pore prediction.\n"
        "- The wall film is named `measured_lossy_wall_film`; no chemical identity is inferred.\n"
        "- Ti substrate uses the existing Rakić model.\n"
        "- Film n,k are used only through narrowband single-wavelength conductivity media.\n"
        "- Flux-derived total absorptance A is the quantitative criterion.\n\n"
        "## Required Answers\n\n"
        "1. Was the measured n,k table read, interpolated, and used in Meep?\n"
        f"   {'Yes, for coated cases.' if not spectra.empty and spectra['n_film'].notna().any() else 'No coated case was completed in this run.'}\n\n"
        "2. What is the actual quantitative evaluation interval?\n"
        f"   {BAND_LO:.4f}-{BAND_HI:.3f} um, constrained by the measured film data overlap.\n\n"
        "3. Which 200/250/300 nm film improves most?\n"
        f"   {best_text}\n\n"
        "4. Absolute and relative enhancement over bare slanted groove:\n"
        f"   {best_text}\n\n"
        "5. Is there a mean_A >= 0.50 high-emissivity candidate?\n"
        f"   {high_text}.\n\n"
        "6. If enhancement remains limited, likely causes include insufficient film n,k loss, "
        "limited coated area, weak geometric coupling, 2D-vs-3D mismatch, or missing multiscale roughness.\n\n"
        "7. Next steps should be selected from angle-resolved validation, period/depth/width scan, "
        "film-thickness optimization, and experimental cross-section/composition validation.\n\n"
        "8. The film is not automatically TiO2 unless the Excel chemistry is independently confirmed.\n\n"
        "## Quality Notes\n\n"
        f"- Field snapshot status: {field_info.get('field_snapshot_status', 'UNKNOWN')}\n"
        f"- Minimum resolution in this run: {int(spectra['resolution'].min()) if not spectra.empty else 'NA'}; "
        f"minimum decay_db: {float(spectra['decay_db'].min()) if not spectra.empty else 'NA'}.\n"
        "- A screen-mode run is only labeled NUMERICAL_SCREENING when resolution>=64 "
        "and decay_db>=60.\n"
        "- Conductivity film field maps are qualitative only; no D04 Lorentz/Drude absorbed-power "
        "formula is used as formal film absorption evidence.\n"
        "- Points with opaque-substrate transmission FAIL are excluded from band metrics.\n"
        "- Wavelengths above Ti Rakić validity are explicitly flagged in the spectra table.\n\n"
        "## Metrics Preview\n\n"
        f"{table}\n\n"
        "## Output Files\n\n"
        f"- Spectra: `{paths['spectra']}`\n"
        f"- Metrics: `{paths['metrics']}`\n"
        f"- Resolution check: `{paths['resolution']}`\n"
        f"- Log: `{paths['log']}`\n",
        encoding="utf-8",
    )


def _markdown_table(df: pd.DataFrame) -> str:
    cols = [
        "case_name", "polarization", "film_thickness_um",
        "mean_A_data_covered_band", "enhancement_over_bare_absolute",
        "emissivity_candidate_class", "enhancement_class",
    ]
    if df.empty:
        return "_No metric rows._"
    df = df[[c for c in cols if c in df.columns]].copy()
    lines = [
        "| " + " | ".join(df.columns) + " |",
        "| " + " | ".join(["---"] * len(df.columns)) + " |",
    ]
    for _, row in df.iterrows():
        vals = []
        for col in df.columns:
            val = row[col]
            vals.append("nan" if pd.isna(val) else (f"{val:.6g}" if isinstance(val, float) else str(val)))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate measured inner-wall lossy film.")
    parser.add_argument("--mode", choices=["smoke", "screen", "refine"], default="smoke")
    parser.add_argument("--nk-csv", type=Path, required=True)
    parser.add_argument("--best-thickness-um", type=float, default=None)
    parser.add_argument("--wavelengths-um", nargs="*", default=None)
    parser.add_argument("--polarizations", nargs="*", choices=["Ez", "Hz"], default=None)
    parser.add_argument("--film-thickness-um", nargs="*", default=None)
    parser.add_argument("--resolution", nargs="*", type=int, default=None)
    parser.add_argument("--period-um", type=float, default=PERIOD_UM)
    parser.add_argument("--top-width-um", type=float, default=TOP_WIDTH_UM)
    parser.add_argument("--bottom-width-um", type=float, default=BOTTOM_WIDTH_UM)
    parser.add_argument("--depth-um", type=float, default=DEPTH_UM)
    parser.add_argument("--tilt-angle-deg", type=float, default=TILT_DEG)
    parser.add_argument("--pml-thickness-um", type=float, default=None)
    parser.add_argument("--substrate-thickness-um", type=float, default=None)
    parser.add_argument("--air-buffer-um", type=float, default=None)
    parser.add_argument("--decay-db", type=float, default=None)
    parser.add_argument("--fwidth-fraction", type=float, default=0.06)
    parser.add_argument("--include-sidewalls-only", action="store_true")
    parser.add_argument("--skip-fields", action="store_true")
    return parser


def _configure_args(args) -> None:
    defaults = _mode_defaults(args.mode)
    nk_table = load_measured_nk_table(args.nk_csv)
    wavelengths = _parse_float_list(args.wavelengths_um)
    if wavelengths is None:
        wavelengths = defaults["wavelengths_um"]
    if wavelengths is None:
        wavelengths = _select_refine_wavelengths(nk_table)
    polarizations = _parse_str_list(args.polarizations)
    if polarizations is None:
        polarizations = defaults["polarizations"]
    thicknesses = _parse_float_list(args.film_thickness_um)
    if thicknesses is None:
        if args.mode == "refine":
            if args.best_thickness_um is None:
                raise ValueError("--best-thickness-um is required for refine mode")
            thicknesses = [0.0, float(args.best_thickness_um)]
        else:
            thicknesses = defaults["film_thicknesses_um"]
    resolutions = args.resolution if args.resolution is not None else defaults["resolution_values"]
    args.wavelengths_um = sorted(set(float(x) for x in wavelengths))
    args.polarizations = polarizations
    args.film_thicknesses_um = sorted(set(float(x) for x in thicknesses))
    args.resolution_values = [int(x) for x in resolutions]
    args.pml_thickness_um = args.pml_thickness_um or defaults["pml_thickness_um"]
    args.substrate_thickness_um = args.substrate_thickness_um or defaults["substrate_thickness_um"]
    args.air_buffer_um = args.air_buffer_um or defaults["air_buffer_um"]
    args.decay_db = args.decay_db if args.decay_db is not None else defaults["decay_db"]


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)
    logger = _setup_logger(paths["log"])
    try:
        defaults = _mode_defaults(args.mode)
        args.pml_thickness_um = args.pml_thickness_um or defaults["pml_thickness_um"]
        args.substrate_thickness_um = (
            args.substrate_thickness_um or defaults["substrate_thickness_um"]
        )
        args.air_buffer_um = args.air_buffer_um or defaults["air_buffer_um"]
        args.decay_db = args.decay_db if args.decay_db is not None else defaults["decay_db"]
        _generate_geometry_plots(paths, args)
        _configure_args(args)
        rows = []
        for resolution in args.resolution_values:
            for pol in args.polarizations:
                for wl in args.wavelengths_um:
                    logger.info(">>> flat_Ti pol=%s wl=%.6g res=%d", pol, wl, resolution)
                    rows.append(_run_one("flat_Ti", 0.0, pol, wl, resolution, args, logger))
                    logger.info(">>> bare_slanted pol=%s wl=%.6g res=%d", pol, wl, resolution)
                    rows.append(_run_one("bare_slanted_groove", 0.0, pol, wl, resolution, args, logger))
                    for thickness in args.film_thicknesses_um:
                        if thickness <= 0:
                            continue
                        logger.info(
                            ">>> coated t=%.4g pol=%s wl=%.6g res=%d",
                            thickness, pol, wl, resolution,
                        )
                        rows.append(
                            _run_one(
                                "inner_wall_film_slanted_groove",
                                thickness,
                                pol,
                                wl,
                                resolution,
                                args,
                                logger,
                            )
                        )
                    if args.include_sidewalls_only and 0.25 in args.film_thicknesses_um:
                        rows.append(
                            _run_one(
                                "sidewalls_only_inner_wall_film_slanted_groove",
                                0.25,
                                pol,
                                wl,
                                resolution,
                                args,
                                logger,
                            )
                        )
        spectra = pd.DataFrame(rows)
        metrics = _compute_metrics(spectra, args)
        resolution = _resolution_check(metrics, args)
        spectra.to_csv(paths["spectra"], index=False)
        metrics.to_csv(paths["metrics"], index=False)
        resolution.to_csv(paths["resolution"], index=False)
        _plot_spectra(spectra, "Ez", paths["ez"])
        _plot_spectra(spectra, "Hz", paths["hz"])
        _plot_proxy(spectra, paths["proxy"])
        _plot_metric_bars(
            metrics,
            paths["mean_vs_t"],
            "mean_A_data_covered_band",
            "D09 mean A vs inner-wall film thickness",
            "mean A",
        )
        _plot_metric_bars(
            metrics,
            paths["enhancement"],
            "enhancement_over_bare_absolute",
            "D09 enhancement over bare slanted groove",
            "absolute enhancement",
        )
        field_info = _run_best_field_snapshot(spectra, metrics, paths, args, logger)
        _write_report(paths, spectra, metrics, resolution, field_info, args)
        logger.info("Wrote %s", paths["spectra"])
        logger.info("Wrote %s", paths["metrics"])
        logger.info("Wrote %s", paths["report"])
    except Exception as exc:
        logger.exception("D09 failed")
        _write_failure_report(paths, exc, args)
        raise


if __name__ == "__main__":
    main()
