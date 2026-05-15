"""D10 planar full-surface measured-film capability screen.

This diagnostics_v2 script tests the material-capability upper bound of the
user-provided measured_lossy_wall_film by covering a flat Ti backplane with a
uniform film.  It is not a final laser-processed microstructure model.
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
import meep as mp
import numpy as np
import pandas as pd

from src.geometry import build_planar_film_on_ti_geometry
from src.io_utils import ensure_dir, project_path, save_figure
from src.materials import (
    TI_RAKIC_VALID_LAMBDA_UM,
    get_measured_lossy_wall_film_medium_single_wavelength,
    get_ti_medium,
    load_measured_nk_table,
)
from src.postprocess import opaque_substrate_transmission_check, wavelength_integrated_average
from src.simulation import run_periodic_2d_metal_single_wavelength


BAND_STRICT_LO = 8.1014
BAND_STRICT_HI = TI_RAKIC_VALID_LAMBDA_UM[1]
BAND_EXTENDED_HI = 12.962
DEFAULT_PERIOD_UM = 10.0
SCREEN_THICKNESSES_UM = [0.00, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00, 3.00, 4.00]


def _paths() -> dict[str, Path]:
    fig_dir = project_path("results", "diagnostics_v2", "figures")
    return {
        "spectra": project_path("results", "diagnostics_v2", "tables", "D10_planar_full_surface_film_spectra.csv"),
        "metrics": project_path("results", "diagnostics_v2", "tables", "D10_planar_full_surface_film_metrics.csv"),
        "resolution": project_path("results", "diagnostics_v2", "tables", "D10_planar_full_surface_film_resolution_check.csv"),
        "report": project_path("results", "diagnostics_v2", "reports", "D10_planar_full_surface_film_capability_report.md"),
        "log": project_path("logs", "diagnostics_v2", "D10_planar_full_surface_film_capability.log"),
        "geom_bare": fig_dir / "D10_geometry_bare_Ti.png",
        "geom_300": fig_dir / "D10_geometry_film_300nm.png",
        "geom_1um": fig_dir / "D10_geometry_film_1um.png",
        "geom_3um": fig_dir / "D10_geometry_film_3um.png",
        "ez": fig_dir / "D10_Ez_thickness_spectra.png",
        "hz": fig_dir / "D10_Hz_thickness_spectra.png",
        "mean": fig_dir / "D10_mean_A_vs_film_thickness.png",
        "best": fig_dir / "D10_best_film_vs_bare_Ti.png",
        "control": fig_dir / "D10_Ez_Hz_planar_control_check.png",
    }


def _setup_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger("D10_planar_full_surface_film")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(path, encoding="utf-8")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def _mode_defaults(mode: str) -> dict:
    if mode == "smoke":
        return dict(
            thicknesses=[0.0, 0.30, 1.00],
            wavelengths=[9.0, 10.0, 12.0],
            polarizations=["Ez"],
            resolutions=[48],
            pml=4.0,
            substrate=8.0,
            air=8.0,
            decay=40.0,
            fwidth=0.06,
        )
    if mode == "screen":
        return dict(
            thicknesses=SCREEN_THICKNESSES_UM,
            wavelengths=None,
            polarizations=["Ez", "Hz"],
            resolutions=[64],
            pml=4.0,
            substrate=8.0,
            air=8.0,
            decay=60.0,
            fwidth=0.06,
        )
    if mode == "refine":
        return dict(
            thicknesses=None,
            wavelengths=None,
            polarizations=["Ez", "Hz"],
            resolutions=[64, 80],
            pml=4.0,
            substrate=8.0,
            air=8.0,
            decay=60.0,
            fwidth=0.06,
        )
    raise ValueError(f"Unknown mode: {mode}")


def _parse_float_list(values: list[str] | None) -> list[float] | None:
    if values is None:
        return None
    out = []
    for item in values:
        for token in item.replace(",", " ").split():
            out.append(float(token))
    return out


def _nearest(table_wl: np.ndarray, target: float, *, upper: float | None = None) -> float:
    wl = table_wl if upper is None else table_wl[table_wl <= upper]
    if wl.size == 0:
        raise ValueError(f"No measured wavelength <= {upper} for target {target}")
    return float(wl[np.argmin(np.abs(wl - target))])


def _select_screen_wavelengths(nk_table: pd.DataFrame) -> list[float]:
    wl = nk_table["wavelength_um"].to_numpy(dtype=float)
    strict = wl[(wl >= BAND_STRICT_LO) & (wl <= BAND_STRICT_HI)]
    if len(strict) < 25:
        raise ValueError("Measured n,k table has fewer than 25 strict-band samples.")
    chosen = set(strict[np.linspace(0, len(strict) - 1, 25).round().astype(int)].astype(float))
    for target in [8.1014, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, BAND_STRICT_HI]:
        chosen.add(_nearest(strict, target, upper=BAND_STRICT_HI))
    for target in [12.5, BAND_EXTENDED_HI]:
        extended = wl[(wl > BAND_STRICT_HI) & (wl <= BAND_EXTENDED_HI)]
        if extended.size:
            chosen.add(_nearest(extended, target, upper=BAND_EXTENDED_HI))
    return sorted(chosen)


def _select_refine_wavelengths(nk_table: pd.DataFrame, max_points: int = 61) -> list[float]:
    strict = nk_table[
        (nk_table["wavelength_um"] >= BAND_STRICT_LO)
        & (nk_table["wavelength_um"] <= BAND_STRICT_HI)
    ]["wavelength_um"].to_numpy(dtype=float)
    if len(strict) <= max_points:
        return strict.astype(float).tolist()
    idx = np.linspace(0, len(strict) - 1, max_points).round().astype(int)
    return strict[idx].astype(float).tolist()


def _read_top3_screen_thicknesses(metrics_path: Path) -> list[float]:
    if not metrics_path.exists():
        raise FileNotFoundError(
            f"Refine requires existing screen metrics: {metrics_path}"
        )
    metrics = pd.read_csv(metrics_path)
    candidates = metrics[
        (metrics["polarization"].isin(["Ez", "Hz"]))
        & (metrics["film_thickness_um"] > 0)
        & metrics["mean_A_8p1014_12p398_strict"].notna()
    ].copy()
    if candidates.empty:
        raise ValueError("No valid coated-film screen metrics available for refine.")
    ranked = (
        candidates.groupby("film_thickness_um")["mean_A_8p1014_12p398_strict"]
        .max()
        .sort_values(ascending=False)
    )
    return [0.0] + [float(x) for x in ranked.head(3).index]


def _configure_args(args) -> None:
    defaults = _mode_defaults(args.mode)
    table = load_measured_nk_table(args.nk_csv)
    args.thicknesses_um = _parse_float_list(args.thickness_um) or defaults["thicknesses"]
    if args.mode == "refine" and args.thicknesses_um is None:
        args.thicknesses_um = _read_top3_screen_thicknesses(_paths()["metrics"])
    args.wavelengths_um = _parse_float_list(args.wavelengths_um) or defaults["wavelengths"]
    if args.wavelengths_um is None:
        args.wavelengths_um = _select_refine_wavelengths(table) if args.mode == "refine" else _select_screen_wavelengths(table)
    args.polarizations = args.polarizations or defaults["polarizations"]
    args.resolutions = args.resolution or defaults["resolutions"]
    args.pml_thickness_um = args.pml_thickness_um or defaults["pml"]
    args.substrate_thickness_um = args.substrate_thickness_um or defaults["substrate"]
    args.air_buffer_um = args.air_buffer_um or defaults["air"]
    args.decay_db = args.decay_db if args.decay_db is not None else defaults["decay"]
    args.fwidth_fraction = args.fwidth_fraction if args.fwidth_fraction is not None else defaults["fwidth"]
    args.thicknesses_um = sorted(set(float(x) for x in args.thicknesses_um))
    args.wavelengths_um = sorted(set(float(x) for x in args.wavelengths_um))
    args.resolutions = [int(x) for x in args.resolutions]
    args.strict_planned_point_count = sum(BAND_STRICT_LO <= wl <= BAND_STRICT_HI for wl in args.wavelengths_um)


def _film_meta_empty() -> dict:
    return {
        "n_film": np.nan,
        "k_film": np.nan,
        "epsilon_real_film": np.nan,
        "epsilon_imag_film": np.nan,
        "D_conductivity_film": np.nan,
        "film_interpolation_flag": "NOT_APPLICABLE",
        "material_validity_flag_film": "NOT_APPLICABLE",
        "film_model_mode": "none",
        "film_data_range_min_um": np.nan,
        "film_data_range_max_um": np.nan,
        "film_warning": "",
    }


def _geometry_factory(thickness: float, film_medium, args):
    ti = get_ti_medium()
    if thickness == 0:
        film_medium = mp.Medium(epsilon=1.0)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        return build_planar_film_on_ti_geometry(
            period_x_um=args.period_um,
            substrate_thickness_um=substrate_thickness_um,
            film_thickness_um=thickness,
            y_surface=y_surface_um,
            medium_substrate=ti,
            medium_film=film_medium,
        )

    return factory


def _transmission_flag(T: float) -> str:
    flag = str(opaque_substrate_transmission_check(np.array([T]))["transmission_quality_flag"][0])
    if flag == "NUMERICAL_PASS":
        return "NUMERICAL_PASS"
    if flag == "WARNING":
        return "WARNING_EXCLUDED_FROM_STRICT_METRIC"
    return "FAIL_EXCLUDED_FROM_ALL_METRICS"


def _run_one(thickness: float, pol: str, wl: float, resolution: int, args, logger) -> dict:
    film_medium = None
    film_meta = _film_meta_empty()
    if thickness > 0:
        film_medium, meta = get_measured_lossy_wall_film_medium_single_wavelength(
            wl, args.nk_csv, allow_extrapolation=False,
        )
        film_meta = {
            "n_film": meta["n"],
            "k_film": meta["k"],
            "epsilon_real_film": meta["epsilon_real"],
            "epsilon_imag_film": meta["epsilon_imag"],
            "D_conductivity_film": meta["D_conductivity"],
            "film_interpolation_flag": meta["interpolation_flag"],
            "material_validity_flag_film": "VALID",
            "film_model_mode": meta["model_mode"],
            "film_data_range_min_um": meta["data_range_um"][0],
            "film_data_range_max_um": meta["data_range_um"][1],
            "film_warning": meta["warning"],
        }

    result = run_periodic_2d_metal_single_wavelength(
        geometry_factory=_geometry_factory(thickness, film_medium, args),
        period_um=args.period_um,
        wavelength_um=wl,
        resolution=resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        decay_db=args.decay_db,
        source_component=pol,
        fwidth_fraction=args.fwidth_fraction,
        solver_version="planar_full_surface_film_single_wavelength_v1",
        source_mode="single_wavelength_narrowband",
        logger=logger,
    )
    ti_flag = "VALID" if wl <= BAND_STRICT_HI else "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING"
    trans_flag = _transmission_flag(result["T"])
    finite = all(np.isfinite(result[x]) for x in ["R", "T", "A", "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"])
    common_valid = (
        finite
        and trans_flag == "NUMERICAL_PASS"
        and result["source_mode"] == "single_wavelength_narrowband"
        and film_meta["material_validity_flag_film"] in {"VALID", "NOT_APPLICABLE"}
    )
    included_strict = common_valid and ti_flag == "VALID" and BAND_STRICT_LO <= wl <= BAND_STRICT_HI
    included_extended = common_valid and BAND_STRICT_LO <= wl <= BAND_EXTENDED_HI
    row = {
        "mode": args.mode,
        "case_name": "bare_flat_Ti" if thickness == 0 else "flat_Ti_with_measured_lossy_wall_film",
        "film_thickness_um": thickness,
        "polarization": pol,
        "wavelength_um": wl,
        "R": result["R"],
        "T": result["T"],
        "A": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "material_validity_flag_ti": ti_flag,
        "transmission_quality_flag": trans_flag,
        "included_in_strict_metric": bool(included_strict),
        "included_in_extended_metric": bool(included_extended),
        "finite_values_flag": "VALID" if finite else "NAN_OR_INF",
        "resolution": result["resolution"],
        "pml_thickness_um": result["pml_thickness_um"],
        "substrate_thickness_um": result["substrate_thickness_um"],
        "air_buffer_um": result["air_buffer_um"],
        "decay_db": result["decay_db"],
        "fwidth_fraction": result["fwidth_fraction"],
        "solver_version": result["solver_version"],
        "source_mode": result["source_mode"],
        "period_um": args.period_um,
        "walltime_s": result["walltime_s"],
        "normalization_note": "D00 convention: R=refl/abs(input), T=-trans/abs(input), A=1-R-T.",
    }
    row.update(film_meta)
    return row


def _route_class(mean_a: float) -> str:
    if not np.isfinite(mean_a):
        return "NOT_QUANTITATIVE"
    if mean_a < 0.30:
        return "MATERIAL_ROUTE_FAIL_FOR_HIGH_EMISSIVITY"
    if mean_a < 0.60:
        return "MATERIAL_HAS_LOSS_BUT_REQUIRES_RESONANT_ARCHITECTURE"
    if mean_a < 0.80:
        return "MATERIAL_ROUTE_WORTH_MICROSTRUCTURE_COUPLING"
    return "HIGH_EMISSIVITY_BASE_LAYER_CANDIDATE"


def _fabrication_class(thickness: float) -> str:
    if thickness <= 0.30:
        return "THIN_MODIFIED_LAYER_RANGE"
    if thickness <= 1.00:
        return "POSSIBLE_THICK_MODIFIED_LAYER_RANGE_NEEDS_EXPERIMENT_CONFIRMATION"
    return "CAPABILITY_UPPER_BOUND_ONLY_NOT_YET_PROCESS_JUSTIFIED"


def _compute_metrics(spectra: pd.DataFrame, args) -> pd.DataFrame:
    rows = []
    group_cols = ["film_thickness_um", "polarization", "resolution"]
    for (thickness, pol, resolution), df in spectra.groupby(group_cols, sort=False):
        strict_mask = df["included_in_strict_metric"].to_numpy(dtype=bool)
        strict_count = int(strict_mask.sum())
        enough = strict_count >= 0.9 * args.strict_planned_point_count
        mean_strict = wavelength_integrated_average(
            df["A"].to_numpy(),
            df["wavelength_um"].to_numpy(),
            BAND_STRICT_LO,
            BAND_STRICT_HI,
            valid_mask=strict_mask,
        )
        strict_status = wavelength_integrated_average.last_status
        if not enough:
            mean_strict = np.nan
            numerical_status = "INCOMPLETE_NUMERICAL_COVERAGE"
        elif strict_status != "valid":
            numerical_status = strict_status.upper()
        else:
            numerical_status = "NUMERICAL_PASS"
        ext_mask = df["included_in_extended_metric"].to_numpy(dtype=bool)
        mean_extended = wavelength_integrated_average(
            df["A"].to_numpy(),
            df["wavelength_um"].to_numpy(),
            BAND_STRICT_LO,
            BAND_EXTENDED_HI,
            valid_mask=ext_mask,
        )
        valid_strict = df[df["included_in_strict_metric"]]
        if valid_strict.empty:
            peak_a = np.nan
            peak_wl = np.nan
        else:
            idx = valid_strict["A"].idxmax()
            peak_a = float(valid_strict.loc[idx, "A"])
            peak_wl = float(valid_strict.loc[idx, "wavelength_um"])
        rows.append({
            "film_thickness_um": thickness,
            "polarization": pol,
            "resolution": resolution,
            "mean_A_8p1014_12p398_strict": mean_strict,
            "mean_A_8p1014_12p962_extended": mean_extended,
            "peak_A_strict": peak_a,
            "peak_wavelength_strict_um": peak_wl,
            "valid_point_count_strict": strict_count,
            "total_planned_point_count_strict": args.strict_planned_point_count,
            "valid_point_count_extended": int(ext_mask.sum()),
            "numerical_status": numerical_status,
            "max_abs_T": float(np.nanmax(np.abs(df["T"].to_numpy()))),
            "has_ti_extended_warning": bool((df["material_validity_flag_ti"] != "VALID").any()),
            "high_emissivity_route_classification": _route_class(mean_strict),
            "fabrication_relevance_classification": _fabrication_class(thickness),
        })
    metrics = pd.DataFrame(rows)
    for idx, row in metrics.iterrows():
        bare = metrics[
            (metrics["film_thickness_um"] == 0)
            & (metrics["polarization"] == row["polarization"])
            & (metrics["resolution"] == row["resolution"])
        ]
        if bare.empty or not np.isfinite(row["mean_A_8p1014_12p398_strict"]):
            continue
        bare_a = float(bare["mean_A_8p1014_12p398_strict"].iloc[0])
        abs_enh = float(row["mean_A_8p1014_12p398_strict"] - bare_a)
        metrics.loc[idx, "enhancement_over_bare_absolute_strict"] = abs_enh
        metrics.loc[idx, "enhancement_over_bare_relative_strict"] = abs_enh / bare_a if bare_a else np.nan
    control_rows = []
    for (thickness, resolution), df in metrics.groupby(["film_thickness_um", "resolution"], sort=False):
        if {"Ez", "Hz"}.issubset(set(df["polarization"])):
            ez = float(df[df["polarization"] == "Ez"]["mean_A_8p1014_12p398_strict"].iloc[0])
            hz = float(df[df["polarization"] == "Hz"]["mean_A_8p1014_12p398_strict"].iloc[0])
            delta = abs(ez - hz) if np.isfinite(ez) and np.isfinite(hz) else np.nan
            control_rows.append((thickness, resolution, delta))
    metrics["ez_hz_mean_A_abs_delta"] = np.nan
    metrics["ez_hz_planar_control_flag"] = "NOT_RUN"
    for thickness, resolution, delta in control_rows:
        mask = (metrics["film_thickness_um"] == thickness) & (metrics["resolution"] == resolution)
        metrics.loc[mask, "ez_hz_mean_A_abs_delta"] = delta
        metrics.loc[mask, "ez_hz_planar_control_flag"] = "NUMERICAL_PASS" if delta < 0.01 else "WARNING"
    return metrics


def _resolution_check(metrics: pd.DataFrame, args) -> pd.DataFrame:
    if args.mode != "refine":
        return pd.DataFrame([{"mode": args.mode, "resolution_check_status": "NOT_RUN"}])
    rows = []
    for (thickness, pol), df in metrics.groupby(["film_thickness_um", "polarization"], sort=False):
        if not {64, 80}.issubset(set(df["resolution"].astype(int))):
            continue
        a64 = float(df[df["resolution"] == 64]["mean_A_8p1014_12p398_strict"].iloc[0])
        a80 = float(df[df["resolution"] == 80]["mean_A_8p1014_12p398_strict"].iloc[0])
        rows.append({
            "film_thickness_um": thickness,
            "polarization": pol,
            "mean_A_res64": a64,
            "mean_A_res80": a80,
            "abs_delta_mean_A": abs(a80 - a64),
            "resolution_check_status": "NUMERICAL_PASS" if abs(a80 - a64) < 0.01 else "WARNING",
        })
    return pd.DataFrame(rows)


def _plot_geometry(thickness: float, out_path: Path, args) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    half = args.period_um / 2.0
    ax.add_patch(plt.Rectangle((-half, -args.substrate_thickness_um), args.period_um, args.substrate_thickness_um, facecolor="#808080", edgecolor="black", label="Ti substrate"))
    if thickness > 0:
        ax.add_patch(plt.Rectangle((-half, 0), args.period_um, thickness, facecolor="#D62728", alpha=0.75, edgecolor="black", label="measured_lossy_wall_film"))
    ax.axhline(0, color="black", lw=1.0)
    ax.set_xlim(-half, half)
    ax.set_ylim(-0.4, max(0.5, thickness + 0.4))
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y relative to Ti surface (um)")
    ax.set_title(f"D10 planar full-surface film, t={thickness:g} um")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.2)
    save_figure(fig, out_path)
    plt.close(fig)


def _make_plots(spectra: pd.DataFrame, metrics: pd.DataFrame, paths: dict[str, Path], args) -> None:
    _plot_geometry(0.0, paths["geom_bare"], args)
    _plot_geometry(0.30, paths["geom_300"], args)
    _plot_geometry(1.00, paths["geom_1um"], args)
    _plot_geometry(3.00, paths["geom_3um"], args)
    for pol, out_path in [("Ez", paths["ez"]), ("Hz", paths["hz"])]:
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        df_pol = spectra[spectra["polarization"] == pol]
        if df_pol.empty:
            ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
            ax.set_axis_off()
        else:
            for thickness, df in df_pol.groupby("film_thickness_um", sort=True):
                ax.plot(df["wavelength_um"], df["A"], marker="o", lw=1.1, label=f"{thickness:g} um")
            ax.axvspan(BAND_STRICT_LO, BAND_STRICT_HI, color="#E8EEF7", alpha=0.45, label="strict")
            ax.axvspan(BAND_STRICT_HI, BAND_EXTENDED_HI, color="#FCE8D5", alpha=0.35, label="extended Ti extrap.")
            ax.set_xlabel("Wavelength (um)")
            ax.set_ylabel("A")
            ax.set_title(f"D10 {pol} spectra by film thickness")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7, ncol=2)
        save_figure(fig, out_path)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for pol, df in metrics.groupby("polarization", sort=False):
        ax.plot(df["film_thickness_um"], df["mean_A_8p1014_12p398_strict"], marker="o", label=pol)
    ax.axhline(0.60, color="#D62728", ls="--", lw=1, label="0.60")
    ax.axhline(0.80, color="#9467BD", ls=":", lw=1, label="0.80")
    ax.set_xlabel("Film thickness (um)")
    ax.set_ylabel("Mean A strict")
    ax.set_title("D10 strict-band mean A vs full-surface film thickness")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    save_figure(fig, paths["mean"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    if metrics.empty:
        ax.text(0.5, 0.5, "No metrics", ha="center", va="center")
        ax.set_axis_off()
    else:
        best = metrics.sort_values("mean_A_8p1014_12p398_strict", ascending=False).head(4)
        ax.bar([f"{r.film_thickness_um:g} um\n{r.polarization}" for r in best.itertuples()], best["mean_A_8p1014_12p398_strict"])
        ax.set_ylabel("Mean A strict")
        ax.set_title("D10 best film thicknesses vs bare Ti")
        ax.grid(axis="y", alpha=0.25)
    save_figure(fig, paths["best"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ctrl = metrics.dropna(subset=["ez_hz_mean_A_abs_delta"]).drop_duplicates(["film_thickness_um", "resolution"])
    if ctrl.empty:
        ax.text(0.5, 0.5, "Ez/Hz paired data not available", ha="center", va="center")
        ax.set_axis_off()
    else:
        ax.plot(ctrl["film_thickness_um"], ctrl["ez_hz_mean_A_abs_delta"], marker="o")
        ax.axhline(0.01, color="red", ls="--", label="0.01 tolerance")
        ax.set_xlabel("Film thickness (um)")
        ax.set_ylabel("|mean A Ez - Hz|")
        ax.set_title("D10 planar Ez/Hz control check")
        ax.grid(alpha=0.25)
        ax.legend()
    save_figure(fig, paths["control"])
    plt.close(fig)


def _md_table(metrics: pd.DataFrame) -> str:
    cols = [
        "film_thickness_um", "polarization", "mean_A_8p1014_12p398_strict",
        "mean_A_8p1014_12p962_extended", "enhancement_over_bare_absolute_strict",
        "high_emissivity_route_classification", "fabrication_relevance_classification",
    ]
    if metrics.empty:
        return "_No metrics._"
    df = metrics[[c for c in cols if c in metrics.columns]].copy()
    lines = ["| " + " | ".join(df.columns) + " |", "| " + " | ".join(["---"] * len(df.columns)) + " |"]
    for _, row in df.iterrows():
        vals = []
        for col in df.columns:
            val = row[col]
            vals.append("nan" if pd.isna(val) else (f"{val:.6g}" if isinstance(val, float) else str(val)))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def _overall_level(metrics: pd.DataFrame, resolution: pd.DataFrame, args) -> str:
    if metrics.empty:
        return "FAIL"
    if args.mode == "smoke":
        return "CODE_PASS"
    if args.mode == "screen":
        formal = min(args.resolutions) >= 64 and args.decay_db >= 60
        return "NUMERICAL_SCREENING" if formal else "WARNING"
    pass_res = not resolution.empty and (resolution.get("resolution_check_status") == "NUMERICAL_PASS").all()
    no_bad = (metrics["numerical_status"] == "NUMERICAL_PASS").all()
    controls = metrics["ez_hz_planar_control_flag"].isin(["NUMERICAL_PASS", "NOT_RUN"]).all()
    return "NUMERICAL_PASS" if pass_res and no_bad and controls else "WARNING"


def _write_report(paths: dict[str, Path], metrics: pd.DataFrame, resolution: pd.DataFrame, args) -> None:
    ensure_dir(paths["report"].parent)
    level = _overall_level(metrics, resolution, args)
    best = metrics[metrics["mean_A_8p1014_12p398_strict"].notna()].sort_values("mean_A_8p1014_12p398_strict", ascending=False).head(1)
    if best.empty:
        best_text = "No strict quantitative mean is available."
        route_text = "No route decision can be made."
        threshold60 = False
        threshold80 = False
    else:
        r = best.iloc[0]
        best_text = f"best strict mean A = {r['mean_A_8p1014_12p398_strict']:.4g} at t={r['film_thickness_um']:g} um, {r['polarization']}."
        threshold60 = bool(r["mean_A_8p1014_12p398_strict"] >= 0.60)
        threshold80 = bool(r["mean_A_8p1014_12p398_strict"] >= 0.80)
        if threshold80:
            route_text = "High-emissivity base-layer candidate found; next add asymmetric microstructure while preserving absorption."
        elif threshold60:
            route_text = "Material route is worth coupling to laser-processable microstructures; do not jump directly to directionality optimization."
        else:
            route_text = "Current measured_lossy_wall_film is not a high-emissivity主体材料 route by itself; consider other layers, composite absorbers, or resonant architectures."
    control_warning = ""
    ctrl = metrics.dropna(subset=["ez_hz_mean_A_abs_delta"])
    if not ctrl.empty:
        max_delta = float(ctrl["ez_hz_mean_A_abs_delta"].max())
        control_warning = f"Max planar Ez/Hz strict mean delta = {max_delta:.4g}."
    report = f"""# D10 Planar Full-Surface Film Capability Screen

Overall result level: **{level}**

## Purpose

This is a material-capability upper-bound test for `measured_lossy_wall_film` on a flat Ti backplane.  It is not the final laser-processed sample geometry and does not rename the film as TiO2.

## Required Answers

1. Was measured n,k used in one-wavelength Meep simulations?
   {'Yes.' if not metrics.empty else 'No completed simulations.'}

2. Best strict-band mean absorptance in 8.1014-12.398 um:
   {best_text}

3. Required thickness for best result:
   {('See best result above.' if not best.empty else 'Not available.')}

4. Is the best thickness experimentally justified?
   Classification is reported per row; t>1 um is only a capability upper bound unless experiments confirm it.

5. Thresholds:
   mean_A >= 0.60: {'yes' if threshold60 else 'no'}; mean_A >= 0.80: {'yes' if threshold80 else 'no'}.

6-8. Route decision:
   {route_text}

9. Strict vs extended:
   Strict conclusions use only 8.1014-12.398 um.  12.398-12.962 um is an extended observation region with `TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING`.

## Quality Notes

- Transmission must satisfy abs(T)<=1e-3 for inclusion in metrics.
- Coverage below 90% of planned strict points is marked `INCOMPLETE_NUMERICAL_COVERAGE`.
- {control_warning}
- Solver source mode must be `single_wavelength_narrowband`; no broadband endpoint normalization is used.

## Metrics Preview

{_md_table(metrics)}

## Outputs

- Spectra: `{paths['spectra']}`
- Metrics: `{paths['metrics']}`
- Resolution check: `{paths['resolution']}`
- Log: `{paths['log']}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="D10 planar full-surface measured-film capability screen.")
    parser.add_argument("--mode", choices=["smoke", "screen", "refine"], default="smoke")
    parser.add_argument("--nk-csv", type=Path, required=True)
    parser.add_argument("--thickness-um", nargs="*", default=None)
    parser.add_argument("--wavelengths-um", nargs="*", default=None)
    parser.add_argument("--polarizations", nargs="*", choices=["Ez", "Hz"], default=None)
    parser.add_argument("--resolution", nargs="*", type=int, default=None)
    parser.add_argument("--period-um", type=float, default=DEFAULT_PERIOD_UM)
    parser.add_argument("--pml-thickness-um", type=float, default=None)
    parser.add_argument("--substrate-thickness-um", type=float, default=None)
    parser.add_argument("--air-buffer-um", type=float, default=None)
    parser.add_argument("--decay-db", type=float, default=None)
    parser.add_argument("--fwidth-fraction", type=float, default=None)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)
    logger = _setup_logger(paths["log"])
    _configure_args(args)
    rows = []
    for resolution in args.resolutions:
        for pol in args.polarizations:
            for thickness in args.thicknesses_um:
                for wl in args.wavelengths_um:
                    logger.info(">>> t=%.4g pol=%s wl=%.6g res=%d", thickness, pol, wl, resolution)
                    rows.append(_run_one(thickness, pol, wl, resolution, args, logger))
    spectra = pd.DataFrame(rows)
    metrics = _compute_metrics(spectra, args)
    resolution = _resolution_check(metrics, args)
    spectra.to_csv(paths["spectra"], index=False)
    metrics.to_csv(paths["metrics"], index=False)
    resolution.to_csv(paths["resolution"], index=False)
    _make_plots(spectra, metrics, paths, args)
    _write_report(paths, metrics, resolution, args)
    logger.info("Wrote %s", paths["spectra"])
    logger.info("Wrote %s", paths["metrics"])
    logger.info("Wrote %s", paths["report"])


if __name__ == "__main__":
    main()
