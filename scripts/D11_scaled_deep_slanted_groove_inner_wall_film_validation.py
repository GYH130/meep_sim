"""D11 scaled deep slanted-groove inner-wall film validation.

This diagnostics_v2 script tests whether a 5x scaled 2D slanted groove with
measured_lossy_wall_film only on the inner sidewalls and bottom can improve
total absorptance in the strict Ti-valid mid-infrared band.
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


BAND_STRICT_LO = 8.1014
BAND_STRICT_HI = TI_RAKIC_VALID_LAMBDA_UM[1]
BAND_EXTENDED_HI = 12.962

DEFAULT_PERIOD_UM = 50.0
DEFAULT_TOP_WIDTH_UM = 20.0
DEFAULT_BOTTOM_WIDTH_UM = 20.0
DEFAULT_DEPTH_UM = 15.0
DEFAULT_TILT_DEG = 20.0


def _paths() -> dict[str, Path]:
    fig_dir = project_path("results", "diagnostics_v2", "figures")
    return {
        "spectra": project_path("results", "diagnostics_v2", "tables", "D11_scaled_deep_groove_inner_wall_film_spectra.csv"),
        "metrics": project_path("results", "diagnostics_v2", "tables", "D11_scaled_deep_groove_inner_wall_film_metrics.csv"),
        "resolution": project_path("results", "diagnostics_v2", "tables", "D11_scaled_deep_groove_resolution_check.csv"),
        "pml": project_path("results", "diagnostics_v2", "tables", "D11_scaled_deep_groove_pml_check.csv"),
        "report": project_path("results", "diagnostics_v2", "reports", "D11_scaled_deep_groove_inner_wall_film_validation_report.md"),
        "log": project_path("logs", "diagnostics_v2", "D11_scaled_deep_groove_inner_wall_film_validation.log"),
        "geom_bare": fig_dir / "D11_geometry_scaled_bare_groove.png",
        "geom_250": fig_dir / "D11_geometry_scaled_inner_wall_film_250nm.png",
        "geom_500": fig_dir / "D11_geometry_scaled_inner_wall_film_500nm.png",
        "ez": fig_dir / "D11_Ez_absorptance_comparison.png",
        "hz": fig_dir / "D11_Hz_absorptance_comparison.png",
        "proxy": fig_dir / "D11_unpolarized_proxy_comparison.png",
        "mean": fig_dir / "D11_mean_A_route_decision.png",
        "enhancement": fig_dir / "D11_enhancement_over_scaled_bare.png",
        "best_e2": fig_dir / "D11_best_case_E2.png",
        "best_h2": fig_dir / "D11_best_case_H2.png",
        "best_overlay": fig_dir / "D11_best_case_geometry_overlay.png",
    }


def _setup_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger("D11_scaled_deep_groove")
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


def _parse_float_list(values: list[str] | None) -> list[float] | None:
    if values is None:
        return None
    out: list[float] = []
    for item in values:
        for token in item.replace(",", " ").split():
            out.append(float(token))
    return out


def _parse_str_list(values: list[str] | None) -> list[str] | None:
    if values is None:
        return None
    out: list[str] = []
    for item in values:
        out.extend(item.replace(",", " ").split())
    return out


def _mode_defaults(mode: str) -> dict:
    if mode == "smoke":
        return dict(
            cases=["scaled_bare_slanted_groove", "scaled_inner_wall_film_500nm"],
            wavelengths=[9.0, 10.0, 12.0],
            polarizations=["Ez", "Hz"],
            resolutions=[32],
            pml_values=[4.0],
            substrate=25.0,
            air=12.0,
            decay=40.0,
            fwidth=0.06,
        )
    if mode == "screen":
        return dict(
            cases=[
                "flat_Ti",
                "scaled_bare_slanted_groove",
                "scaled_inner_wall_film_250nm",
                "scaled_inner_wall_film_500nm",
            ],
            wavelengths=None,
            polarizations=["Ez", "Hz"],
            resolutions=[48],
            pml_values=[4.0],
            substrate=25.0,
            air=12.0,
            decay=60.0,
            fwidth=0.06,
        )
    if mode == "refine":
        return dict(
            cases=["scaled_bare_slanted_groove", "scaled_inner_wall_film_500nm"],
            wavelengths=None,
            polarizations=["Ez", "Hz"],
            resolutions=[64, 80],
            pml_values=[4.0],
            substrate=25.0,
            air=12.0,
            decay=60.0,
            fwidth=0.06,
        )
    raise ValueError(f"Unknown mode: {mode}")


def _nearest(values: np.ndarray, target: float) -> float:
    return float(values[np.argmin(np.abs(values - target))])


def _select_screen_wavelengths(nk_table: pd.DataFrame) -> list[float]:
    wl = nk_table["wavelength_um"].to_numpy(dtype=float)
    strict = wl[(wl >= BAND_STRICT_LO) & (wl <= BAND_STRICT_HI)]
    if len(strict) < 25:
        raise ValueError("Measured n,k table has fewer than 25 strict-band samples.")
    chosen = set(strict[np.linspace(0, len(strict) - 1, 25).round().astype(int)].astype(float))
    for target in [BAND_STRICT_LO, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, BAND_STRICT_HI]:
        chosen.add(_nearest(strict, target))
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


def _case_thickness(case_name: str) -> float:
    if case_name == "scaled_inner_wall_film_250nm":
        return 0.25
    if case_name == "scaled_inner_wall_film_500nm":
        return 0.50
    return 0.0


def _case_label(case_name: str) -> str:
    return {
        "flat_Ti": "flat Ti",
        "scaled_bare_slanted_groove": "scaled bare groove",
        "scaled_inner_wall_film_250nm": "scaled 250 nm film",
        "scaled_inner_wall_film_500nm": "scaled 500 nm film",
    }.get(case_name, case_name)


def _configure_args(args) -> None:
    defaults = _mode_defaults(args.mode)
    table = load_measured_nk_table(args.nk_csv)
    args.cases = _parse_str_list(args.cases) or defaults["cases"]
    args.wavelengths_um = _parse_float_list(args.wavelengths_um) or defaults["wavelengths"]
    if args.wavelengths_um is None:
        args.wavelengths_um = _select_refine_wavelengths(table) if args.mode == "refine" else _select_screen_wavelengths(table)
    args.polarizations = _parse_str_list(args.polarizations) or defaults["polarizations"]
    args.resolutions = [int(x) for x in (args.resolution or defaults["resolutions"])]
    args.pml_values_um = _parse_float_list(args.pml_thickness_um) or defaults["pml_values"]
    if args.mode == "refine" and 6.0 not in args.pml_values_um:
        args.pml_values_um.append(6.0)
    args.substrate_thickness_um = args.substrate_thickness_um or defaults["substrate"]
    args.air_buffer_um = args.air_buffer_um or defaults["air"]
    args.decay_db = args.decay_db if args.decay_db is not None else defaults["decay"]
    args.fwidth_fraction = args.fwidth_fraction if args.fwidth_fraction is not None else defaults["fwidth"]
    args.period_um = args.period_um or DEFAULT_PERIOD_UM
    args.top_width_um = args.top_width_um or DEFAULT_TOP_WIDTH_UM
    args.bottom_width_um = args.bottom_width_um or DEFAULT_BOTTOM_WIDTH_UM
    args.depth_um = args.depth_um or DEFAULT_DEPTH_UM
    args.tilt_angle_deg = args.tilt_angle_deg or DEFAULT_TILT_DEG
    args.wavelengths_um = sorted(set(float(x) for x in args.wavelengths_um))
    args.pml_values_um = sorted(set(float(x) for x in args.pml_values_um))
    args.strict_planned_point_count = sum(BAND_STRICT_LO <= wl <= BAND_STRICT_HI for wl in args.wavelengths_um)
    if args.depth_um >= args.substrate_thickness_um:
        raise ValueError("depth_um must be smaller than substrate_thickness_um.")


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


def _geometry_factory(case_name: str, film_medium, args):
    ti = get_ti_medium()
    air = mp.Medium(epsilon=1.0)
    thickness = _case_thickness(case_name)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        if case_name == "flat_Ti":
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
        if case_name == "scaled_bare_slanted_groove":
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
        if case_name.startswith("scaled_inner_wall_film"):
            if film_medium is None:
                raise ValueError("film_medium is required for coated groove cases.")
            return build_inner_wall_film_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.top_width_um,
                bottom_width_um=args.bottom_width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=args.tilt_angle_deg,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                film_thickness_um=thickness,
                medium_substrate=ti,
                medium_film=film_medium,
                medium_groove=air,
                coating_mode="sidewalls_and_bottom",
            )
        raise ValueError(f"Unknown case_name: {case_name}")

    return factory


def _transmission_flag(T: float) -> str:
    flag = str(opaque_substrate_transmission_check(np.array([T]))["transmission_quality_flag"][0])
    if flag == "NUMERICAL_PASS":
        return "NUMERICAL_PASS"
    if flag == "WARNING":
        return "WARNING_EXCLUDED_FROM_STRICT_METRIC"
    return "FAIL_EXCLUDED_FROM_ALL_METRICS"


def _run_one(case_name: str, pol: str, wl: float, resolution: int, pml_um: float, args, logger) -> dict:
    thickness = _case_thickness(case_name)
    film_medium = None
    film_meta = _film_meta_empty()
    if case_name.startswith("scaled_inner_wall_film"):
        film_medium, meta = get_measured_lossy_wall_film_medium_single_wavelength(
            wl, args.nk_csv, allow_extrapolation=False
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
        geometry_factory=_geometry_factory(case_name, film_medium, args),
        period_um=args.period_um,
        wavelength_um=wl,
        resolution=resolution,
        pml_thickness_um=pml_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        decay_db=args.decay_db,
        source_component=pol,
        fwidth_fraction=args.fwidth_fraction,
        solver_version="scaled_deep_inner_wall_film_single_wavelength_v1",
        source_mode="single_wavelength_narrowband",
        logger=logger,
    )
    ti_flag = "VALID" if wl <= BAND_STRICT_HI else "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING"
    trans_flag = _transmission_flag(float(result["T"]))
    finite = all(np.isfinite(result[k]) for k in ["R", "T", "A", "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"])
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
        "case_name": case_name,
        "film_thickness_um": thickness,
        "coating_mode": "none" if thickness == 0 else "sidewalls_and_bottom",
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
        "top_width_um": args.top_width_um,
        "bottom_width_um": args.bottom_width_um,
        "depth_um": args.depth_um,
        "tilt_angle_deg": args.tilt_angle_deg,
        "film_grid_points": thickness * resolution,
        "normalization_note": "D00 convention: R=refl/abs(input), T=-trans/abs(input), A=1-R-T.",
        "walltime_s": result["walltime_s"],
    }
    row.update(film_meta)
    return row


def _route_class(mean_a: float) -> str:
    if not np.isfinite(mean_a):
        return "NOT_QUANTITATIVE"
    if mean_a < 0.30:
        return "SCALED_INNER_WALL_FILM_ROUTE_FAIL_FOR_HIGH_EMISSIVITY"
    if mean_a < 0.60:
        return "PARTIAL_ABSORPTION_ENHANCEMENT_BUT_BELOW_CONTINUE_THRESHOLD"
    if mean_a < 0.80:
        return "ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION"
    return "HIGH_EMISSIVITY_CANDIDATE_FOR_DIRECTIONAL_FOLLOWUP"


def _compute_metrics(spectra: pd.DataFrame, args) -> pd.DataFrame:
    rows = []
    group_cols = ["case_name", "polarization", "film_thickness_um", "resolution", "pml_thickness_um"]
    for keys, df in spectra.groupby(group_cols, sort=False):
        case, pol, thickness, resolution, pml_um = keys
        strict_mask = df["included_in_strict_metric"].to_numpy(dtype=bool)
        strict_count = int(strict_mask.sum())
        enough = strict_count >= 0.9 * args.strict_planned_point_count
        mean_strict = wavelength_integrated_average(
            df["A"].to_numpy(), df["wavelength_um"].to_numpy(),
            BAND_STRICT_LO, BAND_STRICT_HI, valid_mask=strict_mask,
        )
        strict_status = wavelength_integrated_average.last_status
        if args.mode == "smoke":
            numerical_status = "CODE_PASS" if strict_count > 0 else "FAIL"
        elif not enough:
            mean_strict = np.nan
            numerical_status = "INCOMPLETE_NUMERICAL_COVERAGE"
        elif strict_status != "valid":
            numerical_status = strict_status.upper()
        else:
            numerical_status = "NUMERICAL_PASS" if args.mode == "refine" else "NUMERICAL_SCREENING"
        ext_mask = df["included_in_extended_metric"].to_numpy(dtype=bool)
        mean_extended = wavelength_integrated_average(
            df["A"].to_numpy(), df["wavelength_um"].to_numpy(),
            BAND_STRICT_LO, BAND_EXTENDED_HI, valid_mask=ext_mask,
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
            "case_name": case,
            "polarization": pol,
            "film_thickness_um": thickness,
            "resolution": resolution,
            "pml_thickness_um": pml_um,
            "mean_A_8p1014_12p398_strict": mean_strict,
            "mean_A_8p1014_12p962_extended": mean_extended,
            "peak_A_strict": peak_a,
            "peak_wavelength_strict_um": peak_wl,
            "valid_point_count_strict": strict_count,
            "total_planned_point_count_strict": args.strict_planned_point_count,
            "valid_point_count_extended": int(ext_mask.sum()),
            "numerical_status": numerical_status,
            "max_abs_T": float(np.nanmax(np.abs(df["T"].to_numpy()))),
            "route_decision": _route_class(mean_strict),
        })
    metrics = pd.DataFrame(rows)
    if metrics.empty:
        return metrics

    ref_cols = ["polarization", "resolution", "pml_thickness_um"]
    for idx, row in metrics.iterrows():
        same = metrics[
            (metrics["polarization"] == row["polarization"])
            & (metrics["resolution"] == row["resolution"])
            & (metrics["pml_thickness_um"] == row["pml_thickness_um"])
        ]
        flat = same[same["case_name"] == "flat_Ti"]
        bare = same[same["case_name"] == "scaled_bare_slanted_groove"]
        for ref_name, ref_df in [("flat_Ti", flat), ("scaled_bare", bare)]:
            col_abs = f"enhancement_over_{ref_name}_absolute"
            col_rel = f"enhancement_over_{ref_name}_relative"
            if ref_df.empty or not np.isfinite(row["mean_A_8p1014_12p398_strict"]):
                metrics.loc[idx, col_abs] = np.nan
                metrics.loc[idx, col_rel] = np.nan
                continue
            ref_a = float(ref_df["mean_A_8p1014_12p398_strict"].iloc[0])
            delta = float(row["mean_A_8p1014_12p398_strict"] - ref_a)
            metrics.loc[idx, col_abs] = delta
            metrics.loc[idx, col_rel] = delta / ref_a if ref_a else np.nan

    proxy_rows = []
    for (case, thickness, resolution, pml_um), df in spectra.groupby(
        ["case_name", "film_thickness_um", "resolution", "pml_thickness_um"], sort=False
    ):
        pivot = df.pivot_table(index="wavelength_um", columns="polarization", values=["A", "included_in_strict_metric"], aggfunc="first")
        if not {"Ez", "Hz"}.issubset(set(pivot["A"].columns)):
            continue
        wavelengths = pivot.index.to_numpy(dtype=float)
        proxy_a = 0.5 * (pivot[("A", "Ez")].to_numpy(dtype=float) + pivot[("A", "Hz")].to_numpy(dtype=float))
        valid = pivot[("included_in_strict_metric", "Ez")].to_numpy(dtype=bool) & pivot[("included_in_strict_metric", "Hz")].to_numpy(dtype=bool)
        mean_proxy = wavelength_integrated_average(proxy_a, wavelengths, BAND_STRICT_LO, BAND_STRICT_HI, valid_mask=valid)
        if args.mode != "smoke" and int(valid.sum()) < 0.9 * args.strict_planned_point_count:
            mean_proxy = np.nan
            status = "INCOMPLETE_NUMERICAL_COVERAGE"
        else:
            status = "CODE_PASS" if args.mode == "smoke" else wavelength_integrated_average.last_status.upper()
        proxy_rows.append({
            "case_name": case,
            "polarization": "unpolarized_2D_proxy",
            "film_thickness_um": thickness,
            "resolution": resolution,
            "pml_thickness_um": pml_um,
            "mean_A_8p1014_12p398_strict": mean_proxy,
            "peak_A_strict": float(np.nanmax(proxy_a)) if proxy_a.size else np.nan,
            "peak_wavelength_strict_um": float(wavelengths[np.nanargmax(proxy_a)]) if proxy_a.size else np.nan,
            "valid_point_count_strict": int(valid.sum()),
            "total_planned_point_count_strict": args.strict_planned_point_count,
            "numerical_status": status,
            "route_decision": _route_class(mean_proxy),
            "physics_note": "A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2; not a true 3D sample emissivity.",
        })
    if proxy_rows:
        metrics = pd.concat([metrics, pd.DataFrame(proxy_rows)], ignore_index=True)
    return metrics


def _resolution_check(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (case, pol, pml_um), df in metrics.groupby(["case_name", "polarization", "pml_thickness_um"], sort=False):
        if not {64, 80}.issubset(set(df["resolution"].astype(int))):
            continue
        a64 = float(df[df["resolution"] == 64]["mean_A_8p1014_12p398_strict"].iloc[0])
        a80 = float(df[df["resolution"] == 80]["mean_A_8p1014_12p398_strict"].iloc[0])
        rows.append({
            "case_name": case,
            "polarization": pol,
            "pml_thickness_um": pml_um,
            "mean_A_res64": a64,
            "mean_A_res80": a80,
            "abs_delta_mean_A": abs(a80 - a64),
            "resolution_check_status": "NUMERICAL_PASS" if abs(a80 - a64) < 0.01 else "WARNING",
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame([{"resolution_check_status": "NOT_RUN"}])


def _pml_check(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (case, pol, res), df in metrics.groupby(["case_name", "polarization", "resolution"], sort=False):
        if not {4.0, 6.0}.issubset(set(df["pml_thickness_um"].astype(float))):
            continue
        a4 = float(df[df["pml_thickness_um"] == 4.0]["mean_A_8p1014_12p398_strict"].iloc[0])
        a6 = float(df[df["pml_thickness_um"] == 6.0]["mean_A_8p1014_12p398_strict"].iloc[0])
        rows.append({
            "case_name": case,
            "polarization": pol,
            "resolution": res,
            "mean_A_pml4": a4,
            "mean_A_pml6": a6,
            "abs_delta_mean_A": abs(a6 - a4),
            "pml_check_status": "NUMERICAL_PASS" if abs(a6 - a4) < 0.01 else "WARNING",
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame([{"pml_check_status": "NOT_RUN"}])


def _plot_geometry(thickness: float, out_path: Path, args) -> None:
    y_surface = 0.0
    outer = slanted_groove_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=y_surface,
    )
    half = 0.5 * args.period_um
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.add_patch(plt.Rectangle((-half, -args.substrate_thickness_um), args.period_um, args.substrate_thickness_um, facecolor="#808080", edgecolor="black", label="Ti substrate"))
    if thickness <= 0:
        ax.fill(*zip(*(outer + [outer[0]])), facecolor="white", edgecolor="#1F77B4", lw=1.4, label="air groove")
    else:
        ax.fill(*zip(*(outer + [outer[0]])), facecolor="#D62728", alpha=0.75, edgecolor="black", lw=1.1, label="measured_lossy_wall_film")
        inner = slanted_groove_vertices(
            top_width_um=args.top_width_um - 2.0 * thickness,
            bottom_width_um=args.bottom_width_um - 2.0 * thickness,
            depth_um=args.depth_um - thickness,
            tilt_angle_deg=args.tilt_angle_deg,
            y_surface=y_surface,
        )
        ax.fill(*zip(*(inner + [inner[0]])), facecolor="white", edgecolor="#1F77B4", lw=1.4, label="air core")
    ax.axhline(y_surface, color="black", lw=1.0)
    ax.text(-half + 0.5, y_surface + 0.35, "top Ti land remains uncoated", fontsize=8)
    ax.set_xlim(-half, half)
    ax.set_ylim(-args.depth_um - 1.0, y_surface + 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"D11 scaled groove geometry, inner film t={thickness:g} um")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.2)
    save_figure(fig, out_path)
    plt.close(fig)


def _make_plots(spectra: pd.DataFrame, metrics: pd.DataFrame, paths: dict[str, Path], args) -> None:
    _generate_geometry_plots(paths, args)
    for pol, out in [("Ez", paths["ez"]), ("Hz", paths["hz"])]:
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        df_pol = spectra[spectra["polarization"] == pol]
        if df_pol.empty:
            ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
            ax.set_axis_off()
        else:
            for case, df in df_pol.groupby("case_name", sort=False):
                ax.plot(df["wavelength_um"], df["A"], marker="o", lw=1.2, label=_case_label(case))
            ax.axvspan(BAND_STRICT_LO, BAND_STRICT_HI, color="#E8EEF7", alpha=0.45)
            ax.set_xlabel("Wavelength (um)")
            ax.set_ylabel("A")
            ax.set_title(f"D11 {pol} absorptance comparison")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7)
        save_figure(fig, out)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    pivot = spectra.pivot_table(index=["case_name", "wavelength_um"], columns="polarization", values="A", aggfunc="first").reset_index()
    if {"Ez", "Hz"}.issubset(set(pivot.columns)):
        pivot["A_unpolarized_2D_proxy"] = 0.5 * (pivot["Ez"] + pivot["Hz"])
        for case, df in pivot.groupby("case_name", sort=False):
            ax.plot(df["wavelength_um"], df["A_unpolarized_2D_proxy"], marker="o", lw=1.2, label=_case_label(case))
        ax.axvspan(BAND_STRICT_LO, BAND_STRICT_HI, color="#E8EEF7", alpha=0.45)
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("A proxy")
        ax.set_title("D11 unpolarized 2D proxy")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "No paired Ez/Hz data", ha="center", va="center")
        ax.set_axis_off()
    save_figure(fig, paths["proxy"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    main = metrics[metrics["pml_thickness_um"].fillna(4.0).eq(4.0)]
    if main.empty:
        ax.text(0.5, 0.5, "No metrics", ha="center", va="center")
        ax.set_axis_off()
    else:
        labels = [f"{_case_label(r.case_name)}\n{r.polarization}" for r in main.itertuples()]
        ax.bar(labels, main["mean_A_8p1014_12p398_strict"])
        ax.axhline(0.60, color="red", ls="--", lw=1, label="0.60 continue")
        ax.axhline(0.80, color="purple", ls=":", lw=1, label="0.80 target")
        ax.set_ylabel("Mean A strict")
        ax.set_title("D11 route decision")
        ax.tick_params(axis="x", labelrotation=45)
        ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=8)
    save_figure(fig, paths["mean"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    enh = main[main["case_name"].str.contains("film", na=False)]
    if "enhancement_over_scaled_bare_absolute" in enh.columns and not enh.empty:
        labels = [f"{_case_label(r.case_name)}\n{r.polarization}" for r in enh.itertuples()]
        ax.bar(labels, enh["enhancement_over_scaled_bare_absolute"])
        ax.set_ylabel("Delta mean A vs scaled bare")
        ax.set_title("D11 enhancement over scaled bare groove")
        ax.tick_params(axis="x", labelrotation=45)
        ax.grid(axis="y", alpha=0.25)
    else:
        ax.text(0.5, 0.5, "No enhancement metrics", ha="center", va="center")
        ax.set_axis_off()
    save_figure(fig, paths["enhancement"])
    plt.close(fig)


def _generate_geometry_plots(paths: dict[str, Path], args) -> None:
    _plot_geometry(0.0, paths["geom_bare"], args)
    _plot_geometry(0.25, paths["geom_250"], args)
    _plot_geometry(0.50, paths["geom_500"], args)


def _md_table(metrics: pd.DataFrame) -> str:
    cols = [
        "case_name", "polarization", "mean_A_8p1014_12p398_strict",
        "enhancement_over_scaled_bare_absolute", "enhancement_over_flat_Ti_absolute",
        "route_decision", "numerical_status",
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


def _overall_level(metrics: pd.DataFrame, args) -> str:
    if metrics.empty:
        return "FAIL"
    if args.mode == "smoke":
        return "CODE_PASS"
    if args.mode == "screen":
        return "NUMERICAL_SCREENING"
    return "NUMERICAL_PASS" if (metrics["numerical_status"] == "NUMERICAL_PASS").all() else "WARNING"


def _write_report(paths: dict[str, Path], spectra: pd.DataFrame, metrics: pd.DataFrame, res_check: pd.DataFrame, pml_check: pd.DataFrame, args) -> None:
    ensure_dir(paths["report"].parent)
    level = _overall_level(metrics, args)
    main = metrics[
        (metrics["case_name"] == "scaled_inner_wall_film_500nm")
        & metrics["mean_A_8p1014_12p398_strict"].notna()
    ].copy()
    if main.empty:
        best_text = "No strict quantitative 500 nm metric is available."
        pass60 = False
        pass80 = False
    else:
        best = main.sort_values("mean_A_8p1014_12p398_strict", ascending=False).iloc[0]
        pass60 = bool(best["mean_A_8p1014_12p398_strict"] >= 0.60)
        pass80 = bool(best["mean_A_8p1014_12p398_strict"] >= 0.80)
        best_text = (
            f"Best 500 nm strict mean A = {best['mean_A_8p1014_12p398_strict']:.4g} "
            f"({best['polarization']}, res={int(best['resolution'])}, PML={best['pml_thickness_um']:g} um)."
        )
    if pass80:
        decision = "HIGH_EMISSIVITY_CANDIDATE_FOR_DIRECTIONAL_FOLLOWUP; next use mode decomposition because P=50 um supports many diffraction orders."
    elif pass60:
        decision = "ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION; refine resolution/PML and then optimize geometry/material coverage."
    else:
        decision = "Below the 0.60 continue threshold; do not automatically enter angular directionality optimization."
    report = f"""# D11 Scaled Deep Slanted Groove Inner-Wall Film Validation

Overall result level: **{level}**

## Required Answers

1. Geometry built?
   Yes: P={args.period_um:g} um, top/bottom width={args.top_width_um:g}/{args.bottom_width_um:g} um, depth={args.depth_um:g} um, tilt={args.tilt_angle_deg:g} deg.

2. Is the 500 nm film only on the inner walls and bottom?
   The geometry function uses an outer film prism plus an inner air prism opening at the top surface.  The generated geometry plots mark the top Ti land as uncoated.

3. 500 nm result:
   {best_text}

4. Diffraction-order caveat:
   This script uses total flux monitors over the full period, so A is suitable for total absorption screening.  It does not diagnose directional emission.  If A reaches the continue threshold, use mode decomposition next.

5. Route decision:
   {decision}

6. Model scope:
   This is a 2D equivalent slanted-groove model with `measured_lossy_wall_film`; it is not a true 3D laser-drilled oblique-hole array and the film is not renamed TiO2.

## Metrics Preview

{_md_table(metrics)}

## Numerical Checks

Resolution check rows: {len(res_check)}

PML check rows: {len(pml_check)}

## Outputs

- Spectra: `{paths['spectra']}`
- Metrics: `{paths['metrics']}`
- Resolution check: `{paths['resolution']}`
- PML check: `{paths['pml']}`
- Log: `{paths['log']}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def _write_placeholder_field_figures(paths: dict[str, Path], args, reason: str) -> None:
    for key, title in [
        ("best_e2", "D11 best-case |E|^2 not generated"),
        ("best_h2", "D11 best-case |H|^2 not generated"),
        ("best_overlay", "D11 best-case geometry overlay not generated"),
    ]:
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        ax.text(0.5, 0.5, reason, ha="center", va="center", wrap=True)
        ax.set_title(title)
        ax.set_axis_off()
        save_figure(fig, paths[key])
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["smoke", "screen", "refine"], default="smoke")
    parser.add_argument("--nk-csv", type=Path, required=True)
    parser.add_argument("--cases", nargs="*")
    parser.add_argument("--wavelengths-um", nargs="*")
    parser.add_argument("--polarizations", nargs="*")
    parser.add_argument("--resolution", nargs="*", type=int)
    parser.add_argument("--pml-thickness-um", nargs="*")
    parser.add_argument("--substrate-thickness-um", type=float)
    parser.add_argument("--air-buffer-um", type=float)
    parser.add_argument("--decay-db", type=float)
    parser.add_argument("--fwidth-fraction", type=float)
    parser.add_argument("--period-um", type=float)
    parser.add_argument("--top-width-um", type=float)
    parser.add_argument("--bottom-width-um", type=float)
    parser.add_argument("--depth-um", type=float)
    parser.add_argument("--tilt-angle-deg", type=float)
    args = parser.parse_args()
    _configure_args(args)

    paths = _paths()
    logger = _setup_logger(paths["log"])
    logger.info("Starting D11 mode=%s", args.mode)
    _generate_geometry_plots(paths, args)

    rows = []
    for resolution in args.resolutions:
        for pml_um in args.pml_values_um:
            for case_name in args.cases:
                for pol in args.polarizations:
                    for wl in args.wavelengths_um:
                        logger.info(
                            ">>> case=%s pol=%s wl=%g res=%s pml=%g",
                            case_name, pol, wl, resolution, pml_um,
                        )
                        rows.append(_run_one(case_name, pol, wl, resolution, pml_um, args, logger))

    spectra = pd.DataFrame(rows)
    metrics = _compute_metrics(spectra, args)
    res_check = _resolution_check(metrics)
    pml_check = _pml_check(metrics)

    for path in [paths["spectra"], paths["metrics"], paths["resolution"], paths["pml"]]:
        ensure_dir(path.parent)
    spectra.to_csv(paths["spectra"], index=False)
    metrics.to_csv(paths["metrics"], index=False)
    res_check.to_csv(paths["resolution"], index=False)
    pml_check.to_csv(paths["pml"], index=False)
    _make_plots(spectra, metrics, paths, args)
    _write_placeholder_field_figures(
        paths,
        args,
        "Field maps are only generated after a screen candidate reaches mean_A_strict >= 0.30; this run did not execute the field-map branch.",
    )
    _write_report(paths, spectra, metrics, res_check, pml_check, args)
    logger.info("Wrote %s", paths["spectra"])
    logger.info("Wrote %s", paths["metrics"])
    logger.info("Wrote %s", paths["report"])


if __name__ == "__main__":
    main()
