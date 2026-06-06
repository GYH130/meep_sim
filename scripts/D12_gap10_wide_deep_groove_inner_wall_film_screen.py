"""D12 gap-10 wide/deep straight-groove inner-wall film screen.

This diagnostics_v2 script tests whether widening the P=50 um groove opening
from 20 um to 40 um, reducing the top Ti gap to 10 um, and deepening the groove
to 30 um lets ``measured_lossy_wall_film`` on only the inner sidewalls and
bottom improve total absorptance in the strict Ti-valid mid-infrared band.

Important scope notes:
- One narrowband Meep run is performed per wavelength.
- Total flux monitors span the full period, so R/T/A are total powers summed
  over all propagating diffraction orders.
- Directionality is not evaluated here; it needs a later mode-decomposition
  analysis if total absorption reaches the continuation threshold.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import logging
import multiprocessing
import os
import re
import sys
import time
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib-codex-cache"),
)

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
    interpolate_nk_at_wavelength,
    load_measured_nk_table,
)
from src.postprocess import opaque_substrate_transmission_check, wavelength_integrated_average
from src.simulation import run_periodic_2d_metal_single_wavelength


BAND_STRICT_LO = 8.1014
BAND_STRICT_HI = TI_RAKIC_VALID_LAMBDA_UM[1]
BAND_EXTENDED_HI = 12.962
CONTINUE_THRESHOLD = 0.60
HIGH_EMISSIVITY_TARGET = 0.80
FIELD_MAP_THRESHOLD = 0.30

DEFAULT_PERIOD_UM = 50.0
DEFAULT_TOP_WIDTH_UM = 40.0
DEFAULT_BOTTOM_WIDTH_UM = 40.0
DEFAULT_DEPTH_UM = 30.0
DEFAULT_TILT_DEG = 0.0
DEFAULT_SUBSTRATE_THICKNESS_UM = 45.0
DEFAULT_AIR_BUFFER_UM = 15.0
DEFAULT_PML_UM = 4.0
TARGET_GAP_UM = 10.0

FLAT = "flat_Ti"
BARE = "gap10_bare_wide_deep_straight_groove"
FILM_500 = "gap10_inner_wall_film_500nm"
FILM_1UM = "gap10_inner_wall_film_1um"
FILM_2UM = "gap10_inner_wall_film_2um"
PROXY = "unpolarized_2D_proxy"
ALL_CASES = [FLAT, BARE, FILM_500, FILM_1UM, FILM_2UM]
DYNAMIC_FILM_PREFIX = "gap10_inner_wall_film_"
DYNAMIC_FILM_SUFFIX = "um"


def _paths(output_tag: str | None = None) -> dict[str, Path]:
    fig_dir = project_path("results", "diagnostics_v2", "figures")
    if output_tag:
        tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", output_tag).strip("_")
        if not tag:
            raise ValueError("--output-tag must contain at least one filename-safe character.")
        return {
            "spectra": project_path("results", "diagnostics_v2", "tables", f"{tag}_spectra.csv"),
            "spectra_checkpoint": project_path("results", "diagnostics_v2", "tables", f"{tag}_spectra_checkpoint.csv"),
            "metrics": project_path("results", "diagnostics_v2", "tables", f"{tag}_metrics.csv"),
            "resolution": project_path("results", "diagnostics_v2", "tables", f"{tag}_resolution_check.csv"),
            "pml": project_path("results", "diagnostics_v2", "tables", f"{tag}_pml_check.csv"),
            "report": project_path("results", "diagnostics_v2", "reports", f"{tag}_report.md"),
            "log": project_path("logs", "diagnostics_v2", f"{tag}.log"),
            "geom_bare": fig_dir / f"{tag}_geometry_gap10_bare.png",
            "geom_500": fig_dir / f"{tag}_geometry_gap10_film_500nm.png",
            "geom_1um": fig_dir / f"{tag}_geometry_gap10_film_1um.png",
            "geom_2um": fig_dir / f"{tag}_geometry_gap10_film_2um.png",
            "ez": fig_dir / f"{tag}_Ez_absorptance_comparison.png",
            "hz": fig_dir / f"{tag}_Hz_absorptance_comparison.png",
            "proxy": fig_dir / f"{tag}_unpolarized_proxy_comparison.png",
            "mean": fig_dir / f"{tag}_mean_A_route_decision.png",
            "enhancement": fig_dir / f"{tag}_enhancement_over_bare_gap10.png",
            "gap10_vs_gap30": fig_dir / f"{tag}_gap10_vs_gap30_reference.png",
            "best_e2": fig_dir / f"{tag}_best_case_E2.png",
            "best_h2": fig_dir / f"{tag}_best_case_H2.png",
            "best_overlay": fig_dir / f"{tag}_best_case_geometry_overlay.png",
        }
    return {
        "spectra": project_path(
            "results",
            "diagnostics_v2",
            "tables",
            "D12_gap10_wide_deep_groove_inner_wall_film_spectra.csv",
        ),
        "spectra_checkpoint": project_path(
            "results",
            "diagnostics_v2",
            "tables",
            "D12_gap10_wide_deep_groove_inner_wall_film_spectra_checkpoint.csv",
        ),
        "metrics": project_path(
            "results",
            "diagnostics_v2",
            "tables",
            "D12_gap10_wide_deep_groove_inner_wall_film_metrics.csv",
        ),
        "resolution": project_path(
            "results",
            "diagnostics_v2",
            "tables",
            "D12_gap10_wide_deep_groove_resolution_check.csv",
        ),
        "pml": project_path(
            "results",
            "diagnostics_v2",
            "tables",
            "D12_gap10_wide_deep_groove_pml_check.csv",
        ),
        "report": project_path(
            "results",
            "diagnostics_v2",
            "reports",
            "D12_gap10_wide_deep_groove_inner_wall_film_report.md",
        ),
        "log": project_path(
            "logs",
            "diagnostics_v2",
            "D12_gap10_wide_deep_groove_inner_wall_film.log",
        ),
        "geom_bare": fig_dir / "D12_geometry_gap10_bare.png",
        "geom_500": fig_dir / "D12_geometry_gap10_film_500nm.png",
        "geom_1um": fig_dir / "D12_geometry_gap10_film_1um.png",
        "geom_2um": fig_dir / "D12_geometry_gap10_film_2um.png",
        "ez": fig_dir / "D12_Ez_absorptance_comparison.png",
        "hz": fig_dir / "D12_Hz_absorptance_comparison.png",
        "proxy": fig_dir / "D12_unpolarized_proxy_comparison.png",
        "mean": fig_dir / "D12_mean_A_route_decision.png",
        "enhancement": fig_dir / "D12_enhancement_over_bare_gap10.png",
        "gap10_vs_gap30": fig_dir / "D12_gap10_vs_gap30_reference.png",
        "best_e2": fig_dir / "D12_best_case_E2.png",
        "best_h2": fig_dir / "D12_best_case_H2.png",
        "best_overlay": fig_dir / "D12_best_case_geometry_overlay.png",
    }


def _setup_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger("D12_gap10_wide_deep_groove")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    if not _mpi_am_master():
        logger.addHandler(logging.NullHandler())
        return logger
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


def _mpi_am_master() -> bool:
    try:
        import meep as mp

        return bool(mp.am_master())
    except Exception:
        return True


def _setup_worker_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger(f"D12_gap10_wide_deep_groove.worker.{os.getpid()}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler(path, encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
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
        return {
            "cases": [BARE, FILM_1UM],
            "wavelengths": [9.0, 10.0, 12.0],
            "polarizations": ["Ez", "Hz"],
            "resolutions": [24],
            "pml_values": [DEFAULT_PML_UM],
            "substrate": DEFAULT_SUBSTRATE_THICKNESS_UM,
            "air": DEFAULT_AIR_BUFFER_UM,
            "decay": 40.0,
            "fwidth": 0.06,
        }
    if mode == "screen":
        return {
            "cases": ALL_CASES,
            "wavelengths": None,
            "polarizations": ["Ez", "Hz"],
            "resolutions": [32],
            "pml_values": [DEFAULT_PML_UM],
            "substrate": DEFAULT_SUBSTRATE_THICKNESS_UM,
            "air": DEFAULT_AIR_BUFFER_UM,
            "decay": 60.0,
            "fwidth": 0.06,
        }
    if mode == "refine":
        return {
            "cases": [BARE, FILM_500, FILM_1UM, FILM_2UM],
            "wavelengths": None,
            "polarizations": ["Ez", "Hz"],
            "resolutions": [48, 64],
            "pml_values": [DEFAULT_PML_UM],
            "substrate": DEFAULT_SUBSTRATE_THICKNESS_UM,
            "air": DEFAULT_AIR_BUFFER_UM,
            "decay": 60.0,
            "fwidth": 0.06,
        }
    raise ValueError(f"Unknown mode: {mode}")


def _nearest(values: np.ndarray, target: float) -> float:
    return float(values[np.argmin(np.abs(values - target))])


def _strict_raw_wavelengths(nk_table: pd.DataFrame) -> np.ndarray:
    wl = nk_table["wavelength_um"].to_numpy(dtype=float)
    return wl[(wl >= BAND_STRICT_LO) & (wl <= BAND_STRICT_HI)]


def _select_screen_wavelengths(nk_table: pd.DataFrame) -> list[float]:
    strict = _strict_raw_wavelengths(nk_table)
    if len(strict) < 25:
        raise ValueError("Measured n,k table has fewer than 25 strict-band samples.")
    chosen = set(strict[np.linspace(0, len(strict) - 1, 25).round().astype(int)].astype(float))
    for target in [
        BAND_STRICT_LO,
        8.5,
        9.0,
        9.5,
        10.0,
        10.5,
        11.0,
        11.5,
        12.0,
        BAND_STRICT_HI,
    ]:
        chosen.add(_nearest(strict, target))
    return sorted(chosen)


def _select_refine_wavelengths(nk_table: pd.DataFrame, max_points: int = 61) -> list[float]:
    strict = _strict_raw_wavelengths(nk_table)
    if len(strict) <= max_points:
        return strict.astype(float).tolist()
    idx = np.linspace(0, len(strict) - 1, max_points).round().astype(int)
    return strict[idx].astype(float).tolist()


def _dynamic_film_case_name(thickness_um: float) -> str:
    token = f"{thickness_um:g}".replace(".", "p")
    return f"{DYNAMIC_FILM_PREFIX}{token}{DYNAMIC_FILM_SUFFIX}"


def _dynamic_case_thickness(case_name: str) -> float | None:
    if not case_name.startswith(DYNAMIC_FILM_PREFIX) or not case_name.endswith(DYNAMIC_FILM_SUFFIX):
        return None
    token = case_name[len(DYNAMIC_FILM_PREFIX):-len(DYNAMIC_FILM_SUFFIX)]
    if not token:
        return None
    try:
        thickness = float(token.replace("p", "."))
    except ValueError:
        return None
    return thickness if thickness > 0 else None


def _is_inner_wall_film_case(case_name: str) -> bool:
    return case_name in {FILM_500, FILM_1UM, FILM_2UM} or _dynamic_case_thickness(case_name) is not None


def _case_thickness(case_name: str) -> float:
    if case_name == FILM_500:
        return 0.50
    if case_name == FILM_1UM:
        return 1.0
    if case_name == FILM_2UM:
        return 2.0
    dynamic = _dynamic_case_thickness(case_name)
    if dynamic is not None:
        return dynamic
    return 0.0


def _case_label(case_name: str) -> str:
    dynamic = _dynamic_case_thickness(case_name)
    if dynamic is not None:
        return f"scaled {dynamic:g} um film"
    return {
        FLAT: "flat Ti",
        BARE: "gap10 bare wide/deep straight groove",
        FILM_500: "gap10 500 nm film",
        FILM_1UM: "gap10 1 um film",
        FILM_2UM: "gap10 2 um film",
    }.get(case_name, case_name)


def _configure_args(args) -> None:
    defaults = _mode_defaults(args.mode)
    table = load_measured_nk_table(args.nk_csv)
    requested_cases = _parse_str_list(args.cases)
    args.film_thicknesses_um = _parse_float_list(args.film_thicknesses_um)
    if args.film_thicknesses_um:
        if any(t <= 0 for t in args.film_thicknesses_um):
            raise ValueError("--film-thicknesses-um values must be positive.")
        dynamic_cases = [_dynamic_film_case_name(t) for t in args.film_thicknesses_um]
        args.dynamic_film_thickness_by_case = dict(zip(dynamic_cases, args.film_thicknesses_um))
        args.cases = requested_cases or dynamic_cases
    else:
        args.dynamic_film_thickness_by_case = {}
        args.cases = requested_cases or defaults["cases"]
    unknown = sorted(
        case for case in set(args.cases)
        if case not in ALL_CASES and _dynamic_case_thickness(case) is None
    )
    if unknown:
        raise ValueError(
            "Unknown D12 case(s): "
            f"{unknown}. Valid fixed cases: {ALL_CASES}; dynamic film cases use "
            f"{DYNAMIC_FILM_PREFIX}<thickness>{DYNAMIC_FILM_SUFFIX}, e.g. gap10_inner_wall_film_1um."
        )

    args.wavelengths_um = _parse_float_list(args.wavelengths_um) or defaults["wavelengths"]
    if args.wavelengths_um is None:
        args.wavelengths_um = (
            _select_refine_wavelengths(table)
            if args.mode == "refine"
            else _select_screen_wavelengths(table)
        )
    args.polarizations = _parse_str_list(args.polarizations) or defaults["polarizations"]
    for pol in args.polarizations:
        if pol not in {"Ez", "Hz"}:
            raise ValueError(f"polarization must be Ez or Hz, got {pol!r}")

    args.resolutions = [int(x) for x in (args.resolution or defaults["resolutions"])]
    args.pml_values_um = _parse_float_list(args.pml_thickness_um) or defaults["pml_values"]
    if args.mode == "refine" and 6.0 not in args.pml_values_um:
        args.pml_values_um.append(6.0)

    args.substrate_thickness_um = args.substrate_thickness_um or defaults["substrate"]
    args.air_buffer_um = args.air_buffer_um or defaults["air"]
    args.decay_db = args.decay_db if args.decay_db is not None else defaults["decay"]
    args.fwidth_fraction = args.fwidth_fraction if args.fwidth_fraction is not None else defaults["fwidth"]
    args.courant = args.courant if args.courant is not None else 0.5
    if not (0 < args.courant <= 0.5):
        raise ValueError("--courant must be in the interval (0, 0.5].")
    args.period_um = args.period_um or DEFAULT_PERIOD_UM
    args.top_width_um = args.top_width_um or DEFAULT_TOP_WIDTH_UM
    args.bottom_width_um = args.bottom_width_um or DEFAULT_BOTTOM_WIDTH_UM
    args.depth_um = args.depth_um or DEFAULT_DEPTH_UM
    args.tilt_angle_deg = args.tilt_angle_deg or DEFAULT_TILT_DEG
    args.wavelengths_um = sorted(set(float(x) for x in args.wavelengths_um))
    args.pml_values_um = sorted(set(float(x) for x in args.pml_values_um))
    args.strict_planned_point_count = sum(
        BAND_STRICT_LO <= wl <= BAND_STRICT_HI for wl in args.wavelengths_um
    )
    if args.depth_um >= args.substrate_thickness_um:
        raise ValueError("depth_um must be smaller than substrate_thickness_um.")
    gap_um = args.period_um - args.top_width_um
    if not np.isclose(args.top_width_um + gap_um, args.period_um):
        raise ValueError("top_width_um + gap_um must equal period_um.")
    for case in args.cases:
        thickness = _case_thickness(case)
        if thickness <= 0:
            continue
        if thickness >= 0.5 * args.top_width_um:
            raise ValueError(f"{case}: film thickness collapses the top opening.")
        if thickness >= 0.5 * args.bottom_width_um:
            raise ValueError(f"{case}: film thickness collapses the bottom opening.")
        if thickness >= args.depth_um:
            raise ValueError(f"{case}: film thickness must be smaller than groove depth.")


def _vertices_inside_cell(vertices: list[tuple[float, float]], period_um: float) -> bool:
    half = 0.5 * period_um
    return all(abs(x) <= half + 1e-9 for x, _ in vertices)


def _geometry_checks(args) -> pd.DataFrame:
    rows: list[dict] = []
    outer = slanted_groove_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=0.0,
    )
    rows.append({
        "check_name": "depth_smaller_than_substrate",
        "status": "PASS" if args.depth_um < args.substrate_thickness_um else "FAIL",
        "details": f"depth={args.depth_um:g}, substrate={args.substrate_thickness_um:g}",
    })
    rows.append({
        "check_name": "outer_groove_inside_unit_cell",
        "status": "PASS" if _vertices_inside_cell(outer, args.period_um) else "FAIL",
        "details": f"period={args.period_um:g}, vertices={outer}",
    })
    bottom_shift = args.depth_um * np.tan(np.deg2rad(args.tilt_angle_deg))
    bottom_right = bottom_shift + 0.5 * args.bottom_width_um
    bottom_left = bottom_shift - 0.5 * args.bottom_width_um
    rows.append({
        "check_name": "tilted_bottom_vertices_inside_unit_cell",
        "status": "PASS" if abs(bottom_left) <= 0.5 * args.period_um + 1e-9 and abs(bottom_right) <= 0.5 * args.period_um + 1e-9 else "FAIL",
        "details": (
            f"bottom_center_shift={bottom_shift:.6g} um, bottom_left={bottom_left:.6g} um, "
            f"bottom_right={bottom_right:.6g} um, allowed=[{-0.5 * args.period_um:.6g}, {0.5 * args.period_um:.6g}]"
        ),
    })
    film_check_thicknesses = sorted(
        {0.0, 0.50, 1.0, 2.0}
        | {float(_case_thickness(case)) for case in getattr(args, "cases", []) if _case_thickness(case) > 0}
    )
    for thickness in film_check_thicknesses:
        if thickness == 0:
            rows.append({
                "check_name": "zero_film_degenerates_to_bare_groove",
                "status": "PASS",
                "details": "D12 uses build_slanted_groove_geometry for bare case and the geometry API degenerates film_thickness_um=0 to that same bare geometry.",
            })
            continue
        inner_top = args.top_width_um - 2.0 * thickness
        inner_bottom = args.bottom_width_um - 2.0 * thickness
        inner_depth = args.depth_um - thickness
        inner = slanted_groove_vertices(
            top_width_um=inner_top,
            bottom_width_um=inner_bottom,
            depth_um=inner_depth,
            tilt_angle_deg=args.tilt_angle_deg,
            y_surface=0.0,
        )
        rows.append({
            "check_name": f"air_core_not_collapsed_{int(round(thickness * 1000))}nm",
            "status": "PASS" if inner_top > 0 and inner_bottom > 0 and inner_depth > 0 else "FAIL",
            "details": f"inner_top={inner_top:g}, inner_bottom={inner_bottom:g}, inner_depth={inner_depth:g}",
        })
        rows.append({
            "check_name": f"inner_air_core_inside_unit_cell_{int(round(thickness * 1000))}nm",
            "status": "PASS" if _vertices_inside_cell(inner, args.period_um) else "FAIL",
            "details": f"period={args.period_um:g}, vertices={inner}",
        })
    rows.append({
        "check_name": "gap_is_10um",
        "status": "PASS" if np.isclose(args.period_um - args.top_width_um, TARGET_GAP_UM) else "FAIL",
        "details": f"gap={args.period_um - args.top_width_um:g} um",
    })
    rows.append({
        "check_name": "period_original_D12_straight_slot",
        "status": "PASS",
        "details": (
            f"period={args.period_um:g} um; straight-slot D12 keeps the physical array pitch at P=50 um."
        ),
    })
    rows.append({
        "check_name": "top_width_plus_gap_equals_period",
        "status": "PASS" if np.isclose(args.top_width_um + (args.period_um - args.top_width_um), args.period_um) else "FAIL",
        "details": f"top_width={args.top_width_um:g}, gap={args.period_um - args.top_width_um:g}, period={args.period_um:g}",
    })
    rows.append({
        "check_name": "top_flat_ti_surface_uncoated",
        "status": "PASS",
        "details": "Only an outer groove-cavity film prism and an inner air prism are used; no top-land film blocks are added.",
    })
    return pd.DataFrame(rows)


def _raise_on_failed_geometry_checks(checks: pd.DataFrame) -> None:
    failed = checks[checks["status"] != "PASS"]
    if not failed.empty:
        raise ValueError("D12 geometry checks failed:\n" + failed.to_string(index=False))


def _apply_refine_gate(args, paths: dict[str, Path], logger: logging.Logger) -> None:
    if args.mode != "refine":
        return
    if not paths["metrics"].exists():
        raise RuntimeError(
            "Refine mode requires an existing D12 screen metrics CSV. "
            "Run --mode screen first and only refine if any film case reaches mean_A_strict >= 0.60."
        )
    screen_metrics = pd.read_csv(paths["metrics"])
    required_cols = {"case_name", "polarization", "mean_A_8p1014_12p398_strict"}
    missing = sorted(required_cols - set(screen_metrics.columns))
    if missing:
        raise RuntimeError(f"Existing metrics CSV is missing required refine-gate columns: {missing}")
    screen_films = screen_metrics[
        (screen_metrics["case_name"].isin([FILM_500, FILM_1UM, FILM_2UM]))
        & (screen_metrics["polarization"].isin(["Ez", "Hz", PROXY]))
    ].copy()
    best = float(screen_films["mean_A_8p1014_12p398_strict"].max()) if not screen_films.empty else np.nan
    if not np.isfinite(best) or best < CONTINUE_THRESHOLD:
        raise RuntimeError(
            "Refine is blocked because screen did not show "
            f"any D12 film case reaching mean_A_strict >= {CONTINUE_THRESHOLD:g}; "
            f"best available value is {best}."
        )
    ranked = (
        screen_films.groupby("case_name")["mean_A_8p1014_12p398_strict"]
        .max()
        .sort_values(ascending=False)
    )
    keep = [BARE] + ranked.head(2).index.tolist()
    args.cases = [case for case in keep if case in ALL_CASES]
    logger.info("D12 refine cases selected from screen metrics: %s", args.cases)


def _film_meta_empty() -> dict:
    return {
        "material_name_film": "NOT_APPLICABLE",
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


def _film_medium_for_wavelength(wl: float, args):
    film_medium, meta = get_measured_lossy_wall_film_medium_single_wavelength(
        wl,
        args.nk_csv,
        allow_extrapolation=False,
    )
    return film_medium, meta


def _geometry_factory(case_name: str, film_medium, args):
    import meep as mp

    ti = get_ti_medium()
    air = mp.Medium(epsilon=1.0)
    thickness = _case_thickness(case_name)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        if case_name == FLAT:
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
        if case_name == BARE:
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
        if _is_inner_wall_film_case(case_name):
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
    if _is_inner_wall_film_case(case_name):
        film_medium, meta = _film_medium_for_wavelength(wl, args)
        film_meta = {
            "material_name_film": meta["material_name"],
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
        courant=args.courant,
        solver_version="gap10_wide_deep_inner_wall_film_single_wavelength_v1",
        source_mode="single_wavelength_narrowband",
        logger=logger,
    )
    ti_flag = "VALID" if wl <= BAND_STRICT_HI else "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING"
    trans_flag = _transmission_flag(float(result["T"]))
    finite = all(
        np.isfinite(result[k])
        for k in ["R", "T", "A", "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"]
    )
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


def _failed_task_row(
    case_name: str,
    pol: str,
    wl: float,
    resolution: int,
    pml_um: float,
    args,
    exc: Exception,
) -> dict:
    thickness = _case_thickness(case_name)
    return {
        "mode": args.mode,
        "case_name": case_name,
        "film_thickness_um": thickness,
        "coating_mode": "none" if thickness == 0 else "sidewalls_and_bottom",
        "polarization": pol,
        "wavelength_um": wl,
        "R": np.nan,
        "T": np.nan,
        "A": np.nan,
        "input_flux_raw": np.nan,
        "reflection_flux_raw": np.nan,
        "transmission_flux_raw": np.nan,
        "raw_transmittance": np.nan,
        "signed_transmittance": np.nan,
        "material_validity_flag_ti": "VALID" if wl <= BAND_STRICT_HI else "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING",
        "transmission_quality_flag": "FAIL_EXCLUDED_FROM_ALL_METRICS",
        "included_in_strict_metric": False,
        "included_in_extended_metric": False,
        "finite_values_flag": "TASK_FAILED",
        "resolution": resolution,
        "pml_thickness_um": pml_um,
        "substrate_thickness_um": args.substrate_thickness_um,
        "air_buffer_um": args.air_buffer_um,
        "decay_db": args.decay_db,
        "fwidth_fraction": args.fwidth_fraction,
        "solver_version": "gap10_wide_deep_inner_wall_film_single_wavelength_v1",
        "source_mode": "single_wavelength_narrowband",
        "period_um": args.period_um,
        "top_width_um": args.top_width_um,
        "bottom_width_um": args.bottom_width_um,
        "depth_um": args.depth_um,
        "tilt_angle_deg": args.tilt_angle_deg,
        "film_grid_points": thickness * resolution,
        "normalization_note": "Task failed before valid R/T/A could be computed.",
        "walltime_s": np.nan,
        "task_error_type": type(exc).__name__,
        "task_error": str(exc),
        **_film_meta_empty(),
    }


def _route_class(mean_a: float) -> str:
    if not np.isfinite(mean_a):
        return "NOT_QUANTITATIVE"
    if mean_a < 0.30:
        return "GAP10_WIDE_DEEP_GROOVE_ROUTE_FAIL_FOR_HIGH_EMISSIVITY"
    if mean_a < CONTINUE_THRESHOLD:
        return "PARTIAL_ABSORPTION_ENHANCEMENT_BUT_BELOW_CONTINUE_THRESHOLD"
    if mean_a < HIGH_EMISSIVITY_TARGET:
        return "ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION"
    return "HIGH_EMISSIVITY_CANDIDATE_FOR_DIRECTIONAL_FOLLOWUP"


def _numerical_status(args, strict_count: int, strict_status: str) -> str:
    if args.mode == "smoke":
        return "CODE_PASS" if strict_count > 0 else "FAIL"
    enough = strict_count >= 0.9 * args.strict_planned_point_count
    if not enough:
        return "INCOMPLETE_NUMERICAL_COVERAGE"
    if strict_status != "valid":
        return strict_status.upper()
    return "NUMERICAL_PASS" if args.mode == "refine" else "NUMERICAL_SCREENING"


def _metric_from_group(
    df: pd.DataFrame,
    args,
    *,
    case: str,
    pol: str,
    thickness: float,
    resolution: int,
    pml_um: float,
    proxy_note: str = "",
) -> dict:
    strict_mask = df["included_in_strict_metric"].to_numpy(dtype=bool)
    strict_count = int(strict_mask.sum())
    mean_strict = wavelength_integrated_average(
        df["A"].to_numpy(dtype=float),
        df["wavelength_um"].to_numpy(dtype=float),
        BAND_STRICT_LO,
        BAND_STRICT_HI,
        valid_mask=strict_mask,
    )
    strict_status = wavelength_integrated_average.last_status
    numerical_status = _numerical_status(args, strict_count, strict_status)
    if numerical_status == "INCOMPLETE_NUMERICAL_COVERAGE":
        mean_strict = np.nan

    ext_mask = df["included_in_extended_metric"].to_numpy(dtype=bool)
    mean_extended = wavelength_integrated_average(
        df["A"].to_numpy(dtype=float),
        df["wavelength_um"].to_numpy(dtype=float),
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
    return {
        "case_name": case,
        "polarization": pol,
        "film_thickness_um": thickness,
        "coating_mode": "none" if thickness == 0 else "sidewalls_and_bottom",
        "resolution": resolution,
        "pml_thickness_um": pml_um,
        "mean_A_8p1014_12p398_strict": mean_strict,
        "mean_A_8p1014_12p962_extended": mean_extended,
        "extended_metric_warning": "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING",
        "peak_A_strict": peak_a,
        "peak_wavelength_strict_um": peak_wl,
        "valid_point_count_strict": strict_count,
        "total_planned_point_count_strict": args.strict_planned_point_count,
        "valid_point_count_extended": int(ext_mask.sum()),
        "numerical_status": numerical_status,
        "max_abs_T": float(np.nanmax(np.abs(df["T"].to_numpy(dtype=float)))),
        "route_decision": _route_class(mean_strict),
        "physics_note": proxy_note,
    }


def _compute_metrics(spectra: pd.DataFrame, args) -> pd.DataFrame:
    rows: list[dict] = []
    group_cols = ["case_name", "polarization", "film_thickness_um", "resolution", "pml_thickness_um"]
    for keys, df in spectra.groupby(group_cols, sort=False):
        case, pol, thickness, resolution, pml_um = keys
        rows.append(
            _metric_from_group(
                df,
                args,
                case=case,
                pol=pol,
                thickness=float(thickness),
                resolution=int(resolution),
                pml_um=float(pml_um),
            )
        )

    for (case, thickness, resolution, pml_um), df in spectra.groupby(
        ["case_name", "film_thickness_um", "resolution", "pml_thickness_um"],
        sort=False,
    ):
        pivot = df.pivot_table(
            index="wavelength_um",
            columns="polarization",
            values=["A", "included_in_strict_metric", "included_in_extended_metric", "T"],
            aggfunc="first",
        )
        if not {"Ez", "Hz"}.issubset(set(pivot["A"].columns)):
            continue
        wavelengths = pivot.index.to_numpy(dtype=float)
        proxy_a = 0.5 * (
            pivot[("A", "Ez")].to_numpy(dtype=float)
            + pivot[("A", "Hz")].to_numpy(dtype=float)
        )
        strict_valid = (
            pivot[("included_in_strict_metric", "Ez")].to_numpy(dtype=bool)
            & pivot[("included_in_strict_metric", "Hz")].to_numpy(dtype=bool)
        )
        extended_valid = (
            pivot[("included_in_extended_metric", "Ez")].to_numpy(dtype=bool)
            & pivot[("included_in_extended_metric", "Hz")].to_numpy(dtype=bool)
        )
        proxy_t = 0.5 * (
            pivot[("T", "Ez")].to_numpy(dtype=float)
            + pivot[("T", "Hz")].to_numpy(dtype=float)
        )
        proxy_df = pd.DataFrame(
            {
                "wavelength_um": wavelengths,
                "A": proxy_a,
                "T": proxy_t,
                "included_in_strict_metric": strict_valid,
                "included_in_extended_metric": extended_valid,
            }
        )
        rows.append(
            _metric_from_group(
                proxy_df,
                args,
                case=case,
                pol=PROXY,
                thickness=float(thickness),
                resolution=int(resolution),
                pml_um=float(pml_um),
                proxy_note=(
                    "A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2; this is a 2D "
                    "polarization-average proxy and is not a true 3D non-polarized emissivity."
                ),
            )
        )

    metrics = pd.DataFrame(rows)
    if metrics.empty:
        return metrics
    metrics = _add_enhancements(metrics)
    return metrics


def _append_proxy_spectra(spectra: pd.DataFrame, args) -> pd.DataFrame:
    proxy_rows: list[dict] = []
    for (case, thickness, resolution, pml_um), df in spectra.groupby(
        ["case_name", "film_thickness_um", "resolution", "pml_thickness_um"],
        sort=False,
    ):
        pivot = df.pivot_table(
            index="wavelength_um",
            columns="polarization",
            values=["R", "T", "A", "included_in_strict_metric", "included_in_extended_metric"],
            aggfunc="first",
        )
        if not {"Ez", "Hz"}.issubset(set(pivot["A"].columns)):
            continue
        for wl in pivot.index.to_numpy(dtype=float):
            row = {
                "mode": args.mode,
                "case_name": case,
                "film_thickness_um": float(thickness),
                "coating_mode": "none" if float(thickness) == 0 else "sidewalls_and_bottom",
                "polarization": PROXY,
                "wavelength_um": float(wl),
                "R": 0.5 * (float(pivot.loc[wl, ("R", "Ez")]) + float(pivot.loc[wl, ("R", "Hz")])),
                "T": 0.5 * (float(pivot.loc[wl, ("T", "Ez")]) + float(pivot.loc[wl, ("T", "Hz")])),
                "A": 0.5 * (float(pivot.loc[wl, ("A", "Ez")]) + float(pivot.loc[wl, ("A", "Hz")])),
                "included_in_strict_metric": bool(
                    pivot.loc[wl, ("included_in_strict_metric", "Ez")]
                    and pivot.loc[wl, ("included_in_strict_metric", "Hz")]
                ),
                "included_in_extended_metric": bool(
                    pivot.loc[wl, ("included_in_extended_metric", "Ez")]
                    and pivot.loc[wl, ("included_in_extended_metric", "Hz")]
                ),
                "resolution": int(resolution),
                "pml_thickness_um": float(pml_um),
                "source_mode": "derived_from_Ez_Hz",
                "physics_note": "A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2; not a true 3D emissivity.",
            }
            proxy_rows.append(row)
    if not proxy_rows:
        return spectra
    return pd.concat([spectra, pd.DataFrame(proxy_rows)], ignore_index=True, sort=False)


def _add_enhancements(metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = metrics.copy()
    key_cols = ["polarization", "resolution", "pml_thickness_um"]
    d11_ref = _load_d11_reference_metrics()
    for idx, row in metrics.iterrows():
        same = metrics[
            (metrics["polarization"] == row["polarization"])
            & (metrics["resolution"] == row["resolution"])
            & (metrics["pml_thickness_um"] == row["pml_thickness_um"])
        ]
        for ref_name, out_name in [(FLAT, "flat_Ti"), (BARE, "gap10_bare")]:
            ref = same[same["case_name"] == ref_name]
            abs_col = f"enhancement_over_{out_name}_absolute"
            rel_col = f"enhancement_over_{out_name}_relative"
            if ref.empty or not np.isfinite(row["mean_A_8p1014_12p398_strict"]):
                metrics.loc[idx, abs_col] = np.nan
                metrics.loc[idx, rel_col] = np.nan
                continue
            ref_a = float(ref["mean_A_8p1014_12p398_strict"].iloc[0])
            if not np.isfinite(ref_a):
                metrics.loc[idx, abs_col] = np.nan
                metrics.loc[idx, rel_col] = np.nan
                continue
            delta = float(row["mean_A_8p1014_12p398_strict"] - ref_a)
            metrics.loc[idx, abs_col] = delta
            metrics.loc[idx, rel_col] = delta / ref_a if ref_a != 0 else np.nan
        ref_a = d11_ref.get((float(row["film_thickness_um"]), row["polarization"]))
        if ref_a is None or not np.isfinite(row["mean_A_8p1014_12p398_strict"]):
            metrics.loc[idx, "enhancement_over_gap30_D11_reference_absolute"] = np.nan
        else:
            metrics.loc[idx, "enhancement_over_gap30_D11_reference_absolute"] = (
                float(row["mean_A_8p1014_12p398_strict"]) - ref_a
            )
    return metrics


def _load_d11_reference_metrics() -> dict[tuple[float, str], float]:
    candidates = [
        project_path("results", "diagnostics_v2", "tables", "D11_scaled_deep_groove_inner_wall_film_metrics.csv"),
        project_path("results", "diagnostics_v2", "tables", "D11_thickness_1um_quick_metrics.csv"),
        project_path("results", "diagnostics_v2", "tables", "D11_thickness_2um_quick_metrics.csv"),
        project_path("results", "diagnostics_v2", "tables", "D11_thickness_3um_quick_metrics.csv"),
    ]
    refs: dict[tuple[float, str], float] = {}
    for path in candidates:
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        required = {"film_thickness_um", "polarization", "mean_A_8p1014_12p398_strict"}
        if not required.issubset(df.columns):
            continue
        df = df[df["polarization"].isin(["Ez", "Hz", PROXY])].copy()
        for row in df.itertuples(index=False):
            try:
                key = (float(row.film_thickness_um), str(row.polarization))
                val = float(row.mean_A_8p1014_12p398_strict)
            except (TypeError, ValueError):
                continue
            if np.isfinite(val):
                refs[key] = val
    return refs


def _resolution_check(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame([{"resolution_check_status": "NOT_RUN"}])
    rows = []
    for (case, pol, pml_um), df in metrics.groupby(["case_name", "polarization", "pml_thickness_um"], sort=False):
        if not {48, 64}.issubset(set(df["resolution"].astype(int))):
            continue
        a48 = float(df[df["resolution"] == 48]["mean_A_8p1014_12p398_strict"].iloc[0])
        a64 = float(df[df["resolution"] == 64]["mean_A_8p1014_12p398_strict"].iloc[0])
        delta = abs(a64 - a48)
        rows.append({
            "case_name": case,
            "polarization": pol,
            "pml_thickness_um": pml_um,
            "mean_A_res48": a48,
            "mean_A_res64": a64,
            "abs_delta_mean_A": delta,
            "resolution_check_status": "NUMERICAL_PASS" if np.isfinite(delta) and delta < 0.01 else "WARNING",
        })
    return pd.DataFrame(rows) if rows else pd.DataFrame([{"resolution_check_status": "NOT_RUN"}])


def _pml_check(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame([{"pml_check_status": "NOT_RUN"}])
    rows = []
    for (case, pol, res), df in metrics.groupby(["case_name", "polarization", "resolution"], sort=False):
        if not {4.0, 6.0}.issubset(set(df["pml_thickness_um"].astype(float))):
            continue
        a4 = float(df[df["pml_thickness_um"] == 4.0]["mean_A_8p1014_12p398_strict"].iloc[0])
        a6 = float(df[df["pml_thickness_um"] == 6.0]["mean_A_8p1014_12p398_strict"].iloc[0])
        delta = abs(a6 - a4)
        rows.append({
            "case_name": case,
            "polarization": pol,
            "resolution": res,
            "mean_A_pml4": a4,
            "mean_A_pml6": a6,
            "abs_delta_mean_A": delta,
            "pml_check_status": "NUMERICAL_PASS" if np.isfinite(delta) and delta < 0.01 else "WARNING",
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
    ax.add_patch(
        plt.Rectangle(
            (-half, -args.substrate_thickness_um),
            args.period_um,
            args.substrate_thickness_um,
            facecolor="#808080",
            edgecolor="black",
            label="Ti substrate",
        )
    )
    if thickness <= 0:
        ax.fill(
            *zip(*(outer + [outer[0]])),
            facecolor="white",
            edgecolor="#1F77B4",
            lw=1.4,
            label="air groove core",
        )
    else:
        ax.fill(
            *zip(*(outer + [outer[0]])),
            facecolor="#D62728",
            alpha=0.75,
            edgecolor="black",
            lw=1.1,
            label="measured_lossy_wall_film",
        )
        inner = slanted_groove_vertices(
            top_width_um=args.top_width_um - 2.0 * thickness,
            bottom_width_um=args.bottom_width_um - 2.0 * thickness,
            depth_um=args.depth_um - thickness,
            tilt_angle_deg=args.tilt_angle_deg,
            y_surface=y_surface,
        )
        ax.fill(
            *zip(*(inner + [inner[0]])),
            facecolor="white",
            edgecolor="#1F77B4",
            lw=1.4,
            label="air groove core",
        )
    ax.axhline(y_surface, color="black", lw=1.0)
    ax.text(-half + 0.5, y_surface + 0.35, "top Ti land remains uncoated", fontsize=8)
    ax.text(half - 12.5, y_surface + 0.35, "bare top Ti surface", fontsize=8)
    ax.set_xlim(-half, half)
    ax.set_ylim(-args.depth_um - 1.0, y_surface + 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"D12 gap10 wide/deep groove geometry, inner film t={thickness:g} um")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.2)
    save_figure(fig, out_path)
    plt.close(fig)


def _generate_geometry_plots(paths: dict[str, Path], args) -> None:
    _plot_geometry(0.0, paths["geom_bare"], args)
    _plot_geometry(0.50, paths["geom_500"], args)
    _plot_geometry(1.0, paths["geom_1um"], args)
    _plot_geometry(2.0, paths["geom_2um"], args)


def _plot_spectra(spectra: pd.DataFrame, paths: dict[str, Path], args) -> None:
    for pol, out in [("Ez", paths["ez"]), ("Hz", paths["hz"])]:
        fig, ax = plt.subplots(figsize=(8.2, 5.0))
        df_pol = spectra[spectra["polarization"] == pol]
        if df_pol.empty:
            ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
            ax.set_axis_off()
        else:
            main = df_pol[
                (df_pol["resolution"] == min(args.resolutions))
                & (df_pol["pml_thickness_um"] == DEFAULT_PML_UM)
            ]
            if main.empty:
                main = df_pol
            for case, df in main.groupby("case_name", sort=False):
                ax.plot(df["wavelength_um"], df["A"], marker="o", lw=1.2, label=_case_label(case))
            ax.axvspan(BAND_STRICT_LO, BAND_STRICT_HI, color="#E8EEF7", alpha=0.45)
            ax.set_xlabel("Wavelength (um)")
            ax.set_ylabel("A")
            ax.set_title(f"D12 {pol} absorptance comparison")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7)
        save_figure(fig, out)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    pivot = spectra.pivot_table(
        index=["case_name", "wavelength_um", "resolution", "pml_thickness_um"],
        columns="polarization",
        values="A",
        aggfunc="first",
    ).reset_index()
    if {"Ez", "Hz"}.issubset(set(pivot.columns)):
        pivot["A_unpolarized_2D_proxy"] = 0.5 * (pivot["Ez"] + pivot["Hz"])
        main = pivot[
            (pivot["resolution"] == min(args.resolutions))
            & (pivot["pml_thickness_um"] == DEFAULT_PML_UM)
        ]
        if main.empty:
            main = pivot
        for case, df in main.groupby("case_name", sort=False):
            ax.plot(df["wavelength_um"], df["A_unpolarized_2D_proxy"], marker="o", lw=1.2, label=_case_label(case))
        ax.axvspan(BAND_STRICT_LO, BAND_STRICT_HI, color="#E8EEF7", alpha=0.45)
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("A_unpolarized_2D_proxy")
        ax.set_title("D12 unpolarized 2D proxy")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "No paired Ez/Hz data", ha="center", va="center")
        ax.set_axis_off()
    save_figure(fig, paths["proxy"])
    plt.close(fig)


def _plot_metric_bars(metrics: pd.DataFrame, paths: dict[str, Path]) -> None:
    main = metrics[
        (metrics["pml_thickness_um"].astype(float) == DEFAULT_PML_UM)
        & (metrics["resolution"].astype(int) == metrics["resolution"].astype(int).min())
    ].copy()

    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    if main.empty:
        ax.text(0.5, 0.5, "No metrics", ha="center", va="center")
        ax.set_axis_off()
    else:
        labels = [f"{_case_label(r.case_name)}\n{r.polarization}" for r in main.itertuples()]
        x = np.arange(len(main))
        ax.bar(x, main["mean_A_8p1014_12p398_strict"].astype(float), color="#4C78A8")
        ax.axhline(CONTINUE_THRESHOLD, color="red", ls="--", lw=1, label="0.60 continue")
        ax.axhline(HIGH_EMISSIVITY_TARGET, color="purple", ls=":", lw=1, label="0.80 target")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Mean A strict")
        ax.set_title("D12 route decision")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=8)
    save_figure(fig, paths["mean"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    enh = main[main["case_name"].isin([FILM_500, FILM_1UM, FILM_2UM])].copy()
    if "enhancement_over_gap10_bare_absolute" in enh.columns and not enh.empty:
        labels = [f"{_case_label(r.case_name)}\n{r.polarization}" for r in enh.itertuples()]
        x = np.arange(len(enh))
        ax.bar(x, enh["enhancement_over_gap10_bare_absolute"].astype(float), color="#59A14F")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Delta mean A vs gap10 bare")
        ax.set_title("D12 enhancement over bare gap10 groove")
        ax.grid(axis="y", alpha=0.25)
    else:
        ax.text(0.5, 0.5, "No enhancement metrics", ha="center", va="center")
        ax.set_axis_off()
    save_figure(fig, paths["enhancement"])
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    if "enhancement_over_gap30_D11_reference_absolute" in main.columns and not main.empty:
        refs = main[main["case_name"].isin([FILM_500, FILM_1UM, FILM_2UM])].copy()
        refs = refs[refs["enhancement_over_gap30_D11_reference_absolute"].notna()]
        if not refs.empty:
            labels = [f"{_case_label(r.case_name)}\n{r.polarization}" for r in refs.itertuples()]
            x = np.arange(len(refs))
            ax.bar(x, refs["enhancement_over_gap30_D11_reference_absolute"].astype(float), color="#F28E2B")
            ax.axhline(0, color="black", lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
            ax.set_ylabel("Delta mean A vs D11 gap30 reference")
            ax.set_title("D12 gap10 vs D11 gap30 reference")
            ax.grid(axis="y", alpha=0.25)
        else:
            ax.text(0.5, 0.5, "No matched D11 reference metrics", ha="center", va="center")
            ax.set_axis_off()
    else:
        ax.text(0.5, 0.5, "No D11 reference metrics", ha="center", va="center")
        ax.set_axis_off()
    save_figure(fig, paths["gap10_vs_gap30"])
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


def _draw_geometry_overlay(ax, args, thickness: float) -> None:
    y_surface = 0.0
    outer = slanted_groove_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=y_surface,
    )
    ox, oy = zip(*(outer + [outer[0]]))
    ax.plot(ox, oy, color="white", lw=1.2)
    if thickness > 0:
        inner = slanted_groove_vertices(
            top_width_um=args.top_width_um - 2.0 * thickness,
            bottom_width_um=args.bottom_width_um - 2.0 * thickness,
            depth_um=args.depth_um - thickness,
            tilt_angle_deg=args.tilt_angle_deg,
            y_surface=y_surface,
        )
        ix, iy = zip(*(inner + [inner[0]]))
        ax.plot(ix, iy, color="#00FFFF", lw=1.0)


def _plot_field_map(data: np.ndarray, title: str, path: Path, extent, args, thickness: float) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    vmax = np.nanpercentile(data, 99.5) if np.isfinite(data).any() else None
    im = ax.imshow(
        data.T,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="inferno",
        vmin=0,
        vmax=vmax,
    )
    _draw_geometry_overlay(ax, args, thickness)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    save_figure(fig, path)
    plt.close(fig)


def _write_placeholder_field_figures(paths: dict[str, Path], reason: str, args) -> dict:
    for key, title in [
        ("best_e2", "D12 best-case |E|^2 not generated"),
        ("best_h2", "D12 best-case |H|^2 not generated"),
    ]:
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        ax.text(0.5, 0.5, reason, ha="center", va="center", wrap=True)
        ax.set_title(title)
        ax.set_axis_off()
        save_figure(fig, paths[key])
        plt.close(fig)
    _plot_geometry(0.50, paths["best_overlay"], args)
    return {"field_snapshot_status": "NOT_RUN", "field_snapshot_note": reason}


def _run_best_field_snapshot(spectra: pd.DataFrame, metrics: pd.DataFrame, paths: dict[str, Path], args, logger) -> dict:
    import meep as mp

    if args.mode != "screen":
        return _write_placeholder_field_figures(
            paths,
            "Field maps are only gated from screen-mode metrics.",
            args,
        )
    candidates = metrics[
        (metrics["case_name"].isin([BARE, FILM_500, FILM_1UM, FILM_2UM]))
        & (metrics["polarization"].isin(["Ez", "Hz"]))
        & (metrics["mean_A_8p1014_12p398_strict"].fillna(-1) >= FIELD_MAP_THRESHOLD)
    ].copy()
    if candidates.empty:
        return _write_placeholder_field_figures(
            paths,
            "Screen did not reach mean_A_strict >= 0.30 for any gap10 structure.",
            args,
        )

    best_metric = candidates.sort_values("mean_A_8p1014_12p398_strict", ascending=False).iloc[0]
    pol = str(best_metric["polarization"])
    resolution = int(best_metric["resolution"])
    pml_um = float(best_metric["pml_thickness_um"])
    peak_wl = float(best_metric["peak_wavelength_strict_um"])
    case_name = str(best_metric["case_name"])
    thickness = _case_thickness(case_name)
    if not np.isfinite(peak_wl):
        return _write_placeholder_field_figures(paths, "No valid strict peak wavelength was available.", args)

    film = None
    if thickness > 0:
        film, _ = _film_medium_for_wavelength(peak_wl, args)
    ti = get_ti_medium()
    air = mp.Medium(epsilon=1.0)
    bottom_buffer = pml_um
    cell_y = 2 * pml_um + args.air_buffer_um + args.substrate_thickness_um + bottom_buffer
    y_top_edge = 0.5 * cell_y
    y_top_pml_inner = y_top_edge - pml_um
    y_surface = y_top_pml_inner - args.air_buffer_um
    y_src = y_top_pml_inner - 0.25 * args.air_buffer_um
    y_refl = y_surface + 0.5 * args.air_buffer_um
    y_dft_min = y_surface - args.depth_um - 0.8
    y_dft_max = y_surface + 1.2
    y_center = 0.5 * (y_dft_min + y_dft_max)
    y_size = y_dft_max - y_dft_min
    if thickness > 0:
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
            coating_mode="sidewalls_and_bottom",
        )
    else:
        geometry = build_slanted_groove_geometry(
            period_x_um=args.period_um,
            top_width_um=args.top_width_um,
            bottom_width_um=args.bottom_width_um,
            depth_um=args.depth_um,
            tilt_angle_deg=args.tilt_angle_deg,
            substrate_thickness_um=args.substrate_thickness_um,
            y_surface=y_surface,
            medium_substrate=ti,
            medium_groove=air,
        )
    e_comps, h_comps, src_c = _field_components(pol)
    fcen = 1.0 / peak_wl
    sim = mp.Simulation(
        cell_size=mp.Vector3(args.period_um, cell_y, 0),
        boundary_layers=[mp.PML(thickness=pml_um, direction=mp.Y)],
        sources=[
            mp.Source(
                mp.GaussianSource(
                    frequency=fcen,
                    fwidth=args.fwidth_fraction * fcen,
                    is_integrated=True,
                ),
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
    logger.info(">>> D12 field snapshot case=%s pol=%s wl=%.6g", case_name, pol, peak_wl)
    try:
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
        extent = (-0.5 * args.period_um, 0.5 * args.period_um, y_dft_min - y_surface, y_dft_max - y_surface)
        _plot_field_map(e2, f"D12 |E|^2, {pol}, {peak_wl:g} um", paths["best_e2"], extent, args, thickness)
        _plot_field_map(h2, f"D12 |H|^2, {pol}, {peak_wl:g} um", paths["best_h2"], extent, args, thickness)
        _plot_geometry(thickness, paths["best_overlay"], args)
        return {
            "field_snapshot_status": "COMPLETED",
            "field_snapshot_case": case_name,
            "field_snapshot_polarization": pol,
            "field_snapshot_wavelength_um": peak_wl,
            "field_snapshot_resolution": resolution,
            "field_snapshot_pml_thickness_um": pml_um,
        }
    except Exception as exc:
        logger.exception("D12 field snapshot failed; preserving spectra/metrics.")
        return _write_placeholder_field_figures(paths, f"Field snapshot failed: {type(exc).__name__}: {exc}", args)


def _md_table(metrics: pd.DataFrame) -> str:
    cols = [
        "case_name",
        "polarization",
        "mean_A_8p1014_12p398_strict",
        "peak_A_strict",
        "enhancement_over_gap10_bare_absolute",
        "enhancement_over_gap30_D11_reference_absolute",
        "enhancement_over_flat_Ti_absolute",
        "route_decision",
        "numerical_status",
    ]
    if metrics.empty:
        return "_No metrics._"
    df = metrics[[c for c in cols if c in metrics.columns]].copy()
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


def _overall_level(metrics: pd.DataFrame, args) -> str:
    if metrics.empty:
        return "FAIL"
    if args.mode == "smoke":
        return "CODE_PASS" if (metrics["numerical_status"] == "CODE_PASS").any() else "FAIL"
    if args.mode == "screen":
        return "NUMERICAL_SCREENING"
    status = set(metrics["numerical_status"].dropna().astype(str))
    return "NUMERICAL_PASS" if status == {"NUMERICAL_PASS"} else "WARNING"


def _best_film(metrics: pd.DataFrame) -> pd.Series | None:
    required = {"case_name", "mean_A_8p1014_12p398_strict"}
    if metrics.empty or not required.issubset(metrics.columns):
        return None
    main = metrics[
        (metrics["case_name"].isin([FILM_500, FILM_1UM, FILM_2UM]))
        & metrics["mean_A_8p1014_12p398_strict"].notna()
    ].copy()
    if main.empty:
        return None
    return main.sort_values("mean_A_8p1014_12p398_strict", ascending=False).iloc[0]


def _metric_lookup(metrics: pd.DataFrame, case: str, pol: str) -> pd.Series | None:
    required = {"case_name", "polarization", "pml_thickness_um"}
    if metrics.empty or not required.issubset(metrics.columns):
        return None
    df = metrics[
        (metrics["case_name"] == case)
        & (metrics["polarization"] == pol)
        & (metrics["pml_thickness_um"].astype(float) == DEFAULT_PML_UM)
    ].copy()
    if df.empty:
        return None
    df = df.sort_values("resolution")
    return df.iloc[0]


def _value_text(row: pd.Series | None, col: str) -> str:
    if row is None or col not in row or pd.isna(row[col]):
        return "not available"
    return f"{float(row[col]):.6g}"


def _synergy_text(metrics: pd.DataFrame) -> str:
    pieces = []
    for pol in ["Ez", "Hz", PROXY]:
        flat = _metric_lookup(metrics, FLAT, pol)
        bare = _metric_lookup(metrics, BARE, pol)
        films = [
            _metric_lookup(metrics, case, pol)
            for case in [FILM_500, FILM_1UM, FILM_2UM]
        ]
        films = [row for row in films if row is not None]
        film = None if not films else max(
            films,
            key=lambda r: float(r["mean_A_8p1014_12p398_strict"]) if pd.notna(r["mean_A_8p1014_12p398_strict"]) else -np.inf,
        )
        if flat is None or bare is None or film is None:
            continue
        flat_a = float(flat["mean_A_8p1014_12p398_strict"])
        bare_a = float(bare["mean_A_8p1014_12p398_strict"])
        film_a = float(film["mean_A_8p1014_12p398_strict"])
        if not all(np.isfinite([flat_a, bare_a, film_a])):
            continue
        geom_gain = bare_a - flat_a
        film_gain = film_a - bare_a
        if abs(geom_gain) >= abs(film_gain) * 1.5:
            label = "mainly geometry scaling"
        elif abs(film_gain) >= abs(geom_gain) * 1.5:
            label = "mainly inner-wall film coupling"
        else:
            label = "geometry and film are comparable / synergistic"
        pieces.append(f"{pol}: geometry gain={geom_gain:.4g}, film gain={film_gain:.4g}, {label}")
    return "; ".join(pieces) if pieces else "not enough valid metrics to separate geometry and film contributions"


def _write_report(
    paths: dict[str, Path],
    spectra: pd.DataFrame,
    metrics: pd.DataFrame,
    res_check: pd.DataFrame,
    pml_check: pd.DataFrame,
    geometry_checks: pd.DataFrame,
    field_info: dict,
    args,
) -> None:
    ensure_dir(paths["report"].parent)
    level = _overall_level(metrics, args)
    geometry_failed = bool((not geometry_checks.empty) and (geometry_checks["status"] != "PASS").any())
    geometry_only = bool(getattr(args, "geometry_only", False))
    if geometry_failed:
        level = "FAIL"
    elif geometry_only:
        level = "GEOMETRY_PASS"
    bottom_row = geometry_checks[geometry_checks["check_name"] == "tilted_bottom_vertices_inside_unit_cell"]
    bottom_ok = (not bottom_row.empty) and str(bottom_row["status"].iloc[0]) == "PASS"
    best = _best_film(metrics)
    best_mean = np.nan if best is None else float(best["mean_A_8p1014_12p398_strict"])
    pass60 = bool(np.isfinite(best_mean) and best_mean >= CONTINUE_THRESHOLD)
    pass80 = bool(np.isfinite(best_mean) and best_mean >= HIGH_EMISSIVITY_TARGET)
    best_text = (
        "No strict quantitative film metric is available."
        if best is None
        else (
            f"Best film strict mean A = {best_mean:.6g} "
            f"({best['case_name']}, {best['polarization']}, res={int(best['resolution'])}, "
            f"PML={float(best['pml_thickness_um']):g} um)."
        )
    )
    if geometry_failed:
        decision = (
            "GEOMETRY_INVALID_DO_NOT_RUN_MEEP. The requested gap=10 um, wide/deep, "
            f"depth={args.depth_um:g} um, tilt={args.tilt_angle_deg:g} deg groove leaves the P={args.period_um:g} um periodic unit cell. "
            "Use a straight slot, reduce tilt, reduce bottom width, or change the lateral placement before "
            "using this as a quantitative absorption route."
        )
    elif geometry_only:
        decision = (
            "GEOMETRY_ONLY_PASS. The geometry checks passed, but no Meep absorption "
            "simulation was requested in this geometry-only run."
        )
    elif pass80:
        decision = (
            "HIGH_EMISSIVITY_CANDIDATE_FOR_DIRECTIONAL_FOLLOWUP. Because P=50 um "
            "supports multiple propagating diffraction orders, the next step must "
            "use mode decomposition before making a directional-emission claim."
        )
    elif pass60:
        decision = (
            "ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION. Continue with deep-groove "
            "dimension, film-thickness, and coverage optimization plus stricter resolution/PML checks."
        )
    else:
        decision = (
            "GAP10_WIDE_DEEP_GROOVE_ROUTE_FAIL_FOR_HIGH_EMISSIVITY or below the continuation "
            "threshold. Even with gap=10 um, width=40 um, and depth=30 um, the inner-wall-only "
            "film route should not be treated as the main high-emissivity route if it remains "
            "below 0.60."
        )

    rows_500 = {
        pol: _metric_lookup(metrics, FILM_500, pol)
        for pol in ["Ez", "Hz", PROXY]
    }
    rows_1um = {
        pol: _metric_lookup(metrics, FILM_1UM, pol)
        for pol in ["Ez", "Hz", PROXY]
    }
    rows_2um = {
        pol: _metric_lookup(metrics, FILM_2UM, pol)
        for pol in ["Ez", "Hz", PROXY]
    }
    rows_bare = {
        pol: _metric_lookup(metrics, BARE, pol)
        for pol in ["Ez", "Hz", PROXY]
    }
    rows_flat = {
        pol: _metric_lookup(metrics, FLAT, pol)
        for pol in ["Ez", "Hz", PROXY]
    }
    period_note = "This straight-slot variant keeps the original D12 P=50 um period."
    answers = f"""1. P={args.period_um:g} um, gap={args.period_um - args.top_width_um:g} um, top/bottom width={args.top_width_um:g}/{args.bottom_width_um:g} um, h={args.depth_um:g} um, tilt={args.tilt_angle_deg:g} deg geometry built?
   {'No. The requested nominal geometry fails the periodic-cell checks, so no Meep absorption run was performed.' if geometry_failed else 'Yes. Geometry checks passed.'} {period_note} Substrate={args.substrate_thickness_um:g} um.

2. Does the tilted bottom remain fully inside the unit cell?
   {'Yes.' if bottom_ok else 'No.'} See `outer_groove_inside_unit_cell` and `inner_air_core_inside_unit_cell_*` in the geometry checks. For this nominal geometry, the bottom center offset is about {args.depth_um * np.tan(np.deg2rad(args.tilt_angle_deg)):.4g} um.

3. Is the film only on inner sidewalls and bottom?
   Yes by construction: the script uses `build_inner_wall_film_slanted_groove_geometry(..., coating_mode="sidewalls_and_bottom")`, an outer groove film prism, and an inner air prism that opens at the top surface. No top-land film block is added. The geometry PNGs mark the top Ti land as uncoated.

4. Gap10 500 nm improvement over bare gap10 groove:
   Ez: {_value_text(rows_500['Ez'], 'enhancement_over_gap10_bare_absolute')} absolute; Hz: {_value_text(rows_500['Hz'], 'enhancement_over_gap10_bare_absolute')} absolute; 2D proxy: {_value_text(rows_500[PROXY], 'enhancement_over_gap10_bare_absolute')} absolute.

5. Does the bare gap10 wide/deep groove improve over flat Ti?
   Ez: bare mean {_value_text(rows_bare['Ez'], 'mean_A_8p1014_12p398_strict')} vs flat {_value_text(rows_flat['Ez'], 'mean_A_8p1014_12p398_strict')}; Hz: bare mean {_value_text(rows_bare['Hz'], 'mean_A_8p1014_12p398_strict')} vs flat {_value_text(rows_flat['Hz'], 'mean_A_8p1014_12p398_strict')}; proxy: bare mean {_value_text(rows_bare[PROXY], 'mean_A_8p1014_12p398_strict')} vs flat {_value_text(rows_flat[PROXY], 'mean_A_8p1014_12p398_strict')}.

6. Which film thickness is best among 0.5, 1, and 2 um?
   Best available film row: {best_text}

7. Best result source:
   The best available row identifies whether it came from Ez, Hz, or the 2D proxy: {best_text}

8. Comparison with old D11 gap=30 references:
   500 nm proxy delta: {_value_text(rows_500[PROXY], 'enhancement_over_gap30_D11_reference_absolute')}; 1 um proxy delta: {_value_text(rows_1um[PROXY], 'enhancement_over_gap30_D11_reference_absolute')}; 2 um proxy delta: {_value_text(rows_2um[PROXY], 'enhancement_over_gap30_D11_reference_absolute')}.

9. Enhancement source:
   {_synergy_text(metrics)}.

10. Does it reach mean_A_strict >= 0.60?
   {'Not evaluated because the geometry is invalid' if geometry_failed else ('Yes' if pass60 else 'No')}. {best_text}

11. Does it reach mean_A_strict >= 0.80?
   {'Not evaluated because the geometry is invalid' if geometry_failed else ('Yes' if pass80 else 'No')}.

12. If below 0.60:
   If a valid geometry later remains below 0.60, then even with gap={args.period_um - args.top_width_um:g} um, top width={args.top_width_um:g} um, bottom width={args.bottom_width_um:g} um, and depth={args.depth_um:g} um, the inner-wall-only coating route is not suitable as the main high-emissivity route. The next step should shift to full-surface high-absorption film plus microstructure directionality control.

13. If between 0.60 and 0.80:
   Continue with deep-groove size, film thickness, and coverage optimization, and perform stricter resolution validation.

14. If >=0.80:
   Proceed to directionality testing, but because P={args.period_um:g} um supports multiple diffraction orders in this band, use mode decomposition to analyze order-resolved directionality.

15. 2D scope:
   The current result is a 2D equivalent straight-groove model, not a true 3D laser-processed hole array.
"""

    report = f"""# D12 Gap10 Wide/Deep Groove Inner-Wall Film Screen

Overall result level: **{level}**

## Required Answers

{answers}

## Route Decision

{decision}

## Quantitative Scope

- Formal metric: `mean_A_8p1014_12p398_strict`.
- Extended observation metric: `mean_A_8p1014_12p962_extended`, flagged with `TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING`.
- `A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2` is only a 2D polarization-average proxy and is not the non-polarized thermal emissivity of a true 3D laser-processed hole array.
- Flux monitors cover the full period, so R/T/A are total powers. This is suitable for total absorption screening, not for directional-emission conclusions.

## Geometry Checks

{_md_table_geometry(geometry_checks)}

## Metrics Preview

{_md_table(metrics)}

## Numerical Checks

- Resolution check rows: {len(res_check)}
- PML check rows: {len(pml_check)}
- Field snapshot status: {field_info.get('field_snapshot_status', 'UNKNOWN')}
- Field snapshot note: {field_info.get('field_snapshot_note', '')}

## Outputs

- Spectra: `{paths['spectra']}`
- Metrics: `{paths['metrics']}`
- Resolution check: `{paths['resolution']}`
- PML check: `{paths['pml']}`
- Log: `{paths['log']}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def _md_table_geometry(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No geometry checks._"
    lines = [
        "| check_name | status | details |",
        "| --- | --- | --- |",
    ]
    for row in df.itertuples(index=False):
        details = str(row.details).replace("|", "/")
        lines.append(f"| {row.check_name} | {row.status} | {details} |")
    return "\n".join(lines)


def _write_failure_report(paths: dict[str, Path], exc: Exception, args) -> None:
    ensure_dir(paths["report"].parent)
    paths["report"].write_text(
        "# D12 Gap10 Wide/Deep Groove Inner-Wall Film Screen\n\n"
        "Overall result level: **FAIL**\n\n"
        f"Failure: `{type(exc).__name__}: {exc}`\n\n"
        "No physical conclusion should be drawn from this failed run.\n",
        encoding="utf-8",
    )


def _make_plots(spectra: pd.DataFrame, metrics: pd.DataFrame, paths: dict[str, Path], args) -> None:
    _generate_geometry_plots(paths, args)
    _plot_spectra(spectra, paths, args)
    _plot_metric_bars(metrics, paths)


def _build_tasks(args) -> list[tuple[str, str, float, int, float]]:
    tasks = []
    for resolution in args.resolutions:
        for pml_um in args.pml_values_um:
            for case_name in args.cases:
                for pol in args.polarizations:
                    for wl in args.wavelengths_um:
                        tasks.append((case_name, pol, float(wl), int(resolution), float(pml_um)))
    return tasks


def _resolve_workers(args, task_count: int) -> int:
    if task_count <= 0:
        return 1
    requested = str(args.workers).strip().lower()
    cpu_count = os.cpu_count() or 1
    if requested == "auto":
        if args.mode == "smoke":
            return 1
        if args.mode == "refine":
            return max(1, min(4, task_count, cpu_count))
        return max(1, min(8, task_count, cpu_count))
    try:
        workers = int(requested)
    except ValueError as exc:
        raise ValueError("--workers must be a positive integer or 'auto'.") from exc
    if workers < 1:
        raise ValueError("--workers must be >= 1.")
    return min(workers, task_count)


def _task_sort_key(row: dict) -> tuple:
    case_order = {case: i for i, case in enumerate(ALL_CASES)}
    pol_order = {"Ez": 0, "Hz": 1, PROXY: 2}
    return (
        int(row.get("resolution", 0)),
        float(row.get("pml_thickness_um", 0.0)),
        case_order.get(row.get("case_name"), 99),
        pol_order.get(row.get("polarization"), 99),
        float(row.get("wavelength_um", 0.0)),
    )


def _run_one_task(payload: tuple[int, tuple[str, str, float, int, float], object, str]) -> tuple[int, dict]:
    index, task, args, log_path = payload
    case_name, pol, wl, resolution, pml_um = task
    logger = _setup_worker_logger(Path(log_path))
    logger.info(
        ">>> worker task=%d case=%s pol=%s wl=%g res=%s pml=%g",
        index,
        case_name,
        pol,
        wl,
        resolution,
        pml_um,
    )
    try:
        row = _run_one(case_name, pol, wl, resolution, pml_um, args, logger)
    except Exception as exc:
        logger.exception(
            "Task failed but will be recorded: case=%s pol=%s wl=%g res=%s pml=%g",
            case_name,
            pol,
            wl,
            resolution,
            pml_um,
        )
        row = _failed_task_row(case_name, pol, wl, resolution, pml_um, args, exc)
    return index, row


def _write_spectra_checkpoint(rows: list[dict], path: Path, logger) -> None:
    if not rows or not _mpi_am_master():
        return
    ensure_dir(path.parent)
    checkpoint = pd.DataFrame(sorted(rows, key=_task_sort_key))
    checkpoint.to_csv(path, index=False)
    logger.info("Checkpointed %d completed spectra row(s) to %s", len(checkpoint), path)


def _run_tasks(tasks: list[tuple[str, str, float, int, float]], args, paths: dict[str, Path], logger) -> list[dict]:
    workers = _resolve_workers(args, len(tasks))
    logger.info(
        "Prepared %d D12 single-wavelength tasks; workers=%d; task_timeout_s=%s (not enforced)",
        len(tasks),
        workers,
        args.task_timeout_s,
    )
    if workers == 1:
        rows = []
        for index, task in enumerate(tasks):
            case_name, pol, wl, resolution, pml_um = task
            logger.info(
                ">>> task=%d/%d case=%s pol=%s wl=%g res=%s pml=%g",
                index + 1,
                len(tasks),
                case_name,
                pol,
                wl,
                resolution,
                pml_um,
            )
            try:
                row = _run_one(case_name, pol, wl, resolution, pml_um, args, logger)
            except Exception as exc:
                logger.exception(
                    "Task failed but will be recorded: case=%s pol=%s wl=%g res=%s pml=%g",
                    case_name,
                    pol,
                    wl,
                    resolution,
                    pml_um,
                )
                row = _failed_task_row(case_name, pol, wl, resolution, pml_um, args, exc)
            rows.append(row)
            _write_spectra_checkpoint(rows, paths["spectra_checkpoint"], logger)
        return sorted(rows, key=_task_sort_key)

    payloads = [
        (index, task, args, str(paths["log"]))
        for index, task in enumerate(tasks)
    ]
    rows_by_index: dict[int, dict] = {}
    mp_context = multiprocessing.get_context("spawn")
    executor = concurrent.futures.ProcessPoolExecutor(
        max_workers=workers,
        mp_context=mp_context,
    )
    payload_iter = iter(payloads)
    pending: dict[concurrent.futures.Future, tuple[int, float]] = {}

    def submit_next() -> bool:
        try:
            payload = next(payload_iter)
        except StopIteration:
            return False
        fut = executor.submit(_run_one_task, payload)
        pending[fut] = (payload[0], time.monotonic())
        return True

    for _ in range(workers):
        submit_next()

    shutdown_started = False
    try:
        while pending:
            done, _ = concurrent.futures.wait(
                set(pending),
                timeout=1.0,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            for fut in done:
                pending.pop(fut, None)
                index, row = fut.result()
                rows_by_index[index] = row
                logger.info("Completed task %d/%d", len(rows_by_index), len(tasks))
                _write_spectra_checkpoint(
                    [rows_by_index[i] for i in sorted(rows_by_index)],
                    paths["spectra_checkpoint"],
                    logger,
                )
                submit_next()
    finally:
        if not shutdown_started:
            executor.shutdown(wait=True, cancel_futures=False)
    rows = [rows_by_index[i] for i in sorted(rows_by_index)]
    return sorted(rows, key=_task_sort_key)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["smoke", "screen", "refine"], default="smoke")
    parser.add_argument("--nk-csv", type=Path, required=True)
    parser.add_argument("--cases", nargs="*")
    parser.add_argument(
        "--film-thicknesses-um",
        nargs="*",
        default=None,
        help="Optional dynamic inner-wall film thicknesses. Generates cases such as gap10_inner_wall_film_1um.",
    )
    parser.add_argument("--wavelengths-um", nargs="*")
    parser.add_argument("--polarizations", nargs="*")
    parser.add_argument("--resolution", nargs="*", type=int)
    parser.add_argument("--pml-thickness-um", nargs="*")
    parser.add_argument("--substrate-thickness-um", type=float)
    parser.add_argument("--air-buffer-um", type=float)
    parser.add_argument("--decay-db", type=float)
    parser.add_argument("--fwidth-fraction", type=float)
    parser.add_argument(
        "--courant",
        type=float,
        default=None,
        help="Meep Courant factor. Lower values such as 0.25 can stabilize dispersive Ti runs.",
    )
    parser.add_argument("--period-um", type=float)
    parser.add_argument("--top-width-um", type=float)
    parser.add_argument("--bottom-width-um", type=float)
    parser.add_argument("--depth-um", type=float)
    parser.add_argument("--tilt-angle-deg", type=float)
    parser.add_argument(
        "--workers",
        default="auto",
        help="Number of independent single-wavelength tasks to run in parallel, or 'auto'.",
    )
    parser.add_argument(
        "--task-timeout-s",
        type=float,
        default=None,
        help="Optional timeout in seconds for submitted parallel tasks. Default: no hard timeout.",
    )
    parser.add_argument(
        "--output-tag",
        default=None,
        help="Optional filename prefix for all D12 outputs, used to run independent jobs without overwriting files.",
    )
    parser.add_argument(
        "--geometry-only",
        action="store_true",
        help="Only run geometry checks and write geometry/report artifacts; skip Meep solves.",
    )
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    is_master = _mpi_am_master()
    paths = _paths(args.output_tag)
    for path in paths.values():
        ensure_dir(path.parent)
    logger = _setup_logger(paths["log"])
    try:
        _configure_args(args)
        _apply_refine_gate(args, paths, logger)
        logger.info("Starting D12 mode=%s", args.mode)
        geometry_checks = _geometry_checks(args)
        failed_geometry = geometry_checks[geometry_checks["status"] != "PASS"]
        if not failed_geometry.empty:
            if is_master:
                _generate_geometry_plots(paths, args)
                empty_spectra = pd.DataFrame()
                empty_metrics = pd.DataFrame()
                res_check = pd.DataFrame([{"resolution_check_status": "NOT_RUN_GEOMETRY_FAILED"}])
                pml_check = pd.DataFrame([{"pml_check_status": "NOT_RUN_GEOMETRY_FAILED"}])
                empty_spectra.to_csv(paths["spectra"], index=False)
                empty_metrics.to_csv(paths["metrics"], index=False)
                res_check.to_csv(paths["resolution"], index=False)
                pml_check.to_csv(paths["pml"], index=False)
                _write_report(
                    paths,
                    empty_spectra,
                    empty_metrics,
                    res_check,
                    pml_check,
                    geometry_checks,
                    {
                        "field_snapshot_status": "NOT_RUN_GEOMETRY_FAILED",
                        "field_snapshot_note": "The tilted groove leaves the periodic unit cell.",
                    },
                    args,
                )
                print("FAIL")
            logger.error("D12 geometry checks failed:\n%s", failed_geometry.to_string(index=False))
            return
        if args.geometry_only:
            if is_master:
                _generate_geometry_plots(paths, args)
                empty_spectra = pd.DataFrame()
                empty_metrics = pd.DataFrame()
                res_check = pd.DataFrame([{"resolution_check_status": "NOT_RUN_GEOMETRY_ONLY"}])
                pml_check = pd.DataFrame([{"pml_check_status": "NOT_RUN_GEOMETRY_ONLY"}])
                empty_spectra.to_csv(paths["spectra"], index=False)
                empty_metrics.to_csv(paths["metrics"], index=False)
                res_check.to_csv(paths["resolution"], index=False)
                pml_check.to_csv(paths["pml"], index=False)
                _write_report(
                    paths,
                    empty_spectra,
                    empty_metrics,
                    res_check,
                    pml_check,
                    geometry_checks,
                    {
                        "field_snapshot_status": "NOT_RUN_GEOMETRY_ONLY",
                        "field_snapshot_note": "Geometry-only run skipped Meep field snapshots.",
                    },
                    args,
                )
                print("CODE_PASS")
            logger.info("D12 geometry-only run passed.")
            return
        if is_master:
            _generate_geometry_plots(paths, args)

        tasks = _build_tasks(args)
        rows = _run_tasks(tasks, args, paths, logger)

        spectra = pd.DataFrame(rows)
        metrics = _compute_metrics(spectra, args)
        res_check = _resolution_check(metrics)
        pml_check = _pml_check(metrics)
        spectra_out = _append_proxy_spectra(spectra, args)

        if is_master:
            spectra_out.to_csv(paths["spectra"], index=False)
            metrics.to_csv(paths["metrics"], index=False)
            res_check.to_csv(paths["resolution"], index=False)
            pml_check.to_csv(paths["pml"], index=False)
            _make_plots(spectra, metrics, paths, args)
            field_info = _run_best_field_snapshot(spectra, metrics, paths, args, logger)
        else:
            field_info = {
                "field_snapshot_status": "SKIPPED_ON_NON_MASTER_RANK",
                "field_snapshot_note": "MPI non-master rank does not write output artifacts.",
            }
        if is_master:
            _write_report(paths, spectra, metrics, res_check, pml_check, geometry_checks, field_info, args)
            logger.info("Wrote %s", paths["spectra"])
            logger.info("Wrote %s", paths["metrics"])
            logger.info("Wrote %s", paths["report"])
    except Exception as exc:
        logger.exception("D12 failed")
        if is_master:
            _write_failure_report(paths, exc, args)
        raise


if __name__ == "__main__":
    main()
