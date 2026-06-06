"""D13 wrapped 20-degree gap10 groove with 2 um inner-wall film.

This diagnostics_v2 script keeps the physical D12 gap10 pitch:

    P=50 um, top/bottom groove width=40 um, depth=30 um, tilt=20 deg.

The 20-degree bottom shift makes one slanted groove cross the right unit-cell
boundary.  D13 represents the same periodic structure by clipping the groove
polygon to the cell and wrapping the overflow piece back to the left boundary.
The top opening remains 40 um, so the periodic top Ti gap remains 10 um.

Outputs include normal-incidence total emissivity proxy and true angular
absorptance A(lambda, theta), which is the directional emissivity proxy under
Kirchhoff reciprocity for the same direction and polarization.
"""

from __future__ import annotations

import argparse
import cmath
import concurrent.futures
import logging
import math
import multiprocessing
import os
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
from matplotlib.patches import Polygon as MplPolygon
import numpy as np
import pandas as pd

from src.io_utils import ensure_dir, project_path, save_figure
from src.materials import (
    TI_RAKIC_VALID_LAMBDA_UM,
    get_measured_lossy_wall_film_medium_single_wavelength,
    get_ti_medium,
)
from src.postprocess import opaque_substrate_transmission_check, wavelength_integrated_average
from src.simulation import run_periodic_2d_metal_single_wavelength


TAG = "D13_tilt20_wrapped_gap10_w40_h30_film2um"
PROXY = "unpolarized_2D_proxy"

BAND_STRICT_LO = 8.1014
BAND_STRICT_HI = TI_RAKIC_VALID_LAMBDA_UM[1]

DEFAULT_WAVELENGTHS_UM = [8.1014, 9.0016, 10.009, 10.985, 12.345]
DEFAULT_ANGLES_DEG = list(np.arange(-70.0, 70.0 + 1e-9, 10.0))
DEFAULT_POLARIZATIONS = ["Ez", "Hz"]

DEFAULT_PERIOD_UM = 50.0
DEFAULT_TOP_WIDTH_UM = 40.0
DEFAULT_BOTTOM_WIDTH_UM = 40.0
DEFAULT_DEPTH_UM = 30.0
DEFAULT_TILT_DEG = 20.0
DEFAULT_FILM_UM = 2.0
DEFAULT_SUBSTRATE_UM = 45.0
DEFAULT_AIR_UM = 15.0
DEFAULT_PML_UM = 4.0
DEFAULT_RESOLUTION = 12
DEFAULT_COURANT = 0.1
DEFAULT_RETRY_COURANT = 0.05
DEFAULT_DECAY_DB = 20.0
DEFAULT_FWIDTH_FRACTION = 0.08


def _paths(output_tag: str = TAG) -> dict[str, Path]:
    fig_dir = project_path("results", "diagnostics_v2", "figures")
    tab_dir = project_path("results", "diagnostics_v2", "tables")
    rep_dir = project_path("results", "diagnostics_v2", "reports")
    log_dir = project_path("logs", "diagnostics_v2")
    return {
        "normal_spectra": tab_dir / f"{output_tag}_normal_spectra.csv",
        "normal_metrics": tab_dir / f"{output_tag}_normal_metrics.csv",
        "angle_spectra": tab_dir / f"{output_tag}_angle_spectra.csv",
        "angle_metrics": tab_dir / f"{output_tag}_angle_metrics.csv",
        "normal_checkpoint": tab_dir / f"{output_tag}_normal_spectra_checkpoint.csv",
        "angle_checkpoint": tab_dir / f"{output_tag}_angle_spectra_checkpoint.csv",
        "report": rep_dir / f"{output_tag}_report.md",
        "log": log_dir / f"{output_tag}.log",
        "geometry": fig_dir / f"{output_tag}_geometry.png",
        "normal_spectra_fig": fig_dir / f"{output_tag}_normal_spectra.png",
        "angle_heatmap_ez": fig_dir / f"{output_tag}_angle_heatmap_Ez.png",
        "angle_heatmap_hz": fig_dir / f"{output_tag}_angle_heatmap_Hz.png",
        "angle_proxy": fig_dir / f"{output_tag}_angle_proxy_band_average.png",
    }


def _setup_logger(path: Path, name: str = "D13_tilt20_wrapped") -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger(name)
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


def _setup_worker_logger(path: Path) -> logging.Logger:
    return _setup_logger(path, name=f"D13_tilt20_wrapped.worker.{os.getpid()}")


def _parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return list(default)
    out: list[float] = []
    for item in values:
        for token in item.replace(",", " ").split():
            out.append(float(token))
    return out


def _parse_str_list(values: list[str] | None, default: list[str]) -> list[str]:
    if not values:
        return list(default)
    out: list[str] = []
    for item in values:
        out.extend(item.replace(",", " ").split())
    return out


def _source_component(pol: str):
    import meep as mp

    if pol == "Ez":
        return mp.Ez
    if pol == "Hz":
        return mp.Hz
    raise ValueError(f"Unknown polarization: {pol}")


def _slanted_vertices(
    *,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    y_surface: float,
) -> list[tuple[float, float]]:
    dx = depth_um * math.tan(math.radians(tilt_angle_deg))
    return [
        (-top_width_um / 2.0, y_surface),
        (+top_width_um / 2.0, y_surface),
        (+bottom_width_um / 2.0 + dx, y_surface - depth_um),
        (-bottom_width_um / 2.0 + dx, y_surface - depth_um),
    ]


def _clip_x(poly: list[tuple[float, float]], *, x_min: float, x_max: float) -> list[tuple[float, float]]:
    def clip_side(points: list[tuple[float, float]], inside, intersect):
        if not points:
            return []
        out = []
        prev = points[-1]
        prev_inside = inside(prev)
        for cur in points:
            cur_inside = inside(cur)
            if cur_inside:
                if not prev_inside:
                    out.append(intersect(prev, cur))
                out.append(cur)
            elif prev_inside:
                out.append(intersect(prev, cur))
            prev, prev_inside = cur, cur_inside
        return out

    def left_intersect(a, b):
        ax, ay = a
        bx, by = b
        t = (x_min - ax) / (bx - ax)
        return (x_min, ay + t * (by - ay))

    def right_intersect(a, b):
        ax, ay = a
        bx, by = b
        t = (x_max - ax) / (bx - ax)
        return (x_max, ay + t * (by - ay))

    out = clip_side(poly, lambda p: p[0] >= x_min - 1e-10, left_intersect)
    out = clip_side(out, lambda p: p[0] <= x_max + 1e-10, right_intersect)
    return _clean_polygon(out)


def _clean_polygon(poly: list[tuple[float, float]], tol: float = 1e-9) -> list[tuple[float, float]]:
    if not poly:
        return []
    cleaned: list[tuple[float, float]] = []
    for x, y in poly:
        pt = (0.0 if abs(x) < tol else float(x), 0.0 if abs(y) < tol else float(y))
        if not cleaned or abs(pt[0] - cleaned[-1][0]) > tol or abs(pt[1] - cleaned[-1][1]) > tol:
            cleaned.append(pt)
    if len(cleaned) > 1:
        first, last = cleaned[0], cleaned[-1]
        if abs(first[0] - last[0]) <= tol and abs(first[1] - last[1]) <= tol:
            cleaned.pop()
    if _polygon_area(cleaned) < 1e-8:
        return []
    return cleaned


def _polygon_area(poly: list[tuple[float, float]]) -> float:
    if len(poly) < 3:
        return 0.0
    area = 0.0
    for (x0, y0), (x1, y1) in zip(poly, poly[1:] + poly[:1]):
        area += x0 * y1 - x1 * y0
    return abs(area) * 0.5


def _wrapped_polygons(poly: list[tuple[float, float]], period_um: float) -> list[list[tuple[float, float]]]:
    half = period_um / 2.0
    pieces: list[list[tuple[float, float]]] = []
    for shift_i in range(-2, 3):
        shifted = [(x + shift_i * period_um, y) for x, y in poly]
        clipped = _clip_x(shifted, x_min=-half, x_max=half)
        if clipped:
            pieces.append(clipped)
    return pieces


def _wrapped_inner_wall_film_geometry(
    *,
    y_surface: float,
    period_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    film_thickness_um: float,
    substrate_thickness_um: float,
    medium_substrate,
    medium_film,
    medium_groove,
) -> list:
    import meep as mp

    if film_thickness_um <= 0:
        raise ValueError("D13 expects a positive inner-wall film thickness.")
    if 2.0 * film_thickness_um >= min(top_width_um, bottom_width_um):
        raise ValueError("Film thickness collapses the inner air core.")
    if film_thickness_um >= depth_um:
        raise ValueError("Film thickness must be smaller than groove depth.")
    if depth_um >= substrate_thickness_um:
        raise ValueError("Groove depth must be smaller than substrate thickness.")

    outer = _slanted_vertices(
        top_width_um=top_width_um,
        bottom_width_um=bottom_width_um,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    inner = _slanted_vertices(
        top_width_um=top_width_um - 2.0 * film_thickness_um,
        bottom_width_um=bottom_width_um - 2.0 * film_thickness_um,
        depth_um=depth_um - film_thickness_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    outer_pieces = _wrapped_polygons(outer, period_um)
    inner_pieces = _wrapped_polygons(inner, period_um)
    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_um, substrate_thickness_um, mp.inf),
        )
    ]
    for piece in outer_pieces:
        geom.append(
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in piece],
                height=mp.inf,
                axis=mp.Vector3(0, 0, 1),
                material=medium_film,
            )
        )
    for piece in inner_pieces:
        geom.append(
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in piece],
                height=mp.inf,
                axis=mp.Vector3(0, 0, 1),
                material=medium_groove,
            )
        )
    return geom


def _geometry_factory(film_medium, args):
    import meep as mp

    ti = get_ti_medium()
    air = mp.Medium(epsilon=1.0)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        return _wrapped_inner_wall_film_geometry(
            y_surface=y_surface_um,
            period_um=args.period_um,
            top_width_um=args.top_width_um,
            bottom_width_um=args.bottom_width_um,
            depth_um=args.depth_um,
            tilt_angle_deg=args.tilt_angle_deg,
            film_thickness_um=args.film_thickness_um,
            substrate_thickness_um=substrate_thickness_um,
            medium_substrate=ti,
            medium_film=film_medium,
            medium_groove=air,
        )

    return factory


def _layout(args) -> dict[str, float]:
    bottom_buffer_um = args.pml_thickness_um
    cell_y = (
        2.0 * args.pml_thickness_um
        + args.air_buffer_um
        + args.substrate_thickness_um
        + bottom_buffer_um
    )
    y_top_edge = 0.5 * cell_y
    y_bottom_edge = -0.5 * cell_y
    y_top_pml_inner = y_top_edge - args.pml_thickness_um
    y_bottom_pml_inner = y_bottom_edge + args.pml_thickness_um
    y_surface = y_top_pml_inner - args.air_buffer_um
    y_substrate_bottom = y_surface - args.substrate_thickness_um
    y_src = y_top_pml_inner - 0.25 * args.air_buffer_um
    y_refl = y_surface + 0.5 * args.air_buffer_um
    y_trans = 0.5 * (y_substrate_bottom + y_bottom_pml_inner)
    return {
        "cell_y": cell_y,
        "y_surface": y_surface,
        "y_src": y_src,
        "y_refl": y_refl,
        "y_trans": y_trans,
    }


def _film_medium_for_wavelength(wl: float, args):
    medium, meta = get_measured_lossy_wall_film_medium_single_wavelength(
        wl,
        args.nk_csv,
        allow_extrapolation=False,
    )
    return medium, {
        "material_name_film": meta["material_name"],
        "n_film": meta["n"],
        "k_film": meta["k"],
        "epsilon_real_film": meta["epsilon_real"],
        "epsilon_imag_film": meta["epsilon_imag"],
        "D_conductivity_film": meta["D_conductivity"],
        "film_interpolation_flag": meta["interpolation_flag"],
        "film_model_mode": meta["model_mode"],
        "film_data_range_min_um": meta["data_range_um"][0],
        "film_data_range_max_um": meta["data_range_um"][1],
        "film_warning": meta["warning"],
    }


def _transmission_flag(T: float) -> str:
    flag = str(opaque_substrate_transmission_check(np.array([T]))["transmission_quality_flag"][0])
    if flag == "NUMERICAL_PASS":
        return "NUMERICAL_PASS"
    if flag == "WARNING":
        return "WARNING_EXCLUDED_FROM_STRICT_METRIC"
    return "FAIL_EXCLUDED_FROM_ALL_METRICS"


def _finite_valid(row: dict) -> bool:
    keys = ["R", "T", "A", "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"]
    return bool(all(np.isfinite(row.get(k, np.nan)) for k in keys))


def _decorate_common_row(row: dict, *, wl: float, pol: str, theta: float | None, args, film_meta: dict) -> dict:
    ti_flag = "VALID" if wl <= BAND_STRICT_HI else "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING"
    finite = _finite_valid(row)
    trans_flag = _transmission_flag(float(row["T"])) if finite else "FAIL_EXCLUDED_FROM_ALL_METRICS"
    valid = finite and trans_flag == "NUMERICAL_PASS" and ti_flag == "VALID"
    common = {
        "case_name": "tilt20_wrapped_gap10_film2um",
        "film_thickness_um": args.film_thickness_um,
        "coating_mode": "sidewalls_and_bottom",
        "polarization": pol,
        "wavelength_um": wl,
        "theta_deg": np.nan if theta is None else theta,
        "included_in_strict_metric": bool(valid and BAND_STRICT_LO <= wl <= BAND_STRICT_HI),
        "finite_values_flag": "VALID" if finite else "NAN_OR_INF",
        "material_validity_flag_ti": ti_flag,
        "transmission_quality_flag": trans_flag,
        "period_um": args.period_um,
        "top_width_um": args.top_width_um,
        "bottom_width_um": args.bottom_width_um,
        "depth_um": args.depth_um,
        "tilt_angle_deg": args.tilt_angle_deg,
        "bottom_shift_um": args.depth_um * math.tan(math.radians(args.tilt_angle_deg)),
        "film_grid_points": args.film_thickness_um * args.resolution,
        "normalization_note": "D00 convention: R=refl/abs(input), T=-trans/abs(input), A=1-R-T.",
        "geometry_note": "Wrapped periodic polygon pieces preserve P=50 and gap=10.",
    }
    common.update(row)
    common.update(film_meta)
    return common


def _run_normal_once(pol: str, wl: float, courant: float, args, logger) -> dict:
    film_medium, film_meta = _film_medium_for_wavelength(wl, args)
    result = run_periodic_2d_metal_single_wavelength(
        geometry_factory=_geometry_factory(film_medium, args),
        period_um=args.period_um,
        wavelength_um=wl,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        decay_db=args.decay_db,
        source_component=pol,
        fwidth_fraction=args.fwidth_fraction,
        courant=courant,
        solver_version="D13_wrapped_gap10_normal_v1",
        source_mode="single_wavelength_normal_incidence",
        logger=logger,
    )
    row = {
        "R": result["R"],
        "T": result["T"],
        "A": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "resolution": result["resolution"],
        "pml_thickness_um": result["pml_thickness_um"],
        "substrate_thickness_um": result["substrate_thickness_um"],
        "air_buffer_um": result["air_buffer_um"],
        "decay_db": result["decay_db"],
        "fwidth_fraction": result["fwidth_fraction"],
        "courant": courant,
        "solver_version": result["solver_version"],
        "source_mode": result["source_mode"],
        "walltime_s": result["walltime_s"],
        "retry_level": "primary" if courant == args.courant else "retry_courant",
    }
    return _decorate_common_row(row, wl=wl, pol=pol, theta=None, args=args, film_meta=film_meta)


def _run_angle_once(pol: str, wl: float, theta: float, courant: float, args, logger) -> dict:
    import meep as mp

    film_medium, film_meta = _film_medium_for_wavelength(wl, args)
    ti = get_ti_medium()
    layout = _layout(args)
    fcen = 1.0 / wl
    fwidth = args.fwidth_fraction * fcen
    kx = fcen * math.sin(math.radians(theta))
    src_c = _source_component(pol)

    def amp_func(r):
        return cmath.exp(1j * 2.0 * math.pi * kx * r.x)

    cell = mp.Vector3(args.period_um, layout["cell_y"], 0)
    pml_layers = [mp.PML(thickness=args.pml_thickness_um, direction=mp.Y)]
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=fwidth, is_integrated=True),
            component=src_c,
            center=mp.Vector3(0, layout["y_src"], 0),
            size=mp.Vector3(args.period_um, 0, 0),
            amp_func=amp_func,
        )
    ]
    k_point = mp.Vector3(kx, 0, 0)
    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        Courant=courant,
        k_point=k_point,
        geometry=[],
        dimensions=2,
    )
    refl_ref = sim_ref.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(center=mp.Vector3(0, layout["y_refl"], 0), size=mp.Vector3(args.period_um, 0, 0)),
    )
    t0 = time.time()
    sim_ref.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, layout["y_refl"], 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    input_flux_raw = float(np.array(mp.get_fluxes(refl_ref))[0])
    ref_data = sim_ref.get_flux_data(refl_ref)
    t_ref = time.time() - t0
    if abs(input_flux_raw) == 0:
        raise RuntimeError("Reference input flux is zero.")

    geometry = _wrapped_inner_wall_film_geometry(
        y_surface=layout["y_surface"],
        period_um=args.period_um,
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        film_thickness_um=args.film_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        medium_substrate=ti,
        medium_film=film_medium,
        medium_groove=mp.Medium(epsilon=1.0),
    )
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        Courant=courant,
        k_point=k_point,
        geometry=geometry,
        dimensions=2,
    )
    refl = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(center=mp.Vector3(0, layout["y_refl"], 0), size=mp.Vector3(args.period_um, 0, 0)),
    )
    trans = sim.add_flux(
        fcen,
        0,
        1,
        mp.FluxRegion(center=mp.Vector3(0, layout["y_trans"], 0), size=mp.Vector3(args.period_um, 0, 0)),
    )
    sim.load_minus_flux_data(refl, ref_data)
    t0 = time.time()
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, layout["y_refl"], 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    reflection_flux_raw = float(np.array(mp.get_fluxes(refl))[0])
    transmission_flux_raw = float(np.array(mp.get_fluxes(trans))[0])
    t_struct = time.time() - t0

    R = reflection_flux_raw / abs(input_flux_raw)
    raw_T = transmission_flux_raw / abs(input_flux_raw)
    T = -transmission_flux_raw / abs(input_flux_raw)
    A = 1.0 - R - T
    row = {
        "R": R,
        "T": T,
        "A": A,
        "input_flux_raw": input_flux_raw,
        "reflection_flux_raw": reflection_flux_raw,
        "transmission_flux_raw": transmission_flux_raw,
        "raw_transmittance": raw_T,
        "signed_transmittance": T,
        "frequency_meep": fcen,
        "kx_meep": kx,
        "resolution": args.resolution,
        "pml_thickness_um": args.pml_thickness_um,
        "substrate_thickness_um": args.substrate_thickness_um,
        "air_buffer_um": args.air_buffer_um,
        "decay_db": args.decay_db,
        "fwidth_fraction": args.fwidth_fraction,
        "courant": courant,
        "solver_version": "D13_wrapped_gap10_angle_v1",
        "source_mode": "single_wavelength_oblique_incidence",
        "walltime_s": t_ref + t_struct,
        "retry_level": "primary" if courant == args.courant else "retry_courant",
        "kirchhoff_note": "A(lambda,theta) is the directional emissivity proxy under reciprocity.",
    }
    return _decorate_common_row(row, wl=wl, pol=pol, theta=theta, args=args, film_meta=film_meta)


def _failed_row(kind: str, pol: str, wl: float, theta: float | None, args, exc: Exception) -> dict:
    return {
        "case_name": "tilt20_wrapped_gap10_film2um",
        "film_thickness_um": args.film_thickness_um,
        "coating_mode": "sidewalls_and_bottom",
        "polarization": pol,
        "wavelength_um": wl,
        "theta_deg": np.nan if theta is None else theta,
        "R": np.nan,
        "T": np.nan,
        "A": np.nan,
        "input_flux_raw": np.nan,
        "reflection_flux_raw": np.nan,
        "transmission_flux_raw": np.nan,
        "raw_transmittance": np.nan,
        "signed_transmittance": np.nan,
        "included_in_strict_metric": False,
        "finite_values_flag": "TASK_FAILED",
        "material_validity_flag_ti": "VALID" if wl <= BAND_STRICT_HI else "TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING",
        "transmission_quality_flag": "FAIL_EXCLUDED_FROM_ALL_METRICS",
        "period_um": args.period_um,
        "top_width_um": args.top_width_um,
        "bottom_width_um": args.bottom_width_um,
        "depth_um": args.depth_um,
        "tilt_angle_deg": args.tilt_angle_deg,
        "bottom_shift_um": args.depth_um * math.tan(math.radians(args.tilt_angle_deg)),
        "resolution": args.resolution,
        "pml_thickness_um": args.pml_thickness_um,
        "substrate_thickness_um": args.substrate_thickness_um,
        "air_buffer_um": args.air_buffer_um,
        "decay_db": args.decay_db,
        "fwidth_fraction": args.fwidth_fraction,
        "courant": args.retry_courant,
        "solver_version": f"D13_wrapped_gap10_{kind}_v1",
        "source_mode": "task_failed",
        "walltime_s": np.nan,
        "retry_level": "failed_after_retry",
        "task_error_type": type(exc).__name__,
        "task_error": str(exc),
    }


def _run_with_retry(kind: str, pol: str, wl: float, theta: float | None, args, logger) -> dict:
    for courant in [args.courant, args.retry_courant]:
        try:
            if kind == "normal":
                row = _run_normal_once(pol, wl, courant, args, logger)
            else:
                row = _run_angle_once(pol, wl, float(theta), courant, args, logger)
            if _finite_valid(row):
                return row
            raise RuntimeError("R/T/A or raw flux contains NaN/Inf.")
        except Exception as exc:
            last_exc = exc
            logger.exception(
                "D13 %s task failed at courant=%g: pol=%s wl=%g theta=%s",
                kind,
                courant,
                pol,
                wl,
                theta,
            )
    return _failed_row(kind, pol, wl, theta, args, last_exc)


def _task_worker(payload):
    index, kind, pol, wl, theta, args, log_path = payload
    logger = _setup_worker_logger(Path(log_path))
    logger.info(">>> D13 task %d kind=%s pol=%s wl=%g theta=%s", index, kind, pol, wl, theta)
    row = _run_with_retry(kind, pol, wl, theta, args, logger)
    return index, row


def _resolve_workers(args, task_count: int) -> int:
    if task_count <= 0:
        return 1
    requested = str(args.workers).strip().lower()
    cpu_count = os.cpu_count() or 1
    if requested == "auto":
        return max(1, min(12, cpu_count, task_count))
    workers = int(requested)
    if workers < 1:
        raise ValueError("--workers must be >= 1.")
    return min(workers, task_count)


def _run_tasks(
    kind: str,
    tasks: list[tuple[str, float, float | None]],
    args,
    path: Path,
    checkpoint: Path,
    logger,
    *,
    checkpoint_seed: pd.DataFrame | None = None,
):
    workers = _resolve_workers(args, len(tasks))
    logger.info("Prepared %d D13 %s tasks; workers=%d", len(tasks), kind, workers)
    payloads = [(i, kind, pol, wl, theta, args, str(path)) for i, (pol, wl, theta) in enumerate(tasks)]
    rows_by_index: dict[int, dict] = {}
    if workers == 1:
        for payload in payloads:
            index, row = _task_worker(payload)
            rows_by_index[index] = row
            _write_checkpoint(rows_by_index, checkpoint, logger, checkpoint_seed=checkpoint_seed)
    else:
        mp_context = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=mp_context) as executor:
            future_map = {executor.submit(_task_worker, payload): payload[0] for payload in payloads}
            for fut in concurrent.futures.as_completed(future_map):
                index, row = fut.result()
                rows_by_index[index] = row
                logger.info("Completed %s task %d/%d", kind, len(rows_by_index), len(tasks))
                _write_checkpoint(rows_by_index, checkpoint, logger, checkpoint_seed=checkpoint_seed)
    return [rows_by_index[i] for i in sorted(rows_by_index)]


def _task_key(kind: str, pol: str, wl: float, theta: float | None) -> tuple:
    if kind == "normal":
        return (pol, round(float(wl), 9))
    return (pol, round(float(wl), 9), round(float(theta), 9))


def _row_task_key(kind: str, row: pd.Series) -> tuple:
    theta = None if kind == "normal" else float(row["theta_deg"])
    return _task_key(kind, str(row["polarization"]), float(row["wavelength_um"]), theta)


def _load_resume_rows(kind: str, checkpoint: Path, logger) -> tuple[pd.DataFrame, set[tuple]]:
    if not checkpoint.exists():
        return pd.DataFrame(), set()
    try:
        existing = pd.read_csv(checkpoint)
    except Exception as exc:
        logger.warning("Could not read %s resume checkpoint %s: %s", kind, checkpoint, exc)
        return pd.DataFrame(), set()
    if existing.empty:
        return existing, set()
    required = {"polarization", "wavelength_um", "finite_values_flag"}
    if kind == "angle":
        required.add("theta_deg")
    if not required.issubset(existing.columns):
        logger.warning("%s checkpoint missing required columns for resume; ignoring %s", kind, checkpoint)
        return pd.DataFrame(), set()
    valid = existing[existing["finite_values_flag"].eq("VALID")].copy()
    keys = {_row_task_key(kind, row) for _, row in valid.iterrows()}
    logger.info("Loaded %d valid %s checkpoint row(s) from %s", len(keys), kind, checkpoint)
    return valid, keys


def _filter_resume_tasks(
    kind: str,
    tasks: list[tuple[str, float, float | None]],
    completed_keys: set[tuple],
    logger,
) -> list[tuple[str, float, float | None]]:
    if not completed_keys:
        return tasks
    remaining = [
        task for task in tasks
        if _task_key(kind, task[0], task[1], task[2]) not in completed_keys
    ]
    logger.info(
        "Resume filtered %s tasks: %d completed, %d remaining",
        kind,
        len(tasks) - len(remaining),
        len(remaining),
    )
    return remaining


def _write_checkpoint(
    rows_by_index: dict[int, dict],
    path: Path,
    logger,
    *,
    checkpoint_seed: pd.DataFrame | None = None,
) -> None:
    ensure_dir(path.parent)
    rows = [rows_by_index[i] for i in sorted(rows_by_index)]
    frame = pd.DataFrame(rows)
    if checkpoint_seed is not None and not checkpoint_seed.empty:
        frame = pd.concat([checkpoint_seed, frame], ignore_index=True)
    frame.to_csv(path, index=False)
    logger.info("Checkpointed %d row(s) to %s", len(frame), path)


def _add_proxy_rows(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    group_cols = ["wavelength_um"]
    if "theta_deg" in df.columns and df["theta_deg"].notna().any():
        group_cols.append("theta_deg")
    rows = []
    for key, group in df.groupby(group_cols, dropna=False, sort=False):
        by_pol = group.set_index("polarization")
        if "Ez" not in by_pol.index or "Hz" not in by_pol.index:
            continue
        ez = by_pol.loc["Ez"]
        hz = by_pol.loc["Hz"]
        row = ez.to_dict()
        row["polarization"] = PROXY
        for col in ["R", "T", "A", "raw_transmittance", "signed_transmittance"]:
            row[col] = float(np.nanmean([ez[col], hz[col]]))
        row["input_flux_raw"] = np.nan
        row["reflection_flux_raw"] = np.nan
        row["transmission_flux_raw"] = np.nan
        row["included_in_strict_metric"] = bool(ez["included_in_strict_metric"] and hz["included_in_strict_metric"])
        row["finite_values_flag"] = "VALID" if row["included_in_strict_metric"] else "PROXY_FROM_INCOMPLETE_POLARIZATIONS"
        row["source_mode"] = "proxy_average"
        row["solver_version"] = "D13_proxy_average_v1"
        rows.append(row)
    if not rows:
        return df
    return pd.concat([df, pd.DataFrame(rows)], ignore_index=True)


def _normal_metrics(normal: pd.DataFrame, args) -> pd.DataFrame:
    rows = []
    for pol, group in normal.groupby("polarization", sort=False):
        mask = group["included_in_strict_metric"].to_numpy(dtype=bool)
        mean_a = wavelength_integrated_average(
            group["A"].to_numpy(dtype=float),
            group["wavelength_um"].to_numpy(dtype=float),
            BAND_STRICT_LO,
            BAND_STRICT_HI,
            valid_mask=mask,
        )
        status = wavelength_integrated_average.last_status
        valid_count = int(mask.sum())
        if valid_count == 1 and len(args.wavelengths_um) == 1:
            mean_a = float(group.loc[mask, "A"].iloc[0])
            status = "valid_single_point_smoke"
        if valid_count < len(args.wavelengths_um) or status != "valid":
            if status != "valid_single_point_smoke":
                mean_a = np.nan
        valid = group[mask]
        peak_idx = valid["A"].idxmax() if not valid.empty else None
        rows.append(
            {
                "polarization": pol,
                "film_thickness_um": args.film_thickness_um,
                "mean_A_8p1014_12p398_strict": mean_a,
                "peak_A_strict": float(valid.loc[peak_idx, "A"]) if peak_idx is not None else np.nan,
                "peak_wavelength_strict_um": float(valid.loc[peak_idx, "wavelength_um"]) if peak_idx is not None else np.nan,
                "valid_point_count_strict": valid_count,
                "total_planned_point_count_strict": len(args.wavelengths_um),
                "max_abs_T": float(group["T"].abs().max()),
                "numerical_status": "NUMERICAL_SCREENING" if np.isfinite(mean_a) else "INCOMPLETE_NUMERICAL_COVERAGE",
            }
        )
    return pd.DataFrame(rows)


def _angle_metrics(angle: pd.DataFrame, args) -> pd.DataFrame:
    rows = []
    for pol, group in angle.groupby("polarization", sort=False):
        band_rows = []
        for theta, theta_df in group.groupby("theta_deg", sort=True):
            mask = theta_df["included_in_strict_metric"].to_numpy(dtype=bool)
            mean_a = wavelength_integrated_average(
                theta_df["A"].to_numpy(dtype=float),
                theta_df["wavelength_um"].to_numpy(dtype=float),
                BAND_STRICT_LO,
                BAND_STRICT_HI,
                valid_mask=mask,
            )
            status = wavelength_integrated_average.last_status
            if int(mask.sum()) == 1 and len(args.wavelengths_um) == 1:
                mean_a = float(theta_df.loc[mask, "A"].iloc[0])
                status = "valid_single_point_smoke"
            if int(mask.sum()) < len(args.wavelengths_um) or status not in {"valid", "valid_single_point_smoke"}:
                mean_a = np.nan
            band_rows.append((float(theta), float(mean_a), int(mask.sum())))
        band = pd.DataFrame(band_rows, columns=["theta_deg", "band_mean_A_strict", "valid_point_count_strict"])
        if band["band_mean_A_strict"].notna().any():
            best = band.loc[band["band_mean_A_strict"].idxmax()]
        else:
            best = pd.Series({"theta_deg": np.nan, "band_mean_A_strict": np.nan})
        plus = _interp_angle_metric(band, args.theta0_deg)
        minus = _interp_angle_metric(band, -args.theta0_deg)
        rows.append(
            {
                "metric_scope": "band_average_by_angle",
                "polarization": pol,
                "best_theta_deg": float(best["theta_deg"]),
                "best_band_mean_A_strict": float(best["band_mean_A_strict"]),
                "theta0_deg": args.theta0_deg,
                "A_plus_theta0": plus,
                "A_minus_theta0": minus,
                "asymmetry_plus_minus_theta0": plus - minus if np.all(np.isfinite([plus, minus])) else np.nan,
                "ratio_plus_minus_theta0": plus / minus if np.isfinite(minus) and abs(minus) > 1e-12 else np.nan,
                "mean_over_angles": float(band["band_mean_A_strict"].mean()),
                "valid_angle_count": int(band["band_mean_A_strict"].notna().sum()),
                "total_angle_count": len(args.angles_deg),
                "max_abs_T": float(group["T"].abs().max()),
            }
        )
        for _, row in band.iterrows():
            rows.append(
                {
                    "metric_scope": "per_angle",
                    "polarization": pol,
                    "theta_deg": row["theta_deg"],
                    "band_mean_A_strict": row["band_mean_A_strict"],
                    "valid_point_count_strict": row["valid_point_count_strict"],
                    "total_planned_point_count_strict": len(args.wavelengths_um),
                }
            )
    return pd.DataFrame(rows)


def _interp_angle_metric(band: pd.DataFrame, theta: float) -> float:
    sub = band.dropna(subset=["band_mean_A_strict"]).sort_values("theta_deg")
    if sub.empty or theta < sub["theta_deg"].min() or theta > sub["theta_deg"].max():
        return np.nan
    return float(np.interp(theta, sub["theta_deg"], sub["band_mean_A_strict"]))


def _geometry_checks(args) -> pd.DataFrame:
    y = 0.0
    outer = _slanted_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=y,
    )
    inner = _slanted_vertices(
        top_width_um=args.top_width_um - 2.0 * args.film_thickness_um,
        bottom_width_um=args.bottom_width_um - 2.0 * args.film_thickness_um,
        depth_um=args.depth_um - args.film_thickness_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=y,
    )
    half = args.period_um / 2.0
    outer_pieces = _wrapped_polygons(outer, args.period_um)
    inner_pieces = _wrapped_polygons(inner, args.period_um)
    dx = args.depth_um * math.tan(math.radians(args.tilt_angle_deg))
    rows = [
        {
            "check_name": "raw_bottom_shift_requires_wrap",
            "status": "PASS" if dx + args.bottom_width_um / 2.0 > half else "FAIL",
            "details": f"bottom_shift={dx:.6g} um; raw bottom right={dx + args.bottom_width_um / 2.0:.6g}; half_period={half:g}",
        },
        {
            "check_name": "top_gap_is_10um",
            "status": "PASS" if abs(args.period_um - args.top_width_um - 10.0) < 1e-9 else "FAIL",
            "details": f"period={args.period_um:g}, top_width={args.top_width_um:g}, gap={args.period_um - args.top_width_um:g}",
        },
        {
            "check_name": "inner_air_core_not_collapsed",
            "status": "PASS"
            if args.top_width_um - 2 * args.film_thickness_um > 0
            and args.bottom_width_um - 2 * args.film_thickness_um > 0
            and args.depth_um - args.film_thickness_um > 0
            else "FAIL",
            "details": (
                f"inner_top={args.top_width_um - 2 * args.film_thickness_um:g}, "
                f"inner_bottom={args.bottom_width_um - 2 * args.film_thickness_um:g}, "
                f"inner_depth={args.depth_um - args.film_thickness_um:g}"
            ),
        },
        {
            "check_name": "wrapped_outer_pieces_inside_cell",
            "status": _pieces_inside_status(outer_pieces, half),
            "details": f"outer_piece_count={len(outer_pieces)}, pieces={outer_pieces}",
        },
        {
            "check_name": "wrapped_inner_pieces_inside_cell",
            "status": _pieces_inside_status(inner_pieces, half),
            "details": f"inner_piece_count={len(inner_pieces)}, pieces={inner_pieces}",
        },
    ]
    return pd.DataFrame(rows)


def _pieces_inside_status(pieces: list[list[tuple[float, float]]], half: float) -> str:
    if not pieces:
        return "FAIL"
    for piece in pieces:
        for x, _ in piece:
            if x < -half - 1e-8 or x > half + 1e-8:
                return "FAIL"
    return "PASS"


def _plot_geometry(paths: dict[str, Path], args) -> None:
    outer = _slanted_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=0.0,
    )
    inner = _slanted_vertices(
        top_width_um=args.top_width_um - 2.0 * args.film_thickness_um,
        bottom_width_um=args.bottom_width_um - 2.0 * args.film_thickness_um,
        depth_um=args.depth_um - args.film_thickness_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=0.0,
    )
    outer_pieces = _wrapped_polygons(outer, args.period_um)
    inner_pieces = _wrapped_polygons(inner, args.period_um)
    half = args.period_um / 2.0
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.add_patch(plt.Rectangle((-half, -args.substrate_thickness_um), args.period_um, args.substrate_thickness_um, color="#B8B8B8"))
    for piece in outer_pieces:
        ax.add_patch(MplPolygon(piece, closed=True, facecolor="#F28E2B", edgecolor="black", alpha=0.85))
    for piece in inner_pieces:
        ax.add_patch(MplPolygon(piece, closed=True, facecolor="white", edgecolor="#1F77B4", alpha=1.0))
    ax.axvline(-half, color="0.25", lw=1)
    ax.axvline(half, color="0.25", lw=1)
    ax.axhline(0, color="0.25", lw=1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-half - 3, half + 3)
    ax.set_ylim(-args.depth_um - 5, 5)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("D13 wrapped 20 deg gap10 groove, 2 um inner-wall film")
    save_figure(fig, paths["geometry"])
    plt.close(fig)


def _plot_normal(normal: pd.DataFrame, paths: dict[str, Path]) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for pol, sub in normal.groupby("polarization", sort=False):
        ax.plot(sub["wavelength_um"], sub["A"], marker="o", lw=1.6, label=pol)
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Emissivity proxy A")
    ax.set_ylim(0, 1.05)
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend()
    ax.set_title("D13 normal-incidence emissivity")
    save_figure(fig, paths["normal_spectra_fig"])
    plt.close(fig)


def _plot_heatmap(angle: pd.DataFrame, pol: str, path: Path) -> None:
    sub = angle[angle["polarization"] == pol]
    if sub.empty:
        return
    pivot = sub.pivot_table(index="theta_deg", columns="wavelength_um", values="A", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    im = ax.imshow(
        pivot.to_numpy(),
        origin="lower",
        aspect="auto",
        extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()],
        vmin=0,
        vmax=1,
        cmap="viridis",
    )
    fig.colorbar(im, ax=ax, label="A(lambda, theta)")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Incident angle theta (deg)")
    ax.set_title(f"D13 angle-resolved emissivity proxy, {pol}")
    save_figure(fig, path)
    plt.close(fig)


def _plot_angle_proxy(angle_metrics: pd.DataFrame, path: Path) -> None:
    per_angle = angle_metrics[
        (angle_metrics["metric_scope"] == "per_angle")
        & (angle_metrics["polarization"].isin(["Ez", "Hz", PROXY]))
    ]
    if per_angle.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for pol, sub in per_angle.groupby("polarization", sort=False):
        ax.plot(sub["theta_deg"], sub["band_mean_A_strict"], marker="o", lw=1.6, label=pol)
    ax.set_xlabel("Incident angle theta (deg)")
    ax.set_ylabel("Band mean A")
    ax.set_ylim(0, 1.05)
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend()
    ax.set_title("D13 angle-resolved band mean")
    save_figure(fig, path)
    plt.close(fig)


def _md_table(df: pd.DataFrame, max_rows: int = 24) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows)
    lines = [
        "| " + " | ".join(view.columns) + " |",
        "| " + " | ".join(["---"] * len(view.columns)) + " |",
    ]
    for _, row in view.iterrows():
        vals = []
        for col in view.columns:
            value = row[col]
            if isinstance(value, float):
                vals.append(f"{value:.6g}")
            else:
                vals.append(str(value).replace("|", "/"))
        lines.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        lines.append("| " + " | ".join(["..."] + [f"{len(df) - max_rows} more rows"] + [""] * (len(view.columns) - 2)) + " |")
    return "\n".join(lines)


def _write_report(paths: dict[str, Path], checks: pd.DataFrame, normal_metrics: pd.DataFrame, angle_metrics: pd.DataFrame, args) -> None:
    ensure_dir(paths["report"].parent)
    normal_view = normal_metrics[
        [
            "polarization",
            "mean_A_8p1014_12p398_strict",
            "peak_A_strict",
            "peak_wavelength_strict_um",
            "valid_point_count_strict",
            "numerical_status",
        ]
    ] if not normal_metrics.empty else pd.DataFrame()
    if "metric_scope" in angle_metrics.columns:
        angle_view = angle_metrics[angle_metrics["metric_scope"] == "band_average_by_angle"]
    else:
        angle_view = pd.DataFrame()
    if "polarization" in angle_view.columns:
        best_proxy = angle_view[angle_view["polarization"] == PROXY]
    else:
        best_proxy = pd.DataFrame()
    if not best_proxy.empty:
        best_line = (
            f"Best proxy angle: theta={float(best_proxy.iloc[0]['best_theta_deg']):.6g} deg, "
            f"band mean A={float(best_proxy.iloc[0]['best_band_mean_A_strict']):.6g}."
        )
    else:
        best_line = "Best proxy angle: not available."
    report = f"""# D13 Wrapped Tilt20 Gap10 Film2um Angle Screen

Overall result level: **NUMERICAL_SCREENING**

## Geometry

- Period: {args.period_um:g} um.
- Top/bottom groove width: {args.top_width_um:g}/{args.bottom_width_um:g} um.
- Top Ti gap: {args.period_um - args.top_width_um:g} um.
- Depth: {args.depth_um:g} um.
- Tilt: {args.tilt_angle_deg:g} deg.
- Bottom shift: {args.depth_um * math.tan(math.radians(args.tilt_angle_deg)):.6g} um.
- Film: {args.film_thickness_um:g} um on sidewalls and bottom.
- Representation: clipped/wrapped periodic polygon pieces, not enlarged-period geometry.

## Geometry Checks

{_md_table(checks, max_rows=20)}

## Normal-Incidence Emissivity

{_md_table(normal_view, max_rows=10)}

## Angle-Resolved Summary

{_md_table(angle_view, max_rows=10)}

{best_line}

## Numerical Scope

- A(lambda, theta) is used as the directional emissivity proxy under reciprocity and local thermal equilibrium.
- This is a 2D wrapped periodic slanted-groove model, not a true 3D hole-array model.
- Angle grid: {min(args.angles_deg):g} to {max(args.angles_deg):g} deg, step approximately {args.angle_step_deg:g} deg.
- Resolution={args.resolution}, Courant={args.courant}, retry Courant={args.retry_courant}, decay_db={args.decay_db}, fwidth_fraction={args.fwidth_fraction}.

## Outputs

- Normal spectra: `{paths['normal_spectra']}`
- Normal metrics: `{paths['normal_metrics']}`
- Angle spectra: `{paths['angle_spectra']}`
- Angle metrics: `{paths['angle_metrics']}`
- Geometry figure: `{paths['geometry']}`
- Log: `{paths['log']}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def _prepare_args(args: argparse.Namespace) -> argparse.Namespace:
    args.wavelengths_um = _parse_float_list(args.wavelengths_um, DEFAULT_WAVELENGTHS_UM)
    if args.angles_deg:
        args.angles_deg = _parse_float_list(args.angles_deg, DEFAULT_ANGLES_DEG)
    else:
        args.angles_deg = list(np.arange(args.angle_min_deg, args.angle_max_deg + 1e-9, args.angle_step_deg))
        if not any(abs(x) < 1e-9 for x in args.angles_deg):
            args.angles_deg.append(0.0)
        args.angles_deg = sorted(set(round(float(x), 10) for x in args.angles_deg))
    args.polarizations = _parse_str_list(args.polarizations, DEFAULT_POLARIZATIONS)
    return args


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nk-csv", type=Path, default=project_path("data", "processed", "measured_lossy_wall_film_nk_sheet2.csv"))
    p.add_argument("--period-um", type=float, default=DEFAULT_PERIOD_UM)
    p.add_argument("--top-width-um", type=float, default=DEFAULT_TOP_WIDTH_UM)
    p.add_argument("--bottom-width-um", type=float, default=DEFAULT_BOTTOM_WIDTH_UM)
    p.add_argument("--depth-um", type=float, default=DEFAULT_DEPTH_UM)
    p.add_argument("--tilt-angle-deg", type=float, default=DEFAULT_TILT_DEG)
    p.add_argument("--film-thickness-um", type=float, default=DEFAULT_FILM_UM)
    p.add_argument("--substrate-thickness-um", type=float, default=DEFAULT_SUBSTRATE_UM)
    p.add_argument("--air-buffer-um", type=float, default=DEFAULT_AIR_UM)
    p.add_argument("--pml-thickness-um", type=float, default=DEFAULT_PML_UM)
    p.add_argument("--resolution", type=int, default=DEFAULT_RESOLUTION)
    p.add_argument("--courant", type=float, default=DEFAULT_COURANT)
    p.add_argument("--retry-courant", type=float, default=DEFAULT_RETRY_COURANT)
    p.add_argument("--decay-db", type=float, default=DEFAULT_DECAY_DB)
    p.add_argument("--fwidth-fraction", type=float, default=DEFAULT_FWIDTH_FRACTION)
    p.add_argument("--wavelengths-um", nargs="*")
    p.add_argument("--polarizations", nargs="*")
    p.add_argument("--angle-min-deg", type=float, default=-70.0)
    p.add_argument("--angle-max-deg", type=float, default=70.0)
    p.add_argument("--angle-step-deg", type=float, default=10.0)
    p.add_argument("--angles-deg", nargs="*")
    p.add_argument("--theta0-deg", type=float, default=30.0)
    p.add_argument("--workers", default="auto")
    p.add_argument("--output-tag", default=TAG)
    p.add_argument("--geometry-only", action="store_true")
    p.add_argument("--skip-normal", action="store_true")
    p.add_argument("--skip-angle", action="store_true")
    p.add_argument(
        "--resume",
        action="store_true",
        help="Reuse valid rows from checkpoint CSVs and run only missing tasks.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _prepare_args(build_arg_parser().parse_args(argv))
    paths = _paths(args.output_tag)
    for path in paths.values():
        ensure_dir(path.parent)
    logger = _setup_logger(paths["log"])
    logger.info("Starting D13 wrapped tilt20 diagnostics with args=%s", vars(args))

    checks = _geometry_checks(args)
    _plot_geometry(paths, args)
    if (checks["status"] != "PASS").any():
        _write_report(paths, checks, pd.DataFrame(), pd.DataFrame(), args)
        failed = checks[checks["status"] != "PASS"]
        raise ValueError("D13 geometry checks failed:\n" + failed.to_string(index=False))

    if args.geometry_only:
        _write_report(paths, checks, pd.DataFrame(), pd.DataFrame(), args)
        logger.info("D13 geometry-only run passed.")
        return 0

    normal = pd.DataFrame()
    normal_metrics = pd.DataFrame()
    if not args.skip_normal:
        normal_tasks = [(pol, wl, None) for pol in args.polarizations for wl in args.wavelengths_um]
        resume_normal, completed = _load_resume_rows("normal", paths["normal_checkpoint"], logger) if args.resume else (pd.DataFrame(), set())
        normal_tasks = _filter_resume_tasks("normal", normal_tasks, completed, logger)
        normal_rows = _run_tasks(
            "normal",
            normal_tasks,
            args,
            paths["log"],
            paths["normal_checkpoint"],
            logger,
            checkpoint_seed=resume_normal if args.resume else None,
        )
        normal_raw = pd.concat([resume_normal, pd.DataFrame(normal_rows)], ignore_index=True)
        if not normal_raw.empty:
            normal_raw.to_csv(paths["normal_checkpoint"], index=False)
        normal = _add_proxy_rows(normal_raw)
        normal.to_csv(paths["normal_spectra"], index=False)
        normal_metrics = _normal_metrics(normal, args)
        normal_metrics.to_csv(paths["normal_metrics"], index=False)
        _plot_normal(normal, paths)
    elif paths["normal_metrics"].exists():
        normal_metrics = pd.read_csv(paths["normal_metrics"])

    angle = pd.DataFrame()
    angle_metrics = pd.DataFrame()
    if not args.skip_angle:
        angle_tasks = [
            (pol, wl, theta)
            for pol in args.polarizations
            for wl in args.wavelengths_um
            for theta in args.angles_deg
        ]
        resume_angle, completed = _load_resume_rows("angle", paths["angle_checkpoint"], logger) if args.resume else (pd.DataFrame(), set())
        angle_tasks = _filter_resume_tasks("angle", angle_tasks, completed, logger)
        angle_rows = _run_tasks(
            "angle",
            angle_tasks,
            args,
            paths["log"],
            paths["angle_checkpoint"],
            logger,
            checkpoint_seed=resume_angle if args.resume else None,
        )
        angle_raw = pd.concat([resume_angle, pd.DataFrame(angle_rows)], ignore_index=True)
        if not angle_raw.empty:
            angle_raw.to_csv(paths["angle_checkpoint"], index=False)
        angle = _add_proxy_rows(angle_raw)
        angle.to_csv(paths["angle_spectra"], index=False)
        angle_metrics = _angle_metrics(angle, args)
        angle_metrics.to_csv(paths["angle_metrics"], index=False)
        _plot_heatmap(angle, "Ez", paths["angle_heatmap_ez"])
        _plot_heatmap(angle, "Hz", paths["angle_heatmap_hz"])
        _plot_angle_proxy(angle_metrics, paths["angle_proxy"])

    _write_report(paths, checks, normal_metrics, angle_metrics, args)
    logger.info("Saved D13 report: %s", paths["report"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
