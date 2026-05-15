"""D04_field_and_absorbed_power_diagnostic.py — field/hotspot diagnostic.

研究目的
--------
在代表性波长处输出 2D Ti 表面结构的 epsilon、|E|^2、|H|^2 与局域耗散功率密度，
判断当前斜槽是否形成了把能量局域到有损 Ti 区域的吸收模式。

物理假设
--------
1. 所有长度单位为 um，Meep 频率 f = 1 / wavelength_um；
2. 正入射、2D 周期模型；Ez 与 Hz 是二维模型中的两种独立偏振；
3. Ti 使用 Meep 内置 Rakić Drude-Lorentz 模型，超过 12.398 um 的点属于外推；
4. A_flux 使用 D00/D01 验证后的通量符号：R=refl/|input|, T=-trans/|input|,
   A=1-R-T；
5. absorbed_power_density 参考 Meep absorbed_power_density 示例，对 DFT 场使用
   yee_grid=True，并用 2*pi*f*Im(conj(E).D) 计算色散介质频域耗散。

通过/失败判据
-------------
1. 每个 case/polarization/wavelength 的 flux R/T/A 必须有限；
2. A_flux 与 A_volume 的差异小于 volume_flux_tol 时视为体积分校核通过；
3. 输出必须保留 raw flux、符号约定、归一化说明；
4. hotspot 比例只作为机制诊断，不作为求解器失败条件。
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
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
from matplotlib.path import Path as MplPath

from src.geometry import build_slanted_groove_geometry, slanted_groove_vertices
from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium


DEFAULTS = dict(
    wavelengths_um=[8.0, 10.0, 12.0],
    use_d02_peaks=True,
    resolution=32,
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    field_air_above_um=2.5,
    fwidth_fraction=0.08,
    decay_db=40.0,
    wall_band_um=0.35,
    bottom_band_um=0.35,
    volume_flux_tol=0.12,
    hotspot_ratio_threshold=0.30,
)


@dataclass(frozen=True)
class CaseSpec:
    case_name: str
    label: str
    kind: str
    tilt_deg: float | None


def _paths() -> dict[str, Path]:
    return {
        "integrals_csv": project_path(
            "results", "diagnostics", "tables",
            "D04_absorbed_power_integrals.csv",
        ),
        "hotspot_csv": project_path(
            "results", "diagnostics", "tables", "D04_hotspot_metrics.csv",
        ),
        "fig_dir": project_path("results", "diagnostics", "figures"),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D04_field_absorption_mechanism_report.md",
        ),
    }


def _case_specs() -> dict[str, CaseSpec]:
    return {
        "flat_Ti": CaseSpec("flat_Ti", "Flat Ti", "flat", None),
        "symmetric_groove": CaseSpec(
            "symmetric_groove", "Symmetric groove, tilt=0", "groove", 0.0,
        ),
        "slanted_groove": CaseSpec(
            "slanted_groove", "Slanted groove, tilt=20", "groove", 20.0,
        ),
    }


def _safe_wavelength_label(wavelength_um: float) -> str:
    return f"{wavelength_um:.3f}um".replace(".", "p").replace("-", "m")


def _parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return list(default)
    out = []
    for item in values:
        for token in item.replace(",", " ").split():
            out.append(float(token))
    return out


def _unique_sorted(values: list[float], tol: float = 1e-6) -> list[float]:
    out: list[float] = []
    for val in sorted(float(v) for v in values):
        if not out or abs(val - out[-1]) > tol:
            out.append(val)
    return out


def _d02_peak_wavelengths(enabled: bool, logger) -> list[float]:
    if not enabled:
        return []
    path = project_path(
        "results", "diagnostics", "tables", "D02_polarization_metrics.csv",
    )
    if not path.is_file():
        logger.warning("D02 metrics not found: %s; skip automatic peak wavelengths", path)
        return []
    df = pd.read_csv(path)
    peaks = []
    for pol in ("Ez", "Hz"):
        sub = df[df["polarization"] == pol]
        if sub.empty or "peak_A_8_13um" not in sub:
            continue
        idx = sub["peak_A_8_13um"].astype(float).idxmax()
        peaks.append(float(sub.loc[idx, "peak_wavelength_um"]))
    logger.info("D02 peak wavelengths added: %s", peaks)
    return peaks


def _cell_layout(args: argparse.Namespace) -> dict[str, float]:
    bottom_buffer_um = args.pml_thickness_um
    cell_y = (
        2.0 * args.pml_thickness_um
        + args.air_buffer_um
        + args.substrate_thickness_um
        + bottom_buffer_um
    )
    y_top_edge = +0.5 * cell_y
    y_bottom_edge = -0.5 * cell_y
    y_top_pml_inner = y_top_edge - args.pml_thickness_um
    y_bottom_pml_inner = y_bottom_edge + args.pml_thickness_um
    y_surface = y_top_pml_inner - args.air_buffer_um
    y_substrate_bottom = y_surface - args.substrate_thickness_um
    y_src = y_top_pml_inner - 0.25 * args.air_buffer_um
    y_refl = y_surface + 0.5 * args.air_buffer_um
    y_trans = 0.5 * (y_substrate_bottom + y_bottom_pml_inner)
    y_dft_max = min(y_surface + args.field_air_above_um, y_refl)
    y_dft_min = y_substrate_bottom
    return dict(
        cell_y=cell_y,
        y_surface=y_surface,
        y_substrate_bottom=y_substrate_bottom,
        y_src=y_src,
        y_refl=y_refl,
        y_trans=y_trans,
        y_dft_min=y_dft_min,
        y_dft_max=y_dft_max,
        y_dft_center=0.5 * (y_dft_min + y_dft_max),
        y_dft_size=y_dft_max - y_dft_min,
    )


def _build_geometry(case: CaseSpec, args: argparse.Namespace, ti, y_surface: float) -> list:
    if case.kind == "flat":
        return [
            mp.Block(
                material=ti,
                center=mp.Vector3(0, y_surface - 0.5 * args.substrate_thickness_um, 0),
                size=mp.Vector3(args.period_um, args.substrate_thickness_um, mp.inf),
            )
        ]
    return build_slanted_groove_geometry(
        period_x_um=args.period_um,
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=float(case.tilt_deg),
        substrate_thickness_um=args.substrate_thickness_um,
        y_surface=y_surface,
        medium_substrate=ti,
    )


def _components_for_polarization(polarization: str) -> tuple[list[int], list[int], list[int], int]:
    if polarization == "Ez":
        e_components = [mp.Ez]
        d_components = [mp.Dz]
        h_components = [mp.Hx, mp.Hy]
        source_component = mp.Ez
    elif polarization == "Hz":
        e_components = [mp.Ex, mp.Ey]
        d_components = [mp.Dx, mp.Dy]
        h_components = [mp.Hz]
        source_component = mp.Hz
    else:
        raise ValueError(f"polarization must be Ez or Hz, got {polarization}")
    return e_components, d_components, h_components, source_component


def _add_dft(sim, components: list[int], fcen: float, center, size):
    # Meep's absorbed-power-density example uses yee_grid=True for DFT fields in
    # dispersive media; keep that explicit here for provenance.
    return sim.add_dft_fields(
        components,
        fcen,
        0,
        1,
        center=center,
        size=size,
        yee_grid=True,
    )


def _crop_to_common(arrays: list[np.ndarray]) -> list[np.ndarray]:
    min_shape = tuple(min(a.shape[i] for a in arrays) for i in range(arrays[0].ndim))
    slices = tuple(slice(0, n) for n in min_shape)
    return [np.asarray(a)[slices] for a in arrays]


def _sum_abs2(arrays: list[np.ndarray]) -> np.ndarray:
    cropped = _crop_to_common(arrays)
    out = np.zeros(cropped[0].shape, dtype=float)
    for arr in cropped:
        out += np.abs(arr) ** 2
    return out


def _absorbed_power_density(
    e_arrays: list[np.ndarray],
    d_arrays: list[np.ndarray],
    frequency: float,
) -> np.ndarray:
    e_crop, d_crop = _crop_to_common(e_arrays), _crop_to_common(d_arrays)
    min_shape = tuple(min(a.shape[i] for a in e_crop + d_crop) for i in range(e_crop[0].ndim))
    slices = tuple(slice(0, n) for n in min_shape)
    density = np.zeros(min_shape, dtype=float)
    for e_arr, d_arr in zip(e_crop, d_crop):
        density += np.imag(np.conj(e_arr[slices]) * d_arr[slices])
    return 2.0 * math.pi * frequency * density


def _grid_for_array(arr: np.ndarray, args: argparse.Namespace, layout: dict[str, float]):
    nx, ny = arr.shape[:2]
    x = np.linspace(-0.5 * args.period_um, 0.5 * args.period_um, nx)
    y = np.linspace(layout["y_dft_min"], layout["y_dft_max"], ny)
    return x, y


def _line_distance(px, py, x1, y1, x2, y2):
    vx = x2 - x1
    vy = y2 - y1
    denom = vx * vx + vy * vy
    if denom == 0:
        return np.hypot(px - x1, py - y1)
    t = ((px - x1) * vx + (py - y1) * vy) / denom
    t = np.clip(t, 0.0, 1.0)
    proj_x = x1 + t * vx
    proj_y = y1 + t * vy
    return np.hypot(px - proj_x, py - proj_y)


def _masks(
    case: CaseSpec,
    arr: np.ndarray,
    args: argparse.Namespace,
    layout: dict[str, float],
) -> dict[str, np.ndarray]:
    x, y = _grid_for_array(arr, args, layout)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    substrate_mask = (
        (yy <= layout["y_surface"])
        & (yy >= layout["y_surface"] - args.substrate_thickness_um)
    )
    if case.kind == "flat":
        groove_mask = np.zeros_like(substrate_mask, dtype=bool)
        wall_mask = np.zeros_like(substrate_mask, dtype=bool)
        bottom_mask = np.zeros_like(substrate_mask, dtype=bool)
    else:
        verts = slanted_groove_vertices(
            top_width_um=args.top_width_um,
            bottom_width_um=args.bottom_width_um,
            depth_um=args.depth_um,
            tilt_angle_deg=float(case.tilt_deg),
            y_surface=layout["y_surface"],
        )
        points = np.column_stack([xx.ravel(), yy.ravel()])
        groove_mask = MplPath(verts).contains_points(points).reshape(xx.shape)
        metal_mask = substrate_mask & ~groove_mask

        tl, tr, br, bl = verts
        left_dist = _line_distance(xx, yy, tl[0], tl[1], bl[0], bl[1])
        right_dist = _line_distance(xx, yy, tr[0], tr[1], br[0], br[1])
        bottom_dist = _line_distance(xx, yy, bl[0], bl[1], br[0], br[1])
        wall_mask = metal_mask & (
            (left_dist <= args.wall_band_um) | (right_dist <= args.wall_band_um)
        )
        bottom_mask = metal_mask & (bottom_dist <= args.bottom_band_um)
        return dict(
            substrate=substrate_mask,
            groove=groove_mask,
            metal=metal_mask,
            wall=wall_mask,
            bottom=bottom_mask,
        )

    metal_mask = substrate_mask & ~groove_mask
    return dict(
        substrate=substrate_mask,
        groove=groove_mask,
        metal=metal_mask,
        wall=wall_mask,
        bottom=bottom_mask,
    )


def _integrate_density(density: np.ndarray, args: argparse.Namespace, layout: dict[str, float]) -> float:
    dx = args.period_um / max(density.shape[0] - 1, 1)
    dy = layout["y_dft_size"] / max(density.shape[1] - 1, 1)
    return float(np.sum(density) * dx * dy)


def _plot_scalar(
    data: np.ndarray,
    args: argparse.Namespace,
    layout: dict[str, float],
    title: str,
    cbar_label: str,
    out_path: Path,
    *,
    log_scale: bool = False,
    cmap: str = "magma",
) -> Path:
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    x, y = _grid_for_array(data, args, layout)
    plot_data = np.asarray(data, dtype=float)
    if log_scale:
        positive = plot_data[np.isfinite(plot_data) & (plot_data > 0)]
        floor = np.percentile(positive, 1) if positive.size else 1e-30
        plot_data = np.log10(np.maximum(plot_data, max(floor, 1e-30)))
        cbar_label = f"log10({cbar_label})"
    im = ax.imshow(
        plot_data.T,
        origin="lower",
        extent=[x.min(), x.max(), y.min(), y.max()],
        aspect="auto",
        cmap=cmap,
    )
    ax.axhline(layout["y_surface"], color="white", lw=0.8, ls=":", alpha=0.8)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    save_figure(fig, out_path)
    plt.close(fig)
    return out_path


def _run_one(case: CaseSpec, polarization: str, wavelength_um: float,
             args: argparse.Namespace, logger) -> tuple[dict, dict, list[Path]]:
    ti = get_ti_medium(lambda_min_um=wavelength_um, lambda_max_um=wavelength_um)
    layout = _cell_layout(args)
    fcen = 1.0 / wavelength_um
    fwidth = args.fwidth_fraction * fcen
    cell = mp.Vector3(args.period_um, layout["cell_y"], 0)
    pml_layers = [mp.PML(thickness=args.pml_thickness_um, direction=mp.Y)]
    e_components, d_components, h_components, source_component = _components_for_polarization(
        polarization
    )
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=fwidth, is_integrated=True),
            component=source_component,
            center=mp.Vector3(0, layout["y_src"], 0),
            size=mp.Vector3(args.period_um, 0, 0),
        )
    ]

    logger.info(">>> D04 ref case=%s pol=%s wl=%.4g", case.case_name, polarization, wavelength_um)
    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        k_point=mp.Vector3(),
        geometry=[],
        dimensions=2,
    )
    refl_ref = sim_ref.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(
            center=mp.Vector3(0, layout["y_refl"], 0),
            size=mp.Vector3(args.period_um, 0, 0),
        ),
    )
    t0 = time.time()
    sim_ref.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=source_component,
            pt=mp.Vector3(0, layout["y_refl"], 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    input_flux_raw = float(np.array(mp.get_fluxes(refl_ref))[0])
    ref_data = sim_ref.get_flux_data(refl_ref)
    ref_time = time.time() - t0
    if input_flux_raw == 0:
        raise RuntimeError("reference input flux is zero")

    geometry = _build_geometry(case, args, ti, layout["y_surface"])
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        k_point=mp.Vector3(),
        geometry=geometry,
        dimensions=2,
    )
    refl = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(
            center=mp.Vector3(0, layout["y_refl"], 0),
            size=mp.Vector3(args.period_um, 0, 0),
        ),
    )
    trans = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(
            center=mp.Vector3(0, layout["y_trans"], 0),
            size=mp.Vector3(args.period_um, 0, 0),
        ),
    )
    sim.load_minus_flux_data(refl, ref_data)

    dft_center = mp.Vector3(0, layout["y_dft_center"], 0)
    dft_size = mp.Vector3(args.period_um, layout["y_dft_size"], 0)
    all_components = e_components + d_components + h_components
    dft = _add_dft(sim, all_components, fcen, dft_center, dft_size)

    logger.info(">>> D04 struct case=%s pol=%s wl=%.4g", case.case_name, polarization, wavelength_um)
    t0 = time.time()
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=source_component,
            pt=mp.Vector3(0, layout["y_refl"], 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    struct_time = time.time() - t0

    refl_flux_raw = float(np.array(mp.get_fluxes(refl))[0])
    trans_flux_raw = float(np.array(mp.get_fluxes(trans))[0])
    R = refl_flux_raw / abs(input_flux_raw)
    T_raw = trans_flux_raw / abs(input_flux_raw)
    T = -trans_flux_raw / abs(input_flux_raw)
    A_flux = 1.0 - R - T

    e_arrays = [np.asarray(sim.get_dft_array(dft, comp, 0)) for comp in e_components]
    d_arrays = [np.asarray(sim.get_dft_array(dft, comp, 0)) for comp in d_components]
    h_arrays = [np.asarray(sim.get_dft_array(dft, comp, 0)) for comp in h_components]
    e2 = _sum_abs2(e_arrays)
    h2 = _sum_abs2(h_arrays)
    absorbed_density = _absorbed_power_density(e_arrays, d_arrays, fcen)
    e2, h2, absorbed_density = _crop_to_common([e2, h2, absorbed_density])

    eps = np.asarray(
        sim.get_array(
            component=mp.Dielectric,
            center=dft_center,
            size=dft_size,
        )
    )
    eps = _crop_to_common([eps, absorbed_density])[0]
    absorbed_density = _crop_to_common([eps, absorbed_density])[1]
    e2 = _crop_to_common([eps, e2])[1]
    h2 = _crop_to_common([eps, h2])[1]

    masks = _masks(case, absorbed_density, args, layout)
    positive_density = np.maximum(absorbed_density, 0.0)
    total_absorbed_power = _integrate_density(absorbed_density, args, layout)
    positive_absorbed_power = _integrate_density(positive_density, args, layout)
    A_volume = total_absorbed_power / abs(input_flux_raw)
    A_volume_positive = positive_absorbed_power / abs(input_flux_raw)
    wall_power = _integrate_density(positive_density * masks["wall"], args, layout)
    bottom_power = _integrate_density(positive_density * masks["bottom"], args, layout)
    metal_power = _integrate_density(positive_density * masks["metal"], args, layout)
    denom = positive_absorbed_power if positive_absorbed_power > 0 else np.nan

    metal_e2 = e2[masks["metal"]]
    air_e2 = e2[~masks["substrate"]]
    groove_e2 = e2[masks["groove"]]
    wall_e2 = e2[masks["wall"]]
    bottom_e2 = e2[masks["bottom"]]
    air_reference = float(np.nanmedian(air_e2)) if air_e2.size else np.nan
    if not np.isfinite(air_reference) or air_reference <= 0:
        air_reference = np.nan

    case_label = f"{case.case_name}_{polarization}_{_safe_wavelength_label(wavelength_um)}"
    fig_dir = _paths()["fig_dir"]
    figures = [
        _plot_scalar(
            eps.real,
            args,
            layout,
            f"D04 {case.case_name} {polarization} {wavelength_um:g} um epsilon",
            "epsilon",
            fig_dir / f"D04_{case_label}_epsilon.png",
            cmap="viridis",
        ),
        _plot_scalar(
            e2,
            args,
            layout,
            f"D04 {case.case_name} {polarization} {wavelength_um:g} um |E|^2",
            "|E|^2",
            fig_dir / f"D04_{case_label}_E2.png",
            log_scale=True,
            cmap="magma",
        ),
        _plot_scalar(
            h2,
            args,
            layout,
            f"D04 {case.case_name} {polarization} {wavelength_um:g} um |H|^2",
            "|H|^2",
            fig_dir / f"D04_{case_label}_H2.png",
            log_scale=True,
            cmap="plasma",
        ),
        _plot_scalar(
            positive_density,
            args,
            layout,
            f"D04 {case.case_name} {polarization} {wavelength_um:g} um absorbed power",
            "positive absorbed power density",
            fig_dir / f"D04_{case_label}_absorbed_power.png",
            log_scale=True,
            cmap="inferno",
        ),
    ]

    finite_flux = bool(np.all(np.isfinite([R, T, A_flux, input_flux_raw, refl_flux_raw, trans_flux_raw])))
    volume_flux_abs_diff = abs(A_volume - A_flux)
    integrals = dict(
        case_name=case.case_name,
        case_label=case.label,
        polarization=polarization,
        wavelength_um=wavelength_um,
        frequency_meep=fcen,
        input_flux_raw=input_flux_raw,
        reflection_flux_raw=refl_flux_raw,
        transmission_flux_raw=trans_flux_raw,
        reflectance=R,
        raw_transmittance=T_raw,
        transmittance=T,
        absorptance_flux=A_flux,
        absorbed_power_volume_raw=total_absorbed_power,
        absorbed_power_volume_positive=positive_absorbed_power,
        absorptance_volume=A_volume,
        absorptance_volume_positive=A_volume_positive,
        volume_flux_abs_difference=volume_flux_abs_diff,
        volume_flux_rel_difference=(
            volume_flux_abs_diff / abs(A_flux) if A_flux != 0 else np.nan
        ),
        pass_or_fail="PASS" if finite_flux and volume_flux_abs_diff <= args.volume_flux_tol else "FAIL",
        normalization_note=(
            "R = reflection_flux_raw / abs(input_flux_raw); "
            "T = -transmission_flux_raw / abs(input_flux_raw); "
            "A_flux = 1 - R - T; A_volume = integral(2*pi*f*Im(conj(E).D)) / abs(input_flux_raw)"
        ),
        absorbed_power_note=(
            "DFT fields requested with yee_grid=True; positive integrals clip small negative numerical artifacts "
            "for hotspot ratios only"
        ),
        walltime_s=ref_time + struct_time,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        air_buffer_um=args.air_buffer_um,
        substrate_thickness_um=args.substrate_thickness_um,
        decay_db=args.decay_db,
    )

    hotspots = dict(
        case_name=case.case_name,
        case_label=case.label,
        polarization=polarization,
        wavelength_um=wavelength_um,
        total_positive_absorbed_power=positive_absorbed_power,
        metal_positive_absorbed_power=metal_power,
        power_near_groove_wall=wall_power,
        power_near_groove_bottom=bottom_power,
        wall_power_fraction=wall_power / denom if np.isfinite(denom) else np.nan,
        bottom_power_fraction=bottom_power / denom if np.isfinite(denom) else np.nan,
        metal_power_fraction=metal_power / denom if np.isfinite(denom) else np.nan,
        max_E2_total=float(np.nanmax(e2)),
        max_E2_metal=float(np.nanmax(metal_e2)) if metal_e2.size else np.nan,
        max_E2_groove_air=float(np.nanmax(groove_e2)) if groove_e2.size else np.nan,
        max_E2_wall_band=float(np.nanmax(wall_e2)) if wall_e2.size else np.nan,
        max_E2_bottom_band=float(np.nanmax(bottom_e2)) if bottom_e2.size else np.nan,
        median_E2_air_above_surface=air_reference,
        max_E2_metal_over_air_median=(
            float(np.nanmax(metal_e2)) / air_reference
            if metal_e2.size and np.isfinite(air_reference) else np.nan
        ),
        max_E2_wall_over_air_median=(
            float(np.nanmax(wall_e2)) / air_reference
            if wall_e2.size and np.isfinite(air_reference) else np.nan
        ),
        max_E2_bottom_over_air_median=(
            float(np.nanmax(bottom_e2)) / air_reference
            if bottom_e2.size and np.isfinite(air_reference) else np.nan
        ),
        hotspot_note=(
            "wall/bottom fractions use positive absorbed power and geometry masks around groove boundaries; "
            "flat_Ti has no groove wall/bottom region"
        ),
    )
    return integrals, hotspots, figures


def _write_report(integrals: pd.DataFrame, hotspots: pd.DataFrame,
                  args: argparse.Namespace, figures: list[Path]) -> Path:
    paths = _paths()
    ensure_dir(paths["report"].parent)
    if integrals.empty:
        raise RuntimeError("no D04 results to report")
    all_pass = bool((integrals["pass_or_fail"] == "PASS").all())
    by_case = hotspots.groupby(["case_name", "polarization"], sort=False).agg(
        mean_wall_fraction=("wall_power_fraction", "mean"),
        mean_bottom_fraction=("bottom_power_fraction", "mean"),
        max_metal_e2_enhancement=("max_E2_metal_over_air_median", "max"),
        max_wall_e2_enhancement=("max_E2_wall_over_air_median", "max"),
        max_bottom_e2_enhancement=("max_E2_bottom_over_air_median", "max"),
    ).reset_index()

    slanted = by_case[by_case["case_name"] == "slanted_groove"]
    flat = by_case[by_case["case_name"] == "flat_Ti"]
    hz_row = slanted[slanted["polarization"] == "Hz"]
    ez_row = slanted[slanted["polarization"] == "Ez"]
    if not hz_row.empty and not ez_row.empty:
        hz_more_wall = (
            float(hz_row["mean_wall_fraction"].iloc[0])
            > float(ez_row["mean_wall_fraction"].iloc[0])
        )
        hz_more_bottom = (
            float(hz_row["mean_bottom_fraction"].iloc[0])
            > float(ez_row["mean_bottom_fraction"].iloc[0])
        )
    else:
        hz_more_wall = False
        hz_more_bottom = False

    max_slanted_e2 = (
        float(slanted["max_metal_e2_enhancement"].max()) if not slanted.empty else np.nan
    )
    has_local_field = bool(np.isfinite(max_slanted_e2) and max_slanted_e2 > 5.0)
    max_slanted_wall_frac = (
        float(slanted["mean_wall_fraction"].max()) if not slanted.empty else np.nan
    )
    effective_loss_mode = bool(
        np.isfinite(max_slanted_wall_frac)
        and max_slanted_wall_frac >= args.hotspot_ratio_threshold
    )

    def md_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_No rows._"
        cols = list(df.columns)
        rows = ["| " + " | ".join(cols) + " |",
                "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, row in df.iterrows():
            cells = []
            for col in cols:
                val = row[col]
                if isinstance(val, float):
                    cells.append(f"{val:.6g}")
                else:
                    cells.append(str(val))
            rows.append("| " + " | ".join(cells) + " |")
        return "\n".join(rows)

    lines = [
        "# D04 Field and Absorbed Power Diagnostic",
        "",
        "## Purpose",
        "Check whether representative Ti grooves localize electromagnetic energy into lossy Ti regions.",
        "",
        "## Physical assumptions",
        "- Length unit: um; Meep frequency: f = 1 / wavelength_um.",
        "- 2D periodic, normal incidence, independent Ez/Hz polarizations.",
        "- Ti: Meep built-in Rakić Drude-Lorentz model; wavelengths above "
        f"{TI_RAKIC_VALID_LAMBDA_UM[1]} um are extrapolated.",
        "- Flux convention: R = reflection_flux_raw / abs(input_flux_raw); "
        "T = -transmission_flux_raw / abs(input_flux_raw); A_flux = 1 - R - T.",
        "- Absorbed power density follows Meep's absorbed_power_density example: "
        "DFT fields are requested with yee_grid=True and evaluated as "
        "2*pi*f*Im(conj(E).D). Hotspot fractions use positive-clipped density to avoid "
        "small negative numerical artifacts dominating ratios.",
        "",
        "## Pass/fail criteria",
        f"- A_flux and A_volume absolute difference <= {args.volume_flux_tol:g}.",
        "- R/T/A and raw flux fields finite.",
        f"- Wall/bottom hotspot ratio >= {args.hotspot_ratio_threshold:g} is treated as a "
        "mechanism flag, not a solver failure.",
        "",
        "## Numerical results",
        md_table(integrals[[
            "case_name", "polarization", "wavelength_um", "absorptance_flux",
            "absorptance_volume", "volume_flux_abs_difference", "pass_or_fail",
        ]]),
        "",
        "## Hotspot summary",
        md_table(by_case),
        "",
        "## Required answers",
        f"1. Does the current slanted groove produce clear local field enhancement? "
        f"{'Numerically indicated' if has_local_field else 'Not clearly indicated'} "
        f"(max metal |E|^2 / air-median |E|^2 = {max_slanted_e2:.3g}).",
        f"2. Is the enhanced field located in lossy Ti rather than mainly in air? "
        f"{'Partly yes' if effective_loss_mode else 'Not strongly supported by this run'} "
        f"(largest mean wall absorbed-power fraction = {max_slanted_wall_frac:.3g}).",
        f"3. Does Hz more readily form groove-wall or groove-bottom dissipation than Ez? "
        f"{'Yes, in this numerical run' if (hz_more_wall or hz_more_bottom) else 'Not proven by this run'}.",
        f"4. Is low spectral absorption explained by no effective localized loss mode? "
        f"{'Supported as a numerical interpretation' if not effective_loss_mode else 'Not fully supported; localized loss exists but may be insufficient in total power'}.",
        "",
        "## Verified conclusions",
        f"- Flux/volume consistency overall status: {'PASS' if all_pass else 'FAIL'} for the executed set.",
        "- Generated field maps and integral tables retain raw flux and normalization provenance.",
        "",
        "## Hypotheses",
        "- If wall/bottom absorbed-power fractions remain small in the full run, the simple 2D groove likely fails "
        "to create an efficient localized Ti loss channel in 8-13 um.",
        "- Differences between Ez and Hz remain 2D polarization effects, not direct proof of 3D unpolarized emissivity.",
        "",
        "## Needs higher-fidelity confirmation",
        "- 3D finite grooves, oxide layers, roughness, rounded sidewalls, and measured Ti optical constants.",
        "- More converged D04 runs if A_volume differs strongly from A_flux in any case.",
        "",
        "## Output files",
        f"- `{paths['integrals_csv']}`",
        f"- `{paths['hotspot_csv']}`",
        f"- Figures generated: {len(figures)} PNG files under `{paths['fig_dir']}`",
        "",
        "## Run configuration",
        "```json",
        json.dumps(vars(args), indent=2, ensure_ascii=False),
        "```",
    ]
    paths["report"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths["report"]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="D04 field and absorbed-power-density diagnostic for Ti grooves."
    )
    p.add_argument("--wavelengths_um", nargs="*", default=None,
                   help="Wavelengths in um, e.g. --wavelengths_um 8 10 12")
    p.add_argument("--no_d02_peaks", action="store_true",
                   help="Do not append D02 peak wavelengths.")
    p.add_argument("--cases", nargs="+", default=list(_case_specs().keys()),
                   choices=list(_case_specs().keys()))
    p.add_argument("--polarizations", nargs="+", default=["Ez", "Hz"], choices=["Ez", "Hz"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--period_um", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--top_width_um", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width_um", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--field_air_above_um", type=float,
                   default=DEFAULTS["field_air_above_um"])
    p.add_argument("--fwidth_fraction", type=float, default=DEFAULTS["fwidth_fraction"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--wall_band_um", type=float, default=DEFAULTS["wall_band_um"])
    p.add_argument("--bottom_band_um", type=float, default=DEFAULTS["bottom_band_um"])
    p.add_argument("--volume_flux_tol", type=float, default=DEFAULTS["volume_flux_tol"])
    p.add_argument("--hotspot_ratio_threshold", type=float,
                   default=DEFAULTS["hotspot_ratio_threshold"])
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    logger = setup_logger("D04_field_and_absorbed_power_diagnostic")
    logger.info("=== D04_field_and_absorbed_power_diagnostic ===")
    logger.info("args=%s", vars(args))

    paths = _paths()
    for path in paths.values():
        ensure_dir(path if path.suffix == "" else path.parent)

    wavelengths = _parse_float_list(args.wavelengths_um, DEFAULTS["wavelengths_um"])
    wavelengths += _d02_peak_wavelengths(not args.no_d02_peaks, logger)
    wavelengths = _unique_sorted(wavelengths)
    logger.info("wavelengths_um=%s", wavelengths)

    cases = _case_specs()
    integrals_rows: list[dict] = []
    hotspot_rows: list[dict] = []
    figures: list[Path] = []
    for case_name in args.cases:
        case = cases[case_name]
        for pol in args.polarizations:
            for wl in wavelengths:
                integrals, hotspots, figs = _run_one(case, pol, wl, args, logger)
                integrals_rows.append(integrals)
                hotspot_rows.append(hotspots)
                figures.extend(figs)

    integrals_df = pd.DataFrame(integrals_rows)
    hotspots_df = pd.DataFrame(hotspot_rows)
    ensure_dir(paths["integrals_csv"].parent)
    integrals_df.to_csv(paths["integrals_csv"], index=False)
    hotspots_df.to_csv(paths["hotspot_csv"], index=False)
    report = _write_report(integrals_df, hotspots_df, args, figures)

    logger.info("saved integrals: %s", paths["integrals_csv"])
    logger.info("saved hotspots: %s", paths["hotspot_csv"])
    logger.info("saved report: %s", report)
    logger.info("saved %d figures under %s", len(figures), paths["fig_dir"])
    if not (integrals_df["pass_or_fail"] == "PASS").all():
        logger.warning("D04 finished with FAIL rows; inspect A_flux vs A_volume")
        return 1
    logger.info("D04 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
