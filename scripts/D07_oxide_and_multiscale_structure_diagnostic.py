"""D07 oxide and multiscale structure sensitivity diagnostic.

This script adds a replaceable TiO2 interface and simplified oxidized-geometry
models to test whether missing oxide/multiscale surface effects could explain
the weak emissivity enhancement of the bare Ti slanted groove.

Important limitation
--------------------
The default TiO2 medium is a placeholder/demo lossless dielectric.  Results are
for sensitivity testing only and are not quantitative TiO2 mid-infrared
predictions.
"""

from __future__ import annotations

import argparse
import json
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
from matplotlib.patches import Polygon, Rectangle
import meep as mp
import numpy as np
import pandas as pd

from src.geometry import (
    build_oxidized_slanted_groove_geometry,
    build_slanted_groove_geometry,
    slanted_groove_vertices,
)
from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import (
    TI_RAKIC_VALID_LAMBDA_UM,
    get_ti_medium,
    get_tio2_medium,
)
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
    tilt_deg=20.0,
    oxide_thicknesses_um=[0.0, 0.05, 0.1, 0.2, 0.5, 1.0],
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    decay_db=40.0,
    tio2_demo_index=2.4,
    enhancement_tol=0.03,
    degeneracy_tol=1e-6,
    geometry_check_thickness_um=0.5,
)


@dataclass(frozen=True)
class CaseSpec:
    case_name: str
    case_label: str
    substrate: str
    structure: str
    oxide_mode: str


def _paths() -> dict[str, Path]:
    return {
        "spectra": project_path(
            "results", "diagnostics", "tables", "D07_oxide_thickness_sweep.csv",
        ),
        "metrics": project_path(
            "results", "diagnostics", "tables",
            "D07_oxide_enhancement_metrics.csv",
        ),
        "geometry": project_path(
            "results", "diagnostics", "figures", "D07_oxide_geometry_check.png",
        ),
        "spectra_ez": project_path(
            "results", "diagnostics", "figures",
            "D07_oxide_spectral_comparison_Ez.png",
        ),
        "spectra_hz": project_path(
            "results", "diagnostics", "figures",
            "D07_oxide_spectral_comparison_Hz.png",
        ),
        "mean_a": project_path(
            "results", "diagnostics", "figures",
            "D07_mean_A_vs_oxide_thickness.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D07_oxide_multiscale_diagnostic_report.md",
        ),
        "best_e2": project_path(
            "results", "diagnostics", "figures", "D07_best_case_E2.png",
        ),
        "best_h2": project_path(
            "results", "diagnostics", "figures", "D07_best_case_H2.png",
        ),
        "best_absorbed": project_path(
            "results", "diagnostics", "figures",
            "D07_best_case_absorbed_power.png",
        ),
    }


def _case_specs() -> dict[str, CaseSpec]:
    return {
        "bare_flat_Ti": CaseSpec(
            "bare_flat_Ti", "Bare flat Ti", "Ti", "flat", "none",
        ),
        "oxide_flat_Ti": CaseSpec(
            "oxide_flat_Ti", "Flat Ti + top oxide film", "Ti", "flat",
            "top_film_only",
        ),
        "bare_slanted_groove": CaseSpec(
            "bare_slanted_groove", "Bare slanted groove", "Ti", "slanted",
            "none",
        ),
        "oxide_slanted_groove_top_film_only": CaseSpec(
            "oxide_slanted_groove_top_film_only",
            "Slanted groove + top oxide film", "Ti", "slanted",
            "top_film_only",
        ),
        "oxide_slanted_groove_conformal_approx": CaseSpec(
            "oxide_slanted_groove_conformal_approx",
            "Slanted groove + conformal oxide approx", "Ti", "slanted",
            "conformal_approx",
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


def _band(df: pd.DataFrame, lo: float = 8.0, hi: float = 13.0) -> pd.DataFrame:
    return df[(df["wavelength_um"] >= lo) & (df["wavelength_um"] <= hi)]


def _flat_geometry_factory(spec: CaseSpec, thickness: float, ti, tio2, args):
    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        geom = [
            mp.Block(
                material=ti,
                center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
            )
        ]
        if spec.oxide_mode != "none" and thickness > 0:
            geom.append(
                mp.Block(
                    material=tio2,
                    center=mp.Vector3(0, y_surface_um + 0.5 * thickness, 0),
                    size=mp.Vector3(args.period_um, thickness, mp.inf),
                )
            )
        return geom
    return factory


def _slanted_geometry_factory(spec: CaseSpec, thickness: float, ti, tio2, args):
    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        if spec.oxide_mode == "none":
            return build_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.top_width_um,
                bottom_width_um=args.bottom_width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=args.tilt_deg,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
            )
        return build_oxidized_slanted_groove_geometry(
            period_x_um=args.period_um,
            top_width_um=args.top_width_um,
            bottom_width_um=args.bottom_width_um,
            depth_um=args.depth_um,
            tilt_angle_deg=args.tilt_deg,
            substrate_thickness_um=substrate_thickness_um,
            y_surface=y_surface_um,
            oxide_thickness_um=thickness,
            medium_substrate=ti,
            medium_oxide=tio2,
            oxide_mode=spec.oxide_mode,
        )
    return factory


def _run_one(spec: CaseSpec, pol: str, thickness: float,
             args: argparse.Namespace, logger) -> pd.DataFrame:
    ti = get_ti_medium(args.wavelength_min_um, args.wavelength_max_um)
    tio2 = get_tio2_medium(
        args.wavelength_min_um,
        args.wavelength_max_um,
        demo_index=args.tio2_demo_index,
        allow_placeholder=True,
    )
    if spec.structure == "flat":
        factory = _flat_geometry_factory(spec, thickness, ti, tio2, args)
    else:
        factory = _slanted_geometry_factory(spec, thickness, ti, tio2, args)

    logger.info(
        ">>> case=%s pol=%s oxide=%.4g mode=%s",
        spec.case_name, pol, thickness, spec.oxide_mode,
    )
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
        source_component=pol,
        logger=logger,
    )
    return pd.DataFrame({
        "case_name": spec.case_name,
        "case_label": spec.case_label,
        "polarization": pol,
        "wavelength_um": result["wavelength_um"],
        "reflectance": result["R"],
        "transmittance": result["T"],
        "absorptance": result["A"],
        "emissivity_proxy": result["A"],
        "oxide_thickness_um": thickness,
        "oxide_mode": spec.oxide_mode,
        "tio2_model": "placeholder_demo",
        "tio2_demo_index": args.tio2_demo_index,
        "placeholder_warning": (
            "TiO2 placeholder_demo is lossless/non-dispersive and only for "
            "sensitivity testing; not a quantitative mid-IR TiO2 model."
        ),
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "normalization_note": (
            "R=reflection_flux_raw/abs(input_flux_raw); "
            "T=-transmission_flux_raw/abs(input_flux_raw); A=1-R-T"
        ),
        "period_um": args.period_um,
        "top_width_um": 0.0 if spec.structure == "flat" else args.top_width_um,
        "bottom_width_um": 0.0 if spec.structure == "flat" else args.bottom_width_um,
        "depth_um": 0.0 if spec.structure == "flat" else args.depth_um,
        "tilt_deg": np.nan if spec.structure == "flat" else args.tilt_deg,
        "walltime_s": float(result["walltime_s"]),
    })


def _metrics(spectra: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    base_lookup = {}
    for (case_name, pol, thickness, mode), df in spectra.groupby(
        ["case_name", "polarization", "oxide_thickness_um", "oxide_mode"],
        sort=False,
    ):
        band = _band(df)
        if band.empty:
            raise RuntimeError(f"No 8-13 um samples for {case_name} {pol} t={thickness}")
        peak_idx = band["absorptance"].idxmax()
        row = dict(
            case_name=case_name,
            case_label=str(df["case_label"].iloc[0]),
            polarization=pol,
            oxide_thickness_um=thickness,
            oxide_mode=mode,
            mean_A_8_13um=float(band["absorptance"].mean()),
            peak_A=float(band.loc[peak_idx, "absorptance"]),
            peak_wavelength_um=float(band.loc[peak_idx, "wavelength_um"]),
            min_A=float(df["absorptance"].min()),
            max_A=float(df["absorptance"].max()),
            all_finite=bool(np.all(np.isfinite(df[[
                "reflectance", "transmittance", "absorptance",
                "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw",
            ]].to_numpy()))),
            has_raw_flux=True,
        )
        if case_name == "bare_flat_Ti":
            base_lookup[("flat", pol)] = row["mean_A_8_13um"]
        if case_name == "bare_slanted_groove":
            base_lookup[("slanted", pol)] = row["mean_A_8_13um"]
        rows.append(row)

    metrics = pd.DataFrame(rows)
    enhancements = []
    for _, row in metrics.iterrows():
        key = ("flat", row["polarization"]) if "flat" in row["case_name"] else (
            "slanted", row["polarization"]
        )
        base = base_lookup.get(key, np.nan)
        enhancements.append(row["mean_A_8_13um"] - base)
    metrics["enhancement_over_bare"] = enhancements
    metrics["significant_enhancement"] = metrics["enhancement_over_bare"] >= args.enhancement_tol
    metrics["pass_or_fail"] = np.where(metrics["all_finite"], "PASS", "FAIL")
    return metrics


def _check_zero_thickness_degeneracy(metrics: pd.DataFrame,
                                     args: argparse.Namespace) -> pd.DataFrame:
    metrics = metrics.copy()
    metrics["zero_thickness_degeneracy_delta"] = np.nan
    for pol in metrics["polarization"].unique():
        bare_flat = metrics[
            (metrics["case_name"] == "bare_flat_Ti") & (metrics["polarization"] == pol)
        ]
        bare_slanted = metrics[
            (metrics["case_name"] == "bare_slanted_groove") & (metrics["polarization"] == pol)
        ]
        for idx, row in metrics.iterrows():
            if row["polarization"] != pol or abs(row["oxide_thickness_um"]) > 1e-12:
                continue
            ref = bare_flat if "flat" in row["case_name"] else bare_slanted
            if ref.empty:
                continue
            delta = abs(float(row["mean_A_8_13um"]) - float(ref["mean_A_8_13um"].iloc[0]))
            metrics.loc[idx, "zero_thickness_degeneracy_delta"] = delta
            if delta > args.degeneracy_tol:
                metrics.loc[idx, "pass_or_fail"] = "FAIL"
    return metrics


def _plot_geometry_check(args: argparse.Namespace, out_path: Path) -> None:
    t = args.geometry_check_thickness_um
    y_surface = 0.0
    outer = slanted_groove_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_deg,
        y_surface=y_surface,
    )
    inner = slanted_groove_vertices(
        top_width_um=args.top_width_um - 2 * t,
        bottom_width_um=args.bottom_width_um - 2 * t,
        depth_um=args.depth_um - t,
        tilt_angle_deg=args.tilt_deg,
        y_surface=y_surface,
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    half_p = args.period_um / 2
    ax.add_patch(Rectangle(
        (-half_p, -args.substrate_thickness_um),
        args.period_um,
        args.substrate_thickness_um,
        color="#8C8C8C",
        label="Ti substrate",
    ))
    ax.add_patch(Polygon(outer, closed=True, facecolor="#F2C14E",
                         edgecolor="black", label="TiO2 lining (outer)"))
    ax.add_patch(Polygon(inner, closed=True, facecolor="white",
                         edgecolor="#1f77b4", label="Air core"))
    left_w = max(0.0, -args.top_width_um / 2 + half_p)
    right_w = left_w
    if left_w > 0:
        ax.add_patch(Rectangle(
            (-half_p, 0), left_w, t, color="#F2C14E", alpha=0.9,
            label="Top TiO2 film",
        ))
        ax.add_patch(Rectangle((args.top_width_um / 2, 0), right_w, t,
                               color="#F2C14E", alpha=0.9))
    ax.set_xlim(-half_p, half_p)
    ax.set_ylim(-args.substrate_thickness_um, t + 0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(
        "D07 conformal_approx geometry check "
        "(placeholder TiO2, not quantitative)"
    )
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.2)
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_spectra(spectra: pd.DataFrame, pol: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    sub = spectra[spectra["polarization"] == pol]
    if sub.empty:
        ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
        ax.set_axis_off()
    else:
        for (case, thickness, mode), df in sub.groupby(
            ["case_name", "oxide_thickness_um", "oxide_mode"], sort=False,
        ):
            if thickness not in sorted(sub["oxide_thickness_um"].unique())[::max(1, len(sub["oxide_thickness_um"].unique()) // 4)]:
                if thickness != sub["oxide_thickness_um"].max():
                    continue
            label = f"{case}, t={thickness:g}, {mode}"
            ax.plot(df["wavelength_um"], df["absorptance"], lw=1.4, label=label)
        ax.axvspan(8, 13, color="#E8EEF7", alpha=0.45)
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("Absorptance")
        ax.set_title(f"D07 oxide sensitivity spectra ({pol})")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, ncol=1)
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_mean_a(metrics: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    for ax, pol in zip(axes, ["Ez", "Hz"]):
        sub = metrics[metrics["polarization"] == pol]
        if sub.empty:
            ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
            ax.set_axis_off()
            continue
        for case, df in sub.groupby("case_name", sort=False):
            ax.plot(
                df["oxide_thickness_um"],
                df["mean_A_8_13um"],
                marker="o",
                label=case,
            )
        ax.set_title(pol)
        ax.set_xlabel("Oxide thickness (um)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Mean A, 8-13 um")
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        axes[-1].legend(handles, labels, fontsize=7, loc="best")
    fig.suptitle("D07 mean absorptance vs oxide thickness")
    save_figure(fig, out_path)
    plt.close(fig)


def _field_components(pol: str) -> tuple[list[int], list[int], list[int], int]:
    if pol == "Ez":
        return [mp.Ez], [mp.Dz], [mp.Hx, mp.Hy], mp.Ez
    if pol == "Hz":
        return [mp.Ex, mp.Ey], [mp.Dx, mp.Dy], [mp.Hz], mp.Hz
    raise ValueError(f"Unknown polarization: {pol}")


def _crop_common(arrays: list[np.ndarray]) -> list[np.ndarray]:
    min_shape = tuple(min(a.shape[i] for a in arrays) for i in range(arrays[0].ndim))
    slices = tuple(slice(0, n) for n in min_shape)
    return [np.asarray(a)[slices] for a in arrays]


def _sum_abs2(arrays: list[np.ndarray]) -> np.ndarray:
    cropped = _crop_common(arrays)
    out = np.zeros(cropped[0].shape, dtype=float)
    for arr in cropped:
        out += np.abs(arr) ** 2
    return out


def _absorbed_power_density(
    e_arrays: list[np.ndarray],
    d_arrays: list[np.ndarray],
    frequency: float,
) -> np.ndarray:
    e_crop, d_crop = _crop_common(e_arrays), _crop_common(d_arrays)
    arrays = _crop_common(e_crop + d_crop)
    n = len(arrays) // 2
    absorbed = np.zeros(arrays[0].shape, dtype=float)
    for e_arr, d_arr in zip(arrays[:n], arrays[n:]):
        absorbed += np.imag(np.conj(e_arr) * d_arr)
    return np.maximum(2.0 * np.pi * frequency * absorbed, 0.0)


def _plot_field_map(data: np.ndarray, title: str, out_path: Path,
                    extent: tuple[float, float, float, float]) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    vmax = np.nanpercentile(data, 99.5) if np.any(np.isfinite(data)) else 1.0
    if vmax <= 0:
        vmax = None
    im = ax.imshow(
        data.T,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="inferno",
        vmin=0,
        vmax=vmax,
    )
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    save_figure(fig, out_path)
    plt.close(fig)


def _run_best_field_snapshot(
    metrics: pd.DataFrame,
    args: argparse.Namespace,
    paths: dict[str, Path],
    logger,
) -> dict:
    if args.skip_best_field_snapshot or metrics.empty:
        return {"field_snapshot_status": "skipped"}

    best = metrics.sort_values("enhancement_over_bare", ascending=False).iloc[0]
    spec = _case_specs()[str(best["case_name"])]
    pol = str(best["polarization"])
    wl = float(best["peak_wavelength_um"])
    thickness = float(best["oxide_thickness_um"])
    ti = get_ti_medium(wl, wl)
    tio2 = get_tio2_medium(wl, wl, demo_index=args.tio2_demo_index)

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
    y_substrate_bottom = y_surface - args.substrate_thickness_um
    y_src = y_top_pml_inner - 0.25 * args.air_buffer_um
    y_refl = y_surface + 0.5 * args.air_buffer_um
    y_dft_min = y_substrate_bottom
    y_dft_max = min(y_refl, y_surface + 2.0)
    y_dft_center = 0.5 * (y_dft_min + y_dft_max)
    y_dft_size = y_dft_max - y_dft_min

    if spec.structure == "flat":
        geometry = _flat_geometry_factory(spec, thickness, ti, tio2, args)(
            y_surface, args.substrate_thickness_um,
        )
    else:
        geometry = _slanted_geometry_factory(spec, thickness, ti, tio2, args)(
            y_surface, args.substrate_thickness_um,
        )

    e_components, d_components, h_components, src_c = _field_components(pol)
    components = e_components + d_components + h_components
    fcen = 1.0 / wl
    fwidth = 0.08 * fcen
    sim = mp.Simulation(
        cell_size=mp.Vector3(args.period_um, cell_y, 0),
        boundary_layers=[mp.PML(thickness=args.pml_thickness_um, direction=mp.Y)],
        sources=[
            mp.Source(
                mp.GaussianSource(frequency=fcen, fwidth=fwidth, is_integrated=True),
                component=src_c,
                center=mp.Vector3(0, y_src, 0),
                size=mp.Vector3(args.period_um, 0, 0),
            )
        ],
        resolution=args.resolution,
        k_point=mp.Vector3(),
        geometry=geometry,
        dimensions=2,
    )
    dft = sim.add_dft_fields(
        components,
        fcen,
        0,
        1,
        center=mp.Vector3(0, y_dft_center, 0),
        size=mp.Vector3(args.period_um, y_dft_size, 0),
        yee_grid=True,
    )
    logger.info(
        ">>> best field snapshot case=%s pol=%s wl=%.4g oxide=%.4g",
        spec.case_name, pol, wl, thickness,
    )
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    e_arrays = [sim.get_dft_array(dft, comp, 0) for comp in e_components]
    d_arrays = [sim.get_dft_array(dft, comp, 0) for comp in d_components]
    h_arrays = [sim.get_dft_array(dft, comp, 0) for comp in h_components]
    e2 = _sum_abs2(e_arrays)
    h2 = _sum_abs2(h_arrays)
    absorbed = _absorbed_power_density(e_arrays, d_arrays, fcen)
    extent = (
        -0.5 * args.period_um,
        0.5 * args.period_um,
        y_dft_min,
        y_dft_max,
    )
    label = (
        f"{spec.case_name}, {pol}, {wl:g} um, oxide {thickness:g} um; "
        "TiO2 placeholder"
    )
    _plot_field_map(e2, "|E|^2 " + label, paths["best_e2"], extent)
    _plot_field_map(h2, "|H|^2 " + label, paths["best_h2"], extent)
    _plot_field_map(
        absorbed,
        "Absorbed power density " + label,
        paths["best_absorbed"],
        extent,
    )
    return {
        "field_snapshot_status": "completed",
        "field_snapshot_case": spec.case_name,
        "field_snapshot_pol": pol,
        "field_snapshot_wavelength_um": wl,
        "field_snapshot_oxide_thickness_um": thickness,
    }


def _write_report(paths: dict[str, Path], metrics: pd.DataFrame,
                  args: argparse.Namespace, field_info: dict) -> None:
    ensure_dir(paths["report"].parent)
    passed = bool((metrics["pass_or_fail"] == "PASS").all()) if not metrics.empty else False
    best = metrics.sort_values("enhancement_over_bare", ascending=False).head(1)
    if best.empty:
        best_text = "No completed metric rows."
        flat_vs_conformal = "No completed metric rows."
    else:
        row = best.iloc[0]
        best_text = (
            f"Best sensitivity enhancement: {row['case_name']} {row['polarization']} "
            f"t={row['oxide_thickness_um']:g} um, mean_A={row['mean_A_8_13um']:.4f}, "
            f"enhancement_over_bare={row['enhancement_over_bare']:.4f}."
        )
        flat = metrics[metrics["case_name"] == "oxide_flat_Ti"]
        conformal = metrics[metrics["case_name"] == "oxide_slanted_groove_conformal_approx"]
        if not flat.empty and not conformal.empty:
            flat_best = float(flat["enhancement_over_bare"].max())
            conf_best = float(conformal["enhancement_over_bare"].max())
            flat_vs_conformal = (
                f"Best flat-film enhancement = {flat_best:.4f}; best conformal "
                f"slanted-groove enhancement = {conf_best:.4f}."
            )
        else:
            flat_vs_conformal = "Flat/conformal comparison is incomplete in this run."

    report = f"""# D07 Oxide And Multiscale Structure Diagnostic

## Purpose

Test whether the weak emissivity enhancement of bare Ti slanted grooves may be
caused by missing oxide or composite surface effects.

## Source-Code Changes

- `src.materials.get_tio2_medium()` was added as a replaceable TiO2 interface.
  The current implementation is `placeholder_demo`: lossless, non-dispersive,
  and explicitly not quantitative.
- `src.geometry.build_oxidized_slanted_groove_geometry()` was added for
  `top_film_only` and `conformal_approx` oxide layouts.  `oxide_thickness_um=0`
  returns the existing bare slanted-groove geometry.

These changes are backward compatible with existing `get_ti_medium()`,
`build_slanted_groove_geometry()`, and `run_periodic_2d_metal_spectrum()`.

## Physical Assumptions

- Lengths are in um and Meep frequency is `f = 1 / wavelength_um`.
- Ti is Meep's Rakić Drude-Lorentz model.
- TiO2 uses placeholder_demo with `n={args.tio2_demo_index:g}`.
- The TiO2 placeholder is only a sensitivity model. It has no mid-IR phonon
  absorption and no validated wavelength range.
- A uses the D00/D01 flux convention:
  `R=refl/abs(input)`, `T=-trans/abs(input)`, `A=1-R-T`.

## Pass/Fail Criteria

- All R/T/A and raw flux columns are finite.
- Raw flux columns are preserved.
- Oxide thickness 0 degenerates to the corresponding bare result within
  `{args.degeneracy_tol:g}` in 8-13 um mean absorptance.

Overall status: **{'PASS' if passed else 'FAIL'}**

## Verified Numerical Findings

- {best_text}
- {flat_vs_conformal}

## Hypotheses, Not Final Quantitative Claims

- A higher placeholder-index oxide response would indicate that surface
  dielectric layers can tune optical coupling, not that real TiO2 gives that
  exact absorptance.
- If conformal groove oxide improves more than flat oxide in the placeholder
  sweep, the real sample should be checked for sidewall/bottom oxide coverage.
- The slant angle is still primarily assessed for directionality by D05/D06;
  this D07 normal-incidence sweep only tests total absorptance sensitivity.

## Needs Experiment Or Better Model

- Reliable mid-IR TiO2/TiOx optical constants or a Drude-Lorentz fit.
- Oxide thickness distribution after femtosecond laser processing.
- Cross-section geometry, roughness, and possible nanoparticle/recast layers.
- 3D morphology; the current model is a 2D periodic approximation.

## Required Answers

1. Does oxide significantly improve 8-13 um mean absorptance?
   {best_text}

2. Is improvement mainly from flat oxide or groove coverage?
   {flat_vs_conformal}

3. Does slant mainly affect absorption enhancement or directionality?
   In this script, slant is tested only through normal-incidence total
   absorptance. Directionality remains the domain of D05/D06.

4. What should experiments characterize first?
   Prioritize oxide thickness, oxide composition/phase and mid-IR n,k,
   cross-section groove geometry, then multiscale roughness.

## Output Files

- Spectra CSV: `{paths['spectra']}`
- Metrics CSV: `{paths['metrics']}`
- Geometry check: `{paths['geometry']}`
- Ez spectra: `{paths['spectra_ez']}`
- Hz spectra: `{paths['spectra_hz']}`
- Mean A plot: `{paths['mean_a']}`
- Best-case |E|^2: `{paths['best_e2']}`
- Best-case |H|^2: `{paths['best_h2']}`
- Best-case absorbed power: `{paths['best_absorbed']}`

Best-case field snapshot status: `{json.dumps(field_info, ensure_ascii=False)}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="D07 oxide sensitivity diagnostic.")
    p.add_argument(
        "--cases", nargs="+",
        default=[
            "bare_flat_Ti",
            "oxide_flat_Ti",
            "bare_slanted_groove",
            "oxide_slanted_groove_top_film_only",
            "oxide_slanted_groove_conformal_approx",
        ],
        choices=sorted(_case_specs().keys()),
    )
    p.add_argument("--polarizations", nargs="+", default=["Ez", "Hz"], choices=["Ez", "Hz"])
    p.add_argument("--oxide_thicknesses_um", nargs="*", default=None)
    p.add_argument("--wavelength_min_um", type=float, default=DEFAULTS["wavelength_min_um"])
    p.add_argument("--wavelength_max_um", type=float, default=DEFAULTS["wavelength_max_um"])
    p.add_argument("--nfreq", type=int, default=DEFAULTS["nfreq"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--period_um", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--top_width_um", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width_um", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--tilt_deg", type=float, default=DEFAULTS["tilt_deg"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float, default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--tio2_demo_index", type=float, default=DEFAULTS["tio2_demo_index"])
    p.add_argument("--enhancement_tol", type=float, default=DEFAULTS["enhancement_tol"])
    p.add_argument("--degeneracy_tol", type=float, default=DEFAULTS["degeneracy_tol"])
    p.add_argument("--geometry_check_thickness_um", type=float,
                   default=DEFAULTS["geometry_check_thickness_um"])
    p.add_argument("--skip_best_field_snapshot", action="store_true")
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    logger = setup_logger("D07_oxide_and_multiscale_structure_diagnostic")
    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)

    thicknesses = _parse_float_list(
        args.oxide_thicknesses_um,
        DEFAULTS["oxide_thicknesses_um"],
    )
    if any(t < 0 for t in thicknesses):
        raise ValueError("All oxide thicknesses must be >= 0")
    thicknesses = sorted(set(float(t) for t in thicknesses))

    case_map = _case_specs()
    spectra_parts = []
    for case_name in args.cases:
        spec = case_map[case_name]
        run_thicknesses = [0.0] if spec.oxide_mode == "none" else thicknesses
        for pol in args.polarizations:
            for thickness in run_thicknesses:
                spectra_parts.append(_run_one(spec, pol, thickness, args, logger))

    spectra = pd.concat(spectra_parts, ignore_index=True)
    metrics = _check_zero_thickness_degeneracy(_metrics(spectra, args), args)
    spectra.to_csv(paths["spectra"], index=False)
    metrics.to_csv(paths["metrics"], index=False)

    _plot_geometry_check(args, paths["geometry"])
    _plot_spectra(spectra, "Ez", paths["spectra_ez"])
    _plot_spectra(spectra, "Hz", paths["spectra_hz"])
    _plot_mean_a(metrics, paths["mean_a"])
    field_info = _run_best_field_snapshot(metrics, args, paths, logger)
    _write_report(paths, metrics, args, field_info)

    logger.info("Wrote %s", paths["spectra"])
    logger.info("Wrote %s", paths["metrics"])
    logger.info("Wrote %s", paths["report"])
    if (metrics["pass_or_fail"] == "FAIL").any():
        logger.warning("D07 completed with failed validation rows.")
    else:
        logger.info("D07 completed successfully.")


if __name__ == "__main__":
    main()
