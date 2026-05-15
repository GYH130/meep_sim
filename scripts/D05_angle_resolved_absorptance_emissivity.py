"""D05_angle_resolved_absorptance_emissivity.py — true angular absorptance.

This diagnostic computes A(lambda, theta, polarization) by running one
narrowband oblique-incidence simulation for every wavelength/angle pair.  Under
reciprocity and local thermal equilibrium, this absorptance is the directional
emissivity proxy for the same direction and polarization.

Physical assumptions
--------------------
1. Length unit is um; Meep frequency is f = 1 / wavelength_um.
2. The surface normal is y and the periodic direction is x.
3. For each fixed wavelength and incident angle, kx = f*sin(theta).  The Bloch
   boundary uses k_point=(kx,0,0), and the extended source uses
   amp_func(r)=exp(i*2*pi*kx*x), following Meep's standard oblique plane-wave
   construction.  y is terminated by PML, so ky is not a Bloch component.
4. A uses the D00/D01-validated flux convention: R=refl/|input|,
   T=-trans/|input|, A=1-R-T.
5. scripts/04_angle_resolved_emission.py is a normal-incidence scattering
   direction plot; it is not the same quantity as this script computes.

Pass/fail criteria
------------------
1. All R/T/A and raw flux values must be finite.
2. flat_Ti and symmetric_groove should satisfy A(+theta) ~= A(-theta) within
   symmetry_tol.
3. theta=0 results are compared to D02 where available and reported as warnings
   when the difference exceeds theta0_d02_tol.
"""

from __future__ import annotations

import argparse
import cmath
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

from src.geometry import build_slanted_groove_geometry
from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium


DEFAULTS = dict(
    wavelengths_um=[8.0, 10.0, 12.0],
    full_wavelengths_um=[8.0, 9.0, 10.0, 11.0, 12.0, 13.0],
    angle_min_deg=-70.0,
    angle_max_deg=70.0,
    angle_step_deg=10.0,
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    resolution=32,
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    fwidth_fraction=0.06,
    decay_db=30.0,
    theta0_deg=30.0,
    symmetry_tol=0.03,
    theta0_d02_tol=0.04,
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
            "results", "diagnostics", "tables",
            "D05_angle_resolved_absorptance.csv",
        ),
        "metrics": project_path(
            "results", "diagnostics", "tables", "D05_directionality_metrics.csv",
        ),
        "heatmap_ez": project_path(
            "results", "diagnostics", "figures", "D05_absorptance_heatmap_Ez.png",
        ),
        "heatmap_hz": project_path(
            "results", "diagnostics", "figures", "D05_absorptance_heatmap_Hz.png",
        ),
        "directionality": project_path(
            "results", "diagnostics", "figures",
            "D05_band_averaged_directionality.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D05_angle_resolved_emissivity_report.md",
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
        "mirrored_slanted_groove": CaseSpec(
            "mirrored_slanted_groove", "Mirrored slanted groove, tilt=-20",
            "groove", -20.0,
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


def _angle_values(args: argparse.Namespace) -> list[float]:
    n = int(round((args.angle_max_deg - args.angle_min_deg) / args.angle_step_deg))
    vals = [args.angle_min_deg + i * args.angle_step_deg for i in range(n + 1)]
    if not any(abs(v) < 1e-9 for v in vals):
        vals.append(0.0)
    return sorted(set(round(v, 10) for v in vals))


def _cell_layout(args: argparse.Namespace) -> dict[str, float]:
    bottom_buffer = args.pml_thickness_um
    cell_y = (
        2 * args.pml_thickness_um
        + args.air_buffer_um
        + args.substrate_thickness_um
        + bottom_buffer
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
    return dict(
        cell_y=cell_y,
        y_surface=y_surface,
        y_src=y_src,
        y_refl=y_refl,
        y_trans=y_trans,
    )


def _geometry(case: CaseSpec, args: argparse.Namespace, ti, y_surface: float) -> list:
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


def _source_component(pol: str):
    if pol == "Ez":
        return mp.Ez
    if pol == "Hz":
        return mp.Hz
    raise ValueError(f"Unknown polarization: {pol}")


def _run_one(case: CaseSpec, pol: str, wl: float, theta: float,
             args: argparse.Namespace, logger) -> dict:
    ti = get_ti_medium(lambda_min_um=wl, lambda_max_um=wl)
    layout = _cell_layout(args)
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
    logger.info(
        ">>> case=%s pol=%s wl=%.4g theta=%.3g kx=%.6g",
        case.case_name, pol, wl, theta, kx,
    )

    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        k_point=k_point,
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
            c=src_c,
            pt=mp.Vector3(0, layout["y_refl"], 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    input_flux_raw = float(np.array(mp.get_fluxes(refl_ref))[0])
    ref_data = sim_ref.get_flux_data(refl_ref)
    ref_time = time.time() - t0
    if input_flux_raw == 0:
        raise RuntimeError("Reference input flux is zero")

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        k_point=k_point,
        geometry=_geometry(case, args, ti, layout["y_surface"]),
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
    t0 = time.time()
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, layout["y_refl"], 0),
            decay_by=10 ** (-args.decay_db / 10.0),
        )
    )
    struct_time = time.time() - t0
    refl_raw = float(np.array(mp.get_fluxes(refl))[0])
    trans_raw = float(np.array(mp.get_fluxes(trans))[0])

    R = refl_raw / abs(input_flux_raw)
    T_raw = trans_raw / abs(input_flux_raw)
    T = -trans_raw / abs(input_flux_raw)
    A = 1.0 - R - T
    finite = bool(np.all(np.isfinite([R, T, A, input_flux_raw, refl_raw, trans_raw])))
    return dict(
        case_name=case.case_name,
        case_label=case.case_label,
        polarization=pol,
        wavelength_um=wl,
        frequency_meep=fcen,
        theta_deg=theta,
        kx_meep=kx,
        input_flux_raw=input_flux_raw,
        reflection_flux_raw=refl_raw,
        transmission_flux_raw=trans_raw,
        reflectance_total=R,
        raw_transmittance=T_raw,
        transmittance=T,
        absorptance=A,
        emissivity_proxy=A,
        all_finite=finite,
        pass_or_fail="PASS" if finite else "FAIL",
        normalization_note=(
            "R = reflection_flux_raw / abs(input_flux_raw); "
            "T = -transmission_flux_raw / abs(input_flux_raw); A = 1 - R - T"
        ),
        angle_note=(
            "Separate narrowband simulation for each wavelength/theta; "
            "kx = frequency*sin(theta), amp_func = exp(i*2*pi*kx*x)"
        ),
        kirchhoff_note=(
            "emissivity_proxy=A only under reciprocity and local thermal equilibrium "
            "for the same direction and polarization"
        ),
        walltime_s=ref_time + struct_time,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        air_buffer_um=args.air_buffer_um,
        substrate_thickness_um=args.substrate_thickness_um,
        decay_db=args.decay_db,
    )


def _interp_theta(df: pd.DataFrame, theta: float) -> float:
    sub = df.sort_values("theta_deg")
    xs = sub["theta_deg"].to_numpy(dtype=float)
    ys = sub["absorptance"].to_numpy(dtype=float)
    if theta < xs.min() or theta > xs.max():
        return np.nan
    return float(np.interp(theta, xs, ys))


def _d02_delta(case_name: str, pol: str, wl: float, a0: float) -> float:
    path = project_path(
        "results", "diagnostics", "tables", "D02_polarization_spectra.csv",
    )
    if not path.is_file():
        return np.nan
    mapping = {
        "flat_Ti": "flat_ti",
        "symmetric_groove": "symmetric_groove_tilt0",
        "slanted_groove": "slanted_groove_tilt20",
        "mirrored_slanted_groove": None,
    }
    d02_case = mapping.get(case_name)
    if d02_case is None:
        return np.nan
    d02 = pd.read_csv(path)
    sub = d02[(d02["case_name"] == d02_case) & (d02["polarization"] == pol)]
    if sub.empty:
        return np.nan
    sub = sub.sort_values("wavelength_um")
    xs = sub["wavelength_um"].to_numpy(dtype=float)
    ys = sub["absorptance"].to_numpy(dtype=float)
    if wl < xs.min() or wl > xs.max():
        return np.nan
    return float(a0 - np.interp(wl, xs, ys))


def _metrics(spectra: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    for key, df in spectra.groupby(
        ["case_name", "case_label", "polarization", "wavelength_um"], sort=False
    ):
        case_name, case_label, pol, wl = key
        idx = df["absorptance"].idxmax()
        a_pos = _interp_theta(df, args.theta0_deg)
        a_neg = _interp_theta(df, -args.theta0_deg)
        a0 = _interp_theta(df, 0.0)
        rows.append(dict(
            metric_scope="per_wavelength",
            case_name=case_name,
            case_label=case_label,
            polarization=pol,
            wavelength_um=wl,
            theta_of_max_absorptance=float(df.loc[idx, "theta_deg"]),
            max_absorptance=float(df.loc[idx, "absorptance"]),
            theta0_deg=args.theta0_deg,
            absorptance_plus_theta0=a_pos,
            absorptance_minus_theta0=a_neg,
            asymmetry_at_theta0=a_pos - a_neg if np.all(np.isfinite([a_pos, a_neg])) else np.nan,
            ratio_at_theta0=a_pos / a_neg if np.isfinite(a_neg) and abs(a_neg) > 1e-12 else np.nan,
            theta0_absorptance=a0,
            theta0_delta_vs_D02=_d02_delta(case_name, pol, wl, a0),
            mean_A_theta=np.nan,
            mean_directionality_ratio=np.nan,
            validation_note="per wavelength",
        ))

    for key, df in spectra.groupby(["case_name", "case_label", "polarization"], sort=False):
        case_name, case_label, pol = key
        mean_by_theta = df.groupby("theta_deg", as_index=False)["absorptance"].mean()
        idx = mean_by_theta["absorptance"].idxmax()
        a_pos = _interp_theta(mean_by_theta, args.theta0_deg)
        a_neg = _interp_theta(mean_by_theta, -args.theta0_deg)
        rows.append(dict(
            metric_scope="band_average",
            case_name=case_name,
            case_label=case_label,
            polarization=pol,
            wavelength_um=np.nan,
            theta_of_max_absorptance=float(mean_by_theta.loc[idx, "theta_deg"]),
            max_absorptance=float(mean_by_theta.loc[idx, "absorptance"]),
            theta0_deg=args.theta0_deg,
            absorptance_plus_theta0=a_pos,
            absorptance_minus_theta0=a_neg,
            asymmetry_at_theta0=a_pos - a_neg if np.all(np.isfinite([a_pos, a_neg])) else np.nan,
            ratio_at_theta0=a_pos / a_neg if np.isfinite(a_neg) and abs(a_neg) > 1e-12 else np.nan,
            theta0_absorptance=_interp_theta(mean_by_theta, 0.0),
            theta0_delta_vs_D02=np.nan,
            mean_A_theta=float(mean_by_theta["absorptance"].mean()),
            mean_directionality_ratio=a_pos / a_neg if np.isfinite(a_neg) and abs(a_neg) > 1e-12 else np.nan,
            validation_note="band average over simulated wavelengths",
        ))
    return pd.DataFrame(rows)


def _symmetry_residuals(spectra: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, df in spectra.groupby(["case_name", "polarization", "wavelength_um"]):
        vals = df.set_index("theta_deg")["absorptance"].to_dict()
        diffs = []
        for theta, a_pos in vals.items():
            if theta > 0 and -theta in vals:
                diffs.append(abs(a_pos - vals[-theta]))
        rows.append(dict(
            case_name=key[0],
            polarization=key[1],
            wavelength_um=key[2],
            max_symmetry_residual=max(diffs) if diffs else np.nan,
        ))
    return pd.DataFrame(rows)


def _mirror_residuals(spectra: pd.DataFrame) -> pd.DataFrame:
    if "mirrored_slanted_groove" not in set(spectra["case_name"]):
        return pd.DataFrame()
    rows = []
    slanted = spectra[spectra["case_name"] == "slanted_groove"]
    mirrored = spectra[spectra["case_name"] == "mirrored_slanted_groove"]
    for (pol, wl), df in slanted.groupby(["polarization", "wavelength_um"]):
        mirror = mirrored[
            (mirrored["polarization"] == pol) & (mirrored["wavelength_um"] == wl)
        ]
        if mirror.empty:
            continue
        vals = df.set_index("theta_deg")["absorptance"].to_dict()
        mvals = mirror.set_index("theta_deg")["absorptance"].to_dict()
        diffs = [abs(a - mvals[-theta]) for theta, a in vals.items() if -theta in mvals]
        rows.append(dict(
            polarization=pol,
            wavelength_um=wl,
            max_mirror_residual=max(diffs) if diffs else np.nan,
        ))
    return pd.DataFrame(rows)


def _plot_heatmap(spectra: pd.DataFrame, pol: str, path: Path) -> None:
    sub = spectra[spectra["polarization"] == pol]
    if sub.empty:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, f"No {pol} data", ha="center", va="center")
        ax.axis("off")
        save_figure(fig, path)
        plt.close(fig)
        return
    cases = list(sub["case_name"].drop_duplicates())
    fig, axes = plt.subplots(1, len(cases), figsize=(4.8 * len(cases), 4.2), squeeze=False)
    vmax = max(0.25, float(sub["absorptance"].max()))
    for ax, case_name in zip(axes[0], cases):
        df = sub[sub["case_name"] == case_name]
        pivot = df.pivot_table(
            index="theta_deg", columns="wavelength_um", values="absorptance",
            aggfunc="mean",
        ).sort_index()
        x_min = float(pivot.columns.min())
        x_max = float(pivot.columns.max())
        y_min = float(pivot.index.min())
        y_max = float(pivot.index.max())
        if x_min == x_max:
            x_min -= 0.5
            x_max += 0.5
        if y_min == y_max:
            y_min -= 0.5
            y_max += 0.5
        im = ax.imshow(
            pivot.to_numpy(),
            origin="lower",
            aspect="auto",
            extent=[x_min, x_max, y_min, y_max],
            cmap="magma",
            vmin=0,
            vmax=vmax,
        )
        ax.set_title(f"{case_name} {pol}")
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("Incident angle (deg)")
        fig.colorbar(im, ax=ax, label="A")
    fig.suptitle(f"D05 angle-resolved absorptance, {pol}")
    save_figure(fig, path)
    plt.close(fig)


def _plot_directionality(spectra: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, pol in zip(axes, ["Ez", "Hz"]):
        sub = spectra[spectra["polarization"] == pol]
        for case_name, df in sub.groupby("case_name", sort=False):
            mean_by_theta = df.groupby("theta_deg", as_index=False)["absorptance"].mean()
            ax.plot(mean_by_theta["theta_deg"], mean_by_theta["absorptance"],
                    marker="o", ms=3, lw=1.5, label=case_name)
        ax.axvline(0, color="0.4", ls=":", lw=1)
        ax.set_title(pol)
        ax.set_xlabel("Incident angle (deg)")
        ax.grid(True, ls=":", alpha=0.5)
    axes[0].set_ylabel("Mean A over simulated wavelengths")
    handles, labels = [], []
    for ax in axes:
        h, lab = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(lab)
    if handles:
        axes[-1].legend(handles, labels, fontsize=8)
    fig.suptitle("D05 band-averaged directionality")
    save_figure(fig, path)
    plt.close(fig)


def _md_table(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows)
    cols = list(view.columns)
    rows = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in view.iterrows():
        vals = []
        for col in cols:
            val = row[col]
            vals.append(f"{val:.6g}" if isinstance(val, float) else str(val))
        rows.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        rows.append(
            "| " + " | ".join(["..."] + [f"{len(df) - max_rows} more rows"] + [""] * (len(cols) - 2)) + " |"
        )
    return "\n".join(rows)


def _write_report(spectra: pd.DataFrame, metrics: pd.DataFrame,
                  args: argparse.Namespace) -> Path:
    paths = _paths()
    sym = _symmetry_residuals(spectra)
    mirror = _mirror_residuals(spectra)
    band = metrics[metrics["metric_scope"] == "band_average"]
    flat_mean = band[band["case_name"] == "flat_Ti"]["mean_A_theta"].mean()
    sym_mean = band[band["case_name"] == "symmetric_groove"]["mean_A_theta"].mean()
    sl_mean = band[band["case_name"] == "slanted_groove"]["mean_A_theta"].mean()
    total_boost = np.isfinite(flat_mean) and np.isfinite(sl_mean) and (sl_mean - flat_mean > 0.05)

    sl_band = band[band["case_name"] == "slanted_groove"]
    if sl_band.empty:
        max_asym = np.nan
        main_dir = "Not enough slanted_groove data."
        direction_mod = False
    else:
        row = sl_band.loc[sl_band["asymmetry_at_theta0"].abs().idxmax()]
        max_asym = float(abs(row["asymmetry_at_theta0"]))
        direction_mod = max_asym > args.symmetry_tol
        main_dir = (
            f"polarization={row['polarization']}, "
            f"theta_max={row['theta_of_max_absorptance']:.3g} deg, "
            f"A(+/-theta0) ratio={row['ratio_at_theta0']:.3g}"
        )

    d02_warn = metrics[metrics["theta0_delta_vs_D02"].abs() > args.theta0_d02_tol]
    control_fail = sym[
        sym["case_name"].isin(["flat_Ti", "symmetric_groove"])
        & (sym["max_symmetry_residual"] > args.symmetry_tol)
    ]
    all_finite = bool((spectra["pass_or_fail"] == "PASS").all())

    lines = [
        "# D05 Angle-Resolved Absorptance / Emissivity Proxy",
        "",
        "## Purpose",
        "Compute true angular absorptance by separate oblique-incidence simulations.",
        "",
        "## Physical assumptions",
        "- Length unit: um; Meep frequency: f = 1 / wavelength_um.",
        "- Surface normal is y and periodic direction is x.",
        "- kx = frequency*sin(theta); source phase is exp(i*2*pi*kx*x).",
        "- Each wavelength and angle is simulated separately, not via one broadband fixed-k run.",
        "- emissivity_proxy=A only under reciprocity, local thermal equilibrium, and same direction/polarization Kirchhoff correspondence.",
        f"- Ti Rakić validity upper wavelength is {TI_RAKIC_VALID_LAMBDA_UM[1]} um; 13 um is extrapolated.",
        "",
        "## Pass/fail criteria",
        "- All raw flux and R/T/A values finite.",
        f"- flat_Ti and symmetric_groove symmetry residual <= {args.symmetry_tol:g}.",
        f"- theta=0 D02 comparison warning threshold <= {args.theta0_d02_tol:g}.",
        "",
        "## Required answers",
        f"1. Does the current slanted groove improve total emissivity? {'Yes in this run' if total_boost else 'Not clearly'} "
        f"(flat mean={flat_mean:.4g}, symmetric mean={sym_mean:.4g}, slanted mean={sl_mean:.4g}).",
        f"2. Does it show directional asymmetry? {'Yes' if direction_mod else 'Not clearly'} "
        f"(largest band-averaged |A(+theta0)-A(-theta0)|={max_asym:.4g}).",
        f"3. Main directionality: {main_dir}",
        f"4. Current structure looks more like {'a direction-modulating structure' if direction_mod and not total_boost else 'a low/moderate-emissivity diagnostic structure'}.",
        "",
        "## Verified Numerical Conclusions",
        f"- Finite flux status: {'PASS' if all_finite else 'FAIL'}.",
        f"- Control symmetry failures: {len(control_fail)}.",
        f"- D02 theta=0 warnings: {len(d02_warn)}.",
        "",
        "## Band-Averaged Metrics",
        _md_table(band, max_rows=20),
        "",
        "## Symmetry Residuals",
        _md_table(sym, max_rows=30),
        "",
        "## Mirror Residuals",
        _md_table(mirror, max_rows=30),
        "",
        "## Hypotheses",
        "- If slanted_groove asymmetry exceeds flat/symmetric controls while mean A stays low, the geometry mainly modulates direction rather than raising total emissivity.",
        "- Stronger high-emissivity behavior may require oxide layers, roughness, multiscale features, or different period/depth.",
        "",
        "## Needs Higher-Fidelity Confirmation",
        "- Full 5-degree angle grid, 8-13 um wavelength grid, converged Hz resolution, 3D geometry, oxide/roughness, measured Ti optical constants.",
        "",
        "## Output Files",
        f"- `{paths['spectra']}`",
        f"- `{paths['metrics']}`",
        f"- `{paths['heatmap_ez']}`",
        f"- `{paths['heatmap_hz']}`",
        f"- `{paths['directionality']}`",
        "",
        "## Run Configuration",
        "```json",
        json.dumps(vars(args), indent=2, ensure_ascii=False),
        "```",
    ]
    ensure_dir(paths["report"].parent)
    paths["report"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths["report"]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="D05 angular absorptance/emissivity proxy.")
    p.add_argument("--wavelengths_um", nargs="*", default=None)
    p.add_argument("--full_wavelength_grid", action="store_true")
    p.add_argument("--angle_min_deg", type=float, default=DEFAULTS["angle_min_deg"])
    p.add_argument("--angle_max_deg", type=float, default=DEFAULTS["angle_max_deg"])
    p.add_argument("--angle_step_deg", type=float, default=DEFAULTS["angle_step_deg"])
    p.add_argument(
        "--cases", nargs="+",
        default=["flat_Ti", "symmetric_groove", "slanted_groove"],
        choices=list(_case_specs().keys()),
    )
    p.add_argument("--include_mirror", action="store_true")
    p.add_argument("--polarizations", nargs="+", default=["Ez", "Hz"], choices=["Ez", "Hz"])
    p.add_argument("--period_um", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--top_width_um", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width_um", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float, default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--fwidth_fraction", type=float, default=DEFAULTS["fwidth_fraction"])
    p.add_argument("--decay_db", type=float, default=DEFAULTS["decay_db"])
    p.add_argument("--theta0_deg", type=float, default=DEFAULTS["theta0_deg"])
    p.add_argument("--symmetry_tol", type=float, default=DEFAULTS["symmetry_tol"])
    p.add_argument("--theta0_d02_tol", type=float, default=DEFAULTS["theta0_d02_tol"])
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    logger = setup_logger("D05_angle_resolved_absorptance_emissivity")
    logger.info("=== D05_angle_resolved_absorptance_emissivity ===")
    logger.info("args=%s", vars(args))

    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)

    default_wls = DEFAULTS["full_wavelengths_um"] if args.full_wavelength_grid else DEFAULTS["wavelengths_um"]
    wavelengths = _parse_float_list(args.wavelengths_um, default_wls)
    angles = _angle_values(args)
    cases = list(args.cases)
    if args.include_mirror and "mirrored_slanted_groove" not in cases:
        cases.append("mirrored_slanted_groove")
    specs = _case_specs()
    logger.info(
        "cases=%s polarizations=%s wavelengths=%s angles=%s total=%d",
        cases, args.polarizations, wavelengths, angles,
        len(cases) * len(args.polarizations) * len(wavelengths) * len(angles),
    )

    rows = []
    for case_name in cases:
        for pol in args.polarizations:
            for wl in wavelengths:
                for theta in angles:
                    rows.append(_run_one(specs[case_name], pol, wl, theta, args, logger))

    spectra = pd.DataFrame(rows)
    spectra.to_csv(paths["spectra"], index=False)
    metrics = _metrics(spectra, args)
    metrics.to_csv(paths["metrics"], index=False)
    _plot_heatmap(spectra, "Ez", paths["heatmap_ez"])
    _plot_heatmap(spectra, "Hz", paths["heatmap_hz"])
    _plot_directionality(spectra, paths["directionality"])
    report = _write_report(spectra, metrics, args)

    logger.info("saved spectra: %s", paths["spectra"])
    logger.info("saved metrics: %s", paths["metrics"])
    logger.info("saved report: %s", report)
    all_pass = bool((spectra["pass_or_fail"] == "PASS").all())
    logger.info("D05 %s", "PASS" if all_pass else "FAIL")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
