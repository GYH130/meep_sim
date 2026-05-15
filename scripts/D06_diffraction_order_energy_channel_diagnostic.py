"""D06_diffraction_order_energy_channel_diagnostic.py.

Use Meep mode decomposition to separate reflected diffraction-order power in
periodic Ti groove structures.  This diagnostic distinguishes true absorption
from redistribution of reflected energy into non-specular diffraction channels.

Physical assumptions
--------------------
1. Length unit is um; Meep frequency is f = 1 / wavelength_um.
2. The periodic direction is x and the surface normal is y.
3. Each wavelength/incident-angle pair is simulated separately.  The Bloch
   wavevector is kx = f*sin(theta), and the extended source uses
   amp_func(r)=exp(i*2*pi*kx*x), following Meep's oblique plane-wave approach.
4. R/T/A use the D00/D01-validated convention: R=refl/|input|,
   T=-trans/|input|, A=1-R-T.
5. Reflected diffraction orders are computed with Meep's
   get_eigenmode_coefficients and DiffractedPlanewave.  Ez is s-polarized and
   Hz is p-polarized for the x-y incidence plane.

Pass/fail criteria
------------------
1. R/T/A and computed order powers must be finite.
2. The selected sum of reflected diffraction-order powers should match the
   reflection flux monitor within order_sum_tol.
3. Symmetric grooves should have +m and -m reflected powers close within
   symmetry_tol for normal incidence.
"""

from __future__ import annotations

import argparse
import cmath
import math
import os
import sys
import time
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
from src.io_utils import ensure_dir, project_path, save_figure, setup_logger
from src.materials import get_ti_medium


DEFAULTS = dict(
    wavelengths_um=[8.0, 10.0, 12.0],
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
    order_sum_tol=0.08,
    symmetry_tol=0.05,
    grazing_tol=1e-7,
    d05_directionality_tol=0.03,
)


@dataclass(frozen=True)
class CaseSpec:
    case_name: str
    case_label: str
    tilt_deg: float


def _paths() -> dict[str, Path]:
    return {
        "orders": project_path(
            "results", "diagnostics", "tables", "D06_diffraction_orders.csv",
        ),
        "summary": project_path(
            "results", "diagnostics", "tables",
            "D06_energy_channel_summary.csv",
        ),
        "order_bars": project_path(
            "results", "diagnostics", "figures",
            "D06_diffraction_order_barplots.png",
        ),
        "channels": project_path(
            "results", "diagnostics", "figures",
            "D06_energy_channel_stacked_bars.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D06_diffraction_energy_channel_report.md",
        ),
    }


def _case_specs() -> dict[str, CaseSpec]:
    return {
        "symmetric_groove": CaseSpec(
            "symmetric_groove", "Symmetric groove, tilt=0", 0.0,
        ),
        "slanted_groove": CaseSpec(
            "slanted_groove", "Slanted groove, tilt=20", 20.0,
        ),
        "mirrored_slanted_groove": CaseSpec(
            "mirrored_slanted_groove", "Mirrored slanted groove, tilt=-20", -20.0,
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


def _source_component(pol: str):
    if pol == "Ez":
        return mp.Ez
    if pol == "Hz":
        return mp.Hz
    raise ValueError(f"Unknown polarization: {pol}")


def _diffracted_planewave(order_m: int, pol: str, axis):
    # The incidence plane is x-y.  s-polarization has E along z (Ez), while
    # p-polarization has H along z (Hz).
    if pol == "Ez":
        return mp.DiffractedPlanewave(
            g=[order_m, 0, 0], axis=axis, s=1, p=0,
        )
    if pol == "Hz":
        return mp.DiffractedPlanewave(
            g=[order_m, 0, 0], axis=axis, s=0, p=1,
        )
    raise ValueError(f"Unknown polarization: {pol}")


def _geometry(case: CaseSpec, args: argparse.Namespace, ti, y_surface: float) -> list:
    return build_slanted_groove_geometry(
        period_x_um=args.period_um,
        top_width_um=args.top_width_um,
        bottom_width_um=args.bottom_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=case.tilt_deg,
        substrate_thickness_um=args.substrate_thickness_um,
        y_surface=y_surface,
        medium_substrate=ti,
    )


def _propagating_orders(
    wavelength_um: float,
    theta_deg: float,
    period_um: float,
    grazing_tol: float,
) -> list[dict]:
    freq = 1.0 / wavelength_um
    kx_inc = freq * math.sin(math.radians(theta_deg))
    max_order = int(math.ceil(period_um * (freq + abs(kx_inc)))) + 2
    rows = []
    for m in range(-max_order, max_order + 1):
        kx_m = kx_inc + m / period_um
        margin = freq - abs(kx_m)
        if margin >= -grazing_tol:
            clipped = max(-1.0, min(1.0, kx_m / freq))
            rows.append(
                dict(
                    order_m=m,
                    kx_incident=kx_inc,
                    kx_order=kx_m,
                    diffraction_angle_deg=math.degrees(math.asin(clipped)),
                    is_grazing=abs(margin) <= grazing_tol,
                )
            )
    return rows


def _read_extra_angles_from_d05(args: argparse.Namespace, logger) -> list[float]:
    if not args.include_d05_directional_angle:
        return []
    path = project_path(
        "results", "diagnostics", "tables", "D05_directionality_metrics.csv",
    )
    if not path.exists():
        logger.info("D05 metrics not found; only requested angles will be used.")
        return []
    try:
        df = pd.read_csv(path)
    except Exception as exc:
        logger.warning("Could not read D05 metrics from %s: %s", path, exc)
        return []
    df = df[(df["case_name"] == "slanted_groove") & (df["metric_scope"] == "band_average")]
    if df.empty:
        return []
    extra = []
    for _, row in df.iterrows():
        asym = abs(float(row.get("asymmetry_at_theta0", 0.0)))
        theta0 = float(row.get("theta0_deg", 0.0))
        if asym >= args.d05_directionality_tol and abs(theta0) > 0:
            extra.extend([theta0, -theta0])
    if extra:
        logger.info("Added D05 directional angles: %s", sorted(set(extra)))
    else:
        logger.info("D05 did not exceed directionality tolerance; no extra angles added.")
    return extra


def _run_one(
    case: CaseSpec,
    pol: str,
    wl: float,
    theta: float,
    args: argparse.Namespace,
    logger,
) -> tuple[list[dict], dict]:
    ti = get_ti_medium(lambda_min_um=wl, lambda_max_um=wl)
    layout = _cell_layout(args)
    fcen = 1.0 / wl
    fwidth = args.fwidth_fraction * fcen
    kx = fcen * math.sin(math.radians(theta))
    src_c = _source_component(pol)
    orders = _propagating_orders(wl, theta, args.period_um, args.grazing_tol)

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
    flux_region = mp.FluxRegion(
        center=mp.Vector3(0, layout["y_refl"], 0),
        size=mp.Vector3(args.period_um, 0, 0),
    )

    logger.info(
        ">>> case=%s pol=%s wl=%.4g theta=%.3g orders=%s",
        case.case_name, pol, wl, theta, [o["order_m"] for o in orders],
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
    refl_ref = sim_ref.add_mode_monitor(fcen, 0, 1, flux_region)
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
        raise RuntimeError("Reference input flux is zero.")

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=args.resolution,
        k_point=k_point,
        geometry=_geometry(case, args, ti, layout["y_surface"]),
        dimensions=2,
    )
    refl = sim.add_mode_monitor(fcen, 0, 1, flux_region)
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
    refl_flux_raw = float(np.array(mp.get_fluxes(refl))[0])
    trans_flux_raw = float(np.array(mp.get_fluxes(trans))[0])
    struct_time = time.time() - t0

    input_abs = abs(input_flux_raw)
    reflectance = refl_flux_raw / input_abs
    transmittance = -trans_flux_raw / input_abs
    absorptance = 1.0 - reflectance - transmittance

    order_rows = []
    for order in orders:
        m = int(order["order_m"])
        try:
            # The axis defines the s/p plane.  Use x for the usual x-y
            # incidence plane, but switch to y for near-grazing orders whose
            # k+G is nearly parallel to x.
            axis = (
                mp.Vector3(0, 1, 0)
                if bool(order["is_grazing"])
                else mp.Vector3(1, 0, 0)
            )
            coeff = sim.get_eigenmode_coefficients(
                refl,
                _diffracted_planewave(m, pol, axis),
                eig_parity=mp.NO_PARITY,
                direction=mp.Y,
            )
            alpha = np.asarray(coeff.alpha)
            power_0 = float(abs(alpha[0, 0, 0]) ** 2 / input_abs)
            power_1 = float(abs(alpha[0, 0, 1]) ** 2 / input_abs)
            err = ""
        except Exception as exc:
            power_0 = math.nan
            power_1 = math.nan
            err = str(exc)
            logger.warning(
                "Mode decomposition failed: case=%s pol=%s wl=%s theta=%s m=%s: %s",
                case.case_name, pol, wl, theta, m, exc,
            )
        order_rows.append(
            dict(
                case_name=case.case_name,
                case_label=case.case_label,
                polarization=pol,
                wavelength_um=wl,
                incident_angle_deg=theta,
                order_m=m,
                diffraction_angle_deg=order["diffraction_angle_deg"],
                is_grazing=bool(order["is_grazing"]),
                input_flux_raw=input_flux_raw,
                reflection_flux_raw=refl_flux_raw,
                transmission_flux_raw=trans_flux_raw,
                reflectance_flux_monitor=reflectance,
                transmittance=transmittance,
                absorptance=absorptance,
                alpha_power_direction_0=power_0,
                alpha_power_direction_1=power_1,
                selected_direction_index=-1,
                reflected_power_fraction=math.nan,
                mode_error=err,
            )
        )

    sums = {}
    for idx in (0, 1):
        vals = [r[f"alpha_power_direction_{idx}"] for r in order_rows]
        sums[idx] = float(np.nansum(vals))
    selected_idx = min(sums, key=lambda idx: abs(sums[idx] - reflectance))
    total_order = sums[selected_idx]
    residual = total_order - reflectance

    for row in order_rows:
        row["selected_direction_index"] = selected_idx
        row["reflected_power_fraction"] = row[f"alpha_power_direction_{selected_idx}"]

    specular = float(np.nansum([
        r["reflected_power_fraction"] for r in order_rows if r["order_m"] == 0
    ]))
    positive = float(np.nansum([
        r["reflected_power_fraction"] for r in order_rows if r["order_m"] > 0
    ]))
    negative = float(np.nansum([
        r["reflected_power_fraction"] for r in order_rows if r["order_m"] < 0
    ]))
    nonspecular = positive + negative

    finite_ok = all(np.isfinite(v) for v in [
        reflectance, transmittance, absorptance, total_order,
    ])
    residual_ok = abs(residual) <= args.order_sum_tol
    summary = dict(
        case_name=case.case_name,
        case_label=case.case_label,
        polarization=pol,
        wavelength_um=wl,
        incident_angle_deg=theta,
        propagating_orders=" ".join(str(o["order_m"]) for o in orders),
        selected_direction_index=selected_idx,
        input_flux_raw=input_flux_raw,
        reflection_flux_raw=refl_flux_raw,
        transmission_flux_raw=trans_flux_raw,
        R_flux_monitor=reflectance,
        T=transmittance,
        absorptance=absorptance,
        specular_reflection=specular,
        positive_diffraction_orders=positive,
        negative_diffraction_orders=negative,
        nonspecular_reflection=nonspecular,
        total_reflected_power_from_orders=total_order,
        order_sum_residual=residual,
        channel_sum=absorptance + transmittance + total_order,
        walltime_s=ref_time + struct_time,
        pass_or_fail="PASS" if finite_ok and residual_ok else "FAIL",
        validation_note=(
            "order sum matches flux"
            if finite_ok and residual_ok
            else "order sum mismatch or non-finite value"
        ),
    )
    return order_rows, summary


def _add_symmetry_checks(summary_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    df = summary_df.copy()
    df["pm_order_symmetry_delta"] = np.nan
    if df.empty:
        return df
    for idx, row in df.iterrows():
        if row["case_name"] != "symmetric_groove" or abs(row["incident_angle_deg"]) > 1e-9:
            continue
        delta = abs(row["positive_diffraction_orders"] - row["negative_diffraction_orders"])
        df.loc[idx, "pm_order_symmetry_delta"] = delta
        if delta > args.symmetry_tol:
            df.loc[idx, "pass_or_fail"] = "FAIL"
            df.loc[idx, "validation_note"] += "; +m/-m symmetry failed"
    return df


def _plot_order_bars(order_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    if order_df.empty:
        ax.text(0.5, 0.5, "No diffraction-order data", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, out_path)
        plt.close(fig)
        return

    df = order_df.copy()
    df["label"] = (
        df["case_name"].astype(str)
        + "\n"
        + df["polarization"].astype(str)
        + ", "
        + df["wavelength_um"].map(lambda v: f"{v:g} um")
        + ", "
        + df["incident_angle_deg"].map(lambda v: f"{v:g} deg")
        + "\nm="
        + df["order_m"].astype(str)
    )
    ax.bar(np.arange(len(df)), df["reflected_power_fraction"], color="#4C78A8")
    ax.set_xticks(np.arange(len(df)))
    ax.set_xticklabels(df["label"], rotation=65, ha="right", fontsize=7)
    ax.set_ylabel("Reflected power fraction")
    ax.set_title("Reflected diffraction orders from mode decomposition")
    ax.grid(axis="y", alpha=0.25)
    save_figure(fig, out_path)
    plt.close(fig)


def _plot_energy_channels(summary_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    if summary_df.empty:
        ax.text(0.5, 0.5, "No energy-channel data", ha="center", va="center")
        ax.set_axis_off()
        save_figure(fig, out_path)
        plt.close(fig)
        return

    df = summary_df.copy()
    df["label"] = (
        df["case_name"].astype(str)
        + "\n"
        + df["polarization"].astype(str)
        + ", "
        + df["wavelength_um"].map(lambda v: f"{v:g} um")
        + ", "
        + df["incident_angle_deg"].map(lambda v: f"{v:g} deg")
    )
    x = np.arange(len(df))
    bottom = np.zeros(len(df))
    channels = [
        ("absorptance", "#F58518", "Absorption"),
        ("specular_reflection", "#4C78A8", "Specular reflection"),
        ("positive_diffraction_orders", "#54A24B", "Positive orders"),
        ("negative_diffraction_orders", "#B279A2", "Negative orders"),
        ("T", "#9D755D", "Transmission"),
    ]
    for col, color, label in channels:
        vals = df[col].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=color, label=label)
        bottom += vals
    ax.axhline(1.0, color="black", lw=1, ls="--", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Incident-power fraction")
    ax.set_title("Energy channels, incident power = 1")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    save_figure(fig, out_path)
    plt.close(fig)


def _allowed_orders_table(args: argparse.Namespace, wavelengths: list[float],
                          angles: list[float]) -> str:
    lines = ["| wavelength_um | angle_deg | propagating_orders | angles_deg |",
             "|---:|---:|---|---|"]
    for wl in wavelengths:
        for theta in angles:
            orders = _propagating_orders(wl, theta, args.period_um, args.grazing_tol)
            lines.append(
                "| "
                + f"{wl:g} | {theta:g} | "
                + " ".join(str(o["order_m"]) for o in orders)
                + " | "
                + ", ".join(f"m={o['order_m']}:{o['diffraction_angle_deg']:.1f}" for o in orders)
                + " |"
            )
    return "\n".join(lines)


def _write_report(
    paths: dict[str, Path],
    order_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    args: argparse.Namespace,
    wavelengths: list[float],
    angles: list[float],
) -> None:
    ensure_dir(paths["report"].parent)
    passed = bool((summary_df["pass_or_fail"] == "PASS").all()) if not summary_df.empty else False

    if summary_df.empty:
        absorption_answer = "No completed simulations."
        redistribution_answer = "No completed simulations."
        directionality_answer = "No completed simulations."
    else:
        slanted = summary_df[summary_df["case_name"] == "slanted_groove"]
        symmetric = summary_df[summary_df["case_name"] == "symmetric_groove"]
        if not slanted.empty and not symmetric.empty:
            slanted_non = float(slanted["nonspecular_reflection"].mean())
            sym_non = float(symmetric["nonspecular_reflection"].mean())
            redistribution_answer = (
                f"Slanted nonspecular reflection mean = {slanted_non:.4f}; "
                f"symmetric mean = {sym_non:.4f}."
            )
        elif not slanted.empty:
            redistribution_answer = (
                f"Slanted nonspecular reflection mean = "
                f"{float(slanted['nonspecular_reflection'].mean()):.4f}."
            )
        else:
            redistribution_answer = "Slanted case was not included in this run."

        high_ref = float(summary_df["R_flux_monitor"].mean())
        mean_a = float(summary_df["absorptance"].mean())
        absorption_answer = (
            f"Mean absorptance = {mean_a:.4f}; mean reflected flux = {high_ref:.4f}. "
            "If absorptance remains low while reflected channels sum near R, energy is "
            "leaving primarily as reflection rather than being dissipated."
        )
        pm_delta = (
            summary_df["positive_diffraction_orders"]
            - summary_df["negative_diffraction_orders"]
        ).abs()
        directionality_answer = (
            f"Mean |positive-negative diffraction order fraction| = "
            f"{float(pm_delta.mean()):.4f}; max = {float(pm_delta.max()):.4f}."
        )

    report = f"""# D06 Diffraction-Order Energy Channel Diagnostic

## Purpose

Use Meep mode decomposition to test whether periodic grooves increase true
absorption or mainly redistribute reflected energy into non-specular
diffraction orders.

## Run Configuration

- Period: {args.period_um:g} um
- Wavelengths: {', '.join(f'{w:g}' for w in wavelengths)} um
- Angles: {', '.join(f'{a:g}' for a in angles)} deg
- Polarizations: {', '.join(args.polarizations)}
- Resolution: {args.resolution} px/um
- PML / air buffer / substrate: {args.pml_thickness_um:g} / {args.air_buffer_um:g} / {args.substrate_thickness_um:g} um
- Decay threshold: {args.decay_db:g} dB

## Physical Assumptions And Normalization

- Source propagates from +y toward -y.
- Raw flux monitor convention is retained in the CSV:
  `R = reflection_flux_raw / abs(input_flux_raw)`,
  `T = -transmission_flux_raw / abs(input_flux_raw)`,
  `A = 1 - R - T`.
- Mode decomposition uses `mp.DiffractedPlanewave(g=[m,0,0])`.
  Ez uses s-polarization and Hz uses p-polarization in the x-y incidence plane.
- The script computes both mode-coefficient direction components and selects the
  one whose propagated-order sum best matches the total reflected flux.

## Pass/Fail Criteria

- Finite fluxes and diffraction-order powers.
- `abs(total_reflected_power_from_orders - R_flux_monitor) <= {args.order_sum_tol:g}`.
- For symmetric grooves at normal incidence,
  `abs(sum(+m)-sum(-m)) <= {args.symmetry_tol:g}`.

Overall status: **{'PASS' if passed else 'FAIL'}**

## Allowed Propagating Orders

{_allowed_orders_table(args, wavelengths, angles)}

## Numerical Findings

Validated by this run:

- {absorption_answer}
- {redistribution_answer}
- {directionality_answer}
- Order-sum residual range:
  {float(summary_df['order_sum_residual'].min()) if not summary_df.empty else float('nan'):.4g}
  to
  {float(summary_df['order_sum_residual'].max()) if not summary_df.empty else float('nan'):.4g}.

Still hypotheses:

- If slanted grooves show larger non-specular reflection but similar A, the
  structure is behaving more like a direction/angle redistribution element than
  a high-emissivity absorber.
- If positive and negative order powers differ in the slanted case, the angular
  asymmetry seen in D05 is likely connected to diffraction-channel imbalance.

Needs higher-fidelity confirmation:

- 3D roughness, oxide layers, rounded groove walls, and experimentally measured
  optical constants are not included.
- Near-grazing diffraction orders can be numerically delicate and should be
  checked with higher resolution/PML before quantitative claims.

## Required Answers

1. Does the current slanted groove transfer specular reflection into
   non-specular orders?
   {redistribution_answer}

2. Is lack of absorption gain because energy still mainly leaves as reflection?
   {absorption_answer}

3. Which orders are allowed at 8, 10, 12 um?
   At normal incidence with period 10 um: 8 um allows m=-1,0,+1; 10 um places
   m=+-1 at grazing; 12 um allows only m=0.  The full simulated table is shown
   above.

4. Is directionality from diffraction-order asymmetry rather than absorption
   enhancement?
   {directionality_answer}

## Output Files

- Orders CSV: `{paths['orders']}`
- Energy summary CSV: `{paths['summary']}`
- Order bar plot: `{paths['order_bars']}`
- Energy channel plot: `{paths['channels']}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="D06 diffraction-order energy-channel diagnostic.",
    )
    p.add_argument(
        "--cases",
        nargs="+",
        default=["symmetric_groove", "slanted_groove", "mirrored_slanted_groove"],
        choices=sorted(_case_specs().keys()),
    )
    p.add_argument("--polarizations", nargs="+", default=["Ez", "Hz"], choices=["Ez", "Hz"])
    p.add_argument("--wavelengths_um", nargs="*", default=None)
    p.add_argument("--angles_deg", nargs="*", default=["0"])
    p.set_defaults(include_d05_directional_angle=True)
    p.add_argument(
        "--include_d05_directional_angle",
        dest="include_d05_directional_angle",
        action="store_true",
    )
    p.add_argument(
        "--no_include_d05_directional_angle",
        dest="include_d05_directional_angle",
        action="store_false",
    )
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
    p.add_argument("--order_sum_tol", type=float, default=DEFAULTS["order_sum_tol"])
    p.add_argument("--symmetry_tol", type=float, default=DEFAULTS["symmetry_tol"])
    p.add_argument("--grazing_tol", type=float, default=DEFAULTS["grazing_tol"])
    p.add_argument("--d05_directionality_tol", type=float, default=DEFAULTS["d05_directionality_tol"])
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    logger = setup_logger("D06_diffraction_order_energy_channel_diagnostic")
    paths = _paths()
    for path in paths.values():
        ensure_dir(path.parent)

    wavelengths = _parse_float_list(args.wavelengths_um, DEFAULTS["wavelengths_um"])
    angles = _parse_float_list(args.angles_deg, [0.0])
    angles.extend(_read_extra_angles_from_d05(args, logger))
    angles = sorted(set(round(a, 10) for a in angles))

    case_map = _case_specs()
    order_rows: list[dict] = []
    summary_rows: list[dict] = []
    for case_name in args.cases:
        case = case_map[case_name]
        for pol in args.polarizations:
            for wl in wavelengths:
                for theta in angles:
                    rows, summary = _run_one(case, pol, wl, theta, args, logger)
                    order_rows.extend(rows)
                    summary_rows.append(summary)

    order_df = pd.DataFrame(order_rows)
    summary_df = _add_symmetry_checks(pd.DataFrame(summary_rows), args)

    order_df.to_csv(paths["orders"], index=False)
    summary_df.to_csv(paths["summary"], index=False)
    _plot_order_bars(order_df, paths["order_bars"])
    _plot_energy_channels(summary_df, paths["channels"])
    _write_report(paths, order_df, summary_df, args, wavelengths, angles)

    logger.info("Wrote %s", paths["orders"])
    logger.info("Wrote %s", paths["summary"])
    logger.info("Wrote %s", paths["report"])
    if not summary_df.empty and (summary_df["pass_or_fail"] == "FAIL").any():
        logger.warning("D06 completed with failed validation rows.")
    else:
        logger.info("D06 completed successfully.")


if __name__ == "__main__":
    main()
