"""D03_numerical_convergence_diagnostic.py — 数值收敛诊断

研究目的
--------
对平面 Ti 与代表性 20° 斜槽结构进行单因素数值收敛测试，判断当前 8-13 um
发射率/吸收率结论是否对 resolution、PML、空气缓冲、基底厚度和衰减阈值稳定。

物理假设
--------
1. 所有长度单位为 um，Meep 频率 f = 1 / wavelength_um；
2. 正入射、2D 周期结构；不透明 Ti 基底下 emissivity_proxy = absorptance；
3. Ti 使用 Meep 内置 Rakić Drude-Lorentz 模型，超过 12.398 um 的点属于外推；
4. 当前诊断只判断数值稳定性，不证明真实飞秒加工 3D 粗糙/氧化结构的发射率。

通过/失败判据
-------------
1. 每个仿真结果 R/T/A 和 raw flux 必须有限；
2. 每个输出保留 input_flux_raw、reflection_flux_raw、transmission_flux_raw；
3. 单因素扫描中，与该扫描最高精度 case 的 8-13 um 最大 A 差异小于
   convergence_tol 时视为收敛；
4. 平面 Ti 同时与半无限 Fresnel 解比较。
"""

from __future__ import annotations

import argparse
import cmath
import json
import sys
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
from src.simulation import run_periodic_2d_metal_spectrum


DEFAULTS = dict(
    wavelength_min_um=8.0,
    wavelength_max_um=13.0,
    nfreq=21,
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    baseline_resolution=48,
    baseline_pml_thickness_um=4.0,
    baseline_air_buffer_um=8.0,
    baseline_substrate_thickness_um=4.0,
    baseline_decay_db=60.0,
    resolution_values=[32, 48, 64],
    pml_values=[2.0, 4.0, 6.0],
    air_buffer_values=[8.0, 12.0],
    substrate_values=[2.0, 4.0, 6.0],
    decay_values=[40.0, 60.0],
    convergence_tol=1e-2,
    high_emissivity_threshold=0.20,
)


@dataclass(frozen=True)
class NumericConfig:
    resolution: int
    pml_thickness_um: float
    air_buffer_um: float
    substrate_thickness_um: float
    decay_db: float

    def label(self) -> str:
        return (
            f"res={self.resolution},pml={self.pml_thickness_um:g},"
            f"air={self.air_buffer_um:g},sub={self.substrate_thickness_um:g},"
            f"decay={self.decay_db:g}"
        )


def _paths() -> dict[str, Path]:
    return {
        "resolution_csv": project_path(
            "results", "diagnostics", "tables", "D03_resolution_convergence.csv",
        ),
        "pml_csv": project_path(
            "results", "diagnostics", "tables", "D03_pml_convergence.csv",
        ),
        "other_csv": project_path(
            "results", "diagnostics", "tables",
            "D03_buffer_substrate_decay_convergence.csv",
        ),
        "resolution_png": project_path(
            "results", "diagnostics", "figures", "D03_resolution_spectra.png",
        ),
        "summary_png": project_path(
            "results", "diagnostics", "figures",
            "D03_convergence_metric_summary.png",
        ),
        "report": project_path(
            "results", "diagnostics", "reports",
            "D03_numerical_convergence_report.md",
        ),
    }


def _case_specs() -> list[dict]:
    return [
        dict(case_name="flat_Ti", kind="flat", polarization="Ez",
             label="Flat Ti"),
        dict(case_name="slanted_Ez", kind="slanted", polarization="Ez",
             label="Slanted tilt=20 Ez"),
        dict(case_name="slanted_Hz", kind="slanted", polarization="Hz",
             label="Slanted tilt=20 Hz"),
    ]


def _baseline(args: argparse.Namespace) -> NumericConfig:
    return NumericConfig(
        resolution=args.baseline_resolution,
        pml_thickness_um=args.baseline_pml_thickness_um,
        air_buffer_um=args.baseline_air_buffer_um,
        substrate_thickness_um=args.baseline_substrate_thickness_um,
        decay_db=args.baseline_decay_db,
    )


def _sweep_items(args: argparse.Namespace) -> list[dict]:
    base = _baseline(args)
    items = []
    for val in args.resolution_values:
        items.append(dict(
            scan_name="resolution",
            varied_parameter="resolution",
            parameter_value=float(val),
            config=NumericConfig(int(val), base.pml_thickness_um, base.air_buffer_um,
                                 base.substrate_thickness_um, base.decay_db),
            highest_value=float(max(args.resolution_values)),
        ))
    for val in args.pml_values:
        items.append(dict(
            scan_name="pml",
            varied_parameter="pml_thickness_um",
            parameter_value=float(val),
            config=NumericConfig(base.resolution, float(val), base.air_buffer_um,
                                 base.substrate_thickness_um, base.decay_db),
            highest_value=float(max(args.pml_values)),
        ))
    for val in args.air_buffer_values:
        items.append(dict(
            scan_name="air_buffer",
            varied_parameter="air_buffer_um",
            parameter_value=float(val),
            config=NumericConfig(base.resolution, base.pml_thickness_um, float(val),
                                 base.substrate_thickness_um, base.decay_db),
            highest_value=float(max(args.air_buffer_values)),
        ))
    for val in args.substrate_values:
        items.append(dict(
            scan_name="substrate",
            varied_parameter="substrate_thickness_um",
            parameter_value=float(val),
            config=NumericConfig(base.resolution, base.pml_thickness_um,
                                 base.air_buffer_um, float(val), base.decay_db),
            highest_value=float(max(args.substrate_values)),
        ))
    for val in args.decay_values:
        items.append(dict(
            scan_name="decay",
            varied_parameter="decay_db",
            parameter_value=float(val),
            config=NumericConfig(base.resolution, base.pml_thickness_um,
                                 base.air_buffer_um, base.substrate_thickness_um,
                                 float(val)),
            highest_value=float(max(args.decay_values)),
        ))
    return items


def _ti_fresnel(wavelengths_um: np.ndarray, ti_medium) -> tuple[np.ndarray, np.ndarray]:
    lo_valid, hi_valid = TI_RAKIC_VALID_LAMBDA_UM
    r_values = []
    for wl in wavelengths_um:
        if wl < lo_valid or wl > hi_valid:
            r_values.append(np.nan)
            continue
        eps = ti_medium.epsilon(1.0 / float(wl))[0][0]
        n_ti = cmath.sqrt(eps)
        r_values.append(abs((1.0 - n_ti) / (1.0 + n_ti)) ** 2)
    r = np.asarray(r_values, dtype=float)
    return r, 1.0 - r


def _run_case(case: dict, config: NumericConfig, args: argparse.Namespace, logger) -> pd.DataFrame:
    ti = get_ti_medium(
        lambda_min_um=args.wavelength_min_um,
        lambda_max_um=args.wavelength_max_um,
    )
    if case["kind"] == "flat":
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return [
                mp.Block(
                    material=ti,
                    center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                    size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
                )
            ]
    else:
        def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
            return build_slanted_groove_geometry(
                period_x_um=args.period_um,
                top_width_um=args.top_width_um,
                bottom_width_um=args.bottom_width_um,
                depth_um=args.depth_um,
                tilt_angle_deg=20.0,
                substrate_thickness_um=substrate_thickness_um,
                y_surface=y_surface_um,
                medium_substrate=ti,
            )

    logger.info(">>> case=%s config=%s", case["case_name"], config.label())
    result = run_periodic_2d_metal_spectrum(
        geometry_factory=factory,
        period_um=args.period_um,
        wavelength_min_um=args.wavelength_min_um,
        wavelength_max_um=args.wavelength_max_um,
        resolution=config.resolution,
        pml_thickness_um=config.pml_thickness_um,
        substrate_thickness_um=config.substrate_thickness_um,
        air_buffer_um=config.air_buffer_um,
        nfreq=args.nfreq,
        decay_db=config.decay_db,
        source_component=case["polarization"],
        logger=logger,
    )

    required = ["input_flux_raw", "reflection_flux_raw", "transmission_flux_raw"]
    missing = [k for k in required if k not in result]
    if missing:
        raise RuntimeError(f"{case['case_name']} missing raw flux fields: {missing}")

    df = pd.DataFrame({
        "case_name": case["case_name"],
        "case_label": case["label"],
        "polarization": case["polarization"],
        "wavelength_um": result["wavelength_um"],
        "reflectance": result["R"],
        "transmittance": result["T"],
        "absorptance": result["A"],
        "emissivity_proxy": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        "resolution": config.resolution,
        "pml_thickness_um": config.pml_thickness_um,
        "air_buffer_um": config.air_buffer_um,
        "substrate_thickness_um": config.substrate_thickness_um,
        "decay_db": config.decay_db,
        "config_label": config.label(),
        "walltime_s": float(result["walltime_s"]),
        "normalization_note": (
            "R = reflection_flux_raw / abs(input_flux_raw); "
            "T = -transmission_flux_raw / abs(input_flux_raw); "
            "A = 1 - R - T"
        ),
    })
    return df


def _band(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["wavelength_um"] >= 8.0) & (df["wavelength_um"] <= 13.0)]


def _add_metrics_for_sweep(sweep_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    ti = get_ti_medium(lambda_min_um=args.wavelength_min_um,
                       lambda_max_um=args.wavelength_max_um)
    for (scan_name, case_name), group in sweep_df.groupby(["scan_name", "case_name"], sort=False):
        ref_value = group["highest_accuracy_parameter_value"].iloc[0]
        ref = group[group["parameter_value"] == ref_value]
        if ref.empty:
            raise RuntimeError(f"missing highest accuracy row for {scan_name}/{case_name}")
        ref = ref.sort_values("wavelength_um")
        ref_band = _band(ref)
        for param_value, df in group.groupby("parameter_value", sort=False):
            df = df.sort_values("wavelength_um")
            band = _band(df)
            if not np.allclose(df["wavelength_um"].to_numpy(),
                               ref["wavelength_um"].to_numpy()):
                raise RuntimeError(f"wavelength mismatch for {scan_name}/{case_name}")
            peak_idx = band["absorptance"].idxmax()
            diff = float(np.max(np.abs(
                band["absorptance"].to_numpy()
                - ref_band["absorptance"].to_numpy()
            )))
            flat_fresnel_diff = np.nan
            if case_name == "flat_Ti":
                _r_f, a_f = _ti_fresnel(df["wavelength_um"].to_numpy(), ti)
                mask = np.isfinite(a_f)
                flat_fresnel_diff = float(np.max(np.abs(
                    df["absorptance"].to_numpy()[mask] - a_f[mask]
                ))) if np.any(mask) else np.nan
            finite = bool(np.all(np.isfinite(df[[
                "reflectance", "transmittance", "absorptance",
                "input_flux_raw", "reflection_flux_raw", "transmission_flux_raw",
            ]].to_numpy())))
            rows.append(dict(
                scan_name=scan_name,
                varied_parameter=df["varied_parameter"].iloc[0],
                parameter_value=float(param_value),
                highest_accuracy_parameter_value=float(ref_value),
                case_name=case_name,
                case_label=df["case_label"].iloc[0],
                polarization=df["polarization"].iloc[0],
                resolution=int(df["resolution"].iloc[0]),
                pml_thickness_um=float(df["pml_thickness_um"].iloc[0]),
                air_buffer_um=float(df["air_buffer_um"].iloc[0]),
                substrate_thickness_um=float(df["substrate_thickness_um"].iloc[0]),
                decay_db=float(df["decay_db"].iloc[0]),
                mean_A_8_13um=float(band["absorptance"].mean()),
                peak_A_8_13um=float(band.loc[peak_idx, "absorptance"]),
                peak_wavelength_um=float(band.loc[peak_idx, "wavelength_um"]),
                max_abs_difference_from_highest_accuracy_case=diff,
                max_abs_difference_from_fresnel_flat_Ti=flat_fresnel_diff,
                max_abs_T=float(np.max(np.abs(df["transmittance"]))),
                min_A=float(df["absorptance"].min()),
                max_A=float(df["absorptance"].max()),
                raw_flux_retained=True,
                all_finite=finite,
                pass_or_fail="PASS" if finite and diff <= args.convergence_tol else "FAIL",
                normalization_note=df["normalization_note"].iloc[0],
            ))
    return pd.DataFrame(rows)


def _expand_metric_to_spectra(sweep_df: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    keys = ["scan_name", "case_name", "parameter_value"]
    keep = keys + [
        "mean_A_8_13um", "peak_A_8_13um", "peak_wavelength_um",
        "max_abs_difference_from_highest_accuracy_case",
        "max_abs_difference_from_fresnel_flat_Ti", "pass_or_fail",
    ]
    return sweep_df.merge(metrics[keep], on=keys, how="left")


def _make_sweep_tables(
    spectra_cache: dict[tuple[str, NumericConfig], pd.DataFrame],
    items: list[dict],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames = []
    for item in items:
        for case in _case_specs():
            df = spectra_cache[(case["case_name"], item["config"])].copy()
            df["scan_name"] = item["scan_name"]
            df["varied_parameter"] = item["varied_parameter"]
            df["parameter_value"] = item["parameter_value"]
            df["highest_accuracy_parameter_value"] = item["highest_value"]
            frames.append(df)
    all_sweeps = pd.concat(frames, ignore_index=True)

    metrics = _add_metrics_for_sweep(all_sweeps, args)
    enriched = _expand_metric_to_spectra(all_sweeps, metrics)
    resolution = enriched[enriched["scan_name"] == "resolution"].copy()
    pml = enriched[enriched["scan_name"] == "pml"].copy()
    other = enriched[enriched["scan_name"].isin(["air_buffer", "substrate", "decay"])].copy()
    return resolution, pml, other, metrics


def _plot_resolution(resolution_df: pd.DataFrame, out_path: Path) -> Path:
    cases = list(resolution_df["case_name"].drop_duplicates())
    fig, axes = plt.subplots(len(cases), 1, figsize=(7.8, 3.3 * len(cases)), sharex=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        subset = resolution_df[resolution_df["case_name"] == case]
        for val, df in subset.groupby("parameter_value", sort=True):
            ax.plot(df["wavelength_um"], df["absorptance"], lw=1.5,
                    label=f"res={int(val)}")
        ax.axvspan(8, 13, color="orange", alpha=0.10)
        ax.set_ylabel("A")
        ax.set_title(str(subset["case_label"].iloc[0]))
        ax.grid(True, ls=":", alpha=0.6)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel("Wavelength (um)")
    fig.suptitle("D03 resolution convergence spectra", y=0.995)
    return save_figure(fig, out_path)


def _plot_summary(metrics: pd.DataFrame, out_path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    summary = metrics.copy()
    summary["label"] = (
        summary["scan_name"] + "\n" + summary["case_name"] + "\n"
        + summary["parameter_value"].map(lambda x: f"{x:g}")
    )
    x = np.arange(len(summary))
    ax.bar(x, summary["max_abs_difference_from_highest_accuracy_case"])
    ax.set_xticks(x)
    ax.set_xticklabels(summary["label"], rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("max |delta A| vs highest accuracy, 8-13 um")
    ax.set_title("D03 convergence metric summary")
    ax.grid(True, axis="y", ls=":", alpha=0.6)
    return save_figure(fig, out_path)


def _lookup_diff(metrics: pd.DataFrame, scan: str, case: str, value: float) -> float | None:
    rows = metrics[
        (metrics["scan_name"] == scan)
        & (metrics["case_name"] == case)
        & np.isclose(metrics["parameter_value"], value)
    ]
    if rows.empty:
        return None
    return float(rows["max_abs_difference_from_highest_accuracy_case"].iloc[0])


def _baseline_diffs(metrics: pd.DataFrame, args: argparse.Namespace) -> dict[str, float | None]:
    return {
        "resolution_32": max(
            (_lookup_diff(metrics, "resolution", c["case_name"], 32.0) or 0.0)
            for c in _case_specs()
        ),
        "pml_2": max(
            (_lookup_diff(metrics, "pml", c["case_name"], 2.0) or 0.0)
            for c in _case_specs()
        ),
        "baseline_resolution_48": max(
            (_lookup_diff(metrics, "resolution", c["case_name"],
                          float(args.baseline_resolution)) or 0.0)
            for c in _case_specs()
        ),
        "baseline_pml_4": max(
            (_lookup_diff(metrics, "pml", c["case_name"],
                          float(args.baseline_pml_thickness_um)) or 0.0)
            for c in _case_specs()
        ),
        "baseline_air_8": max(
            (_lookup_diff(metrics, "air_buffer", c["case_name"],
                          float(args.baseline_air_buffer_um)) or 0.0)
            for c in _case_specs()
        ),
        "baseline_substrate_4": max(
            (_lookup_diff(metrics, "substrate", c["case_name"],
                          float(args.baseline_substrate_thickness_um)) or 0.0)
            for c in _case_specs()
        ),
        "baseline_decay_60": max(
            (_lookup_diff(metrics, "decay", c["case_name"],
                          float(args.baseline_decay_db)) or 0.0)
            for c in _case_specs()
        ),
    }


def _write_report(metrics: pd.DataFrame, args: argparse.Namespace,
                  paths: dict[str, Path]) -> Path:
    diffs = _baseline_diffs(metrics, args)
    res32_ok = diffs["resolution_32"] is not None and diffs["resolution_32"] <= args.convergence_tol
    pml2_ok = diffs["pml_2"] is not None and diffs["pml_2"] <= args.convergence_tol
    baseline_keys = [
        "baseline_resolution_48", "baseline_pml_4", "baseline_air_8",
        "baseline_substrate_4", "baseline_decay_60",
    ]
    baseline_ok = all((diffs[k] is not None and diffs[k] <= args.convergence_tol)
                      for k in baseline_keys)

    slanted = metrics[metrics["case_name"].isin(["slanted_Ez", "slanted_Hz"])]
    max_slanted_mean = float(slanted["mean_A_8_13um"].max())
    min_slanted_mean = float(slanted["mean_A_8_13um"].min())
    low_stable = max_slanted_mean < args.high_emissivity_threshold

    recommended = (
        f"baseline resolution={args.baseline_resolution}, pml={args.baseline_pml_thickness_um:g}, "
        f"air_buffer={args.baseline_air_buffer_um:g}, substrate={args.baseline_substrate_thickness_um:g}, "
        f"decay={args.baseline_decay_db:g}"
        if baseline_ok else
        "use the highest tested settings for parameters that exceed tolerance "
        "(resolution=max, pml=max, air_buffer=max, substrate=max, decay=max)"
    )

    lines = [
        "# D03 Numerical Convergence Diagnostic Report",
        "",
        "## Run Configuration",
        "",
        f"- wavelength range: {args.wavelength_min_um:g}-{args.wavelength_max_um:g} um",
        f"- nfreq: {args.nfreq}",
        f"- baseline: resolution={args.baseline_resolution}, pml={args.baseline_pml_thickness_um:g}, "
        f"air={args.baseline_air_buffer_um:g}, substrate={args.baseline_substrate_thickness_um:g}, "
        f"decay={args.baseline_decay_db:g}",
        f"- convergence tolerance: {args.convergence_tol:g}",
        "",
        "## Pass / Fail Criteria",
        "",
        "- all R/T/A and raw flux columns finite",
        "- each row preserves input_flux_raw, reflection_flux_raw, transmission_flux_raw",
        "- max |delta A| in 8-13 um versus the highest-accuracy value in that single-factor scan <= tolerance",
        "",
        "## Required Answers",
        "",
        f"1. Is current default resolution=32 enough for quantitative conclusions? "
        f"{'Yes' if res32_ok else 'No / not proven by this run'}; "
        f"max difference vs highest resolution = {diffs['resolution_32']}.",
        f"2. Does current PML=2 um affect 8-13 um results? "
        f"{'No significant effect within tolerance' if pml2_ok else 'Yes / not proven negligible'}; "
        f"max difference vs highest PML = {diffs['pml_2']}.",
        f"3. Recommended unified baseline for later diagnostics: {recommended}.",
        f"4. Is the low-emissivity conclusion stable to numerical parameters? "
        f"{'Yes' if low_stable else 'No'}; slanted mean A range across scans = "
        f"[{min_slanted_mean:.6f}, {max_slanted_mean:.6f}].",
        "",
        "## Verified Numerical Conclusions",
        "",
    ]
    for _, row in metrics.iterrows():
        lines.append(
            f"- {row['scan_name']} {row['case_name']} {row['varied_parameter']}={row['parameter_value']:g}: "
            f"mean A={row['mean_A_8_13um']:.6f}, "
            f"max |delta A|={row['max_abs_difference_from_highest_accuracy_case']:.3e}, "
            f"{row['pass_or_fail']}."
        )
    lines.extend([
        "",
        "## Hypotheses Still Not Proven",
        "",
        "- D03 does not prove which physical mode controls absorption; it only tests numerical sensitivity.",
        "- D03 does not include oxidation, rounded sidewalls, multiscale roughness, or 3D finite structures.",
        "",
        "## Needs Experiment Or Higher-Fidelity Model",
        "",
        "- Ti optical constants near the Rakić upper bound and laser-induced oxide layers need experimental confirmation.",
        "- If any scan fails tolerance in full runs, later diagnostics should use the stricter numerical settings.",
        "",
        "## Output Files",
        "",
        f"- resolution CSV: `{paths['resolution_csv']}`",
        f"- PML CSV: `{paths['pml_csv']}`",
        f"- buffer/substrate/decay CSV: `{paths['other_csv']}`",
        f"- resolution PNG: `{paths['resolution_png']}`",
        f"- summary PNG: `{paths['summary_png']}`",
        f"- report: `{paths['report']}`",
    ])
    ensure_dir(paths["report"].parent)
    paths["report"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths["report"]


def _float_list(values: list[float]) -> list[float]:
    return [float(v) for v in values]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--wavelength_min_um", type=float, default=DEFAULTS["wavelength_min_um"])
    p.add_argument("--wavelength_max_um", type=float, default=DEFAULTS["wavelength_max_um"])
    p.add_argument("--nfreq", type=int, default=DEFAULTS["nfreq"])
    p.add_argument("--period_um", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--top_width_um", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width_um", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth_um", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--baseline_resolution", type=int, default=DEFAULTS["baseline_resolution"])
    p.add_argument("--baseline_pml_thickness_um", type=float,
                   default=DEFAULTS["baseline_pml_thickness_um"])
    p.add_argument("--baseline_air_buffer_um", type=float,
                   default=DEFAULTS["baseline_air_buffer_um"])
    p.add_argument("--baseline_substrate_thickness_um", type=float,
                   default=DEFAULTS["baseline_substrate_thickness_um"])
    p.add_argument("--baseline_decay_db", type=float, default=DEFAULTS["baseline_decay_db"])
    p.add_argument("--resolution_values", type=int, nargs="+",
                   default=DEFAULTS["resolution_values"])
    p.add_argument("--pml_values", type=float, nargs="+", default=DEFAULTS["pml_values"])
    p.add_argument("--air_buffer_values", type=float, nargs="+",
                   default=DEFAULTS["air_buffer_values"])
    p.add_argument("--substrate_values", type=float, nargs="+",
                   default=DEFAULTS["substrate_values"])
    p.add_argument("--decay_values", type=float, nargs="+", default=DEFAULTS["decay_values"])
    p.add_argument("--convergence_tol", type=float, default=DEFAULTS["convergence_tol"])
    p.add_argument("--high_emissivity_threshold", type=float,
                   default=DEFAULTS["high_emissivity_threshold"])
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.wavelength_min_um <= 0 or args.wavelength_min_um >= args.wavelength_max_um:
        raise ValueError("wavelength_min_um must be positive and < wavelength_max_um")
    if args.depth_um >= min(args.substrate_values + [args.baseline_substrate_thickness_um]):
        raise ValueError("depth_um must be smaller than all tested substrate thickness values")
    if args.nfreq <= 1:
        raise ValueError("nfreq must be > 1")
    if min(args.resolution_values + [args.baseline_resolution]) < 32:
        raise ValueError("Ti Drude-Lorentz stability requires resolution >= 32")

    args.pml_values = _float_list(args.pml_values)
    args.air_buffer_values = _float_list(args.air_buffer_values)
    args.substrate_values = _float_list(args.substrate_values)
    args.decay_values = _float_list(args.decay_values)

    logger = setup_logger("D03_numerical_convergence_diagnostic")
    logger.info("=== D03_numerical_convergence_diagnostic ===")
    logger.info("args = %s", vars(args))
    mp.verbosity(1)

    paths = _paths()
    items = _sweep_items(args)
    needed = {}
    for item in items:
        for case in _case_specs():
            needed[(case["case_name"], item["config"])] = case

    spectra_cache: dict[tuple[str, NumericConfig], pd.DataFrame] = {}
    for (case_name, config), case in needed.items():
        spectra_cache[(case_name, config)] = _run_case(case, config, args, logger)

    resolution_df, pml_df, other_df, metrics = _make_sweep_tables(
        spectra_cache, items, args,
    )

    ensure_dir(paths["resolution_csv"].parent)
    resolution_df.to_csv(paths["resolution_csv"], index=False)
    pml_df.to_csv(paths["pml_csv"], index=False)
    other_df.to_csv(paths["other_csv"], index=False)
    _plot_resolution(resolution_df, paths["resolution_png"])
    _plot_summary(metrics, paths["summary_png"])
    _write_report(metrics, args, paths)

    ok = bool((metrics["pass_or_fail"] == "PASS").all())
    logger.info("resolution CSV → %s", paths["resolution_csv"])
    logger.info("PML CSV → %s", paths["pml_csv"])
    logger.info("other CSV → %s", paths["other_csv"])
    logger.info("resolution PNG → %s", paths["resolution_png"])
    logger.info("summary PNG → %s", paths["summary_png"])
    logger.info("report → %s", paths["report"])
    logger.info("metrics summary: %s", json.dumps(
        metrics.to_dict(orient="records"), ensure_ascii=False, default=str))
    logger.info("overall status: %s", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
