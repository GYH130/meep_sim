"""Validate diagnostics_v2 common numerical metrics.

This script is a smoke/acceptance check for the shared postprocessing and
single-wavelength solver utilities.  It writes only to diagnostics_v2 outputs
and does not overwrite D00-D08 results.
"""

from __future__ import annotations

import argparse
import cmath
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path("/private/tmp") / "matplotlib-codex-cache"),
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import meep as mp
import numpy as np
import pandas as pd

from src.io_utils import ensure_dir, project_path
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium
from src.postprocess import (
    assess_quantitative_validity,
    broadband_source_quality_mask,
    compute_RTA_downward_incidence,
    material_validity_mask,
    opaque_substrate_transmission_check,
    wavelength_integrated_average,
)
from src.simulation import run_periodic_2d_metal_single_wavelength


def _paths() -> dict[str, Path]:
    return {
        "table": project_path(
            "results", "diagnostics_v2", "tables",
            "V00_common_metrics_validation.csv",
        ),
        "report": project_path(
            "results", "diagnostics_v2", "reports",
            "V00_common_metrics_validation_report.md",
        ),
        "log": project_path(
            "logs", "diagnostics_v2", "V00_validate_common_metrics.log",
        ),
    }


def _setup_logger(path: Path) -> logging.Logger:
    ensure_dir(path.parent)
    logger = logging.getLogger("V00_validate_common_metrics")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler(path, encoding="utf-8")
    fh.setFormatter(formatter)
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info("Log file: %s", path)
    return logger


def _status(ok: bool) -> str:
    return "NUMERICAL_PASS" if ok else "FAIL"


def _flat_ti_single(args: argparse.Namespace, logger) -> dict:
    ti = get_ti_medium(args.wavelength_um, args.wavelength_um)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        return [
            mp.Block(
                material=ti,
                center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                size=mp.Vector3(args.period_um, substrate_thickness_um, mp.inf),
            )
        ]

    result = run_periodic_2d_metal_single_wavelength(
        geometry_factory=factory,
        period_um=args.period_um,
        wavelength_um=args.wavelength_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        decay_db=args.decay_db,
        source_component=args.polarization,
        fwidth_fraction=args.fwidth_fraction,
        logger=logger,
    )
    eps = ti.epsilon(1.0 / args.wavelength_um)[0][0]
    n_ti = cmath.sqrt(eps)
    r_fresnel = abs((1.0 - n_ti) / (1.0 + n_ti)) ** 2
    result["fresnel_A"] = 1.0 - r_fresnel
    result["abs_A_minus_fresnel"] = abs(result["A"] - result["fresnel_A"])
    return result


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Validate diagnostics_v2 common metrics.")
    p.add_argument("--wavelength_um", type=float, default=10.0)
    p.add_argument("--period_um", type=float, default=10.0)
    p.add_argument("--resolution", type=int, default=32)
    p.add_argument("--pml_thickness_um", type=float, default=2.0)
    p.add_argument("--substrate_thickness_um", type=float, default=4.0)
    p.add_argument("--air_buffer_um", type=float, default=4.0)
    p.add_argument("--decay_db", type=float, default=20.0)
    p.add_argument("--polarization", choices=["Ez", "Hz"], default="Ez")
    p.add_argument("--fwidth_fraction", type=float, default=0.06)
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    paths = _paths()
    for key in ("table", "report"):
        ensure_dir(paths[key].parent)
    logger = _setup_logger(paths["log"])

    rows = []
    wavelengths = np.array([13.0, 8.0, 9.0])
    values = np.array([1.0, 0.0, 0.0])
    integrated = wavelength_integrated_average(values, wavelengths, 8.0, 13.0)
    simple = float(np.mean(values))
    rows.append({
        "check_name": "wavelength_integrated_average",
        "status": _status(not np.isclose(integrated, simple) and np.isclose(integrated, 0.4)),
        "value": integrated,
        "reference_value": simple,
        "message": "trapezoid wavelength average differs from simple mean",
    })

    rta = compute_RTA_downward_incidence(
        np.array([8.0]), np.array([0.42]), np.array([-10.0]),
    )
    rows.append({
        "check_name": "signed_transmission_definition",
        "status": _status(np.isclose(rta["T"][0], -0.042)),
        "value": float(rta["T"][0]),
        "reference_value": -0.042,
        "message": "downward incidence T=-transmission_flux_raw/abs(input_flux_raw)",
    })

    source = broadband_source_quality_mask(np.array([1.0, 1e-4, 0.5]))
    low_ok = source["source_quality_flag"][1] == "LOW_SOURCE_SNR"
    rows.append({
        "check_name": "low_source_snr_flag",
        "status": _status(bool(low_ok)),
        "value": float(source["source_relative_flux"][1]),
        "reference_value": 1e-3,
        "message": "weak broadband endpoint is excluded from quantitative metrics",
    })

    trans = opaque_substrate_transmission_check(np.array([-0.042]))
    rows.append({
        "check_name": "negative_transmission_failure",
        "status": _status(trans["transmission_quality_flag"][0] == "FAIL"),
        "value": -0.042,
        "reference_value": 5e-3,
        "message": "T=-0.042 is a FAIL for opaque Ti substrate checks",
    })

    sim = _flat_ti_single(args, logger)
    mat_mask = material_validity_mask(
        np.array([sim["wavelength_um"]]),
        TI_RAKIC_VALID_LAMBDA_UM[1],
    )
    quality = assess_quantitative_validity(
        np.array([sim["wavelength_um"]]),
        np.array([sim["T"]]),
        np.array([sim["input_flux_raw"]]),
        mat_mask,
        sim["source_mode"],
    )
    flat_ok = sim["abs_A_minus_fresnel"] < 0.01 and abs(sim["T"]) < 1e-3
    rows.append({
        "check_name": "single_wavelength_flat_ti_10um",
        "status": _status(flat_ok),
        "value": sim["A"],
        "reference_value": sim["fresnel_A"],
        "message": "single-wavelength flat Ti A matches Fresnel and T is opaque",
        "solver_version": sim["solver_version"],
        "source_mode": sim["source_mode"],
        "polarization": sim["polarization"],
        "wavelength_um": sim["wavelength_um"],
        "resolution": sim["resolution"],
        "pml_thickness_um": sim["pml_thickness_um"],
        "substrate_thickness_um": sim["substrate_thickness_um"],
        "air_buffer_um": sim["air_buffer_um"],
        "decay_db": sim["decay_db"],
        "material_model": "Ti_Rakic_meep_builtin",
        "material_validity_flag": quality["material_validity_flag"][0],
        "numerical_quality_flag": quality["quality_flag"][0],
        "input_flux_raw": sim["input_flux_raw"],
        "reflection_flux_raw": sim["reflection_flux_raw"],
        "transmission_flux_raw": sim["transmission_flux_raw"],
        "reflectance": sim["R"],
        "transmittance": sim["T"],
        "absorptance": sim["A"],
        "abs_A_minus_fresnel": sim["abs_A_minus_fresnel"],
    })

    df = pd.DataFrame(rows)
    default_fields = {
        "solver_version": "diagnostics_v2_common_metrics",
        "source_mode": "synthetic",
        "polarization": args.polarization,
        "wavelength_um": np.nan,
        "resolution": args.resolution,
        "pml_thickness_um": args.pml_thickness_um,
        "substrate_thickness_um": args.substrate_thickness_um,
        "air_buffer_um": args.air_buffer_um,
        "decay_db": args.decay_db,
        "material_model": "synthetic",
        "material_validity_flag": "VALID",
        "numerical_quality_flag": "VALID",
    }
    for col, value in default_fields.items():
        if col not in df:
            df[col] = value
        else:
            df[col] = df[col].fillna(value)
    df.to_csv(paths["table"], index=False)

    overall = "NUMERICAL_PASS" if (df["status"] == "NUMERICAL_PASS").all() else "FAIL"
    lines = ["| check_name | status | value | reference_value | message |",
             "|---|---|---:|---:|---|"]
    for _, row in df.iterrows():
        lines.append(
            f"| {row['check_name']} | {row['status']} | "
            f"{row['value']:.6g} | {row['reference_value']:.6g} | "
            f"{row['message']} |"
        )

    report = f"""# V00 Common Metrics Validation

Overall status: **{overall}**

## Purpose

Validate diagnostics_v2 postprocessing and single-wavelength solver utilities
without writing to the original diagnostics directory.

## Acceptance Checks

{chr(10).join(lines)}

## Flat Ti 10 um Single-Wavelength Result

- A_meep: {sim['A']:.6g}
- A_fresnel: {sim['fresnel_A']:.6g}
- abs(A_meep - A_fresnel): {sim['abs_A_minus_fresnel']:.6g}
- T: {sim['T']:.6g}
- quality_flag: {quality['quality_flag'][0]}

## Output Files

- Table: `{paths['table']}`
- Report: `{paths['report']}`
- Log: `{paths['log']}`
"""
    paths["report"].write_text(report, encoding="utf-8")
    logger.info("Wrote %s", paths["table"])
    logger.info("Wrote %s", paths["report"])
    if overall != "NUMERICAL_PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
