"""Audit and evaluate the frozen Ordal-derived Ti model without running Meep.

No fitting is performed here.  Frequencies, oscillator frequencies and damping
rates use Meep's cyclic-frequency convention with one length unit = one µm.
PCHIP targets between source nodes are a specified optical model, not new
measurements.  The source does not validate a particular fabricated sample or
its high-temperature behavior.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def _parameters(model):
    if model.get("model") != "susceptibility":
        raise ValueError("MATERIAL_MODEL_INVALID: susceptibility model required")
    if float(model.get("meep_length_unit_um", math.nan)) != 1.0:
        raise ValueError("MATERIAL_UNIT_MISMATCH: one Meep length unit must equal 1 µm")
    base = float(model["epsilon"])
    if not math.isfinite(base) or base <= 0:
        raise ValueError("MATERIAL_INVALID: epsilon_inf must be finite and positive")
    poles = model.get("susceptibilities", [])
    if not isinstance(poles, list) or not poles:
        raise ValueError("MATERIAL_INVALID: missing susceptibility poles")
    checked = []
    for pole in poles:
        kind = str(pole["type"]).lower()
        f0, gamma, sigma = (float(pole[key]) for key in ("frequency", "gamma", "sigma"))
        if kind not in {"lorentz", "drude"}:
            raise ValueError("MATERIAL_INVALID: unsupported susceptibility")
        if not all(math.isfinite(v) for v in (f0, gamma, sigma)) or f0 <= 0 or gamma < 0 or sigma < 0:
            raise ValueError("MATERIAL_INVALID: invalid passive oscillator parameters")
        checked.append((kind, f0, gamma, sigma))
    band = np.asarray(model.get("wavelength_range_um", []), dtype=float)
    if band.shape != (2,) or not np.isfinite(band).all() or band[0] <= 0 or band[1] <= band[0]:
        raise ValueError("MATERIAL_INVALID: wavelength_range_um must be a positive ascending pair")
    return base, checked, band


def load_model(path):
    """Load a frozen JSON model and check its intrinsic physical parameters."""
    model = json.loads(Path(path).read_text())
    _parameters(model)
    return model


def epsilon(model, f):
    """Return complex relative epsilon at positive cyclic frequency f (1/µm)."""
    base, poles, _ = _parameters(model)
    frequency = np.asarray(f, dtype=float)
    if not np.isfinite(frequency).all() or np.any(frequency <= 0):
        raise ValueError("MATERIAL_INVALID_FREQUENCY")
    result = np.full(frequency.shape, base, dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for kind, f0, gamma, sigma in poles:
            denominator = -frequency * frequency - 1j * gamma * frequency
            if kind == "lorentz":
                denominator = denominator + f0 * f0
            result += sigma * f0 * f0 / denominator
    if not np.isfinite(result).all():
        raise ValueError("MATERIAL_NONFINITE_EPSILON")
    return result


def audit_model(model, resolution, courant, wavelength):
    """Fail closed before solver/cache admission; q<1 is only a rough screen."""
    _, poles, band = _parameters(model)
    if model.get("fit_accepted") is not True or model.get("route_a_ti") is not True:
        raise ValueError("MATERIAL_FIT_NOT_ACCEPTED")
    resolution, courant = float(resolution), float(courant)
    if not math.isfinite(resolution) or resolution <= 0 or not math.isfinite(courant) or not 0 < courant <= .5:
        raise ValueError("MATERIAL_INVALID_TIMESTEP: resolution>0 and 0<Courant<=0.5 required")
    waves = np.asarray(wavelength, dtype=float)
    if waves.size == 0 or not np.isfinite(waves).all() or np.any(waves <= 0):
        raise ValueError("MATERIAL_INVALID_WAVELENGTH")
    if np.any(waves < band[0]) or np.any(waves > band[1]):
        raise ValueError("MATERIAL_OUT_OF_BAND")
    highest = max((f0 for kind, f0, _, _ in poles if kind == "lorentz"), default=0.0)
    dt = courant / resolution
    q = math.pi * highest * dt
    if not math.isfinite(q) or q >= 1:
        raise ValueError(f"MATERIAL_UNSTABLE: q={q!r} must be <1")
    values = epsilon(model, 1 / waves)
    if np.any(values.imag < -1e-10):
        raise ValueError("MATERIAL_NONPASSIVE")
    return {"valid": True, "status": "ADMITTED_ROUGH_SCREEN_ONLY", "q": q,
            "pi_fmax_dt": q, "resolution": resolution, "courant": courant,
            "dt": dt, "highest_lorentz_frequency": highest,
            "rough_stability_passed": True, "numerically_qualified": False,
            "meep_length_unit_um": 1.0, "wavelength_range_um": band.tolist(),
            "checked_wavelength_um": waves.tolist(),
            "passivity_min_imag_epsilon": float(np.min(values.imag)),
            "criterion": "pi * highest_Lorentz_frequency * Courant / resolution < 1",
            "note": "Necessary rough dispersion screen, not time-domain convergence certification"}


def validate_fit(model, raw_csv, fit_report):
    """Recompute the fixed-scale 301-point PCHIP and original-node fit checks.

    Reported metrics are cross-checked against recomputation, and both stored
    source hashes must match the bytes of raw_csv.  Returns an invalid report
    instead of silently accepting stale fit metadata.  Never changes the model.
    """
    result = {"valid": False, "accepted": False, "errors": [],
              "interpretation": "specified Ti optical model; interpolated values are not measurements",
              "physical_sample_validated": False}
    errors = result["errors"]
    try:
        _, _, band = _parameters(model)
        report = (json.loads(Path(fit_report).read_text())
                  if isinstance(fit_report, (str, Path)) else fit_report)
        raw_path = Path(raw_csv)
        digest = hashlib.sha256(raw_path.read_bytes()).hexdigest()
        result["source_sha256"] = digest
        if model.get("source_sha256") != digest or report.get("source_sha256") != digest:
            errors.append("MATERIAL_SOURCE_HASH_MISMATCH")
        if model.get("fit_accepted") is not True or report.get("accepted") is not True:
            errors.append("MATERIAL_FIT_NOT_ACCEPTED")
        if model.get("route_a_ti") is not True:
            errors.append("MATERIAL_ROUTE_A_REQUIRED")
        data = np.genfromtxt(raw_path, delimiter=",", names=True, dtype=float)
        if data.dtype.names != ("wavelength_um", "n", "k"):
            raise ValueError("MATERIAL_SOURCE_COLUMNS_INVALID")
        wavelengths, n, k = (np.atleast_1d(data[name]) for name in data.dtype.names)
        if len(wavelengths) < 2 or not np.isfinite(np.column_stack([wavelengths, n, k])).all():
            raise ValueError("MATERIAL_SOURCE_NONFINITE")
        if np.any(wavelengths <= 0) or np.any(n < 0) or np.any(k < 0):
            raise ValueError("MATERIAL_SOURCE_NONPASSIVE")
        order = np.argsort(wavelengths)
        wavelengths, n, k = wavelengths[order], n[order], k[order]
        if np.any(np.diff(wavelengths) <= 0):
            raise ValueError("MATERIAL_SOURCE_DUPLICATE_WAVELENGTH")
        if band[0] < wavelengths[0] or band[1] > wavelengths[-1]:
            raise ValueError("MATERIAL_SOURCE_COVERAGE_MISSING: extrapolation prohibited")
        grid = np.linspace(band[0], band[1], 301)
        ni = PchipInterpolator(wavelengths, n, extrapolate=False)
        ki = PchipInterpolator(wavelengths, k, extrapolate=False)
        target = (ni(grid) + 1j * ki(grid)) ** 2
        scale = max(float(np.sqrt(np.mean(np.abs(target) ** 2))), 1.0)
        grid_error = np.abs(epsilon(model, 1 / grid) - target) / scale
        original = (wavelengths >= band[0]) & (wavelengths <= band[1])
        if not np.any(original):
            raise ValueError("MATERIAL_ORIGINAL_NODES_MISSING")
        node_target = (n[original] + 1j * k[original]) ** 2
        node_error = np.abs(epsilon(model, 1 / wavelengths[original]) - node_target) / scale
        passivity_grid = np.linspace(1 / band[1], 1 / band[0], 3010)
        metrics = {
            "normalization_scale": scale,
            "normalized_rmse": float(np.sqrt(np.mean(grid_error ** 2))),
            "normalized_max_error": float(grid_error.max()),
            "original_node_rmse": float(np.sqrt(np.mean(node_error ** 2))),
            "original_node_max_error": float(node_error.max()),
            "passivity_min_imag_epsilon": float(epsilon(model, passivity_grid).imag.min()),
            "original_nodes": int(original.sum()), "interpolated_evaluation_points": 301,
        }
        result.update(metrics)
        result.update(band_um=band.tolist(), meep_length_unit_um=1.0,
                      normalization="max(RMS(abs(epsilon_data)),1), fixed over band",
                      interpolation="PCHIP n,k(lambda), no extrapolation")
        for name, value in metrics.items():
            previous = report.get(name)
            if previous is None or not np.isfinite(float(previous)) or not np.isclose(value, float(previous), rtol=1e-8, atol=1e-12):
                errors.append(f"MATERIAL_FIT_REPORT_MISMATCH:{name}")
        if not np.array_equal(np.asarray(report.get("band_um", [])), band):
            errors.append("MATERIAL_FIT_REPORT_MISMATCH:band_um")
        for provenance_key in ("doi", "source_url", "data_file"):
            source = model.get("provenance", {}).get(provenance_key)
            if not source or source != report.get("source", {}).get(provenance_key):
                errors.append(f"MATERIAL_PROVENANCE_MISMATCH:{provenance_key}")
        if metrics["normalized_rmse"] > .02 or metrics["original_node_rmse"] > .02:
            errors.append("MATERIAL_FIT_RMSE_EXCEEDED")
        if metrics["normalized_max_error"] > .05 or metrics["original_node_max_error"] > .05:
            errors.append("MATERIAL_FIT_MAX_ERROR_EXCEEDED")
        if metrics["passivity_min_imag_epsilon"] < -1e-10:
            errors.append("MATERIAL_NONPASSIVE")
        result["valid"] = result["accepted"] = not errors
    except (KeyError, TypeError, ValueError, OSError, OverflowError) as exc:
        errors.append(str(exc))
    result["status"] = "ACCEPTED" if result["valid"] else "MATERIAL_FIT_INVALID"
    return result


def build_medium(model):
    """Create Meep's matching dispersive medium; Meep is imported only here."""
    base, poles, band = _parameters(model)
    import meep as mp

    susceptibilities = []
    for kind, frequency, gamma, sigma in poles:
        cls = mp.LorentzianSusceptibility if kind == "lorentz" else mp.DrudeSusceptibility
        susceptibilities.append(cls(frequency=frequency, gamma=gamma, sigma=sigma))
    return mp.Medium(epsilon=base, E_susceptibilities=susceptibilities,
                     valid_freq_range=mp.FreqRange(min=1 / band[1], max=1 / band[0]))
