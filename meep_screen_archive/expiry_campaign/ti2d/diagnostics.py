"""Independent, Meep-free physical diagnostics for the two-dimensional solver.

Coordinates are Meep's x/y: physical slab-normal z is simulation y.  Frequencies
and wavevectors are cycles per micrometre, rather than radians per micrometre.
DFT power uses Meep's Re(E* x H) convention (no extra factor of one half).

References:
https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/#absorbed-power-density
https://meep.readthedocs.io/en/latest/Mode_Decomposition/
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np


def _real_scalar(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value) or np.iscomplexobj(value):
        raise ValueError(f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _array(value, name, *, real=False):
    array = np.asarray(value)
    if array.size == 0 or array.dtype.kind not in "biufc":
        raise ValueError(f"{name} must be a nonempty numeric array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    if real and np.iscomplexobj(array):
        raise ValueError(f"{name} must be real")
    return array


def _polarization(pol):
    if not isinstance(pol, str) or pol.lower() not in ("s", "p"):
        raise ValueError("pol must be 's' (Ez) or 'p' (Hz)")
    return pol.lower()


def fresnel_slab(eps, lambda_um, theta_deg, pol, thickness_um=60):
    """Finite air/material/air slab reflectance, transmittance and absorption.

    ``eps`` is the passive relative permittivity for exp(-i omega t).  Only a
    decaying propagation exponential is evaluated, including optically thick
    metal.  ``r`` and ``t`` are Ez amplitudes for s and Hz amplitudes for p.
    Incident and exit media are identical, hence R=|r|² and T=|t|².
    """
    wavelength = _real_scalar(lambda_um, "lambda_um", positive=True)
    thickness = _real_scalar(thickness_um, "thickness_um", nonnegative=True)
    theta = _real_scalar(theta_deg, "theta_deg")
    pol = _polarization(pol)
    if abs(theta) >= 90:
        raise ValueError("theta_deg must lie strictly between -90 and 90")
    if not np.isscalar(eps):
        raise ValueError("eps must be a finite complex scalar")
    eps = complex(eps)
    if not math.isfinite(eps.real) or not math.isfinite(eps.imag):
        raise ValueError("eps must be finite")
    if eps.imag < 0:
        raise ValueError("passive eps must have nonnegative imaginary part")
    angle = math.radians(theta)
    q0 = math.cos(angle)
    q1 = complex(np.sqrt(complex(eps - math.sin(angle) ** 2)))
    if q1.imag < 0 or (q1.imag == 0 and q1.real < 0):
        q1 = -q1
    phase_scale = 2 * math.pi * thickness / wavelength
    if thickness == 0:
        r, t = 0j, 1 + 0j
    elif abs(q1) < 1e-12:
        # Exact zero-normal-wavevector limit avoids a removable 0/0.
        if pol == "p" and abs(eps) < 1e-24:
            denominator = 2 * q0 - 1j * phase_scale
            r = 1j * phase_scale / denominator
            t = 2 * q0 / denominator
        else:
            admittance = q0 if pol == "s" else eps * q0
            denominator = 2 - 1j * phase_scale * admittance
            r = -1j * phase_scale * admittance / denominator
            t = 2 / denominator
    else:
        a = q0 if pol == "s" else eps * q0
        reflection = (a - q1) / (a + q1)
        one_way = complex(np.exp(1j * phase_scale * q1))
        round_trip = one_way * one_way
        denominator = 1 - reflection * reflection * round_trip
        r = reflection * (1 - round_trip) / denominator
        t = (1 - reflection * reflection) * one_way / denominator
    reflectance = float(abs(r) ** 2)
    transmittance = float(abs(t) ** 2)
    absorption = 1 - reflectance - transmittance
    # Remove floating-point roundoff only, never clip a material imbalance.
    if abs(absorption) < 32 * np.finfo(float).eps:
        absorption = 0.0
    if not all(math.isfinite(v) for v in (reflectance, transmittance, absorption)):
        raise ValueError("Fresnel result is not finite")
    return {
        "R": reflectance, "T": transmittance, "A": float(absorption),
        "r": complex(r), "t": complex(t), "pol": pol,
        "lambda_um": wavelength, "theta_deg": theta,
        "thickness_um": thickness,
    }


def periodic_projection(x, fields, period, kx, f, pol):
    """Decompose a full air-plane period into upward/downward Fourier orders.

    The incident Bloch component is normally ``kx=-f*sin(theta)``. Input fields
    retain their Bloch phase: each order is projected against
    exp(-2*pi*i*(kx+m/period)*x).  The arrays must be collocated at the same
    plane. Subtract reference *complex fields* before calling when appropriate.

    Both N equal-weight samples covering a period (including midpoint samples)
    and N+1 endpoint-inclusive samples are accepted. A duplicate endpoint is
    removed only after checking its Bloch phase. Near-cutoff orders are flagged
    and assigned no directional amplitude; no division by zero is attempted.
    """
    period = _real_scalar(period, "period", positive=True)
    kx = _real_scalar(kx, "kx")
    f = _real_scalar(f, "f", positive=True)
    pol = _polarization(pol)
    x = _array(x, "x", real=True).astype(float, copy=False)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("x must be a one-dimensional array with at least 2 samples")
    spacing = np.diff(x)
    dx = float(spacing[0])
    tolerance = max(period * 1e-10, 1e-13)
    if dx <= 0 or not np.allclose(spacing, dx, rtol=1e-9, atol=tolerance):
        raise ValueError("x must be strictly increasing and uniformly spaced")
    if not isinstance(fields, Mapping):
        raise ValueError("fields must map component names to collocated arrays")
    names = ("Ex", "Hz") if pol == "p" else ("Ez", "Hx")
    if any(name not in fields for name in names):
        raise ValueError(f"{pol} projection requires fields {names}")
    arrays = {name: _array(value, f"fields[{name}]") for name, value in fields.items()}
    if any(value.shape != x.shape for value in arrays.values()):
        raise ValueError("every field array must exactly match x.shape; no squeezing/broadcasting")
    duplicate_endpoint = math.isclose((x.size - 1) * dx, period, rel_tol=1e-9, abs_tol=tolerance)
    if duplicate_endpoint:
        phase = np.exp(2j * math.pi * kx * period)
        for name, values in arrays.items():
            scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
            if abs(values[-1] - phase * values[0]) > 1e-7 * scale:
                raise ValueError(f"{name} endpoints do not satisfy the specified Bloch phase")
        x = x[:-1]
        arrays = {name: values[:-1] for name, values in arrays.items()}
    elif not math.isclose(x.size * dx, period, rel_tol=1e-9, abs_tol=tolerance):
        raise ValueError("x must cover exactly one period without gaps or endpoint double counting")
    if x.size < 2:
        raise ValueError("at least two distinct periodic samples are required")
    n = x.size
    mode_numbers = np.arange(-(n // 2), (n - 1) // 2 + 1, dtype=int)
    needed_min = math.ceil((-f - kx) * period - 1e-10)
    needed_max = math.floor((f - kx) * period + 1e-10)
    if needed_min < int(mode_numbers[0]) or needed_max > int(mode_numbers[-1]):
        raise ValueError("x grid cannot resolve all propagating and cutoff diffraction orders")
    orders = []
    for m in mode_numbers:
        order_kx = kx + int(m) / period
        phase = np.exp(-2j * math.pi * order_kx * x)
        coefficient = {name: complex(np.mean(values * phase)) for name, values in arrays.items()}
        q_squared = 1 - (order_kx / f) ** 2
        cutoff = abs(q_squared) <= 1e-10
        propagating = q_squared > 1e-10
        q = complex(np.sqrt(complex(q_squared)))
        if cutoff:
            up, down = None, None
            up_power, down_power = 0.0, 0.0
        else:
            if pol == "p":
                up = (coefficient["Hz"] - coefficient["Ex"] / q) / 2
                down = (coefficient["Hz"] + coefficient["Ex"] / q) / 2
            else:
                up = (coefficient["Ez"] + coefficient["Hx"] / q) / 2
                down = (coefficient["Ez"] - coefficient["Hx"] / q) / 2
            up_power = period * q.real * abs(up) ** 2 if propagating else 0.0
            down_power = period * q.real * abs(down) ** 2 if propagating else 0.0
        orders.append({
            "m": int(m), "kx": float(order_kx), "q": float(q.real),
            "q_imag": float(q.imag), "propagating": bool(propagating),
            "cutoff": bool(cutoff), "up_amplitude": up, "down_amplitude": down,
            "up_power": float(up_power), "down_power": float(down_power),
        })
    up_power = float(sum(order["up_power"] for order in orders))
    down_power = float(sum(order["down_power"] for order in orders))
    if pol == "p":
        direct_flux = -period * np.mean(np.real(np.conj(arrays["Ex"]) * arrays["Hz"]))
    else:
        direct_flux = period * np.mean(np.real(np.conj(arrays["Ez"]) * arrays["Hx"]))
    return {
        "orders": orders, "up_power": up_power, "down_power": down_power,
        "net_up_power": up_power - down_power,
        "direct_net_up_power": float(direct_flux),
        "period": period, "samples": int(n), "dx": period / n,
        "duplicate_endpoint_removed": bool(duplicate_endpoint), "pol": pol,
        "cutoff_orders": [order["m"] for order in orders if order["cutoff"]],
        "normalization": "Meep DFT: period*q*abs(amplitude)^2; no 1/2",
    }


def _check_coordinates(coordinates, shape, name):
    if isinstance(coordinates, Mapping):
        values = list(coordinates.values())
        axes = list(coordinates.keys())
    elif isinstance(coordinates, (tuple, list)):
        values = list(coordinates)
        axes = [str(i) for i in range(len(values))]
    else:
        raise ValueError(f"{name} coordinates must contain one vector or grid per array axis")
    if len(values) != len(shape):
        raise ValueError(f"{name} coordinate count must equal field rank")
    output = {}
    for axis, (label, value) in enumerate(zip(axes, values)):
        coordinate = _array(value, f"{name} coordinate {label}", real=True)
        if coordinate.shape not in ((shape[axis],), shape):
            raise ValueError(f"{name} coordinate {label} does not match field shape")
        if coordinate.ndim == 1 and coordinate.size > 1 and not np.all(np.diff(coordinate) > 0):
            raise ValueError(f"{name} coordinate vectors must be strictly increasing")
        output[str(label)] = {"shape": list(coordinate.shape),
                              "min": float(np.min(coordinate)), "max": float(np.max(coordinate))}
    return output


def integrate_absorption(pairs, f, input_flux):
    """Integrate omega*Im(conj(E)*D) on explicitly collocated weighted grids.

    Each pair supplies name, E, D, weights, ti_mask, and coordinates (or coords).
    No interpolation, reshaping, broadcasting, calibration, sign correction or
    clipping is performed here. Weights already contain the physical area.
    ``input_flux`` is the strictly positive incident-power magnitude. ``A`` is
    the full integration-region absorption; ``A_ti`` and ``A_outside_ti`` expose
    the material-mask split for checking vacuum/interface contamination.
    """
    f = _real_scalar(f, "f", positive=True)
    input_flux = _real_scalar(input_flux, "input_flux", positive=True)
    if not isinstance(pairs, (list, tuple)) or not pairs:
        raise ValueError("pairs must be a nonempty list of collocated component pairs")
    omega = 2 * math.pi * f
    components = []
    seen = set()
    for pair in pairs:
        if not isinstance(pair, Mapping):
            raise ValueError("each absorption pair must be a mapping")
        required = ("name", "E", "D", "weights", "ti_mask")
        if any(key not in pair for key in required):
            raise ValueError(f"absorption pair must contain {required}")
        name = str(pair["name"])
        if not name or name in seen:
            raise ValueError("absorption component names must be nonempty and unique")
        seen.add(name)
        E = _array(pair["E"], f"{name}.E")
        D = _array(pair["D"], f"{name}.D")
        weights = _array(pair["weights"], f"{name}.weights", real=True)
        mask = _array(pair["ti_mask"], f"{name}.ti_mask", real=True)
        if E.ndim == 0 or any(value.shape != E.shape for value in (D, weights, mask)):
            raise ValueError(f"{name}: E,D,weights,ti_mask shapes must match exactly")
        if np.any(weights < 0) or not np.any(weights > 0):
            raise ValueError(f"{name}: weights must be nonnegative with positive total")
        if not np.all((mask == 0) | (mask == 1)):
            raise ValueError(f"{name}: ti_mask must be boolean or exactly 0/1")
        coordinate_value = pair.get("coordinates", pair.get("coords"))
        metadata = _check_coordinates(coordinate_value, E.shape, name)
        density = omega * np.imag(np.conj(E) * D)
        if not np.all(np.isfinite(density)):
            raise ValueError(f"{name}: absorption density overflowed")
        weighted = density * weights
        if not np.all(np.isfinite(weighted)):
            raise ValueError(f"{name}: weighted absorption density overflowed")
        ti_mask = mask.astype(bool)
        total = float(np.sum(weighted))
        ti = float(np.sum(weighted[ti_mask]))
        outside = float(np.sum(weighted[~ti_mask]))
        if not all(math.isfinite(value) for value in (total, ti, outside, float(np.sum(weights)))):
            raise ValueError(f"{name}: absorption integral overflowed")
        components.append({
            "name": name, "shape": list(E.shape), "coordinates": metadata,
            "weight_sum": float(np.sum(weights)), "ti_weight_sum": float(np.sum(weights[ti_mask])),
            "absorbed_power": total, "ti_absorbed_power": ti,
            "outside_ti_absorbed_power": outside,
            "negative_absorbed_power": float(np.sum(weighted[weighted < 0])),
            "positive_absorbed_power": float(np.sum(weighted[weighted > 0])),
            "density_min": float(np.min(density)), "density_max": float(np.max(density)),
            "A": total / input_flux,
        })
    total = float(sum(item["absorbed_power"] for item in components))
    ti = float(sum(item["ti_absorbed_power"] for item in components))
    outside = float(sum(item["outside_ti_absorbed_power"] for item in components))
    if not all(math.isfinite(value) for value in (total, ti, outside, total / input_flux,
                                                ti / input_flux, outside / input_flux)):
        raise ValueError("normalized absorption integral overflowed")
    return {
        "absorbed_power": total, "A": total / input_flux,
        "ti_absorbed_power": ti, "A_ti": ti / input_flux,
        "outside_ti_absorbed_power": outside, "A_outside_ti": outside / input_flux,
        "total_imag_E_conj_D": total / omega,
        "input_flux": input_flux, "frequency": f, "omega": omega,
        "components": components,
        "formula": "omega * sum(Im(conj(E)*D) * quadrature_weights) / input_flux",
        "empirical_calibration": False,
    }


class StopTracker:
    """Require decay and stable R/T in consecutive complete post-source windows.

    Probe and non-PML envelope inputs are intensities (sum of |E_component|²),
    not amplitudes or integrated fluxes. Each denominator is the larger of the
    channel's historical peak and a reference-excitation floor. By default the
    floor equals the reference intensity peak, so a probe that was never
    appreciably excited is still checked against the incident intensity budget;
    it is never dropped or divided by its own near-zero peak.

    Supply ``reference_intensity_peak`` from a qualified reference excitation
    and describe its origin in ``reference_provenance``. Without an external
    reference (e.g. while running that reference), the floor uses the historical
    non-PML envelope peak sampled at or before ``source_end`` and freezes after
    cutoff. A post-source tail cannot establish or increase this floor.
    ``reference_floor_fraction`` can tighten the floor, but cannot enlarge it
    beyond the supplied excitation. All intensities must have matching units
    and source scaling. A weak probe's allowed residual is bounded by
    ``decay_threshold * reference_floor_fraction * reference_intensity_peak``.

    R/T stability uses an absolute range, never a tail-relative drift. Callers
    enforce a time limit separately: it cannot set convergence.
    Samples should be at most half a window apart. Missing R/T observations,
    absent excitation, long gaps, and incomplete windows cannot pass.
    """

    def __init__(self, source_end=None, *, window=20, consecutive_windows=3,
                 decay_threshold=1e-4, rt_tolerance=1e-3,
                 reference_intensity_peak=None, reference_provenance=None,
                 reference_floor_fraction=1.0):
        self.source_end = None if source_end is None else _real_scalar(source_end, "source_end", nonnegative=True)
        self.window = _real_scalar(window, "window", positive=True)
        count = _real_scalar(consecutive_windows, "consecutive_windows", positive=True)
        if not count.is_integer():
            raise ValueError("consecutive_windows must be a positive integer")
        self.consecutive_windows = int(count)
        self.decay_threshold = _real_scalar(decay_threshold, "decay_threshold", positive=True)
        if self.decay_threshold >= 1:
            raise ValueError("decay_threshold must be less than one")
        self.rt_tolerance = _real_scalar(rt_tolerance, "rt_tolerance", positive=True)
        self.reference_floor_fraction = _real_scalar(
            reference_floor_fraction, "reference_floor_fraction", positive=True)
        if self.reference_floor_fraction > 1:
            raise ValueError("reference_floor_fraction must be at most one")
        self.reference_intensity_peak = None if reference_intensity_peak is None else _real_scalar(
            reference_intensity_peak, "reference_intensity_peak", positive=True)
        if self.reference_intensity_peak is not None:
            if not isinstance(reference_provenance, str) or not reference_provenance.strip():
                raise ValueError("an external reference intensity requires reference_provenance")
            _real_scalar(self.reference_intensity_peak * self.reference_floor_fraction,
                         "reference intensity denominator floor", positive=True)
        elif reference_provenance is not None:
            raise ValueError("reference_provenance requires reference_intensity_peak")
        self.reference_provenance = reference_provenance
        self.trace = []
        self._peaks = {}
        self._envelope_peak = 0.0
        self._source_on_envelope_peak = 0.0
        self._probe_names = None
        self._windows = []
        self.converged = False

    def _normalization(self):
        external = self.reference_intensity_peak is not None
        peak = self.reference_intensity_peak if external else self._source_on_envelope_peak
        floor = peak * self.reference_floor_fraction
        return {
            "method": "historical_intensity_peak_with_reference_excitation_floor",
            "quantity": "sum of |E_component|^2; not amplitude or integrated flux",
            "reference_intensity_peak": peak,
            "reference_peak_source": ("external_reference_excitation" if external
                                      else "source_on_non_pml_envelope"),
            "reference_provenance": (self.reference_provenance if external else
                                     "this trace: non-PML envelope samples at time <= source_end"),
            "reference_floor_fraction": self.reference_floor_fraction,
            "denominator_floor": floor,
            "weak_probe_residual_bound": self.decay_threshold * floor,
            "reference_ready": bool(floor > 0),
            "reference_frozen": bool(external or (self.trace and
                                      self.trace[-1]["time"] > self.source_end)),
            "all_probes_required": True,
        }

    def sample(self, time, probe_power, envelope_power, R=None, T=None, source_end=None):
        time = _real_scalar(time, "time", nonnegative=True)
        if self.trace and time <= self.trace[-1]["time"]:
            raise ValueError("sample times must be strictly increasing")
        if source_end is not None:
            candidate = _real_scalar(source_end, "source_end", nonnegative=True)
            if self.source_end is not None and candidate != self.source_end:
                raise ValueError("source_end cannot change during a trace")
            self.source_end = candidate
        if self.source_end is None:
            raise ValueError("source_end must be supplied before sampling")
        if not isinstance(probe_power, Mapping) or not probe_power:
            raise ValueError("probe_power must be a nonempty mapping of |E|² values")
        if any(not isinstance(name, str) or not name for name in probe_power):
            raise ValueError("probe names must be nonempty strings")
        names = frozenset(probe_power)
        if self._probe_names is not None and names != self._probe_names:
            raise ValueError("probe names must be unchanged throughout the trace")
        powers = {name: _real_scalar(value, f"probe_power[{name}]", nonnegative=True)
                  for name, value in probe_power.items()}
        envelope = _real_scalar(envelope_power, "envelope_power", nonnegative=True)
        R = None if R is None else _real_scalar(R, "R")
        T = None if T is None else _real_scalar(T, "T")
        self._probe_names = names
        for name, power in powers.items():
            self._peaks[name] = max(self._peaks.get(name, 0.0), power)
        self._envelope_peak = max(self._envelope_peak, envelope)
        if time <= self.source_end:
            self._source_on_envelope_peak = max(self._source_on_envelope_peak, envelope)
        normalization = self._normalization()
        normalization["reference_frozen"] = bool(self.reference_intensity_peak is not None or
                                                  time > self.source_end)
        floor = normalization["denominator_floor"]
        denominators = {name: max(peak, floor) for name, peak in self._peaks.items()}
        ratios = {name: power / denominators[name] if normalization["reference_ready"] else None
                  for name, power in powers.items()}
        # These labels describe observed excitation, not an assertion that a
        # small field is numerical noise. Every classified probe remains gated.
        classification = {
            name: ("zero_observed" if peak == 0 else
                   "reference_unavailable" if floor == 0 else
                   "below_reference_decay_budget" if peak <= self.decay_threshold * floor else
                   "resolved_excitation")
            for name, peak in self._peaks.items()
        }
        envelope_denominator = max(self._envelope_peak, floor)
        envelope_ratio = envelope / envelope_denominator if normalization["reference_ready"] else None
        record = {"time": time, "probe_power": powers, "probe_peak": dict(self._peaks),
                  "probe_ratio": ratios, "probe_normalization_denominator": denominators,
                  "probe_classification": classification,
                  "probe_normalization_source": {
                      name: "reference_excitation_floor" if peak < floor else "historical_probe_peak"
                      for name, peak in self._peaks.items()},
                  "normalization": normalization, "envelope_power": envelope,
                  "envelope_peak": self._envelope_peak, "envelope_ratio": envelope_ratio,
                  "envelope_normalization_denominator": envelope_denominator,
                  "R": R, "T": T,
                  "decay_pass": bool(normalization["reference_ready"] and self._envelope_peak > 0 and
                                     max([envelope_ratio, *ratios.values()]) <= self.decay_threshold)}
        self.trace.append(record)
        # Closed windows are evaluated once; historic-peak ratios are retained
        # as sampled, so a later burst cannot retroactively make old tails pass.
        complete = max(0, int(math.floor((time - self.source_end) / self.window + 1e-12)))
        while len(self._windows) < complete:
            index = len(self._windows)
            left = self.source_end + index * self.window
            right = left + self.window
            samples = [item for item in self.trace if left < item["time"] <= right + 1e-10]
            times = [left, *[item["time"] for item in samples], right]
            covered = len(samples) >= 2 and max(np.diff(times)) <= self.window / 2 + 1e-10
            r_values = [item["R"] for item in samples if item["R"] is not None]
            t_values = [item["T"] for item in samples if item["T"] is not None]
            r_range = float(max(r_values) - min(r_values)) if len(r_values) >= 2 else None
            t_range = float(max(t_values) - min(t_values)) if len(t_values) >= 2 else None
            have_rt = len(r_values) == len(samples) and len(t_values) == len(samples) and len(samples) >= 2
            decay = bool(samples and all(item["decay_pass"] for item in samples))
            stable = bool(have_rt and r_range < self.rt_tolerance and t_range < self.rt_tolerance)
            self._windows.append({"start": left, "end": right, "samples": len(samples),
                                  "covered": bool(covered), "decay_pass": decay,
                                  "R_range": r_range, "T_range": t_range,
                                  "rt_stable": stable, "passed": bool(covered and decay and stable)})
        last = self._windows[-self.consecutive_windows:]
        # Check the full three-window tail too: small per-window monotone drift
        # is not permission for a larger accumulated change in R or T.
        tail_start = last[0]["start"] if last else time
        tail = [item for item in self.trace if item["time"] > tail_start]
        all_rt = bool(tail and all(item["R"] is not None and item["T"] is not None for item in tail))
        tail_stable = bool(all_rt and
                           max(item["R"] for item in tail) - min(item["R"] for item in tail) < self.rt_tolerance and
                           max(item["T"] for item in tail) - min(item["T"] for item in tail) < self.rt_tolerance)
        self.converged = bool(len(last) == self.consecutive_windows and
                              all(item["passed"] for item in last) and tail_stable and
                              all(item["decay_pass"] for item in tail))
        record["converged"] = self.converged
        return self.converged

    def summary(self):
        consecutive = 0
        for window in reversed(self._windows):
            if not window["passed"]:
                break
            consecutive += 1
        return {
            "converged": self.converged,
            "status": "converged" if self.converged else "not_converged",
            "source_end": self.source_end, "last_time": self.trace[-1]["time"] if self.trace else None,
            "window": self.window, "required_windows": self.consecutive_windows,
            "consecutive_passed_windows": consecutive, "decay_threshold": self.decay_threshold,
            "rt_tolerance": self.rt_tolerance, "probe_peaks": dict(self._peaks),
            "non_pml_envelope_peak": self._envelope_peak,
            "source_on_non_pml_envelope_peak": self._source_on_envelope_peak,
            "windows": [dict(item) for item in self._windows],
            "last_sample": dict(self.trace[-1]) if self.trace else None,
            "normalization": self._normalization(),
            "timeout_is_convergence": False,
        }


def direction_metrics(eplus, eminus, errorplus=None, errorminus=None, *, coarse=False):
    """Signed directional contrast with an explicit five-error signal floor.

    D=(eplus-eminus)/(eplus+eminus). Bounds are deterministic absolute error
    estimates, not assumed independent statistical standard deviations. A
    missing error estimate or unresolved denominator leaves D unverified.
    FWHM cannot be inferred from a pair of samples and is always returned as
    None; coarse=True explicitly records an unresolved coarse-grid width.
    """
    eplus = _real_scalar(eplus, "eplus", nonnegative=True)
    eminus = _real_scalar(eminus, "eminus", nonnegative=True)
    denominator = eplus + eminus
    raw = (eplus - eminus) / denominator if denominator > 0 else None
    result = {
        "eplus": eplus, "eminus": eminus, "sum": denominator,
        "raw_D": raw, "D": None, "error_D": None,
        "status": "not_verified", "verified": False,
        "denominator_floor": None, "FWHM": None,
        "FWHM_status": "unresolved_coarse_grid" if coarse else "not_evaluated",
        "unresolved_FWHM": bool(coarse),
    }
    if errorplus is None or errorminus is None:
        result["reason"] = "independent_error_estimates_required"
        return result
    errorplus = _real_scalar(errorplus, "errorplus", nonnegative=True)
    errorminus = _real_scalar(errorminus, "errorminus", nonnegative=True)
    error_sum = errorplus + errorminus
    floor = 5 * error_sum
    result.update(errorplus=errorplus, errorminus=errorminus, denominator_floor=floor)
    if denominator <= floor:
        result["reason"] = "denominator_not_above_five_error_floor"
        return result
    # Conservative finite-error bound; denominator-error coupling retained.
    uncertainty = 2 * (eminus * errorplus + eplus * errorminus) / (denominator * (denominator - error_sum))
    result.update(D=raw, error_D=float(uncertainty), status="verified", verified=True,
                  reason="independent_error_floor_passed")
    return result
