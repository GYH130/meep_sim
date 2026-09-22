"""Offline TM Yee half-space predictions, not FDTD runs or convergence evidence.

Assumes a flat interface on an Ex node, air above, and homogeneous passive
material below. It does not qualify a finite film or structured geometry.
Meep averages epsilon_inf only: https://meep.readthedocs.io/en/stable/Subpixel_Smoothing/
ADE recurrence: https://github.com/NanoComp/meep/blob/master/src/susceptibility.cpp
"""
import cmath
import math


def _real(value, name, minimum=0.0, strict=True):
    value = float(value)
    if not math.isfinite(value) or (value <= minimum if strict else value < minimum):
        raise ValueError(f"Invalid {name}")
    return value


def discrete_epsilon(model, f, dt):
    """Independent centered-ADE epsilon; dt=0 gives the continuum pole model."""
    f, dt = _real(f, "frequency"), _real(dt, "dt", strict=False)
    base = _real(model["epsilon"], "epsilon_inf")
    if model.get("model", "susceptibility") != "susceptibility":
        raise ValueError("Expected susceptibility model")
    if model.get("meep_length_unit_um", 1.0) != 1.0:
        raise ValueError("Expected one length unit = one micrometer")
    w = 2 * math.pi * f
    if w * dt >= math.pi:
        raise ValueError("Frequency exceeds temporal Nyquist limit")
    omega = 2 * math.sin(w * dt / 2) / dt if dt else w
    damping_frequency = math.sin(w * dt) / dt if dt else w
    result = complex(base)
    for pole in model.get("susceptibilities", []):
        kind = str(pole["type"]).lower()
        if kind not in ("lorentz", "drude"):
            raise ValueError("Unsupported susceptibility")
        w0 = 2 * math.pi * _real(pole["frequency"], "pole frequency")
        gamma = 2 * math.pi * _real(pole["gamma"], "gamma", strict=False)
        sigma = _real(pole["sigma"], "sigma", strict=False)
        if kind == "lorentz" and w0 * dt >= 2:
            raise ValueError("Unstable Lorentz ADE timestep")
        denominator = (w0 * w0 if kind == "lorentz" else 0) - omega**2
        denominator -= 1j * gamma * damping_frequency
        if denominator == 0:
            raise ValueError("Singular undamped material resonance")
        result += sigma * w0 * w0 / denominator
    if not all(math.isfinite(v) for v in (result.real, result.imag)) or result.imag < 0:
        raise ValueError("Expected finite passive epsilon")
    return result


def planar_tm_prediction(model, wavelength_um, theta_deg, resolution, courant=.25,
                         interface_mode="half_instantaneous", interface_epsilon=None):
    """Predict R from discrete equations, independently of saved fields by default.

    Ex_j and Hz_(j+1/2) use depth into the metal and exp(-i omega t).
    half_instantaneous represents the interface in the inspected Meep setup;
    bulk disables instantaneous averaging, and full_half averages all epsilon.
    An explicit interface_epsilon overrides that constitutive value. Predictions
    at r16/r20/r24 or other grids are not FDTD results or convergence certificates.
    """
    wave = _real(wavelength_um, "wavelength")
    h = 1 / _real(resolution, "resolution")
    courant = _real(courant, "Courant")
    theta = float(theta_deg)
    if courant > .5 or not math.isfinite(theta) or abs(theta) >= 90:
        raise ValueError("Require Courant <= 0.5 and abs(theta) < 90 degrees")
    dt, w = courant * h, 2 * math.pi / wave
    eps = discrete_epsilon(model, 1 / wave, dt)
    omega = 2 * math.sin(w * dt / 2) / dt
    kx = 2 * math.sin(-w * math.sin(math.radians(theta)) * h / 2) / h
    vacuum_argument = h * h * (omega * omega - kx * kx) / 4
    if not 0 < vacuum_argument < 1:
        raise ValueError("Incident vacuum mode is not propagating on this grid")
    a = cmath.exp(2j * math.asin(math.sqrt(vacuum_argument)))
    p = cmath.exp(2j * cmath.asin(h * cmath.sqrt(omega**2 * eps - kx**2) / 2))
    if abs(p) > 1 + 1e-12:
        raise ValueError("Growing branch is not a passive half-space solution")
    modes = {"half_instantaneous": eps - (model["epsilon"] - 1) / 2,
             "bulk": eps, "full_half": (eps + 1) / 2}
    if interface_mode not in modes:
        raise ValueError("Unknown interface mode")
    eps_i = modes[interface_mode] if interface_epsilon is None else complex(interface_epsilon)
    if not all(math.isfinite(v) for v in (eps_i.real, eps_i.imag)) or eps_i.imag < 0:
        raise ValueError("Expected finite passive interface epsilon")
    if eps == 0 or omega**2 * eps == kx**2:
        raise ValueError("Degenerate material admittance")
    bm, ba = omega - kx**2 / (omega * eps), omega - kx**2 / omega
    boundary = (p - 1) / bm + omega * h * h * eps_i
    r = (1 - 1 / a - ba * boundary) / (ba * boundary - (1 - a))
    return dict(R=abs(r)**2, r=r, depth_factor=p, epsilon_discrete=eps,
                epsilon_interface=eps_i, vacuum_depth_factor=a,
                omega_discrete=omega, kx_discrete=kx, prediction_only=True,
                fdtd_run=False, convergence_certified=False)
