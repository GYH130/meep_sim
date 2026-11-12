#!/usr/bin/env python3
"""Score saved ray summaries using normal and off-axis collection cones.

The sampled angular fields remain in the output. Interpolation and integration
are post-processing, with no additional ray solves or spectral assumptions.
Only the Python standard library is needed.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
import json
import math
import os
from pathlib import Path
from tempfile import NamedTemporaryFile


NORMAL_ANGLES = tuple(range(5, 86))
OFF_AXIS_ANGLES = (5, 10, 15, 20, 25, 30, 45, 60)
TARGET_AXES = tuple((theta, phi) for theta in (15, 30, 45) for phi in (0, 90, 180, 270))
PERIOD = 2 * math.pi
HORIZON = math.pi / 2
EPS = 1e-10
GL8 = (
    (-0.9602898564975363, 0.1012285362903763),
    (-0.7966664774136267, 0.2223810344533745),
    (-0.5255324099163290, 0.3137066458778873),
    (-0.1834346424956498, 0.3626837833783620),
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def projected_coordinate(angle_deg: float) -> float:
    return math.sin(math.radians(angle_deg)) ** 2


@dataclass
class RingIntegral:
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        self.step = PERIOD / len(self.values)
        self.closed = self.values + self.values[:1]
        self.prefix = [0.0]
        for y0, y1 in zip(self.closed, self.closed[1:]):
            self.prefix.append(self.prefix[-1] + self.step * (y0 + y1) / 2)

    def primitive(self, phi: float) -> float:
        """Integral of the periodic piecewise-linear field from 0 to phi."""
        turns = math.floor(phi / PERIOD)
        within = phi - turns * PERIOD
        index = min(int(within / self.step), len(self.values) - 1)
        dx = within - index * self.step
        y0, y1 = self.closed[index:index + 2]
        return turns * self.prefix[-1] + self.prefix[index] + dx * (
            y0 + dx * (y1 - y0) / (2 * self.step)
        )

    def arc(self, phi_low: float, phi_high: float) -> float:
        return self.primitive(phi_high) - self.primitive(phi_low)

    @property
    def mean(self) -> float:
        return math.fsum(self.values) / len(self.values)


def normal_integrals(theta: tuple[float, ...], rings: tuple[RingIntegral, ...],
                     radii: tuple[float, ...]) -> dict[float, float]:
    """Exact integral of azimuth means linearly interpolated in sin²(theta)."""
    u = tuple(map(projected_coordinate, theta))
    means = tuple(ring.mean for ring in rings)
    if theta[-1] < 90:
        u += (1.0,)
        means += means[-1:]
    prefix = [0.0]
    for x0, x1, y0, y1 in zip(u, u[1:], means, means[1:]):
        prefix.append(prefix[-1] + (x1 - x0) * (y0 + y1) / 2)
    result = {}
    for alpha in radii:
        x = projected_coordinate(alpha)
        index = min(bisect_right(u, x) - 1, len(u) - 2)
        dx = x - u[index]
        slope = (means[index + 1] - means[index]) / (u[index + 1] - u[index])
        result[alpha] = prefix[index] + dx * (means[index] + dx * slope / 2)
    return result


@lru_cache(maxsize=4096)
def cap_quadrature(theta: tuple[float, ...], phi_count: int, axis_theta: float,
                   axis_phi: float, radius: float) -> tuple[tuple, float]:
    """Return interpolation/arc terms and projected cap weight in the upper hemisphere.

    The spherical cap has a contiguous periodic azimuth interval at each theta.
    Resolve its edges analytically, including where they cross the input phi
    nodes. This avoids assigning whole grid cells to one side of a cap edge.
    """
    require(0 < axis_theta <= 90 and 0 < radius <= 90, "Invalid off-axis cap")
    samples = tuple(map(math.radians, theta))
    theta_c, phi_c, alpha = map(math.radians, (axis_theta, axis_phi, radius))
    cos_c, sin_c, cos_a = math.cos(theta_c), math.sin(theta_c), math.cos(alpha)
    low = max(0.0, theta_c - alpha)
    high = min(HORIZON, theta_c + alpha)
    borders = [low, high, *(t for t in samples if low < t < high)]
    if low < alpha - theta_c < high:
        borders.append(alpha - theta_c)
    for index in range(phi_count):
        phi = index * PERIOD / phi_count
        b = sin_c * math.cos(phi - phi_c)
        magnitude = math.hypot(cos_c, b)
        if magnitude <= EPS or cos_a >= magnitude:
            continue
        phase = math.atan2(b, cos_c)
        offset = math.acos(cos_a / magnitude)
        borders.extend(t for t in (phase - offset, phase + offset) if low < t < high)
    borders = sorted(set(borders))

    terms = []
    projected_weights = []
    for lo, hi in zip(borders, borders[1:]):
        if hi - lo <= 1e-12:
            continue
        center, half = (lo + hi) / 2, (hi - lo) / 2
        # A sine transformation smooths the square-root arc endpoints.
        for node, weight in GL8:
            eta = HORIZON * node
            angle = center + half * math.sin(eta)
            sin_t, cos_t = math.sin(angle), math.cos(angle)
            limit = (cos_a - cos_t * cos_c) / (sin_t * sin_c)
            if limit >= 1:
                continue
            arc_half = math.pi if limit <= -1 else math.acos(limit)
            w = weight * half * HORIZON * math.cos(eta) * sin_t * cos_t / math.pi
            lower = min(max(0, bisect_right(samples, angle) - 1), len(samples) - 1)
            upper = min(lower + 1, len(samples) - 1)
            if lower == upper:
                upper_fraction = 0.0
            else:
                u_low = math.sin(samples[lower]) ** 2
                u_high = math.sin(samples[upper]) ** 2
                upper_fraction = (sin_t ** 2 - u_low) / (u_high - u_low)
            terms.append((lower, upper, w * (1 - upper_fraction), w * upper_fraction,
                          phi_c - arc_half, phi_c + arc_half))
            projected_weights.append(w * 2 * arc_half)
    projected_weight = math.fsum(projected_weights)
    if theta_c + alpha <= HORIZON + 1e-12:
        analytic = cos_c * math.sin(alpha) ** 2
        require(math.isclose(projected_weight, analytic, rel_tol=5e-7, abs_tol=1e-9),
                f"Spherical-cap quadrature disagrees with analytic projected area: {axis_theta}, {axis_phi}, {radius}")
    return tuple(terms), projected_weight


def cap_integral(rings: tuple[RingIntegral, ...], terms: tuple) -> float:
    return math.fsum(w0 * rings[i0].arc(phi0, phi1) + w1 * rings[i1].arc(phi0, phi1)
                     for i0, i1, w0, w1, phi0, phi1 in terms)


def cone_metrics(power: float, power_upper: float, hemisphere: float,
                 hemisphere_upper: float, projected_weight: float) -> dict:
    require(0 < projected_weight < 1 + EPS, "Invalid projected cap weight")
    require(-EPS <= power <= hemisphere + EPS and power <= power_upper + EPS,
            "Cone power outside its physical bounds")
    outside = max(0.0, hemisphere - power)
    outside_upper = max(outside, hemisphere_upper - power_upper)
    inside_mean = power / projected_weight
    outside_mean = outside / (1 - projected_weight) if projected_weight < 1 - EPS else None
    eta = power / hemisphere if hemisphere > 0 else None
    lower_denominator = power + outside_upper
    upper_denominator = power_upper + outside
    return {
        "projected_blackbody_weight": projected_weight,
        "cone_power_over_full_hemisphere_blackbody": power,
        "cone_average_ideal_emissivity": inside_mean,
        "cone_fraction_of_emitted_power": eta,
        "cone_outside_power_over_full_hemisphere_blackbody": outside,
        "cone_outside_fraction_of_emitted_power": outside / hemisphere if hemisphere > 0 else None,
        "cone_outside_mean_ideal_emissivity": outside_mean,
        "cone_to_outside_mean_ratio": inside_mean / outside_mean if outside_mean and outside_mean > 1e-6 else None,
        "unresolved_only_power_interval": [power, power_upper],
        "unresolved_only_emitted_fraction_interval": [
            power / lower_denominator if lower_denominator > 0 else None,
            power_upper / upper_denominator if upper_denominator > 0 else None,
        ],
    }


def validate_summary(data: dict, path: Path) -> tuple[str, tuple[str, ...]]:
    require(str(data.get("state", "")).startswith("COMPLETE") and not data.get("failed"),
            f"{path}: only complete, failure-free summaries are scored")
    require(not data.get("not_dispatched"), f"{path}: missing undispatched directions")
    rows = data["geometry_summaries"]
    require(len(rows) == data["completed_geometries"] == data["planned_geometry_count"],
            f"{path}: geometry counts do not match")
    count = sum(row["completed_angles"] for row in rows)
    require(count == data["completed_conditions"] == data["condition_count"],
            f"{path}: direction counts do not match")
    columns = tuple(data["azimuth_sample_columns"])
    require(len(columns) == len(set(columns)) and {"phi_deg", "A_total", "residual"} <= set(columns),
            f"{path}: invalid azimuth sample columns")
    return str(data["phase"]), columns


def unpack_row(row: dict, columns: tuple[str, ...]) -> tuple:
    shape = row["shape"]
    name = shape["id"]
    angles = sorted(row["angles"], key=lambda ring: ring["theta_deg"])
    theta = tuple(float(ring["theta_deg"]) for ring in angles)
    require(len(theta) >= 3 and theta[0] == 0 and 89 <= theta[-1] <= 90
            and all(a < b for a, b in zip(theta, theta[1:])), f"{name}: incomplete polar grid")
    azimuth = [list(ring["azimuth_samples"]) for ring in angles]
    phi_column = columns.index("phi_deg")
    phi_count = len(azimuth[1])
    require(phi_count >= 3, f"{name}: incomplete azimuth grid")
    a_index, r_index = columns.index("A_total"), columns.index("residual")
    field, upper_field = [], []
    partition = {"A_active", "A_sidewall", "A_land"}
    identity_columns = {key: columns.index(key) for key in sorted(partition) if key in columns}
    for index, (ring, samples) in enumerate(zip(angles, azimuth)):
        expected = [0.0] if index == 0 else [k * 360 / phi_count for k in range(phi_count)]
        require(len(samples) == len(expected) == ring["azimuth_count"], f"{name}: missing azimuth samples")
        absorption, upper = [], []
        for sample, phi in zip(samples, expected):
            require(len(sample) == len(columns) and all(math.isfinite(v) for v in sample)
                    and math.isclose(sample[phi_column], phi, abs_tol=1e-8), f"{name}: invalid azimuth sample")
            a, residual = sample[a_index], sample[r_index]
            require(-EPS <= a <= 1 + EPS and -EPS <= residual and a + residual <= 1 + EPS,
                    f"{name}: directional energy bounds invalid")
            if identity_columns.keys() == partition:
                values = [sample[i] for i in identity_columns.values()]
                require(all(-EPS <= v <= 1 + EPS for v in values)
                        and math.isclose(a, math.fsum(values), rel_tol=0, abs_tol=1e-8),
                        f"{name}: active/sidewall/land absorptance does not sum to total")
            absorption.append(a)
            upper.append(a + residual)
        require(math.isclose(math.fsum(absorption) / len(absorption), ring["mean_A"], abs_tol=1e-9),
                f"{name}: directional mean does not match retained samples")
        if index == 0:
            absorption *= phi_count
            upper *= phi_count
        field.append(RingIntegral(tuple(absorption)))
        upper_field.append(RingIntegral(tuple(upper)))
    expected_count = 1 + (len(theta) - 1) * phi_count
    require(row["completed_angles"] == row["expected_angles"] == expected_count,
            f"{name}: incomplete geometry")
    return theta, tuple(field), tuple(upper_field), azimuth, identity_columns


def score_row(row: dict, columns: tuple[str, ...], radii: tuple[float, ...] = NORMAL_ANGLES,
              off_axis_radii: tuple[float, ...] = OFF_AXIS_ANGLES) -> dict:
    shape = row["shape"]
    theta, field, upper, azimuth, regions = unpack_row(row, columns)
    phi_count = len(field[0].values)
    total = normal_integrals(theta, field, (*radii, 90))
    total_upper = normal_integrals(theta, upper, (*radii, 90))
    hemi, hemi_upper = total[90], total_upper[90]
    if all(a.values == b.values for a, b in zip(field, upper)):
        upper = field

    def sample_peak(projected: bool) -> dict:
        value, negative_polar, negative_index = max(
            (a * (math.cos(math.radians(t)) if projected else 1), -t, -i)
            for t, ring in zip(theta, field) for i, a in enumerate(ring.values))
        polar, azimuth_index = -negative_polar, -negative_index
        return {"value": value, "theta_deg": polar,
                "phi_deg": azimuth_index * 360 / phi_count if polar else None}

    regional_metrics = {}
    for channel, index in regions.items():
        rings = tuple(RingIntegral(tuple(s[index] for s in ring) * (phi_count if i == 0 else 1))
                      for i, ring in enumerate(azimuth))
        integrated = normal_integrals(theta, rings, (*radii, 90))
        regional_metrics[channel] = {
            "hemisphere_power_over_full_hemisphere_blackbody": integrated[90],
            "normal_cone_power_over_full_hemisphere_blackbody": {
                f"{angle:g}": integrated[angle] for angle in radii},
        }

    off_axis = []
    for axis_theta, axis_phi in TARGET_AXES:
        values = {}
        for radius in off_axis_radii:
            terms, weight = cap_quadrature(theta, phi_count, axis_theta, axis_phi, radius)
            power = cap_integral(field, terms)
            power_upper = power if upper is field else cap_integral(upper, terms)
            values[f"{radius:g}"] = cone_metrics(power, power_upper, hemi, hemi_upper, weight)
        if values:
            off_axis.append({"axis_theta_deg": axis_theta, "axis_phi_deg": axis_phi, "cones": values})

    q_lo = next((a for a in radii if hemi > 0 and total[a] >= 0.5 * hemi), None)
    q_hi = next((a for a in radii if hemi > 0 and total[a] >= 0.8 * hemi), None)
    return {
        "shape": shape,
        "sampled_theta_deg": theta,
        "azimuth_mean": [ring.mean for ring in field],
        "azimuth_samples": azimuth,
        "azimuth_sample_columns": columns,
        "hemisphere_power_over_full_hemisphere_blackbody": hemi,
        "unresolved_only_hemisphere_power_interval": [hemi, hemi_upper],
        "unsampled_horizon_absolute_power_bound": 1 - projected_coordinate(theta[-1]),
        "sampled_radiance_peak": sample_peak(False),
        "sampled_projected_power_peak": sample_peak(True),
        "first_reported_normal_cone_radius_at_50_pct_deg": q_lo,
        "first_reported_normal_cone_radius_at_80_pct_deg": q_hi,
        "regional_power": regional_metrics,
        "normal_cones": {f"{a:g}": cone_metrics(total[a], total_upper[a], hemi, hemi_upper,
                                                projected_coordinate(a)) for a in radii},
        "off_axis_target_cones": off_axis,
    }


def compare_surfaces(records: list[dict]) -> list[dict]:
    """Compare distinct surface profiles only within one identical geometry."""
    by_shape = {}
    for record in records:
        if record["surface_profile"] is not None:
            by_shape.setdefault(record["shape"]["id"], []).append(record)
    comparisons = []
    for sid, group in by_shape.items():
        if len(group) < 2:
            continue
        baseline = next((row for row in group if row["phase"] == "ideal"), None)
        if baseline is None:
            continue
        for row in group:
            if row is baseline:
                continue
            require(row["shape"] == baseline["shape"], f"{sid}: surface runs have different geometry")
            comparisons.append({"shape_id": sid, "baseline_phase": baseline["phase"],
                                "comparison_phase": row["phase"],
                                "hemisphere_power_delta": row["hemisphere_power_over_full_hemisphere_blackbody"]
                                - baseline["hemisphere_power_over_full_hemisphere_blackbody"],
                                "normal_cone_power_delta": {
                                    key: row["normal_cones"][key]["cone_power_over_full_hemisphere_blackbody"]
                                    - base["cone_power_over_full_hemisphere_blackbody"]
                                    for key, base in baseline["normal_cones"].items()}})
    return comparisons


def make_report(paths: list[Path], off_axis_radii: tuple[float, ...]) -> dict:
    records, sources, case_keys = [], [], set()
    for path in paths:
        raw = path.read_bytes()
        data = json.loads(raw)
        phase, columns = validate_summary(data, path)
        source_index = len(sources)
        sources.append({"path": str(path.resolve()), "sha256": sha256(raw).hexdigest(),
                        "phase": phase, "cases": data["completed_geometries"],
                        "directions": data["completed_conditions"],
                        "rays_per_direction": data["rays_per_condition"],
                        "surface_profile": data.get("surface_profile")})
        for row in data["geometry_summaries"]:
            key = (phase, row["shape"]["id"])
            require(key not in case_keys, f"Duplicate geometry/phase: {key}")
            case_keys.add(key)
            scored = score_row(row, columns, off_axis_radii=off_axis_radii)
            scored.update(phase=phase, source_index=source_index, surface_profile=data.get("surface_profile"))
            records.append(scored)
    return {
        "kind": "SAVED_GEOMETRIC_OPTICS_ANGULAR_POSTPROCESSING",
        "scope": "Ideal or prescribed gray geometric optics. No wavelength, temperature, physical TA2 or full-wave prediction.",
        "methods": {
            "normal": "Piecewise linear in sin²(theta) of full-period uniform azimuth means, analytically integrated.",
            "off_axis": "Spherical cap clipped to upper hemisphere. Exact periodic piecewise-linear azimuth arc integrals, interpolation in sin²(theta), and Gaussian polar quadrature with the actual surface normal projection.",
            "normalization": "Power divided by a unit-cell footprint blackbody emitting into the whole upper hemisphere. Cone mean divides this by the same cone's projected blackbody fraction. Cone concentration instead divides by that structure's hemisphere output.",
            "horizon": "Hold the last sampled response to 90 degrees. Absolute normalized horizon uncertainty is bounded separately by its projected area fraction.",
            "uncertainty": "Unresolved-only intervals do not include Monte Carlo, angular sampling, geometric facets, material or interpolation errors. Near-zero outside mean (<=1e-6) is reported through leakage, not a ratio.",
            "sampling": "Off-axis metrics are quadratures of interpolants of retained samples, not additional ray evaluations. Sampled peak coordinates do not locate a unique optimum; tied plateaus can have multiple directions.",
        },
        "normal_half_angles_deg": NORMAL_ANGLES,
        "off_axis_half_angles_deg": off_axis_radii,
        "sources": sources,
        "cases": records,
        "matched_surface_comparisons": compare_surfaces(records),
    }


def write_json(path: Path, content: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile("w", dir=path.parent, encoding="utf-8", delete=False) as file:
            temporary = Path(file.name)
            json.dump(content, file, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
            file.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summaries", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--off-axis-half-angles", nargs="+", type=float, default=OFF_AXIS_ANGLES)
    parser.add_argument("--normal-only", action="store_true", help="Skip off-axis quadratures; retain all azimuth samples")
    args = parser.parse_args()
    radii = () if args.normal_only else tuple(sorted(set(args.off_axis_half_angles)))
    if any(not math.isfinite(a) or not 0 < a <= 90 for a in radii):
        parser.error("off-axis half-angles must be finite and lie in (0, 90] degrees")
    report = make_report(args.summaries, radii)
    write_json(args.output, report)
    print(json.dumps({"state": "COMPLETE_ANGULAR_POSTPROCESSING", "cases": len(report["cases"]),
                      "sources": len(report["sources"]), "output": str(args.output.resolve())}))


if __name__ == "__main__":
    main()
