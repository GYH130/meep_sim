#!/usr/bin/env python3
"""Reduce completed angular post-processing to matched family comparisons."""
from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
import json
import math
from pathlib import Path

from analyze_results import RingIntegral, cap_integral, cap_quadrature, require, write_json


PHASES = ("ideal", "floor90", "wall04")
ANCHORS = ("G23", "G99", "MID")
RADII = (15, 30, 45)
POWER = "cone_power_over_full_hemisphere_blackbody"
FRACTION = "cone_fraction_of_emitted_power"
MEAN = "cone_average_ideal_emissivity"
HEMISPHERE = "hemisphere_power_over_full_hemisphere_blackbody"


def row_index(report: dict) -> dict[tuple[str, str], dict]:
    records = report["cases"]
    indexed = {(row["phase"], row["shape"]["id"]): row for row in records}
    ids = {row["shape"]["id"] for row in records}
    require(len(ids) == 24 and len(records) == len(indexed) == 72
            and indexed.keys() == set(product(PHASES, ids)), "Expected 24 matched geometries in each phase")
    for sid in sorted(ids):
        group = [indexed[phase, sid] for phase in PHASES]
        require(all(row["shape"] == group[0]["shape"] for row in group[1:]),
                f"Geometry changed between surface phases: {sid}")
    return indexed


def sample_fields(row: dict, column: str) -> tuple[tuple[float, ...], ...]:
    index = row["azimuth_sample_columns"].index(column)
    return tuple(tuple(sample[index] for sample in ring) for ring in row["azimuth_samples"])


def cone_axis_x(row: dict, radius: int = 15) -> dict:
    """Optimize receiver power on an interpolated signed x-z angular cut only."""
    theta = tuple(row["sampled_theta_deg"])
    samples = sample_fields(row, "A_total")
    count = len(samples[1])
    rings = (RingIntegral(samples[0] * count), *(RingIntegral(values) for values in samples[1:]))
    power_by_tilt = {}
    for tilt in range(-35, 36):
        if tilt == 0:
            power = row["normal_cones"][str(radius)][POWER]
        else:
            terms, _ = cap_quadrature(theta, count, abs(tilt), 0 if tilt > 0 else 180, radius)
            power = cap_integral(rings, terms)
        power_by_tilt[tilt] = power
    best = max(power_by_tilt, key=lambda tilt: (power_by_tilt[tilt], -abs(tilt), tilt))
    normal = power_by_tilt[0]
    return {
        "half_angle_deg": radius,
        "best_signed_receiver_tilt_deg": best,
        "best_axis_theta_deg": abs(best),
        "best_axis_phi_deg": 0 if best > 0 else (180 if best < 0 else None),
        "power_at_best_axis": power_by_tilt[best],
        "fraction_at_best_axis": power_by_tilt[best] / row[HEMISPHERE] if row[HEMISPHERE] else None,
        "normal_axis_power": normal,
        "gain_vs_normal_relative": power_by_tilt[best] / normal - 1 if normal else None,
        "power_by_signed_receiver_tilt": power_by_tilt,
    }


def ranking(indexed: dict[tuple[str, str], dict]) -> dict:
    groups = defaultdict(list)
    for (phase, _), row in indexed.items():
        groups[phase, row["shape"]["anchor"]].append(row)
    report = {}
    for (phase, anchor), group in sorted(groups.items()):
        require(anchor in ANCHORS and len(group) == 8, f"Missing family variants: {phase}/{anchor}")
        result = {}
        for alpha in RADII:
            def get(row):
                return row["normal_cones"][str(alpha)]

            def rank(key):
                return [row["shape"]["variant"] for row in sorted(
                    group, key=lambda row: (-round(get(row)[key], 12), -get(row)[POWER], row["shape"]["variant"]))]

            result[str(alpha)] = {"q_rank": rank(POWER), "eta_rank": rank(FRACTION)}
        report.setdefault(phase, {})[anchor] = result
    return report


def groove_invariant_check(row: dict) -> dict:
    """Independently integrate a y-invariant, x-symmetric gray groove over beta.

    With tangential direction cosine s_y, dΩ = dβ ds_y and cosθ equals
    sqrt(1-s_y²) cosβ. The integral in s_y is π/2. Reflection symmetry in x
    gives the normalized hemisphere integral ∫₀¹ A(arcsin u) du.
    """
    require(row["shape"]["variant"] == "T0" and row["shape"]["y_invariant"],
            "Invariant quadrature requires a symmetric infinite groove")
    indices = {name: i for i, name in enumerate(row["azimuth_sample_columns"])}
    values = defaultdict(list)
    for theta_deg, samples in zip(row["sampled_theta_deg"], row["azimuth_samples"]):
        theta = math.radians(theta_deg)
        for sample in samples:
            phi = math.radians(sample[indices["phi_deg"]])
            beta = math.atan2(math.sin(theta) * math.cos(phi), math.cos(theta))
            coordinate = round(abs(math.sin(beta)), 10)
            values[coordinate].append(sample[indices["A_total"]])
    nodes = sorted((u, math.fsum(samples) / len(samples)) for u, samples in values.items())
    nodes.append((1.0, nodes[-1][1]))
    invariant = math.fsum((u1 - u0) * (a0 + a1) / 2 for (u0, a0), (u1, a1) in zip(nodes, nodes[1:]))
    original = row[HEMISPHERE]
    return {"angular_grid_hemisphere_q": original, "invariant_hemisphere_q": invariant,
            "invariant_minus_grid": invariant - original, "floor_area_fraction": row["shape"]["floor_area_fraction"],
            "max_duplicate_coordinate_spread": max(max(samples) - min(samples) for samples in values.values()),
            "beta_nodes_before_horizon_extension": len(values)}


def build_summary(report: dict) -> dict:
    indexed = row_index(report)
    cases = {}
    checks = defaultdict(dict)
    centered = {(phase, anchor): cone_axis_x(indexed[phase, f"{anchor}_A0"])
                for phase, anchor in product(PHASES, ANCHORS)}
    for (phase, sid), row in sorted(indexed.items()):
        shape = row["shape"]
        record = {
            "shape": shape,
            "normal_direction_emissivity": row["azimuth_mean"][0],
            "hemisphere_q": row[HEMISPHERE],
            "normal_cones": {str(alpha): {key: row["normal_cones"][str(alpha)][source]
                                          for key, source in (("q", POWER), ("eta", FRACTION), ("epsilon_inside", MEAN))}
                             for alpha in RADII},
            "radiance_peak_sampled": row["sampled_radiance_peak"],
            "projected_power_peak_sampled": row["sampled_projected_power_peak"],
            "hemisphere_region_q": {region: values[HEMISPHERE] for region, values in row["regional_power"].items()},
        }
        if "dxplus" in sid or "dxminus" in sid:
            sweep = cone_axis_x(row)
            reference = centered[phase, shape["anchor"]]["power_by_signed_receiver_tilt"][sweep["best_signed_receiver_tilt_deg"]]
            sweep["matched_centered_A0_q_at_same_receiver_axis"] = reference
            sweep["q_gain_vs_matched_centered_A0_at_same_axis"] = sweep["power_at_best_axis"] - reference
            record["x_axis_receiver_sweep"] = sweep
        if shape["variant"] == "T0":
            checks["groove_invariant_quadrature"][f"{phase}/{sid}"] = groove_invariant_check(row)
        cases.setdefault(phase, {})[sid] = record

    for phase, anchor in product(PHASES, ANCHORS):
        a = indexed[phase, f"{anchor}_A0_dxplus"]
        b = indexed[phase, f"{anchor}_A0_dxminus"]
        fields_a = sample_fields(a, "A_total")
        fields_b = sample_fields(b, "A_total")
        # Reflect x: phi -> 180 - phi (mod 360), with invariant theta.
        n = len(fields_a[1])
        mirror_deltas = [abs(fields_a[0][0] - fields_b[0][0])]
        for ring_a, ring_b in zip(fields_a[1:], fields_b[1:]):
            mirror_deltas.extend(abs(value - ring_b[(n // 2 - i) % n]) for i, value in enumerate(ring_a))
        checks["signed_offset_mirror"][f"{phase}/{anchor}"] = {
            "max_sample_abs_error": max(mirror_deltas),
            "mean_sample_abs_error": math.fsum(mirror_deltas) / len(mirror_deltas),
        }

    for sid in cases["ideal"]:
        a = indexed["ideal", sid]
        b = indexed["floor90", sid]
        c = indexed["wall04", sid]
        checks["floor_change"][sid] = {
            "ideal_hemisphere_minus_floor_area_fraction": a[HEMISPHERE] - a["shape"]["floor_area_fraction"],
            "max_abs_floor90_minus_point9_ideal": max(
                abs(g - .9 * f) for ring_a, ring_b in zip(sample_fields(a, "A_total"), sample_fields(b, "A_total"))
                for f, g in zip(ring_a, ring_b)),
        }
        checks["wall_effect"][sid] = {
            "hemisphere_q_delta": c[HEMISPHERE] - b[HEMISPHERE],
            "normal_cone_q_delta": {str(alpha): c["normal_cones"][str(alpha)][POWER] - b["normal_cones"][str(alpha)][POWER]
                                    for alpha in RADII},
            "normal_cone_eta_delta": {str(alpha): c["normal_cones"][str(alpha)][FRACTION] - b["normal_cones"][str(alpha)][FRACTION]
                                      for alpha in RADII},
            "sidewall_and_land_hemisphere_q": sum(c["regional_power"][region][HEMISPHERE]
                                                  for region in ("A_sidewall", "A_land")),
        }

    return {
        "scope": "Same-anchor projected area, opening, floor area and depth matched. Specular geometric optics and assigned scalar absorptance, no measured TA2, film, wavelength or temperature.",
        "definitions": {
            "q": "Power relative to a blackbody of the same unit-cell footprint emitting into the full hemisphere.",
            "eta": "Cone share of the structure's own hemispherical output.",
            "x_axis_receiver_sweep": "Fixed 15° half-angle on interpolated x–z cut; positive axis points along +x (phi=0), negative along −x (phi=180). The 1° tilt increment is post-processing, not 1° ray sampling or a 2D search for an optimal receiver axis.",
            "q_rank": "Same-anchor comparisons at fixed cone and surface. Ties within ray/angular error are not resolved.",
            "groove_invariant_quadrature": "Independent hemisphere diagnostic for the y-invariant, x-symmetric T0 groove. Collapse saved directions to beta=atan2(sin(theta)*cos(phi),cos(theta)); average repeats in u=|sin(beta)| and integrate in u. Does not overwrite sampled-grid cone ranks or add solves.",
        },
        "sources": report["sources"],
        "cases": cases,
        "within_anchor_normal_cone_ranks": ranking(indexed),
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build_summary(json.loads(args.metrics.read_text()))
    write_json(args.output, result)
    print(json.dumps({"state": "COMPLETE_MATCHED_FAMILY_SUMMARY", "cases": sum(map(len, result["cases"].values())),
                      "output": str(args.output.resolve())}))


if __name__ == "__main__":
    main()
