"""Validate both coarse result sets and select twelve paired M0 refinements."""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


ANGLES = (15, 30)
DIMENSIONS = ("height_um", "bottom_width_um", "upper_width_um", "pitch_um")
SHAPE_KEYS = ("id", "family", *DIMENSIONS, "curve_segments")
LANES = ("15_strength", "30_strength", "15_concentration", "30_concentration",
         "15_balance", "30_balance")
SAMPLE_COLUMNS = ["phi_deg", "A_total", "residual", "se_spatial_only"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    data = path.read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def shape_identity(shape):
    return tuple(shape[key] for key in SHAPE_KEYS)


def phase_grid(phase):
    theta = sorted(phase["theta_deg"])
    phi = sorted(phase["phi_deg"])
    require(theta[0] == 0 and theta[-1] == 89 and len(set(theta)) == len(theta),
            "Expected polar grid from 0 to 89 degrees")
    require(phi == list(range(0, 360, 15)), "Expected full 15-degree azimuth grid")
    return theta, phi, 1 + (len(theta) - 1) * len(phi)


def integrate(theta, samples, limit):
    """Integrate azimuth means in normalized projected solid angle sin²(theta)."""
    x = [math.sin(math.radians(angle)) ** 2 for angle in theta] + [1.0]
    samples = [*samples, samples[-1]]
    stop = math.sin(math.radians(limit)) ** 2
    terms = []
    for x0, x1, y0, y1 in zip(x, x[1:], samples, samples[1:]):
        if x0 >= stop:
            break
        dx = min(x1, stop) - x0
        terms.append(dx * (y0 + dx * (y1 - y0) / (2 * (x1 - x0))))
    return math.fsum(terms)


def score_shape(row, expected_shape, theta, phi, angle_count):
    shape = row["shape"]
    sid = shape["id"]
    require(shape_identity(shape) == shape_identity(expected_shape), f"{sid}: geometry mismatch")
    require(row["completed_angles"] == row["expected_angles"] == angle_count,
            f"{sid}: incomplete direction count")
    rings = sorted(row["angles"], key=lambda ring: ring["theta_deg"])
    require([ring["theta_deg"] for ring in rings] == theta, f"{sid}: polar grid mismatch")
    means, uppers, residuals = [], [], []
    for ring in rings:
        expected_phi = phi if ring["theta_deg"] else [0]
        samples = ring["azimuth_samples"]
        require(ring["azimuth_count"] == len(samples) == len(expected_phi)
                and [sample[0] for sample in samples] == expected_phi,
                f"{sid}: azimuth grid mismatch at {ring['theta_deg']} degrees")
        require(all(len(sample) == 4 and all(math.isfinite(value) for value in sample)
                    and -1e-10 <= sample[1] <= 1 + 1e-10 and 0 <= sample[2] <= 1
                    and sample[1] + sample[2] <= 1 + 1e-10 and sample[3] >= 0 for sample in samples),
                f"{sid}: invalid absorptance, residual or sampling error")
        mean = math.fsum(sample[1] for sample in samples) / len(samples)
        upper = math.fsum(sample[1] + sample[2] for sample in samples) / len(samples)
        residual = max(sample[2] for sample in samples)
        require(math.isclose(mean, ring["mean_A"], abs_tol=1e-10, rel_tol=0)
                and math.isclose(upper, ring["mean_A_upper_bound"], abs_tol=1e-10, rel_tol=0)
                and math.isclose(residual, ring["max_residual"], abs_tol=1e-10, rel_tol=0),
                f"{sid}: azimuth summary disagrees with samples")
        means.append(max(0.0, mean))
        uppers.append(max(0.0, upper))
        residuals.append(residual)

    total = integrate(theta, means, 90)
    total_upper = integrate(theta, uppers, 90)
    require(total > 0, f"{sid}: zero hemisphere power cannot be ranked")
    cones = {}
    for angle in ANGLES:
        weight = math.sin(math.radians(angle)) ** 2
        power = integrate(theta, means, angle)
        power_upper = integrate(theta, uppers, angle)
        out = max(0.0, total - power)
        out_upper = max(out, total_upper - power_upper)
        cones[str(angle)] = {
            "mean_ideal_emissivity_inside": power / weight,
            "power_fraction": min(1.0, power / total),
            "cell_normalized_power": power,
            "mean_ideal_emissivity_interval": [power / weight, power_upper / weight],
            "power_fraction_interval": [
                power / (power + out_upper) if power + out_upper else 0.0,
                power_upper / (power_upper + out) if power_upper + out else 0.0,
            ],
        }
    h, b, d, p = (shape[key] for key in DIMENSIONS)
    return {
        "shape": expected_shape,
        "cones": cones,
        "hemispherical_ideal_emissivity": total,
        "hemispherical_ideal_emissivity_upper_bound": total_upper,
        "floor_area_fraction_of_cell": (b / p) ** 2,
        "hemi_minus_floor_area_fraction": total - (b / p) ** 2,
        "land_width_along_neighbor_centers_um": p - d,
        "opening_to_pitch_ratio": d / p,
        "height_to_opening_ratio": h / d,
        "bottom_to_opening_ratio": b / d,
        "max_residual": max(residuals),
    }


def score_phase(summary, manifest, manifest_sha256, phase_name):
    phase = manifest["phases"][phase_name]
    theta, phi, angles = phase_grid(phase)
    shape_ids = phase["shape_ids"]
    shapes = {shape["id"]: shape for shape in manifest["shapes"]}
    require(len(shapes) == len(manifest["shapes"]) and len(set(shape_ids)) == len(shape_ids)
            and set(shape_ids) <= shapes.keys(), "Invalid phase geometry inventory")
    require(summary["kind"] == manifest["kind"] and summary["phase"] == phase_name
            and summary["state"] == "COMPLETE_IDEAL_PILOT"
            and not summary["failed"] and summary["not_dispatched"] == 0,
            f"{phase_name}: incomplete or failed run")
    require(summary["manifest_sha256"] == manifest_sha256
            and summary["engine_source_sha256"] == manifest["engine_source_sha256"]
            and summary["azimuth_sample_columns"] == SAMPLE_COLUMNS,
            f"{phase_name}: source or sample column mismatch")
    require(summary["condition_count"] == summary["completed_conditions"] == len(shape_ids) * angles
            and summary["completed_geometries"] == summary["planned_geometry_count"] == len(shape_ids)
            and summary["angles_per_shape"] == angles
            and summary["rays_per_condition"] == phase["rays_per_condition"],
            f"{phase_name}: solver counts or ray settings mismatch")
    rows = summary["geometry_summaries"]
    require(len(rows) == len(shape_ids) and {row["shape"]["id"] for row in rows} == set(shape_ids),
            f"{phase_name}: result geometry inventory mismatch")
    return {
        row["shape"]["id"]: score_shape(row, shapes[row["shape"]["id"]], theta, phi, angles)
        for row in rows
    }


def matched_pairs(metrics):
    groups = {}
    for sid, record in metrics.items():
        match = re.fullmatch(r"(A[03])_(G\d+)", sid)
        require(match is not None and match[1] == record["shape"]["family"],
                f"{sid}: unrecognized matched-pair ID")
        groups.setdefault(match[2], {})[match[1]] = record
    dimensions = set()
    for gid, group in groups.items():
        require(set(group) == {"A0", "A3"}, f"{gid}: missing matched wall profile")
        curved, straight = (
            tuple(group[family]["shape"][key] for key in DIMENSIONS) for family in ("A0", "A3")
        )
        require(curved == straight and curved not in dimensions,
                f"{gid}: mismatched or duplicate dimensions")
        dimensions.add(curved)
    return groups


def lane_scores(group, lane):
    angle, metric = lane.split("_")
    scores = []
    for family, row in group.items():
        cone = row["cones"][angle]
        strength = cone["mean_ideal_emissivity_inside"]
        fraction = cone["power_fraction"]
        primary = {"strength": strength, "concentration": fraction,
                   "balance": strength * fraction}[metric]
        # Fraction can saturate at one. Use actual in-cone output as its tie break.
        scores.append((round(primary, 12), strength, fraction, family))
    return max(scores)


def shortlist(pairs, reused_pairs):
    ranks = {
        lane: sorted(pairs, key=lambda gid: (
            *(-value for value in lane_scores(pairs[gid], lane)[:3]), gid))
        for lane in LANES
    }
    positions = {lane: {gid: i for i, gid in enumerate(order, 1)} for lane, order in ranks.items()}
    selected = {}

    def include(gid, lane):
        score, strength, fraction, family = lane_scores(pairs[gid], lane)
        selected[gid] = {
            "pair_id": gid, "reason": lane, "driving_family": family, "lane_score": score,
            "driving_strength": strength, "driving_fraction": fraction,
            "coarse_ranks_among_100_pairs": {name: order[gid] for name, order in positions.items()},
        }

    # Retain the previous near-leader only when it still ranks in the top 12 of a strength lane.
    if "G05" in pairs and "G05" not in reused_pairs:
        strength_lane = min(("15_strength", "30_strength"), key=lambda lane: positions[lane]["G05"])
        if positions[strength_lane]["G05"] <= 12:
            include("G05", strength_lane)

    # Round robin across both target cones and strength, concentration and balance.
    while len(selected) < 12:
        for lane in LANES:
            gid = next(gid for gid in ranks[lane] if gid not in selected and gid not in reused_pairs)
            include(gid, lane)
            if len(selected) == 12:
                break
    return list(selected.values()), ranks


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "previous-coarse", "previous-refine", "coarse",
                 "output-manifest", "output-selection"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    prior = args.previous_coarse.parent
    require(args.previous_refine.parent == prior, "Prior files must share their pinned receipt")
    receipt, receipt_sha = read_json(prior / "receipt.json")
    previous, previous_sha = read_json(prior / "manifest.json")
    old_coarse, old_coarse_sha = read_json(args.previous_coarse)
    old_refine, old_refine_sha = read_json(args.previous_refine)
    manifest, manifest_sha = read_json(args.manifest)
    coarse, coarse_sha = read_json(args.coarse)
    for name, sha in (("manifest.json", previous_sha), (args.previous_coarse.name, old_coarse_sha),
                      (args.previous_refine.name, old_refine_sha)):
        require(receipt["files"][name] == sha, f"Pinned prior hash differs: {name}")
    require(previous["schema"] == manifest["schema"] == 2
            and manifest["source_geometry_manifest_sha256"] == previous_sha
            and all(previous[key] == manifest[key] for key in
                    ("kind", "engine_source", "engine_source_sha256", "seed_base", "max_bounces")),
            "Previous and new manifests must use the same engine and ray settings")
    new_phase, old_phase = manifest["phases"]["coarse"], previous["phases"]["coarse"]
    require(phase_grid(new_phase) == phase_grid(old_phase)
            and new_phase["rays_per_condition"] == old_phase["rays_per_condition"],
            "Previous and new coarse angular sampling must match")

    old_metrics = score_phase(old_coarse, previous, previous_sha, "coarse")
    old_refined = score_phase(old_refine, previous, previous_sha, "refine")
    new_metrics = score_phase(coarse, manifest, manifest_sha, "coarse")
    shapes = {shape["id"]: shape for shape in manifest["shapes"]}
    require(len(old_metrics) == 20 and len(new_metrics) == 180 and len(shapes) == 200
            and not old_metrics.keys() & new_metrics.keys()
            and shapes.keys() == old_metrics.keys() | new_metrics.keys(),
            "Expected twenty reused and 180 new coarse geometries")
    require(all(shape_identity(row["shape"]) == shape_identity(shapes[sid])
                for sid, row in old_metrics.items()), "Reused dimensions differ in the new manifest")
    require(len(old_refined) == 6 and old_refined.keys() <= old_metrics.keys(),
            "Expected six previous fine geometries with matching coarse inputs")
    pairs = matched_pairs({**old_metrics, **new_metrics})
    reused_pairs = matched_pairs(old_refined)
    selections, ranks = shortlist(pairs, reused_pairs)
    pair_ids = [row["pair_id"] for row in selections]
    shape_ids = [f"{family}_{gid}" for gid in pair_ids for family in ("A0", "A3")]

    phase = previous["phases"]["refine"].copy()
    angles = phase_grid(phase)[2]
    phase.update(shape_ids=shape_ids, shape_count=len(shape_ids), angles_per_shape=angles,
                 condition_count=len(shape_ids) * angles)
    refine_manifest = {
        key: manifest[key] for key in
        ("schema", "kind", "scope", "engine_source", "engine_source_sha256", "seed_base", "max_bounces")
    }
    refine_manifest.update(
        source_geometry_manifest_sha256=manifest_sha, shapes=[shapes[sid] for sid in shape_ids],
        shape_count=len(shape_ids), condition_count=phase["condition_count"], phases={"refine": phase},
    )
    manifest_data = json_bytes(refine_manifest)
    selection = {
        "kind": "MATCHED_M0_REFINEMENT_SELECTION",
        "source_coarse_shape_count": len(shapes),
        "new_coarse_condition_count": coarse["completed_conditions"],
        "reused_coarse_condition_count": old_coarse["completed_conditions"],
        "reused_refine_shape_ids": sorted(old_refined),
        "selected_pair_ids": pair_ids, "shape_ids": shape_ids,
        "selected_shape_count": len(shape_ids), "angles_per_shape": angles,
        "rays_per_condition": phase["rays_per_condition"], "condition_count": phase["condition_count"],
        "refine_manifest_sha256": hashlib.sha256(manifest_data).hexdigest(),
        "inputs": {
            "coarse_manifest_sha256": manifest_sha, "coarse_summary_sha256": coarse_sha,
            "previous_manifest_sha256": previous_sha, "previous_coarse_summary_sha256": old_coarse_sha,
            "previous_refine_summary_sha256": old_refine_sha, "previous_receipt_sha256": receipt_sha,
        },
        "method": {
            "quadrature": "Azimuth mean; integrate in sin²(theta), linearly between samples, "
                          "hold the 89-degree mean through 90 degrees.",
            "lanes": "Select matched pairs round robin across 15/30-degree cone strength, fraction, "
                     "and their product. Each lane uses the better single wall profile of each pair.",
            "g05": "Include G05 as a comparison point if it remains top 12 in either cone strength lane.",
            "limits": "M0 ideal black floor and mirror walls. Engineering ratios are descriptors, "
                      "not fabrication pass criteria. Residual intervals exclude angular and spatial sampling error.",
        },
        "selections": selections, "pair_ranks": ranks,
        "metrics": {**old_metrics, **new_metrics}, "previous_refined_metrics": old_refined,
    }
    require(not args.output_manifest.exists() and not args.output_selection.exists(),
            "Refinement manifest and selection outputs are create-only")
    with args.output_manifest.open("xb") as output:
        output.write(manifest_data)
    with args.output_selection.open("xb") as output:
        output.write(json_bytes(selection))
    print(json.dumps({key: selection[key] for key in
                      ("selected_pair_ids", "selected_shape_count", "condition_count", "refine_manifest_sha256")}))


if __name__ == "__main__":
    main()
