"""Run native and CLI geometry-family checks without third-party libraries."""

import argparse
import json
import math
from pathlib import Path
import subprocess


GEOMETRY = {
    "A0": ("hex_flat_to_flat", 6, True),
    "A3": ("hex_flat_to_flat", 6, True),
    "R0": ("polygon_vertex_diameter", 128, True),
    "S0": ("square_side", 4, False),
    "T0": ("trench_clear_x_width", 0, False),
}
PITCH, HEIGHT, BOTTOM, OPENING = 80.0, 82.0, 25.0, 66.0
RAYS = 2048
ROUND_AREA_FACTOR = 128 * math.sin(2 * math.pi / 128) / 8


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def close(actual, expected, tolerance, message):
    require(abs(actual - expected) <= tolerance, f"{message}: {actual} != {expected}")


def command(binary, family, surface, theta=38, phi=31, offset=0, floor=None, wall=None):
    argv = [
        str(binary),
        "--structure", family, "--surface", surface,
        "--theta", str(theta), "--phi", str(phi), "--offset-x-um", str(offset),
        "--pitch-um", str(PITCH), "--height-um", str(HEIGHT),
        "--bottom-width-um", str(BOTTOM), "--upper-width-um", str(OPENING),
        "--rays", str(RAYS), "--seed", "28021996", "--max-bounces", "256",
        "--segments", "128" if family == "A0" else "1",
    ]
    if floor is not None:
        argv += ["--floor-alpha", str(floor)]
    if wall is not None:
        argv += ["--wall-alpha", str(wall)]
    return argv


def run_json(argv):
    process = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=30)
    if process.returncode:
        raise RuntimeError(f"native check failed ({process.returncode}): {process.stderr.strip()}")
    return json.loads(process.stdout)


def geometry_areas(family):
    _, _, triangular = GEOMETRY[family]
    cell = PITCH**2 * (math.sqrt(3) / 2 if triangular else 1)
    if family == "T0":
        return cell, PITCH * BOTTOM, PITCH * OPENING
    factor = math.sqrt(3) / 2 if family in ("A0", "A3") else ROUND_AREA_FACTOR if family == "R0" else 1
    return cell, factor * BOTTOM**2, factor * OPENING**2


def check_geometry(row, family):
    definition, sides, _ = GEOMETRY[family]
    cell, floor, opening = geometry_areas(family)
    metadata = row["geometry_meta"]
    require(row["width_definition"] == definition, f"{family} width convention")
    require(metadata["opening_sides"] == sides, f"{family} aperture topology")
    require(metadata["y_invariant"] == (family == "T0"), f"{family} invariance")
    close(row["unit_cell_area_um2"], cell, 1e-8, f"{family} footprint")
    close(row["active_area_um2"], floor, 1e-8, f"{family} floor")
    close(row["opening_area_um2"], opening, 1e-8, f"{family} opening")
    close(metadata["floor_area_fraction"], floor / cell, 1e-11, f"{family} floor fraction")
    close(metadata["opening_area_fraction"], opening / cell, 1e-11, f"{family} opening fraction")


def check_energy(row, family, surface):
    require(row["surface"] == surface and row["structure"] == family, "native result binding")
    require(row["geometry_unresolved_rays"] == 0, f"{family} unresolved rays")
    require(row["max_ray_closure_error"] <= 1e-11, f"{family} per-ray closure")
    for channel in ("s", "p", "unpolarized"):
        result = row["channels"][channel]
        close(result["A_total"] + result["escape"] + result["residual"], 1, 1e-11, f"{family} closure")
        close(result["A_active"] + result["A_wall"], result["A_total"], 1e-11, f"{family} A partition")
        close(result["A_sidewall"] + result["A_land"], result["A_wall"], 1e-11, f"{family} wall partition")


def test_family(binary, family):
    ideal = run_json(command(binary, family, "M0"))
    check_geometry(ideal, family)
    check_energy(ideal, family, "M0")
    require((ideal["floor_alpha"], ideal["wall_alpha"]) == (1, 0), f"{family} M0 parameters")
    matched_grey = run_json(command(binary, family, "MG", floor=1, wall=0))
    grey = run_json(command(binary, family, "MG", floor=0.9, wall=0.04))
    check_energy(matched_grey, family, "MG")
    check_energy(grey, family, "MG")
    require((grey["floor_alpha"], grey["wall_alpha"]) == (0.9, 0.04), f"{family} MG parameters")
    require(grey["channels"]["unpolarized"]["A_land"] > 0, f"{family} grey land absorbs")
    for key in ("A_total", "A_active", "A_sidewall", "A_land", "escape", "residual"):
        close(ideal["channels"]["unpolarized"][key], matched_grey["channels"]["unpolarized"][key],
              1e-13, f"{family} M0/MG equality {key}")

    # Black floor and wall catch every ray on its first surface. At normal incidence
    # only the projected floor contributes to A_active. Finite spatial sampling
    # therefore needs a binomial tolerance, unlike deterministic area arithmetic.
    black = run_json(command(binary, family, "MG", theta=0, floor=1, wall=1))
    check_energy(black, family, "MG")
    cell, floor, _ = geometry_areas(family)
    expected = floor / cell
    tolerance = 6 * math.sqrt(expected * (1 - expected) / RAYS) + 1 / RAYS
    close(black["channels"]["unpolarized"]["A_active"], expected, tolerance,
          f"{family} direct normal floor catch")
    close(black["channels"]["unpolarized"]["A_total"], 1, 1e-13, f"{family} black-surface catch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    native = run_json([str(binary), "--self-test"])
    require(native["self_test"] == "PASS", "native self-test")
    for family in GEOMETRY:
        test_family(binary, family)
    for family in ("A0", "A3"):
        offset = run_json(command(binary, family, "MG", offset=-5, floor=0.9, wall=0.04))
        require(offset["offset_x_um"] == -5, f"{family} offset output")
        check_energy(offset, family, "MG")

    normal = run_json(command(binary, "T0", "MG", theta=0, phi=90, floor=0.9, wall=0.04))
    along = run_json(command(binary, "T0", "MG", theta=75, phi=90, floor=0.9, wall=0.04))
    check_energy(along, "T0", "MG")
    require(along["periodic_wraps"] > RAYS, "trench rays must cross y-normalization seams")
    for key in ("A_total", "A_active", "A_sidewall", "A_land", "escape"):
        close(along["channels"]["unpolarized"][key], normal["channels"]["unpolarized"][key],
              1e-11, f"trench along-y invariance {key}")
    print(json.dumps({
        "status": "PASS",
        "native_checks": native["checks"],
        "cli_families": list(GEOMETRY),
        "cli_checks": ["footprint_and_aperture_area", "m0_mg_endpoint_agreement", "grey_energy_partition",
                       "direct_normal_floor_catch", "offset_binding", "trench_y_invariance"],
        "rays_per_cli_fixture": RAYS,
    }, separators=(",", ":")))


if __name__ == "__main__":
    main()
