"""D14 left-apex triangular wedge groove with 2 um inner-wall film.

This diagnostics_v2 script reuses the D13 narrowband absorption workflow, but
replaces the wrapped parallelogram groove with a left-apex triangular wedge:

    P=50 um, top opening=40 um, depth=30 um, apex shift=30*tan(20 deg).

The 2 um measured-loss film remains on the groove sidewalls/apex shell using
the same outer-film plus inner-air representation as D13.  Outputs use an
independent D14 tag so existing D13 checkpoints are not overwritten.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

_D13_PATH = Path(__file__).with_name("D13_tilt20_wrapped_gap10_film2um_angle_screen.py")
_SPEC = importlib.util.spec_from_file_location("d13_angle_screen_base", _D13_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"Could not load D13 base script from {_D13_PATH}")
d13 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = d13
_SPEC.loader.exec_module(d13)
_D13_PREPARE_ARGS = d13._prepare_args


TAG = "D14_left_wedge_gap10_w40_h30_film2um"

d13.TAG = TAG
d13.DEFAULT_BOTTOM_WIDTH_UM = 0.0


def _left_wedge_vertices(
    *,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    y_surface: float,
) -> list[tuple[float, float]]:
    del bottom_width_um
    dx = depth_um * math.tan(math.radians(tilt_angle_deg))
    return [
        (-top_width_um / 2.0, y_surface),
        (+top_width_um / 2.0, y_surface),
        (-top_width_um / 2.0 - dx, y_surface - depth_um),
    ]


def _shifted_line_inward(
    p: tuple[float, float],
    q: tuple[float, float],
    distance_um: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    px, py = p
    qx, qy = q
    dx = qx - px
    dy = qy - py
    length = math.hypot(dx, dy)
    if length <= 0:
        raise ValueError("Cannot offset a zero-length wedge edge.")
    # The D14 vertices are clockwise.  The interior is on the right-hand side
    # of each directed edge, so this normal points into the air wedge.
    nx = dy / length
    ny = -dx / length
    return (px + distance_um * nx, py + distance_um * ny), (
        qx + distance_um * nx,
        qy + distance_um * ny,
    )


def _line_intersection(
    a0: tuple[float, float],
    a1: tuple[float, float],
    b0: tuple[float, float],
    b1: tuple[float, float],
) -> tuple[float, float]:
    x1, y1 = a0
    x2, y2 = a1
    x3, y3 = b0
    x4, y4 = b1
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(den) < 1e-12:
        raise ValueError("Offset wedge lines are parallel.")
    det_a = x1 * y2 - y1 * x2
    det_b = x3 * y4 - y3 * x4
    return (
        (det_a * (x3 - x4) - (x1 - x2) * det_b) / den,
        (det_a * (y3 - y4) - (y1 - y2) * det_b) / den,
    )


def _line_at_y(
    p: tuple[float, float],
    q: tuple[float, float],
    y_target: float,
) -> tuple[float, float]:
    px, py = p
    qx, qy = q
    if abs(qy - py) < 1e-12:
        raise ValueError("Cannot intersect a horizontal line with y=y_surface.")
    u = (y_target - py) / (qy - py)
    return px + u * (qx - px), y_target


def _point_line_distance(
    point: tuple[float, float],
    a: tuple[float, float],
    b: tuple[float, float],
) -> float:
    px, py = point
    ax, ay = a
    bx, by = b
    return abs((bx - ax) * (ay - py) - (ax - px) * (by - ay)) / math.hypot(bx - ax, by - ay)


def _left_wedge_inner_vertices(
    *,
    top_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    film_thickness_um: float,
    y_surface: float,
) -> list[tuple[float, float]]:
    outer = _left_wedge_vertices(
        top_width_um=top_width_um,
        bottom_width_um=0.0,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    left_top, right_top, apex = outer
    right_offset = _shifted_line_inward(right_top, apex, film_thickness_um)
    left_offset = _shifted_line_inward(apex, left_top, film_thickness_um)
    return [
        _line_at_y(*left_offset, y_surface),
        _line_at_y(*right_offset, y_surface),
        _line_intersection(*right_offset, *left_offset),
    ]


def _left_wedge_inner_wall_film_geometry(
    *,
    y_surface: float,
    period_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    film_thickness_um: float,
    substrate_thickness_um: float,
    medium_substrate,
    medium_film,
    medium_groove,
) -> list:
    import meep as mp

    del bottom_width_um
    if film_thickness_um <= 0:
        raise ValueError("D14 expects a positive inner-wall film thickness.")
    if 2.0 * film_thickness_um >= top_width_um:
        raise ValueError("Film thickness collapses the wedge air core.")
    if film_thickness_um >= depth_um:
        raise ValueError("Film thickness must be smaller than wedge depth.")
    if depth_um >= substrate_thickness_um:
        raise ValueError("Wedge depth must be smaller than substrate thickness.")

    outer = _left_wedge_vertices(
        top_width_um=top_width_um,
        bottom_width_um=0.0,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    inner = _left_wedge_inner_vertices(
        top_width_um=top_width_um,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        film_thickness_um=film_thickness_um,
        y_surface=y_surface,
    )
    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_um, substrate_thickness_um, mp.inf),
        )
    ]
    for piece in d13._wrapped_polygons(outer, period_um):
        geom.append(
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in piece],
                height=mp.inf,
                axis=mp.Vector3(0, 0, 1),
                material=medium_film,
            )
        )
    for piece in d13._wrapped_polygons(inner, period_um):
        geom.append(
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in piece],
                height=mp.inf,
                axis=mp.Vector3(0, 0, 1),
                material=medium_groove,
            )
        )
    return geom


def _geometry_factory(film_medium, args):
    import meep as mp

    ti = d13.get_ti_medium()
    air = mp.Medium(epsilon=1.0)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        return _left_wedge_inner_wall_film_geometry(
            y_surface=y_surface_um,
            period_um=args.period_um,
            top_width_um=args.top_width_um,
            bottom_width_um=0.0,
            depth_um=args.depth_um,
            tilt_angle_deg=args.tilt_angle_deg,
            film_thickness_um=args.film_thickness_um,
            substrate_thickness_um=substrate_thickness_um,
            medium_substrate=ti,
            medium_film=film_medium,
            medium_groove=air,
        )

    return factory


def _geometry_checks(args):
    half = args.period_um / 2.0
    dx = args.depth_um * math.tan(math.radians(args.tilt_angle_deg))
    apex_x = -args.top_width_um / 2.0 - dx
    outer = _left_wedge_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=0.0,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=0.0,
    )
    inner = _left_wedge_inner_vertices(
        top_width_um=args.top_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        film_thickness_um=args.film_thickness_um,
        y_surface=0.0,
    )
    inner_top = inner[1][0] - inner[0][0]
    inner_depth = -inner[2][1]
    left_distance = _point_line_distance(inner[0], outer[2], outer[0])
    right_distance = _point_line_distance(inner[1], outer[1], outer[2])
    apex_left_distance = _point_line_distance(inner[2], outer[2], outer[0])
    apex_right_distance = _point_line_distance(inner[2], outer[1], outer[2])
    return d13.pd.DataFrame(
        [
            {
                "check_name": "top_gap10_preserved",
                "status": "PASS" if abs(args.period_um - args.top_width_um - 10.0) < 1e-9 else "FAIL",
                "details": f"period={args.period_um:g}, top_width={args.top_width_um:g}, gap={args.period_um - args.top_width_um:g}",
            },
            {
                "check_name": "left_apex_wraps_periodically",
                "status": "PASS" if apex_x < -half else "FAIL",
                "details": f"apex_x={apex_x:.6g} um, half_period={half:g} um, shift={dx:.6g} um",
            },
            {
                "check_name": "film_shell_positive",
                "status": "PASS" if inner_top > 0 and inner_depth > 0 else "FAIL",
                "details": f"inner_top={inner_top:g} um, inner_depth={inner_depth:g} um",
            },
            {
                "check_name": "wall_normal_film_offset",
                "status": "PASS"
                if max(
                    abs(left_distance - args.film_thickness_um),
                    abs(right_distance - args.film_thickness_um),
                    abs(apex_left_distance - args.film_thickness_um),
                    abs(apex_right_distance - args.film_thickness_um),
                )
                < 1e-8
                else "FAIL",
                "details": (
                    f"left_top={left_distance:.6g}, right_top={right_distance:.6g}, "
                    f"apex_left={apex_left_distance:.6g}, apex_right={apex_right_distance:.6g} um"
                ),
            },
            {
                "check_name": "depth_below_substrate",
                "status": "PASS" if args.depth_um < args.substrate_thickness_um else "FAIL",
                "details": f"depth={args.depth_um:g} um, substrate={args.substrate_thickness_um:g} um",
            },
            {
                "check_name": "wrapped_polygon_pieces_exist",
                "status": "PASS"
                if d13._wrapped_polygons(outer, args.period_um) and d13._wrapped_polygons(inner, args.period_um)
                else "FAIL",
                "details": (
                    f"outer_pieces={len(d13._wrapped_polygons(outer, args.period_um))}, "
                    f"inner_pieces={len(d13._wrapped_polygons(inner, args.period_um))}"
                ),
            },
        ]
    )


def _add_wedge_patches(ax, outer, inner, *, show_outer_outline: bool = True) -> None:
    ax.add_patch(d13.MplPolygon(outer, closed=True, facecolor="#F39C34", edgecolor="#1f77b4", lw=1.2))
    ax.add_patch(d13.MplPolygon(inner, closed=True, facecolor="white", edgecolor="#1f77b4", lw=1.2))
    if show_outer_outline:
        ax.add_patch(d13.MplPolygon(outer, closed=True, facecolor="none", edgecolor="black", lw=1.2))


def _plot_continuous_geometry(path: Path, args, outer, inner) -> None:
    xs = [pt[0] for pt in outer + inner]
    fig, ax = d13.plt.subplots(figsize=(7.0, 4.8))
    pad = 3.0
    x_min = min(xs) - pad
    x_max = max(xs) + pad
    ax.add_patch(
        d13.plt.Rectangle(
            (x_min, -args.substrate_thickness_um),
            x_max - x_min,
            args.substrate_thickness_um,
            color="#B8B8B8",
        )
    )
    _add_wedge_patches(ax, outer, inner)
    ax.axhline(0, color="#444444", lw=1.2)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-args.depth_um - 5, 5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("D14 one complete continuous left-apex wedge, 2 um normal-offset film")
    d13.save_figure(fig, path)
    d13.plt.close(fig)


def _plot_geometry(paths: dict[str, Path], args) -> None:
    outer = _left_wedge_vertices(
        top_width_um=args.top_width_um,
        bottom_width_um=0.0,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        y_surface=0.0,
    )
    inner = _left_wedge_inner_vertices(
        top_width_um=args.top_width_um,
        depth_um=args.depth_um,
        tilt_angle_deg=args.tilt_angle_deg,
        film_thickness_um=args.film_thickness_um,
        y_surface=0.0,
    )
    outer_pieces = d13._wrapped_polygons(outer, args.period_um)
    inner_pieces = d13._wrapped_polygons(inner, args.period_um)
    half = args.period_um / 2.0
    continuous_path = paths["geometry"].with_name(paths["geometry"].stem + "_continuous.png")

    fig, ax = d13.plt.subplots(figsize=(6.4, 5.0))
    ax.add_patch(
        d13.plt.Rectangle(
            (-half, -args.substrate_thickness_um),
            args.period_um,
            args.substrate_thickness_um,
            color="#B8B8B8",
        )
    )
    for piece in outer_pieces:
        ax.add_patch(d13.MplPolygon(piece, closed=True, facecolor="#F39C34", edgecolor="#1f77b4", lw=1.2))
    for piece in inner_pieces:
        ax.add_patch(d13.MplPolygon(piece, closed=True, facecolor="white", edgecolor="#1f77b4", lw=1.2))
    for piece in outer_pieces:
        ax.add_patch(d13.MplPolygon(piece, closed=True, facecolor="none", edgecolor="black", lw=1.2))

    ax.axhline(0, color="#444444", lw=1.2)
    ax.axvline(-half, color="#555555", lw=1.0)
    ax.axvline(half, color="#555555", lw=1.0)
    ax.set_xlim(-half - 3, half + 3)
    ax.set_ylim(-args.depth_um - 5, 5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("D14 left-apex triangular wedge, gap10, 2 um film")
    d13.save_figure(fig, paths["geometry"])
    d13.plt.close(fig)
    _plot_continuous_geometry(continuous_path, args, outer, inner)


def _write_report(paths, checks, normal_metrics, angle_metrics, args) -> None:
    dx = args.depth_um * math.tan(math.radians(args.tilt_angle_deg))
    apex_x = -args.top_width_um / 2.0 - dx
    continuous_path = paths["geometry"].with_name(paths["geometry"].stem + "_continuous.png")
    check_lines = ["| check_name | status | details |", "| --- | --- | --- |"]
    for row in checks.itertuples(index=False):
        check_lines.append(f"| {row.check_name} | {row.status} | {row.details} |")
    report = f"""# D14 Left-Apex Triangular Wedge Angle Screen

## Geometry

- Period: {args.period_um:g} um.
- Top wedge opening: {args.top_width_um:g} um.
- Top Ti gap: {args.period_um - args.top_width_um:g} um.
- Depth: {args.depth_um:g} um.
- Left apex: x={apex_x:.6g} um, y={-args.depth_um:g} um.
- Apex shift magnitude: {dx:.6g} um from tan({args.tilt_angle_deg:g} deg).
- Film: {args.film_thickness_um:g} um normal-offset shell on sidewalls/apex.
- Representation: clipped/wrapped periodic triangular wedge pieces.

## Geometry Checks

{chr(10).join(check_lines)}

## Outputs

- Normal spectra: `{paths['normal_spectra']}`
- Normal metrics: `{paths['normal_metrics']}`
- Angle spectra: `{paths['angle_spectra']}`
- Angle metrics: `{paths['angle_metrics']}`
- Geometry figure: `{paths['geometry']}`
- Continuous wedge geometry figure: `{continuous_path}`
- Log: `{paths['log']}`
"""
    paths["report"].write_text(report, encoding="utf-8")


def _prepare_args(args):
    args = _D13_PREPARE_ARGS(args)
    args.bottom_width_um = 0.0
    return args


d13._slanted_vertices = _left_wedge_vertices
d13._wrapped_inner_wall_film_geometry = _left_wedge_inner_wall_film_geometry
d13._geometry_factory = _geometry_factory
d13._geometry_checks = _geometry_checks
d13._plot_geometry = _plot_geometry
d13._write_report = _write_report
d13._prepare_args = _prepare_args


if __name__ == "__main__":
    raise SystemExit(d13.main())
