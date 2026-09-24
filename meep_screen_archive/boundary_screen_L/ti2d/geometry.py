"""Pure geometry for a periodic, finite, tilted blind slot (all lengths in µm).

The top is physical z=0 and ``depth=-z``.  In (x, depth), the axis is
``(sin(alpha), cos(alpha))``.  A cavity is the intersection of
``abs(x*cos(alpha)-depth*sin(alpha)) <= W/2``,
``x*sin(alpha)+depth*cos(alpha) <= L``, and ``depth >= 0``.
There is deliberately no axis-coordinate >= 0 cap at the entrance.
"""

from __future__ import annotations

import math
from itertools import combinations

import numpy as np


def _polygon(width, length, sine, cosine):
    """Intersect the four closed half-planes independently of vertex formulas."""
    normals = np.array([[cosine, -sine], [-cosine, sine],
                        [sine, cosine], [0.0, -1.0]])
    bounds = np.array([width / 2, width / 2, length, 0.0])
    tolerance = 1e-11 * max(1.0, width, length)
    vertices = []
    for i, j in combinations(range(4), 2):
        pair = normals[[i, j]]
        if abs(np.linalg.det(pair)) < 1e-14:
            continue
        point = np.linalg.solve(pair, bounds[[i, j]])
        if np.all(normals @ point <= bounds + tolerance):
            if not any(np.linalg.norm(point - old) <= tolerance for old in vertices):
                point[np.abs(point) < 1e-14] = 0.0
                vertices.append(point)
    if len(vertices) < 3:
        raise ValueError("GEOMETRY_DEGENERATE: fewer than three vertices")
    points = np.asarray(vertices)
    center = points.mean(axis=0)
    angle = np.arctan2(points[:, 1] - center[1], points[:, 0] - center[0])
    return points[np.argsort(angle)]


def _point_segment_distance(point, start, end):
    delta = end - start
    denominator = float(delta @ delta)
    if denominator == 0:
        return float(np.linalg.norm(point - start))
    fraction = np.clip(float((point - start) @ delta) / denominator, 0.0, 1.0)
    return float(np.linalg.norm(point - start - fraction * delta))


def polygon_distance(first, second):
    """Exact Euclidean distance of two convex closed polygons, including overlap.

    The separating-axis test includes containment and collinear edge contact.
    If disjoint, a closest pair contains a vertex and a finite edge.
    """
    first, second = np.asarray(first, float), np.asarray(second, float)
    separated = False
    for polygon in (first, second):
        for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
            edge = end - start
            normal = np.array([-edge[1], edge[0]])
            a, b = first @ normal, second @ normal
            if a.max() < b.min() or b.max() < a.min():
                separated = True
                break
        if separated:
            break
    if not separated:
        return 0.0
    return min(
        _point_segment_distance(point, start, end)
        for points, edges in ((first, second), (second, first))
        for point in points
        for start, end in zip(edges, np.roll(edges, -1, axis=0))
    )


def validate_geometry(g):
    """Return a JSON-safe report; invalid inputs never qualify a simulation.

    ``periodic_shifts`` are integer cell indices k, so a copied vertex is
    ``(x + k*period_um, depth)``.  They include every cavity whose *finite*
    polygon intersects the canonical x cell [-P/2, P/2].  ``min_wall_um`` and
    ``residual_um`` report achieved distances, not requested thresholds.
    The residual gate is strict: thickness > maximum_depth + residual_required.
    """
    report = {"valid": False, "status": "GEOMETRY_INVALID", "errors": [],
              "vertices_x_depth": [], "periodic_shifts": [],
              "width_normal_um": None, "surface_opening_um": None,
              "vertical_center_depth_um": None, "maximum_depth_um": None,
              "min_wall_um": None, "residual_um": None}
    try:
        values = {key: float(g[key]) for key in
                  ("period_um", "surface_fill", "axis_length_um", "tilt_deg")}
        values.update({key: float(g.get(key, default)) for key, default in
                       (("thickness_um", 60.0), ("min_wall_um", 1.0),
                        ("residual_um", 10.0))})
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("GEOMETRY_NONFINITE")
        period, fill, length, tilt = (values[key] for key in
                                     ("period_um", "surface_fill", "axis_length_um", "tilt_deg"))
        thickness = values["thickness_um"]
        required_wall, required_residual = values["min_wall_um"], values["residual_um"]
        if period <= 0 or length <= 0 or thickness <= 0:
            raise ValueError("GEOMETRY_NONPOSITIVE_LENGTH")
        if not 0 < fill < 1:
            raise ValueError("GEOMETRY_SURFACE_FILL: require 0 < F < 1")
        if not abs(tilt) < 90:
            raise ValueError("GEOMETRY_TILT: require abs(tilt_deg) < 90")
        if required_wall < 1.0 or required_residual < 10.0:
            raise ValueError("GEOMETRY_REQUIREMENTS: wall >= 1 µm and residual >= 10 µm")
        angle = math.radians(tilt)
        sine, cosine = math.sin(angle), math.cos(angle)
        width = period * fill * cosine
        bottom_top = length * cosine - width * abs(sine) / 2
        if bottom_top <= 0:
            raise ValueError("GEOMETRY_SURFACE_BOTTOM_REGIME: bottom intersects surface")
        vertices = _polygon(width, length, sine, cosine)
        maximum_depth = float(vertices[:, 1].max())
        span = float(np.ptp(vertices[:, 0]))
        wall = math.inf
        checked = []
        k = 1
        # Any omitted neighbor has horizontal separation k*P-span, and strip
        # separation k*P*cos(alpha)-W.  Their maximum is a rigorous lower
        # bound, so a bounding box can skip work but never reject geometry.
        while True:
            lower_bound = max(0.0, k * period - span, k * period * cosine - width)
            if lower_bound > wall:
                break
            for sign in (-1, 1):
                shifted = vertices + np.array([sign * k * period, 0.0])
                wall = min(wall, polygon_distance(vertices, shifted))
                checked.append(sign * k)
            if wall == 0:
                break
            k += 1
        x_min, x_max = float(vertices[:, 0].min()), float(vertices[:, 0].max())
        first = math.ceil((-period / 2 - x_max) / period - 1e-12)
        last = math.floor((period / 2 - x_min) / period + 1e-12)
        residual = thickness - maximum_depth
        report.update(values)
        report.update(vertices_x_depth=vertices.tolist(), periodic_shifts=list(range(first, last + 1)),
                      width_normal_um=width, surface_opening_um=width / cosine,
                      vertical_center_depth_um=length * cosine, maximum_depth_um=maximum_depth,
                      bottom_minimum_depth_um=bottom_top,
                      min_wall_um=wall, residual_um=residual,
                      required_min_wall_um=required_wall, required_residual_um=required_residual,
                      wall_lower_bound_um=max(0.0, period * cosine - width),
                      neighbor_shifts_checked=checked, x_span_um=span,
                      surface_bottom_regime_valid=True)
        errors = report["errors"]
        if wall < required_wall:
            errors.append("GEOMETRY_WALL_TOO_THIN")
        if not thickness > maximum_depth + required_residual:
            errors.append("GEOMETRY_RESIDUAL_TOO_THIN")
        report["valid"] = not errors
        report["status"] = "VALID" if not errors else errors[0]
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        report["errors"].append(str(exc))
        report["status"] = str(exc).split(":", 1)[0]
    return report


def periodic_slot_mask(x, d, g):
    """Membership in the union of finite periodic cavities; arrays broadcast.

    Canonicalizing x alone is insufficient for a tilted finite cavity.  After
    canonicalization, all finite neighboring polygons crossing that cell are
    tested with their own axis and bottom planes.
    """
    report = validate_geometry(g)
    if not report["valid"]:
        raise ValueError(report["status"])
    x, d = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(d, dtype=float))
    if not np.isfinite(x).all() or not np.isfinite(d).all():
        raise ValueError("GEOMETRY_NONFINITE_COORDINATE")
    period = float(g["period_um"])
    wrapped = (x + period / 2) % period - period / 2
    angle = math.radians(float(g["tilt_deg"]))
    sine, cosine = math.sin(angle), math.cos(angle)
    width, length = report["width_normal_um"], float(g["axis_length_um"])
    mask = np.zeros(x.shape, dtype=bool)
    tolerance = 1e-12 * max(1.0, period, length)
    for k in report["periodic_shifts"]:
        local_x = wrapped - k * period
        mask |= ((np.abs(local_x * cosine - d * sine) <= width / 2 + tolerance)
                 & (local_x * sine + d * cosine <= length + tolerance)
                 & (d >= 0))
    return mask
