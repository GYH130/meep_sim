"""Audited array adapters; this module never imports or initializes Meep.

Volume monitors must be created with yee_grid=True.  Each E/D pair retains its
own component grid, including the complete sampled y extent.  Plane monitors
use Meep's default centered grid with a zero-height integration region.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

from .common import canonical_period
from .geometry import periodic_slot_mask


def _grid_parameters(resolution, period, kx):
    resolution, period, kx = float(resolution), float(period), float(kx)
    if not all(math.isfinite(v) for v in (resolution, period, kx)):
        raise ValueError("ARRAY_NONFINITE_GRID_PARAMETER")
    if resolution <= 0 or period <= 0:
        raise ValueError("ARRAY_NONPOSITIVE_GRID_PARAMETER")
    return resolution, period, kx


def _field(sim, monitor, component, name):
    value = np.asarray(sim.get_dft_array(monitor, component, 0))
    if value.size == 0 or not np.isfinite(value).all():
        raise ValueError(f"ARRAY_EMPTY_OR_NONFINITE:{name}")
    return value


def _native_component(sim, monitor, component, name, resolution):
    value = _field(sim, monitor, component, name)
    dims, lower, upper = sim.get_array_slice_dimensions(component, vol=monitor.where)
    dims = np.asarray(dims)
    if dims.shape != (3,) or not np.isfinite(dims).all() or np.any(dims < 1) or np.any(dims != np.floor(dims)):
        raise ValueError(f"NATIVE_GRID_DIMENSIONS_INVALID:{name}")
    dims = dims.astype(int)
    if dims[2] != 1 or tuple(value.shape) != tuple(dims[:2]):
        raise ValueError(f"NATIVE_GRID_SHAPE_MISMATCH:{name}:{value.shape}:{tuple(dims)}")
    lo = np.asarray([lower.x, lower.y, lower.z], dtype=float)
    hi = np.asarray([upper.x, upper.y, upper.z], dtype=float)
    if not np.isfinite(lo).all() or not np.isfinite(hi).all():
        raise ValueError(f"NATIVE_GRID_NONFINITE_CORNER:{name}")
    tolerance = 1e-9 / max(resolution, 1.0)
    if abs(hi[2] - lo[2]) > tolerance:
        raise ValueError(f"NATIVE_GRID_NOT_TWO_DIMENSIONAL:{name}")
    coordinates = []
    for axis in (0, 1):
        coordinate = lo[axis] + np.arange(dims[axis], dtype=float) / resolution
        if not np.isclose(coordinate[-1], hi[axis], rtol=0, atol=tolerance):
            raise ValueError(f"NATIVE_GRID_CORNER_SPACING_MISMATCH:{name}:{axis}")
        coordinates.append(coordinate)
    return value, coordinates[0], coordinates[1], lo, hi


def native_volume_pairs(sim, monitor, pairs, resolution, period, kx, geometry, surface_y):
    """Return absorption pair records, complete raw arrays, and audit metadata.

    ``pairs`` maps a stable name to (E_component, D_component).  The simulation
    must already have initialized the DFT monitor.  Native component corners
    come from get_array_slice_dimensions, not centered get_array_metadata.
    Only redundant x-period guards are removed, through canonical_period's
    Bloch-copy verification; no y samples or shape mismatches are cropped.
    """
    resolution, period, kx = _grid_parameters(resolution, period, kx)
    if not isinstance(pairs, Mapping) or not pairs:
        raise ValueError("NATIVE_COMPONENT_PAIRS_REQUIRED")
    surface_y = float(surface_y)
    thickness = float(geometry.get("thickness_um", 60.0))
    if not math.isfinite(surface_y) or not math.isfinite(thickness) or thickness <= 0:
        raise ValueError("NATIVE_MATERIAL_REGION_INVALID")
    records, raw, component_metadata = [], {}, []
    for name, components in pairs.items():
        if not isinstance(name, str) or not name or len(components) != 2:
            raise ValueError("NATIVE_COMPONENT_PAIR_INVALID")
        e_component, d_component = components
        E, ex, ey, elo, ehi = _native_component(sim, monitor, e_component, name + ".E", resolution)
        D, dx, dy, dlo, dhi = _native_component(sim, monitor, d_component, name + ".D", resolution)
        if E.shape != D.shape or not np.array_equal(ex, dx) or not np.array_equal(ey, dy):
            raise ValueError(f"NATIVE_ED_GRID_NOT_COLLOCATED:{name}")
        x, selected, keep = canonical_period(ex, {"E": E, "D": D}, period, kx, resolution)
        X, Y = np.meshgrid(x, ey, indexing="ij")
        depth = surface_y - Y
        slab = (depth >= 0) & (depth <= thickness)
        if geometry.get("case") == "flat":
            ti_mask = slab
        else:
            ti_mask = slab & ~periodic_slot_mask(X, depth, geometry)
        weights = np.full(selected["E"].shape, resolution ** -2, dtype=float)
        records.append({"name": name, "E": selected["E"], "D": selected["D"],
                        "weights": weights, "ti_mask": ti_mask,
                        "coordinates": {"x": x, "y": ey}})
        prefix = name + "_"
        raw.update({prefix + "E_full": E, prefix + "D_full": D,
                    prefix + "x_full": ex, prefix + "y_full": ey,
                    prefix + "D_x_full": dx, prefix + "D_y_full": dy,
                    prefix + "E_lower_corner": elo, prefix + "E_upper_corner": ehi,
                    prefix + "D_lower_corner": dlo, prefix + "D_upper_corner": dhi,
                    prefix + "x_selected_indices": np.flatnonzero(keep),
                    prefix + "x_keep": keep, prefix + "x_canonical": x,
                    prefix + "y_canonical": ey, prefix + "E_canonical": selected["E"],
                    prefix + "D_canonical": selected["D"],
                    prefix + "weights_canonical": weights, prefix + "ti_mask_canonical": ti_mask})
        component_metadata.append({"name": name, "full_shape": list(E.shape),
                                   "canonical_shape": list(selected["E"].shape),
                                   "x_guards_removed": int((~keep).sum()),
                                   "x_selected_indices": np.flatnonzero(keep).tolist(),
                                   "x_min_um": float(x[0]), "x_max_um": float(x[-1]),
                                   "y_min_um": float(ey[0]), "y_max_um": float(ey[-1]),
                                   "minimum_depth_um": float(depth.min()),
                                   "maximum_depth_um": float(depth.max()),
                                   "weight_sum_um2": float(weights.sum()),
                                   "ti_weight_sum_um2": float(weights[ti_mask].sum()),
                                   "ed_collocated": True, "bloch_guards_verified": True})
    metadata = {"coordinate_source": "component-specific get_array_slice_dimensions",
                "grid": "native Yee; component-specific E/D collocation",
                "quadrature": "uniform (1/resolution)^2 at every retained native y sample",
                "centered_metadata_used": False, "y_samples_removed": 0,
                "period_um": period, "resolution": resolution, "kx": kx,
                "surface_y_um": surface_y, "thickness_um": thickness,
                "components": component_metadata}
    return records, raw, metadata


def plane_fields(sim, monitor, components, P, kx, res):
    """Read centered fields on a zero-height plane and audit x-period guards.

    Full centered cubature weights are preserved as evidence.  The returned
    canonical field samples have uniform spacing for the Floquet transform;
    their uniform period weights are also saved separately.
    """
    res, P, kx = _grid_parameters(res, P, kx)
    if not isinstance(components, Mapping) or not components:
        raise ValueError("PLANE_COMPONENTS_REQUIRED")
    x, y, z, weights = (np.asarray(value) for value in sim.get_array_metadata(dft_cell=monitor))
    if x.ndim != 1 or y.shape != (1,) or z.shape != (1,) or weights.shape != x.shape:
        raise ValueError("PLANE_METADATA_SHAPE_MISMATCH")
    if not all(np.isfinite(value).all() for value in (x, y, z, weights)):
        raise ValueError("PLANE_METADATA_NONFINITE")
    if np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError("PLANE_METADATA_WEIGHTS_INVALID")
    if not np.isclose(float(weights.sum()), P, rtol=1e-9, atol=1e-9):
        raise ValueError("PLANE_METADATA_WEIGHT_SUM_MISMATCH")
    arrays = {}
    for name, component in components.items():
        if not isinstance(name, str) or not name:
            raise ValueError("PLANE_COMPONENT_NAME_INVALID")
        value = _field(sim, monitor, component, name)
        if value.ndim != 1 or value.shape != x.shape:
            raise ValueError(f"PLANE_COMPONENT_SHAPE_MISMATCH:{name}")
        arrays[name] = value
    canonical_x, selected, keep = canonical_period(x, arrays, P, kx, res)
    raw = {"x_full": x, "y_full": y, "z_full": z, "weights_full": weights,
           "x_selected_indices": np.flatnonzero(keep), "x_keep": keep,
           "x_canonical": canonical_x, "weights_selected_from_metadata": weights[keep],
           "weights_canonical_uniform": np.full(canonical_x.shape, 1 / res)}
    for name, value in arrays.items():
        raw[name + "_full"] = value
        raw[name + "_canonical"] = selected[name]
    return canonical_x, selected, raw
