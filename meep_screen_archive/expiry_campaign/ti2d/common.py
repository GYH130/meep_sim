"""Small atomic, hashed artifacts; no network access and no legacy writes."""
from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return {'real': float(value.real), 'imag': float(value.imag)}
    if isinstance(value, np.generic):
        return value.item()
    return value


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def digest_object(value):
    return hashlib.sha256(json.dumps(jsonable(value), sort_keys=True, allow_nan=False,
                                     separators=(',', ':')).encode()).hexdigest()


def implementation_digest():
    root = Path(__file__).parent
    return digest_object({p.name: digest_file(p) for p in sorted(root.glob('*.py'))})


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(jsonable(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n'
    fd, temp = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as out:
            out.write(payload)
            out.flush()
            os.fsync(out.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def atomic_npz(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as out:
            np.savez_compressed(out, **arrays)
            out.flush()
            os.fsync(out.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def canonical_period(x, arrays, period, kx, resolution):
    """Explicitly select one canonical period and audit redundant Bloch samples.

    Meep centered/native arrays can include one or two interpolation guards.
    This is a declared periodic-domain restriction, never shape-mismatch repair.
    The original arrays and selection are saved as evidence by callers.
    """
    try:
        if any(np.iscomplexobj(v) for v in (period, kx, resolution)):
            raise ValueError('PERIOD_SCALAR_COMPLEX')
        period, kx, resolution = float(period), float(kx), float(resolution)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError('PERIOD_SCALAR_INVALID') from exc
    if not all(math.isfinite(v) for v in (period, kx, resolution)) or min(period, resolution) <= 0:
        raise ValueError('PERIOD_SCALAR_INVALID')
    if np.iscomplexobj(x):
        raise ValueError('PERIOD_GRID_COMPLEX')
    x = np.asarray(x, float)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError('PERIOD_GRID_INVALID')
    if not np.isfinite(x).all():
        raise ValueError('PERIOD_GRID_NONFINITE')
    dx = 1.0 / resolution
    if not np.allclose(np.diff(x), dx, atol=1e-10, rtol=1e-10):
        raise ValueError('PERIOD_SPACING_MISMATCH')
    n = int(round(period * resolution))
    if abs(n * dx - period) > 1e-9:
        raise ValueError('PERIOD_NOT_INTEGER_GRID')
    keep = (x >= -period / 2 - 1e-9) & (x < period / 2 - 1e-9)
    if int(keep.sum()) != n:
        raise ValueError('PERIOD_COVERAGE_MISMATCH')
    selected = {}
    for name, a in arrays.items():
        a = np.asarray(a)
        if a.ndim == 0 or a.shape[0] != len(x):
            raise ValueError(f'PERIOD_COMPONENT_SHAPE_MISMATCH:{name}')
        if not np.isfinite(a).all():
            raise ValueError(f'PERIOD_COMPONENT_NONFINITE:{name}')
        scale = max(float(np.max(np.abs(a))), 1e-280)
        # All guards must be periodic copies of a retained point.
        for i in np.flatnonzero(~keep):
            xi = ((x[i] + period / 2) % period) - period / 2
            j = int(np.argmin(abs(x - xi)))
            if not keep[j] or abs(x[j] - xi) > 1e-8:
                raise ValueError('PERIOD_GUARD_HAS_NO_PARTNER')
            phase = np.exp(2j * np.pi * kx * (x[i] - x[j]))
            if np.max(np.abs(a[i] - a[j] * phase)) > 1e-6 * scale:
                raise ValueError(f'BLOCH_GUARD_MISMATCH:{name}')
        selected[name] = a[keep]
    return x[keep], selected, keep
