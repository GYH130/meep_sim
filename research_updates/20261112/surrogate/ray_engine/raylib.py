"""Thin Python driver around the M3-extended ray tracer (./ray_material).

FEM geometry mapping (floor-fill SiO2 of thickness t in an A3 frustum of height H measured
from the Ti floor, floor across-flats d, opening across-flats 69.5455842631299 um, pitch 80):
the air-filled part of the cavity is itself an A3 frustum of height H - t whose bottom is the
flat film top (apothem a_t = d/2 + (34.77279213156495 - d/2) t / H) and whose top is the
unchanged opening. Its floor is modelled as a planar SiO2(t)/Ti stack (surface M3).
"""
import hashlib
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import materials as M

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "ray_material")
CACHE = os.path.join(HERE, "runs", "cache")
OPEN_ACROSS = 69.5455842631299
OPEN_APOTHEM = 34.77279213156495
PITCH = 80.0


def film_top_apothem(H, t, d_bottom=15.0):
    a0 = d_bottom / 2.0
    return a0 + (OPEN_APOTHEM - a0) * t / H


def build_args(H, t, lam, theta, phi, rays, seed, d_bottom=15.0, surface="M3", max_bounces=128):
    nti, kti = M.ti_nk(lam)
    if surface == "M3":
        nf, kf = M.sio2_nk(lam)
        a_t = film_top_apothem(H, t, d_bottom)
        height = H - t
        bottom = 2 * a_t
    elif surface == "M1":
        nf = kf = 0.0
        height = H
        bottom = d_bottom
    else:
        raise ValueError(surface)
    args = ["--structure", "A3", "--surface", surface, "--height-um", repr(float(height)),
            "--bottom-width-um", repr(float(bottom)), "--upper-width-um", repr(OPEN_ACROSS), "--pitch-um", repr(PITCH),
            "--theta", repr(float(theta)), "--phi", repr(float(phi)), "--n", repr(nti), "--k", repr(kti),
            "--wavelength-um", repr(float(lam)), "--rays", str(int(rays)), "--seed", str(int(seed)),
            "--max-bounces", str(int(max_bounces))]
    if surface == "M3":
        args += ["--film-n", repr(nf), "--film-k", repr(kf), "--film-thickness-um", repr(float(t))]
    return args


def run_one(H, t, lam, theta, phi, rays, seed, d_bottom=15.0, surface="M3", max_bounces=128, use_cache=True):
    args = build_args(H, t, lam, theta, phi, rays, seed, d_bottom, surface, max_bounces)
    key = hashlib.sha256((" ".join(args) + "|" + _exe_hash()).encode()).hexdigest()[:24]
    path = os.path.join(CACHE, key + ".json")
    if use_cache and os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    out = subprocess.run([EXE] + args, check=True, capture_output=True, text=True).stdout
    d = json.loads(out)
    d["_args"] = args
    d["_design"] = dict(H=H, t=t, d_bottom=d_bottom, lam=lam, theta=theta, phi=phi, rays=rays, seed=seed, surface=surface)
    os.makedirs(CACHE, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(d, fh)
    return d


_EXE_HASH = None


def _exe_hash():
    global _EXE_HASH
    if _EXE_HASH is None:
        with open(os.path.join(HERE, "ray_material.cpp"), "rb") as fh:
            _EXE_HASH = hashlib.sha256(fh.read()).hexdigest()[:16]
    return _EXE_HASH


def run_many(jobs, workers=8):
    """jobs: list of dicts of run_one kwargs. Returns list of results in order."""
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_one, **j) for j in jobs]
        return [f.result() for f in futs]


def channel(res, pol):
    """pol in {'s','p','unpolarized'}; returns dict of tallies (fractions) incl. Ti total."""
    c = res["channels"][pol]
    out = {k: c[k] for k in ("A_total", "A_active", "A_wall", "escape", "residual", "A_film", "A_Ti_floor", "A_land")}
    out["A_Ti_total"] = c["A_wall"] + c["A_Ti_floor"]
    out["A_wall_only"] = c["A_wall"] - c["A_land"]
    out["se_A_total"] = c["se"]["A_total"]
    out["se_A_film"] = c["se"]["A_film"]
    out["se_A_Ti_total_approx"] = (c["se"]["A_wall"] ** 2 + c["se"]["A_Ti_floor"] ** 2) ** 0.5
    return out
