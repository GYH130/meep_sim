"""Optical constants used by the ray/FEM comparison.

Ti  : Ordal table (material.json "source_table"), piecewise-linear n,k vs wavelength.
SiO2: Kischkat (refractiveindex.info YAML, "tabulated nk"), piecewise-linear n,k vs wavelength,
      valid only 1.53846-14.28571 um (no extrapolation).
Convention: exp(-i omega t), N = n + i k (k >= 0 passive).
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TI_JSON = os.path.join(HERE, "inputs", "Ti_Ordal_material.json")
SIO2_YML = os.path.join(HERE, "inputs", "SiO2_Kischkat.yml")

# Values quoted for the FEM runs (used as regression checks of the interpolation).
FEM_TI = {9.0: (3.6382215, 17.731697), 10.5: (4.308219763636363, 20.63042690909091)}
FEM_SIO2 = {9.0: (0.5947997319778178, 1.7674602002464557), 10.5: (2.2447212969283283, 0.06982232081911283)}


def _load_ti():
    d = json.load(open(TI_JSON))
    rows = d["source_table"]
    wl = np.array([r["wavelength_um"] for r in rows], float)
    n = np.array([r["n"] for r in rows], float)
    k = np.array([r["k"] for r in rows], float)
    o = np.argsort(wl)
    return wl[o], n[o], k[o]


def _load_sio2():
    rows = []
    in_data = False
    for line in open(SIO2_YML):
        s = line.strip()
        if s.startswith("data:"):
            in_data = True
            continue
        if in_data:
            parts = s.split()
            if len(parts) == 3:
                try:
                    rows.append([float(p) for p in parts])
                    continue
                except ValueError:
                    pass
            if rows:
                break
    a = np.array(rows, float)
    o = np.argsort(a[:, 0])
    return a[o, 0], a[o, 1], a[o, 2]


TI_WL, TI_N, TI_K = _load_ti()
SIO2_WL, SIO2_N, SIO2_K = _load_sio2()


def _interp(wl_tab, n_tab, k_tab, lam, name):
    lam = float(lam)
    if lam < wl_tab[0] - 1e-12 or lam > wl_tab[-1] + 1e-12:
        raise ValueError(f"{name}: wavelength {lam} um outside table range {wl_tab[0]}-{wl_tab[-1]} (no extrapolation)")
    return float(np.interp(lam, wl_tab, n_tab)), float(np.interp(lam, wl_tab, k_tab))


def ti_nk(lam_um):
    return _interp(TI_WL, TI_N, TI_K, lam_um, "Ti(Ordal)")


def sio2_nk(lam_um):
    return _interp(SIO2_WL, SIO2_N, SIO2_K, lam_um, "SiO2(Kischkat)")


def check_against_fem(verbose=True):
    out = []
    for lam, (n, k) in FEM_TI.items():
        nn, kk = ti_nk(lam)
        out.append(("Ti", lam, n, k, nn, kk, abs(nn - n), abs(kk - k)))
    for lam, (n, k) in FEM_SIO2.items():
        nn, kk = sio2_nk(lam)
        out.append(("SiO2", lam, n, k, nn, kk, abs(nn - n), abs(kk - k)))
    if verbose:
        for r in out:
            print("%-5s %5.2f um  FEM n=%.15g k=%.15g | interp n=%.15g k=%.15g | |dn|=%.2e |dk|=%.2e" % r)
    return out


if __name__ == "__main__":
    print("Ti table rows:", len(TI_WL), "range", TI_WL[0], TI_WL[-1])
    print("SiO2 table rows:", len(SIO2_WL), "range", SIO2_WL[0], SIO2_WL[-1])
    res = check_against_fem()
    worst = max(max(r[6], r[7]) for r in res)
    print("worst |delta| vs FEM-quoted values:", worst)
    for lam in (8.0, 9.0, 10.5, 11.5, 12.5, 13.0):
        print(f"lambda={lam:5.2f}  Ti n,k={ti_nk(lam)}  SiO2 n,k={sio2_nk(lam)}")
