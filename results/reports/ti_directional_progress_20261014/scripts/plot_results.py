#!/usr/bin/env python3
"""Reproduce the current scientific figures from the sanitized JSON release.

Requires Python 3, NumPy, Matplotlib. Run from any working directory:
    python scripts/plot_results.py

No solver is invoked. Absorption estimates are not uncertainty intervals.
Sparse wavelength samples are deliberately not connected into a spectrum.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
import numpy as np

NAVY = "#1D4168"
TEAL = "#138A88"
ORANGE = "#D87835"
GRAY = "#6E7780"
LIGHT = "#D7DEE4"


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def setup():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10.5,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#7A858F",
        "axes.linewidth": 0.7,
        "grid.color": "#E3E8ED",
        "grid.linewidth": 0.6,
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
    })


def style(ax, title):
    ax.set_title(title, loc="left", pad=11)
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)


def header(fig, title, subtitle):
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=13.5, fontweight="bold")
    fig.text(0.08, 0.875, subtitle, ha="left", fontsize=10, color=GRAY)


def save(fig, folder, name):
    folder.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "svg"):
        fig.savefig(folder / f"{name}.{suffix}", dpi=240)
    plt.close(fig)


def metric(row, name):
    return float(row["metrics"][name])


def select(records, *, fill=None, wavelength=None, theta=None):
    return [r for r in records
            if (fill is None or np.isclose(float(r["fill"]), fill))
            and (wavelength is None or np.isclose(float(r["wavelength_um"]), wavelength))
            and (theta is None or np.isclose(float(r["theta_deg"]), theta))]


def one(records, **kwargs):
    found = select(records, **kwargs)
    if len(found) != 1:
        raise ValueError(f"Expected one point for {kwargs}, found {len(found)}")
    return found[0]


def fill_screen(records, out):
    fills = np.array([0.50, 0.65, 0.80])
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.2))
    fig.subplots_adjust(left=.095, right=.975, top=.73, bottom=.20, wspace=.34)
    header(fig, "Fill-fraction screening at 10.5 μm",
           "p polarization  |  incidence ±30°  |  r = 20  |  exploratory")
    ax, bx = axes
    style(ax, "(a) Direction-resolved absorption")
    for theta, color, lab in ((30, NAVY, "+30°"), (-30, TEAL, "−30°")):
        rows = [one(records, fill=f, wavelength=10.5, theta=theta) for f in fills]
        af = np.array([metric(r, "A_flux") for r in rows]) * 100
        av = np.array([metric(r, "A_vol") for r in rows]) * 100
        ax.plot(fills, af, "o--", lw=1, color=color, ms=6, label=f"{lab}, flux")
        ax.plot(fills, av, "s", mfc="white", mec=color, mew=1.4, ms=5.5,
                label=f"{lab}, volume", zorder=4)
    ax.set(xlabel="Air opening fraction F", ylabel="Absorptance (%)", ylim=(0, 25))
    ax.set_xticks(fills, ["0.50", "0.65", "0.80"])
    ax.legend(loc="lower left", ncol=2,
              frameon=False, handlelength=1.6, columnspacing=1)
    style(bx, "(b) Two-angle contrast")
    for name, marker, face, label in (("A_flux", "o", ORANGE, "Flux"),
                                      ("A_vol", "s", "white", "Volume")):
        delta = [100 * (metric(one(records, fill=f, wavelength=10.5, theta=30), name)
                        - metric(one(records, fill=f, wavelength=10.5, theta=-30), name))
                 for f in fills]
        bx.plot(fills, delta, marker=marker, ls="--" if name == "A_flux" else "none",
                color=ORANGE, mfc=face, mec=ORANGE, mew=1.3, ms=6, lw=1, label=label)
        if name == "A_flux":
            for x, y in zip(fills, delta):
                bx.annotate(f"{y:.2f}", (x, y), xytext=(0, 9), textcoords="offset points",
                            ha="center", fontsize=10, color=ORANGE)
    bx.set(xlabel="Air opening fraction F", ylabel="A(+30°) − A(−30°) (pp)", ylim=(0, 13.5))
    bx.set_xticks(fills, ["0.50", "0.65", "0.80"])
    bx.legend(loc="lower right", ncol=2, frameon=False)
    fig.text(.095, .026, "Dashed segments guide the eye; three samples do not establish a global optimum.",
             fontsize=9, color=GRAY)
    save(fig, out, "01_fill_comparison")


def wavelength_samples(records, out):
    waves = np.array([8., 10.5, 13.])
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 4.2))
    fig.subplots_adjust(left=.095, right=.975, top=.73, bottom=.20, wspace=.35)
    header(fig, "Three wavelength samples at F = 0.65",
           "p polarization  |  incidence ±30°  |  r = 20  |  exploratory")
    style(ax, "(a) Absorption samples")
    for theta, color, label in ((30, NAVY, "+30°"), (-30, TEAL, "−30°")):
        rows = [one(records, fill=.65, wavelength=w, theta=theta) for w in waves]
        for key, marker, face, text in (("A_flux", "o", color, "flux"),
                                        ("A_vol", "s", "white", "volume")):
            y = [100 * metric(r, key) for r in rows]
            ax.plot(waves, y, marker=marker, ls="none", color=color, mfc=face,
                    mec=color, mew=1.3, ms=6, label=f"{label}, {text}")
    ax.set(xlabel="Wavelength (μm)", ylabel="Absorptance (%)", ylim=(0, 27), xlim=(7.4, 13.6))
    ax.set_xticks(waves, ["8", "10.5", "13"])
    ax.legend(loc="lower left", ncol=2,
              frameon=False, columnspacing=1, handlelength=1)
    style(bx, "(b) Two-angle contrast samples")
    for key, marker, face, label in (("A_flux", "o", ORANGE, "Flux"),
                                     ("A_vol", "s", "white", "Volume")):
        y = [100 * (metric(one(records, fill=.65, wavelength=w, theta=30), key)
                    - metric(one(records, fill=.65, wavelength=w, theta=-30), key))
             for w in waves]
        bx.plot(waves, y, marker=marker, ls="none", color=ORANGE,
                mfc=face, mec=ORANGE, mew=1.3, ms=6, label=label)
        if key == "A_flux":
            for x, val in zip(waves, y):
                bx.annotate(f"{val:.2f}", (x, val), xytext=(0, 10), textcoords="offset points",
                            ha="center", color=ORANGE, fontsize=10)
    bx.set(xlabel="Wavelength (μm)", ylabel="A(+30°) − A(−30°) (pp)",
           ylim=(0, 14), xlim=(7.4, 13.6))
    bx.set_xticks(waves, ["8", "10.5", "13"])
    bx.legend(loc="lower right", ncol=2, frameon=False)
    fig.text(.095, .026, "Discrete points only: not a continuous 8–13 μm spectrum or a band-averaged emissivity.",
             fontsize=9, color=GRAY)
    save(fig, out, "02_wavelength_samples")


def residuals(records, out):
    rows = sorted(records, key=lambda r: (float(r["fill"]), float(r["wavelength_um"]),
                                          -float(r["theta_deg"])))
    x = np.arange(len(rows))
    colors = [NAVY if r["theta_deg"] > 0 else TEAL for r in rows]
    labels = [f"{float(r['fill']):.2f}\n{float(r['wavelength_um']):g}\n{float(r['theta_deg']):+g}"
              for r in rows]
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.4, 5.3), sharex=True)
    fig.subplots_adjust(left=.12, right=.975, top=.75, bottom=.245, hspace=.50)
    header(fig, "Residual diagnostics for the 10 current structured points",
           "Same material model  |  p polarization  |  r = 20  |  exploratory")
    style(ax, "(a) Flux–volume absorption difference")
    vals = np.array([100*(metric(r, "A_flux")-metric(r, "A_vol")) for r in rows])
    ax.bar(x, vals, width=.60, color=colors, zorder=3)
    ax.axhline(0, color=GRAY, lw=.8)
    ax.set(ylabel="Aflux − Avol (pp)", ylim=(min(vals)*1.5, max(vals)*1.6))
    for xx, val in zip(x, vals):
        ax.text(xx, val+(.025 if val > 0 else -.025), f"{val:+.3f}", ha="center",
                va="bottom" if val > 0 else "top", fontsize=8.5)
    handles = [Line2D([], [], marker="s", ls="none", color=NAVY, label="+30°"),
               Line2D([], [], marker="s", ls="none", color=TEAL, label="−30°")]
    ax.legend(handles=handles, ncol=2, frameon=False, loc="lower right", fontsize=9)
    style(bx, "(b) Signed net transmission")
    tv = np.array([100*metric(r, "T") for r in rows])
    bx.axhline(0, color=GRAY, lw=.8)
    bx.bar(x, tv, width=.60, color=colors, zorder=3)
    bx.set(ylabel="T (% of incident power)", ylim=(min(tv)*1.32, .012))
    bx.set_xticks(x, labels)
    bx.set_xlabel("Air opening fraction F / wavelength (μm) / incidence angle (°)", labelpad=10)
    fig.text(.12, .037, "Residuals are not error bars. Negative T is retained, not clipped or treated as physical gain.",
             fontsize=9, color=GRAY)
    save(fig, out, "03_numerical_residuals")


def source_diagnostic(data, out):
    rows = data["records"]
    rows = sorted(rows, key=lambda r: (0 if r["source_type"].lower().startswith("e") else 1,
                                       0 if r["polarization"].lower() == "s" else 1))
    labels = [f"{'Electric' if r['source_type'].lower().startswith('e') else 'Magnetic'}\n"
              f"{r['polarization']} ({r['source_component']})" for r in rows]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 4.2))
    fig.subplots_adjust(left=.095, right=.975, top=.73, bottom=.27, wspace=.36)
    header(fig, "Flat-interface source-type controls",
           "10.5 μm  |  normal incidence  |  r = 20  |  diagnostic, not a slot result")
    x = np.arange(len(rows))
    style(ax, "(a) Absorption by source convention")
    ax.plot(x, [100*metric(r, "A_flux") for r in rows], "o", ms=6.5, color=NAVY, label="Flux")
    ax.plot(x, [100*metric(r, "A_vol") for r in rows], "s", ms=6, mfc="white",
            mec=ORANGE, mew=1.5, label="Volume")
    ax.set_xticks(x, labels)
    ax.set(ylabel="Absorptance (%)", ylim=(3.1, 3.95), xlim=(-.5, 3.5))
    ax.legend(frameon=False, ncol=2, loc="upper right")
    style(bx, "(b) Flux–volume residual")
    rv = [100*(metric(r, "A_flux")-metric(r, "A_vol")) for r in rows]
    bx.bar(x, rv, width=.58, color=[TEAL, TEAL, ORANGE, ORANGE])
    bx.axhline(0, lw=.8, color=GRAY)
    bx.set_xticks(x, labels)
    bx.set(ylabel="Aflux − Avol (pp)", ylim=(-.34, .28), xlim=(-.5, 3.5))
    for xx, val in zip(x, rv):
        bx.text(xx, val + (.025 if val > 0 else -.025), f"{val:+.3f}",
                ha="center", va="bottom" if val > 0 else "top", fontsize=9)
    fig.text(.095, .075, "At each source type, normal-incidence s and p agree to numerical precision.",
             fontsize=9.5, color=GRAY)
    fig.text(.095, .028, "Electric vs magnetic Aflux differs by 0.371 pp; this is not a universal correction for slots.",
             fontsize=9, color=GRAY)
    save(fig, out, "04_source_type_diagnostic")


def geometry(data, out):
    rows = data["geometry"]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 4.7), sharey=True)
    fig.subplots_adjust(left=.08, right=.975, top=.73, bottom=.22, wspace=.17)
    header(fig, "Canonical 2D periodic blind-slot geometry",
           "Period P = 50 μm  |  axis length L = 40 μm  |  tilt = 30°  |  Ti thickness = 60 μm")
    for i, (ax, row) in enumerate(zip(axes, rows)):
        period = row["period_um"]
        thickness = row["thickness_um"]
        ax.add_patch(Rectangle((-period/2, 0), period, thickness,
                               facecolor="#C1CDD6", edgecolor=GRAY, lw=.7))
        vertices = np.array(row["vertices_x_depth"])
        # All finite cavity polygons intersecting the canonical periodic cell.
        # The JSON vertices are the output of the pure half-plane geometry validator.
        for shift in (-1, 0, 1):
            ax.add_patch(Polygon(vertices + [shift*period, 0], closed=True,
                                 facecolor="white", edgecolor=NAVY, lw=1.2))
        half = row["surface_opening_um"]/2
        ax.annotate("", xy=(-half, -4), xytext=(half, -4),
                    arrowprops={"arrowstyle": "<->", "color": TEAL, "lw": .8})
        ax.text(0, -7, f"w = {row['surface_opening_um']:g} μm", ha="center", va="bottom",
                color=TEAL, fontsize=9.5)
        ax.plot([0, 20], [0, 40*np.cos(np.pi/6)], "--", color=ORANGE, lw=.8)
        ax.axvline(-period/2, color=GRAY, ls=":", lw=.8)
        ax.axvline(period/2, color=GRAY, ls=":", lw=.8)
        ax.text(0, 53.5, f"Ti base\n{row['residual_um']:.2f} μm", ha="center", va="center",
                fontsize=9.5, color=NAVY)
        ax.set(xlim=(-period/2, period/2), ylim=(60, -13), xlabel="x (μm)")
        ax.set_aspect("equal")
        ax.set_xticks([-25, 0, 25])
        ax.set_yticks([0, 20, 40, 60])
        ax.set_title(f"({chr(97+i)}) F = {row['surface_fill']:.2f}", loc="left", pad=11)
    axes[0].set_ylabel("Depth into Ti (μm)")
    fig.text(.08, .074, "F = w/P is the air opening fraction. White regions are air cavities; colored regions are Ti.",
             fontsize=9, color=GRAY)
    fig.text(.08, .029, "At fixed tilted-axis length, changing F also changes maximum depth and residual base thickness.",
             fontsize=9, color=GRAY)
    save(fig, out, "05_geometry")


def flat_validation(data, out):
    rows = data["records"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
    fig.subplots_adjust(left=.10, right=.975, top=.73, bottom=.27, wspace=.35)
    header(fig, "Flat-interface Fresnel reflection checks",
           "Current material model  |  p: r = 20, s: r = 28  |  sampled wavelengths only")
    style(ax, "(a) p polarization, ±30°, r = 20")
    for theta, color, marker in ((30, NAVY, "o"), (-30, TEAL, "s")):
        group = sorted([r for r in rows if r["polarization"] == "p" and r["theta_deg"] == theta],
                       key=lambda r: r["wavelength_um"])
        ax.plot([r["wavelength_um"] for r in group], [r["R_error_pp"] for r in group],
                marker=marker, ls="none", ms=6, color=color,
                mfc=color if theta > 0 else "white", mew=1.3, label=f"{theta:+g}°")
    ax.axhline(.25, color=ORANGE, ls="--", lw=1)
    ax.text(19.8, .285, "0.25 pp gate", color=ORANGE, fontsize=9, ha="right")
    ax.set(xlabel="Wavelength (μm)", ylabel="|R − RFresnel| (pp)", ylim=(0, 1.4), xlim=(1.5, 21.5))
    ax.set_xticks([3, 8, 10.5, 13, 20], ["3", "8", "10.5", "13", "20"])
    ax.legend(frameon=False, ncol=2, loc="upper right")
    style(bx, "(b) s, +30°, r = 28")
    sg = [r for r in rows if r["polarization"] == "s"]
    bx.bar([0, 1], [r["R_error_pp"] for r in sg], color=TEAL, width=.5)
    bx.axhline(.25, color=ORANGE, ls="--", lw=1)
    bx.set(ylim=(0, .75), ylabel="|R − RFresnel| (pp)")
    bx.set_xticks([0, 1], ["0.20", "0.10"])
    bx.set_xlabel("Source fwidth fraction\nat 10.5 μm")
    for i, r in enumerate(sg):
        bx.text(i, r["R_error_pp"]+.03, f"{r['R_error_pp']:.3f}", ha="center", fontsize=10)
    fig.text(.10, .075, "4 of 12 controls pass their strict flat gate; 8 remain diagnostic failures.",
             fontsize=9.5, color=GRAY)
    fig.text(.10, .028, "A completed solve is not an accuracy certificate; no continuous-band qualification is claimed.",
             fontsize=9, color=GRAY)
    save(fig, out, "06_flat_validation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    setup()
    data_dir = args.root / "data"
    out_dir = args.root / "figures" / "current"
    structured = load(data_dir / "current_structured_results.json")
    records = structured["records"]
    if len(records) != 10:
        raise ValueError("The release figure set expects exactly 10 current structured points")
    if any(r["polarization"].lower() != "p" for r in records):
        raise ValueError("Current screening figures only apply to p polarization")
    fill_screen(records, out_dir)
    wavelength_samples(records, out_dir)
    residuals(records, out_dir)
    source_diagnostic(load(data_dir / "source_type_diagnostic.json"), out_dir)
    geometry(load(data_dir / "geometry_and_settings.json"), out_dir)
    flat_validation(load(data_dir / "flat_validation_results.json"), out_dir)
    print(f"Created 6 figures in PNG and SVG formats: {out_dir}")


if __name__ == "__main__":
    main()
