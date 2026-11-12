#!/usr/bin/env python3
"""Summarize and plot the re-scored saved M0 geometry screen."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import PercentFormatter
import numpy as np


POWER = "cone_power_over_full_hemisphere_blackbody"
FRACTION = "cone_fraction_of_emitted_power"
MEAN = "cone_average_ideal_emissivity"
HEMISPHERE = "hemisphere_power_over_full_hemisphere_blackbody"
REFERENCE = (10, 15, 20, 25, 30, 40, 45, 60)
EXAMPLES = {"A0_G23": "#0F7884", "A0_G99": "#C45934", "A0_G26": "#654F94"}


def summaries(report: dict) -> dict:
    coarse = {row["shape"]["id"]: row for row in report["cases"] if row["phase"] == "coarse"}
    fine = {row["shape"]["id"]: row for row in report["cases"] if row["phase"] == "refine"}
    require_counts = (len(coarse), len(fine))
    if require_counts != (200, 30):
        raise ValueError(f"Expected 200 coarse and 30 refined cases, got {require_counts}")

    def angle_rows(collection, radius):
        return {sid: record["normal_cones"][str(radius)] for sid, record in collection.items()}

    angle_ranking = {}
    for radius in REFERENCE:
        values = angle_rows(coarse, radius)
        strength = sorted(values, key=lambda key: (-values[key][POWER], key))
        concentration = sorted(values, key=lambda key: (-round(values[key][FRACTION], 12), -values[key][POWER], key))
        heuristic = sorted(values, key=lambda key: (-values[key][MEAN] * values[key][FRACTION], key))
        front = []
        for key, v in values.items():
            if not any(r[POWER] >= v[POWER] - 1e-12 and r[FRACTION] >= v[FRACTION] - 1e-12
                       and (r[POWER] > v[POWER] + 1e-12 or r[FRACTION] > v[FRACTION] + 1e-12)
                       for r in values.values()):
                front.append(key)
        angle_ranking[str(radius)] = {
            "coarse_top_normalized_power": strength[:3],
            "coarse_top_concentration_with_power_tiebreak": concentration[:3],
            "coarse_top_old_product_heuristic": heuristic[:3],
            "nondominated_among_coarse": sorted(front),
            "nondominated_without_refinement": sorted(set(front) - fine.keys()),
        }

    quality = {}
    for radius in (15, 30):
        coarse_at = angle_rows(coarse, radius)
        fine_at = angle_rows(fine, radius)
        quality[str(radius)] = {
            "max_abs_coarse_fine_shift_percentage_points": {
                key: 100 * max(abs(row[key] - coarse_at[sid][key]) for sid, row in fine_at.items())
                for key in (MEAN, FRACTION, POWER)
            },
            "near_saturated_concentration_count": {
                "coarse": sum(row[FRACTION] >= 0.9999 for row in coarse_at.values()),
                "refine": sum(row[FRACTION] >= 0.9999 for row in fine_at.values()),
            },
        }

    example_metrics = {}
    for sid in EXAMPLES:
        row = fine[sid]
        example_metrics[sid] = {
            "shape": row["shape"],
            HEMISPHERE: row[HEMISPHERE],
            "normal_cones": {str(radius): row["normal_cones"][str(radius)] for radius in (15, 30, 40, 60)},
            "off_axis_30deg_around_phi0": {
                str(axis["axis_theta_deg"]): axis["cones"]["30"]
                for axis in row["off_axis_target_cones"] if axis["axis_phi_deg"] == 0},
        }
    return {
        "scope": "Saved M0 ideal ray data. Post-processing only. 200 coarse cases and 30 previously computed refinements, including 6 earlier refinements. No spectral or temperature data.",
        "angle_ranking_sensitivity": angle_ranking,
        "paired_coarse_fine_sensitivity": quality,
        "illustrative_fine_cases": example_metrics,
        "sources": report["sources"],
    }


def font_style() -> bool:
    cjk = Path(os.environ.get("MEEP_CJK_FONT", "Arial Unicode.ttf"))
    if cjk.is_file():
        font_manager.fontManager.addfont(str(cjk))
    plt.rcParams.update({
        "font.family": "Arial Unicode MS" if cjk.is_file() else "DejaVu Sans",
        "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
        "axes.labelcolor": "#1D3448", "axes.edgecolor": "#A1B0BE", "grid.color": "#E3E9EE",
        "text.color": "#1D3448", "xtick.color": "#4D6275", "ytick.color": "#4D6275",
        "axes.unicode_minus": False, "svg.fonttype": "none", "legend.frameon": False,
    })
    return cjk.is_file()


def save(fig, path: Path) -> None:
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight", dpi=220)
    plt.close(fig)


def normal_plot(report: dict, output: Path, cjk: bool) -> None:
    fine = {row["shape"]["id"]: row for row in report["cases"] if row["phase"] == "refine"}
    radii = np.array(report["normal_half_angles_deg"])
    fig, (power_ax, fraction_ax) = plt.subplots(1, 2, figsize=(11.2, 4.3), layout="constrained")
    fig.suptitle("同一结构在不同法向角锥内的输出与集中度" if cjk else "Normal cone throughput and concentration", x=.02, ha="left", fontsize=17)
    for sid, color in EXAMPLES.items():
        values = fine[sid]["normal_cones"]
        label = sid.replace("_", "–")
        for ax, key in ((power_ax, POWER), (fraction_ax, FRACTION)):
            ax.plot(radii, [values[str(a)][key] for a in radii], color=color, linewidth=2.1, label=label)
    fraction_ax.plot(radii, np.sin(np.radians(radii)) ** 2, color="#778A9B", linestyle=":", linewidth=1.7,
                     label="朗伯面" if cjk else "Lambertian")
    for ax in (power_ax, fraction_ax):
        for radius in (15, 30):
            ax.axvline(radius, color="#9CAFBF", linestyle="--", linewidth=.85)
        ax.set(xlim=(5, 85), xlabel="法向角锥半角 (°)" if cjk else "Normal cone half-angle (°)")
        ax.set_xticks((5, 15, 30, 45, 60, 75, 85))
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.grid(True)
        ax.set_axisbelow(True)
    power_ax.set_ylabel("锥内功率 / 同面积黑体半球功率" if cjk else "Cone power / footprint blackbody hemisphere")
    fraction_ax.set_ylabel("锥内功率占半球比例" if cjk else "Cone power / emitted hemisphere")
    power_ax.set(ylim=(0, .205), yticks=np.arange(0, .201, .05))
    fraction_ax.set_ylim(0, 1.03)
    power_ax.legend(loc="lower right")
    fraction_ax.legend(loc="lower right", ncols=2)
    save(fig, output / "normal_cone_sensitivity")


def off_axis_plot(report: dict, output: Path, cjk: bool) -> None:
    data = next(row for row in report["cases"] if row["phase"] == "refine" and row["shape"]["id"] == "A0_G99")
    angle = np.array(data["sampled_theta_deg"])
    indices = {name: i for i, name in enumerate(data["azimuth_sample_columns"])}
    phi = np.array([sample[indices["phi_deg"]] for sample in data["azimuth_samples"][1]])
    field = np.array([[sample[indices["A_total"]] for sample in ring] * (len(phi) if i == 0 else 1)
                      for i, ring in enumerate(data["azimuth_samples"])])
    fig = plt.figure(figsize=(11.2, 4.2), layout="constrained")
    axes = fig.add_gridspec(1, 3, width_ratios=(1, 1, 1.55))
    fig.suptitle("A0–G99 的二维角域和离轴接收锥" if cjk else "A0–G99 angular field and off-axis receiver cones",
                 x=.02, ha="left", fontsize=17)
    for index, (values, title) in enumerate(((field, "方向发射率" if cjk else "Directional emissivity"),
                                             (field * np.cos(np.radians(angle))[:, None],
                                              "单位立体角投影输出" if cjk else "Projected output per solid angle"))):
        ax = fig.add_subplot(axes[0, index], projection="polar")
        step = 360 / len(phi)
        edges = np.radians(np.arange(len(phi) + 1) * step - step / 2)
        r_edges = np.r_[0, (angle[1:] + angle[:-1]) / 2, 90]
        ax.grid(False)
        image = ax.pcolormesh(edges, r_edges, values, cmap="cividis", vmin=0, vmax=1)
        ax.set_rgrids([15, 30, 60], angle=45, color="white", fontsize=9)
        ax.set_thetagrids([0, 90, 180, 270])
        ax.set_ylim(0, 90)
        ax.set_title(title, fontsize=11, pad=12)
        ax.grid(True, color="white", linewidth=.4, alpha=.65)
    fig.colorbar(image, ax=fig.axes[:2], orientation="horizontal", shrink=.65, pad=.08,
                 format=PercentFormatter(1, decimals=0), ticks=(0, .25, .5, .75, 1))
    ax = fig.add_subplot(axes[0, 2])
    for phi_center, color in zip((0, 90, 180, 270), ("#0F7884", "#C45934", "#654F94", "#427D57")):
        profile = [axis for axis in data["off_axis_target_cones"] if axis["axis_phi_deg"] == phi_center]
        x = [0, *(axis["axis_theta_deg"] for axis in profile)]
        y = [data["normal_cones"]["30"][POWER], *(axis["cones"]["30"][POWER] for axis in profile)]
        ax.plot(x, y, marker="o", linewidth=1.4, markersize=4, color=color, label=f"φ₀={phi_center}°")
    ax.set(xticks=(0, 15, 30, 45), xlabel="接收轴偏离法线 (°)" if cjk else "Receiver-axis tilt (°)",
           ylabel="锥内归一化功率" if cjk else "Normalized cone power", ylim=(0, .2))
    ax.set_title("接收锥半角固定为 30°" if cjk else "30° receiver-cone half-angle", fontsize=11)
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.grid(True)
    ax.legend(loc="upper right", fontsize=9)
    save(fig, output / "angular_map_receiver_axes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics", type=Path)
    args = parser.parse_args()
    output = args.metrics.parent
    report = json.loads(args.metrics.read_text())
    compact = summaries(report)
    (output / "key_sensitivity.json").write_text(json.dumps(compact, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    cjk = font_style()
    normal_plot(report, output, cjk)
    off_axis_plot(report, output, cjk)


if __name__ == "__main__":
    main()
