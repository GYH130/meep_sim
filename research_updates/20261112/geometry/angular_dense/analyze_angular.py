"""Turn complete M0 angle-resolved ray summaries into plots and cone metrics."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import os

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np


NAVY = "#17314B"
MUTED = "#536A7C"
GRID = "#E5ECF1"
WALL = "#DAE5EE"
FAMILIES = {"A0": "#077C87", "A3": "#BE5B31"}
PAIR_IDS = ("G01", "G06", "G09")
CONES = (15, 30)


def setup_style() -> bool:
    cjk_font = Path(os.environ.get("MEEP_CJK_FONT", "Arial Unicode.ttf"))
    has_cjk = cjk_font.is_file()
    if has_cjk:
        font_manager.fontManager.addfont(str(cjk_font))
    plt.rcParams.update({
        "font.family": "Arial Unicode MS" if has_cjk else "DejaVu Sans",
        "font.size": 10, "pdf.fonttype": 42, "axes.unicode_minus": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#B7C5CF", "axes.labelcolor": NAVY,
        "axes.titlecolor": NAVY, "axes.titlesize": 12, "axes.titlepad": 9,
        "xtick.color": MUTED, "ytick.color": MUTED, "grid.color": GRID,
        "legend.frameon": False, "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    return has_cjk


def read_phase(path: Path, phase: str) -> tuple[dict, dict[str, "AngularGeometry"]]:
    data = json.loads(path.read_text())
    if data.get("phase", phase) != phase:
        raise ValueError(f"{path} is not the {phase} phase")
    if not str(data.get("state", "")).startswith("COMPLETE") or data.get("failed"):
        raise ValueError(f"{path} is not a completed, failure-free phase")
    rows = data["geometry_summaries"]
    count = sum(row["completed_angles"] for row in rows)
    expected_count = data["condition_count"]
    if count != expected_count or data["completed_conditions"] != expected_count:
        raise ValueError(f"{path} has missing condition results")
    if len(rows) != data["completed_geometries"] or len(rows) != data["planned_geometry_count"]:
        raise ValueError(f"{path} has missing geometry results")
    geometries = {}
    for row in rows:
        result = AngularGeometry.from_row(row)
        if result.name in geometries:
            raise ValueError(f"Duplicate geometry {result.name}")
        geometries[result.name] = result
    return data, geometries


@dataclass
class AngularGeometry:
    shape: dict
    theta: np.ndarray
    average: np.ndarray
    min_a: np.ndarray
    max_a: np.ndarray
    residual: np.ndarray
    azimuth: np.ndarray | None
    azimuth_values: np.ndarray | None

    @property
    def name(self) -> str:
        return self.shape["id"]

    @property
    def family(self) -> str:
        return self.shape["family"]

    @classmethod
    def from_row(cls, row: dict) -> "AngularGeometry":
        angles = sorted(row["angles"], key=lambda sample: sample["theta_deg"])
        if len(angles) < 3:
            raise ValueError("Need at least three polar angle rings")
        theta = np.asarray([sample["theta_deg"] for sample in angles], dtype=float)
        if theta[0] != 0 or np.any(np.diff(theta) <= 0) or theta[-1] < 89 or theta[-1] > 90:
            raise ValueError(f"Incomplete or invalid polar hemisphere for {row['shape']['id']}")
        counts = [sample["azimuth_count"] for sample in angles]
        if (counts[0] != 1 or len(set(counts[1:])) != 1 or counts[1] < 3
                or sum(counts) != row["completed_angles"]
                or row["completed_angles"] != row["expected_angles"]):
            raise ValueError(f"Incomplete angle rings for {row['shape']['id']}")
        average, min_a, max_a, residual = (
            np.asarray([sample[key] for sample in angles], dtype=float)
            for key in ("mean_A", "min_A", "max_A", "max_residual")
        )
        arrays = [theta, average, min_a, max_a, residual]
        if any(not np.all(np.isfinite(values)) for values in arrays):
            raise ValueError(f"Nonfinite angular data for {row['shape']['id']}")
        if (np.any(min_a < -1e-10) or np.any(max_a > 1 + 1e-10)
                or np.any(average < min_a - 1e-10) or np.any(average > max_a + 1e-10)
                or np.any(residual < 0) or np.any(residual > 1)):
            raise ValueError(f"Invalid absorptance or residual for {row['shape']['id']}")

        azimuth = azimuth_values = None
        if any("azimuth_samples" in sample for sample in angles):
            azimuth = np.arange(counts[1]) * (360 / counts[1])
            matrix = []
            for i, sample in enumerate(angles):
                ring = np.asarray(sample["azimuth_samples"], dtype=float)
                expected_phi = [0] if i == 0 else azimuth
                if (ring.shape != (counts[i], 4)
                        or not np.all(np.isfinite(ring))
                        or not np.allclose(ring[:, 0], expected_phi)
                        or not np.isclose(ring[:, 1].mean(), average[i], atol=1e-10, rtol=0)):
                    raise ValueError(f"Azimuth summary disagrees with samples for {row['shape']['id']}")
                matrix.append(np.repeat(ring[0, 1], counts[1]) if i == 0 else ring[:, 1])
            azimuth_values = np.asarray(matrix)
        return cls(row["shape"], theta, average, min_a, max_a, residual, azimuth, azimuth_values)

    def cumulative(self, half_angles) -> np.ndarray:
        # The normalized projected-solid-angle element is d(sin²θ).
        # Trapezoidal integration in this coordinate preserves a Lambertian reference exactly.
        u = np.sin(np.deg2rad(self.theta)) ** 2
        a = self.average
        if u[-1] < 1:
            u = np.append(u, 1.0)
            a = np.append(a, a[-1])
        integral = np.concatenate(([0.0], np.cumsum(np.diff(u) * (a[:-1] + a[1:]) / 2)))
        x = np.sin(np.deg2rad(np.asarray(half_angles, dtype=float))) ** 2
        idx = np.clip(np.searchsorted(u, x, side="right") - 1, 0, len(u) - 2)
        dx = x - u[idx]
        slope = (a[idx + 1] - a[idx]) / (u[idx + 1] - u[idx])
        return integral[idx] + dx * (a[idx] + slope * dx / 2)

    def metrics(self) -> dict:
        hemisphere = float(self.cumulative(90))
        shape = self.shape
        floor_fraction = (shape["bottom_width_um"] / shape["pitch_um"]) ** 2
        metrics = {
            "shape": shape,
            "hemispherical_ideal_emissivity": hemisphere,
            "floor_area_fraction_of_cell": floor_fraction,
            "hemi_minus_floor_area_fraction": hemisphere - floor_fraction,
            "opening_area_fraction_of_cell": (shape["upper_width_um"] / shape["pitch_um"]) ** 2,
            "land_width_along_neighbor_centers_um": shape["pitch_um"] - shape["upper_width_um"],
            "height_to_opening_ratio": shape["height_um"] / shape["upper_width_um"],
            "max_residual": float(np.max(self.residual)),
            "max_absorptance_azimuth_range": float(np.max(self.max_a - self.min_a)),
            "cones": {},
        }
        for angle in CONES:
            x = np.sin(np.deg2rad(angle)) ** 2
            power = float(self.cumulative(angle))
            outside = max(0.0, (hemisphere - power) / (1 - x))
            inside = power / x
            metrics["cones"][str(angle)] = {
                "mean_ideal_emissivity_inside": inside,
                "power_fraction": power / hemisphere if hemisphere > 0 else None,
                "mean_ideal_emissivity_outside": outside,
                "mean_inside_to_outside_ratio": inside / outside if outside > 1e-14 else None,
                "cell_normalized_power": power,
            }
        return metrics


def check_refinement(coarse: dict, refined: dict) -> None:
    if not set(refined).issubset(coarse):
        raise ValueError("Refinement contains unmatched shapes")
    for name, finer in refined.items():
        if finer.shape != coarse[name].shape:
            raise ValueError(f"Shape parameters differ across phases for {name}")


def phase_comparison(coarse: dict, refined: dict) -> dict:
    result = {}
    for angle in CONES:
        key = str(angle)
        result[key] = {}
        for metric in ("mean_ideal_emissivity_inside", "power_fraction"):
            phases = [coarse, refined]
            comparable = [
                name for name in sorted(refined)
                if all(phase[name]["cones"][key][metric] is not None for phase in phases)
            ]
            ranks = [
                sorted(comparable, key=lambda name: (-phase[name]["cones"][key][metric], name))
                for phase in phases
            ]
            delta = {
                name: 100 * (refined[name]["cones"][key][metric] - coarse[name]["cones"][key][metric])
                for name in comparable
            }
            result[key][metric] = {
                "coarse_rank_among_refined_shapes": ranks[0],
                "refined_rank": ranks[1],
                "rank_unchanged": ranks[0] == ranks[1],
                "refine_minus_coarse_percentage_points": delta,
                "max_abs_shift_pp": max((abs(value) for value in delta.values()), default=None),
            }
    return result


def cpc_width(z: np.ndarray) -> np.ndarray:
    angle = np.pi / 12
    sin, cos = np.sin(angle), np.cos(angle)
    k = 20 * (1 + sin)
    return (-sin * (k + cos * z) + np.sqrt(k * k + 2 * k * cos * z)) / (cos * cos) - 10


def profile(result: AngularGeometry) -> tuple[np.ndarray, np.ndarray]:
    shape = result.shape
    bottom, opening, height = (shape[key] for key in ("bottom_width_um", "upper_width_um", "height_um"))
    segments = shape["curve_segments"] if result.family == "A0" else 1
    t = np.linspace(0, 1, segments + 1)
    if result.family == "A0":
        relative_width = (cpc_width(90 * t) - 10) / (cpc_width(90) - 10)
    else:
        relative_width = t
    return (bottom + (opening - bottom) * relative_width) / 2, height * t


def letter(ax, label: str, x: float = -0.10):
    ax.text(x, 1.045, f"({label})", transform=ax.transAxes, color=MUTED,
            fontsize=9, va="bottom")


def setup_plot(ax):
    ax.grid(True, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3)


def save_figure(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=220)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def wall_panel(ax, curved: AngularGeometry, straight: AngularGeometry, cjk: bool) -> None:
    shape = curved.shape
    height, bottom, opening, pitch = (
        shape[key] for key in ("height_um", "bottom_width_um", "upper_width_um", "pitch_um")
    )
    ax.fill_betweenx([0, height], [-pitch / 2, -pitch / 2], [pitch / 2, pitch / 2],
                     color=WALL, alpha=0.5, linewidth=0)
    r, z = profile(curved)
    ax.fill_betweenx(z, -r, r, color="white", linewidth=0)
    for result in (curved, straight):
        r, z = profile(result)
        style = "-" if result.family == "A0" else "--"
        for sign in (-1, 1):
            ax.plot(sign * r, z, style, lw=1.9, color=FAMILIES[result.family], zorder=3)
    ax.plot([-pitch / 2, -opening / 2], [height, height], lw=1, color=MUTED)
    ax.plot([opening / 2, pitch / 2], [height, height], lw=1, color=MUTED)
    ax.plot([-bottom / 2, bottom / 2], [0, 0], color=NAVY, lw=3, solid_capstyle="butt")
    ax.annotate("", xy=(-opening / 2, height + 5), xytext=(opening / 2, height + 5),
                arrowprops={"arrowstyle": "<->", "color": MUTED, "lw": 0.85})
    ax.text(0, height + 8, f"D = {opening:.1f} μm", ha="center", va="bottom", color=NAVY, fontsize=9)
    ax.text(0, -5, f"b = {bottom:.1f} μm  ·  H = {height:.1f} μm",
            ha="center", va="top", color=NAVY, fontsize=9)
    ax.set(xlim=(-pitch / 2 - 4, pitch / 2 + 4), ylim=(-16, height + 24), aspect="equal")
    ax.set_title(f"{curved.name.split('_', 1)[1]}  {'同尺寸截面' if cjk else 'matched cross-section'}",
                 loc="left", fontsize=11)
    ax.axis("off")


def paired_plot(coarse: dict, refined: dict, out: Path, cjk: bool) -> None:
    pairs = []
    for pair_id in PAIR_IDS:
        names = [f"{family}_{pair_id}" for family in FAMILIES]
        if all(name in coarse for name in names):
            pairs.append([refined.get(name, coarse[name]) for name in names])
    if not pairs:
        return
    fig, axes = plt.subplots(len(pairs), 3, figsize=(14, 3.2 * len(pairs) + 0.75), squeeze=False)
    fig.subplots_adjust(left=0.055, right=0.975, top=0.855, bottom=0.105, wspace=0.36, hspace=0.53)
    fig.text(0.055, 0.965, "同尺寸腔体的角域对比" if cjk else "Angular response of matched cavities",
             size=20, color=NAVY, weight="bold")
    title = (
        "M0 理想射线模型  ·  黑底、镜面侧壁  ·  全方位角平均  ·  按原胞面积归一化"
        if cjk else "M0 ideal rays  ·  black floor and mirror walls  ·  full azimuth average  ·  unit-cell area"
    )
    fig.text(0.055, 0.927, title, size=11, color=MUTED)
    if cjk:
        legends = ["A0 曲壁", "A3 直壁", "朗伯面"]
    else:
        legends = ["A0 curved", "A3 straight", "Lambertian"]
    handles = [
        plt.Line2D([0], [0], color=FAMILIES["A0"], lw=2),
        plt.Line2D([0], [0], color=FAMILIES["A3"], lw=2, linestyle="--"),
        plt.Line2D([0], [0], color="#7E8E9C", lw=1.4, linestyle=":"),
    ]
    fig.legend(handles, legends, loc="upper right", bbox_to_anchor=(0.98, 0.93), ncol=3)
    for i, (curved, straight) in enumerate(pairs):
        geom, angular, cone = axes[i]
        wall_panel(geom, curved, straight, cjk)
        for ax in (angular, cone):
            setup_plot(ax)
            for half_angle in CONES:
                ax.axvline(half_angle, color="#B1C0CC", linestyle=":", lw=0.9, zorder=1)
            ax.set_xlim(0, 90)
            ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
            ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        polar = np.linspace(0, 90, 361)
        for result in (curved, straight):
            color = FAMILIES[result.family]
            style = "-" if result.family == "A0" else "--"
            angular.fill_between(result.theta, result.min_a, result.max_a, color=color, alpha=0.13, lw=0)
            angular.plot(result.theta, result.average, linestyle=style, color=color, lw=2)
            if result.name in refined:
                rough = coarse[result.name]
                angular.scatter(rough.theta, rough.average, s=9, color=color, alpha=0.40, edgecolors="none")
            power = result.cumulative(polar)
            if power[-1] > 0:
                cone.plot(polar, power / power[-1], linestyle=style, color=color, lw=2)
                x = np.asarray(CONES)
                cone.scatter(x, result.cumulative(x) / power[-1], color=color, s=24, zorder=4,
                             marker="o" if result.family == "A0" else "s")
        cone.plot(polar, np.sin(np.deg2rad(polar)) ** 2, color="#7E8E9C", linestyle=":", lw=1.5)
        angular.set_ylim(0, 1)
        cone.set_ylim(0, 1)
        angular.set_xlabel("极角 θ (°)" if cjk else "Polar angle θ (°)")
        cone.set_xlabel("角锥半角 (°)" if cjk else "Cone half-angle (°)")
        angular.set_ylabel("理想方向发射率" if cjk else "Ideal directional emissivity")
        cone.set_ylabel("锥内功率占半球比例" if cjk else "Fraction of hemispheric power")
        if i == 0:
            angular.set_title("各方向有多强" if cjk else "Emission strength by angle", loc="left")
            cone.set_title("有多少功率集中在角锥内" if cjk else "How much power lies in the cone", loc="left")
        for j, ax in enumerate((geom, angular, cone)):
            letter(ax, chr(ord("a") + i * 3 + j), x=-0.18 if j == 0 else -0.10)
    text = (
        "浅色带为不同方位角的取值范围。曲线为加密结果，小圆点为粗采样。"
        if cjk else "Shading spans the azimuth minimum to maximum. Lines show refined data where available; points show coarse samples."
    ) if refined else (
        "浅色带为不同方位角的取值范围。曲线由粗采样结果连线。"
        if cjk else "Shading spans the azimuth minimum to maximum. Lines connect coarse samples."
    )
    fig.text(0.055, 0.041, text, fontsize=10, color=MUTED)
    save_figure(fig, out / "angular_response_pairs")


def label_points(ax, points: list[tuple[float, float, str, str]], marker_positions) -> None:
    """Place short labels near points while keeping labels inside their axes."""
    ax.figure.canvas.draw()
    to_pixel = ax.transData.transform
    boundary = ax.get_window_extent()
    offsets = [(6, 7), (6, -9), (-6, 7), (-6, -9), (9, 15), (-9, 15),
               (11, -17), (-11, -17), (15, 0), (-15, 0), (10, 23), (-10, -24),
               (25, 18), (-25, -18)]
    scale = ax.figure.dpi / 72
    clearance = 3.5 * scale
    occupied = [
        (px - clearance, py - clearance, px + clearance, py + clearance)
        for px, py in to_pixel(marker_positions)
    ]
    for x, y, name, color in sorted(points, key=lambda point: (-point[1], point[0])):
        px, py = to_pixel((x, y))
        for dx, dy in offsets:
            horizontal = "left" if dx >= 0 else "right"
            width, height = (4.4 * len(name) + 2) * scale, 9 * scale
            x0 = px + dx * scale - (width if dx < 0 else 0)
            y0 = py + dy * scale - height / 2
            x1, y1 = x0 + width, y0 + height
            if not (boundary.x0 + 2 <= x0 < x1 <= boundary.x1 - 2
                    and boundary.y0 + 2 <= y0 < y1 <= boundary.y1 - 2):
                continue
            if any(x0 < q[2] and x1 > q[0] and y0 < q[3] and y1 > q[1] for q in occupied):
                continue
            occupied.append((x0, y0, x1, y1))
            arrow = {"arrowstyle": "-", "lw": 0.45, "color": "#A9B8C5"} if abs(dy) > 12 or abs(dx) > 12 else None
            ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                        fontsize=8, color=color, ha=horizontal, va="center", arrowprops=arrow)
            break


def tradeoff_plot(coarse: dict, refined: dict, out: Path, cjk: bool) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.76, bottom=0.20, wspace=0.24)
    fig.text(0.09, 0.943, "角锥内有多强、能集中多少" if cjk else "Emission strength and angular concentration",
             size=20, color=NAVY, weight="bold")
    scope = (
        "M0 理想模型。实心点为粗采样，空心点为加密复算。"
        if cjk else "M0 ideal rays. Filled markers show coarse estimates and hollow markers show refinement."
    )
    fig.text(0.09, 0.891, scope, size=11, color=MUTED)
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none", color=FAMILIES["A0"]),
        plt.Line2D([0], [0], marker="s", linestyle="none", color=FAMILIES["A3"]),
        plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white",
                   markeredgecolor=NAVY, color=NAVY),
        plt.Line2D([0], [0], linestyle=":", color="#7E8E9C"),
    ]
    labels = ["A0 曲壁", "A3 直壁", "加密复算", "朗伯面"] if cjk else ["A0 curved", "A3 straight", "Refined", "Lambertian"]
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.97, 0.887), ncol=4)
    for i, (ax, angle) in enumerate(zip(axes, CONES)):
        setup_plot(ax)
        key = str(angle)
        points = []
        xs, ys = [], []
        reference = np.sin(np.deg2rad(angle)) ** 2
        ax.axvline(reference, color="#7E8E9C", linestyle=":", lw=1.5)
        for name, record in coarse.items():
            family = record["shape"]["family"]
            color = FAMILIES[family]
            marker = "o" if family == "A0" else "s"
            data = record["cones"][key]
            x, y = data["power_fraction"], data["mean_ideal_emissivity_inside"]
            if x is None:
                continue
            xs.append(x)
            ys.append(y)
            ax.scatter(x, y, color=color, marker=marker, s=40, zorder=3,
                       edgecolor="white", linewidth=0.65, alpha=0.95)
            label = name.split("_", 1)[1]
            label_x, label_y = x, y
            if name in refined:
                fine = refined[name]["cones"][key]
                xn, yn = fine["power_fraction"], fine["mean_ideal_emissivity_inside"]
                if xn is not None:
                    xs.append(xn)
                    ys.append(yn)
                    ax.plot([x, xn], [y, yn], lw=1, color=color, alpha=0.55)
                    ax.scatter(xn, yn, facecolor="white", edgecolor=color, marker=marker,
                               s=31, lw=1.25, zorder=4)
                    label_x, label_y = xn, yn
            if label in (*PAIR_IDS, "G05"):
                points.append((label_x, label_y, label, color))
        xmax = max(xs, default=1.0)
        ymin, ymax = min(ys, default=0), max(ys, default=1)
        ax.set_xlim(0, min(1.045, xmax + 0.07))
        ax.set_ylim(max(0, ymin - 0.08), min(1.03, ymax + 0.12))
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.set_xlabel("锥内功率占半球比例" if cjk else "Fraction of hemispheric power")
        ax.set_ylabel("锥内平均理想发射率" if cjk else "Cone-average ideal emissivity")
        ax.set_title(f"半角 {angle}°" if cjk else f"{angle}° cone half-angle", loc="left")
        letter(ax, chr(ord("a") + i))
        label_points(ax, points, list(zip(xs, ys)))
    footer = (
        "越靠上表示锥内发射率越高，越靠右表示功率越集中。图中标注三组复算结构和候选 G05。"
        if cjk else "Higher means stronger cone emission. Farther right means greater concentration. Labels mark the three refined pairs and G05."
    )
    fig.text(0.09, 0.067, footer, fontsize=10, color=MUTED)
    save_figure(fig, out / "cone_tradeoff")


def azimuth_plot(geometries: dict, out: Path, cjk: bool) -> None:
    results = [geometries.get(f"{family}_G06") for family in FAMILIES]
    if any(result is None or result.azimuth is None for result in results):
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 6), subplot_kw={"projection": "polar"})
    fig.subplots_adjust(left=0.055, right=0.855, top=0.78, bottom=0.16, wspace=0.30)
    fig.text(0.065, 0.945, "G06 的全方位角响应" if cjk else "Full azimuth response for G06",
             size=20, color=NAVY, weight="bold")
    fig.text(0.065, 0.885, "M0 理想模型  ·  圆心对应法向，半径表示极角" if cjk else
             "M0 ideal rays  ·  Center at the surface normal, radius represents polar angle.", size=11, color=MUTED)
    for ax, result in zip(axes, results):
        step = 360 / len(result.azimuth)
        phi_edges = np.arange(len(result.azimuth) + 1) * step - step / 2
        theta_edges = np.concatenate(([0], (result.theta[:-1] + result.theta[1:]) / 2, [90]))
        ax.grid(False)
        mesh = ax.pcolormesh(np.deg2rad(phi_edges), theta_edges, result.azimuth_values,
                             cmap="cividis", vmin=0, vmax=1, shading="flat")
        ax.set_ylim(0, 90)
        ax.set_rgrids([15, 30, 60], ["15°", "30°", "60°"], angle=45, color="white", fontsize=8)
        ax.set_xticks(np.deg2rad(np.arange(0, 360, 60)))
        ax.grid(True, lw=0.55, alpha=0.6, color="white")
        title = {"A0": "A0 曲壁", "A3": "A3 直壁"}[result.family] if cjk else result.name
        ax.set_title(title, pad=19, color=FAMILIES[result.family])
        ax.spines["polar"].set_color("#AFBEC8")
    cax = fig.add_axes([0.89, 0.21, 0.018, 0.52])
    colorbar = fig.colorbar(mesh, cax=cax, ticks=np.linspace(0, 1, 6), format=PercentFormatter(1, decimals=0))
    colorbar.set_label("理想方向发射率" if cjk else "Ideal directional emissivity", labelpad=8)
    fig.text(0.065, 0.05, "两图使用同一色标。色块按采样方向的发射率着色。" if cjk else
             "Shared color scale. Cells use the emissivity at sampled directions.", fontsize=10, color=MUTED)
    save_figure(fig, out / "g06_azimuth_maps")


def main() -> None:
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", type=Path, default=directory / "fetched/coarse_summary.json")
    parser.add_argument("--refine", type=Path)
    parser.add_argument("--out", type=Path, default=directory / "figures")
    args = parser.parse_args()
    coarse_meta, coarse = read_phase(args.coarse, "coarse")
    refine_path = args.refine or directory / "fetched/refine_summary.json"
    refined_meta, refined = read_phase(refine_path, "refine") if refine_path.is_file() else ({}, {})
    if args.refine and not refine_path.is_file():
        parser.error(f"Missing refinement summary: {refine_path}")
    check_refinement(coarse, refined)
    coarse_metrics = {name: geometry.metrics() for name, geometry in coarse.items()}
    refined_metrics = {name: geometry.metrics() for name, geometry in refined.items()}
    comparison = phase_comparison(coarse_metrics, refined_metrics) if refined else {}
    result = {
        "scope": "M0 ideal ray transport with a black floor and mirror walls, normalized by unit-cell area",
        "quadrature": "Azimuth mean followed by piecewise linear integration in u=sin²(theta); hold final angle value to 90 degrees",
        "unsampled_horizon_tail_weight": max(
            1 - float(np.sin(np.deg2rad(geometry.theta[-1])) ** 2)
            for geometry in (*coarse.values(), *refined.values())
        ),
        "floor_area_ratio_diagnostic": "Reciprocal ideal black floor area is a hemispheric emission upper bound; difference includes numerical quadrature and Monte Carlo effects and is not a pass gate",
        "comparison_scope": "Differences combine angular densification and changes in rays per direction",
        "coarse_source": str(args.coarse.resolve()),
        "refined_source": str(refine_path.resolve()) if refined else None,
        "coarse_conditions": coarse_meta["completed_conditions"],
        "refined_conditions": refined_meta.get("completed_conditions", 0),
        "coarse": coarse_metrics,
        "refined": refined_metrics,
        "phase_comparison": comparison,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "cone_metrics.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    cjk = setup_style()
    paired_plot(coarse, refined, args.out, cjk)
    tradeoff_plot(coarse_metrics, refined_metrics, args.out, cjk)
    azimuth_plot({**coarse, **refined}, args.out, cjk)
    print(json.dumps({"coarse_shapes": len(coarse), "refined_shapes": len(refined), "output": str(args.out.resolve())}))


if __name__ == "__main__":
    main()
