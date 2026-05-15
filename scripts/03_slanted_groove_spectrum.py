"""03_slanted_groove_spectrum.py — Ti 基底 2D 周期斜槽光谱与倾角扫描

研究目的
--------
在 02 (对称矩形槽) 的基础上引入 **非对称几何** —— 单胞内一个倾斜四边形槽
(顶宽 ≠ 底宽，或顶/底中心在 x 方向错位)。研究倾角 α 对 5–15 μm 光谱响应
(R, A, ε_proxy) 的影响，为后续 “角分辨发射 / 定向发射” 建模打下基础。

物理任务
--------
1. Ti 基底（不透明近似，T ≈ 0）+ 顶部周期性斜槽 (period_um);
2. 几何参数：top_width_um, bottom_width_um, depth_um, tilt_angle_deg；
3. 几何工厂位于 `src/geometry.build_slanted_groove_geometry`：
     - tilt_angle = 0 + top == bottom → 对称矩形（退化）；
     - tilt_angle ≠ 0 → 顶底中心错位 depth·tan(α)，整体非对称；
   这与“把对称矩形整体旋转”是**完全不同**的：本函数在固定单胞内改变轮廓，
   保留了表面周期结构的物理含义，而坐标旋转只是参考系变换、不改光谱；
4. 仿真核心复用 `src/simulation.run_periodic_2d_metal_spectrum`，因此
   监视器布局 / 源 / PML 与 01-02 完全一致；
5. 单偏振 Ez (TE w.r.t. 光栅，E 沿槽轴 z) —— 与 02 保持一致，便于直接对照；
6. 批量比较 α ∈ {0°, 10°, 20°, 30°}。

为什么 2D 等效模型够用作 “先看趋势”
------------------------------------
- 真实飞秒激光斜烧蚀产生的是 3D 倾斜圆/椭圆截面孔，有限深度、侧壁圆滑、
  常带 TiO/TiO₂ 氧化层；
- 当前 2D 模型保留了 “斜壁 + 周期性 + 一个特征宽度” 的最少必要要素，
  足以在正入射下看出非对称几何如何改变 R/A 谱的形状与积分发射率；
- **绝对量** 不可直接与 3D 实验对比；**趋势** (随 α 单调或非单调变化、
  共振峰位置漂移) 可用作机器学习数据的预探针。

后续扩展（按复杂度递增）
------------------------
1. 引入圆角：把 4 顶点 prism 换成 8 顶点（顶/底各两个小圆弧），近似真实
   烧蚀槽口的曲率半径；
2. 加氧化层：在 prism 槽壁内嵌一层 ε_TiO₂ ≈ 6 (中红外) 的薄 Medium；
3. 升 3D 仿真：把 axis=z 的 prism 换成沿斜轴方向的 mp.Prism / mp.Cone，
   配合 3D 周期边界（x、z 双向 Bloch）；
4. 角分辨：在本脚本基础上扫描入射 / 出射 θ，得到 ε(λ, θ)。

输入 / 输出
-----------
输入：仅命令行参数 (--period, --top_width, --bottom_width, --depth,
       --tilt_angles ...)。
输出：
  - results/tables/slanted_groove_spectra.csv     # 长表 (α, λ, R, T, A, ε)
  - results/figures/slanted_groove_geometry.png    # 4 个 α 下的几何轮廓
  - results/figures/slanted_groove_angle_comparison.png  # 谱比较
  - logs/03_slanted_groove_spectrum_<时间戳>.log

物理假设与近似（必须随结果一同记录）
----------------------------------
1. 不透明基底 → T ≈ 0，emissivity_proxy ≡ absorptance；
2. 2D 等效（槽沿 z 无穷长），非真实 3D 斜孔；
3. 单偏振 Ez；
4. Ti Rakić 1998 在 13–15 μm 为 Drude 尾部外推；
5. 室温、无氧化层、表面理想光滑、侧壁无圆角。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import meep as mp
import numpy as np
import pandas as pd

from src.geometry import (
    build_rectangular_groove_geometry,
    build_slanted_groove_geometry,
    slanted_groove_vertices,
)
from src.io_utils import (
    make_run_tag,
    project_path,
    save_figure,
    save_spectrum_csv,
    setup_logger,
)
from src.materials import TI_RAKIC_VALID_LAMBDA_UM, get_ti_medium
from src.postprocess import energy_conservation_check
from src.simulation import run_periodic_2d_metal_spectrum


# ---------------------------------------------------------------------------
# 默认参数（全部可被 CLI 覆盖）
# ---------------------------------------------------------------------------

DEFAULTS = dict(
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    tilt_angles_deg=[0.0, 10.0, 20.0, 30.0],
    wavelength_min_um=5.0,
    wavelength_max_um=15.0,
    resolution=32,             # Ti Drude–Lorentz 稳定下限
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    nfreq=121,
    decay_db=40.0,
)

# Ti 稳定下限 (memory: feedback_ti_resolution)
TI_MAX_LORENTZ_FREQ = 15.67
MIN_SAFE_RES = int(2 * TI_MAX_LORENTZ_FREQ) + 1  # 32


# ---------------------------------------------------------------------------
# 对外 API
# ---------------------------------------------------------------------------

def simulate_slanted_groove(
    period_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    wavelength_min_um: float = DEFAULTS["wavelength_min_um"],
    wavelength_max_um: float = DEFAULTS["wavelength_max_um"],
    resolution: int = DEFAULTS["resolution"],
    *,
    pml_thickness_um: float = DEFAULTS["pml_thickness_um"],
    substrate_thickness_um: float = DEFAULTS["substrate_thickness_um"],
    air_buffer_um: float = DEFAULTS["air_buffer_um"],
    nfreq: int = DEFAULTS["nfreq"],
    decay_db: float = DEFAULTS["decay_db"],
    logger=None,
) -> pd.DataFrame:
    """2D 周期斜槽 Ti 表面正入射光谱仿真。

    Returns
    -------
    pandas.DataFrame
        Columns: wavelength_um, reflectance, transmittance, absorptance,
        emissivity_proxy。emissivity_proxy ≡ absorptance (不透明基底)。

    Raises
    ------
    ValueError
        参数非法（含 resolution < 32、宽度越出单胞、顶点越界等）。
    """
    if resolution < MIN_SAFE_RES:
        raise ValueError(
            f"resolution={resolution} 低于 Ti Drude–Lorentz 稳定下限 {MIN_SAFE_RES}。"
        )

    medium_ti = get_ti_medium(lambda_min_um=wavelength_min_um,
                              lambda_max_um=wavelength_max_um)

    def factory(y_surface_um: float, substrate_thickness: float) -> list:
        return build_slanted_groove_geometry(
            period_x_um=period_um,
            top_width_um=top_width_um,
            bottom_width_um=bottom_width_um,
            depth_um=depth_um,
            tilt_angle_deg=tilt_angle_deg,
            substrate_thickness_um=substrate_thickness,
            y_surface=y_surface_um,
            medium_substrate=medium_ti,
        )

    result = run_periodic_2d_metal_spectrum(
        geometry_factory=factory,
        period_um=period_um,
        wavelength_min_um=wavelength_min_um,
        wavelength_max_um=wavelength_max_um,
        resolution=resolution,
        pml_thickness_um=pml_thickness_um,
        substrate_thickness_um=substrate_thickness_um,
        air_buffer_um=air_buffer_um,
        nfreq=nfreq,
        decay_db=decay_db,
        source_component="Ez",
        logger=logger,
    )

    df = pd.DataFrame({
        "wavelength_um": result["wavelength_um"],
        "reflectance": result["R"],
        "transmittance": result["T"],
        "absorptance": result["A"],
        "emissivity_proxy": result["A"],  # 不透明基底 + Kirchhoff
    })
    df.attrs.update(
        period_um=period_um,
        top_width_um=top_width_um,
        bottom_width_um=bottom_width_um,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        resolution=resolution,
        walltime_s=result["walltime_s"],
        geometry_y=result["geometry_y"],
    )
    return df


def simulate_rectangular_control(
    period_um: float,
    groove_width_um: float,
    groove_depth_um: float,
    wavelength_min_um: float,
    wavelength_max_um: float,
    resolution: int,
    *,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    air_buffer_um: float,
    nfreq: int,
    decay_db: float,
    logger=None,
) -> pd.DataFrame:
    """退化对照：用矩形槽几何工厂跑一次 (用于验证 tilt=0 + top==bottom 的斜槽
    给出相同光谱)。"""
    medium_ti = get_ti_medium(lambda_min_um=wavelength_min_um,
                              lambda_max_um=wavelength_max_um)

    def factory(y_surface_um: float, substrate_thickness: float) -> list:
        return build_rectangular_groove_geometry(
            period_x_um=period_um,
            groove_width_um=groove_width_um,
            groove_depth_um=groove_depth_um,
            substrate_thickness_um=substrate_thickness,
            y_surface=y_surface_um,
            medium_substrate=medium_ti,
        )

    result = run_periodic_2d_metal_spectrum(
        geometry_factory=factory,
        period_um=period_um,
        wavelength_min_um=wavelength_min_um,
        wavelength_max_um=wavelength_max_um,
        resolution=resolution,
        pml_thickness_um=pml_thickness_um,
        substrate_thickness_um=substrate_thickness_um,
        air_buffer_um=air_buffer_um,
        nfreq=nfreq,
        decay_db=decay_db,
        source_component="Ez",
        logger=logger,
    )
    return pd.DataFrame({
        "wavelength_um": result["wavelength_um"],
        "reflectance": result["R"],
        "transmittance": result["T"],
        "absorptance": result["A"],
        "emissivity_proxy": result["A"],
    })


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------

def geometry_self_checks(
    *,
    period_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angles_deg: list[float],
    logger,
) -> dict:
    """对一组倾角做几何自检：
       (1) 所有顶点 x ∈ [-P/2, P/2]；
       (2) 不同 α 的轮廓 x 坐标确实不同 (排除 “只是坐标旋转” 的伪非对称)。
    """
    half_P = period_um / 2.0
    vertices_by_alpha: dict[float, list[tuple[float, float]]] = {}
    for alpha in tilt_angles_deg:
        v = slanted_groove_vertices(
            top_width_um=top_width_um,
            bottom_width_um=bottom_width_um,
            depth_um=depth_um,
            tilt_angle_deg=alpha,
        )
        for (x, _y) in v:
            if abs(x) > half_P + 1e-9:
                raise ValueError(
                    f"α={alpha}° 顶点 x={x:.3f} 越出单胞 [-{half_P}, {half_P}] μm"
                )
        vertices_by_alpha[alpha] = v

    # 检查不同 α 产生不同轮廓（仅当 α≠0 时应与 α=0 有差异）
    base = np.array(vertices_by_alpha[tilt_angles_deg[0]])
    distinct_pairs = 0
    for alpha in tilt_angles_deg[1:]:
        diff = np.max(np.abs(np.array(vertices_by_alpha[alpha]) - base))
        if alpha != tilt_angles_deg[0] and diff < 1e-9:
            raise RuntimeError(
                f"α={alpha}° 的顶点与 α={tilt_angles_deg[0]}° 完全一致；"
                f"几何并未随倾角改变。"
            )
        distinct_pairs += int(diff > 1e-9)

    logger.info("geometry self-check: 所有顶点在单胞内；%d/%d 个 α 与基准产生几何差异。",
                distinct_pairs, len(tilt_angles_deg) - 1)
    return vertices_by_alpha


def spectrum_self_checks(
    spectra_by_alpha: dict[float, pd.DataFrame],
    *,
    logger,
    energy_tol: float = 0.05,
    T_opaque_tol: float = 1e-2,
    range_tol: float = 1e-3,
) -> dict:
    """能量守恒 + 不透明性 + 不同 α 是否真的改变光谱。"""
    summary = {}
    wl_ref = None
    ti_valid_hi = TI_RAKIC_VALID_LAMBDA_UM[1]
    for alpha, df in spectra_by_alpha.items():
        wl = df["wavelength_um"].to_numpy()
        if wl_ref is None:
            wl_ref = wl
        elif not np.allclose(wl, wl_ref):
            raise RuntimeError(f"α={alpha} 的波长轴与其他不一致。")
        if not np.all(np.diff(wl) > 0):
            raise RuntimeError(f"α={alpha} 的波长未升序。")

        rta = {"R": df["reflectance"].to_numpy(),
               "T": df["transmittance"].to_numpy(),
               "A": df["absorptance"].to_numpy()}
        rep = energy_conservation_check(rta, atol=energy_tol)
        Tmax = float(np.max(np.abs(rta["T"])))
        range_violations = {
            "reflectance": int(np.sum((rta["R"] < -range_tol) | (rta["R"] > 1 + range_tol))),
            "transmittance": int(np.sum((rta["T"] < -range_tol) | (rta["T"] > 1 + range_tol))),
            "absorptance": int(np.sum((rta["A"] < -range_tol) | (rta["A"] > 1 + range_tol))),
        }
        valid_8_13 = (wl >= 8) & (wl <= min(13, ti_valid_hi))
        summary[alpha] = dict(
            max_energy_err=rep["max_abs_dev"],
            mean_energy_err=rep["mean_abs_dev"],
            Tmax=Tmax,
            range_violations=range_violations,
            avg_A_8_13=float(np.mean(rta["A"][(wl >= 8) & (wl <= 13)])),
            avg_A_8_valid=float(np.mean(rta["A"][valid_8_13])),
        )
        if Tmax > T_opaque_tol:
            logger.warning("α=%g° Tmax=%.3e 超过不透明阈值。", alpha, T_opaque_tol)
        if any(range_violations.values()):
            logger.warning(
                "α=%g° 存在 R/T/A 超出物理范围的频点: %s。"
                "常见原因是材料外推波段、源边缘信噪比或收敛不足。",
                alpha, range_violations,
            )
        logger.info(
            "self-check α=%g°  max|R+T+A-1|=%.2e  Tmax=%.2e  "
            "range=%s  <A>_{8-13}=%.3f  <A>_{8-%.3f valid}=%.3f",
            alpha, rep["max_abs_dev"], Tmax, range_violations,
            summary[alpha]["avg_A_8_13"], min(13, ti_valid_hi),
            summary[alpha]["avg_A_8_valid"],
        )

    if wl_ref is not None and float(np.max(wl_ref)) > ti_valid_hi:
        logger.warning(
            "当前最大波长 %.3f μm 超出 Ti Rakić 模型标定上限 %.3f μm；"
            "建议把定量结论限制在有效波段，或替换为覆盖目标波段的光学常数。",
            float(np.max(wl_ref)), ti_valid_hi,
        )

    # 不同 α 的谱必须有非平凡差异（否则 prism 没生效）
    alphas_sorted = sorted(spectra_by_alpha.keys())
    if len(alphas_sorted) >= 2:
        base = spectra_by_alpha[alphas_sorted[0]]["reflectance"].to_numpy()
        max_dR_across = 0.0
        for alpha in alphas_sorted[1:]:
            diff = float(np.max(np.abs(
                spectra_by_alpha[alpha]["reflectance"].to_numpy() - base
            )))
            max_dR_across = max(max_dR_across, diff)
        summary["max_dR_across_alpha"] = max_dR_across
        if max_dR_across < 1e-3:
            logger.warning(
                "所有 α 的 R 谱差异 <1e-3，疑似几何未真正变化或网格太粗。"
            )
        else:
            logger.info("spectrum self-check: max|ΔR| across α = %.3e", max_dR_across)
    return summary


def regression_check_alpha_zero(
    df_alpha0: pd.DataFrame,
    df_rectangular: pd.DataFrame,
    *,
    logger,
    tol: float = 5e-2,
) -> dict:
    """tilt=0 且 top==bottom 的斜槽应与对称矩形槽给出相同光谱。

    两者几何在 Meep 中分别用 mp.Prism 和 mp.Block 表示，离散化栅格差异会
    导致 R 谱有 O(1/res) 级别的小偏差，所以 tol 默认 5%。
    """
    wl = df_alpha0["wavelength_um"].to_numpy()
    if not np.allclose(wl, df_rectangular["wavelength_um"].to_numpy()):
        raise RuntimeError("回归对照的波长轴不一致。")
    dR = df_alpha0["reflectance"].to_numpy() - df_rectangular["reflectance"].to_numpy()
    dA = df_alpha0["absorptance"].to_numpy() - df_rectangular["absorptance"].to_numpy()
    max_dR = float(np.max(np.abs(dR)))
    max_dA = float(np.max(np.abs(dA)))
    logger.info(
        "regression α=0° vs rectangular: max|ΔR|=%.3e  max|ΔA|=%.3e  (tol=%.2g)",
        max_dR, max_dA, tol,
    )
    if max_dR > tol or max_dA > tol:
        logger.warning(
            "α=0° 与矩形槽偏差超过容差 %.2g：max|ΔR|=%.3e max|ΔA|=%.3e；"
            "可能 prism 与 block 离散化在斜壁上栅格不同；可考虑提高 resolution。",
            tol, max_dR, max_dA,
        )
    return dict(max_dR=max_dR, max_dA=max_dA, tol=tol)


# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------

def plot_geometry_overview(
    *,
    period_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    substrate_thickness_um: float,
    tilt_angles_deg: list[float],
    out_path: Path,
) -> Path:
    """画 4 个 α 下的几何轮廓 (xy 截面)。"""
    n = len(tilt_angles_deg)
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 3.4), sharey=True)
    if n == 1:
        axes = [axes]
    half_P = period_um / 2.0
    y_surface = 0.0
    y_sub_bot = -substrate_thickness_um
    for ax, alpha in zip(axes, tilt_angles_deg):
        # 基底矩形 (Ti)
        ax.add_patch(Rectangle(
            (-half_P, y_sub_bot), period_um, substrate_thickness_um,
            facecolor="#c9c9c9", edgecolor="k", lw=0.8,
        ))
        # 槽 (air)
        verts = slanted_groove_vertices(
            top_width_um=top_width_um,
            bottom_width_um=bottom_width_um,
            depth_um=depth_um,
            tilt_angle_deg=alpha,
            y_surface=y_surface,
        )
        ax.add_patch(Polygon(verts, closed=True,
                              facecolor="white", edgecolor="C3", lw=1.6))
        # 单胞边界
        ax.axvline(-half_P, color="C0", ls=":", lw=0.8)
        ax.axvline(+half_P, color="C0", ls=":", lw=0.8)
        ax.axhline(y_surface, color="0.4", ls=":", lw=0.5)

        ax.set_xlim(-half_P * 1.15, half_P * 1.15)
        ax.set_ylim(y_sub_bot - 0.2, depth_um * 0.6)
        ax.set_aspect("equal")
        ax.set_title(f"α = {alpha:g}°")
        ax.set_xlabel("x (μm)")
        ax.grid(True, ls=":", alpha=0.4)
    axes[0].set_ylabel("y (μm)")
    fig.suptitle(
        f"Slanted groove geometry: P={period_um:g}, top_w={top_width_um:g}, "
        f"bot_w={bottom_width_um:g}, D={depth_um:g} μm",
        fontsize=11,
    )
    fig.tight_layout()
    return save_figure(fig, out_path)


def plot_angle_comparison(
    spectra_by_alpha: dict[float, pd.DataFrame],
    out_path: Path,
    *,
    title: str = "",
) -> Path:
    """同一坐标轴上画所有 α 的 R(λ) 和 A(λ)。"""
    fig, (ax_R, ax_A) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True)
    cmap = plt.get_cmap("viridis")
    alphas_sorted = sorted(spectra_by_alpha.keys())
    n = len(alphas_sorted)
    for i, alpha in enumerate(alphas_sorted):
        df = spectra_by_alpha[alpha]
        color = cmap(0.15 + 0.7 * (i / max(n - 1, 1)))
        ax_R.plot(df["wavelength_um"], df["reflectance"], lw=1.7,
                  color=color, label=f"α = {alpha:g}°")
        ax_A.plot(df["wavelength_um"], df["absorptance"], lw=1.7,
                  color=color, label=f"α = {alpha:g}°  (ε_proxy)")
    for ax in (ax_R, ax_A):
        ax.axvspan(8, 13, color="orange", alpha=0.08, label="8–13 μm window")
        ax.grid(True, ls=":", alpha=0.6)
        ax.legend(loc="best", framealpha=0.9, fontsize=9)
    ax_R.set_ylabel("Reflectance R")
    ax_R.set_ylim(0, 1.05)
    ax_R.set_title(title or "Slanted groove: angle sweep, 2D, Ez")
    ax_A.set_xlabel("Wavelength λ (μm)")
    ax_A.set_ylabel("Absorptance A ≈ ε (opaque)")
    return save_figure(fig, out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--period", type=float, default=DEFAULTS["period_um"],
                   help="x 方向周期 (μm)")
    p.add_argument("--top_width", type=float, default=DEFAULTS["top_width_um"],
                   help="槽顶宽 (μm)，必须 < period")
    p.add_argument("--bottom_width", type=float, default=DEFAULTS["bottom_width_um"],
                   help="槽底宽 (μm)，必须 < period")
    p.add_argument("--depth", type=float, default=DEFAULTS["depth_um"],
                   help="槽深 (μm)，必须 > 0 且 < substrate")
    p.add_argument("--tilt_angles", type=float, nargs="+",
                   default=DEFAULTS["tilt_angles_deg"],
                   help="倾角列表 (°)，默认 [0, 10, 20, 30]")
    p.add_argument("--wavelength_min_um", type=float,
                   default=DEFAULTS["wavelength_min_um"])
    p.add_argument("--wavelength_max_um", type=float,
                   default=DEFAULTS["wavelength_max_um"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--pml_thickness_um", type=float,
                   default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float,
                   default=DEFAULTS["air_buffer_um"])
    p.add_argument("--nfreq", type=int, default=DEFAULTS["nfreq"])
    p.add_argument("--no_regression", action="store_true",
                   help="跳过 α=0° vs 对称矩形槽 的回归对照")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    run_tag = make_run_tag("03_slanted_groove_spectrum")
    logger = setup_logger(run_tag)
    logger.info("=== 03_slanted_groove_spectrum ===")
    logger.info("args = %s", vars(args))

    mp.verbosity(1)

    # ---- 0) 几何自检（无需 Meep） ----
    geometry_self_checks(
        period_um=args.period,
        top_width_um=args.top_width,
        bottom_width_um=args.bottom_width,
        depth_um=args.depth,
        tilt_angles_deg=args.tilt_angles,
        logger=logger,
    )

    # ---- 1) 几何可视化 ----
    geom_fig = plot_geometry_overview(
        period_um=args.period,
        top_width_um=args.top_width,
        bottom_width_um=args.bottom_width,
        depth_um=args.depth,
        substrate_thickness_um=args.substrate_thickness_um,
        tilt_angles_deg=sorted(args.tilt_angles),
        out_path=project_path("results", "figures", "slanted_groove_geometry.png"),
    )
    logger.info("geometry PNG → %s", geom_fig)

    # ---- 2) 批量倾角扫描 ----
    spectra_by_alpha: dict[float, pd.DataFrame] = {}
    for alpha in args.tilt_angles:
        logger.info(">>> simulate α = %g°", alpha)
        spectra_by_alpha[alpha] = simulate_slanted_groove(
            period_um=args.period,
            top_width_um=args.top_width,
            bottom_width_um=args.bottom_width,
            depth_um=args.depth,
            tilt_angle_deg=alpha,
            wavelength_min_um=args.wavelength_min_um,
            wavelength_max_um=args.wavelength_max_um,
            resolution=args.resolution,
            pml_thickness_um=args.pml_thickness_um,
            substrate_thickness_um=args.substrate_thickness_um,
            air_buffer_um=args.air_buffer_um,
            nfreq=args.nfreq,
            logger=logger,
        )

    # ---- 3) 光谱自检 ----
    spec_summary = spectrum_self_checks(spectra_by_alpha, logger=logger)

    # ---- 4) α=0 回归 vs 对称矩形槽（同 res/网格） ----
    reg_summary = None
    if (not args.no_regression
            and 0.0 in args.tilt_angles
            and abs(args.top_width - args.bottom_width) < 1e-9):
        logger.info(">>> regression: α=0° (prism) vs rectangular (block)")
        df_rect = simulate_rectangular_control(
            period_um=args.period,
            groove_width_um=args.top_width,
            groove_depth_um=args.depth,
            wavelength_min_um=args.wavelength_min_um,
            wavelength_max_um=args.wavelength_max_um,
            resolution=args.resolution,
            pml_thickness_um=args.pml_thickness_um,
            substrate_thickness_um=args.substrate_thickness_um,
            air_buffer_um=args.air_buffer_um,
            nfreq=args.nfreq,
            decay_db=DEFAULTS["decay_db"],
            logger=logger,
        )
        reg_summary = regression_check_alpha_zero(
            spectra_by_alpha[0.0], df_rect, logger=logger,
        )
    elif not args.no_regression:
        logger.info("跳过回归对照：tilt_angles 中无 0° 或 top_width != bottom_width。")

    # ---- 5) CSV：长表 ----
    long_rows = []
    for alpha, df in spectra_by_alpha.items():
        for _, row in df.iterrows():
            long_rows.append({
                "period_um": args.period,
                "top_width_um": args.top_width,
                "bottom_width_um": args.bottom_width,
                "depth_um": args.depth,
                "tilt_angle_deg": alpha,
                "wavelength_um": row["wavelength_um"],
                "reflectance": row["reflectance"],
                "transmittance": row["transmittance"],
                "absorptance": row["absorptance"],
                "emissivity_proxy": row["emissivity_proxy"],
            })
    long_df = pd.DataFrame(long_rows)
    csv_out = project_path("results", "tables", "slanted_groove_spectra.csv")
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(csv_out, index=False)
    logger.info("CSV → %s  (rows=%d)", csv_out, len(long_df))

    # ---- 6) 谱对比图 ----
    title = (f"Slanted groove sweep: P={args.period:g}, top={args.top_width:g}, "
             f"bot={args.bottom_width:g}, D={args.depth:g} μm  "
             f"(res={args.resolution}/μm, Ez)")
    fig_out = plot_angle_comparison(
        spectra_by_alpha,
        project_path("results", "figures", "slanted_groove_angle_comparison.png"),
        title=title,
    )
    logger.info("comparison PNG → %s", fig_out)

    # ---- 7) 总结 ----
    total_walltime = sum(df.attrs["walltime_s"] for df in spectra_by_alpha.values())
    logger.info("=== done; total walltime = %.1f s ===", total_walltime)
    logger.info("spectrum summary: %s", spec_summary)
    if reg_summary is not None:
        logger.info("regression summary: %s", reg_summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
