"""02_periodic_groove_spectrum.py — Ti 基底 1D 周期矩形槽光谱正问题

研究目的
--------
建立最基本的 “微结构参数 → 光谱响应” 正问题函数 `simulate_periodic_groove`，
为后续批量参数扫描与机器学习数据生成奠定基础。本脚本同时跑：
  - 平板 Ti（无结构）作为基线；
  - 1D 周期矩形槽 Ti（period_um / groove_width_um / groove_depth_um）；
并把两者的 R/T/A/ε_proxy 谱画在同一张图上，直接看到结构带来的差异。

物理任务
--------
1. Ti 基底（不透明近似，T 应 ≈ 0）；
2. 2D 仿真，x 方向周期 = period_um，y 方向 PML；
3. 槽 = 在 Ti 上表面挖一个矩形空气块，宽 groove_width_um，深 groove_depth_um；
4. 波长 5–15 μm；
5. 单偏振（详见下面 “偏振约定”）。

偏振约定
--------
- 仿真平面为 xy；
- 槽周期方向沿 x，槽长度方向沿 z（out-of-plane，无穷长）；
- 入射沿 -y；
- 本脚本只跑 **Ez 偏振**（E 沿 z，即沿槽方向）。在光栅文献中这对应
  **TE / s-偏振 w.r.t. 光栅**（E 平行于槽线）。
- 选 Ez 的理由：与 01_flat_ti_benchmark.py 保持一致，便于直接对照基线。
- 局限：TM/p（E 垂直于槽线，Ex 分量）会激发 cavity/spoof-SPP 模式，与
  发射率增强关系更密切；待 01–02 流程稳定后再补一个 TM 版本，最终
  非偏振发射率取 TE/TM 平均。

为什么先做矩形槽，而不是直接做斜孔
----------------------------------
- 矩形槽是 1D 周期，几何只有 3 个参数（P, W, D），网格 + 周期边界都最简单，
  数值收敛性、能量守恒最容易验证；
- 解析极限清楚：W→0 或 D→0 退化为平板，等效介质理论在 W ≪ λ 时给出近似；
- 飞秒激光实际烧蚀的微孔阵列是 2D 周期且非矩形，几何复杂度高 1–2 个数量级；
  先把 1D 矩形槽这个 “最干净” 的物理模型跑通，后续 2D 锥孔/斜孔脚本可以
  直接复用本脚本的 simulate_* 框架，只换 geometry。

simulate_periodic_groove 未来如何被批量调用
-------------------------------------------
- 该函数签名 (period, width, depth, λ_range, resolution) → DataFrame，是纯函数，
  无全局状态，可直接被 itertools.product / joblib.Parallel 包起来扫参；
- 返回 DataFrame 而非 ndarray，便于直接 concat 出 (P, W, D, λ, R, T, A, ε) 的
  长表，作为 ML 训练集；
- 同一接口未来替换 geometry kind 即可扩展到 2D 孔阵 / 斜孔。

输入 / 输出
-----------
输入：仅命令行参数（不依赖外部 YAML）。
输出：
  - results/tables/periodic_groove_spectrum.csv        # 槽 + 平板对照
  - results/figures/periodic_groove_vs_flat_ti.png
  - logs/02_periodic_groove_spectrum_<时间戳>.log
DataFrame 列：
  wavelength_um, reflectance, transmittance, absorptance, emissivity_proxy

物理假设与近似（必须随结果一同记录）
----------------------------------
1. **不透明基底近似**：substrate_thickness ≫ skin depth，故 T ≈ 0，定义
   emissivity_proxy = absorptance（Kirchhoff 定律 ε = A，前提是局域热平衡）。
   脚本会打印 T 的最大值确认这一近似成立。
2. **2D 等效**：槽沿 z 方向无穷长；真实激光烧蚀结构是有限长 2D 阵列，
   端效应未建模。
3. **单偏振 (Ez)**：非偏振发射率应取 TE+TM 平均，本脚本只给 TE 部分。
4. **Ti Drude–Lorentz (Rakić 1998)** 标定范围 0.248–12.4 μm；13–15 μm 是
   Drude 尾部外推。
5. **室温、无氧化层、表面理想光滑**：实际飞秒激光改性表面有 TiO/TiO₂ 和
   亚波长粗糙度，本模型不含。
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# 让 `python scripts/02_...py` 能 import src/
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meep as mp
import numpy as np
import pandas as pd

from src.io_utils import (
    make_run_tag,
    project_path,
    save_figure,
    save_spectrum_csv,
    setup_logger,
)
from src.materials import (
    TI_RAKIC_VALID_LAMBDA_UM,
    freq_range_for_band,
    get_ti_medium,
    meep_freq_to_wavelength_um,
)
from src.postprocess import energy_conservation_check


# ---------------------------------------------------------------------------
# 默认参数
# ---------------------------------------------------------------------------

DEFAULTS = dict(
    wavelength_min_um=5.0,
    wavelength_max_um=15.0,
    resolution=32,             # Ti Drude–Lorentz 稳定下限，见 memory
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,  # 远超中红外 Ti 趋肤深度；同时给槽留足深度余量
    air_buffer_um=8.0,
    nfreq=201,
    decay_db=40.0,
)

# Ti 稳定性下限（参见 memory: feedback_ti_resolution）
TI_MAX_LORENTZ_FREQ = 15.67
MIN_SAFE_RES = int(2 * TI_MAX_LORENTZ_FREQ) + 1  # 32


# ---------------------------------------------------------------------------
# 几何 + Meep 仿真核心（私有）
# ---------------------------------------------------------------------------

def _build_groove_geometry(
    *,
    period_x_um: float,
    groove_width_um: float,
    groove_depth_um: float,
    substrate_thickness_um: float,
    y_surface: float,
    medium_ti,
) -> list:
    """构造 Ti 基底 + 顶部矩形槽（air 块覆盖）。

    Meep 中后加入的 geometry 覆盖先加入的，因此先放 Ti slab，再放 air block
    “挖” 出槽。groove_width_um=0 或 groove_depth_um=0 时不添加 air，等价于
    平板基线。
    """
    geom = [
        mp.Block(
            material=medium_ti,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_x_um, substrate_thickness_um, mp.inf),
        )
    ]
    if groove_width_um > 0 and groove_depth_um > 0:
        geom.append(
            mp.Block(
                material=mp.Medium(epsilon=1.0),  # 空气
                center=mp.Vector3(0, y_surface - groove_depth_um / 2.0, 0),
                size=mp.Vector3(groove_width_um, groove_depth_um, mp.inf),
            )
        )
    return geom


def _run_meep_sim(
    *,
    period_um: float,
    groove_width_um: float,
    groove_depth_um: float,
    wavelength_min_um: float,
    wavelength_max_um: float,
    resolution: int,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    air_buffer_um: float,
    nfreq: int,
    decay_db: float,
    logger=None,
) -> dict:
    """一次 “参考 run + 结构 run”，返回 R/T/A 与几何元数据。

    与 01_flat_ti_benchmark.py 中 `run_flat_ti_sim` 的布局完全一致，仅在
    geometry 上换成 “Ti slab + 顶部矩形槽”。布局沿 y 自上而下：
        [+y] top PML / 空气 buffer / 源 / 反射 monitor / 空气 /
             Ti 上表面 (y=y_surface) / [槽 air 块 嵌在最上层 depth μm] /
             Ti 体 / 透射 monitor / bottom PML  [-y]
    """
    log = logger.info if logger else print

    f_min, f_max, fcen = freq_range_for_band(wavelength_min_um, wavelength_max_um)
    df = f_max - f_min

    # cell 尺寸
    bottom_buffer_um = pml_thickness_um
    cell_y = (
        2 * pml_thickness_um
        + air_buffer_um
        + substrate_thickness_um
        + bottom_buffer_um
    )
    cell = mp.Vector3(period_um, cell_y, 0)

    # 关键 y 坐标
    y_top_edge = +0.5 * cell_y
    y_bottom_edge = -0.5 * cell_y
    y_top_pml_inner = y_top_edge - pml_thickness_um
    y_bottom_pml_inner = y_bottom_edge + pml_thickness_um
    y_surface = y_top_pml_inner - air_buffer_um
    y_substrate_bottom = y_surface - substrate_thickness_um
    y_src = y_top_pml_inner - 0.25 * air_buffer_um
    y_refl = y_surface + 0.5 * air_buffer_um
    y_trans = 0.5 * (y_substrate_bottom + y_bottom_pml_inner)

    if not (y_surface < y_refl < y_src < y_top_pml_inner):
        raise RuntimeError("源/反射 monitor 的 y 坐标顺序异常。")
    if not (y_bottom_pml_inner < y_trans < y_substrate_bottom):
        raise RuntimeError("透射 monitor 未位于 Ti 底面与 bottom PML 之间。")
    if groove_depth_um >= substrate_thickness_um:
        raise ValueError(
            f"groove_depth_um ({groove_depth_um}) 必须 < substrate_thickness_um "
            f"({substrate_thickness_um})，否则槽穿透基底。"
        )

    pml_layers = [mp.PML(thickness=pml_thickness_um, direction=mp.Y)]

    # 单偏振：Ez，详见文件顶部说明
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=df, is_integrated=True),
            component=mp.Ez,
            center=mp.Vector3(0, y_src, 0),
            size=mp.Vector3(period_um, 0, 0),
        )
    ]

    medium_ti = get_ti_medium(lambda_min_um=wavelength_min_um,
                              lambda_max_um=wavelength_max_um)
    geometry = _build_groove_geometry(
        period_x_um=period_um,
        groove_width_um=groove_width_um,
        groove_depth_um=groove_depth_um,
        substrate_thickness_um=substrate_thickness_um,
        y_surface=y_surface,
        medium_ti=medium_ti,
    )

    # --- 参考 run（空气） ---
    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=resolution,
        k_point=mp.Vector3(),
        geometry=[],
        dimensions=2,
    )
    refl_fr_ref = sim_ref.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0, y_refl, 0),
                      size=mp.Vector3(period_um, 0, 0)),
    )
    t0 = time.time()
    log("[ref run] res=%d, cell=(%.2f, %.2f) μm, nfreq=%d, fcen=%.4f, df=%.4f",
        resolution, cell.x, cell.y, nfreq, fcen, df)
    sim_ref.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20, c=mp.Ez,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-decay_db / 10.0),
        )
    )
    input_flux = np.array(mp.get_fluxes(refl_fr_ref))
    freqs = np.array(mp.get_flux_freqs(refl_fr_ref))
    ref_data = sim_ref.get_flux_data(refl_fr_ref)
    t_ref = time.time() - t0
    log("[ref run] done in %.2f s, |input_flux| ∈ [%.3e, %.3e]",
        t_ref, np.min(np.abs(input_flux)), np.max(np.abs(input_flux)))

    if np.any(np.abs(input_flux) == 0):
        raise RuntimeError("参考 run 中 input_flux 出现 0，源或 monitor 位置可能有误。")

    # --- 结构 run ---
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=resolution,
        k_point=mp.Vector3(),
        geometry=geometry,
        dimensions=2,
    )
    refl_fr = sim.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0, y_refl, 0),
                      size=mp.Vector3(period_um, 0, 0)),
    )
    trans_fr = sim.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0, y_trans, 0),
                      size=mp.Vector3(period_um, 0, 0)),
    )
    sim.load_minus_flux_data(refl_fr, ref_data)

    t0 = time.time()
    log("[struct run] P=%.3f W=%.3f D=%.3f μm  (W=0 或 D=0 表示平板基线)",
        period_um, groove_width_um, groove_depth_um)
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20, c=mp.Ez,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-decay_db / 10.0),
        )
    )
    refl_flux = np.array(mp.get_fluxes(refl_fr))
    trans_flux = np.array(mp.get_fluxes(trans_fr))
    t_struct = time.time() - t0
    log("[struct run] done in %.2f s", t_struct)

    # R / T / A.  水平 flux 面正方向为 +y；透射波沿 -y 传播，
    # 因此物理透射率使用 signed_transmittance。
    R = refl_flux / np.abs(input_flux)
    raw_transmittance = trans_flux / np.abs(input_flux)
    signed_transmittance = -trans_flux / np.abs(input_flux)
    T = signed_transmittance
    A = 1.0 - R - T
    wavelengths_um = np.array([meep_freq_to_wavelength_um(f) for f in freqs])
    order = np.argsort(wavelengths_um)

    return dict(
        wavelength_um=wavelengths_um[order],
        R=R[order],
        T=T[order],
        A=A[order],
        input_flux_raw=input_flux[order],
        reflection_flux_raw=refl_flux[order],
        transmission_flux_raw=trans_flux[order],
        raw_trans_flux=trans_flux[order],
        raw_transmittance=raw_transmittance[order],
        signed_transmittance=signed_transmittance[order],
        walltime_s=t_ref + t_struct,
        geometry_y=dict(
            cell_y=cell_y,
            y_surface=y_surface,
            y_substrate_bottom=y_substrate_bottom,
            y_src=y_src, y_refl=y_refl, y_trans=y_trans,
        ),
    )


# ---------------------------------------------------------------------------
# 对外 API
# ---------------------------------------------------------------------------

def simulate_periodic_groove(
    period_um: float,
    groove_width_um: float,
    groove_depth_um: float,
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
    """1D 周期矩形槽 Ti 表面正入射光谱仿真。

    Parameters
    ----------
    period_um : float
        x 方向周期 (μm)。
    groove_width_um : float
        槽宽 (μm)，必须 < period_um。
    groove_depth_um : float
        槽深 (μm)，必须 > 0 且 < substrate_thickness_um。
    wavelength_min_um, wavelength_max_um : float
        研究波段端点 (μm)。
    resolution : int
        Meep points/μm。Ti 稳定下限 32。

    Returns
    -------
    pandas.DataFrame
        Columns: wavelength_um, reflectance, transmittance, absorptance,
        emissivity_proxy。 emissivity_proxy ≡ absorptance（不透明基底 +
        Kirchhoff 近似）。

    Raises
    ------
    ValueError
        参数非法（W>=P, D<=0, resolution < 32, 波段非法等）。

    Notes
    -----
    本函数是纯函数（除日志外无副作用），适合批量参数扫描。
    """
    # ---- 参数自检 ----
    if period_um <= 0:
        raise ValueError(f"period_um 必须 > 0, 收到 {period_um}")
    if groove_width_um < 0:
        raise ValueError(f"groove_width_um 必须 >= 0, 收到 {groove_width_um}")
    if groove_depth_um < 0:
        raise ValueError(f"groove_depth_um 必须 >= 0, 收到 {groove_depth_um}")
    if groove_width_um >= period_um:
        raise ValueError(
            f"groove_width_um ({groove_width_um}) 必须 < period_um ({period_um})"
        )
    if wavelength_min_um <= 0 or wavelength_min_um >= wavelength_max_um:
        raise ValueError(
            f"波段非法: [{wavelength_min_um}, {wavelength_max_um}] μm"
        )
    if resolution < MIN_SAFE_RES:
        raise ValueError(
            f"resolution={resolution} 低于 Ti Drude–Lorentz 稳定下限 "
            f"{MIN_SAFE_RES} points/μm。"
        )
    if nfreq <= 1:
        raise ValueError(f"nfreq 必须 > 1, 收到 {nfreq}")

    result = _run_meep_sim(
        period_um=period_um,
        groove_width_um=groove_width_um,
        groove_depth_um=groove_depth_um,
        wavelength_min_um=wavelength_min_um,
        wavelength_max_um=wavelength_max_um,
        resolution=resolution,
        pml_thickness_um=pml_thickness_um,
        substrate_thickness_um=substrate_thickness_um,
        air_buffer_um=air_buffer_um,
        nfreq=nfreq,
        decay_db=decay_db,
        logger=logger,
    )

    df = pd.DataFrame({
        "wavelength_um": result["wavelength_um"],
        "reflectance": result["R"],
        "transmittance": result["T"],
        "absorptance": result["A"],
        "input_flux_raw": result["input_flux_raw"],
        "reflection_flux_raw": result["reflection_flux_raw"],
        "transmission_flux_raw": result["transmission_flux_raw"],
        "raw_transmittance": result["raw_transmittance"],
        "signed_transmittance": result["signed_transmittance"],
        # 不透明基底近似下，ε ≈ A（Kirchhoff）；非透明基底场景需另外处理
        "emissivity_proxy": result["A"],
    })
    df.attrs["period_um"] = period_um
    df.attrs["groove_width_um"] = groove_width_um
    df.attrs["groove_depth_um"] = groove_depth_um
    df.attrs["resolution"] = resolution
    df.attrs["walltime_s"] = result["walltime_s"]
    df.attrs["geometry_y"] = result["geometry_y"]
    return df


# ---------------------------------------------------------------------------
# 自检 + 对照
# ---------------------------------------------------------------------------

def self_checks(df_groove: pd.DataFrame, df_flat: pd.DataFrame, logger,
                *, energy_tol: float = 0.05, T_opaque_tol: float = 1e-2,
                range_tol: float = 1e-3) -> dict:
    """物理 + 数值自检。返回汇总 dict（也写入日志）。"""
    wl = df_groove["wavelength_um"].to_numpy()
    if not np.all(np.diff(wl) > 0):
        raise RuntimeError("波长数组未严格升序。")
    if not np.allclose(wl, df_flat["wavelength_um"].to_numpy()):
        raise RuntimeError("槽与平板的波长轴不一致，无法直接对照。")

    # 能量守恒（复用 src/postprocess）
    rta = {
        "R": df_groove["reflectance"].to_numpy(),
        "T": df_groove["transmittance"].to_numpy(),
        "A": df_groove["absorptance"].to_numpy(),
    }
    energy_report = energy_conservation_check(rta, atol=energy_tol)

    # 不透明假设
    Tmax_groove = float(df_groove["transmittance"].abs().max())
    Tmax_flat = float(df_flat["transmittance"].abs().max())

    # R/T/A 物理范围。A 是由 1-R-T 定义的，所以能量闭合残差本身不能证明结果可靠。
    def _range_report(df: pd.DataFrame) -> dict[str, int]:
        report = {}
        for col in ("reflectance", "transmittance", "absorptance"):
            vals = df[col].to_numpy()
            report[col] = int(np.sum((vals < -range_tol) | (vals > 1 + range_tol)))
        return report

    range_groove = _range_report(df_groove)
    range_flat = _range_report(df_flat)

    # 结构 vs 平板差异
    dR = df_groove["reflectance"].to_numpy() - df_flat["reflectance"].to_numpy()
    dA = df_groove["absorptance"].to_numpy() - df_flat["absorptance"].to_numpy()
    idx_max_dA = int(np.argmax(np.abs(dA)))
    ti_valid_lo, ti_valid_hi = TI_RAKIC_VALID_LAMBDA_UM
    valid_8_13 = (wl >= 8) & (wl <= min(13, ti_valid_hi))

    summary = dict(
        max_energy_error=energy_report["max_abs_dev"],
        mean_energy_error=energy_report["mean_abs_dev"],
        Tmax_groove=Tmax_groove,
        Tmax_flat=Tmax_flat,
        range_violations_groove=range_groove,
        range_violations_flat=range_flat,
        max_abs_dR=float(np.max(np.abs(dR))),
        max_abs_dA=float(np.max(np.abs(dA))),
        lambda_um_of_max_dA=float(wl[idx_max_dA]),
        groove_avg_A_8_13=float(df_groove.loc[
            (wl >= 8) & (wl <= 13), "absorptance"].mean()),
        flat_avg_A_8_13=float(df_flat.loc[
            (wl >= 8) & (wl <= 13), "absorptance"].mean()),
        groove_avg_A_8_valid=float(df_groove.loc[
            valid_8_13, "absorptance"].mean()),
        flat_avg_A_8_valid=float(df_flat.loc[
            valid_8_13, "absorptance"].mean()),
    )

    log = logger.info
    log("self-check: max|R+T+A-1| = %.3e   mean = %.3e",
        summary["max_energy_error"], summary["mean_energy_error"])
    log("self-check: Tmax (groove) = %.3e   Tmax (flat) = %.3e   (应 ≪ 1)",
        Tmax_groove, Tmax_flat)
    log("self-check: range violations outside [-%.1e, 1+%.1e]   groove=%s   flat=%s",
        range_tol, range_tol, range_groove, range_flat)
    log("self-check: max|ΔR| = %.3e   max|ΔA| = %.3e  @ λ = %.3f μm",
        summary["max_abs_dR"], summary["max_abs_dA"],
        summary["lambda_um_of_max_dA"])
    log("self-check: 8–13 μm 平均吸收率   groove=%.3f   flat=%.3f",
        summary["groove_avg_A_8_13"], summary["flat_avg_A_8_13"])
    log("self-check: 8–%.3f μm (Ti valid) 平均吸收率   groove=%.3f   flat=%.3f",
        min(13, ti_valid_hi),
        summary["groove_avg_A_8_valid"], summary["flat_avg_A_8_valid"])

    if Tmax_groove > T_opaque_tol or Tmax_flat > T_opaque_tol:
        logger.warning(
            "T 最大值超过不透明阈值 %.1e；substrate_thickness 可能不够。",
            T_opaque_tol,
        )
    if any(range_groove.values()) or any(range_flat.values()):
        logger.warning(
            "存在 R/T/A 超出物理范围的频点；通常来自材料外推波段、源边缘信噪比或收敛不足。"
        )
    if wl.max() > ti_valid_hi:
        logger.warning(
            "当前最大波长 %.3f μm 超出 Ti Rakić 模型标定上限 %.3f μm；"
            "建议把定量结论限制在有效波段，或替换为覆盖目标波段的光学常数。",
            float(wl.max()), ti_valid_hi,
        )
    if summary["max_abs_dA"] < 1e-3:
        logger.warning(
            "槽与平板几乎无差异 (max|ΔA|<1e-3)，可能 W 太小、D 太浅或网格不够。"
        )
    return summary


# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------

def plot_compare(df_groove: pd.DataFrame, df_flat: pd.DataFrame,
                 out_path: Path, *, title: str = "") -> Path:
    """槽 vs 平板：R / A 谱对照（T 在不透明假设下 ≈ 0 不画）。"""
    wl = df_groove["wavelength_um"].to_numpy()
    fig, (ax_R, ax_A) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True)

    ax_R.plot(wl, df_flat["reflectance"], lw=1.6, color="C0",
              label="flat Ti", ls="--")
    ax_R.plot(wl, df_groove["reflectance"], lw=1.8, color="C3",
              label="periodic groove")
    ax_R.set_ylabel("Reflectance R")
    ax_R.set_ylim(0, 1.05)
    ax_R.axvspan(8, 13, color="orange", alpha=0.08,
                 label="8–13 μm window")
    ax_R.legend(loc="lower right", framealpha=0.9)
    ax_R.grid(True, ls=":", alpha=0.6)
    ax_R.set_title(title or "Ti slab vs periodic rectangular groove (2D, Ez)")

    ax_A.plot(wl, df_flat["absorptance"], lw=1.6, color="C0",
              label="flat Ti  (ε_proxy)", ls="--")
    ax_A.plot(wl, df_groove["absorptance"], lw=1.8, color="C3",
              label="periodic groove  (ε_proxy)")
    ax_A.axvspan(8, 13, color="orange", alpha=0.08)
    ax_A.set_ylabel("Absorptance A ≈ ε  (opaque)")
    ax_A.set_xlabel("Wavelength λ (μm)")
    ax_A.set_ylim(0, max(0.05, 1.05 * float(df_groove["absorptance"].max())))
    ax_A.legend(loc="upper right", framealpha=0.9)
    ax_A.grid(True, ls=":", alpha=0.6)

    return save_figure(fig, out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--period", type=float, required=True,
                   help="x 方向周期 (μm)")
    p.add_argument("--width", type=float, required=True,
                   help="槽宽 (μm)，必须 < period")
    p.add_argument("--depth", type=float, required=True,
                   help="槽深 (μm)，必须 > 0 且 < substrate_thickness")
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
    return p.parse_args()


def main() -> int:
    args = parse_args()
    run_tag = make_run_tag("02_periodic_groove_spectrum")
    logger = setup_logger(run_tag)
    logger.info("=== 02_periodic_groove_spectrum ===")
    logger.info("args = %s", vars(args))

    mp.verbosity(1)

    # ---- 1) 平板基线 ----
    logger.info(">>> Step A: 平板 Ti 基线 (W=0, D=0)")
    df_flat = simulate_periodic_groove(
        period_um=args.period,
        groove_width_um=0.0,
        groove_depth_um=0.0,
        wavelength_min_um=args.wavelength_min_um,
        wavelength_max_um=args.wavelength_max_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        nfreq=args.nfreq,
        logger=logger,
    )

    # ---- 2) 槽结构 ----
    logger.info(">>> Step B: 周期矩形槽 P=%.3f W=%.3f D=%.3f μm",
                args.period, args.width, args.depth)
    df_groove = simulate_periodic_groove(
        period_um=args.period,
        groove_width_um=args.width,
        groove_depth_um=args.depth,
        wavelength_min_um=args.wavelength_min_um,
        wavelength_max_um=args.wavelength_max_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        nfreq=args.nfreq,
        logger=logger,
    )

    # ---- 3) 自检 ----
    summary = self_checks(df_groove, df_flat, logger)

    # ---- 4) CSV：把槽与平板的谱写到同一表里，便于直接对照 ----
    wl = df_groove["wavelength_um"].to_numpy()
    csv_path = save_spectrum_csv(
        wl,
        {
            # 槽
            "reflectance_groove": df_groove["reflectance"].to_numpy(),
            "transmittance_groove": df_groove["transmittance"].to_numpy(),
            "absorptance_groove": df_groove["absorptance"].to_numpy(),
            "emissivity_proxy_groove": df_groove["emissivity_proxy"].to_numpy(),
            "raw_transmittance_groove": df_groove["raw_transmittance"].to_numpy(),
            "signed_transmittance_groove": df_groove["signed_transmittance"].to_numpy(),
            # 平板
            "reflectance_flat": df_flat["reflectance"].to_numpy(),
            "transmittance_flat": df_flat["transmittance"].to_numpy(),
            "absorptance_flat": df_flat["absorptance"].to_numpy(),
            "emissivity_proxy_flat": df_flat["emissivity_proxy"].to_numpy(),
            "raw_transmittance_flat": df_flat["raw_transmittance"].to_numpy(),
            "signed_transmittance_flat": df_flat["signed_transmittance"].to_numpy(),
            # 差异
            "delta_R": df_groove["reflectance"].to_numpy() - df_flat["reflectance"].to_numpy(),
            "delta_A": df_groove["absorptance"].to_numpy() - df_flat["absorptance"].to_numpy(),
            # 能量守恒残差
            "energy_error_groove": np.abs(
                1.0 - df_groove["reflectance"].to_numpy()
                - df_groove["transmittance"].to_numpy()
                - df_groove["absorptance"].to_numpy()
            ),
        },
        project_path("results", "tables", "periodic_groove_spectrum.csv"),
    )
    logger.info("CSV → %s", csv_path)

    # ---- 5) PNG ----
    title = (f"Ti groove vs flat: P={args.period:g} W={args.width:g} "
             f"D={args.depth:g} μm  (res={args.resolution}/μm, Ez)")
    fig_path = plot_compare(
        df_groove, df_flat,
        project_path("results", "figures", "periodic_groove_vs_flat_ti.png"),
        title=title,
    )
    logger.info("PNG → %s", fig_path)

    logger.info(
        "=== done in flat=%.1fs + groove=%.1fs s ===",
        df_flat.attrs["walltime_s"], df_groove.attrs["walltime_s"],
    )
    logger.info("SUMMARY: %s", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
