"""01_flat_ti_benchmark.py — 平面 Ti 中红外光谱基准 (5–15 μm)

研究目的
--------
为后续金属微结构发射率仿真建立 **平板 Ti** 基准：在引入任何几何复杂度之前，
先在 2D 域中跑通 “参考归一化 run + 结构 run” 流程，验证以下要素都正确：
  1. 材料调用：meep.materials.Ti (Rakić 1998 Drude–Lorentz) 的中红外行为；
  2. 边界条件：x 方向周期 + y 方向 PML；
  3. 反射率归一化：先空场 run 记录入射 flux，再用 `load_minus_flux_data`
     扣除入射场；
  4. 后处理 / 单位换算 / IO 全链路（CSV、PNG、日志）。

为什么先做平板 Ti
------------------
- 平板有解析极限可做物理直觉对照：金属在中红外应 **高反射、几乎零透射、
  低吸收 (≈1−R)**；任何明显违反这一点（例如 R<0 或 A>1）都说明流程错。
- 在没有几何细节的情况下，**任何能量守恒偏差都来自数值误差**：分辨率、
  PML、buffer 距离、源带宽。这是后续微结构脚本最干净的“失败检测器”。

为什么反射谱要做参考 run
------------------------
Meep 的源不是纯入射波 —— 它在源平面同时向 ±y 注入场。直接用 monitor 读到
的 flux 既含入射也含反射。因此标准做法是：
  step 1 (reference)：移除结构 (全空气)，跑相同源，记录 monitor 上的
    “净入射 flux” 并保存到磁盘；
  step 2 (with structure)：加上 Ti，载入 step 1 的入射 flux 并 **取反相加**
    (`load_minus_flux_data`)，monitor 上剩下的就是反射场。
这样得到的 R(λ) 才是物理反射率。

输入 / 输出
-----------
输入：仅命令行参数（见 argparse；不依赖外部 YAML，便于单文件复现）。
输出：
  - results/tables/flat_ti_spectrum.csv
  - results/figures/flat_ti_spectrum.png
  - logs/01_flat_ti_benchmark.log
  - 若 --convergence：
      results/tables/flat_ti_convergence.csv
      results/figures/flat_ti_convergence.png

适用范围 / 局限
---------------
- 仅 2D，Ez 偏振 (out-of-plane)；对各向同性 Ti + 正入射，等价于全偏振。
  非正入射 / 微结构破坏对称时需要分别仿 TE/TM 并平均；
- Ti 模型标定波长 0.248–12.4 μm；13–15 μm 是 Drude 尾部外推，定性可信
  但定量需 FTIR 校核；
- 不含氧化层 (TiO/TiO₂)、不含温度依赖、不含表面粗糙度；
- 默认 substrate_thickness_um=2 μm 在中红外足以使 T≈0（趋肤深度 ~10⁻²μm），
  脚本会显式打印 T 的最大值确认这一点。
"""

from __future__ import annotations

import argparse
import cmath
import os
import sys
import time
import warnings
from pathlib import Path

# 让 `python scripts/01_flat_ti_benchmark.py` 能 import src/
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")  # 无显示环境也能保存图
import matplotlib.pyplot as plt
import meep as mp
import numpy as np

from src.io_utils import (
    ensure_dir,
    project_path,
    save_figure,
    save_spectrum_csv,
    setup_logger,
)
from src.materials import (
    freq_range_for_band,
    get_ti_medium,
    meep_freq_to_wavelength_um,
)


# ---------------------------------------------------------------------------
# 默认参数（也可由 argparse 覆盖）
# ---------------------------------------------------------------------------

DEFAULTS = dict(
    wavelength_min_um=5.0,
    wavelength_max_um=15.0,
    resolution=30,          # points / μm，convergence 测试会扫 20/40/60
    pml_thickness_um=2.0,   # PML 厚度（与最长波长同量级即可在正入射下达到 >30 dB 吸收）
    substrate_thickness_um=2.0,  # 不透明假设：中红外 Ti 的趋肤深度 ~10⁻² μm，2 μm 极冗余
    nfreq=201,
    air_buffer_um=8.0,      # 源 / monitor 与结构之间的空气缓冲
    period_x_um=1.0,        # x 方向 (周期方向) cell 宽度，对平板无影响
    decay_db=40.0,          # 衰减阈值
)


# ---------------------------------------------------------------------------
# Meep 单次仿真：返回 (wavelengths_um, R, T, A)
# ---------------------------------------------------------------------------

def _build_geometry(substrate_thickness_um: float, cell_x_um: float,
                    medium_ti) -> list:
    """构造一个上表面在 y=0、向下延伸的 Ti 平板。"""
    return [
        mp.Block(
            material=medium_ti,
            center=mp.Vector3(0, -substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(cell_x_um, substrate_thickness_um, mp.inf),
        )
    ]


def run_flat_ti_sim(
    *,
    wavelength_min_um: float,
    wavelength_max_um: float,
    resolution: int,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    nfreq: int,
    air_buffer_um: float = DEFAULTS["air_buffer_um"],
    period_x_um: float = DEFAULTS["period_x_um"],
    decay_db: float = DEFAULTS["decay_db"],
    logger=None,
) -> dict:
    """对单一 resolution 跑一次 “参考 + 结构” 仿真。

    Returns
    -------
    dict with keys: wavelength_um, R, T, A, freqs, ref_input_flux, walltime_s

    Notes
    -----
    几何布局（沿 y, 从顶到底）：
        [+y]  top PML  (厚 pml_thickness_um)
              空气 buffer
              源平面 (y = y_src)               <-- 平面波，沿 -y 方向传播
              反射 monitor (y = y_refl)        <-- 位于源与结构之间
              空气
              Ti 上表面 (y = 0)
              Ti 基底 (厚 substrate_thickness_um)
              透射 monitor (y = y_trans)       <-- 在 Ti 底面 (PML 边界正上方)
              bottom PML
        [-y]
    """
    if wavelength_min_um >= wavelength_max_um:
        raise ValueError("wavelength_min_um 必须小于 wavelength_max_um。")
    if resolution <= 0 or nfreq <= 0:
        raise ValueError("resolution 与 nfreq 必须为正整数。")
    if pml_thickness_um <= 0 or substrate_thickness_um <= 0 or air_buffer_um <= 0:
        raise ValueError("PML、substrate、air_buffer 厚度必须为正数。")

    log = logger.info if logger else print

    # --- 频率 ---
    f_min, f_max, fcen = freq_range_for_band(wavelength_min_um, wavelength_max_um)
    df = f_max - f_min

    # --- cell 尺寸 ---
    # Meep 的 cell 默认以原点为中心；显式从 cell 边界推导所有 y 坐标，避免
    # source/monitor 被放到仿真域外。Ti 上方留 air_buffer，下方留一段空气缓冲
    # 用于透射 monitor 与 bottom PML 解耦。
    bottom_buffer_um = pml_thickness_um
    cell_y = (
        2 * pml_thickness_um
        + air_buffer_um
        + substrate_thickness_um
        + bottom_buffer_um
    )
    cell = mp.Vector3(period_x_um, cell_y, 0)

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

    # --- 边界 ---
    pml_layers = [mp.PML(thickness=pml_thickness_um, direction=mp.Y)]
    # x 方向不写 PML 即默认周期 (Bloch k=0)

    # --- 源 ---
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=df, is_integrated=True),
            component=mp.Ez,  # 2D 中 Ez = out-of-plane；各向同性 + 正入射对偏振不敏感
            center=mp.Vector3(0, y_src, 0),
            size=mp.Vector3(period_x_um, 0, 0),
        )
    ]

    # --- Ti ---
    medium_ti = get_ti_medium(lambda_min_um=wavelength_min_um,
                              lambda_max_um=wavelength_max_um)
    geometry = _build_geometry(substrate_thickness_um, period_x_um, medium_ti)
    geometry[0].center = mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0)

    # =====================================================================
    # Step 1: 参考 run (空场)
    # =====================================================================
    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=resolution,
        k_point=mp.Vector3(),  # 周期边界 + Γ 点 (正入射)
        geometry=[],            # 关键：空气 only
        dimensions=2,
    )
    refl_fr_ref = sim_ref.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0, y_refl, 0),
                      size=mp.Vector3(period_x_um, 0, 0)),
    )

    t0 = time.time()
    log("[ref run] resolution=%d  cell=(%.2f, %.2f) μm  nfreq=%d  fcen=%.4f  df=%.4f",
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
        raise RuntimeError(
            "参考 run 中 input_flux 出现 0，源或 monitor 位置可能有误。"
        )

    # =====================================================================
    # Step 2: 结构 run (带 Ti 基底)
    # =====================================================================
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
                      size=mp.Vector3(period_x_um, 0, 0)),
    )
    trans_fr = sim.add_flux(
        fcen, df, nfreq,
        mp.FluxRegion(center=mp.Vector3(0, y_trans, 0),
                      size=mp.Vector3(period_x_um, 0, 0)),
    )

    # 关键：把参考 run 的入射场以相反符号载入 refl monitor，使 monitor 只剩反射场
    sim.load_minus_flux_data(refl_fr, ref_data)

    t0 = time.time()
    log("[struct run] adding Ti substrate (thickness=%.2f μm)", substrate_thickness_um)
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

    # =====================================================================
    # R / T / A
    # =====================================================================
    # Meep 约定：refl monitor 法线沿 +y。源沿 -y 传播 → 入射 flux 为负。
    # 减去参考入射场后，refl_flux 留下的是 “反射场”（沿 +y），符号为正。
    # R = refl_flux / |input_flux|.  透射波继续沿 -y 传播，因此 trans_flux 的
    # raw 符号为负；物理 T 取 -trans_flux / |input_flux|。
    R = refl_flux / np.abs(input_flux)
    raw_transmittance = trans_flux / np.abs(input_flux)
    signed_transmittance = -trans_flux / np.abs(input_flux)
    T = signed_transmittance
    A = 1.0 - R - T

    wavelengths_um = np.array([meep_freq_to_wavelength_um(f) for f in freqs])

    # Meep 返回的频率是 fcen 附近降序还是升序取决于版本；统一按波长升序排序
    order = np.argsort(wavelengths_um)
    return dict(
        wavelength_um=wavelengths_um[order],
        R=R[order],
        T=T[order],
        A=A[order],
        freqs=freqs[order],
        ref_input_flux=np.abs(input_flux)[order],
        input_flux_raw=input_flux[order],
        reflection_flux_raw=refl_flux[order],
        transmission_flux_raw=trans_flux[order],
        raw_trans_flux=trans_flux[order],
        raw_transmittance=raw_transmittance[order],
        signed_transmittance=signed_transmittance[order],
        walltime_s=t_ref + t_struct,
        geometry_y=dict(
            cell_y=cell_y,
            top_pml_inner=y_top_pml_inner,
            source=y_src,
            reflection_monitor=y_refl,
            ti_surface=y_surface,
            ti_bottom=y_substrate_bottom,
            transmission_monitor=y_trans,
            bottom_pml_inner=y_bottom_pml_inner,
        ),
    )


def analytic_flat_reflectance_ti(wavelengths_um: np.ndarray) -> np.ndarray:
    """同一 Meep Ti 材料在空气/Ti 半无限界面的正入射 Fresnel 反射率。"""
    from meep.materials import Ti

    values = []
    for wl in wavelengths_um:
        try:
            eps = Ti.epsilon(1.0 / float(wl))[0][0]
        except ValueError:
            values.append(np.nan)
            continue
        n_ti = cmath.sqrt(eps)
        values.append(abs((n_ti - 1.0) / (n_ti + 1.0)) ** 2)
    return np.asarray(values, dtype=float)


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------

def self_checks(result: dict, logger, *, band_lo_um=8.0, band_hi_um=13.0,
                in_range_tol: float = 0.02, energy_tol: float = 0.05) -> dict:
    """对一次仿真结果做物理自检。

    Returns
    -------
    dict with summary fields written into the log.
    """
    wl = result["wavelength_um"]
    R, T, A = result["R"], result["T"], result["A"]

    # 1) 波长升序
    if not np.all(np.diff(wl) > 0):
        raise RuntimeError("波长数组未严格升序。")

    # 2) R/T/A 是否主要落在 [0, 1]
    def _frac_out(arr):
        return float(np.mean((arr < -in_range_tol) | (arr > 1 + in_range_tol)))
    frac_out = {k: _frac_out(v) for k, v in zip("RTA", (R, T, A))}

    # 3) 能量守恒
    energy_error = np.abs(1.0 - (R + T + A))
    max_err = float(np.max(energy_error))
    mean_err = float(np.mean(energy_error))

    # 4) 8–13 μm 平均吸收率
    band_mask = (wl >= band_lo_um) & (wl <= band_hi_um)
    avg_A_band = float(np.mean(A[band_mask])) if band_mask.any() else float("nan")
    avg_R_band = float(np.mean(R[band_mask])) if band_mask.any() else float("nan")
    avg_T_band = float(np.mean(T[band_mask])) if band_mask.any() else float("nan")

    R_analytic = analytic_flat_reflectance_ti(wl)
    valid = np.isfinite(R_analytic)
    max_R_vs_analytic = (
        float(np.max(np.abs(R[valid] - R_analytic[valid])))
        if np.any(valid) else float("nan")
    )

    log = logger.info
    log("self-check: max |R+T+A-1|     = %.3e   mean = %.3e", max_err, mean_err)
    log("self-check: fraction outside [-%g, 1+%g]  R=%.2f%%  T=%.2f%%  A=%.2f%%",
        in_range_tol, in_range_tol,
        100 * frac_out["R"], 100 * frac_out["T"], 100 * frac_out["A"])
    log("self-check: 8-13 μm averages   <R>=%.3f  <T>=%.3e  <A>=%.3f",
        avg_R_band, avg_T_band, avg_A_band)
    log("self-check: T range over full band   [%.3e, %.3e]  (应远 << 1 验证不透明)",
        float(T.min()), float(T.max()))
    log("self-check: max |R - R_Fresnel| over Ti valid band = %.3e",
        max_R_vs_analytic)
    log("self-check: y layout = %s", result.get("geometry_y", {}))

    if max_err > energy_tol:
        logger.warning(
            "能量守恒误差 max=%.3e 超过阈值 %.3e，请考虑提高 resolution "
            "或增大 air_buffer / pml_thickness。", max_err, energy_tol,
        )
    if any(v > 0.05 for v in frac_out.values()):
        logger.warning(
            "R/T/A 超出 [0,1] 的样点比例 > 5%%，疑似存在数值伪影，"
            "建议核查源/monitor 距离与 PML 厚度。"
        )

    return dict(
        max_energy_error=max_err,
        mean_energy_error=mean_err,
        frac_out_of_range=frac_out,
        avg_A_8_13um=avg_A_band,
        avg_R_8_13um=avg_R_band,
        avg_T_8_13um=avg_T_band,
        T_max=float(T.max()),
        max_R_vs_analytic=max_R_vs_analytic,
    )


# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------

def plot_spectrum(result: dict, out_path: Path, *, title: str = "") -> Path:
    wl = result["wavelength_um"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(wl, result["R"], label="Reflectance R", lw=1.8)
    ax.plot(wl, result["T"], label="Transmittance T", lw=1.8)
    ax.plot(wl, result["A"], label="Absorptance A = 1 - R - T", lw=1.8)
    ax.axvspan(8, 13, color="orange", alpha=0.08, label="8-13 μm window")
    ax.set_xlabel("Wavelength λ (μm)")
    ax.set_ylabel("Spectral power fraction")
    ax.set_ylim(-0.05, 1.10)
    ax.set_xlim(wl.min(), wl.max())
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(loc="center right", framealpha=0.9)
    ax.set_title(title or "Flat Ti slab, 2D normal-incidence benchmark")
    return save_figure(fig, out_path)


def plot_convergence(results_by_res: dict[int, dict], out_path: Path) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True)
    ax_R, ax_err = axes
    best_res = sorted(results_by_res.keys())[-1]
    best_R = results_by_res[best_res]["R"]

    for res, r in sorted(results_by_res.items()):
        ax_R.plot(r["wavelength_um"], r["R"], lw=1.5, label=f"res={res}/μm")
        if res == best_res:
            continue
        diff = np.abs(r["R"] - best_R)
        ax_err.semilogy(r["wavelength_um"], diff, lw=1.2,
                        label=f"|R{res}-R{best_res}|")

    for ax in axes:
        ax.grid(True, ls=":", alpha=0.6)
        ax.legend()
    ax_R.set_ylabel("Reflectance R")
    ax_R.set_title("Flat Ti — resolution convergence")
    ax_R.set_ylim(0, 1.05)
    ax_err.set_ylabel("Reflectance difference")
    ax_err.set_xlabel("Wavelength λ (μm)")
    return save_figure(fig, out_path)


# ---------------------------------------------------------------------------
# 可选：分辨率收敛测试
# ---------------------------------------------------------------------------

def run_resolution_convergence_test(
    resolutions: list[int],
    *,
    wavelength_min_um: float,
    wavelength_max_um: float,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    nfreq: int,
    logger,
) -> dict[int, dict]:
    """跑多个 resolution，并保存对比 CSV + PNG。

    评估指标：
      - 不同分辨率下 R(λ) 是否重合；
      - 能量守恒误差是否随 resolution 增加而下降。

    Notes
    -----
    Ti (Rakić 1998) 的最高 Lorentz 共振频率为 15.67 Meep 单位 (a=1μm)。
    Meep 的 Drude–Lorentz 辅助微分方程要求 Nyquist 频率 > 共振频率，即：
        resolution > 2 × 15.67 ≈ 32 points/μm
    低于此值会导致 NaN/Inf 不稳定。默认最低安全分辨率为 30/μm（实测临界稳定）。
    """
    if len(resolutions) < 2:
        raise ValueError("至少要 2 个分辨率才能做收敛对比。")

    # 最高 Lorentz 共振频率 (Meep 单位)；低于此 Nyquist 会不稳定
    TI_MAX_LORENTZ_FREQ = 15.67
    MIN_SAFE_RES = int(2 * TI_MAX_LORENTZ_FREQ) + 1  # 32
    for res in resolutions:
        if res < MIN_SAFE_RES:
            raise ValueError(
                f"resolution={res} 低于 Ti Drude–Lorentz 稳定下限 {MIN_SAFE_RES} points/μm "
                f"(Ti 最高 Lorentz 共振 {TI_MAX_LORENTZ_FREQ:.1f} Meep 单位，"
                f"Nyquist={res/2:.1f} < {TI_MAX_LORENTZ_FREQ:.1f})。"
                f"请使用 resolution >= {MIN_SAFE_RES}。"
            )

    results = {}
    for res in resolutions:
        logger.info("=== convergence: resolution = %d points/μm ===", res)
        results[res] = run_flat_ti_sim(
            wavelength_min_um=wavelength_min_um,
            wavelength_max_um=wavelength_max_um,
            resolution=res,
            pml_thickness_um=pml_thickness_um,
            substrate_thickness_um=substrate_thickness_um,
            nfreq=nfreq,
            logger=logger,
        )

    # CSV: 每个 resolution 占一组列
    wl_ref = results[resolutions[0]]["wavelength_um"]
    columns = {}
    for res, r in results.items():
        # 不同 res 的频率轴是同一个 (Meep flux 频点由 fcen/df/nfreq 决定，与分辨率无关)
        if not np.allclose(r["wavelength_um"], wl_ref):
            raise RuntimeError(f"resolution={res} 的波长轴与基准不一致，无法叠加。")
        columns[f"R_res{res}"] = r["R"]
        columns[f"T_res{res}"] = r["T"]
        columns[f"A_res{res}"] = r["A"]
        columns[f"dR_vs_res{resolutions[-1]}_res{res}"] = np.abs(
            r["R"] - results[resolutions[-1]]["R"]
        )

    csv_path = save_spectrum_csv(
        wl_ref, columns,
        project_path("results", "tables", "flat_ti_convergence.csv"),
    )
    fig_path = plot_convergence(
        results, project_path("results", "figures", "flat_ti_convergence.png"),
    )
    logger.info("convergence CSV → %s", csv_path)
    logger.info("convergence PNG → %s", fig_path)

    # 收敛指标：与最高 resolution 比较
    res_sorted = sorted(results.keys())
    best = results[res_sorted[-1]]
    for res in res_sorted[:-1]:
        diff = np.max(np.abs(results[res]["R"] - best["R"]))
        logger.info("max|R(res=%d) - R(res=%d)| = %.3e", res, res_sorted[-1], diff)

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--wavelength_min_um", type=float, default=DEFAULTS["wavelength_min_um"])
    p.add_argument("--wavelength_max_um", type=float, default=DEFAULTS["wavelength_max_um"])
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--pml_thickness_um", type=float, default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--nfreq", type=int, default=DEFAULTS["nfreq"])
    p.add_argument("--air_buffer_um", type=float, default=DEFAULTS["air_buffer_um"])
    p.add_argument("--convergence", action="store_true",
                   help="额外跑 resolution 收敛测试 (默认 [32, 50])")
    p.add_argument("--convergence_resolutions", type=int, nargs="+",
                   default=[32, 50],
                   help="收敛测试使用的 resolution 列表 (points/μm)；"
                        "Ti 稳定下限 32/μm，低于此值会 NaN/Inf")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    logger = setup_logger("01_flat_ti_benchmark")
    logger.info("=== 01_flat_ti_benchmark ===")
    logger.info("args = %s", vars(args))

    # Meep 安静一点（默认会刷很多 init 信息；保留即可便于调试）
    mp.verbosity(1)

    # ---------- 主仿真 ----------
    result = run_flat_ti_sim(
        wavelength_min_um=args.wavelength_min_um,
        wavelength_max_um=args.wavelength_max_um,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        nfreq=args.nfreq,
        air_buffer_um=args.air_buffer_um,
        logger=logger,
    )

    # ---------- 自检 ----------
    summary = self_checks(result, logger)

    # ---------- 保存 CSV ----------
    energy_error = np.abs(1.0 - (result["R"] + result["T"] + result["A"]))
    R_analytic = analytic_flat_reflectance_ti(result["wavelength_um"])
    csv_path = save_spectrum_csv(
        result["wavelength_um"],
        {
            "reflectance": result["R"],
            "transmittance": result["T"],
            "absorptance": result["A"],
            "input_flux_raw": result["input_flux_raw"],
            "reflection_flux_raw": result["reflection_flux_raw"],
            "transmission_flux_raw": result["transmission_flux_raw"],
            "raw_transmittance": result["raw_transmittance"],
            "signed_transmittance": result["signed_transmittance"],
            "fresnel_reflectance": R_analytic,
            "fresnel_absorptance": 1.0 - R_analytic,
            "energy_error": energy_error,
        },
        project_path("results", "tables", "flat_ti_spectrum.csv"),
    )
    logger.info("CSV → %s", csv_path)

    # ---------- 保存图 ----------
    title = (f"Flat Ti slab benchmark, 2D, normal incidence  "
             f"(res={args.resolution}/μm, "
             f"substrate={args.substrate_thickness_um} μm)")
    fig_path = plot_spectrum(
        result,
        project_path("results", "figures", "flat_ti_spectrum.png"),
        title=title,
    )
    logger.info("PNG → %s", fig_path)

    # ---------- 收敛测试 (可选) ----------
    if args.convergence:
        run_resolution_convergence_test(
            args.convergence_resolutions,
            wavelength_min_um=args.wavelength_min_um,
            wavelength_max_um=args.wavelength_max_um,
            pml_thickness_um=args.pml_thickness_um,
            substrate_thickness_um=args.substrate_thickness_um,
            nfreq=args.nfreq,
            logger=logger,
        )

    logger.info("=== done in %.2f s ===", result["walltime_s"])
    logger.info("SUMMARY: %s", summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
