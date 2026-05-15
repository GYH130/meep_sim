"""04_angle_resolved_emission.py — 斜槽结构的正入射散射远场方向图

研究目的
--------
在 03 (斜槽光谱) 的基础上引入 **角度维度** —— 对单波长正入射条件下，计算
非对称斜槽结构散射到上半空间各方向的远场强度 I(θ)，并提取方向性
指标 (峰值角度、左右积分比、Bragg 级强度、半功率束宽)。这是定向发射建模的
先导探针，而不是严格的角分辨热发射率 ε(θ)。

near-to-far 的三步逻辑
----------------------
Meep 的 near-to-far-field (NTFF) 等价于在近场记录电磁切向分量，然后用 2D
格林函数 (Hankel 函数渐近形式) 把它"投射"到远场任意点。完整流程是三步：

  Step 1 (参考 run)：空气单元，加 NTFF monitor 在 y = y_n2f (位于源与结构
    之间)。运行后保存 NTFF 数据 → 这是"纯入射场"的近场记录。
  Step 2 (结构 run)：放入斜槽几何，同样的源与 NTFF monitor。仿真开始前用
    `sim.load_minus_near2far_data(n2f, ref_data)` 把入射场的 NTFF 数据"减"
    掉 → monitor 在仿真过程中累积的将只是 **散射场** 的 NTFF 数据。
  Step 3 (远场扫描)：对每个观测角 θ，在 r_far 距离处的远场点
    p = r_far · (sin θ, cos θ) 处调用 `sim.get_farfield(n2f, p)`，得到该点
    的 (Ex, Ey, Ez, Hx, Hy, Hz)。Ez 偏振下 |Ez|² 即为角分辨辐射强度。

与 "定向发射" 的关系和边界
---------------------------
- 飞秒激光改性表面的工程目标不仅是 8–13 μm **全空间** 平均高发射率，更要在
  特定方向 (例如对热源探测器方位) 集中辐射，这要求 ε(θ_emit) 在某个方向显著
  高于其它方向；
- 互易 / Kirchhoff 论证严格成立的前提是：
    (a) 材料各向同性、局域热平衡、结构互易；
    (b) 探测发射方向 θ_emit 等价于从同一方向 θ_in 入射并计算吸收
        A(θ_in)=1-R_specular-ΣR_diffracted_orders-T。
- 本脚本做的是 **固定 θ_in=0 正入射，扫描 θ_out 散射方向**。这与 ε(θ_emit)
  是互易关系中的反向问题，不应把本脚本输出直接命名为角分辨发射率；
- 因此本脚本是"定向发射先导探针"：用于判断几何是否会把正入射散射重分配到
  非对称方向。真正的 ε(θ_emit) 角度图应在后续脚本中扫描 θ_in，并对每个 θ_in
  计算 1 - R_specular - R_diffracted_orders (不透明基底下 T≈0)。

Bragg 条件对方向性的限制
------------------------
周期结构的远场方向主要由可传播 Bragg 阶控制，正入射时
    sin(theta_m) = m λ / P。
若 P < λ，只有 m=0 specular 阶可传播，方向图会被法向主峰锁住；若 P≈λ，
±1 阶处在掠射/近临界，数值上容易变窄并且角采样要求高。要看到明显方向性，
通常需要让目标波长满足 P ≥ 1.5 λ_target（±1 阶进入中等角度），或者有意把
P 调到略大于 λ_target，使 ±1 阶处于近临界并精细采样。

2D 远场和真实 3D 样品的差异
---------------------------
1. **波前类型**：2D 远场是柱面波 E ∝ A(θ)/√r，3D 是球面波 E ∝ A(θ,φ)/r。
   |E(θ)|² 在 2D 中表示单位极角的功率密度，3D 中表示单位立体角的功率密度。
   两者的"角分布"概念一致，但绝对单位不同。
2. **周期化**：本脚本仿真单胞 + Bloch 周期 + nperiods 个等价副本求和，模拟
   有限 N (= 2·nperiods + 1) 周期阵列。N → ∞ 时角分布退化为离散 Bragg 锐峰
   (δ 函数)；N 有限时是有限宽度的 sinc-状 lobes。
3. **方位无关**：2D 中无 φ 维，等价于沿 z 方向无穷长样品；真实样品有限
   φ-arc，且斜壁在 3D 中是椭圆而非矩形。
4. **偏振**：2D Ez 仅捕捉一种偏振，3D 中 TE+TM 平均才接近非偏振发射。

输入 / 输出
-----------
输入：CLI (--tilt_angles 0 20, --wavelengths 8 10 12, --period, --top_width, …)
输出：
  - results/tables/angle_resolved_emission.csv  (兼容旧名；内容是散射远场长表)
  - results/tables/angle_resolved_metrics.csv    (tilt, λ, θ_max, 左右积分比, Bragg 指标)
  - results/figures/angle_resolved_cartesian.png (3 panel × 2 line)
  - results/figures/angle_resolved_polar.png     (3 polar panel)
  - logs/04_angle_resolved_emission.log

物理假设与近似（必须随结果一同记录）
----------------------------------
1. 不透明基底 (T ≈ 0)；
2. 2D 等效，槽沿 z 无穷长；本脚本输出的是"散射方向图"，与发射方向图
   只在弱角度依赖区间内可以互换；
3. 单偏振 Ez；
4. nperiods 有限 → 角分布是有限孔径的; 真正无穷周期面只在 Bragg 角处发光；
5. r_far 取有限值 (默认 1000 μm)，远 >> λ，仍是渐近近似；
6. Ti Rakić 1998 在 12 μm 附近已贴近标定上界 12.4 μm；
7. 室温、无氧化层、无圆角。
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meep as mp
import numpy as np
import pandas as pd

from src.geometry import build_slanted_groove_geometry
from src.io_utils import (
    ensure_dir,
    project_path,
    save_figure,
    setup_logger,
)
from src.materials import get_ti_medium


# ---------------------------------------------------------------------------
# 默认参数
# ---------------------------------------------------------------------------

DEFAULTS = dict(
    period_um=10.0,
    top_width_um=4.0,
    bottom_width_um=4.0,
    depth_um=3.0,
    tilt_angles_deg=[0.0, 20.0],          # 对称 + 非对称
    wavelengths_um=[8.0, 10.0, 12.0],
    angle_min_deg=-90.0,
    angle_max_deg=+90.0,
    angle_step_deg=0.5,
    resolution=32,
    pml_thickness_um=2.0,
    substrate_thickness_um=4.0,
    air_buffer_um=8.0,
    nperiods=10,                          # 有限阵列宽度 (2*N+1 = 21 periods)
    r_far_um=1000.0,
    src_fractional_bw=0.05,               # 窄带 (Gaussian 5% bandwidth)
    decay_db=40.0,
)

# Ti 稳定下限 (memory: feedback_ti_resolution)
TI_MAX_LORENTZ_FREQ = 15.67
MIN_SAFE_RES = int(2 * TI_MAX_LORENTZ_FREQ) + 1  # 32


# ---------------------------------------------------------------------------
# Meep 仿真装配（参考 + 结构共用同一函数）
# ---------------------------------------------------------------------------

def _build_sim(
    *,
    wavelength_um: float,
    fractional_bw: float,
    period_um: float,
    resolution: int,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    air_buffer_um: float,
    nperiods: int,
    geometry_factory,  # callable(y_surface, substrate_thickness) -> list, 或 None
) -> tuple:
    """构造 meep.Simulation + 一个 NTFF monitor + 一个 flux monitor。

    几何布局沿 y 自上而下：
        [+y] top PML / 空气 buffer / 源 / NTFF & flux monitor /
             空气 / 结构 / bottom PML [-y]
    NTFF monitor 放在源与结构之间，weight=+1 → 把上半空间 (+y) 作为远场观测域。

    若 geometry_factory is None，则空气单元（参考 run）。
    """
    if wavelength_um <= 0:
        raise ValueError(f"wavelength_um 必须 > 0, 收到 {wavelength_um}")
    if fractional_bw <= 0 or fractional_bw >= 1:
        raise ValueError(f"fractional_bw 必须 ∈ (0,1), 收到 {fractional_bw}")

    fcen = 1.0 / wavelength_um
    df = fractional_bw * fcen

    bottom_buffer_um = pml_thickness_um
    cell_y = (
        2 * pml_thickness_um
        + air_buffer_um
        + substrate_thickness_um
        + bottom_buffer_um
    )
    cell = mp.Vector3(period_um, cell_y, 0)

    y_top_edge = +0.5 * cell_y
    y_bottom_edge = -0.5 * cell_y
    y_top_pml_inner = y_top_edge - pml_thickness_um
    y_bottom_pml_inner = y_bottom_edge + pml_thickness_um
    y_surface = y_top_pml_inner - air_buffer_um
    y_src = y_top_pml_inner - 0.25 * air_buffer_um
    y_n2f = y_surface + 0.5 * air_buffer_um  # 在源与结构之间

    if not (y_surface < y_n2f < y_src < y_top_pml_inner):
        raise RuntimeError("源 / NTFF monitor 的 y 坐标顺序异常。")

    pml_layers = [mp.PML(thickness=pml_thickness_um, direction=mp.Y)]
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=df, is_integrated=True),
            component=mp.Ez,
            center=mp.Vector3(0, y_src, 0),
            size=mp.Vector3(period_um, 0, 0),
        )
    ]

    geometry = [] if geometry_factory is None else geometry_factory(y_surface, substrate_thickness_um)

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=resolution,
        k_point=mp.Vector3(),  # 正入射
        geometry=geometry,
        dimensions=2,
    )
    n2f = sim.add_near2far(
        fcen, 0, 1,
        mp.Near2FarRegion(
            center=mp.Vector3(0, y_n2f, 0),
            size=mp.Vector3(period_um, 0, 0),
            weight=+1,  # +y 为外向法线
        ),
        nperiods=nperiods,
    )
    # flux monitor 用于检查总散射通量（散场能量 vs 近场入射）
    flux = sim.add_flux(
        fcen, df, 1,
        mp.FluxRegion(center=mp.Vector3(0, y_n2f, 0),
                      size=mp.Vector3(period_um, 0, 0)),
    )
    meta = dict(fcen=fcen, df=df, y_n2f=y_n2f, y_surface=y_surface,
                cell_y=cell_y)
    return sim, n2f, flux, meta


def _run_until_decayed(sim, *, y_probe: float, decay_db: float, label: str,
                       logger) -> float:
    t0 = time.time()
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20, c=mp.Ez,
            pt=mp.Vector3(0, y_probe, 0),
            decay_by=10 ** (-decay_db / 10.0),
        )
    )
    dt = time.time() - t0
    logger.info("[%s] done in %.2f s", label, dt)
    return dt


def _scan_far_field(sim, n2f, *, angles_deg: np.ndarray,
                    r_far_um: float) -> np.ndarray:
    """在以原点为中心、半径 r_far 的圆弧上扫角度，返回 |Ez|².

    角度约定：θ = 0° 对应 +y (法线)，θ > 0 朝 +x，θ < 0 朝 -x。
    """
    angles_rad = np.radians(angles_deg)
    intensities = np.empty(angles_deg.shape, dtype=float)
    for i, th in enumerate(angles_rad):
        pt = mp.Vector3(r_far_um * np.sin(th), r_far_um * np.cos(th), 0)
        far = sim.get_farfield(n2f, pt)
        # far 是长度 6 的 (Ex, Ey, Ez, Hx, Hy, Hz) 复数；2D Ez 源时 Ez 主导。
        Ez = far[2]
        intensities[i] = float(np.abs(Ez) ** 2)
    return intensities


# ---------------------------------------------------------------------------
# 对外 API: 一个 (tilt, wavelength) 的远场扫描；并按波长复用参考 run
# ---------------------------------------------------------------------------

def run_angle_resolved_sweep(
    *,
    period_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilts_deg: list[float],
    wavelengths_um: list[float],
    angles_deg: np.ndarray,
    resolution: int = DEFAULTS["resolution"],
    pml_thickness_um: float = DEFAULTS["pml_thickness_um"],
    substrate_thickness_um: float = DEFAULTS["substrate_thickness_um"],
    air_buffer_um: float = DEFAULTS["air_buffer_um"],
    nperiods: int = DEFAULTS["nperiods"],
    r_far_um: float = DEFAULTS["r_far_um"],
    src_fractional_bw: float = DEFAULTS["src_fractional_bw"],
    decay_db: float = DEFAULTS["decay_db"],
    logger=None,
) -> dict[tuple[float, float], dict]:
    """对每个 (tilt, λ) 跑一次散射场角分辨远场。

    Returns
    -------
    dict[(tilt, λ), {intensity, R_scalar, scattered_flux, ref_input_flux,
                     walltime_s}]
    """
    if resolution < MIN_SAFE_RES:
        raise ValueError(
            f"resolution={resolution} 低于 Ti Drude–Lorentz 稳定下限 {MIN_SAFE_RES}。"
        )
    if r_far_um <= max(wavelengths_um) * 5:
        raise ValueError(
            f"r_far_um={r_far_um} 太小，应远 >> λ_max={max(wavelengths_um)} (建议 ≥ 5·λ_max)。"
        )
    if not (np.all(angles_deg >= -90.0) and np.all(angles_deg <= 90.0)):
        raise ValueError("angles_deg 必须 ∈ [-90, 90]。")
    if logger is None:
        raise ValueError("logger 必填，便于追踪长时间运行的进度。")

    results: dict[tuple[float, float], dict] = {}

    for wl in wavelengths_um:
        logger.info("============ wavelength = %.2f μm ============", wl)
        medium_ti = get_ti_medium(lambda_min_um=wl, lambda_max_um=wl)

        # --- Step 1: 参考 run (空气) ---
        sim_ref, n2f_ref, flux_ref, meta = _build_sim(
            wavelength_um=wl,
            fractional_bw=src_fractional_bw,
            period_um=period_um,
            resolution=resolution,
            pml_thickness_um=pml_thickness_um,
            substrate_thickness_um=substrate_thickness_um,
            air_buffer_um=air_buffer_um,
            nperiods=nperiods,
            geometry_factory=None,
        )
        logger.info("[ref] fcen=%.4f, df=%.4f, cell_y=%.2f μm",
                    meta["fcen"], meta["df"], meta["cell_y"])
        _run_until_decayed(sim_ref, y_probe=meta["y_n2f"],
                           decay_db=decay_db, label=f"ref λ={wl}μm",
                           logger=logger)
        ref_n2f_data = sim_ref.get_near2far_data(n2f_ref)
        ref_flux_data = sim_ref.get_flux_data(flux_ref)
        ref_input_flux = float(np.array(mp.get_fluxes(flux_ref))[0])
        if ref_input_flux == 0:
            raise RuntimeError(f"λ={wl}: ref_input_flux == 0; 源/monitor 异常。")
        # 入射波沿 -y 传播 → flux 通过 +y 法向面应为负
        if ref_input_flux > 0:
            logger.warning("λ=%g: ref_input_flux > 0 (=%.3e), 与预期 (-y 传播) 相反。",
                           wl, ref_input_flux)

        # --- Step 2 & 3: 每个 tilt 跑一次结构 run + angle scan ---
        for tilt in tilts_deg:
            t_case0 = time.time()

            def factory(y_surface_um, substrate_thickness,
                        _tilt=tilt, _ti=medium_ti):
                return build_slanted_groove_geometry(
                    period_x_um=period_um,
                    top_width_um=top_width_um,
                    bottom_width_um=bottom_width_um,
                    depth_um=depth_um,
                    tilt_angle_deg=_tilt,
                    substrate_thickness_um=substrate_thickness,
                    y_surface=y_surface_um,
                    medium_substrate=_ti,
                )

            sim, n2f, flux, _ = _build_sim(
                wavelength_um=wl,
                fractional_bw=src_fractional_bw,
                period_um=period_um,
                resolution=resolution,
                pml_thickness_um=pml_thickness_um,
                substrate_thickness_um=substrate_thickness_um,
                air_buffer_um=air_buffer_um,
                nperiods=nperiods,
                geometry_factory=factory,
            )
            # 减去入射场的近场记录 → 仿真累积的将仅是散射场
            sim.load_minus_near2far_data(n2f, ref_n2f_data)
            sim.load_minus_flux_data(flux, ref_flux_data)

            _run_until_decayed(sim, y_probe=meta["y_n2f"],
                               decay_db=decay_db,
                               label=f"struct λ={wl}μm α={tilt}°",
                               logger=logger)

            scattered_flux = float(np.array(mp.get_fluxes(flux))[0])
            R_scalar = scattered_flux / abs(ref_input_flux)
            logger.info("λ=%g α=%g: scattered_flux=%.3e, R=%.3f",
                        wl, tilt, scattered_flux, R_scalar)

            intensity = _scan_far_field(sim, n2f,
                                        angles_deg=angles_deg,
                                        r_far_um=r_far_um)
            t_case = time.time() - t_case0
            results[(tilt, wl)] = dict(
                intensity=intensity,
                R_scalar=R_scalar,
                scattered_flux=scattered_flux,
                ref_input_flux=ref_input_flux,
                walltime_s=t_case,
            )
            # 显式释放：Meep Simulation 持仓内存
            del sim

        del sim_ref

    return results


# ---------------------------------------------------------------------------
# 方向性指标
# ---------------------------------------------------------------------------

def diffraction_orders(period_um: float, wavelength_um: float) -> list[tuple[int, float]]:
    """返回正入射下可传播 Bragg 衍射级 (m, theta_deg)，sin(theta)=m*lambda/P。"""
    if period_um <= 0 or wavelength_um <= 0:
        raise ValueError("period_um 与 wavelength_um 必须为正。")
    m_max = int(np.floor(period_um / wavelength_um))
    orders = []
    for m in range(-m_max, m_max + 1):
        s = m * wavelength_um / period_um
        if abs(s) <= 1.0:
            orders.append((m, float(np.degrees(np.arcsin(s)))))
    return orders


def recommended_angle_step_deg(
    *,
    period_um: float,
    wavelength_um: float,
    nperiods: int,
    samples_per_lobe: float = 4.0,
) -> float:
    """有限阵列主瓣采样建议。

    有限阵列宽度 L=(2*nperiods+1)P 时，衍射瓣角宽量级约 λ/L（弧度）。
    为避免漏峰，建议每个主瓣至少采样 samples_per_lobe 个点。
    """
    aperture_um = max(1, 2 * nperiods + 1) * period_um
    lobe_width_deg = float(np.degrees(wavelength_um / aperture_um))
    return max(0.05, lobe_width_deg / samples_per_lobe)


def bragg_regime(period_um: float, wavelength_um: float) -> dict:
    """诊断周期/波长组合是否可能产生明显非 specular 方向性。"""
    orders = diffraction_orders(period_um, wavelength_um)
    nonzero = [(m, th) for m, th in orders if m != 0]
    ratio = period_um / wavelength_um
    if not nonzero:
        regime = "specular_only"
        guidance = "P < λ，仅 m=0 阶传播；倾角对远场方向性的影响通常很弱。"
    elif all(abs(th) >= 75.0 for _m, th in nonzero):
        regime = "near_grazing"
        guidance = "±1 阶接近掠射/临界；峰窄且角采样敏感。"
    elif ratio >= 1.5:
        regime = "multi_order_moderate_angles"
        guidance = "±1 阶进入中等角度，更有机会观察到几何非对称带来的方向性。"
    else:
        regime = "multi_order_high_angles"
        guidance = "存在非零 Bragg 阶，但角度偏大；方向性可能有限，需检查 Bragg 级强度比。"
    return dict(
        period_over_wavelength=ratio,
        orders=orders,
        nonzero_order_count=len(nonzero),
        regime=regime,
        guidance=guidance,
    )


def _order_label(m: int) -> str:
    return f"p{m}" if m >= 0 else f"m{abs(m)}"


def compute_directionality_metrics(angles_deg: np.ndarray,
                                   intensity: np.ndarray,
                                   *,
                                   period_um: float | None = None,
                                   wavelength_um: float | None = None) -> dict:
    """返回 {theta_max_deg, front_back_ratio, half_power_beamwidth_deg}.

    front_back_ratio = I(θ_max) / I(-θ_max)。tilt=0 + 对称几何 → ≈ 1。
    right_left_integral_ratio = ∫_{θ>0} I dθ / ∫_{θ<0} I dθ，更适合衡量左右偏向。
    HPBW: 主峰附近 I 下降到 0.5·I_max 的全宽 (degrees)；若主峰一侧未下降到
          半高即触及 [-90, +90] 边界，HPBW 用边界值替代并记为下限估计。
    """
    if angles_deg.shape != intensity.shape:
        raise ValueError("angles_deg 与 intensity 长度不一致。")
    if intensity.size < 3:
        raise ValueError("intensity 点数 < 3，无法计算 HPBW。")

    i_max = int(np.argmax(intensity))
    theta_max = float(angles_deg[i_max])
    I_max = float(intensity[i_max])

    if I_max <= 0:
        return dict(theta_max_deg=theta_max,
                    front_back_ratio=float("nan"),
                    right_left_integral_ratio=float("nan"),
                    right_left_asymmetry=float("nan"),
                    half_power_beamwidth_deg=float("nan"))

    # front_back_ratio: 线性插值 I(-θ_max)
    I_back = float(np.interp(-theta_max, angles_deg, intensity))
    fbr = I_max / I_back if I_back > 0 else float("inf")

    left_mask = angles_deg < 0
    right_mask = angles_deg > 0
    left_int = float(np.trapezoid(intensity[left_mask], angles_deg[left_mask])) if np.any(left_mask) else 0.0
    right_int = float(np.trapezoid(intensity[right_mask], angles_deg[right_mask])) if np.any(right_mask) else 0.0
    rl_ratio = right_int / left_int if left_int > 0 else float("inf")
    rl_asym = ((right_int - left_int) / (right_int + left_int)
               if (right_int + left_int) > 0 else float("nan"))

    # HPBW: 主峰附近 FWHM
    half = I_max * 0.5

    # 向左 (索引减小方向)
    j = i_max
    while j > 0 and intensity[j] > half:
        j -= 1
    if intensity[j] > half:
        theta_left = float(angles_deg[0])  # 触底
        hpbw_truncated_left = True
    else:
        # 在 j 和 j+1 之间插值
        a, b = float(angles_deg[j]), float(angles_deg[j + 1])
        ia, ib = float(intensity[j]), float(intensity[j + 1])
        if ib == ia:
            theta_left = a
        else:
            theta_left = a + (half - ia) / (ib - ia) * (b - a)
        hpbw_truncated_left = False

    # 向右
    j = i_max
    while j < intensity.size - 1 and intensity[j] > half:
        j += 1
    if intensity[j] > half:
        theta_right = float(angles_deg[-1])
        hpbw_truncated_right = True
    else:
        a, b = float(angles_deg[j - 1]), float(angles_deg[j])
        ia, ib = float(intensity[j - 1]), float(intensity[j])
        if ib == ia:
            theta_right = b
        else:
            theta_right = a + (ia - half) / (ia - ib) * (b - a)
        hpbw_truncated_right = False

    hpbw = float(theta_right - theta_left)
    metrics = dict(
        theta_max_deg=theta_max,
        front_back_ratio=fbr,
        right_left_integral_ratio=rl_ratio,
        right_left_asymmetry=rl_asym,
        half_power_beamwidth_deg=hpbw,
        hpbw_truncated=bool(hpbw_truncated_left or hpbw_truncated_right),
    )

    if period_um is not None and wavelength_um is not None:
        orders = diffraction_orders(period_um, wavelength_um)
        metrics["bragg_orders"] = ";".join(f"{m}:{theta:.3f}" for m, theta in orders)
        regime = bragg_regime(period_um, wavelength_um)
        metrics["period_over_wavelength"] = regime["period_over_wavelength"]
        metrics["nonzero_bragg_order_count"] = regime["nonzero_order_count"]
        metrics["bragg_regime"] = regime["regime"]
        metrics["bragg_guidance"] = regime["guidance"]
        for m, theta in orders:
            label = _order_label(m)
            metrics[f"bragg_{label}_theta_deg"] = theta
            metrics[f"bragg_{label}_I_norm"] = float(np.interp(theta, angles_deg, intensity) / I_max)
        if any(m == -1 for m, _ in orders) and any(m == 1 for m, _ in orders):
            I_minus = float(np.interp(dict(orders)[-1], angles_deg, intensity))
            I_plus = float(np.interp(dict(orders)[1], angles_deg, intensity))
            metrics["bragg_plus_minus_1_ratio"] = I_plus / I_minus if I_minus > 0 else float("inf")
        else:
            metrics["bragg_plus_minus_1_ratio"] = float("nan")
    return metrics


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------

def self_checks(
    results: dict[tuple[float, float], dict],
    angles_deg: np.ndarray,
    *,
    logger,
    period_um: float,
    nperiods: int,
    symmetry_tol: float = 0.10,
    smoothness_tol: float = 0.50,
) -> dict[tuple[float, float], dict]:
    """对每个 (tilt, λ) 做自检：
       1) 归一化 max == 1；
       2) 角分布相邻点跳变 < smoothness_tol；
       3) tilt = 0 时 I(θ) ≈ I(-θ)；tilt ≠ 0 不强制对称；
       4) 远场角积分 vs 近场散射通量比值 (仅日志，不阻塞)。
    """
    angles = np.asarray(angles_deg)
    dtheta_rad = float(np.median(np.diff(np.radians(angles))))
    summary = {}
    for (tilt, wl), case in results.items():
        I = case["intensity"]
        if not np.all(np.isfinite(I)):
            raise RuntimeError(f"(tilt={tilt}, λ={wl}) intensity 含非有限值。")
        Imax = float(np.max(I))
        if Imax <= 0:
            raise RuntimeError(f"(tilt={tilt}, λ={wl}) intensity 全为 0。")
        I_norm = I / Imax
        if abs(np.max(I_norm) - 1.0) > 1e-9:
            raise RuntimeError(f"(tilt={tilt}, λ={wl}) normalized max != 1。")

        angle_step = float(np.median(np.diff(angles)))
        recommended_step = recommended_angle_step_deg(
            period_um=period_um,
            wavelength_um=wl,
            nperiods=nperiods,
        )
        if angle_step > recommended_step:
            logger.warning(
                "(tilt=%g°, λ=%g): angle_step=%.3f° 粗于有限阵列衍射瓣建议 %.3f°；"
                "峰值角、HPBW 与 Bragg 级强度可能欠采样。",
                tilt, wl, angle_step, recommended_step,
            )

        # 平滑性: 相邻点跳变占 max 的比例
        max_jump = float(np.max(np.abs(np.diff(I_norm))))
        if max_jump > smoothness_tol:
            logger.warning(
                "(tilt=%g°, λ=%g): max|ΔI_norm|=%.3f 超过平滑阈值 %.2f，"
                "可能 nperiods 太大产生窄 lobes 或角度采样不足。",
                tilt, wl, max_jump, smoothness_tol,
            )

        # 对称性: 仅 tilt=0 强制
        I_flip = np.interp(-angles, angles, I)
        sym_err = float(np.max(np.abs(I - I_flip)) / Imax)
        if abs(tilt) < 1e-9:
            if sym_err > symmetry_tol:
                logger.warning(
                    "(tilt=0, λ=%g): 对称误差 %.3e 超过阈值 %.2e；"
                    "对称几何下应近似 I(θ)=I(-θ)。",
                    wl, sym_err, symmetry_tol,
                )
            else:
                logger.info("self-check (tilt=0, λ=%g): symmetry err = %.3e ✓",
                            wl, sym_err)
        else:
            logger.info(
                "self-check (tilt=%g°, λ=%g): asymmetry err = %.3e (允许 ≠ 0)",
                tilt, wl, sym_err,
            )

        # 远场积分 vs 近场散射通量 (一致性 sanity，仅记录)
        far_int_arb = float(np.sum(I) * dtheta_rad)
        ratio = far_int_arb / abs(case["scattered_flux"]) if case["scattered_flux"] != 0 else float("nan")
        logger.info(
            "(tilt=%g°, λ=%g): ∫I dθ = %.3e (a.u.) ; "
            "near-field |scattered_flux| = %.3e ; ratio = %.3e",
            tilt, wl, far_int_arb, abs(case["scattered_flux"]), ratio,
        )

        summary[(tilt, wl)] = dict(
            max_jump_normalized=max_jump,
            symmetry_err=sym_err,
            far_field_integral_au=far_int_arb,
            far_to_near_ratio_au=ratio,
            recommended_angle_step_deg=recommended_step,
        )
    return summary


# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------

def _color_for_tilt(tilt: float, tilts_sorted: list[float]) -> tuple:
    cmap = plt.get_cmap("plasma")
    n = max(len(tilts_sorted) - 1, 1)
    return cmap(0.15 + 0.7 * (tilts_sorted.index(tilt) / n))


def plot_cartesian(
    results: dict[tuple[float, float], dict],
    angles_deg: np.ndarray,
    wavelengths: list[float],
    tilts: list[float],
    out_path: Path,
    *,
    title: str = "",
) -> Path:
    fig, axes = plt.subplots(1, len(wavelengths),
                             figsize=(4.0 * len(wavelengths), 4.0),
                             sharey=True)
    if len(wavelengths) == 1:
        axes = [axes]
    tilts_sorted = sorted(tilts)
    for ax, wl in zip(axes, wavelengths):
        for tilt in tilts_sorted:
            I = results[(tilt, wl)]["intensity"]
            I_norm = I / np.max(I)
            ax.plot(angles_deg, I_norm,
                    color=_color_for_tilt(tilt, tilts_sorted),
                    lw=1.6, label=f"α = {tilt:g}°")
        ax.set_xlabel("θ (deg, from +y normal)")
        ax.set_title(f"λ = {wl:g} μm")
        ax.set_xlim(angles_deg.min(), angles_deg.max())
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, ls=":", alpha=0.6)
        ax.axvline(0, color="0.5", ls=":", lw=0.8)
        ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    axes[0].set_ylabel("Normalized far-field intensity |Ez|²")
    fig.suptitle(title or "Angle-resolved scattered far field (2D, Ez)",
                 fontsize=11)
    fig.tight_layout()
    return save_figure(fig, out_path)


def plot_polar(
    results: dict[tuple[float, float], dict],
    angles_deg: np.ndarray,
    wavelengths: list[float],
    tilts: list[float],
    out_path: Path,
    *,
    title: str = "",
) -> Path:
    fig = plt.figure(figsize=(4.0 * len(wavelengths), 4.0))
    tilts_sorted = sorted(tilts)
    for k, wl in enumerate(wavelengths):
        ax = fig.add_subplot(1, len(wavelengths), k + 1, projection="polar")
        ax.set_theta_zero_location("N")  # 0° 朝上 (+y 法线)
        ax.set_theta_direction(-1)       # 顺时针为 +θ → 与 sin θ 朝 +x 一致
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        for tilt in tilts_sorted:
            I = results[(tilt, wl)]["intensity"]
            I_norm = I / np.max(I)
            ax.plot(np.radians(angles_deg), I_norm,
                    color=_color_for_tilt(tilt, tilts_sorted),
                    lw=1.6, label=f"α = {tilt:g}°")
        ax.set_ylim(0, 1.05)
        ax.set_title(f"λ = {wl:g} μm", pad=12)
        ax.grid(True, ls=":", alpha=0.6)
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.15),
                  fontsize=9, ncol=len(tilts_sorted))
    fig.suptitle(title or "Polar far-field |Ez|² (upper hemisphere, normalized)",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    return save_figure(fig, out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--period", type=float, default=DEFAULTS["period_um"])
    p.add_argument("--top_width", type=float, default=DEFAULTS["top_width_um"])
    p.add_argument("--bottom_width", type=float, default=DEFAULTS["bottom_width_um"])
    p.add_argument("--depth", type=float, default=DEFAULTS["depth_um"])
    p.add_argument("--tilt_angles", type=float, nargs="+",
                   default=DEFAULTS["tilt_angles_deg"],
                   help="倾角列表 (°)，默认 [0, 20] = 一个对称 + 一个非对称")
    p.add_argument("--wavelengths", type=float, nargs="+",
                   default=DEFAULTS["wavelengths_um"],
                   help="波长列表 (μm)，默认 [8, 10, 12]")
    p.add_argument("--angle_min", type=float, default=DEFAULTS["angle_min_deg"])
    p.add_argument("--angle_max", type=float, default=DEFAULTS["angle_max_deg"])
    p.add_argument("--angle_step", type=float, default=DEFAULTS["angle_step_deg"],
                   help="远场角度采样步长 (deg)。有限阵列衍射瓣通常建议 0.25–0.5°。")
    p.add_argument("--resolution", type=int, default=DEFAULTS["resolution"])
    p.add_argument("--pml_thickness_um", type=float,
                   default=DEFAULTS["pml_thickness_um"])
    p.add_argument("--substrate_thickness_um", type=float,
                   default=DEFAULTS["substrate_thickness_um"])
    p.add_argument("--air_buffer_um", type=float,
                   default=DEFAULTS["air_buffer_um"])
    p.add_argument("--nperiods", type=int, default=DEFAULTS["nperiods"])
    p.add_argument("--r_far_um", type=float, default=DEFAULTS["r_far_um"])
    p.add_argument("--src_fractional_bw", type=float,
                   default=DEFAULTS["src_fractional_bw"])
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logger = setup_logger("04_angle_resolved_emission")
    logger.info("=== 04_angle_resolved_scattering ===")
    logger.info("args = %s", vars(args))
    logger.warning(
        "本脚本输出正入射散射远场方向图，不是严格角分辨热发射率 ε(θ)。"
    )
    logger.warning(
        "若目标是 ε(θ_emit)，应扫 θ_in 并计算 A(θ_in)=1-R_specular-ΣR_diffracted_orders。"
    )

    mp.verbosity(1)

    angles_deg = np.arange(args.angle_min,
                           args.angle_max + 0.5 * args.angle_step,
                           args.angle_step)
    if angles_deg[-1] > args.angle_max + 1e-6:
        angles_deg = angles_deg[:-1]
    logger.info("angles: %d 个点 ∈ [%.1f, %.1f] step %.1f°",
                angles_deg.size, angles_deg[0], angles_deg[-1], args.angle_step)

    weak_directionality_cases = 0
    for wl in args.wavelengths:
        regime = bragg_regime(args.period, wl)
        orders_txt = ", ".join(f"m={m}: θ={th:.2f}°" for m, th in regime["orders"])
        logger.info(
            "Bragg diagnostic λ=%g μm: P/λ=%.3f, orders=[%s], regime=%s",
            wl, regime["period_over_wavelength"], orders_txt, regime["regime"],
        )
        logger.info("Bragg guidance λ=%g μm: %s", wl, regime["guidance"])
        recommended_step = recommended_angle_step_deg(
            period_um=args.period,
            wavelength_um=wl,
            nperiods=args.nperiods,
        )
        if args.angle_step > recommended_step:
            logger.warning(
                "λ=%g μm: 当前 angle_step=%.3f° 粗于建议 %.3f°。",
                wl, args.angle_step, recommended_step,
            )
        if regime["regime"] in {"specular_only", "near_grazing"}:
            weak_directionality_cases += 1
    if weak_directionality_cases >= max(1, len(args.wavelengths) // 2):
        logger.warning(
            "多数目标波长处于 specular-only 或 near-grazing Bragg 条件；"
            "当前 P=%.3g μm 下倾角方向性预计较弱。建议尝试 P≈1.5λ_target 到 2λ_target "
            "（例如 12–18 μm），或把 P 调到略大于目标 λ 并用更细角度采样研究近临界 ±1 阶。",
            args.period,
        )

    # --- 主仿真扫描 ---
    results = run_angle_resolved_sweep(
        period_um=args.period,
        top_width_um=args.top_width,
        bottom_width_um=args.bottom_width,
        depth_um=args.depth,
        tilts_deg=args.tilt_angles,
        wavelengths_um=args.wavelengths,
        angles_deg=angles_deg,
        resolution=args.resolution,
        pml_thickness_um=args.pml_thickness_um,
        substrate_thickness_um=args.substrate_thickness_um,
        air_buffer_um=args.air_buffer_um,
        nperiods=args.nperiods,
        r_far_um=args.r_far_um,
        src_fractional_bw=args.src_fractional_bw,
        logger=logger,
    )

    # --- 自检 ---
    check_summary = self_checks(
        results, angles_deg, logger=logger,
        period_um=args.period, nperiods=args.nperiods,
    )

    # --- 方向性指标 ---
    metrics_rows = []
    for (tilt, wl), case in results.items():
        m = compute_directionality_metrics(
            angles_deg, case["intensity"],
            period_um=args.period,
            wavelength_um=wl,
        )
        row = {
            "period_um": args.period,
            "top_width_um": args.top_width,
            "bottom_width_um": args.bottom_width,
            "depth_um": args.depth,
            "nperiods": args.nperiods,
            "tilt_angle_deg": tilt,
            "wavelength_um": wl,
            "theta_max_deg": m["theta_max_deg"],
            "front_back_ratio": m["front_back_ratio"],
            "right_left_integral_ratio": m["right_left_integral_ratio"],
            "right_left_asymmetry": m["right_left_asymmetry"],
            "bragg_orders": m.get("bragg_orders", ""),
            "bragg_plus_minus_1_ratio": m.get("bragg_plus_minus_1_ratio", float("nan")),
            "half_power_beamwidth_deg": m["half_power_beamwidth_deg"],
            "hpbw_truncated": m["hpbw_truncated"],
            "near_field_R": case["R_scalar"],
            "walltime_s": case["walltime_s"],
        }
        for key, value in m.items():
            if key.startswith("bragg_") and key not in row:
                row[key] = value
        metrics_rows.append(row)
        logger.info(
            "metrics: tilt=%g°  λ=%g μm  θ_max=%.2f°  F/B=%.3f  "
            "R/L_int=%.3f  asym=%.3f  HPBW=%.2f°%s",
            tilt, wl, m["theta_max_deg"], m["front_back_ratio"],
            m["right_left_integral_ratio"], m["right_left_asymmetry"],
            m["half_power_beamwidth_deg"],
            "  [truncated]" if m["hpbw_truncated"] else "",
        )
    metrics_df = pd.DataFrame(metrics_rows)

    # --- CSV 输出 ---
    long_rows = []
    for (tilt, wl), case in results.items():
        I = case["intensity"]
        I_norm = I / np.max(I)
        for theta, i_val, i_norm in zip(angles_deg, I, I_norm):
            long_rows.append({
                "period_um": args.period,
                "top_width_um": args.top_width,
                "bottom_width_um": args.bottom_width,
                "depth_um": args.depth,
                "nperiods": args.nperiods,
                "tilt_angle_deg": tilt,
                "wavelength_um": wl,
                "angle_deg": float(theta),
                "far_field_intensity": float(i_val),
                "normalized_intensity": float(i_norm),
            })
    long_df = pd.DataFrame(long_rows)

    tables_dir = ensure_dir(project_path("results", "tables"))
    emission_csv = tables_dir / "angle_resolved_emission.csv"
    metrics_csv = tables_dir / "angle_resolved_metrics.csv"
    long_df.to_csv(emission_csv, index=False)
    metrics_df.to_csv(metrics_csv, index=False)
    logger.info("emission CSV → %s  (rows=%d)", emission_csv, len(long_df))
    logger.info("metrics CSV  → %s  (rows=%d)", metrics_csv, len(metrics_df))

    # --- 图 ---
    title_root = (f"P={args.period:g}, top={args.top_width:g}, "
                  f"bot={args.bottom_width:g}, D={args.depth:g} μm, "
                  f"nperiods={args.nperiods}, res={args.resolution}/μm")
    cart_png = plot_cartesian(
        results, angles_deg, args.wavelengths, args.tilt_angles,
        project_path("results", "figures", "angle_resolved_cartesian.png"),
        title=f"Scattered far field: {title_root}",
    )
    polar_png = plot_polar(
        results, angles_deg, args.wavelengths, args.tilt_angles,
        project_path("results", "figures", "angle_resolved_polar.png"),
        title=f"Polar scattered field: {title_root}",
    )
    logger.info("cartesian PNG → %s", cart_png)
    logger.info("polar PNG     → %s", polar_png)

    # --- 总结 ---
    total = sum(c["walltime_s"] for c in results.values())
    logger.info("=== done; total Meep walltime = %.1f s ===", total)
    logger.info("self-check summary: %s", check_summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
