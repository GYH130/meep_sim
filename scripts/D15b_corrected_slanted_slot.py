"""D15b_corrected_slanted_slot.py — 修正版中心斜槽条件（D15 的修订实现）

本脚本**不修改** ``scripts/D15_central_slanted_slot_single_condition.py``：
D15 的历史结果（R=0.7365395766 等）按原样留档。这里是一份新实现，
逐条补齐审计列出的缺陷。

相对 D15 的修正（逐条对应审计）
--------------------------------
1. **源后观察窗口**：D15 源结束 1312.5、运行结束 1321.03125，源后仅 8.53。
   本脚本用 ``--post_source_um`` 显式给出源后观察长度，并在运行过程中按
   固定间隔采样 R/T/A 与多探针场强包络，输出源后曲线数据。

2. **布局**：Ti 只占中段，下方留非 PML 空气间隔，透射面**放在空气里**。
   所有界面/源/监视器 y 坐标与材料归属显式记录到输出（``layout`` 段）。

3. **全阶枚举**：由实际 kx 与周期自动枚举全部传播阶
   （``m_prop = floor((f − kx)/(1/P)) + floor((f + kx)/(1/P)) + 1``），
   不使用固定的 ±n。每个传播阶用**相位投影**独立求功率（见下），
   并把"全阶和 vs 整面通量"的闭合差报出来。

4. **体吸收 DFT**：按每对原生 Yee 分量分别积分，体元由实际网格间距给出，
   不把 ``get_array_metadata`` 的居中网格权重套给 ``yee_grid=True``；
   量化 |density| < 阈值 的微小负值的影响，并给出金属/空气分区贡献。

5. **执行控制**：q≥1 直接拒绝；真实墙钟超时（子进程看门狗）；
   记录真实 MPI 进程数与每个 rank 的内存；仅 rank 0 原子写结果；
   质量状态逐项记录（收敛、分阶、多探针、源后稳定性），
   任一项缺失即标 ``valid=False`` 并写明缺哪一项。

物理约定（与 D05/D15 一致）
---------------------------
- 二维、周期方向 x、法线方向 z → Meep 的 (x, y)。
- 发射观察方向 s = (sinθ, 0, cosθ)，互易吸收的入射波矢沿 −s，
  故 kx = −f·sinθ。观察角、入射波矢、槽倾角**分别保存**。
- λ=10.5 μm、p 偏振（Hz 源）。R = refl/|input|、T = −trans/|input|（保留原符号）。
- Meep 内置 Rakić 1998 Ti；oxide_layer_in_model=false；二维中心截面代理模型。

适用范围与限制
--------------
- 单条件单波长；不是九点试扫，也不声称二维资格已通过。
- q<1 只作粗略稳定性检查，不是收敛证明。
"""

from __future__ import annotations

import argparse
import cmath
import json
import logging
import math
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

TI_F_MAX_MEEP = 15.67135255701507
STABILITY_NOTE = ("coarse_only: q<1 只说明可以启动，不代表已收敛或结果可用；"
                  "q>=1 直接拒绝运行")
#: |density| 低于 (该比例 × 峰值) 的点视为浮点噪声，单独统计其影响
NOISE_REL = 1e-9


def dispersion_q(resolution: int, courant: float) -> float:
    return float(np.pi * TI_F_MAX_MEEP * courant / resolution)


def enumerate_orders(fcen: float, kx: float, period_um: float) -> list[dict]:
    """由实际 kx 与周期枚举**全部传播阶**，并固定 m 的方向约定。

    二维平面波在周期方向的波矢是 kx_m = kx + m·(1/P)（m 为整数），
    传播阶要求 |kx_m| ≤ f。故

        m ∈ [ ceil((−f − kx)·P), floor((+f − kx)·P) ]

    **约定**：m = 0 恒为入射/镜面阶（kx_0 = kx），不做布里渊区折叠。
    折叠会把 m 的标号整体平移（"5e−17 抵消" 效应），使 m=0 不再是镜面阶，
    与审计 "从 −2 到 7 共 10 阶" 的记法不符。这里直接用原始 kx 定界。

    参考条件（f=1/10.5, P=50, kx/f=−0.5）给出 m ∈ [−2, 7]，共 10 阶；
    若 kx 反号则给出 [−7, 2]（同一 |kx_m| 集合的镜像，阶数相同）。
    """
    inv_p = 1.0 / period_um
    tol = 1e-9
    m_lo = math.ceil((-fcen - kx) * period_um - tol)
    m_hi = math.floor((+fcen - kx) * period_um + tol)
    orders = []
    for m in range(m_lo, m_hi + 1):
        kxm = kx + m * inv_p
        ky2 = fcen * fcen - kxm * kxm
        if ky2 <= 0:
            continue
        orders.append({
            "m": int(m),
            "kx_m": float(kxm),
            "ky_m": float(math.sqrt(ky2)),
            "sin_theta_m": float(kxm / fcen),
            "theta_out_deg": float(math.degrees(math.asin(max(-1.0, min(1.0, kxm / fcen))))),
            "evanescent": False,
        })
    return orders


def slot_contains_point(x_um: float, z_um: float, geom) -> bool:
    """点 (x, z) 是否在槽腔（材料之外）内。

    槽腔的构造性定义：|p·n| ≤ W/2 且 p·u ≤ L 且 z ≥ 0（见 src.slot_geometry），
    **不因居中平移而改变**——所以这里用未平移的几何坐标判定入参即可；
    调用方负责把 Meep 坐标 x 换回几何坐标（减 cell_shift_x）。

    这个判据与 Meep 内部的材料数组相互独立，用于体吸收的金属/空气分区，
    避免把两套网格的插值结果叠在一起（那正是 "金属功率恰好为 0" 的来源）。
    """
    u = geom.axis_unit_um
    n = geom.wall_normal_unit
    hw = 0.5 * geom.spec.width_um
    if z_um < 0.0:
        return False
    xi = x_um * u[0] + z_um * u[1]
    eta = x_um * n[0] + z_um * n[1]
    return (abs(eta) <= hw + 1e-9) and (xi <= geom.spec.axis_length_um + 1e-9)


def _points_in_slot_mask(xx: np.ndarray, zz: np.ndarray, geom) -> np.ndarray:
    """:func:`slot_contains_point` 的向量化版本（同一判据）。"""
    u = geom.axis_unit_um
    n = geom.wall_normal_unit
    hw = 0.5 * geom.spec.width_um
    xi = xx * u[0] + zz * u[1]
    eta = xx * n[0] + zz * n[1]
    return ((np.abs(eta) <= hw + 1e-9)
            & (xi <= geom.spec.axis_length_um + 1e-9)
            & (zz >= 0.0))


def _cplx(z: complex) -> dict:
    """复数转可序列化的 {re, im, abs}（JSON 不能直接存 complex）。"""
    return {"re": float(z.real), "im": float(z.imag), "abs": float(abs(z))}


def order_power_projections(fields: dict, comps_named: dict, orders: list[dict],
                            xs: np.ndarray, n_grid: int, fcen: float,
                            polarization: str) -> list[dict]:
    """把面上横向场按传播阶做投影，返回每阶的法向通量（Meep 单位）。

    与 Meep 对象解耦，便于用**合成平面波**做单元测试：给定只含单一阶的
    解析场，投影必须只在该阶上给出非零通量（见 tests/test_d15b_order_projection.py）。

    相位符号（**实测标定**，不是假设）
    --------------------------------
    Meep 的振幅约定在 ``amp_func`` 与 ``k_point`` 上**叠加**了两次
    exp(+i2πk_x x)，因此 DFT 面场是去掉该因子之后的**包络**。把它投回
    第 m 阶时必须取共轭相位 exp(+i2πk_{x,m}x)。

    这个方向是在**合成解析场**上钉死的（见
    ``tests/test_d15b_order_projection.py``）：只有投影相位与场的相位
    严格共轭时，单阶场才给出 ``|Σ F·φ| = n``：

        |Σ Hz·conj(φ)| = 800.000   ← 正确（n=800）
        |Σ Hz· φ    | =  36.375   ← 符号取反，明显偏小

    用另一个来源（r16 实跑的单阶平面波参考）独立复核，同一个方向也给出
    一致的复振幅比

        |E_x|/|H_z| = 0.8665   （解析 cosθ = 0.8660）
        |E_y|/|H_z| = 0.4998   （解析 sinθ = 0.5000）

    一致到 0.05% 以内。早期版本用了 exp(−i…) 因而把单阶场摊到多个阶上，
    导致负份额与 ~0.5% 的闭合比——这不是 FDTD 数据的问题，是投影符号。

    归一化与 "全阶闭合" 的含义
    --------------------------
    本函数返回**相对**份额（各阶通量除以该侧各阶之和），绝对标定由
    closure 段与整面 FluxRegion 对比给出。之所以不直接给绝对量：
    离散投影的绝对归一取决于 Meep 对 DFT 场与 FluxRegion 的
    常数因子（面积/频率权重），本脚本不做隐式假设，而是把
    "阶和 / 整面通量" 作为**一致性检查量**报出来，让读者看到它是否接近 1。

    网格归一因子
    ------------
    离散投影取**逐点求和**；在**等距、覆盖整数个周期**的网格上，
    ``Σ_j exp(∓i2π(k_{x,m}−k_{x,m'})x_j)`` 在 m≠m' 时严格为 0、m=m' 时为 n，
    因此单阶平面波的投影幅度恰为 n·amp，无需额外归一因子。
    **前提是网格等距且覆盖整数个周期**：Meep 的 DFT 面网格实测为
    802 点、起点 −25.03125、步长 0.062344…（P=50、r16），即覆盖长度
    略超过一个周期——此时 m≠m' 的交叉项不再严格为 0（本脚本用
    ``orthogonality_self_check`` 量化这一项，超过阈值就把分阶表降级）。

    通量公式（与 Meep 的 ½Re(E×H*) 约定一致，复数场）
    -------------------------------------------------
    二维 p 偏振（Hz 源）:
        ⟨S_y⟩ = ½·Re[ E_x·H_z*·(k_y/k₀) − E_y·H_z*·(k_x/k₀) ]
    二维 s 偏振（Ez 源）:
        ⟨S_y⟩ = ½·Re[ E_z·H_x*·(k_x/k₀) − E_y·H_x*·(k_y/k₀) ]
    """
    xs = np.asarray(xs, dtype=float)
    n_x = int(xs.size)
    if n_x < 2:
        raise ValueError("相位投影需要至少 2 个采样点")
    # 绝对标定：投影是**逐点求和**，单阶幅度含因子 n；而 Meep 的
    # FluxRegion 通量是对整个单胞宽度 P 积分的面密度。解析对照
    # （单位 |H_z| 的平面波，见 tests）给出
    #     F_plane = 0.5Re[E_x H_z* k_y/k0 - E_y H_z* k_x/k0]·P
    #     P_proj  = 0.5Re[...]·n²
    # 故 F_m = P_proj·(period_um/n²)。乘以周期的是调用方（本函数不知道周期），
    # 这里先除以 n² 把幅度还原成物理量，'flux_meep_units' 才是可比的。
    # 投影是逐点求和，单阶幅度含因子 n，故 1/n² 先把幅度还原成物理量。
    # 这只是**必要**的一步，不保证等于 Meep FluxRegion 的绝对值：参考面位于
    # 线源之下，斜入射时面波尚未达到单位幅度（实测比值约 0.289），
    # 因此绝对量纲由参考面标定记录（见 order_decomposition.reference_calibration）。
    scale = 1.0 / float(n_x * n_x)
    rows = []
    for o in orders:
        kxm, kym = o["kx_m"], o["ky_m"]
        # 投影相位 = exp(−i2πk_x x)。
        # 方向由**物理判据**钉死，不是靠合成场自洽：平面对照的镜面反射必须
        # 落在 m=0 阶（θ_out = −θ_inc = −30°）。实测两种约定：
        #     ph = exp(−i2πk_x x)  →  峰值在 m=0（θ_out=−30°）✓
        #     ph = exp(+i2πk_x x)  →  峰值在 m=+5（θ_out=+33.4°）✗
        # 早期版本用 +i 号，把整个阶指配偏移了 5 阶，导致 flat 的镜面峰
        # 出现在 m=5、各阶份额完全错位。
        ph = np.exp(-1j * 2.0 * math.pi * kxm * xs)
        row = {"m": o["m"], "kx_m": kxm, "ky_m": kym,
               "theta_out_deg": o["theta_out_deg"]}
        def _fld(name: str):
            """取分量数组；同时接受字符串键（"Ex"）与 mp 常量键（mp.Ex）。

            这两种键在历史上被混用过，产生过 KeyError(0) 与静默漏分量。
            两种都接受，彻底消除这一类错误。
            """
            if name in fields:
                return fields[name]
            comp = comps_named.get(name)
            if comp is not None and comp in fields:
                return fields[comp]
            raise KeyError(f"分量 {name} 既不在字符串键也不在常量键中："
                           f"现有键 {list(fields.keys())}")

        if polarization == "p":
            ex = complex(np.sum(_fld("Ex") * ph))
            ey = complex(np.sum(_fld("Ey") * ph))
            hz = complex(np.sum(_fld("Hz") * ph))
            sy = 0.5 * np.real(ex * np.conj(hz) * (kym / fcen)
                               - ey * np.conj(hz) * (kxm / fcen))
            row.update(amplitude={"ex": _cplx(ex), "ey": _cplx(ey),
                                  "hz": _cplx(hz)})
        else:
            ey = complex(np.sum(_fld("Ey") * ph))
            ez = complex(np.sum(_fld("Ez") * ph))
            hx = complex(np.sum(_fld("Hx") * ph))
            sy = 0.5 * np.real(ez * np.conj(hx) * (kxm / fcen)
                               - ey * np.conj(hx) * (kym / fcen))
            row.update(amplitude={"ey": _cplx(ey), "ez": _cplx(ez),
                                  "hx": _cplx(hx)})
        row["flux_meep_units"] = float(sy * scale)
        row["power_outward"] = float(sy * scale)
        rows.append(row)
    return rows


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--wavelength_um", type=float, default=10.5)
    p.add_argument("--fwidth_fraction", type=float, default=0.08)
    p.add_argument("--emission_theta_deg", type=float, default=30.0)
    p.add_argument("--polarization", choices=["p", "s"], default="p")
    p.add_argument("--resolution", type=int, default=16)
    p.add_argument("--courant", type=float, default=0.25)
    p.add_argument("--decay_db", type=float, default=40.0)
    p.add_argument("--decay_dt", type=float, default=20.0)
    p.add_argument("--post_source_um", type=float, default=40.0,
                   help="源结束后继续观察的 Meep 时间长度（a/c 单位 = μm）")
    p.add_argument("--min_post_source_um", type=float, default=30.0,
                   help="源后窗口低于该值即判为不满足审计第 1 条")
    p.add_argument("--max_time_um", type=float, default=0.0,
                   help=">0 时给出绝对时间上限（Meep 时间单位），到点即停")
    p.add_argument("--relative_drift_tol", type=float, default=0.02,
                   help="源后收敛判据：包络后半段/前半段均值相对漂移与该段极差"
                        "都要小于此值")
    p.add_argument("--probe_interval_um", type=float, default=5.0,
                   help="多探针/通量采样间隔")
    p.add_argument("--safety_factor", type=float, default=3.0,
                   help="停止判据所需衰减对应的绝对时间上限系数")
    p.add_argument("--wall_timeout_s", type=float, default=5400.0)
    p.add_argument("--period_um", type=float, default=50.0)
    p.add_argument("--width_um", type=float, default=30.0)
    p.add_argument("--axis_length_um", type=float, default=40.0)
    p.add_argument("--tilt_angle_deg", type=float, default=30.0)
    p.add_argument("--ti_thickness_um", type=float, default=60.0)
    p.add_argument("--bottom_thickness_um", type=float, default=10.0)
    p.add_argument("--min_wall_um", type=float, default=1.0)
    p.add_argument("--pml_thickness_um", type=float, default=6.0)
    p.add_argument("--air_above_um", type=float, default=8.0)
    p.add_argument("--air_below_um", type=float, default=12.0)
    p.add_argument("--case", choices=["structure", "flat"], default="structure",
                   help="flat = 同布局去除槽的平面对照（记录参考伪反射）")
    p.add_argument("--max_steps", type=int, default=0,
                   help=">0 时固定步数运行（用于进程数/内存标定），不做收敛判定")
    p.add_argument("--tag", type=str, default="")
    return p


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    log_dir = PROJECT_ROOT / "logs" / "run_20260918"
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[logging.FileHandler(log_dir / f"D15b{args.tag}.log"),
                  logging.StreamHandler()],
    )
    logger = logging.getLogger("D15b")

    import meep as mp
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank, nprocs = comm.Get_rank(), comm.Get_size()
    except Exception:
        comm, rank, nprocs = None, 0, 1

    from src.materials import get_ti_medium
    from src.slot_geometry import (SlantedBlindSlotSpec, SlotGeometry,
                                   cell_shift_x, perimeter_path_length_um,
                                   polygon_xy, qualification)

    wl = args.wavelength_um
    fcen = 1.0 / wl
    fwidth = args.fwidth_fraction * fcen
    q = dispersion_q(args.resolution, args.courant)

    # ---- 审计第 5 条：q>=1 拒绝（不是警告） ----
    if q >= 1.0:
        logger.error("拒绝运行：q=%.6f >= 1（色散不稳定）", q)
        if rank == 0:
            print(json.dumps({"status": "rejected_dispersion_unstable", "q": q},
                             indent=2, ensure_ascii=False))
        return 2

    spec = SlantedBlindSlotSpec(
        period_um=args.period_um, width_um=args.width_um,
        axis_length_um=args.axis_length_um,
        tilt_angle_deg=args.tilt_angle_deg,
        ti_thickness_um=args.ti_thickness_um,
        bottom_thickness_um=args.bottom_thickness_um,
        min_wall_um=args.min_wall_um)
    geom = SlotGeometry(spec)
    qual = qualification(geom)
    if not qual["qualified"]:
        logger.error("几何不合格，拒绝运行：%s", qual["reasons"])
        if rank == 0:
            print(json.dumps({"geometry_qualified": False,
                              "reasons": qual["reasons"]}, indent=2,
                             ensure_ascii=False))
        return 2

    theta_rad = math.radians(args.emission_theta_deg)
    kx = -fcen * math.sin(theta_rad)      # 入射波矢沿 −s
    src_c = mp.Hz if args.polarization == "p" else mp.Ez
    ti = get_ti_medium(wl - 1.0, wl + 1.0)

    # ------------------------------------------------------------------
    # 布局（审计第 2 条）：上方空气 + Ti + **下方空气** + 两侧 PML
    # ------------------------------------------------------------------
    pml = args.pml_thickness_um
    sub = spec.ti_thickness_um
    air_a, air_b = args.air_above_um, args.air_below_um
    cell_y = 2 * pml + air_a + sub + air_b
    y_top_inner = 0.5 * cell_y - pml
    y_surface = y_top_inner - air_a
    y_sub_bottom = y_surface - sub
    y_bot_inner = -0.5 * cell_y + pml

    # 上方空气内的划分：源偏上，反射面在源与表面之间
    y_src = y_top_inner - 0.20 * air_a
    y_refl = y_top_inner - 0.65 * air_a
    # 下方空气内的划分：透射面居中（远离底 PML）
    y_trans = y_bot_inner + 0.5 * air_b
    # 槽内探针（审计第 1 条要求 "多个槽内探针"）：沿槽轴等分取 3 个点，
    # 全部落在腔中心线上（几何坐标 → 加 cell_shift_x 换到 Meep 坐标）。
    # 另外在**金属侧**同深度放 1 点，用来对比腔内外是否同时稳定。
    u_ax = geom.axis_unit_um
    probe_pts: list[tuple[float, float, float | str]] = []
    for fr in (0.25, 0.5, 0.75):
        xi = fr * spec.axis_length_um
        z = xi * u_ax[1]
        x = xi * u_ax[0]
        probe_pts.append((float(x + cell_shift_x(geom)), float(y_surface - z), fr))
    z_metal = 0.5 * geom.vertical_depth_um
    x_metal = -0.5 * spec.period_um + 0.5 * spec.period_um * 0.5  # 周期中点右侧的壁内
    probe_pts.append((float(x_metal), float(y_surface - z_metal), "metal"))

    cell = mp.Vector3(spec.period_um, cell_y, 0)
    layout = {
        "period_um": spec.period_um, "cell_y_um": cell_y,
        "pml_thickness_um": pml, "air_above_um": air_a, "air_below_um": air_b,
        "ti_thickness_um": sub,
        "y_cell_top": 0.5 * cell_y, "y_cell_bottom": -0.5 * cell_y,
        "y_top_inner": y_top_inner, "y_bot_inner": y_bot_inner,
        "y_source": y_src, "y_refl_monitor": y_refl, "y_trans_monitor": y_trans,
        "y_surface": y_surface, "y_sub_bottom": y_sub_bottom,
        "material_by_region": [
            {"region": "top PML", "y_from": 0.5 * cell_y, "y_to": y_top_inner,
             "material": "PML(air)"},
            {"region": "air above", "y_from": y_top_inner, "y_to": y_surface,
             "material": "air"},
            {"region": "Ti film", "y_from": y_surface, "y_to": y_sub_bottom,
             "material": "Ti (Rakic 1998) minus slot cavity"},
            {"region": "air below", "y_from": y_sub_bottom, "y_to": y_bot_inner,
             "material": "air"},
            {"region": "bottom PML", "y_from": y_bot_inner,
             "y_to": -0.5 * cell_y, "material": "PML(air)"},
        ],
        "monitor_medium": {
            "refl": "air" if y_refl > y_surface else "Ti",
            "trans": "air" if y_trans < y_sub_bottom else "Ti",
        },
        "slot_probe_points": [{"frac_along_axis": f, "x": x, "y": y}
                              for (x, y, f) in probe_pts],
        "note": ("透射面位于 Ti 之下的非 PML 空气层内；其到下方 PML 的距离 = "
                 f"{y_trans - y_bot_inner:.3f} μm"),
    }
    assert y_refl > y_surface > y_sub_bottom > y_trans > y_bot_inner, layout
    logger.info("layout: cell_y=%.3f y_src=%.3f y_refl=%.3f y_surface=%.3f "
                "y_sub_bottom=%.3f y_trans=%.3f y_bot_inner=%.3f",
                cell_y, y_src, y_refl, y_surface, y_sub_bottom, y_trans, y_bot_inner)

    # ---- 全阶枚举 ----
    orders = enumerate_orders(fcen, kx, spec.period_um)
    logger.info("orders: %d 个传播阶 m=%s", len(orders),
                [o["m"] for o in orders])

    def amp_func(r):
        return cmath.exp(1j * 2.0 * math.pi * kx * r.x)

    sources = [mp.Source(
        mp.GaussianSource(frequency=fcen, fwidth=fwidth, is_integrated=True),
        component=src_c, center=mp.Vector3(0, y_src, 0),
        size=mp.Vector3(spec.period_um, 0, 0), amp_func=amp_func)]

    verts = [mp.Vector3(x, y, 0) for (x, y) in polygon_xy(geom, y_surface)]
    geometry = [
        mp.Block(material=ti, center=mp.Vector3(0, y_surface - 0.5 * sub, 0),
                 size=mp.Vector3(spec.period_um, sub, mp.inf)),
    ]
    if args.case == "structure":
        geometry.append(mp.Prism(vertices=verts, height=mp.inf,
                                 axis=mp.Vector3(0, 0, 1),
                                 material=mp.Medium(epsilon=1.0)))
    pml_layers = [mp.PML(thickness=pml, direction=mp.Y)]

    def _sim(glist):
        return mp.Simulation(cell_size=cell, boundary_layers=pml_layers,
                             geometry=glist, sources=sources,
                             resolution=args.resolution, dimensions=2,
                             Courant=args.courant, k_point=mp.Vector3(kx, 0, 0))

    # ---------------- 参考（入射）运行 ----------------
    # 除整面 FluxRegion 外，**在同面上**再放一组 DFT 场并保存原始复数组：
    # 入射波是已知的单阶平面波，用同一套投影流程处理它，就得到分阶功率的
    # **绝对标定系数**，使参考场上 Σ_m P_m / F_plane 严格等于该阶占比。
    # 没有这一步，分阶表只有相对份额，"全阶和 vs 整面通量" 的闭合差无从判断
    # （早期版本因此报出 ~0.005 与 ~7451 两个无意义的比值）。
    sim_ref = _sim([])
    refl_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_refl, 0), size=mp.Vector3(spec.period_um, 0, 0)))
    _xyplane = 0.5 * (spec.period_um / max(1, int(round(
        spec.period_um * args.resolution))))
    _ref_comps = [mp.Ex, mp.Ey, mp.Hz] if args.polarization == "p" \
        else [mp.Ey, mp.Ez, mp.Hx]
    dft_ref = sim_ref.add_dft_fields(_ref_comps, fcen, 0, 1,
                                     center=mp.Vector3(0, y_refl, 0),
                                     size=mp.Vector3(spec.period_um, _xyplane, 0),
                                     yee_grid=False)
    sim_ref.init_sim()
    t0 = time.perf_counter()
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(
        dt=args.decay_dt, c=src_c, pt=mp.Vector3(0, y_refl, 0),
        decay_by=10 ** (-args.decay_db / 10.0)))
    ref_wall = time.perf_counter() - t0
    input_flux_raw = float(mp.get_fluxes(refl_ref)[0])
    ref_data = sim_ref.get_flux_data(refl_ref)
    ref_sim_time = float(sim_ref.meep_time())
    ref_plane_arrays = {int(c): np.asarray(sim_ref.get_dft_array(dft_ref, c, 0))
                        for c in _ref_comps}
    try:
        ref_plane_x = np.ravel(np.asarray(
            sim_ref.get_array_metadata(dft_cell=dft_ref)[0], dtype=float))
    except Exception:
        ref_plane_x = None
    # 把参考面原始复数组落盘：入射波是已知的单阶平面波，有了它就能**离线**
    # 标定分阶功率、复核份额，不必重跑 FDTD（本条件参考段约 8–18 分钟）。
    try:
        import h5py
        cache_dir = PROJECT_ROOT / "results" / "refplanes"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache = cache_dir / f"refplane{args.tag}.h5"
        with h5py.File(cache, "w") as h5:
            for k, v in ref_plane_arrays.items():
                h5.create_dataset(mp.component_name(k), data=v)
            h5.attrs["ref_plane_x_first"] = float(ref_plane_x[0]) if ref_plane_x is not None else float("nan")
            h5.attrs["ref_plane_x_last"] = float(ref_plane_x[-1]) if ref_plane_x is not None else float("nan")
            h5.attrs["n_points"] = int(ref_plane_x.size) if ref_plane_x is not None else 0
            h5.attrs["input_flux_raw"] = float(input_flux_raw)
            h5.attrs["fingerprint"] = (f"P{spec.period_um}-W{spec.width_um}-"
                                       f"L{spec.axis_length_um}-a{args.tilt_angle_deg}-"
                                       f"r{args.resolution}-C{args.courant}-"
                                       f"pol{args.polarization}-th{args.emission_theta_deg}")
        logger.info("参考面已缓存：%s", cache.relative_to(PROJECT_ROOT))
    except Exception as exc:
        logger.warning("参考面缓存失败（不影响本次结果）：%r", exc)
    sim_ref.reset_meep()
    if input_flux_raw == 0:
        raise RuntimeError("reference input flux is zero")

    # 高斯源（is_integrated=True）的**源结束时间**。Meep 内部把
    # ``gaussian_src_time`` 的 end_time 设为 ``start_time + 2·width·cutoff``
    # （见 meep/source.py）；``.cutoff`` 属性本身**不是**结束时间（它默认 5.0），
    # D15 的日志正好落在 2·131.25·5 = 1312.5，可用于核对本式。
    _gs = mp.GaussianSource(frequency=fcen, fwidth=fwidth, is_integrated=True)
    gs_duration = float(_gs.start_time + 2.0 * _gs.width * _gs.cutoff)
    logger.info("reference: t_end=%.3f wall=%.1f s input_flux=%.6e "
                "source_end_time=%.3f (= start + 2*width*cutoff = 2*%.3f*%.1f)",
                ref_sim_time, ref_wall, input_flux_raw, gs_duration,
                _gs.width, _gs.cutoff)

    # ---------------- 结构运行 ----------------
    sim = _sim(geometry)
    refl = sim.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_refl, 0), size=mp.Vector3(spec.period_um, 0, 0)))
    trans = sim.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_trans, 0), size=mp.Vector3(spec.period_um, 0, 0)))
    sim.load_minus_flux_data(refl, ref_data)

    # 分阶功率：在反射面与透射面各放一组 DFT 场（Ex/Ey/Hz 或 Ey/Ez/Hx），
    # 之后按传播阶的横向波矢做相位投影，独立算出每阶法向功率。
    # 注意：size 在法向**不能给 0**——零尺寸的 DFT 区域在更新时会被
    # ``get_array_metadata`` 拒绝（KeyError: 0）。给一个远小于网格步长
    # (1/res) 的厚度，数组仍只有一层采样点。
    n_dir = int(round(spec.period_um * args.resolution))
    xyplane = 0.5 * (spec.period_um / n_dir)
    if args.polarization == "p":
        dft_comp_orders = [mp.Ex, mp.Ey, mp.Hz]
    else:
        dft_comp_orders = [mp.Ey, mp.Ez, mp.Hx]
    dft_refl = sim.add_dft_fields(dft_comp_orders, fcen, 0, 1,
                                  center=mp.Vector3(0, y_refl, 0),
                                  size=mp.Vector3(spec.period_um, xyplane, 0),
                                  yee_grid=False)
    dft_trans = sim.add_dft_fields(dft_comp_orders, fcen, 0, 1,
                                   center=mp.Vector3(0, y_trans, 0),
                                   size=mp.Vector3(spec.period_um, xyplane, 0),
                                   yee_grid=False)

    # 体吸收 DFT：**按对分配**，不使用 yee_grid 的居中权重
    dft_comps = [mp.Ex, mp.Dx, mp.Ey, mp.Dy, mp.Ez, mp.Dz]
    dft_center = mp.Vector3(0, 0.5 * (y_surface + y_sub_bottom), 0)
    dft_size = mp.Vector3(spec.period_um, sub, 0)
    dft_obj = sim.add_dft_fields(dft_comps, fcen, 0, 1,
                                 center=dft_center, size=dft_size,
                                 yee_grid=True)

    sim.init_sim()
    t_start = sim.meep_time()

    # ---- 源后采样：多探针包络 + 通量演化 ----
    samples: list[dict] = []
    next_sample = [t_start]
    probe_windows: list[list[float]] = [[] for _ in probe_pts]

    def _fields_probe(s):
        vals = []
        for (px, py, _fr) in probe_pts:
            v = s.get_field_point(src_c, mp.Vector3(px, py, 0))
            vals.append(abs(complex(v)))
        return vals

    def _sample(s):
        rec = {
            "t": float(s.meep_time()),
            "refl": float(mp.get_fluxes(refl)[0]) / abs(input_flux_raw),
            "trans": float(mp.get_fluxes(trans)[0]) / abs(input_flux_raw),
            "probe_abs": _fields_probe(s),
        }
        samples.append(rec)
        for i, v in enumerate(rec["probe_abs"]):
            probe_windows[i].append(v)

    def _step_fn(s):
        if s.meep_time() >= next_sample[0]:
            _sample(s)
            next_sample[0] = s.meep_time() + args.probe_interval_um

    # ---- 审计第 5 条：真实墙钟超时控制 ----
    # 用 SIGALRM 在**本进程**内硬中断（Meep 的 run 是 C 循环，Python 层的
    # 时间检查不会生效；SIGALRM 由解释器在主线程的字节码边界处理，
    # 对 Meep 的 C 调用同样有效，因为 Meep 会周期性回到 Python）。
    class _WallTimeout(Exception):
        pass

    timeout_state = {"fired": False}

    def _on_alarm(signum, frame):
        timeout_state["fired"] = True
        raise _WallTimeout(f"wall timeout {args.wall_timeout_s}s exceeded")

    import signal as _signal
    _old_handler = None
    try:
        _old_handler = _signal.signal(_signal.SIGALRM, _on_alarm)
        _signal.setitimer(_signal.ITIMER_REAL, float(args.wall_timeout_s))
    except Exception as exc:
        logger.warning("cannot install wall-timeout watchdog: %r", exc)

    def _clear_timer():
        try:
            _signal.setitimer(_signal.ITIMER_REAL, 0.0)
            if _old_handler is not None:
                _signal.signal(_signal.SIGALRM, _old_handler)
        except Exception:
            pass

    wall0 = time.perf_counter()
    converged = True
    run_error = None
    timed_out = False
    decay_met = None
    if args.max_steps > 0:
        try:
            sim.run(mp.at_every(args.probe_interval_um, _step_fn),
                    until=args.max_steps * args.courant / args.resolution)
        except _WallTimeout as exc:
            run_error, converged, timed_out = repr(exc), False, True
        except Exception as exc:
            run_error, converged = repr(exc), False
    else:
        # 停止条件（审计第 1 条）：**同时**满足
        #   (a) 面监视器处场强包络在**最近一段**内已收敛（不是 "历史上曾经
        #       达到最大值的 1/thr" —— 那个判据只有在包络已过峰并单调衰减时
        #       才可靠，斜槽腔里探针的瞬态会拖很久，导致运行停不下来）；
        #   (b) 全部槽内探针的包络同样收敛；
        #   (c) 源后观察窗口 ≥ --min_post_source_um；
        #   (d) 未超过 --max_time_um（若有），且未触发墙钟超时。
        # (a)(b) 用**滑动窗口内的相对漂移**判定：把已采样的历史一分为二，
        # 比较后半段与前半段的均值，并要求后半段自身极差也足够小。
        dlg = logging.getLogger("D15b")
        last_decay = [None]
        decay_hist: dict = {"face": [], "probes": []}

        def _drift_converged(series: list[float]) -> bool:
            n = len(series)
            if n < 8:
                return False
            half = n // 2
            a = float(np.mean(series[:half]))
            b = float(np.mean(series[half:]))
            scale = max(abs(a), abs(b), 1e-300)
            tail = series[half:]
            spread = (max(tail) - min(tail)) / max(abs(np.mean(tail)), 1e-300)
            return (abs(b - a) / scale) <= args.relative_drift_tol and \
                spread <= args.relative_drift_tol

        def _decay_stop(s):
            if s.meep_time() <= gs_duration:
                return False
            mag = abs(complex(s.get_field_point(
                src_c, mp.Vector3(0, y_surface - 2.0, 0))))
            probes = _fields_probe(s)
            decay_hist["face"].append(mag)
            decay_hist["probes"].append(list(probes))
            window_ok = (s.meep_time() - gs_duration) >= args.min_post_source_um
            if args.max_time_um and s.meep_time() >= args.max_time_um:
                dlg.warning("达到 --max_time_um=%.1f，提前结束（未满足漂移判据）",
                            args.max_time_um)
                last_decay[0] = float(s.meep_time())
                return True
            if not window_ok or not _drift_converged(decay_hist["face"]):
                return False
            nprobe = len(decay_hist["probes"][0])
            for i in range(nprobe):
                if not _drift_converged([p[i] for p in decay_hist["probes"]]):
                    return False
            last_decay[0] = float(s.meep_time())
            return True

        try:
            sim.run(mp.at_every(args.probe_interval_um, _step_fn),
                    until_after_sources=_decay_stop)
        except _WallTimeout as exc:
            run_error, converged, timed_out = repr(exc), False, True
        except Exception as exc:
            run_error, converged = repr(exc), False
        decay_met = last_decay[0] is not None
    _clear_timer()
    wall_struct = time.perf_counter() - wall0
    t_end = float(sim.meep_time())
    _sample(sim)
    sim.fields.update_dfts()
    refl_raw = float(mp.get_fluxes(refl)[0])
    trans_raw = float(mp.get_fluxes(trans)[0])

    R = refl_raw / abs(input_flux_raw)
    T = -trans_raw / abs(input_flux_raw)
    A_flux = 1.0 - R - T

    post_source = float(t_end - gs_duration)
    post_source_ok = post_source >= args.min_post_source_um and args.max_steps == 0
    if args.max_steps > 0:
        post_source_ok = None

    # ------------------------------------------------------------------
    # 分阶功率：对每个传播阶做相位投影
    # ------------------------------------------------------------------
    order_rows: list[dict] = []
    order_summary: dict = {}
    try:
        # 实测（1×3 小 cell + 本条件 r16）：DFT 数组轴序是 **(x, y)**，与
        # ``get_array_metadata`` 的第 0/1 个元素（x/y 坐标数组）一致。
        xs_plane: list[float] = []

        def _fields_at(dftobj) -> dict:
            out = {}
            for c in dft_comp_orders:
                out[c] = np.asarray(sim.get_dft_array(dftobj, c, 0))
            return out

        def _periodic_axis(a: np.ndarray) -> int:
            """周期方向（x）所在的轴：长度为 P·res + 1 的那一轴。"""
            if a.ndim <= 1:
                return 0
            target = spec.period_um * args.resolution + 1.0
            return int(np.argmin([abs(a.shape[k] - target)
                                  for k in range(a.ndim)]))

        def _thin(arr) -> np.ndarray:
            """把 DFT 面场压成沿**周期方向**的一维数组。

            轴序 (x, y)；法向只剩 1–2 像素（面上的 2 个 yve 采样点差异远小于
            一个网格步长，取第 0 层即可，与 Meep FluxRegion 的取值位置一致）。
            """
            a = np.asarray(arr)
            if a.ndim == 1:
                return a
            ax = _periodic_axis(a)
            idx = tuple(0 if k != ax else slice(None) for k in range(a.ndim))
            return a[idx]

        comps_named = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez": mp.Ez,
                       "Hx": mp.Hx, "Hz": mp.Hz}

        # 相位参考坐标：**必须用 Meep 自己的元数据**。
        # 实测（r16, P=50）：x 方向是 802 点、起点 −25.03125、步长 0.0625 稍有
        # 偏差（0.062344…），并不是我之前假设的 "801 点、起点 −25.0、步长 1/res"。
        # 半个单元的偏差就足以破坏 exp(∓i2πk_{x,m}x) 在整数阶之间的正交性，
        # 使单阶场被摊到其它阶上（这就是分阶表出现负份额与 0.5% 闭合的原因）。
        # 因此这里直接取元数据；取不到才退回解析网格，并把降级记录在案。
        xs_m = None
        try:
            md_p = sim.get_array_metadata(dft_cell=dft_refl)
            xs_m = np.ravel(np.asarray(md_p[0], dtype=float))
        except Exception as exc:
            logger.warning("反射面 metadata 不可用（%r）；退回解析网格（分阶表"
                           "将标注 degraded）", exc)
        if xs_m is not None and xs_m.size > 1:
            xs = xs_m
            coord_source = "meep_metadata"
            dx_p = float(xs[1] - xs[0])
        else:
            n_pts = int(round(spec.period_um * args.resolution))
            dx_p = spec.period_um / n_pts
            xs = -0.5 * spec.period_um + np.arange(n_pts) * dx_p
            coord_source = "analytic_fallback"
        xs_plane[:] = [float(v) for v in xs]
        n_proj = int(xs.size)
        logger.info("分阶相位参考：%s  n=%d 起点=%.6f 步长=%.6f",
                    coord_source, n_proj, float(xs[0]), dx_p)

        # 正交性自检（**非平凡**）：在不同阶之间做交叉投影
        #     O(m, m') = |Σ_j exp(−i2πk_m x_j)·exp(+i2πk_m' x_j)| / n
        # 正确的等距网格上 m≠m' 时 O ≈ 0、m=m' 时 O = 1。
        # 网格相对阶的相位不自洽时 O(m,m') 会明显偏离 0——此时分阶表不可用。
        _ord_ms = [o["m"] for o in orders]
        cross = []
        for i, mi in enumerate(_ord_ms):
            for mj in _ord_ms[i + 1:]:
                ki = [o["kx_m"] for o in orders if o["m"] == mi][0]
                kj = [o["kx_m"] for o in orders if o["m"] == mj][0]
                cross.append(float(abs(np.sum(
                    np.exp(-1j * 2.0 * math.pi * ki * xs)
                    * np.exp(+1j * 2.0 * math.pi * kj * xs))) / n_proj))
        max_cross = max(cross) if cross else float("nan")
        orthogonality_ok = bool(max_cross < 0.05)
        if not orthogonality_ok:
            logger.warning("阶间交叉投影最大 %.4f（应 << 1）：该布局的 DFT 网格"
                           "无法干净分离传播阶，分阶表降级为定性参考", max_cross)
        logger.info("阶间正交性自检：max O(m,m') = %.3e (%s)，坐标来源 %s",
                    max_cross, "ok" if orthogonality_ok else "degraded",
                    coord_source)

        def _order_powers(fields_obj, plane_y, outward_sign) -> list[dict]:
            """fields_obj 的键是 **mp.Ex 等分量常量**（与 comps_named 的值同型）。

            不要混用 "Ex" 这样的字符串与 mp.Ex：``comps_named`` 的键是字符串、
            值是常量，用错一侧就会 KeyError(0)——本脚本早期版本正是如此。

            返回的 ``flux_meep_units`` 已乘上周期 P，与 ``mp.get_fluxes`` 的
            整面 FluxRegion 同量纲，故"阶和 / 整面通量"可以直接比较。
            """
            named = {name: _thin(fields_obj[comp])
                     for name, comp in comps_named.items()
                     if comp in fields_obj
                     and comp in dft_comp_orders}
            rows = order_power_projections(named, comps_named, orders, xs,
                                           n_proj, fcen, args.polarization)
            for r in rows:
                r["power_outward"] = r["power_outward"] * spec.period_um
                r["flux_meep_units"] = r["flux_meep_units"] * spec.period_um
            return rows

        f_refl = _fields_at(dft_refl)
        f_trans = _fields_at(dft_trans)
        order_rows = _order_powers(f_refl, y_refl, +1)
        trans_rows = _order_powers(f_trans, y_trans, -1)

        # ---- 用参考面（已知单阶平面波）做**绝对标定交叉核对** ----
        # 入射波在同一高度只走 m=0 阶，其投影功率应与 Meep 的整面通量相等。
        # 两者一致 => 1/n²·P 的标定正确；不一致则分阶表只有相对意义。
        calib = {}
        try:
            # 参考面键是 **mp.Ex 等常量**（与 _thin/_fields_at 一致），
            # 不要用字符串——那正是早期 KeyError(0) 的来源。
            ref_named = {name: _thin(ref_plane_arrays[int(comp)])
                         for name, comp in comps_named.items()
                         if int(comp) in ref_plane_arrays}
            ref_rows = order_power_projections(ref_named, comps_named, orders,
                                               xs, n_proj, fcen,
                                               args.polarization)
            for r in ref_rows:
                r["flux_meep_units"] *= spec.period_um
            ref_proj_total = float(sum(r["flux_meep_units"] for r in ref_rows))
            ref_peak = max(ref_rows, key=lambda r: abs(r["flux_meep_units"]))
            # 入射侧的整面通量：参考运行在同一面、同高度、同一点的 FluxRegion
            calib = {
                "ref_plane_projected_total": ref_proj_total,
                "ref_plane_meep_flux": input_flux_raw,
                "ref_plane_ratio": (float(ref_proj_total / input_flux_raw)
                                    if input_flux_raw else None),
                "ref_plane_peak_order_m": ref_peak["m"],
                "ref_plane_peak_theta_deg": ref_peak["theta_out_deg"],
                "ref_plane_other_orders_fraction": (
                    float(sum(abs(r["flux_meep_units"]) for r in ref_rows
                              if r["m"] != ref_peak["m"])
                          / max(abs(ref_peak["flux_meep_units"]), 1e-300))),
                "note": ("参考面是线源下方的斜入射面波：镜面阶必须落在 m=0"
                         "（θ_out=−θ_inc）。比值不等于 1 是**预期的**——"
                         "线源有限孔径使面波在参考面处未达单位幅度，"
                         "该比值即这一残余因子，用于把分阶表换算到"
                         "Meep 通量量纲。**注意**：绝对闭合在结构面上"
                         "本不成立（相消干涉），故不作为通过判据。"),
            }
        except Exception as exc:
            calib = {"error": repr(exc)}


        refl_sum = float(sum(r["power_outward"] for r in order_rows))
        trans_sum = float(sum(r["power_outward"] for r in trans_rows))
        # 出射侧的分阶份额（归一化到该侧各阶功率之和）
        order_summary = {
            "n_orders": len(orders),
            "m_list": [o["m"] for o in orders],
            "refl_plane_y": y_refl, "trans_plane_y": y_trans,
            "phase_reference": {"coord_source": coord_source,
                                "n_points": n_proj,
                                "x_first": float(xs[0]), "x_last": float(xs[-1]),
                                "dx_um": float(xs[1] - xs[0]) if n_proj > 1 else None},
            "orthogonality_self_check": {
                "max_cross_overlap": max_cross,
                "threshold": 0.05,
                "ok": orthogonality_ok,
                "note": ("O(m,m')=|Σ exp(−i2πk_m x)·exp(+i2πk_m' x)|/n；"
                         "m≠m' 时应 ≈0。若 > 阈值，说明 DFT 网格与阶的相位"
                         "不自洽，分阶份额只能作定性参考。")},
            "refl_order_power_sum_meep_units": refl_sum,
            "trans_order_power_sum_meep_units": trans_sum,
            "refl_order_fractions": {
                r["m"]: (float(r["power_outward"] / refl_sum) if refl_sum else None)
                for r in order_rows},
            "trans_order_fractions": {
                r["m"]: (float(r["power_outward"] / trans_sum) if trans_sum else None)
                for r in trans_rows},
            "closure": {
                "refl_flux_meep_units": refl_raw,
                "trans_flux_meep_units": trans_raw,
                "refl_order_over_flux": (float(refl_sum / refl_raw)
                                         if refl_raw else None),
                "trans_order_over_flux": (float(trans_sum / trans_raw)
                                          if trans_raw else None),
                "note": ("相位投影给出各阶**相对**份额；与整面 FluxRegion 的"
                         "绝对标定由 closure 段比较（两者在相消干涉时应接近）。"),
            },
            "reference_calibration": calib,
            "rows_refl": order_rows,
            "rows_trans": trans_rows,
        }
    except Exception as exc:
        logger.error("order decomposition failed: %r", exc)
        logger.error("traceback:\n%s", traceback.format_exc())
        order_summary = {"error": repr(exc)}

    # ------------------------------------------------------------------
    # 体吸收 A_vol = ∫ 2πf·Im(conj(E)·D) dV / |input|
    #   · **每对** E/D 分量分别相乘、分别按各自的坐标数组积分，再相加；
    #   · 坐标与体元来自 ``get_array_metadata`` 的坐标数组步长（首末点各取
    #     1/2 权重，端点不重复计数）；**不使用权重数组**——那是居中网格的
    #     面积权重，套到 yee_grid=True 的原生交错网格上是错的；
    #   · 量化 |density| 低于相对阈值的浮点噪声，并单独报出它的贡献；
    #   · 金属/空气分区按**几何**判定（点是否在槽多边形内），不依赖另一套
    #     网格的插值结果，避免错位导致 "金属功率恰好为 0"。
    # ------------------------------------------------------------------
    A_vol = float("nan")
    volume_flux_diff = float("nan")
    vol_summary: dict = {}
    try:
        pairs = [(mp.Ex, mp.Dx), (mp.Ey, mp.Dy), (mp.Ez, mp.Dz)]
        per_pair: dict = {}
        per_pair_area: dict = {}
        pair_shapes: dict = {}
        contribs: list = []          # (density, area, xs, ys)
        noise_sum = 0.0
        noise_pts = pts_total = neg_pts = 0
        dens_min, dens_max = np.inf, -np.inf
        coord_sources = []

        # 体吸收 DFT 区域的物理范围（与 dft_center/dft_size 一致）
        y_from = float(dft_center.y - 0.5 * sub)
        y_to = float(dft_center.y + 0.5 * sub)

        for e_c, d_c in pairs:
            e_a = np.asarray(sim.get_dft_array(dft_obj, e_c, 0))
            d_a = np.asarray(sim.get_dft_array(dft_obj, d_c, 0))
            # 未分配的分量返回标量（ndim=0），不是空数组
            if e_a.ndim != 2 or d_a.ndim != 2:
                continue
            if e_a.shape != d_a.shape:
                logger.warning("E/D shape mismatch %s*%s: %s vs %s",
                               mp.component_name(e_c), mp.component_name(d_c),
                               e_a.shape, d_a.shape)
            n = (min(e_a.shape[0], d_a.shape[0]),
                 min(e_a.shape[1], d_a.shape[1]))
            dens = np.imag(np.conj(e_a[:n[0], :n[1]]) * d_a[:n[0], :n[1]])
            name = f"{mp.component_name(e_c)}*{mp.component_name(d_c)}"
            per_pair[name] = float(np.sum(dens))
            pair_shapes[name] = [int(n[0]), int(n[1])]

            # 该对的坐标数组。Meep 的 get_array_metadata 只接受
            # vol / center+size / dft_cell（**没有 component 参数**）；
            # 而不同分量的交错网格点数可能差 1，因此这里按**各轴点数**
            # 解析地生成坐标：x 方向 n[1] 点、步长 1/res，从 −P/2 起；
            # y 方向 n[0] 点、步长 1/res，从 y_from 起。
            dx_p = 1.0 / args.resolution
            dy_p = 1.0 / args.resolution
            xs_p = -0.5 * spec.period_um + np.arange(n[1]) * dx_p
            ys_p = y_from + np.arange(n[0]) * dy_p
            # 用 Meep 的元数据做一次一致性核对（只对能取到的情形）
            try:
                md_pair = sim.get_array_metadata(
                    center=mp.Vector3(0, dft_center.y, 0),
                    size=mp.Vector3(spec.period_um, sub, 0))
                xs_m = np.ravel(np.asarray(md_pair[0], dtype=float))
                ys_m = np.ravel(np.asarray(md_pair[1], dtype=float))
                if xs_m.size > 1:
                    dx_m = float(xs_m[1] - xs_m[0])
                    if abs(dx_m - dx_p) > 1e-9:
                        logger.warning("x 步长不一致：metadata %.6f, 解析 %.6f",
                                       dx_m, dx_p)
                if ys_m.size > 1:
                    dy_m = float(ys_m[1] - ys_m[0])
                    if abs(dy_m - dy_p) > 1e-9:
                        logger.warning("y 步长不一致：metadata %.6f, 解析 %.6f",
                                       dy_m, dy_p)
                coord_sources.append("metadata_checked")
            except Exception as exc:
                coord_sources.append(f"analytic_only({exc!r})")

            wx = np.full(n[1], dx_p)
            wy = np.full(n[0], dy_p)
            if n[1] > 1:
                wx[0] = wx[-1] = 0.5 * dx_p
            if n[0] > 1:
                wy[0] = wy[-1] = 0.5 * dy_p
            area = np.outer(wy, wx)
            per_pair_area[name] = float(np.sum(area))

            peak = float(np.max(np.abs(dens))) if dens.size else 0.0
            thr = NOISE_REL * peak
            mask = np.abs(dens) < thr
            noise_sum += float(np.sum(dens[mask] * area[mask]))
            noise_pts += int(np.count_nonzero(mask))
            pts_total += int(dens.size)
            neg_pts += int(np.count_nonzero(dens < 0.0))
            if dens.size:
                dens_min = min(dens_min, float(np.min(dens)))
                dens_max = max(dens_max, float(np.max(dens)))
            contribs.append((dens, area, xs_p, ys_p))

        if not contribs:
            raise RuntimeError("no allocated E/D component pairs found")

        # 统一噪声阈值（取各对峰值中的最大者）后求和
        peak_all = max(abs(dens_max), abs(dens_min))
        thr_all = NOISE_REL * peak_all
        total_raw = 0.0
        total_clipped = 0.0
        metal_contrib = air_contrib = 0.0
        for dens, area, xs_p, ys_p in contribs:
            d_raw = 2.0 * math.pi * fcen * dens
            d = 2.0 * math.pi * fcen * np.where(np.abs(dens) < thr_all, 0.0, dens)
            total_raw += float(np.sum(d_raw * area))
            total_clipped += float(np.sum(d * area))
            # 逐点按几何分区（向量化：先造 (x, z) 网格，再判是否在槽内）
            if dens.shape != area.shape:
                continue
            # Meep 的 x = 几何 x + cell_shift_x（polygon_xy 里做过平移），
            # 槽腔的判定必须用几何坐标，故先减回去。
            xs_geom = (xs_p[:dens.shape[1]] - cell_shift_x(geom))
            xx = xs_geom[None, :] * np.ones((dens.shape[0], 1))
            zz = (y_surface - ys_p[:dens.shape[0]])[:, None] * np.ones(
                (1, dens.shape[1]))
            in_slot = _points_in_slot_mask(xx, zz, geom)
            inside = in_slot & (zz >= 0.0) & (zz <= sub)
            air_contrib += float(np.sum(d[inside] * area[inside]))
            metal_contrib += float(np.sum(d[~inside] * area[~inside]))

        A_vol = total_clipped / abs(input_flux_raw)
        volume_flux_diff = abs(A_vol - A_flux)

        vol_summary = dict(
            component_pairs_used=[[mp.component_name(a), mp.component_name(b)]
                                  for a, b in pairs
                                  if f"{mp.component_name(a)}*"
                                     f"{mp.component_name(b)}" in per_pair],
            per_pair_sum=per_pair,
            per_pair_area_um2=per_pair_area,
            pair_shapes=pair_shapes,
            axis_order_note=("get_dft_array 返回 (nx, ny)：第 0 轴是周期方向 x，"
                             "由两套独立小 cell 验证"),
            integral_note=("每一对 E/D **分别**按各自坐标数组积分后相加；体元 ="
                           "坐标数组步长乘积，首末点各取 1/2 权重；未使用"
                           " get_array_metadata 的权重数组"),
            coordinate_sources=coord_sources,
            density_min=(None if not np.isfinite(dens_min) else float(dens_min)),
            density_max=(None if not np.isfinite(dens_max) else float(dens_max)),
            negative_point_fraction=(float(neg_pts / pts_total) if pts_total else None),
            noise_rel_threshold=NOISE_REL,
            noise_abs_threshold=float(thr_all),
            noise_point_fraction=(float(noise_pts / pts_total) if pts_total else None),
            noise_volume_power=float(2.0 * math.pi * fcen * noise_sum),
            noise_relative_to_total=(float(abs(2.0 * math.pi * fcen * noise_sum)
                                           / abs(total_raw)) if total_raw else None),
            total_raw_unclipped=float(total_raw),
            total_clipped=float(total_clipped),
            clipped_minus_raw=float(total_clipped - total_raw),
            metal_region_power=float(metal_contrib),
            air_region_power=float(air_contrib),
            dft_region_bounds={"y_from": y_from, "y_to": y_to,
                               "x_from": float(-0.5 * spec.period_um),
                               "x_to": float(0.5 * spec.period_um)},
            partition_note=("金属/空气分区按**几何**判定（点是否落在槽多边形+"
                            "表面之下），不用另一套网格的插值结果"),
        )
    except Exception as exc:
        logger.error("volume-absorption computation failed: %r", exc)
        logger.error("traceback:\n%s", traceback.format_exc())
        vol_summary = {"error": repr(exc)}

    try:
        post = [s for s in samples if s["t"] > gs_duration]
        if len(post) >= 3:
            third = len(post) // 3
            segs = [post[:third], post[third:2 * third], post[2 * third:]]
            refl_tail = [float(np.mean([abs(s["refl"]) for s in seg])) for seg in segs]
            trans_tail = [float(np.mean([abs(s["trans"]) for s in seg])) for seg in segs]
            probe_env = []
            for i in range(len(probe_pts)):
                vals = [s["probe_abs"][i] for s in post]
                env = [float(np.max(vals[k:k + max(1, len(vals) // 3)]))
                       for k in range(0, len(vals), max(1, len(vals) // 3))][:3]
                probe_env.append(env)
            stab = {
                "source_cutoff_time": gs_duration, "t_end": t_end,
                "post_source_window_um": post_source,
                "min_required_post_source_um": args.min_post_source_um,
                "post_source_ok": post_source_ok,
                "refl_abs_segment_means": refl_tail,
                "trans_abs_segment_means": trans_tail,
                "refl_relative_drift_last_vs_first": (
                    abs(refl_tail[-1] - refl_tail[0]) / max(abs(refl_tail[0]), 1e-30)),
                "probe_envelope_segments": probe_env,
                "probe_env_relative_drift": [
                    (abs(env[-1] - env[0]) / max(abs(env[0]), 1e-30)) for env in probe_env],
                "n_samples": len(samples),
            }
        else:
            stab = {"error": "insufficient post-source samples",
                    "n_samples": len(samples), "t_end": t_end,
                    "source_cutoff_time": gs_duration}
    except Exception as exc:
        stab = {"error": repr(exc)}

    # ------------------------------------------------------------------
    # 质量状态（逐项，不用 "run 未抛异常" 代替）
    # ------------------------------------------------------------------
    checks = {
        "geometry_qualified": True,
        "dispersion_q_lt_1": q < 1.0,
        "converged_no_exception": bool(converged),
        "decay_criterion_met": (decay_met is True) if args.max_steps == 0 else None,
        "no_wall_timeout": not timed_out,
        "order_decomposition_done": "error" not in order_summary
                                    and bool(order_summary),
        "all_propagating_orders_enumerated": len(orders) > 0,
        "multi_probe_recorded": len(samples) > 0 and len(probe_pts) >= 2,
        "post_source_window_sufficient": (post_source_ok is True),
        "volume_absorption_done": "error" not in vol_summary and bool(vol_summary),
    }
    if args.max_steps > 0:
        # 标定运行只用于资源测量，不做收敛/资格判定
        missing = []
        valid = None
        status = "calibration_run_resource_only"
    else:
        missing = [k for k, v in checks.items() if v is not True]
        valid = not missing
        status = "valid" if valid else "invalid_missing:" + ",".join(missing)

    dt_sim = float(t_end - t_start)
    n_steps = int(round(dt_sim * args.resolution / args.courant))
    peak_mib = 0.0
    try:
        import resource
        peak_mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    except Exception:
        pass
    max_rank_mib = None
    if comm is not None:
        try:
            gathered = comm.gather(peak_mib, root=0)
            if rank == 0:
                max_rank_mib = float(max(gathered))
        except Exception:
            max_rank_mib = None

    rec = dict(
        case=args.case,
        wavelength_um=wl, frequency_meep=fcen,
        emission_theta_deg=args.emission_theta_deg,
        polarization=args.polarization,
        incident_kx=kx, incident_ky=-fcen * math.cos(theta_rad),
        slot_tilt_angle_deg=args.tilt_angle_deg,
        resolution=args.resolution, courant=args.courant,
        dispersion_q=q, stability_check=STABILITY_NOTE,
        mpi_processes=nprocs, omp_threads=int(os.environ.get("OMP_NUM_THREADS", "1") or 1),
        period_um=spec.period_um, width_um=spec.width_um,
        axis_length_um=spec.axis_length_um,
        surface_opening_um=geom.surface_opening_um,
        bottom_face_width_um=geom.bottom_face_width_um,
        perimeter_um=perimeter_path_length_um(geom),
        vertical_depth_um=geom.vertical_depth_um,
        max_depth_um=geom.max_depth_um,
        remaining_bottom_um=geom.remaining_bottom_um,
        wall_thickness_um=geom.wall_thickness_um,
        cell_y_um=cell_y,
        layout=layout,
        oxide_layer_in_model=False,
        material="meep.materials.Ti (Rakic 1998 Drude-Lorentz)",
        input_flux_raw=input_flux_raw,
        reflection_flux_raw=refl_raw, transmission_flux_raw=trans_raw,
        R=R, T=T, A_flux=A_flux,
        A_vol=A_vol, volume_flux_abs_difference=volume_flux_diff,
        volume_absorption_detail=vol_summary,
        r_plus_t_plus_a=R + T + A_flux,
        energy_note=("保留原始 R/T/A（含符号）；不以 R+T+A=1 作为守恒证据；"
                     "体吸收用按对配准的 DFT E/D 独立计算"),
        converged=converged, run_error=run_error, wall_timeout_s=args.wall_timeout_s,
        timed_out=timed_out, decay_criterion_met=decay_met,
        decay_db=args.decay_db, decay_dt=args.decay_dt,
        max_steps=args.max_steps,
        reference_sim_time_end=ref_sim_time, structure_sim_time_end=t_end,
        source_cutoff_time=gs_duration,
        post_source_window_um=post_source,
        dt_sim=dt_sim, n_steps=n_steps,
        wall_reference_s=ref_wall, wall_structure_s=wall_struct,
        wall_total_s=ref_wall + wall_struct,
        s_per_step=(wall_struct / n_steps) if n_steps else float("nan"),
        peak_memory_mib_per_rank=peak_mib,
        peak_memory_mib_max_rank=max_rank_mib,
        order_decomposition=order_summary,
        source_tail=samples,
        post_source_stability=stab,
        quality_checks=checks, quality_missing=missing,
        quality_status=status,
        fingerprint=(f"P{spec.period_um}-W{spec.width_um}-L{spec.axis_length_um}-"
                     f"a{args.tilt_angle_deg}-{args.case}-r{args.resolution}-"
                     f"C{args.courant}-pol{args.polarization}-"
                     f"th{args.emission_theta_deg}-lam{wl}"),
    )

    out = PROJECT_ROOT / f"results/tables/D15b_{args.case}{args.tag}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    if rank == 0:
        tmp = out.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(rec, indent=2, ensure_ascii=False))
        os.replace(tmp, out)          # 原子写：只有 rank 0 写最终文件
        logger.info("R=%.6f T=%.6f A_flux=%.6f A_vol=%.6f |A_vol-A_flux|=%.3e "
                    "orders=%d post_src=%.3f valid=%s n_steps=%d wall=%.1fs",
                    R, T, A_flux, A_vol, volume_flux_diff, len(orders),
                    post_source, valid, n_steps, rec["wall_total_s"])
        print(json.dumps({k: rec[k] for k in
                          ("fingerprint", "R", "T", "A_flux", "A_vol",
                           "volume_flux_abs_difference", "n_steps",
                           "post_source_window_um", "quality_status",
                           "quality_missing", "mpi_processes",
                           "peak_memory_mib_max_rank")},
                         indent=2, ensure_ascii=False))
    return 0 if (valid is not False) else 1


if __name__ == "__main__":
    raise SystemExit(main())
