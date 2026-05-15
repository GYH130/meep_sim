"""Meep 仿真装配与运行。

职责
----
- 根据几何 + 材料 + 源 + 监视器配置构造 `meep.Simulation`；
- 提供常用源（宽带平面波、Gaussian 脉冲）和监视器（flux plane、
  近场 DFT）的封装；
- 运行至能量衰减阈值，并保存 raw 输出到 data/raw/。

设计原则
--------
- **不要** 在脚本里直接拼 Meep 对象；统一通过本模块的工厂函数，
  以便在 postprocess 中能用同一份元数据解释结果。
- 本模块只负责 “搭建并跑”，不做物理后处理（那是 postprocess.py 的事）。

物理假设
--------
- 2D 表面模型中 x 为周期方向，y 为表面法向；
- 默认仿真域沿 y 上下用 PML，x 方向使用周期/Bloch 边界；
- 入射方向默认从 +y 指向 -y（从上方打到结构上表面）；
- 角度入射通过设置 Bloch k 矢量实现（同一频率单角度 → 多频率需多次仿真，
  在脚本里循环 / MPI 并行处理）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class SimulationSetup:
    """一次仿真所需的全部参数（与具体 Meep 对象解耦，便于序列化）。

    Attributes
    ----------
    cell_size_um : tuple[float, float, float]
        Meep cell 尺寸 (x, y, z)，μm。2D 仿真时 z = 0。
    pml_thickness_um : float
        PML 厚度 (μm)，建议 >= 最大波长。
    resolution : int
        points per μm。
    lambda_min_um, lambda_max_um : float
        研究波段。
    n_freq : int
        flux 监视器采样频率数。
    incidence_theta_deg : float
        入射极角（从 z 轴起算），度。
    incidence_phi_deg : float
        入射方位角，度。
    polarization : str
        "TE" 或 "TM"（仅 2D / 平面入射有意义）。
    runtime_decay_db : float
        Meep `until_after_sources` 衰减阈值 (dB)。
    """

    cell_size_um: tuple[float, float, float]
    pml_thickness_um: float
    resolution: int
    lambda_min_um: float
    lambda_max_um: float
    n_freq: int
    incidence_theta_deg: float = 0.0
    incidence_phi_deg: float = 0.0
    polarization: str = "TE"
    runtime_decay_db: float = 30.0


def build_simulation(setup: SimulationSetup, geometry: list, sources: list,
                     boundary_layers: list | None = None) -> Any:
    """构造 `meep.Simulation` 对象（骨架占位）。"""
    raise NotImplementedError("Step 0 骨架：在首个具体仿真脚本中实现。")


def make_planewave_source(setup: SimulationSetup, z_position_um: float) -> list:
    """生成宽带平面波源（骨架占位）。

    参数包含中心频率 / 带宽 / Bloch-k 等，从 `setup` 推导。
    """
    raise NotImplementedError("Step 0 骨架。")


def add_flux_monitors(sim: Any, setup: SimulationSetup,
                      z_refl_um: float, z_trans_um: float) -> dict:
    """在反射面和透射面分别加 flux 监视器。

    Returns
    -------
    dict
        {"refl": meep.FluxRegion-like, "trans": ...}，
        交由 postprocess 模块根据 setup 中的 n_freq、波段恢复物理量。
    """
    raise NotImplementedError("Step 0 骨架。")


def run_until_decayed(sim: Any, setup: SimulationSetup) -> None:
    """按 `runtime_decay_db` 衰减阈值运行仿真，并打印进度日志。"""
    raise NotImplementedError("Step 0 骨架。")


# ---------------------------------------------------------------------------
# 通用 2D 周期金属表面 R/T/A 仿真
# ---------------------------------------------------------------------------
# 这是 01/02 脚本中 “参考 run + 结构 run” 流程的通用版本：把 geometry 通过
# 工厂函数注入，仿真域 (cell、PML、源、monitor) 布局保持一致，便于不同结构
# 直接对照。01_flat_ti_benchmark.py 是历史脚本，保留原实现；02_periodic_groove
# 可酌情在下次重构时切到本接口；03_slanted_groove_spectrum 起统一用本接口。


def run_periodic_2d_metal_spectrum(
    *,
    geometry_factory,
    period_um: float,
    wavelength_min_um: float,
    wavelength_max_um: float,
    resolution: int,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    air_buffer_um: float,
    nfreq: int,
    decay_db: float = 40.0,
    source_component: str = "Ez",
    logger: Any = None,
) -> dict:
    """2D 周期金属表面 R/T/A 谱仿真。

    Parameters
    ----------
    geometry_factory : Callable[[float, float], list]
        签名 ``f(y_surface_um, substrate_thickness_um) -> list[meep.GeometricObject]``。
        必须返回包含 “基底 + 表面结构” 的完整几何列表；本函数不会再额外加基底。
        把 y_surface 作为参数注入是为了让结构正确锚在仿真域内（cell 中心在原点）。
    period_um : float
        x 方向周期 (μm)，决定 cell 的 x 宽度。
    wavelength_min_um, wavelength_max_um : float
        研究波段端点 (μm)。
    resolution : int
        Meep points/μm。
    pml_thickness_um, substrate_thickness_um, air_buffer_um : float
        仿真域几何参数 (μm)。
    nfreq : int
        flux monitor 采样频率数。
    decay_db : float
        `stop_when_fields_decayed` 衰减阈值 (dB)。
    source_component : {"Ez", "Hz"}
        2D 中选 Ez (TE w.r.t. 光栅) 或 Hz (TM)。
    logger : logging.Logger | None
        若提供则用 info 输出进度。

    Returns
    -------
    dict
        {wavelength_um, R, T, A, walltime_s, geometry_y}，均为已按波长升序排序的
        ``numpy.ndarray``（标量项除外）。

    Raises
    ------
    ValueError
        参数非法；监视器 y 坐标顺序异常时 RuntimeError。

    Notes
    -----
    仿真域沿 y 自上而下：
        [+y] top PML / 空气 buffer / 源 / 反射 monitor /
             空气 / 表面 (y=y_surface) / 基底 / 透射 monitor / bottom PML [-y]
    监视器 y 坐标完全由 cell 尺寸推导，避免 magic number。
    """
    import time
    import meep as mp
    import numpy as np
    from .materials import freq_range_for_band, meep_freq_to_wavelength_um

    # ---- 参数自检 ----
    if period_um <= 0:
        raise ValueError(f"period_um 必须 > 0, 收到 {period_um}")
    if resolution <= 0 or nfreq <= 1:
        raise ValueError("resolution 必须 > 0，nfreq 必须 > 1。")
    if pml_thickness_um <= 0 or substrate_thickness_um <= 0 or air_buffer_um <= 0:
        raise ValueError("PML / substrate / air_buffer 厚度必须为正。")
    if wavelength_min_um <= 0 or wavelength_min_um >= wavelength_max_um:
        raise ValueError(
            f"波段非法: [{wavelength_min_um}, {wavelength_max_um}] μm"
        )

    component_map = {"Ez": mp.Ez, "Hz": mp.Hz}
    if source_component not in component_map:
        raise ValueError(
            f"source_component 必须是 {list(component_map)}, 收到 {source_component}"
        )
    src_c = component_map[source_component]

    log = logger.info if logger is not None else (lambda *a, **k: None)

    # ---- 频率 ----
    f_min, f_max, fcen = freq_range_for_band(wavelength_min_um, wavelength_max_um)
    df = f_max - f_min

    # ---- cell ----
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
    y_substrate_bottom = y_surface - substrate_thickness_um
    y_src = y_top_pml_inner - 0.25 * air_buffer_um
    y_refl = y_surface + 0.5 * air_buffer_um
    y_trans = 0.5 * (y_substrate_bottom + y_bottom_pml_inner)

    if not (y_surface < y_refl < y_src < y_top_pml_inner):
        raise RuntimeError("源/反射 monitor 的 y 坐标顺序异常。")
    if not (y_bottom_pml_inner < y_trans < y_substrate_bottom):
        raise RuntimeError("透射 monitor 未位于基底底面与 bottom PML 之间。")

    # ---- 几何 (由调用方注入) ----
    geometry = geometry_factory(y_surface, substrate_thickness_um)
    if not isinstance(geometry, list) or len(geometry) == 0:
        raise ValueError("geometry_factory 必须返回非空 list[meep.GeometricObject]。")

    pml_layers = [mp.PML(thickness=pml_thickness_um, direction=mp.Y)]
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=df, is_integrated=True),
            component=src_c,
            center=mp.Vector3(0, y_src, 0),
            size=mp.Vector3(period_um, 0, 0),
        )
    ]

    # ---- 参考 run (空气) ----
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
    log("[ref run] res=%d, cell=(%.2f, %.2f) μm, nfreq=%d, fcen=%.4f, df=%.4f, comp=%s",
        resolution, cell.x, cell.y, nfreq, fcen, df, source_component)
    sim_ref.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20, c=src_c,
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

    # ---- 结构 run ----
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
    log("[struct run] %d objects in geometry (注入自 geometry_factory)", len(geometry))
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20, c=src_c,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-decay_db / 10.0),
        )
    )
    refl_flux = np.array(mp.get_fluxes(refl_fr))
    trans_flux = np.array(mp.get_fluxes(trans_fr))
    t_struct = time.time() - t0
    log("[struct run] done in %.2f s", t_struct)

    # Meep 约定下，水平 flux 面的正方向为 +y；本项目源沿 -y 入射。
    # 因此透射到下方的功率 raw flux 为负，物理透射率需取相反号。
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


def run_periodic_2d_metal_single_wavelength(
    *,
    geometry_factory,
    period_um: float,
    wavelength_um: float,
    resolution: int,
    pml_thickness_um: float,
    substrate_thickness_um: float,
    air_buffer_um: float,
    decay_db: float,
    source_component: str,
    fwidth_fraction: float = 0.06,
    solver_version: str = "diagnostics_v2_signed_flux_single_wavelength",
    source_mode: str = "single_wavelength",
    logger: Any = None,
) -> dict:
    """Run one narrowband normal-incidence 2D periodic metal simulation.

    This diagnostics_v2 interface avoids broad Gaussian endpoint quality
    problems by running one target wavelength at a time.  The source frequency
    is exactly ``1 / wavelength_um`` and each flux monitor records only that
    frequency.  The reference and structure runs use identical source
    parameters.

    Returns scalar raw fluxes and D00-signed R/T/A values plus provenance fields.
    ``solver_version`` and ``source_mode`` are keyword-overridable so downstream
    diagnostics can record their own provenance labels while reusing the same
    verified flux layout.
    """
    import time
    import meep as mp
    import numpy as np
    from .postprocess import compute_RTA_downward_incidence

    if period_um <= 0:
        raise ValueError(f"period_um 必须 > 0, 收到 {period_um}")
    if wavelength_um <= 0:
        raise ValueError(f"wavelength_um 必须 > 0, 收到 {wavelength_um}")
    if resolution <= 0:
        raise ValueError("resolution 必须 > 0。")
    if pml_thickness_um <= 0 or substrate_thickness_um <= 0 or air_buffer_um <= 0:
        raise ValueError("PML / substrate / air_buffer 厚度必须为正。")
    if fwidth_fraction <= 0:
        raise ValueError("fwidth_fraction 必须为正。")

    component_map = {"Ez": mp.Ez, "Hz": mp.Hz}
    if source_component not in component_map:
        raise ValueError(
            f"source_component 必须是 {list(component_map)}, 收到 {source_component}"
        )
    src_c = component_map[source_component]
    log = logger.info if logger is not None else (lambda *a, **k: None)

    fcen = 1.0 / wavelength_um
    fwidth = fwidth_fraction * fcen

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
    y_substrate_bottom = y_surface - substrate_thickness_um
    y_src = y_top_pml_inner - 0.25 * air_buffer_um
    y_refl = y_surface + 0.5 * air_buffer_um
    y_trans = 0.5 * (y_substrate_bottom + y_bottom_pml_inner)
    if not (y_surface < y_refl < y_src < y_top_pml_inner):
        raise RuntimeError("源/反射 monitor 的 y 坐标顺序异常。")
    if not (y_bottom_pml_inner < y_trans < y_substrate_bottom):
        raise RuntimeError("透射 monitor 未位于基底底面与 bottom PML 之间。")

    geometry = geometry_factory(y_surface, substrate_thickness_um)
    if not isinstance(geometry, list) or len(geometry) == 0:
        raise ValueError("geometry_factory 必须返回非空 list[meep.GeometricObject]。")

    pml_layers = [mp.PML(thickness=pml_thickness_um, direction=mp.Y)]
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=fcen, fwidth=fwidth, is_integrated=True),
            component=src_c,
            center=mp.Vector3(0, y_src, 0),
            size=mp.Vector3(period_um, 0, 0),
        )
    ]

    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=resolution,
        k_point=mp.Vector3(),
        geometry=[],
        dimensions=2,
    )
    refl_ref = sim_ref.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(
            center=mp.Vector3(0, y_refl, 0),
            size=mp.Vector3(period_um, 0, 0),
        ),
    )
    t0 = time.time()
    log(
        "[single ref] wl=%.6g um, f=%.6g, res=%d, comp=%s",
        wavelength_um, fcen, resolution, source_component,
    )
    sim_ref.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-decay_db / 10.0),
        )
    )
    input_flux_raw = float(np.array(mp.get_fluxes(refl_ref))[0])
    ref_data = sim_ref.get_flux_data(refl_ref)
    t_ref = time.time() - t0
    if abs(input_flux_raw) == 0:
        raise RuntimeError("参考 run 中 input_flux_raw 为 0，无法归一化。")

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        sources=sources,
        resolution=resolution,
        k_point=mp.Vector3(),
        geometry=geometry,
        dimensions=2,
    )
    refl = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(
            center=mp.Vector3(0, y_refl, 0),
            size=mp.Vector3(period_um, 0, 0),
        ),
    )
    trans = sim.add_flux(
        fcen, 0, 1,
        mp.FluxRegion(
            center=mp.Vector3(0, y_trans, 0),
            size=mp.Vector3(period_um, 0, 0),
        ),
    )
    sim.load_minus_flux_data(refl, ref_data)
    t0 = time.time()
    log("[single struct] %d geometry objects", len(geometry))
    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            dt=20,
            c=src_c,
            pt=mp.Vector3(0, y_refl, 0),
            decay_by=10 ** (-decay_db / 10.0),
        )
    )
    reflection_flux_raw = float(np.array(mp.get_fluxes(refl))[0])
    transmission_flux_raw = float(np.array(mp.get_fluxes(trans))[0])
    t_struct = time.time() - t0

    rta = compute_RTA_downward_incidence(
        np.array([reflection_flux_raw]),
        np.array([transmission_flux_raw]),
        np.array([input_flux_raw]),
    )
    return dict(
        wavelength_um=float(wavelength_um),
        R=float(rta["R"][0]),
        T=float(rta["T"][0]),
        A=float(rta["A"][0]),
        input_flux_raw=input_flux_raw,
        reflection_flux_raw=reflection_flux_raw,
        transmission_flux_raw=transmission_flux_raw,
        raw_trans_flux=transmission_flux_raw,
        raw_transmittance=float(rta["raw_transmittance"][0]),
        signed_transmittance=float(rta["signed_transmittance"][0]),
        source_mode=source_mode,
        solver_version=solver_version,
        polarization=source_component,
        resolution=resolution,
        pml_thickness_um=pml_thickness_um,
        substrate_thickness_um=substrate_thickness_um,
        air_buffer_um=air_buffer_um,
        decay_db=decay_db,
        fwidth_fraction=fwidth_fraction,
        material_model="caller_supplied_geometry_materials",
        walltime_s=t_ref + t_struct,
        geometry_y=dict(
            cell_y=cell_y,
            y_surface=y_surface,
            y_substrate_bottom=y_substrate_bottom,
            y_src=y_src,
            y_refl=y_refl,
            y_trans=y_trans,
        ),
    )
