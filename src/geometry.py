"""微结构几何定义模块。

职责
----
- 根据参数生成 Meep `geometry` 列表（光栅、微柱、凹坑、多层等）；
- 计算与几何耦合的仿真元数据：周期、最高特征频率（用于 resolution 估计）、
  PML/缓冲层厚度推荐值。

单位
----
所有几何长度统一使用 **μm**（与 Meep 单位 a = 1 μm 一致）。

约定
----
- 表面法向取 +z；
- 周期方向若为 1D 光栅默认沿 x；2D 光栅沿 x、y；
- 基底从 z = 0 向下延伸（z < 0），结构特征位于 z >= 0。

注意
----
本模块只描述几何形状，**不包含材料**。材料在 `materials.py` 中定义，
并在 `simulation.py` 中与几何组合成 Meep `Block`/`Cylinder` 等对象。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class GeometrySpec:
    """通用几何参数描述（与具体生成函数解耦）。

    Attributes
    ----------
    kind : str
        几何类型，例如 "flat" / "grating_1d" / "pillar_2d" / "hole_2d"。
    period_x, period_y : float
        周期 (μm)；非周期方向填 0。
    feature_size : float
        最小特征尺寸 (μm)，用于估算所需 Meep resolution。
    height : float
        结构高度 (μm)。
    extra : dict
        各类型特有参数（占空比、半径、倒角等）。
    """

    kind: str
    period_x: float
    period_y: float
    feature_size: float
    height: float
    extra: dict[str, Any]


def build_geometry(spec: GeometrySpec, material_substrate: Any, material_feature: Any) -> list:
    """根据 GeometrySpec 构建 Meep 几何对象列表。

    Parameters
    ----------
    spec : GeometrySpec
        几何参数。
    material_substrate : meep.Medium
        基底材料（来自 `materials.build_medium`）。
    material_feature : meep.Medium
        结构特征部分的材料（同上）。可与基底相同。

    Returns
    -------
    list[meep.GeometricObject]
        Meep 可直接使用的几何列表。

    Raises
    ------
    ValueError
        当 spec.kind 未识别或参数非法时。
    NotImplementedError
        未来骨架占位；具体生成函数将在后续脚本中按需实现。
    """
    raise NotImplementedError("Step 0 骨架：在后续脚本中按需实现具体几何生成。")


def estimate_resolution(spec: GeometrySpec, wavelength_min_um: float,
                        points_per_wavelength: int = 20,
                        points_per_feature: int = 10) -> int:
    """根据最短波长和最小特征给出建议 resolution（points per μm）。

    取两者要求的较大值，保证既能分辨电磁波长，又能分辨几何细节。
    返回值是整数 points/μm，可直接用于 `meep.Simulation(resolution=...)`。
    """
    if wavelength_min_um <= 0 or spec.feature_size <= 0:
        raise ValueError("wavelength_min_um 与 feature_size 必须为正数。")
    res_wave = points_per_wavelength / wavelength_min_um
    res_feat = points_per_feature / spec.feature_size
    return int(max(res_wave, res_feat) + 0.5)


# ---------------------------------------------------------------------------
# 具体几何工厂
# ---------------------------------------------------------------------------
# 所有 build_* 函数返回 list[meep.GeometricObject]，约定：
#   - 单胞沿 x 居中于原点 (x ∈ [-period_x/2, +period_x/2])；
#   - 基底上表面位于 y = y_surface，基底向 -y 延伸 substrate_thickness_um；
#   - Meep 后加入的几何覆盖先加入的，因此 “基底实体 + 在其上方/内部加 air” 即可
#     表示挖槽 / 挖孔。
# 用 mp.Prism 表示任意四边形：axis 沿 z (out-of-plane) → 2D 仿真中等效为
# 沿 z 方向无穷长的棱柱。


def build_rectangular_groove_geometry(
    *,
    period_x_um: float,
    groove_width_um: float,
    groove_depth_um: float,
    substrate_thickness_um: float,
    y_surface: float,
    medium_substrate: Any,
    medium_groove: Any = None,
) -> list:
    """金属基底 + 顶部对称矩形槽。

    Parameters
    ----------
    medium_groove : meep.Medium, optional
        槽内介质，默认空气 (epsilon=1)。

    Notes
    -----
    groove_width_um == 0 或 groove_depth_um == 0 时退化为无结构平板基底。
    主要用于：
      (1) `02_periodic_groove_spectrum.py` 之外的脚本复用；
      (2) `03_slanted_groove_spectrum.py` 的退化对照（tilt=0, top=bottom 时
          斜槽应与本函数构造的矩形槽给出相同光谱）。
    """
    import meep as mp

    if period_x_um <= 0:
        raise ValueError(f"period_x_um 必须 > 0, 收到 {period_x_um}")
    if groove_width_um < 0 or groove_depth_um < 0:
        raise ValueError("groove_width_um / groove_depth_um 必须 >= 0")
    if groove_width_um >= period_x_um:
        raise ValueError(
            f"groove_width_um ({groove_width_um}) 必须 < period_x_um ({period_x_um})"
        )
    if groove_depth_um >= substrate_thickness_um:
        raise ValueError(
            f"groove_depth_um ({groove_depth_um}) 必须 < substrate_thickness_um "
            f"({substrate_thickness_um})"
        )
    if medium_groove is None:
        medium_groove = mp.Medium(epsilon=1.0)

    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_x_um, substrate_thickness_um, mp.inf),
        )
    ]
    if groove_width_um > 0 and groove_depth_um > 0:
        geom.append(
            mp.Block(
                material=medium_groove,
                center=mp.Vector3(0, y_surface - groove_depth_um / 2.0, 0),
                size=mp.Vector3(groove_width_um, groove_depth_um, mp.inf),
            )
        )
    return geom


def slanted_groove_vertices(
    *,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    y_surface: float = 0.0,
) -> list[tuple[float, float]]:
    """返回斜槽轮廓的四个顶点 (xy 平面，按逆时针)。

    约定
    ----
    - 顶部边在 y = y_surface, 宽 top_width_um, 中心 x = 0；
    - 底部边在 y = y_surface - depth_um, 宽 bottom_width_um,
      中心 x = depth_um · tan(tilt_angle_deg) （即沿 +x 平移）。

    退化检查
    --------
    - tilt_angle_deg == 0 且 top_width == bottom_width → 对称矩形；
    - tilt_angle_deg == 0 且 top_width ≠ bottom_width → 对称梯形（V / 倒 V）；
    - tilt_angle_deg ≠ 0 → 顶底中心错开，整体是非对称四边形 ——
      这正是真实斜烧蚀槽产生的“表面非对称”几何含义，而**不是**把
      对称矩形整体旋转（那只是坐标系变换，物理上仍对称）。

    Returns
    -------
    list of (x, y) tuples，长度 4，顺序：
      top-left, top-right, bottom-right, bottom-left（逆时针）。
    """
    import math

    if top_width_um < 0 or bottom_width_um < 0:
        raise ValueError("top_width_um / bottom_width_um 必须 >= 0")
    if depth_um <= 0:
        raise ValueError(f"depth_um 必须 > 0, 收到 {depth_um}")
    if abs(tilt_angle_deg) >= 89.0:
        raise ValueError(
            f"|tilt_angle_deg| = {abs(tilt_angle_deg)} 必须 < 89°；"
            f"过大时槽几乎水平，不再代表 “斜槽” 的物理含义。"
        )

    dx = depth_um * math.tan(math.radians(tilt_angle_deg))
    y_top = y_surface
    y_bot = y_surface - depth_um
    return [
        (-top_width_um / 2.0, y_top),
        (+top_width_um / 2.0, y_top),
        (+bottom_width_um / 2.0 + dx, y_bot),
        (-bottom_width_um / 2.0 + dx, y_bot),
    ]


def build_slanted_groove_geometry(
    *,
    period_x_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    substrate_thickness_um: float,
    y_surface: float,
    medium_substrate: Any,
    medium_groove: Any = None,
) -> list:
    """金属基底 + 顶部 2D 斜槽 (单胞内一个非对称四边形空气区)。

    几何含义
    --------
    - 槽用 mp.Prism (任意多边形棱柱) 表示，axis 沿 z (out-of-plane)；
    - tilt_angle = 0 + top_width == bottom_width → 对称矩形槽（退化检验点，
      应给出与 `build_rectangular_groove_geometry` 完全相同的几何 → 光谱）；
    - tilt_angle ≠ 0 → 槽底相对槽顶沿 +x 偏移 depth·tan(α)，产生左右非对称
      的侧壁，正入射下也会破坏 ±x 对称性（这是后续研究定向发射的几何根源）。

    2D 等效模型 vs 真实 3D 斜孔
    ---------------------------
    本函数构建的是 “沿 z 无穷长 + xy 平面内斜壁四边形” 的 2D 等效结构：
      - 真实飞秒激光斜烧蚀产生的是 3D 倾斜圆/椭圆截面孔，有限深度、侧壁
        圆滑且常带氧化层 (TiO/TiO₂)；
      - 当前 2D 模型只保留了 “斜壁 + 周期性 + 一个特征宽度” 的最少必要要素，
        足以研究角度对正问题光谱的 **趋势**，但绝对量不可直接与 3D 实验对比；
      - 后续扩展路径：
          (a) 3D 仿真域 + mp.Prism / mp.Cone (沿斜轴) 替换本函数；
          (b) 在 prism 顶/底引入圆角 (改成 8 顶点多边形)；
          (c) 在槽壁内加一层薄 TiO₂ Medium (用第二个嵌套 prism)。

    单胞边界自检
    ------------
    要求所有顶点 x 落在 [-period_x_um/2, +period_x_um/2] 内，否则槽穿出单胞，
    Bloch 周期边界下相邻单胞的槽会相互重叠 / 错位，结果不可解释。
    """
    import meep as mp

    if period_x_um <= 0:
        raise ValueError(f"period_x_um 必须 > 0, 收到 {period_x_um}")
    if depth_um >= substrate_thickness_um:
        raise ValueError(
            f"depth_um ({depth_um}) 必须 < substrate_thickness_um "
            f"({substrate_thickness_um})"
        )
    if medium_groove is None:
        medium_groove = mp.Medium(epsilon=1.0)

    verts_xy = slanted_groove_vertices(
        top_width_um=top_width_um,
        bottom_width_um=bottom_width_um,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    half_P = period_x_um / 2.0
    out_of_cell = [(x, y) for (x, y) in verts_xy if abs(x) > half_P + 1e-9]
    if out_of_cell:
        raise ValueError(
            f"斜槽顶点越出单胞 [-P/2, +P/2]: P={period_x_um}, 越界顶点={out_of_cell}, "
            f"建议增大 period 或减小 top/bottom width / tilt_angle。"
        )

    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_x_um, substrate_thickness_um, mp.inf),
        ),
        mp.Prism(
            vertices=[mp.Vector3(x, y, 0) for (x, y) in verts_xy],
            height=mp.inf,
            axis=mp.Vector3(0, 0, 1),
            material=medium_groove,
        ),
    ]
    return geom


def build_oxidized_slanted_groove_geometry(
    *,
    period_x_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    substrate_thickness_um: float,
    y_surface: float,
    oxide_thickness_um: float,
    medium_substrate: Any,
    medium_oxide: Any,
    medium_groove: Any = None,
    oxide_mode: str = "conformal_approx",
) -> list:
    """Approximate oxidized 2D slanted-groove geometry.

    Modes
    -----
    ``top_film_only``
        Ti substrate plus air groove, with TiO2 blocks only on the top land
        regions outside the groove opening.

    ``conformal_approx``
        Ti substrate, an outer TiO2 groove-lining prism, and a smaller inner
        air prism.  This approximates oxide on sidewalls and the groove bottom
        without simply filling the whole groove with TiO2.  Top land regions
        also receive a TiO2 film.

    Notes
    -----
    The conformal mode is a geometric sensitivity model, not a true normal-
    offset surface mesh.  It intentionally keeps an air core in the groove.
    ``oxide_thickness_um == 0`` returns the bare slanted-groove geometry.
    """
    import meep as mp

    if oxide_thickness_um < 0:
        raise ValueError("oxide_thickness_um must be >= 0")
    if oxide_mode not in {"top_film_only", "conformal_approx"}:
        raise ValueError(
            "oxide_mode must be 'top_film_only' or 'conformal_approx', "
            f"got {oxide_mode!r}"
        )
    if medium_groove is None:
        medium_groove = mp.Medium(epsilon=1.0)
    if oxide_thickness_um == 0:
        return build_slanted_groove_geometry(
            period_x_um=period_x_um,
            top_width_um=top_width_um,
            bottom_width_um=bottom_width_um,
            depth_um=depth_um,
            tilt_angle_deg=tilt_angle_deg,
            substrate_thickness_um=substrate_thickness_um,
            y_surface=y_surface,
            medium_substrate=medium_substrate,
            medium_groove=medium_groove,
        )

    if oxide_thickness_um >= depth_um:
        raise ValueError(
            "oxide_thickness_um must be smaller than groove depth for the "
            "approximate conformal model."
        )
    if 2.0 * oxide_thickness_um >= min(top_width_um, bottom_width_um):
        raise ValueError(
            "oxide_thickness_um is too large for the current groove width; "
            "the inner air core would collapse."
        )

    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_x_um, substrate_thickness_um, mp.inf),
        )
    ]

    top_open_left = -top_width_um / 2.0
    top_open_right = top_width_um / 2.0
    half_p = period_x_um / 2.0
    left_land_width = max(0.0, top_open_left - (-half_p))
    right_land_width = max(0.0, half_p - top_open_right)
    if left_land_width > 0:
        geom.append(
            mp.Block(
                material=medium_oxide,
                center=mp.Vector3(
                    -half_p + 0.5 * left_land_width,
                    y_surface + 0.5 * oxide_thickness_um,
                    0,
                ),
                size=mp.Vector3(left_land_width, oxide_thickness_um, mp.inf),
            )
        )
    if right_land_width > 0:
        geom.append(
            mp.Block(
                material=medium_oxide,
                center=mp.Vector3(
                    top_open_right + 0.5 * right_land_width,
                    y_surface + 0.5 * oxide_thickness_um,
                    0,
                ),
                size=mp.Vector3(right_land_width, oxide_thickness_um, mp.inf),
            )
        )

    outer = slanted_groove_vertices(
        top_width_um=top_width_um,
        bottom_width_um=bottom_width_um,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    if oxide_mode == "top_film_only":
        geom.append(
            mp.Prism(
                vertices=[mp.Vector3(x, y, 0) for x, y in outer],
                height=mp.inf,
                axis=mp.Vector3(0, 0, 1),
                material=medium_groove,
            )
        )
        return geom

    geom.append(
        mp.Prism(
            vertices=[mp.Vector3(x, y, 0) for x, y in outer],
            height=mp.inf,
            axis=mp.Vector3(0, 0, 1),
            material=medium_oxide,
        )
    )
    inner = slanted_groove_vertices(
        top_width_um=top_width_um - 2.0 * oxide_thickness_um,
        bottom_width_um=bottom_width_um - 2.0 * oxide_thickness_um,
        depth_um=depth_um - oxide_thickness_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    geom.append(
        mp.Prism(
            vertices=[mp.Vector3(x, y, 0) for x, y in inner],
            height=mp.inf,
            axis=mp.Vector3(0, 0, 1),
            material=medium_groove,
        )
    )
    return geom


def build_inner_wall_film_slanted_groove_geometry(
    *,
    period_x_um: float,
    top_width_um: float,
    bottom_width_um: float,
    depth_um: float,
    tilt_angle_deg: float,
    substrate_thickness_um: float,
    y_surface: float,
    film_thickness_um: float,
    medium_substrate: Any,
    medium_film: Any,
    medium_groove: Any = None,
    coating_mode: str = "sidewalls_and_bottom",
) -> list:
    """Build a slanted groove with measured lossy film only inside the groove.

    The default ``sidewalls_and_bottom`` mode approximates the coating as an
    outer film prism occupying the original groove cavity plus a smaller inner
    air prism.  The inner air prism still opens at ``y_surface`` so the external
    flat top Ti land is not coated.  This is a geometric approximation, not a
    strict normal-offset constant-thickness surface mesh.

    ``sidewalls_only`` is a mechanism-control approximation: it creates a film
    shell on the groove sidewalls while the groove bottom remains air.  It is
    not used as the default quantitative conclusion.

    ``film_thickness_um == 0`` returns the same bare slanted-groove geometry as
    :func:`build_slanted_groove_geometry`.
    """
    import meep as mp

    if film_thickness_um < 0:
        raise ValueError("film_thickness_um must be >= 0")
    if coating_mode not in {"sidewalls_and_bottom", "sidewalls_only"}:
        raise ValueError(
            "coating_mode must be 'sidewalls_and_bottom' or 'sidewalls_only', "
            f"got {coating_mode!r}"
        )
    if medium_groove is None:
        medium_groove = mp.Medium(epsilon=1.0)
    if film_thickness_um == 0:
        return build_slanted_groove_geometry(
            period_x_um=period_x_um,
            top_width_um=top_width_um,
            bottom_width_um=bottom_width_um,
            depth_um=depth_um,
            tilt_angle_deg=tilt_angle_deg,
            substrate_thickness_um=substrate_thickness_um,
            y_surface=y_surface,
            medium_substrate=medium_substrate,
            medium_groove=medium_groove,
        )
    if film_thickness_um >= depth_um:
        raise ValueError("film_thickness_um must be smaller than depth_um")
    min_width = min(top_width_um, bottom_width_um)
    if 2.0 * film_thickness_um >= min_width:
        raise ValueError(
            "film_thickness_um is too large: the inner air core would collapse."
        )
    if depth_um >= substrate_thickness_um:
        raise ValueError("depth_um must be smaller than substrate_thickness_um")

    outer = slanted_groove_vertices(
        top_width_um=top_width_um,
        bottom_width_um=bottom_width_um,
        depth_um=depth_um,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    half_p = period_x_um / 2.0
    out_of_cell = [(x, y) for (x, y) in outer if abs(x) > half_p + 1e-9]
    if out_of_cell:
        raise ValueError(
            f"Outer groove vertices leave the unit cell: {out_of_cell}"
        )

    if coating_mode == "sidewalls_and_bottom":
        inner_depth = depth_um - film_thickness_um
        inner_bottom_width = bottom_width_um - 2.0 * film_thickness_um
    else:
        inner_depth = depth_um
        inner_bottom_width = bottom_width_um - 2.0 * film_thickness_um
    inner_top_width = top_width_um - 2.0 * film_thickness_um
    if inner_depth <= 0 or inner_top_width <= 0 or inner_bottom_width <= 0:
        raise ValueError("Inner air core collapsed; reduce film_thickness_um.")

    inner = slanted_groove_vertices(
        top_width_um=inner_top_width,
        bottom_width_um=inner_bottom_width,
        depth_um=inner_depth,
        tilt_angle_deg=tilt_angle_deg,
        y_surface=y_surface,
    )
    out_of_cell = [(x, y) for (x, y) in inner if abs(x) > half_p + 1e-9]
    if out_of_cell:
        raise ValueError(f"Inner groove vertices leave the unit cell: {out_of_cell}")

    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_x_um, substrate_thickness_um, mp.inf),
        ),
        mp.Prism(
            vertices=[mp.Vector3(x, y, 0) for x, y in outer],
            height=mp.inf,
            axis=mp.Vector3(0, 0, 1),
            material=medium_film,
        ),
        mp.Prism(
            vertices=[mp.Vector3(x, y, 0) for x, y in inner],
            height=mp.inf,
            axis=mp.Vector3(0, 0, 1),
            material=medium_groove,
        ),
    ]
    return geom


def build_planar_film_on_ti_geometry(
    *,
    period_x_um: float,
    substrate_thickness_um: float,
    film_thickness_um: float,
    y_surface: float,
    medium_substrate: Any,
    medium_film: Any,
) -> list:
    """Planar Ti backplane with optional full-surface measured lossy film.

    This geometry is a material-capability upper-bound test, not a microstructure
    model.  The Ti top surface is at ``y_surface``.  When
    ``film_thickness_um > 0``, the film uniformly covers the full x-period above
    the Ti surface.  ``film_thickness_um == 0`` strictly degenerates to bare
    planar Ti.
    """
    import meep as mp

    if period_x_um <= 0:
        raise ValueError(f"period_x_um must be > 0, got {period_x_um}")
    if substrate_thickness_um <= 0:
        raise ValueError(
            f"substrate_thickness_um must be > 0, got {substrate_thickness_um}"
        )
    if film_thickness_um < 0:
        raise ValueError(f"film_thickness_um must be >= 0, got {film_thickness_um}")

    geom = [
        mp.Block(
            material=medium_substrate,
            center=mp.Vector3(0, y_surface - substrate_thickness_um / 2.0, 0),
            size=mp.Vector3(period_x_um, substrate_thickness_um, mp.inf),
        )
    ]
    if film_thickness_um > 0:
        geom.append(
            mp.Block(
                material=medium_film,
                center=mp.Vector3(0, y_surface + film_thickness_um / 2.0, 0),
                size=mp.Vector3(period_x_um, film_thickness_um, mp.inf),
            )
        )
    return geom
