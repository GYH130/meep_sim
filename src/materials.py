"""金属色散模型与单位换算。

职责
----
- 提供 8–13 μm 中红外波段常用金属（Au、Ag、Cu、W、Ni、Ti 等）的色散模型；
- 提供 Drude / Lorentz–Drude 参数 → `meep.Medium` 的转换；
- 提供波长 ↔ Meep 频率换算（基于 a = 1 μm）。

单位约定
--------
- 真空光速 c = 1（Meep 内部）；
- 特征长度 a = 1 μm；
- 波长 λ [μm] 与 Meep 频率 f 的换算：
      f = a / λ = 1 / λ[μm]
- 角频率 ω = 2π f。

色散数据来源建议
----------------
- Rakić et al., "Optical properties of metallic films for vertical-cavity
  optoelectronic devices", Appl. Opt. 37, 5271 (1998).
- Palik, "Handbook of Optical Constants of Solids".
具体参数加入时请在函数 docstring 中注明引用与适用波长范围。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# 单位换算
# ---------------------------------------------------------------------------

def wavelength_um_to_meep_freq(wavelength_um: float) -> float:
    """λ[μm] → Meep 频率 (a = 1 μm)."""
    if wavelength_um <= 0:
        raise ValueError(f"wavelength_um 必须为正数，收到 {wavelength_um}")
    return 1.0 / wavelength_um


def meep_freq_to_wavelength_um(freq_meep: float) -> float:
    """Meep 频率 → λ[μm] (a = 1 μm)."""
    if freq_meep <= 0:
        raise ValueError(f"freq_meep 必须为正数，收到 {freq_meep}")
    return 1.0 / freq_meep


def freq_range_for_band(lambda_min_um: float, lambda_max_um: float) -> tuple[float, float, float]:
    """给定波长范围返回 (f_min, f_max, f_center) 用于宽带源。"""
    if lambda_min_um <= 0 or lambda_max_um <= 0:
        raise ValueError("波长必须为正数。")
    if lambda_min_um >= lambda_max_um:
        raise ValueError("lambda_min_um 必须小于 lambda_max_um。")
    f_max = wavelength_um_to_meep_freq(lambda_min_um)  # 短波长 → 高频
    f_min = wavelength_um_to_meep_freq(lambda_max_um)
    f_center = 0.5 * (f_min + f_max)
    return f_min, f_max, f_center


# ---------------------------------------------------------------------------
# 色散参数容器
# ---------------------------------------------------------------------------

@dataclass
class DrudeLorentzParams:
    """Drude–Lorentz 色散参数（Meep 内部单位）。

    所有频率均为 Meep 单位（a = 1 μm 时即 1/μm）。

    Attributes
    ----------
    name : str
        材料名（例如 "Au_Rakic"）；用于日志与文件名。
    epsilon_inf : float
        高频介电常数 ε∞。
    drude_omega_p : float
        Drude 项等效等离子体频率（Meep 单位）。
    drude_gamma : float
        Drude 项阻尼。
    lorentz_terms : list[dict]
        每一项 = {"sigma": ..., "omega_0": ..., "gamma": ...}，
        Meep `LorentzianSusceptibility` 直接参数化。
    valid_lambda_um : tuple[float, float] | None
        参数标定有效波长范围 (μm)，仅用于自检；超出会发警告。
    reference : str
        文献引用。
    """

    name: str
    epsilon_inf: float
    drude_omega_p: float
    drude_gamma: float
    lorentz_terms: list[dict] = field(default_factory=list)
    valid_lambda_um: tuple[float, float] | None = None
    reference: str = ""


# ---------------------------------------------------------------------------
# 材料库（占位；在后续脚本中按需补全数值）
# ---------------------------------------------------------------------------

MATERIAL_LIBRARY: dict[str, DrudeLorentzParams] = {}
"""注册的材料字典，键为材料名。

在后续 Step 中通过 `register_material` 加入具体金属参数（如 Au_Rakic 等）。
"""


def register_material(params: DrudeLorentzParams) -> None:
    """将一个 DrudeLorentzParams 注册到全局材料库。"""
    if params.name in MATERIAL_LIBRARY:
        raise ValueError(f"材料 {params.name} 已注册；如需替换请先删除。")
    MATERIAL_LIBRARY[params.name] = params


def build_medium(name: str) -> Any:
    """根据材料名构造 `meep.Medium` 对象。

    Parameters
    ----------
    name : str
        在 MATERIAL_LIBRARY 中的键。

    Returns
    -------
    meep.Medium

    Notes
    -----
    Step 0 骨架：实现延后到加入第一个具体材料时。
    届时在此处 `import meep as mp` 并组装 Drude + Lorentz 项。
    """
    if name not in MATERIAL_LIBRARY:
        raise KeyError(f"未注册的材料：{name}。已注册：{list(MATERIAL_LIBRARY)}")
    raise NotImplementedError("Step 0 骨架：实现延后到第一个具体材料脚本。")


def check_band_within_validity(name: str, lambda_min_um: float, lambda_max_um: float) -> None:
    """自检：所选波段是否在材料色散模型的标定有效范围内。"""
    params = MATERIAL_LIBRARY.get(name)
    if params is None or params.valid_lambda_um is None:
        return
    lo, hi = params.valid_lambda_um
    if lambda_min_um < lo or lambda_max_um > hi:
        import warnings
        warnings.warn(
            f"材料 {name} 标定范围 [{lo}, {hi}] μm，"
            f"当前研究波段 [{lambda_min_um}, {lambda_max_um}] μm 超出，结果需谨慎。"
        )


# ---------------------------------------------------------------------------
# 直接复用 meep.materials 的内置金属
# ---------------------------------------------------------------------------
# 设计取舍：Meep 自带 `meep.materials` 提供了多种金属的 Drude–Lorentz 拟合
# （Rakić 1998, Appl. Opt. 37, 5271）。在没有更可靠的本地实验数据之前，
# 把这些封装成 get_xxx_medium() 函数，统一进出口，便于：
#   1) 在脚本里只 import 一个名字，不必每个脚本都 `from meep.materials import ...`；
#   2) 在此处做波段有效性自检并打日志；
#   3) 将来若改用自家拟合参数，只需替换函数体，调用方不动。


# Rakić 1998 Drude–Lorentz Ti 模型的标定波长范围（μm），与 meep.materials.metal_range 一致。
TI_RAKIC_VALID_LAMBDA_UM: tuple[float, float] = (0.24797, 12.398)
TIO2_PLACEHOLDER_VALID_LAMBDA_UM: tuple[float, float] | None = None
"""No quantitative validity range is claimed for the TiO2 placeholder model."""
MEASURED_LOSSY_WALL_FILM_SOURCE_FILE = "data/raw/measured_lossy_wall_film_nk.xlsx"


def get_ti_medium(lambda_min_um: float | None = None,
                  lambda_max_um: float | None = None):
    """返回 Meep 自带的 Ti 介质（Rakić 1998 Drude–Lorentz 拟合）。

    Parameters
    ----------
    lambda_min_um, lambda_max_um : float, optional
        若提供，则对研究波段做一次有效性自检；超出 0.248–12.4 μm
        部分（属于 Drude 尾部外推）会发出 warning 但不会阻止仿真。

    Returns
    -------
    meep.Medium

    Notes
    -----
    - 该色散模型在中红外大气窗口 8–13 μm 内仅 12.4 μm 以下严格有效，
      13 μm 附近是 Drude 尾部外推，定性可信但定量需以 FTIR 实验校核；
    - 拟合不含温度依赖、不含氧化层 (TiO/TiO₂)，对实际飞秒激光改性的
      Ti 表面只是“干净金属基底”理想化；
    - 各向同性 (sigma_diag 三向相等)，所以 2D 计算中 TE/TM 在正入射时
      给出相同 R/T。
    """
    try:
        from meep.materials import Ti as _Ti  # 延迟导入：让 src 不强依赖 Meep
    except ImportError as exc:
        raise ImportError(
            "未能导入 meep.materials.Ti；请确认已激活 meep_env 或安装 pymeep。"
        ) from exc

    if lambda_min_um is not None and lambda_max_um is not None:
        lo_valid, hi_valid = TI_RAKIC_VALID_LAMBDA_UM
        if lambda_min_um < lo_valid or lambda_max_um > hi_valid:
            import warnings
            warnings.warn(
                f"Ti (Rakić 1998) 标定范围 [{lo_valid}, {hi_valid}] μm，"
                f"当前研究波段 [{lambda_min_um}, {lambda_max_um}] μm 超出，"
                f"超出部分为 Drude 尾部外推，结果定性可信，定量请用 FTIR 校核。"
            )
    return _Ti


def get_tio2_medium(
    lambda_min_um: float | None = None,
    lambda_max_um: float | None = None,
    *,
    model: str = "placeholder_demo",
    demo_index: float = 2.4,
    allow_placeholder: bool = True,
):
    """Return a replaceable TiO2 medium interface.

    This project does not yet include a vetted mid-infrared TiO2 optical-constant
    dataset or Drude-Lorentz fit.  The default ``placeholder_demo`` mode is a
    **lossless, non-dispersive sensitivity-test medium** with user-controlled
    refractive index.  It is useful for checking whether adding a dielectric
    layer could plausibly matter, but it must not be used for final quantitative
    comparison to experiment.

    Parameters
    ----------
    lambda_min_um, lambda_max_um : float, optional
        Accepted for API symmetry and warning messages.  No quantitative
        validity is claimed for the placeholder.
    model : {"placeholder_demo"}
        Future TiO2 fits should be added as new model names without changing
        callers.
    demo_index : float
        Real refractive index used by the placeholder demo medium.
    allow_placeholder : bool
        If False, raise an error instead of returning the placeholder.
    """
    import warnings

    if model != "placeholder_demo":
        raise ValueError(
            f"Unknown TiO2 model {model!r}. Available: 'placeholder_demo'."
        )
    if not allow_placeholder:
        raise ValueError(
            "No vetted mid-IR TiO2 model is available in this project yet. "
            "Pass allow_placeholder=True only for sensitivity tests."
        )
    if demo_index <= 0:
        raise ValueError(f"demo_index must be positive, got {demo_index}")

    warnings.warn(
        "Using TiO2 placeholder_demo: lossless, non-dispersive n="
        f"{demo_index:g}. This is only for sensitivity testing and is not a "
        "quantitative TiO2 mid-IR optical model."
    )
    try:
        import meep as mp
    except ImportError as exc:
        raise ImportError(
            "未能导入 meep；请确认已激活 meep_env 或安装 pymeep。"
        ) from exc
    return mp.Medium(index=demo_index)


def load_measured_nk_table(csv_path: str | Path):
    """Load a cleaned measured n,k table for ``measured_lossy_wall_film``.

    The CSV must contain ``wavelength_um``, ``n``, and ``k`` columns.  The table
    is sorted by wavelength and validated for monotonic wavelength, finite
    positive ``n``, and non-negative ``k``.  This function intentionally does
    not infer chemistry from the filename; the material is treated only as the
    user-provided measured lossy wall film.
    """
    import numpy as np
    import pandas as pd

    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Measured n,k CSV not found: {path}")
    df = pd.read_csv(path)
    required = ["wavelength_um", "n", "k"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Measured n,k CSV is missing columns: {missing}")
    df = df[required].copy()
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df.isna().any().any():
        raise ValueError("Measured n,k CSV contains missing/non-numeric values.")
    df = df.sort_values("wavelength_um").reset_index(drop=True)
    if not (df["wavelength_um"] > 0).all():
        raise ValueError("Measured n,k wavelengths must be positive.")
    if not (df["n"] > 0).all():
        raise ValueError("Measured n values must be positive.")
    if not (df["k"] >= 0).all():
        raise ValueError("Measured k values must be non-negative.")
    wl = df["wavelength_um"].to_numpy(dtype=float)
    if np.any(np.diff(wl) <= 0):
        raise ValueError("Measured n,k wavelengths must be strictly increasing.")
    return df


def interpolate_nk_at_wavelength(
    nk_table,
    wavelength_um: float,
    allow_extrapolation: bool = False,
) -> dict:
    """Linearly interpolate measured n,k at one wavelength.

    Extrapolation is forbidden by default.  The returned complex permittivity is
    computed as ``epsilon_complex = (n + 1j*k)**2``.
    """
    import numpy as np

    if wavelength_um <= 0:
        raise ValueError(f"wavelength_um must be positive, got {wavelength_um}")
    df = nk_table.copy()
    for col in ["wavelength_um", "n", "k"]:
        if col not in df.columns:
            raise ValueError(f"nk_table is missing required column {col!r}")
    wl = df["wavelength_um"].to_numpy(dtype=float)
    n_arr = df["n"].to_numpy(dtype=float)
    k_arr = df["k"].to_numpy(dtype=float)
    data_min = float(np.min(wl))
    data_max = float(np.max(wl))
    if not allow_extrapolation and not (data_min <= wavelength_um <= data_max):
        raise ValueError(
            f"wavelength_um={wavelength_um} is outside measured n,k range "
            f"[{data_min}, {data_max}] um; extrapolation is disabled."
        )
    n = float(np.interp(wavelength_um, wl, n_arr))
    k = float(np.interp(wavelength_um, wl, k_arr))
    epsilon_complex = (n + 1j * k) ** 2
    return {
        "wavelength_um": float(wavelength_um),
        "n": n,
        "k": k,
        "epsilon_real": float(np.real(epsilon_complex)),
        "epsilon_imag": float(np.imag(epsilon_complex)),
        "data_lambda_min_um": data_min,
        "data_lambda_max_um": data_max,
        "interpolation_flag": (
            "INTERPOLATED_WITHIN_DATA_RANGE"
            if data_min <= wavelength_um <= data_max
            else "EXTRAPOLATED_OUTSIDE_DATA_RANGE"
        ),
    }


def get_measured_lossy_wall_film_medium_single_wavelength(
    wavelength_um: float,
    nk_csv_path: str | Path,
    allow_extrapolation: bool = False,
):
    """Return a narrowband Meep medium from user-provided measured n,k data.

    This interface is valid only for one wavelength at a time.  It maps the
    measured complex permittivity to a non-dispersive Meep medium with
    ``D_conductivity`` at the target frequency.  It is not a broadband
    dispersive fit and must not be used for a single broadband FDTD run.
    """
    import math
    import meep as mp

    table = load_measured_nk_table(nk_csv_path)
    nk = interpolate_nk_at_wavelength(
        table,
        wavelength_um,
        allow_extrapolation=allow_extrapolation,
    )
    epsilon_real = nk["epsilon_real"]
    epsilon_imag = nk["epsilon_imag"]
    if epsilon_real <= 0:
        raise ValueError(
            "Conductivity medium requires positive epsilon_real; "
            f"got {epsilon_real} at {wavelength_um} um."
        )
    frequency_meep = 1.0 / wavelength_um
    d_conductivity = 2.0 * math.pi * frequency_meep * epsilon_imag / epsilon_real
    medium = mp.Medium(epsilon=epsilon_real, D_conductivity=d_conductivity)
    metadata = {
        "material_name": "measured_lossy_wall_film",
        "model_mode": "nk_interpolated_single_wavelength_conductivity",
        "wavelength_um": float(wavelength_um),
        "n": nk["n"],
        "k": nk["k"],
        "epsilon_real": epsilon_real,
        "epsilon_imag": epsilon_imag,
        "D_conductivity": float(d_conductivity),
        "quantitative_ready": True,
        "data_range_um": [nk["data_lambda_min_um"], nk["data_lambda_max_um"]],
        "source_file": MEASURED_LOSSY_WALL_FILM_SOURCE_FILE,
        "interpolation_flag": nk["interpolation_flag"],
        "warning": (
            "Valid only for narrowband single-wavelength simulation; not a "
            "broadband dispersive fit."
        ),
    }
    return medium, metadata


def fit_measured_nk_to_lorentz_medium(*args, **kwargs):
    """Placeholder for a future broadband dispersive fit.

    If a future workflow needs a single broadband FDTD calculation, the measured
    n,k data must first be fitted to a causal Lorentz or Drude-Lorentz
    dispersive model.  The current implementation uses one-wavelength
    conductivity media and is not suitable for broadband dispersive simulation.
    """
    raise NotImplementedError(
        "Broadband Lorentz/Drude-Lorentz fitting for measured_lossy_wall_film "
        "has not been implemented."
    )
