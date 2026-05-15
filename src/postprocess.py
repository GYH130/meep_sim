"""Meep 输出 → 物理量（R, T, A, ε）的后处理。

职责
----
- 从 flux 监视器原始数据计算反射率 R、透射率 T、吸收率 A；
- 在 **不透明基底** 假设下，按 Kirchhoff 定律 ε(λ, θ) = A(λ, θ)
  得到光谱发射率与角分辨发射率；
- 提供自检函数（能量守恒、参考解对比等）；
- 输出整理后的 CSV 到 data/processed/，图到 results/figures/。

物理假设（必须随结果一同记录）
------------------------------
1. 局域热平衡 + 互易 → 可用吸收率代替发射率；
2. 金属层足够厚使得透射 T ≈ 0；若不满足，必须在结构下方加 flux 监视器
   并显式计算 T；
3. 非偏振发射率取 TE/TM 平均（这要求两次仿真）。
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


def compute_RTA_downward_incidence(
    reflection_flux_raw: np.ndarray,
    transmission_flux_raw: np.ndarray,
    input_flux_raw: np.ndarray,
) -> dict[str, np.ndarray]:
    """Compute R/T/A for this project's D00-validated downward incidence layout.

    Convention
    ----------
    - The incident plane wave propagates from +y toward -y.
    - The reflected wave propagates toward +y.
    - A horizontal Meep flux monitor reports positive flux along +y, so raw
      transmitted flux below an opaque/downward-illuminated structure is
      negative for physical downward power.

    Therefore:
        R = reflection_flux_raw / abs(input_flux_raw)
        raw_transmittance = transmission_flux_raw / abs(input_flux_raw)
        T = -transmission_flux_raw / abs(input_flux_raw)
        A = 1 - R - T
    """
    reflection_flux_raw = np.asarray(reflection_flux_raw, dtype=float)
    transmission_flux_raw = np.asarray(transmission_flux_raw, dtype=float)
    input_flux_raw = np.asarray(input_flux_raw, dtype=float)
    if not (
        reflection_flux_raw.shape
        == transmission_flux_raw.shape
        == input_flux_raw.shape
    ):
        raise ValueError(
            "flux 数组长度不一致："
            f"refl={reflection_flux_raw.shape}, "
            f"trans={transmission_flux_raw.shape}, "
            f"input={input_flux_raw.shape}"
        )
    if np.any(np.abs(input_flux_raw) == 0):
        raise ValueError("input_flux_raw 中存在 0，无法归一化。")

    denom = np.abs(input_flux_raw)
    R = reflection_flux_raw / denom
    raw_transmittance = transmission_flux_raw / denom
    signed_transmittance = -transmission_flux_raw / denom
    T = signed_transmittance
    A = 1.0 - R - T
    return {
        "R": R,
        "T": T,
        "A": A,
        "raw_transmittance": raw_transmittance,
        "signed_transmittance": signed_transmittance,
    }


def compute_RTA(refl_flux: np.ndarray, trans_flux: np.ndarray,
                ref_flux: np.ndarray) -> dict[str, np.ndarray]:
    """根据归一化 flux 计算 R, T, A（legacy 符号约定）。

    This function is retained for old scripts/tests that use the original
    convention: reflected flux is expected to be negative relative to a positive
    reference flux, hence ``R = -refl_flux / ref_flux`` and
    ``T = trans_flux / ref_flux``.

    New diagnostics for the current 2D Ti surface layout should call
    :func:`compute_RTA_downward_incidence`, whose sign convention was validated
    by D00 for source propagation from +y to -y.

    Parameters
    ----------
    refl_flux : np.ndarray
        结构存在时反射面记录的（已减去入射场参考的）反射 flux，shape (n_freq,)。
    trans_flux : np.ndarray
        透射面 flux，shape (n_freq,)。透明假设下应可置零。
    ref_flux : np.ndarray
        无结构（参考）仿真中入射面记录的入射 flux，shape (n_freq,)。

    Returns
    -------
    dict
        {"R": ..., "T": ..., "A": ...}，每个均为 (n_freq,) 数组。

    Raises
    ------
    ValueError
        三个数组长度不一致或 ref_flux 中存在 0。
    """
    refl_flux = np.asarray(refl_flux)
    trans_flux = np.asarray(trans_flux)
    ref_flux = np.asarray(ref_flux)
    if not (refl_flux.shape == trans_flux.shape == ref_flux.shape):
        raise ValueError(
            f"flux 数组长度不一致：refl={refl_flux.shape}, "
            f"trans={trans_flux.shape}, ref={ref_flux.shape}"
        )
    if np.any(ref_flux == 0):
        raise ValueError("ref_flux 中存在 0，无法归一化。")
    R = -refl_flux / ref_flux  # Meep 反射 flux 符号约定：减号
    T = trans_flux / ref_flux
    A = 1.0 - R - T
    return {"R": R, "T": T, "A": A}


def emissivity_from_absorptance(A: np.ndarray, *, opaque: bool = True) -> np.ndarray:
    """ε(λ, θ) = A(λ, θ)，前提是不透明基底 + 局域热平衡。

    Parameters
    ----------
    A : np.ndarray
        吸收率谱。
    opaque : bool
        若为 False，表示用户明确承认基底非不透明，仍允许通过，
        但调用者需要自行处理热辐射出射方向。
    """
    if not opaque:
        import warnings
        warnings.warn(
            "opaque=False：A ≈ ε 仅在不透明基底假设下严格成立，"
            "请确认结果解释方式。"
        )
    return np.asarray(A)


def energy_conservation_check(rta: dict[str, np.ndarray],
                              atol: float = 1e-3) -> dict[str, float]:
    """自检：R + T + A 是否 ≈ 1，并报告偏差统计。"""
    total = rta["R"] + rta["T"] + rta["A"]
    diff = total - 1.0
    report = {
        "max_abs_dev": float(np.max(np.abs(diff))),
        "mean_abs_dev": float(np.mean(np.abs(diff))),
    }
    if report["max_abs_dev"] > atol:
        import warnings
        warnings.warn(
            f"能量守恒偏差超过阈值 {atol}: max|R+T+A-1| = {report['max_abs_dev']:.3e}"
        )
    return report


def band_average(values: np.ndarray, wavelengths_um: np.ndarray,
                 lambda_lo_um: float, lambda_hi_um: float) -> float:
    """在指定波段对量（如 ε）做波长积分平均。

    旧版本采用简单算术平均；diagnostics_v2 起改为
    :func:`wavelength_integrated_average`，以避免等频率采样转换到不等间隔
    波长轴后产生偏置。
    """
    avg = wavelength_integrated_average(
        values, wavelengths_um, lambda_lo_um, lambda_hi_um,
    )
    if np.isnan(avg):
        raise ValueError(
            f"指定波段 [{lambda_lo_um}, {lambda_hi_um}] μm 内有效点少于 3 个。"
        )
    return avg


def wavelength_integrated_average(
    values: np.ndarray,
    wavelengths_um: np.ndarray,
    lambda_lo_um: float,
    lambda_hi_um: float,
    valid_mask: np.ndarray | None = None,
) -> float:
    """Return trapezoidal wavelength average over a band.

    The data are sorted by wavelength before integration.  Only points inside
    ``[lambda_lo_um, lambda_hi_um]`` and with ``valid_mask=True`` are used.
    If fewer than three valid samples remain, this function returns ``NaN`` and
    sets ``wavelength_integrated_average.last_status`` to
    ``"insufficient_samples"``.
    """
    wavelengths_um = np.asarray(wavelengths_um)
    values = np.asarray(values)
    if wavelengths_um.shape != values.shape:
        raise ValueError("wavelengths_um 与 values 长度不一致。")
    if lambda_lo_um >= lambda_hi_um:
        raise ValueError("lambda_lo_um 必须小于 lambda_hi_um。")
    if valid_mask is None:
        valid_mask = np.ones_like(wavelengths_um, dtype=bool)
    else:
        valid_mask = np.asarray(valid_mask, dtype=bool)
        if valid_mask.shape != wavelengths_um.shape:
            raise ValueError("valid_mask 与 wavelengths_um 长度不一致。")

    finite = np.isfinite(wavelengths_um) & np.isfinite(values)
    mask = (
        finite
        & valid_mask
        & (wavelengths_um >= lambda_lo_um)
        & (wavelengths_um <= lambda_hi_um)
    )
    if np.count_nonzero(mask) < 3:
        wavelength_integrated_average.last_status = "insufficient_samples"
        return float("nan")

    wl = wavelengths_um[mask]
    vals = values[mask]
    order = np.argsort(wl)
    wl = wl[order]
    vals = vals[order]
    span = float(wl[-1] - wl[0])
    if span <= 0:
        wavelength_integrated_average.last_status = "insufficient_span"
        return float("nan")
    wavelength_integrated_average.last_status = "valid"
    return float(np.trapezoid(vals, wl) / span)


wavelength_integrated_average.last_status = "not_called"


def material_validity_mask(
    wavelengths_um: np.ndarray,
    lambda_max_valid_um: float,
) -> np.ndarray:
    """Mask wavelengths within a material model's upper validity limit."""
    wavelengths_um = np.asarray(wavelengths_um, dtype=float)
    if lambda_max_valid_um <= 0:
        raise ValueError("lambda_max_valid_um must be positive.")
    return np.isfinite(wavelengths_um) & (wavelengths_um <= lambda_max_valid_um)


def broadband_source_quality_mask(
    input_flux_raw: np.ndarray,
    relative_threshold: float = 1e-3,
) -> dict[str, np.ndarray]:
    """Flag broadband source samples with weak incident power.

    ``source_relative_flux = abs(input_flux_raw) / max(abs(input_flux_raw))``.
    Points below ``relative_threshold`` are marked ``LOW_SOURCE_SNR`` and should
    be excluded from quantitative band metrics.
    """
    input_flux_raw = np.asarray(input_flux_raw, dtype=float)
    if relative_threshold < 0:
        raise ValueError("relative_threshold must be non-negative.")
    abs_flux = np.abs(input_flux_raw)
    max_flux = np.nanmax(abs_flux) if abs_flux.size else np.nan
    if not np.isfinite(max_flux) or max_flux <= 0:
        rel = np.full_like(abs_flux, np.nan, dtype=float)
        valid = np.zeros_like(abs_flux, dtype=bool)
    else:
        rel = abs_flux / max_flux
        valid = np.isfinite(rel) & (rel >= relative_threshold)
    flags = np.where(valid, "VALID", "LOW_SOURCE_SNR")
    return {
        "source_relative_flux": rel,
        "source_quality_flag": flags,
        "valid_for_quantitative_metric": valid,
    }


def opaque_substrate_transmission_check(
    T: np.ndarray,
    quantitative_tol: float = 1e-3,
    warning_tol: float = 5e-3,
) -> dict[str, np.ndarray]:
    """Assess transmission for the current thick opaque Ti substrate model.

    This check is not appropriate for transparent D00 lossless-interface tests.
    """
    if quantitative_tol < 0 or warning_tol < quantitative_tol:
        raise ValueError("Require 0 <= quantitative_tol <= warning_tol.")
    T = np.asarray(T, dtype=float)
    abs_T = np.abs(T)
    flags = np.full(T.shape, "FAIL", dtype=object)
    flags[abs_T <= warning_tol] = "WARNING"
    flags[abs_T <= quantitative_tol] = "NUMERICAL_PASS"
    valid = flags == "NUMERICAL_PASS"
    return {
        "abs_T": abs_T,
        "transmission_quality_flag": flags,
        "valid_for_opaque_substrate_metric": valid,
    }


def assess_quantitative_validity(
    wavelengths_um: np.ndarray,
    T: np.ndarray,
    input_flux_raw: np.ndarray,
    material_validity_mask: np.ndarray,
    source_mode: str,
) -> dict[str, np.ndarray]:
    """Combine material, source, and opaque-substrate quality flags per point."""
    wavelengths_um = np.asarray(wavelengths_um, dtype=float)
    T = np.asarray(T, dtype=float)
    input_flux_raw = np.asarray(input_flux_raw, dtype=float)
    material_validity_mask = np.asarray(material_validity_mask, dtype=bool)
    if not (
        wavelengths_um.shape
        == T.shape
        == input_flux_raw.shape
        == material_validity_mask.shape
    ):
        raise ValueError("wavelengths, T, input_flux, material mask shapes differ.")

    source = broadband_source_quality_mask(input_flux_raw)
    trans = opaque_substrate_transmission_check(T)
    issues: list[list[str]] = []
    flags = []
    for i in range(wavelengths_um.size):
        item_issues = []
        if not material_validity_mask[i]:
            item_issues.append("MATERIAL_EXTRAPOLATION")
        if source_mode != "single_wavelength" and not source["valid_for_quantitative_metric"][i]:
            item_issues.append("LOW_SOURCE_SNR")
        if trans["transmission_quality_flag"][i] == "FAIL":
            item_issues.append("TRANSMISSION_FAILURE")
        if len(item_issues) == 0:
            flags.append("VALID")
        elif len(item_issues) == 1:
            flags.append(item_issues[0])
        else:
            flags.append("MULTIPLE_WARNINGS")
        issues.append(item_issues)

    quantitative_valid = np.array([flag == "VALID" for flag in flags], dtype=bool)
    return {
        "quality_flag": np.asarray(flags, dtype=object),
        "valid_for_quantitative_metric": quantitative_valid,
        "source_relative_flux": source["source_relative_flux"],
        "source_quality_flag": source["source_quality_flag"],
        "transmission_quality_flag": trans["transmission_quality_flag"],
        "material_validity_flag": np.where(
            material_validity_mask, "VALID", "MATERIAL_EXTRAPOLATION",
        ),
        "quality_issues": np.asarray([";".join(x) for x in issues], dtype=object),
    }


def stack_angular_spectra(spectra_by_angle: dict[float, np.ndarray],
                          wavelengths_um: Iterable[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """把若干角度的光谱整理为 (theta, lambda) 二维矩阵，便于热图绘制。

    Parameters
    ----------
    spectra_by_angle : dict[float, np.ndarray]
        键为角度 (deg)，值为该角度下的 ε(λ) 数组。
    wavelengths_um : array-like
        所有角度共用的波长轴。

    Returns
    -------
    thetas, wavelengths, matrix
        matrix.shape == (n_theta, n_wavelength)
    """
    wavelengths = np.asarray(list(wavelengths_um))
    thetas = np.array(sorted(spectra_by_angle.keys()))
    matrix = np.empty((thetas.size, wavelengths.size), dtype=float)
    for i, th in enumerate(thetas):
        spec = np.asarray(spectra_by_angle[th])
        if spec.shape != wavelengths.shape:
            raise ValueError(
                f"角度 {th} 的光谱长度 {spec.shape} 与波长轴 {wavelengths.shape} 不一致。"
            )
        matrix[i] = spec
    return thetas, wavelengths, matrix
