"""配置加载、日志、CSV / 图像输出。

约定
----
- 所有可调参数集中在 YAML 配置文件（configs/*.yaml）中，禁止散落代码；
- 每次脚本运行：
    1. 把使用的配置文件副本写到 results/reports/<run_tag>.yaml；
    2. 日志写入 logs/<run_tag>.log；
    3. CSV 输出到 data/processed/ 或 results/tables/；
    4. PNG 输出到 results/figures/。
- run_tag 形如  "<script_name>_<YYYYmmdd-HHMMSS>"，由 `make_run_tag` 生成。
"""

from __future__ import annotations

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Run tag / 路径辅助
# ---------------------------------------------------------------------------

def make_run_tag(script_name: str) -> str:
    """生成 '<script>_<YYYYmmdd-HHMMSS>' 风格的运行标签。"""
    stem = Path(script_name).stem
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return f"{stem}_{stamp}"


def ensure_dir(path: str | os.PathLike) -> Path:
    """确保目录存在，返回 Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def project_path(*parts: str) -> Path:
    """以项目根目录为基准拼路径。"""
    return PROJECT_ROOT.joinpath(*parts)


# ---------------------------------------------------------------------------
# 配置加载
# ---------------------------------------------------------------------------

def load_config(path: str | os.PathLike) -> dict[str, Any]:
    """加载 YAML / JSON 配置文件为 dict。"""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"配置文件不存在: {p}")
    suffix = p.suffix.lower()
    if suffix in (".yaml", ".yml"):
        import yaml  # 延迟导入，避免无 yaml 时影响其他工具
        with p.open("r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    if suffix == ".json":
        with p.open("r", encoding="utf-8") as f:
            return json.load(f)
    raise ValueError(f"不支持的配置文件扩展名: {suffix}")


def archive_config(config: dict[str, Any], run_tag: str) -> Path:
    """把本次运行使用的配置存档到 results/reports/<run_tag>.json。"""
    out = ensure_dir(project_path("results", "reports")) / f"{run_tag}.json"
    with out.open("w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)
    return out


# ---------------------------------------------------------------------------
# 日志
# ---------------------------------------------------------------------------

def setup_logger(run_tag: str, level: int = logging.INFO) -> logging.Logger:
    """同时输出到控制台和 logs/<run_tag>.log."""
    log_dir = ensure_dir(project_path("logs"))
    log_file = log_dir / f"{run_tag}.log"

    logger = logging.getLogger(run_tag)
    logger.setLevel(level)
    logger.propagate = False
    if logger.handlers:  # 避免重复添加
        return logger

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info("日志文件: %s", log_file)
    return logger


# ---------------------------------------------------------------------------
# CSV / 图像输出
# ---------------------------------------------------------------------------

def save_spectrum_csv(wavelengths_um: np.ndarray, columns: dict[str, np.ndarray],
                      out_path: str | os.PathLike) -> Path:
    """保存光谱表 (wavelength_um, col1, col2, ...) 到 CSV.

    所有 value 数组长度必须等于 wavelengths_um.
    """
    wavelengths_um = np.asarray(wavelengths_um)
    data = {"wavelength_um": wavelengths_um}
    for name, arr in columns.items():
        arr = np.asarray(arr)
        if arr.shape != wavelengths_um.shape:
            raise ValueError(
                f"列 '{name}' 长度 {arr.shape} 与波长轴 {wavelengths_um.shape} 不一致。"
            )
        data[name] = arr
    df = pd.DataFrame(data)
    p = Path(out_path)
    ensure_dir(p.parent)
    df.to_csv(p, index=False)
    return p


def save_figure(fig, out_path: str | os.PathLike, dpi: int = 200) -> Path:
    """统一保存 matplotlib 图：自动建目录、固定 dpi、tight bbox。"""
    p = Path(out_path)
    ensure_dir(p.parent)
    fig.savefig(p, dpi=dpi, bbox_inches="tight")
    return p
