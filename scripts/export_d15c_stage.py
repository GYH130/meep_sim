#!/usr/bin/env python3
"""Generate a concise D15c status report from existing JSON/logs.

This script never imports or runs Meep. It only reads completed result files.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import subprocess
from typing import Optional


def git_value(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), *args], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unavailable"


def locate_table_dir(root: Path) -> Path:
    for candidate in (root / "results" / "tables", root / "tables",
                      root / "review" / "tables"):
        if any(candidate.glob("D15b_*.json")):
            return candidate
    raise FileNotFoundError("No D15b_*.json found under results/tables, tables, or review/tables")


def matching_log(root: Path, json_name: str) -> Optional[Path]:
    if json_name.startswith("D15b_structure_"):
        name = json_name.replace("D15b_structure_", "D15b_", 1).replace(".json", ".log")
    elif json_name.startswith("D15b_flat_flat_"):
        name = json_name.replace("D15b_flat_flat_", "D15b_flat_", 1).replace(".json", ".log")
    else:
        name = json_name.replace(".json", ".log")
    for candidate in (root / "logs" / name, root / "review" / "logs" / name):
        if candidate.is_file():
            return candidate
    return None


def fmt(value, digits=8) -> str:
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return "n/a"


def build_report(root: Path) -> str:
    table_dir = locate_table_dir(root)
    records = []
    for path in sorted(table_dir.glob("D15b_*.json")):
        data = json.loads(path.read_text())
        if data.get("max_steps"):
            continue
        log_path = matching_log(root, path.name)
        log_text = log_path.read_text(errors="replace") if log_path else ""
        hit_time_cap = "提前结束（未满足漂移判据）" in log_text
        audit_status = "TIME_LIMIT_UNQUALIFIED" if hit_time_cap else str(
            data.get("quality_status", "unknown"))
        records.append((path, data, log_path, audit_status))

    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    branch = git_value(root, "branch", "--show-current")
    commit = git_value(root, "rev-parse", "--short", "HEAD")
    lines = [
        "# D15c 当前结果简报",
        "",
        f"生成时间：{now}",
        f"Git：`{branch}` / `{commit}`",
        "",
        "## 结论",
        "",
        "现有文件证明 Meep 已实际运行并输出结果；当前只有一个斜槽单点的两个网格和一个平面对照。",
        "这些结果提示单点吸收增强，但分阶功率、体吸收分区和停止状态仍需修订，不能用于尺寸规律或正式方向性结论。",
        "",
        "## 已完成算例",
        "",
        "| 条件 | r | R | T | A_flux | 报告 A_vol | 步数 | 墙钟(s) | 审计状态 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for path, data, _log, audit_status in records:
        lines.append(
            "| {case} | {res} | {r} | {t} | {a} | {av} | {steps} | {wall} | {status} |".format(
                case=data.get("case", path.stem), res=data.get("resolution", "n/a"),
                r=fmt(data.get("R")), t=fmt(data.get("T")), a=fmt(data.get("A_flux")),
                av=fmt(data.get("A_vol")), steps=data.get("n_steps", "n/a"),
                wall=fmt(data.get("wall_total_s"), 7), status=audit_status,
            )
        )
    try:
        display_table_dir = table_dir.relative_to(root)
    except ValueError:
        display_table_dir = table_dir
    lines += [
        "",
        "## 资格状态",
        "",
        "- Meep 实际运行：已确认（时间步、墙钟、R/T/A 输出存在）。",
        "- 网格收敛：仅 r16/r20 两点，r24 未完成。",
        "- 停止条件：若日志触及 `max_time_um`，本报告标为 `TIME_LIMIT_UNQUALIFIED`。",
        "- 体吸收：原 A_vol 后处理坐标/权重在独立审计中发现问题，表中仅保留原值供追溯。",
        "- 分阶功率：原公式及入射/反射分离存在问题，旧的负阶份额不用于物理结论。",
        "- s 偏振、正负角、镜像及尺寸曲线：未完成。",
        "",
        "## 解释边界",
        "",
        "结果仅表示 specified Ti optical model 下的二维斜槽单点行为。`oxide_layer_in_model=false` 不代表真实样品没有氧化层。",
        "",
        f"数据目录：`{display_table_dir}`",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    text = build_report(root)
    if args.output:
        output = args.output if args.output.is_absolute() else root / args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_name(output.name + ".tmp")
        temporary.write_text(text + "\n")
        os.replace(temporary, output)
        print(output)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
