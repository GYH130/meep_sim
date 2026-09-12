# 合入原 meep_sim 项目的说明

这是文档与轻量证据包，不包含求解任务、缓存大文件、环境安装或生产代码替换。复制`results/`到原项目同名目录即可阅读。建议提交主题：`docs: document Ti validation progress and unresolved target-slot qualification`。

可以在原README中添加：

```markdown
## Ti验证进展

参见[2026-09-12对照报告](results/reports/ti_progress_20260912/REPORT.md)。
新增材料审计、拟合与稳定性诊断；目标斜槽尺寸扫描尚未通过数值资格。
```

## 后续代码移植映射（建议，不是本包已经实施的补丁）

| 当前工程模块 | 原项目建议位置 | 兼容检查 |
|---|---|---|
| `meep_screen/src/meep_screen/material_audit.py` | 新`src/material_validation.py` | 核对sigma与omega/f单位；禁止移除旧材料接口 |
| `ti_material.py`、`material_fit.py`及原始材料/拟合文件 | `src/materials.py`旁的新拟合模块与`data/` | 保留全部数据来源、温度未知和固定尺度误差定义；依赖须单独审查 |
| `blind_slot.py` | 新`src/blind_slot.py` | 本代码几何变量用x-z命名，但Meep二维内部通过x-y放置；原项目本来用x-y，不可机械重命名字段 |
| `material_audit.py`中的角度/contrast校验 | `src/postprocess.py`旁的新校验模块 | 发射角与入射kx反号；不重用旧angle字段而改变原意 |
| `recovery.py`、`safety.py`及workflow缓存逻辑 | 新运行管理模块 | 当前函数有包内依赖，不能只复制单文件；成功缓存须核验材料/设置/代码指纹 |
| `tests/test_sampled_decay.py`及几何/材料测试 | `tests/` | 合成停止序列测试不是FDTD收敛验收 |

避免整体覆盖`src/simulation.py`：旧项目已有D00/D04/D06能力，当前目标求解循环没有证明更快或整体更准。旧D06“选择最吻合方向”的做法应在单独补丁中改为预定义方向，并用已知传播方向对照验证，而不是只改一个阈值。

建议先创建文档PR，再创建模块级PR；每个代码PR需在原环境重新跑适配对照。不要自动运行`run_qualified_pilot.py`或任何扫描入口。不要把79项软件测试的通过写成物理认证通过。

## 包内容与复核

- `results/reports/ti_progress_20260912/REPORT.md`：可直接在GitHub阅读的中文对照报告。
- `evidence/*.json`：旧表格摘录、新结果/状态及原始文件指纹；保留失败记录，移除本机绝对路径。
- `evidence/manifest.json`：归档SHA-256、ZIP注释及各证据来源。原GitHub HEAD、许可证和冲突尚未核验。

本包不重新授权原数据/代码的许可证。发布前请按仓库现有许可及材料数据库数据来源要求复核；拟合记录中的来源说明是已保存元数据，不是法律意见。
