# TA2 钛基微腔定向热发射项目评估入口

本页索引 `GYH130/meep_sim` 已公开的固定提交，供独立复核研究目标、物理仿真、几何筛选和代理模型。项目当前主波段为 **3–20 μm**，最高工作温度 **523 K**。15° 和 30° 是用于比较的角锥半角。目标角锥、最低发射率和加工窗口仍需确认。仓库默认分支根 README 记录早期 8–13 μm 的 Meep 研究。

| 评估主题 | 首读 | 数据或实现 |
| --- | --- | --- |
| 目标与边界 | [当前研究规格](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/objectives/specification.json) | [11 月 12 日增量说明](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/README.md) |
| Meep 与 COMSOL 基线 | [11 月 4 日报告和复现说明](https://github.com/GYH130/meep_sim/blob/1e890af1f4ded5e18d046a9dae47cfcc4a76cfb3/research_updates/20261104/README.md) | [88 个全波条件和质量标签](https://github.com/GYH130/meep_sim/blob/1e890af1f4ded5e18d046a9dae47cfcc4a76cfb3/research_updates/20261104/data/research_audit/current_archive_88_conditions.csv) · [补样结果](https://github.com/GYH130/meep_sim/blob/1e890af1f4ded5e18d046a9dae47cfcc4a76cfb3/research_updates/20261104/data/latest_ml/all9_scalar_results.csv) |
| 几何与角域 | [控制变量筛选及数量口径](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/geometry/screen200/report/metrics.json) | [结构族试点](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/geometry/family_pilot/report/family_summary.json) · [方法复核](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/geometry/method_review/assessment.json) |
| 代理模型 | [整套几何留出指标](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/surrogate/metrics.json) | [射线与 FEM 逐点对照](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/surrogate/inputs/ray_vs_fem.csv) · [WFT3 冻结预测](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/surrogate/WFT3_preregistered_predictions.json) |
| 最新公开全波锚点 | [WFT3 求解标量及警告](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/fem/WFT3_scalar_result.json) | [G99 曲壁与直壁对照的计划及状态](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/fem/g99/README.md) |
| 来源与复现 | [11 月 4 日来源清单](https://github.com/GYH130/meep_sim/blob/1e890af1f4ded5e18d046a9dae47cfcc4a76cfb3/research_updates/20261104/SOURCE_PROVENANCE.json) | [11 月 12 日归档清单](https://github.com/GYH130/meep_sim/blob/89e81a6129705751f9f16a795bb88e052cc3f6b3/research_updates/20261112/MANIFEST_SHA256.json) |

## 复核重点

1. 从同一条件比较全部传播衍射级的端口功率、材料损耗积分与原始净通量。将功率闭合、网格收敛、材料可信度分别评估。
2. 复算 88 个全波条件的几何分组和质量标签。把理想射线粗扫、细化与结构族试点各自的几何数和角度计算数列明。
3. 从三维角度和光谱覆盖审查锥内平均发射率、锥内功率占比以及每单位原胞面积的锥内功率。理想边界、光学钛/SiO₂ 仿真与实测各用独立证据。
4. 比较均值基线、直接高斯过程、射线和射线残差校正。检查整套结构留出的适用范围，并单独报告 WFT3 新几何的误差。
5. 给出最少的新增验证点、解决的分歧、验收依据和加工试验。优先说明哪些证据会改变路径选择。

请在评价中标注文件、条件和原始数值。提交链接固定到具体版本，便于复查。
