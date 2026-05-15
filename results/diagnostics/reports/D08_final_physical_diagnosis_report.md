# D08 Final Physical Diagnosis Report

## Executive Diagnosis

D00 已通过，但 D03 显示数值未完全收敛，结构定量结论暂缓。

## Stage Status

| Stage | Status | Summary |
|---|---|---|
| D00 Flux Sign And Fresnel Validation | **PASS** | 透射符号选择为 ['signed']；无损界面最大 R/T 误差 1.81e-05/0.000676；平面 Ti 最大 A 误差 0.00166。 |
| D01 Solver Regression | **PASS** | 回归表状态 PASS；通量修正前后平均 A 最大变化 0.000262。 |
| D02 Polarization Diagnostic | **PASS** | Hz-Ez 的 8-13 um 平均 A 差值范围 [-2.03e-07, 0.0152]；未发现 Hz 明显高于 Ez。 |
| D03 Numerical Convergence | **FAIL** | 最大谱差 0.0834；失败对象 ['slanted_Hz']。 |
| D04 Field And Absorbed Power | **WARNING** | Flux/volume 最大差 0.000225；最大 metal |E|^2 增强 0.154，最大槽壁耗散占比 0.676。 |
| D05 Angle-Resolved Absorptance | **PASS** | 最大角平均 A 约 0.0949；最大方向比 1.03；最大 A(+θ)-A(-θ) 0.00242。 |
| D06 Diffraction Energy Channels | **PASS** | 阶功率求和残差最大 1.02e-08；平均 R=0.865，A=0.135，非镜面反射=0.453。 |
| D07 Oxide And Multiscale Sensitivity | **PASS** | 最大 placeholder 氧化层增强 0.0159；未超过显著阈值。 TiO2 仍为 placeholder，不能作实验定量预测。 |

## Confirmed Conclusions

- D00 确认了透射通量应使用 signed transmittance，平面 Ti 与 Fresnel 基准基本一致。
- D02 未显示 Hz 相对 Ez 的强烈、稳健跃迁，固定 Ez 结论不是唯一问题来源。
- D04 支持当前裸 Ti 几何没有形成强局域场增强。
- D05 当前结果未显示强方向性；角分辨结论仍受 smoke/采样范围限制。
- D06 的衍射阶功率和总反射 flux 闭合，说明反射能量通道分析可用。
- D07 placeholder smoke 未显示超过阈值的氧化层增强。

## Unconfirmed Hypotheses

- 裸 Ti 斜槽提升有限可能来自缺少有效局域损耗模式，也可能来自 2D 理想几何过于简单。
- 真实飞秒加工样品的高发射若存在，可能依赖氧化层、粗糙度、再凝固层或多尺度结构。
- 斜槽几何更可能负责方向性/衍射通道调制，而不是单独提供宽带高吸收。

## Next Computation To Spend Resources On

- 优先重跑 D03 中失败的 slanted_Hz resolution 收敛，用 resolution>=48/64 确定统一基准。
- 用 D03 推荐的收敛设置重跑 D02/D05/D06 的关键点，避免 smoke 参数主导结论。
- 替换 D07 的 TiO2 placeholder，加入可靠 mid-IR TiO2/TiOx n,k 或 Drude-Lorentz 拟合。
- 在氧化层模型通过后，再做几何参数扫描和数据集生成。

## Next Experimental Data To Collect

- 氧化层厚度和空间分布，尤其槽壁/槽底覆盖。
- TiO2/TiOx 成分、相态和 8-13 um 光学常数。
- 截面几何：槽深、开口、侧壁角、圆角、周期分布。
- 多尺度粗糙度和再凝固颗粒，用于决定是否需要 3D/有效介质模型。

## Machine-Learning Dataset Readiness

**不建议进入机器学习数据集生成阶段。**

D03 存在失败项，当前不应将这些结果作为机器学习训练标签。

## Important Caveats

- D08 不重新运行 Meep，只汇总已存在 CSV；若上游 D00-D07 是 smoke test，本报告也只是 smoke 级别。
- 如果 D00 或 D03 未通过，不应把当前结构吸收率作为机器学习标签。
- D07 当前 TiO2 是 `placeholder_demo`，不能作为实验定量预测。

## Outputs

- Dashboard: `/Users/luckydog/meep_sim/results/diagnostics/figures/D08_diagnosis_dashboard.png`
- Report: `/Users/luckydog/meep_sim/results/diagnostics/reports/D08_final_physical_diagnosis_report.md`

## Thresholds Used

```json
{
  "hz_delta_tol": 0.02,
  "local_e2_tol": 1.0,
  "low_absorption_threshold": 0.2,
  "directionality_ratio_tol": 1.2,
  "oxide_enhancement_tol": 0.03
}
```
