# 钛斜槽定向热发射计算阶段结果

截至 2026 年 10 月 14 日，同一稳定 Ti 光学模型下已获得 10 个二维斜盲槽探索性结果点。当前最明确的结论是：在已采样的三个空气开口占空比中，F = 0.65 的 10.5 μm、p 偏振、+30° 吸收响应及 ±30° 两角差最大；这不是全局最优，也不是完整带内方向性结论。

研究目标是 8–13 μm 波段的方向性热发射，并在后续拓展至 3–20 μm。当前数据仍是离散波长、离散角度和 p 偏振下的吸收率计算，不代表已完成连续发射率光谱或真实 TA2 样品验证。

## 报告与图集

- [可编辑 Word 报告](report/钛斜槽定向热发射计算阶段报告_20261014.docx)
- [PDF 报告](report/钛斜槽定向热发射计算阶段报告_20261014.pdf)
- [完整图集及图注](FIGURE_GALLERY.md)：6 张当前结果与诊断图、10 张历史诊断与几何图。
- [历史图目录与原图 SHA256](data/existing_figure_catalog.json)

## 当前关键结果

共同条件：二维周期斜盲槽，周期 50 μm，槽轴长 40 μm，倾角 30°，Ti 厚度 60 μm，分辨率 r = 20。F 表示表面空气开口宽度与周期之比，不是金属体积分数。

以下为 10.5 μm、p 偏振的通量吸收率 A_flux；两角差使用百分点，不是相对百分比。

|空气开口占空比 F|A_flux at +30°|A_flux at −30°|两角差|
|---|---:|---:|---:|
|0.50|13.328%|8.510%|4.818 个百分点|
|0.65|19.974%|8.935%|11.039 个百分点|
|0.80|16.104%|8.289%|7.815 个百分点|

独立体吸收 A_vol 给出相同的三个采样点趋势。F = 0.65 处的两角差为 11.038 个百分点，与通量计算接近；两种吸收方法的差异不是完整误差条。−30° 各占空比之间的小差别应谨慎解释。

F = 0.65 的 p 偏振已覆盖 8、10.5、13 μm 三个离散波长。在这三个点，A_flux 的 ±30° 两角差依次为 11.417、11.039、8.557 个百分点。三点不构成连续 8–13 μm 光谱，不能据此计算可靠的带内平均发射率。

![当前占空比筛选](figures/current/01_fill_comparison.png)

## 数值状态和适用边界

- 10 个结构结果保留为探索性结果，不升级为正式资格通过。单点通量与体吸收差约为 0.27–0.38 个百分点，原始小负透射值保持不裁剪。
- 同材料平面控制中，12 项有 4 项通过对应严格平面精度门槛，8 项保留为精度未通过的诊断。完成求解不等于所有波长和偏振均已验证。
- 正入射平面控制显示，同一种电源或磁源口径下 s/p 一致；电源与磁源的 A_flux 相差约 0.371 个百分点。这是诊断范围内的源口径差异，不是可施加于全部斜槽结果的通用修正量。
- 当前尚不能给出完整角谱、峰值角、半高宽、带内方向性积分、非偏振性能或三维圆孔的定量性能。
- 改变 F 时，最深槽底和剩余底厚也随几何改变；本组扫描不是固定垂直深度的单独孔径试验。
- 旧材料模型、旧停止判据和软件示例图均与当前结果分开，不能混合成一条参数曲线。

## 数据文件

|文件|内容|
|---|---|
|[current_structured_results.json](data/current_structured_results.json)|10 个当前结构点的 R、T、A_flux、A_vol、状态、检查项和证据指纹|
|[direction_pairs.json](data/direction_pairs.json)|5 组 ±30° 配对差值与比值；不是完整角域方向性因子|
|[geometry_and_settings.json](data/geometry_and_settings.json)|三种占空比的几何与计算设置|
|[current_model.json](data/current_model.json)|当前 Ti 色散模型及来源信息|
|[material_fit_summary.json](data/material_fit_summary.json)|材料频域拟合摘要及范围限制|
|[flat_validation_results.json](data/flat_validation_results.json)|12 项平面控制结果，保留未通过项|
|[source_type_diagnostic.json](data/source_type_diagnostic.json)|正入射源口径控制与保存场诊断|
|[stability_diagnosis_summary.json](data/stability_diagnosis_summary.json)|被替代宽带拟合的早期发散诊断，与当前稳定模型区分|
|[existing_figure_catalog.json](data/existing_figure_catalog.json)|10 张历史原图的图注、用途、来源标识和 SHA256|

R、T、A 的 JSON 数值均为无量纲比例；绘图中的 `%` 表示乘以 100，`pp` 表示百分点。负 T 和未通过状态有意保留，不能通过裁剪或重命名消除。

## 重绘图表

安装 Python 3、NumPy 和 Matplotlib 后，在本目录运行：

```bash
python scripts/plot_results.py
```

该命令只读取随包 JSON 并生成 6 张当前图的 PNG 与 SVG，不调用 Meep，也不启动计算。

[build_report.py](scripts/build_report.py) 使用这些 JSON 和图表重建 DOCX，需要 python-docx、Pillow 及中文字体。PDF 导出和逐页排版检查是独立步骤。

[extract_scientific_data.py](scripts/extract_scientific_data.py) 保留了从原始工作区证据摘要抽取公开数据的逻辑；运行它需要未随包分发的原始证据文件。本包支持结果查看与图表重绘，**不是包含全部原始场、缓存、运行环境和求解输入的完整 FDTD 复现包**。

## 文件范围

公开内容包括科学指标、模型参数、图表、报告和相应证据指纹。不包含远程连接凭据、私人机器标识、完整运行日志、大体积场文件或虚拟磁盘。证据指纹用于识别既有来源，不等同于原始证据已公开分发。
