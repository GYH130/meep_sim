# 钛基微腔定向热发射：截至 2026-11-04 的研究数据与代码

本补充包汇总 MEEP 二维探索、COMSOL 三维空腔/SiO₂ 底填充研究、最新 9 项机器学习补样、功率独立复核、代理模型验证与加工路线调研。目标为 TA2/钛基激光微结构在 8–13 μm 的发射率及角域调控；现有单点或稀疏数据不能代替完整带内/半球积分。

## 阅读顺序与证据边界

1. `docs/`：本轮总报告（如有）与关键图表说明。
2. `data/history/`：2026-10-14/30 历史科学汇总，含结果、失败/警告和实际覆盖范围。不同阶段材料、几何、方法不同，不能混成一条性能曲线。
3. `data/latest_ml/all9_scalar_results.csv` 与同目录 JSON：最新 9 项物理解和保存场功率方法复核。它们对应 10.5 μm、H80/H100 与 SiO₂ 底填充厚度 2/3 μm 的受控补样；不是 9 个独立尺寸。
4. `data/ml/dataset.json`：89 条档案记录、去重后 86 条物理条件、8 套独立 H/t 设计。`baseline_rows` 为 26 个单点代理候选条件。保留原始净通量警告和质量标签，按整套结构留出，避免将角度/偏振当新增独立尺寸。
5. `data/fabrication/`：可加工性、SiO₂ 材料/填充方法的文献与几何分析，属于建议，不是已加工或实测结果。

总研究档案与 ML 池口径不同：`data/research_audit/` 补入两项此前未进入 ML 池的 H100/t3、9 μm、法向 s/p 警告结果，因此研究底填充档案为 **91 条求解记录、88 个去重条件、8 套 H/t 设计**，其中保存场方法通过37、失败14、未审核37。ML 池仍为89/86，未因整理报告自动重训或静默加入失败点。各阶段报告重叠，不可相加成“全项目总数”。

## 本次阶段报告

- [Word 可编辑报告](docs/TA2钛基微结构定向热发射_阶段研究报告_20261104.docx)
- [PDF 阅读版](docs/TA2钛基微结构定向热发射_阶段研究报告_20261104.pdf)

15页、11幅原创图、11个可编辑公式；包含研究目标、已做工作、结果边界、工艺建议和下一步计划。

**状态必须分开：** 9/9 有实际求解且完成独立标量复核；9/9 保存场功率方法复核通过；原始净通量接受标记为 0/9；全部 `formal_qualification=false`，尚未建立整体网格收敛。严格 0.002 与探索绝对 0.005、相对 2% 均保留。`A_loss_total`、`A_port`、`A_flux` 分开保留，不用“修正值”掩盖差异。

## 代码与复现

- 本包不启动远程任务，不含 SSH、凭据、进程/资源快照、大型 MPH 或三维场文件。
- `src/ml/reproduce_validation.py`：可直接在本包根目录运行的本地代理复核，输出到 `reproduced/`。它按 H/t 整套结构留出；均值、Ridge 与 Matérn GP 分别比较，保留质量池。运行：`python -m pip install -r requirements-analysis.txt`，然后 `python src/ml/reproduce_validation.py`。
- `src/ml/model_helpers.py`：原有模型函数的冻结科学源代码；直接执行其原始 `main()` 仍要求旧工作区目录，不属于本包入口。请使用上述便携入口。
- `src/comsol/`：实际三维几何、材料与功率重建的冻结 Java/Python 科学核心，需 COMSOL Multiphysics/Wave Optics 6.3 授权。`CavityFloorFilm.java` 被原控制器调用，单独编译不等于重启全部作业。原控制器/主机准入与资源策略不随公开包发布；`data/latest_ml/frozen_scientific_plan.json` 保存模型设置与源绑定约定。某些冻结源码旧注释含 48 GiB 缓存描述，后续实际作业使用 24 GiB，缓存不改变物理材料/几何。公开包未包含原始场，因此独立保存场重建需要从原档案取得对应模型/导出字段。
- `src/meep/ti2d/`：2026-10-13 稳定二维桥接版本核心；材料在 `assets/material/`。这是二维历史模型，不能直接复用为后续三维空腔结果。`environment-ti.yml` 是历史单进程配方；当时 MPI 环境与并行参数必须另行一致核对。
- `figures/`：原创科学示意、真实数据图的 PNG 与可编辑 SVG。线段在稀疏点间仅作视觉引导，不表明连续角域或连续光谱已求解。
- `src/figures/build_figures.py` 和 `build_history_figure.py` 为公开路径适配入口，输出到 `reproduced/teacher_figures/`；`original_*.py` 保留原文件。中文字体优先使用 macOS Arial Unicode，其他系统需安装可用 CJK 字体。历史覆盖面板读取冻结 ML 池，须与 `research_audit` 的完整88条件台账区分。
- `src/figures/v2_build_*.py` 与 `fabrication_build_figures.py` 重绘报告使用的既有图组，输出到 `reproduced/v2_figures/`。`src/figures/build_report.py` 使用包内图和数据生成新Word文件到 `reproduced/report/final/`，不会覆盖交付版；依赖 python-docx 及可用中文字体。PDF需另经Word/LibreOffice渲染并逐页核对；生成Word不等于完成渲染验收。冻结报告源的个人根路径已用 `<WORKSPACE_ROOT>` 代替，原始SHA仍在来源清单。

ML 目前用于小样本代理可用性检验与下一步补点设计，尚未证明优于简单均值基线，不能称“已找到最优结构”。正常/斜向吸收通过互易条件映射到相应外向发射方向，两个方向不代表完整角域。

## 数据来源、单位与引用

长度 μm；角度 degree；R/A 为入射功率的无量纲比值；“百分点”=比例差×100。光学钛为文献 n,k 的阻抗边界，并非 TA2 实测光学常数；SiO₂ 为 Kischkat 文献数据的底填充，不是沿侧壁共形镀膜。原始文献出处、材料表与复核指纹在各 JSON 中。

`SOURCE_PROVENANCE.json` 记录每个归档源文件的工作区相对定位符及原始 SHA256，和公开导出后 SHA256。`workspace/`、`server_archive/` 是归档定位符，不是本包可打开的路径；公开数据仅规范化个人路径与非科学主机元数据，未改变数值。`MANIFEST_SHA256.json` 包含公开全部文件哈希（自身除外）。历史根仓库文件不被本补充目录覆盖。未另行推定原始代码或材料数据的新许可证；复用应遵循根仓库和原文献/数据库条件。
