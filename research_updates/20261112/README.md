# TA2 钛基微腔研究迁移补充包 · 2026-11-12

这是 11 月 4 日研究发布后的增量。此前的 Meep、COMSOL、ML 数据和报告保存在 [2026-11-04 发布](../20261104/README.md)（已核对原发布提交 `1e890af1f4ded5e18d046a9dae47cfcc4a76cfb3`）。本包可在本地读取和复算。冻结的原环境、服务器任务记录、原生 COMSOL 模型及大规模场文件应从私有迁移副本恢复。

| 目录 | 内容及状态 |
| --- | --- |
| `objectives` | 主波段 3–20 μm，温度最高 523 K，参考膜层、工艺可行性和角域目标。规格文件保留待定事项。 |
| `geometry/angular_dense` | 20 个理想几何的粗角域数据、6 个细扫几何及角域分析脚本。 |
| `geometry/screen200` | 200 个理想几何及 24 个细扫几何的射线筛选，含原胞尺寸、面形、每角采样结果和配对选择依据。 |
| `geometry/family_pilot` | 三组预设表面条件，各 24 个几何、10,392 个角度条件。数据用于检查面形与表面吸收的敏感性。 |
| `geometry/method_review` | 角度权重、接收轴定义、角锥积分、理想表面假设及其适用范围。 |
| `surrogate` | 8 个独立 H/b/t 几何的 88 个历史条件、整套几何留出验证、两组额外尺寸和 WFT3 冻结预测。可从源代码重训。 |
| `fem` | WFT3 的已保存单点吸收与审核字段，以及 G99 的两种面形参数、物理模型、功率检查和新的部署模板。G99 此处提供输入与方法。 |

## 复算现有数据

进入本目录。用于图表复算的 `numpy`、`matplotlib` 可以单独安装。选择明确的新输出路径以保留已有图。

```sh
python -m pip install numpy matplotlib
python geometry/angular_dense/analyze_angular.py \
  --coarse geometry/angular_dense/fetched/coarse_summary.json \
  --refine geometry/angular_dense/fetched/refine_summary.json \
  --out geometry/angular_dense/replayed
python geometry/family_pilot/analyze_results.py \
  geometry/family_pilot/fetched/ideal_summary.json \
  geometry/family_pilot/fetched/floor90_summary.json \
  geometry/family_pilot/fetched/wall04_summary.json \
  --output geometry/family_pilot/angle_metrics_replayed.json
python geometry/family_pilot/summarize_family_pilot.py \
  geometry/family_pilot/angle_metrics_replayed.json \
  --output geometry/family_pilot/family_summary_replayed.json
```

`analyze_results.py` 还可读取 `screen200/fetched` 和 `screen200/prior` 的四份角点结果，然后将输出传给 `method_review/rescored/summarize_angles.py`，复绘角锥灵敏度图。聚合时保留粗扫和细扫来源。中文重绘用 `MEEP_CJK_FONT` 指定本机 CJK 字体文件。理想黑底和镜面壁的射线模型忽略波动效应。预设灰体表面用于敏感性比较，其输出不直接给出真实 TA2 或膜层的 3–20 μm 光谱。采用 `cos(theta) sin(theta)` 投影权重，可按指定半角和接收轴积分。

射线核使用 C++17 编译。三个几何目录各有输入 `manifest.json`、`run_pilot.py` 及 `resource_guard.py`。可用 `validate-plan` 检查案例和核的哈希。重新运行原 Windows/WSL 后台保护时，在新主机设置 `MEEP_RAY_STAGE_ROOT` 和 `MEEP_WINDOWS_POWERSHELL`，核实新磁盘和资源；使用 `commands` 查看生成命令。保持每份冻结几何输入、采样数、方向、种子和权重一致后再对比。

## 代理模型

在 `surrogate/` 下，复制目录到独立工作位置再重训，以保留原输出及冻结预测。使用能满足 `requirements.txt` 的 Python 环境。

```sh
python -m pip install -r requirements.txt
c++ -O3 -std=c++17 ray_engine/ray_material.cpp -o ray_engine/ray_material
python train.py
python evaluate.py
python predict.py --help
```

预测目标为给定波长、方向和偏振的 `A_loss`。训练宽度为 15 μm，整套结构留出用于检验局部插值。WFT1/WFT2 单点保持在训练集外，WFT3 的三种预测冻结在其实际计算结果出现之前。WFT3 的端口法吸收率为 62.646%、损耗积分为 62.727%，原净通量为 67.150%；原始净通量警告保留，保存场功率核算未获验证，正式合格标记为 false。

## G99 在新服务器恢复

见 [G99 模型与迁移说明](fem/g99/README.md)。保留基准和 8 面分段形的源文件、几何、材料及对照输入。新主机需要 Windows、持有效许可的 COMSOL 6.3 环境、同批两个网格和独立的功率检查。GitHub 包里没有可直接续算的 COMSOL 检查点，当前服务器的真实任务与结果由私有迁移资料记录。

本目录的 `MANIFEST_SHA256.json` 记录每个文件的发布哈希与大小、原始项目来源和可移植处理。哈希不同的文本经过路径配置、脱敏或 JSON 紧凑序列化。历史摘要里用于标识原试验输入的 SHA-256 仍指原件。
