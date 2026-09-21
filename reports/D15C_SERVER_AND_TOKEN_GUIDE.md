# D15c 服务器结果、GitHub 保存与低 token 操作说明

## 当前判断

服务器确实运行了 Meep。现有求解脚本导入 `meep`、建立 `mp.Simulation` 并调用 `sim.run`；日志记录了 92,800/116,000 个时间步、墙钟时间及 R/T/A。

Meep 与 Claude Code 是两个独立过程。Meep 在服务器 CPU 上运行，不调用 Claude API。token 费用来自 Claude Code 在每次搜索、命令、日志读取和判断时重新调用模型。缓存是正常的提示词缓存机制，但本次上下文增长到几十万 token，并被 1,471 次请求反复读取，因此产生高额费用。

现有数据只支持一个 10.5 μm、+30°、p 偏振斜槽单点的初步吸收增强线索。分阶功率、体吸收分区和停止状态仍需修订，不能作为尺寸规律、方向性或 FWHM 结论。

## 在 AutoDL 服务器生成当前报告

先把 `export_d15c_stage.py` 上传到服务器仓库的 `scripts/` 目录，然后在终端运行：

```bash
cd /root/autodl-tmp/meep_sim
source /root/miniconda3/etc/profile.d/conda.sh
conda activate meep_sim

python -c "import meep as mp; print('Meep version:', mp.__version__)"
python scripts/D15c_build_d15c_tables.py | tee results/reports/D15c_table_rebuild.log
python scripts/export_d15c_stage.py \
  --project-root . \
  --output results/reports/D15c_current_status.md

sed -n '1,220p' results/reports/D15c_current_status.md
```

上述命令只重建表格、读取已有 JSON/日志并生成报告，不启动新的 FDTD。

如暂时不上传新脚本，可直接查看现有报告：

```bash
cd /root/autodl-tmp/meep_sim
sed -n '1,260p' review/D15c_stage_report.md
sed -n '1,220p' review/MANIFEST.md
```

建议再生成一个下载包：

```bash
cd /root/autodl-tmp/meep_sim
mkdir -p results/exports
tar -czf results/exports/D15c_current_report.tar.gz \
  results/reports/D15c_current_status.md \
  results/reports/D15c_table_rebuild.log \
  results/tables/D15b_*.json \
  results/tables/D15c_*.csv \
  review/D15c_stage_report.md \
  review/MANIFEST.md
sha256sum results/exports/D15c_current_report.tar.gz
```

## 保存到 GitHub

先确认当前分支和准备提交的文件：

```bash
cd /root/autodl-tmp/meep_sim
git status --short --branch
git branch --show-current

git add \
  scripts/export_d15c_stage.py \
  results/reports/D15c_current_status.md \
  results/reports/D15c_table_rebuild.log \
  results/tables/D15b_*.json \
  results/tables/D15c_*.csv \
  review/D15c_stage_report.md \
  review/MANIFEST.md

git diff --cached --check
git diff --cached --stat
git diff --cached --name-only
```

确认清单没有 API 密钥、凭据或超大 HDF5 文件后提交：

```bash
git commit -m "Add D15c Meep results and qualification report"
git push -u origin "$(git branch --show-current)"
```

不要提交 `results/refplanes/*.h5`、模型密钥、`.env` 或 Claude 配置。若 GitHub 要求认证，使用 SSH key 或 GitHub fine-grained token；不要把 token 写进远程 URL、脚本或聊天。

## 后续低 token 规则

1. 每个阶段开一个新会话，只提供 `HANDOFF.md` 和必要文件路径。
2. 固定 Sonnet 处理编码、命令和日志；Opus 只用于一次设计审查和最终复核。
3. Claude 只负责生成并启动脚本，随后结束会话。Meep 用 `screen` 或 `nohup` 独立运行。
4. 禁止 Claude 轮询日志。求解器只写 `status.json`、结果 JSON 和最终日志尾部。
5. 每轮限制最多 8–12 次工具调用，完成既定产物后停止。
6. 不把完整日志、历史聊天或大段代码反复输入模型；只读取错误行、摘要和指定函数。
7. 在任何新仿真前先输出条件数、单例时长、预计 CPU-h 和费用上限。

推荐给 Claude Code 的固定提示：

> 本轮只处理明确列出的文件和任务，最多 10 次工具调用。禁止轮询后台进程，禁止启动未批准的仿真，禁止读取完整日志；只读取最后 100 行和错误匹配。生成规定产物后立即停止并给出文件清单。若计算仍在运行，写明进程号和结果路径后退出，不等待、不持续检查。
