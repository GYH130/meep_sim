# J：有预算上限的尾场边界诊断

## 当前证据与问题

I 原始计算达到 2.5 小时上限，仿真时间 914.625；仍为 TIME_LIMIT_UNQUALIFIED。
最后一次额外采样间隔为 0.125，其余为 0.25。离线恢复仅将该终点从 FFT 排除，
全部原始文件、状态和旧报告保留，新报告在 I/reports_recovered/。
I 的背面主要空间分量为 m=7；其截止频率约 0.09238/µm，与尾场谱峰在现有频率分辨率下接近。
这是待验证的近掠射尾场线索，不是已证明的 PML 原因或合格光学结果。

## 本轮锁定范围

- 仅两项串行诊断：pml6_control、pml12_test；4 MPI rank、每 rank 单线程。
- 同一冻结 Ordal Route A、二维纯 Ti、无氧化层、r32/C0.25、原激励、原探针及孔槽。
- 仿真终点 t1200；对照最高 4h，PML12 最高 5h；本轮最高 9h，原累计 48h 不变。
- 依据实测 I 8955s/t914.625 估计对照约 3.27h；PML12 网格高度增加 12.5%，估算约 3.68h。
  估计不是完成承诺，超时不加时、不自动重试。
- I 无可恢复的场状态，因此 J 从零场启动，不称作续算。
- 第二项只在第一项完成且与 H 比较通过原 >=t1000 门槛后启动。
- 唯一物理变量是 PML 厚度 6→12µm；物理区域、源、槽和监视面坐标不变。
- 不扩展尺寸扫描；所有结果 qualification=false；仅诊断边界厚度敏感性。
- 今天 16:52 的租赁截止假设已撤销。具体明天到期时刻未知，不将等待窗口当作租赁到期。

## 修复及运行方式

新 probe_v2 将终点样本单独保存并校验，周期采样数据维持等间隔。
恢复分析拒绝内部缺点、非有限数值、文件损坏或不同的配对采样时刻；不插值、不放宽门槛。
计划及依赖指纹在启动前冻结；使用共享预算锁，退出时确认所属 MPI 子进程清理。
现有 Git 改动及相关依赖备份到新运行 inputs/，不改旧分支、不关旧 Claude。

服务器项目中的入口（使用已有 meep_sim 环境）：

```sh
python -B -m unittest -v test_recover_tail_analysis test_tail_j
python -B recover_tail_analysis.py runs/tail_20260923_I
python -B run_tail_j.py prepare
python -B run_tail_j.py run
```

上述 prepare/run 是单次入口；已启动后不要重复运行。
状态：runs/tail_20260923_J/status.json；结果：reports_recovered/。
原始压缩包：exports/J_raw_evidence.tar.gz；完成清单：exports/J_export.json。
本机 collect_J_once.py 等待完成后校验下载原始证据，再推送脱敏精简快照到既有专用归档分支。
Mac 保持唤醒和联网才能自动下载/上传；服务器队列不依赖本机连接，后台不调用 AI。
AI 实际 token 与人民币账单无法获取，不编造成本；计算资源由共享账本记录。
