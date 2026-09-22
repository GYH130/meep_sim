# F：斜槽早期稳定性单变量诊断

本阶段只回答：保留相同材料、网格、时间步与几何时，关闭界面平滑是否消除已观察到的早期数值增长。它不是二维资格验证或尺寸扫描。

## 已知证据与冻结边界

- E 基准斜槽在 t≈78 的非 PML 电场强度达到约 1.26e302，随后溢出，状态 CRASHED；不是超时、内存限制或 API 故障。
- D 同材料 r32 平面反射率误差约 0.0020984，通过 0.0025 门槛；平面通过不代表斜槽通过。
- Ordal Route A 材料、baseline 几何、r32、Courant=0.25、p 偏振、+30°观察角、10.5 μm 保持不变。
- 不修改 ti2d；core SHA256 为 24d09c50dcf6c0cb9271947dbcb04f03b62d9d0500bd7bace34e6909d54834c1。
- 两臂依次为 eps_averaging=True / False。均省去被动 DFT/通量/体积分监视器和参考 DFT 载入；不使用生产缓存。因此保留平滑臂必须先复现早期增长，才能支持平滑相关归因。

## 预算和停止

- 最多两个任务，串行、各 4 MPI ranks、每 rank 单线程；单臂最多 1200 秒，两臂合计最多 1800 秒，包括退出等待。
- 使用原累计 48 小时账本；开始基数 30952.50267100334 秒。进程树内存限制为启动时容器可用内存的 70%。
- 仿真时间仅至 100；源约在 525 才结束，短测未发散不证明完整稳定或收敛。
- 每 2 时间单位记录非 PML |E|² 包络及局部探针。连续三个样本超过参考峰值的 1e12 倍且相邻强度增长均超过 10 倍，或发现非有限场值，立即提前停止。
- 非物理性的时间/内存/进程错误停止后续臂；不自动重试或扩展。

## 证据与命令

服务器根目录：`/root/autodl-tmp/meep_sim/meep_screen`。
结果：`repair_diagnostics/slot_stability_20260922_F/`。

- `plan.json` / `authorization.json` / `inputs/`：冻结计划、指纹、原输入及账本快照。
- `status.json`：短状态；`report.md` / `report.json`：结束后两臂结论及资源账本。
- 各臂 `trace.json`：包络、峰值位置、探针和停止过程；`peak_patch_*.npz`：最多四份小区域原始复场及几何掩膜；`sample_grid.npz`：坐标。
- `process.json` / `stdout.log`：资源、退出码和原始输出。所有结论均保持 qualification=false。

只读状态命令：

```bash
cd /root/autodl-tmp/meep_sim/meep_screen
cat repair_diagnostics/slot_stability_20260922_F/status.json
```

冻结和一次启动入口（已启动后不要再执行）：

```bash
/root/miniconda3/envs/meep_sim/bin/python -B run_slot_stability_diagnostic.py --prepare
/root/miniconda3/envs/meep_sim/bin/python -B run_slot_stability_diagnostic.py --run
```

监管器在服务器后台独立运行。AI 完成启动检查后结束本轮，不等待求解或反复读取日志。下一次只检查短汇总和有必要的证据；如果两臂均不稳定或原问题未复现，不宣称修复成功。即使关闭平滑臂短测稳定，仍需另行批准完整运行及平面/斜槽精度复核；不得直接开启 20 条件扫描。

不自动 Git 推送；旧结果和用户文件保持不变。字段中的“几何区域”是理想几何分类，不是有效介电张量的直接测量。
